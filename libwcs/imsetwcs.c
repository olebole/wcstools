/* File libwcs/imsetwcs.c
 * February 21, 1996
 * By Elwood Downey, revised by Doug Mink
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fitshead.h"
#include "wcs.h"


#define MINSTARS	4	/* min stars we require from GSC and image */
#define MINBIN		4	/* min coincident hits required */
#define MAXSTARS	25	/* max star pairs we need to try using */
#define MAXGSC		100	/* max Guide Stars to find in image region */

extern int findStars ();
extern int findRegistration ();
extern int findCoords ();

static int getNominalPos ();
static void sky2im ();
static void setFITSWCS ();
void sortStars ();
static void sortGuideStars ();
static int starstatSortF ();

/* set the C* WCS fields in header based on the gsc limiting values and dimmer.
 * do it by finding stars in image and in GSC and finding the rotation and
 *   offsets which result in a best-fit.
 * verbose generates extra info on stdout.
 * try using deeper gsc searches if we have trouble.
 * return 1 if all ok, else 0
 */

static int tolerance = 20;	/* +/- this many pixels is a hit */
void settolerance (tol)
int tol;
{
    tolerance = tol;
    return;
}

static double gsclim = 17;	/* initial GSC search limiting magnitude */
void setgsclim (lim)
double lim;
{
    gsclim = lim;
    return;
}

static double rot0 = 0.0;	/* Initial image rotation */
void setrot (rot)
double rot;
{
    rot0 = rot;
    return;
}

static double secpix0 = 0.0;	/* Set image scale--override image header */
void setsecpix (secpix)
double secpix;
{
    secpix0 = secpix;
    return;
}

static double flip = -1.0;	/* Flip around N-S axis */
void setflip ()
{
    flip = 1.0;
    return;
}

static double frac = 1.0;	/* Additional catalog or image stars to use */
void setfrac (frac0)
double frac0;
{
    frac = frac0;
    return;
}

static int trimatch = 0;	/* Star matching method */
void setmatch ()
{
    trimatch = 1;
    return;
}


/* set the C* WCS fields in fip based on the given limiting GSC mag.
 * do it by finding stars in fip and in GSC down to gsclim and finding the
 *   angle and offsets which result in a best-fit.
 * verbose generates extra info on stdout.
 * return 0 if all ok, else -1
 */

int
setWCSFITS (header, image, verbose)

char	*header;	/* FITS header */
char	*image;		/* Image pixels */
int	verbose;

{
    double *gnum=0;		/* GSC star numbers */
    double *gra=0, *gdec=0;	/* GSC stars, ra and dec, rads */
    double *gm=0;		/* GCS magnitudes */
    int ng;			/* n GSC stars */
    int nbg;		/* n brighttest GCS stars we actually use */
    double *sx=0, *sy=0;	/* image stars, pixels */
    double *sb=0;		/* image star brightesses */
    int ns;			/* n image stars */
    int nbs;		/* n brightest image stars we actually use */
    double *gx, *gy;	/* GSC stars projected onto plane, pixels */
    double ra0, dec0;	/* nominal center from RA/DEC FITS fields */
    double pixsz;		/* pixel size, rads */
    int imw, imh;		/* image size, pixels */
    double dra, ddec;	/* errors, in rads */
    double rot, dx, dy;
    double ra1,ra2,dec1,dec2,mag1,mag2;
    double secpix;
    double x0, y0;
    int ngmax;
    int classd;
    int nbin, nbytes;
    int ret = 0;
    int i;
    int gscread();

    /* get nominal position and scale */
    if (getNominalPos (header,verbose,&ra0,&dec0,&pixsz,&imw,&imh) < 0)
	goto out;

    /* Set the RA and Dec limits in degrees for GSCREAD */
    ddec = pixsz * imh * 0.48;
    dra = (pixsz * imw * 0.48) / cos (dec0);
    ra1 = raddeg (ra0 - dra);
    ra2 = raddeg (ra0 + dra);
    dec1 = raddeg (dec0 - ddec);
    dec2 = raddeg (dec0 + ddec);
    mag1 = -2.0;
    mag2 = gsclim;

    ngmax = MAXGSC;
    nbytes = MAXGSC * sizeof (double);
    gnum = (double *) malloc (nbytes);
    gra = (double *) malloc (nbytes);
    gdec = (double *) malloc (nbytes);
    gm = (double *) malloc (nbytes);
    gx = (double *) malloc (nbytes);
    gy = (double *) malloc (nbytes);

    /* Find the nearby GSC stars, in ra/dec */
    ng = gscread (ra1,ra2,dec1,dec2,mag1,mag2,classd,ngmax,gnum,gra,gdec,gm);
    if (ng < MINSTARS) {
	if (ng < 0)
	fprintf (stderr, "Error getting GSC stars: %d\n", ng);
	else if (ng == 0)
	fprintf (stderr,"Need at least %d GSC stars but found none\n",
							MINSTARS);
	else
	fprintf (stderr, "Need at least %d GSC stars but only found %d\n",
							MINSTARS, ng);
	goto out;
	}

    /* Convert ra and dec to radians */
    for (i = 0; i < ng; i++) {
	gra[i] = degrad (gra[i]);
	gdec[i] = degrad (gdec[i]);
	}

    /* project the GSC stars onto a plane nominally at ra0/dec0 to
     * convert them into pixels
     */
    if (!gx || !gy) {
	fprintf (stderr, "Could not malloc temp space of %d bytes\n",
					    ng*sizeof(double)*2);
	goto out;
	}
    sky2im (imw, imh, pixsz, ra0, dec0, ng, gra, gdec, gx, gy);

    /* Sort Guide Stars by brightness */
    sortGuideStars (gnum, gra, gdec, gx, gy, gm, ng);

    /* need only use the brightest MAXSTARS Guide Stars */
    if (ng > MAXSTARS) {
	nbg = MAXSTARS;
	if (verbose)
	    printf ("using %d / %d GSC stars brighter than %.1f\n",
		    gm[nbg-1], nbg,ng);
	}
    else {
	nbg = ng;
	if (verbose) {
	    if (gsclim > 0.0)
		printf ("using all %d GSC stars brighter than %.1f\n",
			ng,gsclim);
	    else
		printf ("using all %d GSC stars\n", ng);
	    }
	}

    for (i = 0; i < ng; i++) {
	char rstr[64], dstr[64];
	ra2str (rstr, raddeg (gra[i]), 3);
	dec2str (dstr, raddeg (gdec[i]), 2);
	printf ("GS: %9.4f %s %s %5.2f %6.1f %6.1f\n",
		gnum[i],rstr,dstr,gm[i],gx[i],gy[i]);
	}

    /* Discover star-like things in the image, in pixels */
    ns = findStars (header, image, &sx, &sy, &sb);
    if (ns < MINSTARS) {
	fprintf (stderr, "Need at least %d image stars but only found %d\n",
							MINSTARS, ns);
	goto out;
	}

    /* Sort star-like objects in image by brightness */
    sortStars (sx, sy, sb, ns);

    /* Use only as many star-like objects as Guide Stars */
    if (ns > nbg) {
	nbs = nbg * frac;
	if (nbs > ns)
	   nbs = ns;
	if (verbose) {
	    printf ("using brightest %d / %d GSC stars\n", nbg, ng);
	    printf ("using brightest %d / %d image stars\n", nbs,ns);
	    }
	}
    else {
	nbs = ns;
	nbg = ng * frac;
	if (nbg > ng)
	    nbg = ng;
	if (verbose)
	    printf ("using all %d image stars\n", ns);
	}

    if (verbose) {
	for (i = 0; i < nbs; i++)
	printf ("Im: %6.1f %6.1f %8.1f\n", sx[i],sy[i],sb[i]);
	}

    /* Find offset, scale, and rotation to im x,y which best matches ref x,y */
    rot = rot0;
    dx = 0.0;
    dy = 0.0;

    /* Find offset, scale, and rotation to image x,y which best matches
	reference ra,dec */
    if (trimatch) {
	/* secpix = 3600.0 * raddeg (pixsz); */
	secpix = 0.0;
	nbin = findCoords (ra0, dec0, nbg, gra, gdec, gm, gx, gy, nbs, sx, sy,
			   &secpix, &rot);
	pixsz = degrad (secpix / 3600.0);
	}

    /* Find offset and rotation to image x,y which best matches reference x,y */
    else {
	nbin = findRegistration (sx, sy, nbs, gx, gy, nbg, imw, imh, MINBIN,
				 tolerance, &rot, &dx, &dy);
	rot = -rot;
	}
    if (nbin < 0) {
	fprintf (stderr, "Star registration failed.\n");
	goto out;
	}
    else if (nbin < MINBIN) {
	fprintf (stderr, "Require %d bin hits but only found %d\n", MINBIN,nbin);
	goto out;
	}
    else if (verbose)
	printf ("%d bin hits\n", nbin);

    dra = dx * pixsz / cos(dec0);
    ddec = dy * pixsz;
    ra0 = ra0 - dra;
    dec0 = dec0 - ddec;

    if (verbose) {
	char str[64];

	ra2str (str, raddeg (dra), 3);
	printf ("Delta  RA: %s = %4g pixels\n", str, dx);
	dec2str (str, raddeg (ddec), 2);
	printf ("Delta Dec: %s = %4g pixels\n", str, dy);
	printf ("Rotation:  %g degrees\n", raddeg(rot));
	}

    /* compute and set (or replace) WCS fields in header */
    x0 = (double) imw / 2.0;
    y0 = (double) imh / 2.0;
    setFITSWCS (header, rot, ra0, dec0, x0, y0, pixsz);

    ret = 1;

    out:
    if (gra)
	free ((char *)gra);
    if (gdec)
	free ((char *)gdec);
    if (gm)
	free ((char *)gm);
    if (gnum)
	free ((char *)gnum);
    if (gx)
	free ((char *)gx);
    if (gy)
	free ((char *)gy);

    if (sx)
	free ((char *)sx);
    if (sy)
	free ((char *)sy);
    if (sb)
	free ((char *)sb);

    return (ret);
}


/* find nominal image center and scale of given image from local header info.
 * If center does not use FK5 (J2000) equinox, convert it
 * return 0 if ok else return -1
 */

static int
getNominalPos (header, verbose, rap, decp, radperpixp, wp, hp)

char	*header;
int	verbose;
double	*rap;
double	*decp;
double	*radperpixp;
int	*wp;
int	*hp;
{
    int nax;
    double secpix, dra, ddec, eq;

    /* find nominal center from RA and DEC fields */
    dra = 0.0;
    if (hgetra (header, "RA", &dra) == 0) {
	fprintf (stderr, "No RA field\n");
	return (-1);
	}
    ddec = 0.0;
    if (hgetdec (header, "DEC", &ddec) == 0) {
	fprintf (stderr, "No DEC field\n");
	return (-1);
	}

    /* Equinox of coordinates */
    if (hgetr8 (header, "EPOCH", &eq) == 0) {
	if (hgetr8 (header, "EQUINOX", &eq) == 0)
	    eq = 1950.0;
	}

    /* If not J2000, convert */
    if (eq == 1950.0)
	fk425 (&dra, &ddec);

    *rap = degrad (dra);
    *decp = degrad (ddec);

    /* Plate scale from SECPIX header parameter */
    if (secpix0 == 0.0) {
	secpix = 0.5;
	if (hgetr8 (header, "SECPIX", &secpix) == 0) {
	    fprintf (stderr, "Cannot find SECPIX in header\n");
	    return (-1);
	    }
	}
    else
	secpix = secpix0;
    *radperpixp = degrad (secpix / 3600.0);
    /* printf ("SECPIX is %g -> %g\n",secpix,*radperpixp); */

    /* Image size */
    nax = 0;
    if (hgeti4 (header,"NAXIS",&nax) < 1)
	return (-1);
    else {
	if (hgeti4 (header,"NAXIS1",wp) < 1)
	    return (-1);
	else {
	    if (hgeti4 (header,"NAXIS2",hp) < 1)
		return (-1);
	    }
	}

    if (verbose) {
	char rstr[64], dstr[64];
	ra2str (rstr, raddeg (*rap), 3);
	dec2str (dstr, raddeg (*decp), 2);
	printf ("RA=%s DEC=%s W=%d H=%d ArcSecs/Pixel=%g\n", rstr, dstr, 
				*wp, *hp, 3600.0*raddeg(*radperpixp));
	}

    return (0);
}


/* Project rap/decp onto plane at xp/yp using nominal ra0, dec0 and
 * build up a minimal FITS header so we can use RADectoXY().
 */

static void
sky2im (imw, imh, pixsz, ra0, dec0, n, rap, decp, xp, yp)

int	imw;
int	imh;
double	pixsz;
double	ra0;
double	dec0;
int	n;
double	rap[];
double	decp[];
double	xp[];
double	yp[];
{
    char *header;
    int i, headlen, off;
    struct WorldCoor *wcs, *wcsset();

/*  headlen = 2 * FITSBLOCK;
    header = malloc (headlen);
    strcpy (header, "END");
    for (i = 3; i < headlen; i++)
	header[i] = ' ';

    bare-bones header
    hputl (header, "SIMPLE", 1);
    hputi4 (header, "BITPIX", 16);
    hputi4 (header, "NAXIS", 2);
    hputi4 (header, "NAXIS1", imw);
    hputi4 (header, "NAXIS2", imh);

    add WCS header based on nominal position, no rotation
    x0 = (double) imw / 2.0;
    y0 = (double) imh / 2.0;
    setFITSWCS (header, rot0, ra0, dec0, x0, y0, pixsz);
 */
    wcs = wcsset (raddeg(ra0), raddeg(dec0), 3600.8*raddeg(pixsz),
		  imw, imh, rot0, 2000);

    /* use the nominal WCS info to find x/y on image */
    for (i = 0; i < n; i++) {
	xp[i] = 0.0;
	yp[i] = 0.0;
/*	RADec2xy (header, rap[i], decp[i], &xp[i], &yp[i]); */
	wcs2pix (wcs, raddeg(rap[i]), raddeg(decp[i]), &xp[i], &yp[i], &off);
	}

/*  free (header); */
    free (wcs);
    return;
}


/* Set FITS C* fields, assuming ra/dec refers to the center pixel */

static void
setFITSWCS (header, rot, ra, dec, x0, y0, pixsz)

char	*header;
double	rot;
double	ra;
double	dec;
double	x0;
double	y0;
double	pixsz;

{
    hputr8 (header, "EQUINOX", 2000.0);
    hputs  (header, "CTYPE1", "RA---TAN");
    hputr8 (header, "CRVAL1", raddeg(ra));
    hputr8 (header, "CDELT1", flip * raddeg(pixsz));
    hputr8 (header, "CRPIX1", x0);
    hputr8 (header, "CROTA1", raddeg(rot));

    hputs  (header, "CTYPE2", "DEC--TAN");
    hputr8 (header, "CRVAL2", raddeg(dec));
    hputr8 (header, "CDELT2", raddeg(pixsz));
    hputr8 (header, "CRPIX2", y0);
    hputr8 (header, "CROTA2", 0.0);
    return;
}

/* structure we use just for building the global star lists needed for sorting
 */
typedef struct {
    double n, ra, dec, x, y, b;
} _StarInfo;

/* Sort Guide Stars by increasing magnitude */

static void
sortGuideStars (n, ra, dec, sx, sy, sb, ns)

double	*n;
double	*ra;
double	*dec;
double	*sx;
double	*sy;
double	*sb;
int	ns;

{
    _StarInfo *sip;
    int i;

    sip = (_StarInfo *) malloc (ns * sizeof(_StarInfo));

    for (i = 0; i < ns; i++) {
	sip[i].n = n[i];
	sip[i].ra = ra[i];
	sip[i].dec = dec[i];
	sip[i].x = sx[i];
	sip[i].y = sy[i];
	sip[i].b = -sb[i];
    }

    qsort ((void *)sip, ns, sizeof(_StarInfo), starstatSortF);

    for (i = 0; i < ns; i++) {
	n[i] = sip[i].n;
	ra[i] = sip[i].ra;
	dec[i] = sip[i].dec;
	sx[i] = sip[i].x;
	sy[i] = sip[i].y;
	sb[i] = -sip[i].b;
    }

    free ((char *)sip);
    return;
}

/* Sort stars by decreasing brightness */

void
sortStars (sx, sy, sb, ns)

double	*sx;
double	*sy;
double	*sb;
int	ns;

{
    _StarInfo *sip;
    int i;

    sip = (_StarInfo *) malloc (ns * sizeof(_StarInfo));

    for (i = 0; i < ns; i++) {
	sip[i].x = sx[i];
	sip[i].y = sy[i];
	sip[i].b = sb[i];
    }

    qsort ((void *)sip, ns, sizeof(_StarInfo), starstatSortF);

    for (i = 0; i < ns; i++) {
	sx[i] = sip[i].x;
	sy[i] = sip[i].y;
	sb[i] = sip[i].b;
    }

    free ((char *)sip);
    return;
}

/* Sort in *decreasing* order */

static int
starstatSortF (ssp1, ssp2)

void *ssp1, *ssp2;

{
    double d = ((_StarInfo *)ssp2)->b - ((_StarInfo *)ssp1)->b;

    if (d < 0)
	return (-1);
    if (d > 0)
	return (1);
    return (0);
}
