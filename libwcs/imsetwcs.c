/* File libwcs/imsetwcs.c
 * September 4, 1996
 * By Doug Mink, based on UIowa code
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fitshead.h"
#include "wcs.h"
#include "lwcs.h"

#define GSC		1	/* refcat value for HST Guide Star Catalog */
#define UJC		2	/* refcat value for USNO UJ Star Catalog */

extern int FindStars ();
extern int TriMatch ();
extern int FocasMatch ();
extern int StarMatch ();
extern int findStars ();
extern int gscread();
extern int ujcread();
extern int tabread();
extern void MagSortStars ();
extern void FluxSortStars ();

extern void fk425(), fk425e(), fk524e();

static struct WorldCoor *GetFITSWCS ();
static void SetFITSWCS ();

/* set the C* WCS fields in  a FITS header based on a reference catalog
 * do it by finding stars in the image and in the reference catalog and
 * finding the rotation and offsets which result in a best-fit.
 * verbose generates extra info on stderr.
 * try using deeper reference star catalog searches if there is trouble.
 * return 1 if all ok, else 0
 */

/* These parameters can be set on the command line */
static double tolerance = PIXDIFF;	/* +/- this many pixels is a hit */
static double reflim = MAGLIM;		/* reference catalog magnitude limit */
static double secpix0 = PSCALE;		/* Set image scale--override header */
static int refcat = GSC;		/* reference catalog switch */
static char refcatname[32]="GSC";	/* reference catalog name */
static double rot0 = 0.0;		/* Initial image rotation */
static int fk4 = 0;			/* Command line center is FK4 */
static int classd = -1;			/* Guide Star Catalog object classes */
static int uplate = 0;			/* UJ Catalog plate number to use */
static double frac = 1.0;		/* Additional catalog/image stars */
static double ra0 = -99.0;		/* Initial center RA in degrees */
static double dec0 = -99.0;		/* Initial center Dec in degrees */
static char wcstype[8]="TAN";		/* WCS projection name */
static int nfit0 = 0;			/* Number of parameters to fit
					   (0=1 per match */


/* set the C* WCS fields in the input image header based on the given limiting
 * reference mag.
 * Finding stars in the input image and in the reference catalog down to
 * reflim and compute the angle and offsets which result in the best fit.
 * verbose generates extra info on stdout.
 * return 0 if all ok, else -1
 */

int
SetWCSFITS (header, image, verbose)

char	*header;	/* FITS header */
char	*image;		/* Image pixels */
int	verbose;

{
    double *gnum=0;	/* Reference star numbers */
    double *gra=0;	/* Reference star right ascensions in degrees */
    double *gdec=0;	/* Reference star declinations in degrees */
    double *gm=0;	/* Reference star magnitudes */
    double *gx=0;	/* Reference star image X-coordinates in pixels */
    double *gy=0;	/* Reference star image Y-coordinates in pixels */
    int *gc=0;		/* Reference object types */
    int *goff=0;	/* Reference star offscale flags */
    int ng;		/* Number of reference stars in image */
    int nbg;		/* Number of brightest reference stars actually used */
    double *sx=0;	/* Image star image X-coordinates in pixels */
    double *sy=0;	/* Image star image X-coordinates in pixels */
    double sra;		/* Image star right ascension in degrees */
    double sdec;	/* Image star declination in degrees */
    double *sb=0;	/* Image star integrated fluxes */
    int *sp=0;		/* Image star peak fluxes in counts */
    int ns;		/* Number of image stars */
    int nbs;		/* Number of brightest image stars actually used */
    double cra, cdec;	/* Nominal center in degrees from RA/DEC FITS fields */
    double dra, ddec;	/* Image half-widths in degrees */
    double secpix;	/* Pixel size in arcseconds */
    int imw, imh;	/* Image size, pixels */
    double ra1,ra2,dec1,dec2,mag1,mag2;
    int minstars;
    int ngmax;
    int nbin, nbytes;
    int ret = 0;
    int is, ig;
    char rstr[32], dstr[32];
    struct WorldCoor *wcs=0;	/* WCS structure */

    /* get nominal position and scale */
    wcs = GetFITSWCS (header,verbose,&cra,&cdec,&dra,&ddec,&secpix,&imw,&imh,2000);
    if (nowcs (wcs)) {
	ret = 0;
	goto out;
	}

    /* Set the RA and Dec limits in degrees for reference star search */
    ra1 = cra - dra;
    ra2 = cra + dra;
    dec1 = cdec - ddec;
    dec2 = cdec + ddec;
    if (reflim > 0.0) {
	mag1 = -2.0;
	mag2 = reflim;
	}
    else {
	mag1 = 0.0;
	mag2 = 0.0;
	}

    /* Allocate arrays for results of reference star search */
    ngmax = MAXREF;
    nbytes = MAXREF * sizeof (double);
    gnum = (double *) malloc (nbytes);
    gra = (double *) malloc (nbytes);
    gdec = (double *) malloc (nbytes);
    gm = (double *) malloc (nbytes);
    gc = (int *) malloc (nbytes);

    /* Find the nearby reference stars, in ra/dec */
    if (refcat == UJC)
	ng = ujcread (ra1,ra2,dec1,dec2,mag1,mag2,uplate,ngmax,
		      gnum,gra,gdec,gm,gc,verbose);
    else if (refcat == GSC)
	ng = gscread (ra1,ra2,dec1,dec2,mag1,mag2,classd,ngmax,
		      gnum,gra,gdec,gm,gc,verbose*100);
    else if (refcatname[0] > 0)
	ng = tabread (refcatname,ra1,ra2,dec1,dec2,mag1,mag2,ngmax,
		      gnum,gra,gdec,gm,gc,verbose);
    else {
	fprintf (stderr,"No reference star catalog specified\n");
	ret = 0;
	goto out;
	}

    if (nfit0 > 1)
	minstars = nfit0;
    else if (nfit0 == 1)
	minstars = 2;
    else
	minstars = MINSTARS;
	
    if (ng < minstars) {
	if (ng < 0)
	    fprintf (stderr, "Error getting reference stars: %d\n", ng);
	else if (ng == 0)
	    fprintf (stderr,"Need >= %d reference stars; found none\n",
							minstars);
	else
	    fprintf (stderr, "Need >= %d reference stars; found only %d\n",
							minstars, ng);
	ret = 0;
	goto out;
	}

    /* Project the reference stars into pixels on a plane at ra0/dec0 */
    gx = (double *) malloc (nbytes);
    gy = (double *) malloc (nbytes);
    goff = (int *) malloc (nbytes);
    if (!gx || !gy) {
	fprintf (stderr, "Could not malloc temp space of %d bytes\n",
					    ng*sizeof(double)*2);
	ret = 0;
	goto out;
	}

    /* use the nominal WCS info to find x/y on image */
    for (ig = 0; ig < ng; ig++) {
	gx[ig] = 0.0;
	gy[ig] = 0.0;
	wcs2pix (wcs, gra[ig], gdec[ig], &gx[ig], &gy[ig], &goff[ig]);
	}

    /* Sort reference stars by brightness (magnitude) */
    MagSortStars (gnum, gra, gdec, gx, gy, gm, gc, ng);

    /* Use only the brightest MAXSTARS reference stars */
    nbg = MAXSTARS;
    if (ng > nbg) {
	if (verbose)
	    printf ("using %d / %d reference stars brighter than %.1f\n",
		     nbg, ng, gm[nbg-1]);
	}
    else {
	nbg = ng;
	if (verbose) {
	    if (reflim > 0.0)
		printf ("using all %d reference stars brighter than %.1f\n",
			ng,reflim);
	    else
		printf ("using all %d reference stars\n", ng);
	    }
	}

    printf ("%s:\n",refcatname);
    for (ig = 0; ig < ng; ig++) {
	ra2str (rstr, gra[ig], 3);
	dec2str (dstr, gdec[ig], 2);
	printf (" %9.4f %s %s %5.2f %6.1f %6.1f\n",
		gnum[ig],rstr,dstr,gm[ig],gx[ig],gy[ig]);
	}

    /* Discover star-like things in the image, in pixels */
    ns = FindStars (header, image, &sx, &sy, &sb, &sp, verbose);
    if (ns < minstars) {
	fprintf (stderr, "Need at least %d image stars but only found %d\n",
							minstars, ns);
	ret = 0;
	goto out;
	}

    /* Sort star-like objects in image by brightness */
    FluxSortStars (sx, sy, sb, sp, ns);

    /* Use only as many star-like objects as reference stars */
    /* (Actually use frac * the number of GS if more than image stars or */
    /*  frac * the number of image stars if more than reference stars) */
    if (ns > nbg) {
	nbs = nbg * frac;
	if (nbs > ns)
	   nbs = ns;
	if (verbose) {
	    printf ("using brightest %d / %d reference stars\n", nbg, ng);
	    printf ("using brightest %d / %d image stars\n", nbs,ns);
	    }
	}
    else {
	nbs = ns;
	nbg = nbs * frac;
	if (nbg > ng)
	    nbg = ng;
	if (verbose)
	    printf ("using all %d image stars\n", ns);
	}

    if (verbose) {
	for (is = 0; is < nbs; is++)
	printf ("Im: %6.1f %6.1f %8.1f %d\n", sx[is],sy[is],sb[is],sp[is]);
	}

    /* Check offsets between all pairs of image stars and reference stars */
    nbin = StarMatch (nbs,sx,sy, nbg,gra,gdec,gx,gy, tolerance, wcs, nfit0, verbose);

    if (nbin < 0) {
	fprintf (stderr, "Star registration failed.\n");
	ret = 0;
	goto out;
	}
    else if (verbose)
	printf ("%d bin hits\n", nbin);

    SetFITSWCS (header, wcs);

    if (verbose) {
	double ra,dec;
	printf ("Rotation:  %g degrees  Arcsec/pixel: %f %f\n",
		wcs->rot, 3600.0*wcs->xinc, 3600.0*wcs->yinc);
	ra2str (rstr, wcs->xref, 3);
	dec2str (dstr, wcs->yref, 2);
	printf ("New center: %s  %s (FK5)\n",rstr,dstr);
	ra = wcs->xref;
	dec = wcs->yref;
	(void)fk524e (&ra,&dec,wcs->epoch);
	ra2str (rstr, ra, 3);
	dec2str (dstr, dec, 2);
	printf ("New center: %s  %s (FK4)\n",rstr,dstr);
	}

    /* Find star matches for this offset and print them */
    if (verbose) {
	double x, y, dx, dy, dxy, sep, xysep, tol2;
	double dxsum=0.0, dysum=0.0, dxysum=0.0, sepsum=0.0;
	int offscl, nmatch=0;
	/* wcs = wcsinit (header); */
	tol2 = tolerance * tolerance;
	printf ("Matches in %s:\n", refcatname);
	for (ig = 0; ig < ng; ig++) {
	    wcs2pix (wcs, gra[ig], gdec[ig], &x, &y, &offscl);
	    if (!offscl) {
		for (is = 0; is < ns; is++) {
		    dx = x - sx[is];
		    dy = y - sy[is];
		    dxy = (dx * dx) + (dy * dy);
		    if (dxy < tol2) {
			nmatch++;
			dxsum = dxsum + dx;
			dysum = dysum + dy;
			dxysum = dxysum + sqrt (dxy);
			pix2wcs (wcs, sx[is], sy[is], &sra, &sdec);
			xysep = sqrt (dxy);
			sep = 3600.0 * wcsdist (gra[ig], gdec[ig], sra, sdec);
			sepsum = sepsum + sep;
			ra2str (rstr, gra[ig], 3);
			dec2str (dstr, gdec[ig], 2);
			printf ("%9.4f %s %s %5.2f",
			    gnum[ig], rstr, dstr, gm[ig]);
			printf ("%6.1f %6.1f %6d %5.2f %5.2f\n",
			    sx[is], sy[is], sp[is], xysep, sep);
			}
		    }
		}
	    }
	if (nmatch > 0) {
	    dx = dxsum / (double)nmatch;
	    dy = dysum / (double)nmatch;
	    dxy = dxysum / (double)nmatch;
	    sep = sepsum / (double)nmatch;
	    printf ("Mean dx= %.2f  dy= %.2f  dxy= %.2f  sep= %.2f\n",
		    dx, dy, dxy, sep);
	    }
	else
	    fprintf (stderr,"SetWCSFITS: Error in WCS, no matches found\n");
	}

    ret = 1;

    out:
    if (gra) free ((char *)gra);
    if (gdec) free ((char *)gdec);
    if (gm) free ((char *)gm);
    if (gnum) free ((char *)gnum);
    if (gx) free ((char *)gx);
    if (gy) free ((char *)gy);
    if (gc) free ((char *)gc);

    if (sx) free ((char *)sx);
    if (sy) free ((char *)sy);
    if (sb) free ((char *)sb);
    if (sp) free ((char *)sp);
    if (wcs) free (wcs);

    return (ret);
}


/* Set a nominal world coordinate system from image header info.
 * If the image center is not FK5 (J2000) equinox, convert it
 * Return a WCS structure if OK, else return NULL
 */

static struct WorldCoor *
GetFITSWCS (header, verbose, cra, cdec, dra, ddec, secpix, wp, hp, eqref)

char	*header;	/* Image FITS header */
int	verbose;	/* Extra printing if =1 */
double	*cra;		/* Center right ascension in degrees (returned) */
double	*cdec;		/* Center declination in degrees (returned) */
double	*dra;		/* Right ascension half-width in degrees (returned) */
double	*ddec;		/* Declination half-width in degrees (returned) */
double	*secpix;	/* Arcseconds per pixel (returned) */
int	*wp;		/* Image width in pixels (returned) */
int	*hp;		/* Image height in pixels (returned) */
int	eqref;		/* Equinox of reference catalog */
{
    int nax;
    int equinox, eqcoor;
    double epoch;
    struct WorldCoor *wcs;

    /* Set image dimensions */
    nax = 0;
    if (hgeti4 (header,"NAXIS",&nax) < 1)
	return (NULL);
    else {
	if (hgeti4 (header,"NAXIS1",wp) < 1)
	    return (NULL);
	else {
	    if (hgeti4 (header,"NAXIS2",hp) < 1)
		return (NULL);
	    }
	}

    /* Find center for pre-existing WCS, if there is one */
    wcs = wcsinit (header);
    if (iswcs (wcs)) {
	wcssize (wcs, cra, cdec, dra, ddec);
	if (wcs->xref == 0.0 && wcs->yref == 0.0) {
	    wcs->xref = *cra;
	    wcs->yref = *cdec;
	    wcs->xrefpix = (double) wcs->nxpix * 0.5;
	    wcs->yrefpix = (double) wcs->nypix * 0.5;
	    wcs->xinc = *dra / wcs->xrefpix;
	    wcs->yinc = *ddec / wcs->yrefpix;
	    hchange (header,"PLTRAH","PLT0RAH");
	    wcs->plate_fit = 0;
	    }
	if (strcmp (wcs->radecsys, "FK4") == 0) {
	    fk425e (cra, cdec, wcs->epoch);
	    wcsshift (wcs, *cra, *cdec, "FK5");
	    }
	}

    /* Otherwise use nominal center from RA and DEC fields */
    else {
	*cra = 0.0;
	*cdec = 0.0;
	if (hgetra (header, "RA", cra) == 0) {
	    fprintf (stderr, "No RA field\n");
	    return (NULL);
	    }
	if (hgetdec (header, "DEC", cdec) == 0) {
	    fprintf (stderr, "No DEC field\n");
	    return (NULL);
	    }

	/* Equinox of coordinates */
	if (hgeti4 (header, "EPOCH", &equinox) == 0) {
	    if (hgeti4 (header, "EQUINOX", &equinox) == 0)
		equinox = 1950;
	    }

	/* Epoch of image (observation date) */
	if (hgetdate (header," OBS-DATE", &epoch) == 0) {
	    if (hgetdate (header," DATE-OBS", &epoch) == 0)
		epoch = (double) equinox;
	    }

	/* If coordinate equinox not reference catalog equinox, convert */
	if (equinox != eqref) {
	    if (eqref == 2000)
		fk425e (cra, cdec, epoch);
	    else
		fk524e (cra, cdec, epoch);
	    }
	}

    /* Set plate scale from command line, if it is there */
    if (secpix0 > 0.0) {
	*dra = (secpix0 * *wp * 0.5 / 3600.0) / cos (degrad(*cdec));
	*ddec = secpix0 * *hp * 0.5 / 3600.0;
	*secpix = secpix0;
	hputnr8 (header,"SECPIX",5,*secpix);
	wcs->yinc = secpix0 / 3600.0;
	wcs->xinc = -wcs->yinc;
	}

    /* Otherwise set plate scale from FITS header */
    else {

	/* Plate scale from WCS if it is present */
	if (iswcs (wcs))
	    *secpix = 3600.0 * 2.0 * *ddec / (double) *hp;

	/* Plate scale from SECPIX header parameter */
	else {
	    *secpix = 0.0;
	    if (hgetr8 (header, "SECPIX", secpix) == 0) {
		if (hgetr8 (header, "SECPIX1", secpix) == 0) {
		    fprintf (stderr, "Cannot find SECPIX in header\n");
		    return (NULL);
		    }
		}
	    }
	}

    /* Set WCS structure if it has not already been set */
    if (nowcs (wcs))
	wcs = wcsset (*cra, *cdec, *secpix, *wp, *hp, rot0, eqref, epoch);

    /* Reset to center position from the command line, if there is one */
    if (ra0 > -99.0 && dec0 > -99.0) {
	char rstr[32],dstr[32];

    /* Reset center for reference star search */
	*cra = ra0;
	*cdec = dec0;

	/* If coordinate equinox not reference catalog equinox, convert */
	if (fk4)
	    eqcoor = 1950;
	else
	    eqcoor = 2000;
	if (eqcoor != eqref) {
	    if (eqref == 2000)
		fk425e (cra, cdec, epoch);
	    else
		fk524e (cra, cdec, epoch);
	    }
	if (fk4)
	    wcsshift (wcs, *cra, *cdec, "FK4");
	else
	    wcsshift (wcs, *cra, *cdec, "FK5");

	ra2str (rstr, ra0, 3);
        dec2str (dstr, dec0, 2);
	hputs (header,"RA",rstr);
	hputs (header,"DEC",dstr);
	hputi4 (header,"EQUINOX",eqcoor);

	if (verbose) {
	    if (eqcoor != eqref) {
		if (fk4)
		    printf ("Center reset to RA=%s DEC=%s (FK4)\n", rstr, dstr);
		else
		    printf ("Center reset to RA=%s DEC=%s (FK5)\n", rstr, dstr);
		ra2str (rstr, *cra, 3);
        	dec2str (dstr, *cdec, 2);
		}
	    if (eqref == 2000)
		printf ("Center reset to RA=%s DEC=%s (FK5)\n", rstr, dstr);
	    else
		printf ("Center reset to RA=%s DEC=%s (FK4)\n", rstr, dstr);
	    }
	}

    /* Image size from header */

    if (verbose) {
	char rstr[64], dstr[64];
	ra2str (rstr, *cra, 3);
	dec2str (dstr, *cdec, 2);
	printf ("RA=%s DEC=%s W=%d H=%d ArcSecs/Pixel=%g\n", rstr, dstr, 
				*wp, *hp, *secpix);
	}

    return (wcs);
}


/* Set FITS C* fields, assuming ra/dec refers to the center pixel */

static void
SetFITSWCS (header, wcs)

char	*header;	/* Image FITS header */
struct WorldCoor *wcs;	/* WCS structure */

{
    double ep, ra, dec;
    char wcstemp[16];

    /* Rename old center coordinates */
    if (hgetra (header,"RA", &ra))
	hchange (header,"RA","WRA");
    if (hgetdec (header,"DEC", &dec))
	hchange (header,"DEC","WDEC");

    if (hgetr8 (header, "EQUINOX", &ep))
	hchange (header, "EQUINOX", "WEQUINOX");

    /* Only change EPOCH if it is used instead of EQUINOX */
    else if (hgetr8 (header, "EPOCH", &ep))
	hchange (header, "EPOCH", "WEPOCH");

    /* Set new center coordinates */
    hputra (header,"RA",wcs->xref);
    hputdec (header,"DEC",wcs->yref);
    hputr8 (header, "EQUINOX", wcs->equinox);
    hputs (header, "RADECSYS", wcs->radecsys);

    strcpy (wcstemp, "RA---");
    strcat (wcstemp, wcstype);
    hputs  (header, "CTYPE1", wcstemp);
    hputnr8 (header, "CRVAL1", 9, wcs->xref);
    hputnr8 (header, "CRPIX1", 4, wcs->xrefpix);
    hputnr8 (header, "CDELT1", 9, wcs->xinc);
    hputnr8 (header, "CROTA1", 3, wcs->rot);

    strcpy (wcstemp, "DEC--");
    strcat (wcstemp, wcstype);
    hputs  (header, "CTYPE2", wcstemp);
    hputnr8 (header, "CRVAL2", 9, wcs->yref);
    hputnr8 (header, "CRPIX2", 4, wcs->yrefpix);
    hputnr8 (header, "CDELT2", 9, wcs->yinc);
    hputnr8 (header, "CROTA2", 3, 0.0);
    return;
}

void settolerance (tol)
double tol;
{ tolerance = tol; return; }

void setrefcat (cat)
char *cat;
{  if (cat[0]=='G' || cat[0]=='g' || strcmp(cat,"gsc")==0 || strcmp(cat,"GSC")==0)
	refcat = GSC;
   else if (cat[0]=='U' || cat[0]=='u' || strcmp(cat,"ujc")==0 || strcmp(cat,"UJC")==0)
	refcat = UJC;
    else
	refcat = 0;
    strcpy (refcatname, cat); return; }

void setwcstype (type)
char *type;
{ strcpy (wcstype, type); return; }


void setreflim (lim)
double lim;
{ reflim = lim; return; }

void setrot (rot)
double rot;
{ rot0 = rot; return; }

void setsecpix (secpix)
double secpix;
{ secpix0 = secpix; return; }

void setfk4 ()
{ fk4 = 1; return; }

void setclass (class)
int class;
{ classd = class; return; }

void setplate (plate)
int plate;
{ uplate = plate; return; }

void setfrac (frac0)
double frac0;
{ if (frac0 < 1.0) frac = 1.0 + frac0;
    else frac = frac0; return; }

void setcenter (rastr, decstr)
char *rastr, *decstr;
{ ra0 = str2ra (rastr); dec0 = str2dec (decstr); return; }

void setnfit (nfit)
int nfit;
{ if (nfit < 2) nfit = 2;
  else if (nfit > 5) nfit = 5;
  nfit0 = nfit; return; }


/* Feb 29 1996	New program
 * Apr 30 1996	Add FOCAS-style catalog matching
 * May  1 1996	Add initial image center from command line
 * May  2 1996	Set up four star matching modes
 * May 15 1996	Pass verbose flag; allow other reference catalogs
 * May 16 1996	Remove sorting to separate file sortstar.c
 * May 17 1996	Add class and verbose arguments
 * May 22 1996  Allow for named reference catalog
 * May 23 1996	Use pre-existing WCS for center, if it is present
 * May 29 1996	Simplify program by always using WCS structure
 * May 30 1996	Make reference/image pair matching the default method
 * Jun 11 1996  Number and zero positions of image stars
 * Jun 12 1996	Be more careful with nominal WCS setting
 * Jun 14 1996	Add residual table
 * Jun 28 1996	Set FITS header from WCS
 * Jul  3 1996	Set epoch from old equinox if not already set
 * Jul 19 1996	Declare tabread
 * Jul 19 1996	Set image center in WCS if DSS WCS
 * Jul 22 1996	Debug tab table reading
 * Aug  5 1996	Add option to change WCS projection
 * Aug  5 1996	Check for SECPIX1 as well as SECPIX
 * Aug  5 1996	Set number of parameters to fit here
 * Aug  7 1996	Save specified number of decimal places in header parameters
 * Aug  7 1996	Rename old center parameters
 * Aug 26 1996	Decrease default pixel tolerance from 20 to 10
 * Aug 26 1996	Drop unused variable EQ in SETFITSWCS
 * Aug 28 1996	Improve output format for matched stars
 * Sep  1 1996	Set some defaults in lwcs.h
 * Sep  3 1996	Fix bug to set plate scale from command line
 * Sep  4 1996	Print reference catalog name on separate line from entries
 * Sep 17 1996	Clean up code
 */
