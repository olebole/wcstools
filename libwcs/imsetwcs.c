/* File libwcs/imsetwcs.c
 * December 1, 1998
 * By Doug Mink, based on UIowa code
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fitshead.h"
#include "wcs.h"
#include "lwcs.h"
#include "wcscat.h"

extern int FindStars();
extern int TriMatch();
extern int FocasMatch();
extern int StarMatch();
extern int FitPlate();
extern struct WorldCoor *GetFITSWCS ();
extern char *getimcat();
extern void SetFITSWCS();
extern void fk524e();


/* Set the C* WCS fields in a FITS header based on a reference catalog
 * by finding stars in the image and in the reference catalog and
 * fitting the scaling, rotation, and offsets.
 * verbose generates extra info on stderr.
 * Try using deeper reference star catalog searches if there is trouble.
 * Return 1 if all ok, else 0
 */

/* These parameters can be set on the command line */
static double tolerance = PIXDIFF;	/* +/- this many pixels is a hit */
static double refmag1 = MAGLIM1;	/* reference catalog magnitude limit */
static double refmag2 = MAGLIM2;	/* reference catalog magnitude limit */
static int refcat = GSC;		/* reference catalog switch */
static char refcatname[32]="GSC";	/* reference catalog name */
static int classd = -1;			/* Guide Star Catalog object classes */
static int uplate = 0;			/* UJ Catalog plate number to use */
static double frac = 1.0;		/* Additional catalog/image stars */
static int nofit = 0;			/* if =1, do not fit WCS */
static int maxcat = MAXSTARS;		/* Maximum number of catalog stars to use */
static int fitwcs = 1;			/* If 1, fit WCS, else use current WCS */
static int fitplate = 0;		/* If 1, fit polynomial, else do not */
static double imfrac0 = 0.0;		/* If > 0.0, multiply image dimensions
					   by this for search */
static int iterate0 = 0;		/* If 1, search field again */
static int recenter0 = 0;		/* If 1, search again with new center*/
static char matchcat[32]="";		/* Match catalog name */
static void PrintRes();
extern void SetFITSPlate();

/* Set the C* WCS fields in the input image header based on the given limiting
 * reference mag.
 * Finding stars in the input image and in the reference catalog between
 * refmag1 and refmag2 and compute the angle and offsets which result in the best fit.
 * verbose generates extra info on stdout.
 * return 0 if all ok, else -1
 */

int
SetWCSFITS (filename, header, image, refcatname, verbose)

char	*filename;	/* image file name */
char	*header;	/* FITS header */
char	*image;		/* Image pixels */
char	*refcatname;	/* Name of reference catalog */
int	verbose;

{
    double *gnum=0;	/* Reference star numbers */
    double *gra=0;	/* Reference star right ascensions in degrees */
    double *gdec=0;	/* Reference star declinations in degrees */
    double *gm=0;	/* Reference star magnitudes */
    double *gmb=0;	/* Reference star blue magnitudes */
    double *gx=0;	/* Reference star image X-coordinates in pixels */
    double *gy=0;	/* Reference star image Y-coordinates in pixels */
    int *gc=0;		/* Reference object types */
    int *goff=0;	/* Reference star offscale flags */
    int ng;		/* Number of reference stars in image */
    int nbg;		/* Number of brightest reference stars actually used */
    double *sx=0;	/* Image star image X-coordinates in pixels */
    double *sy=0;	/* Image star image X-coordinates in pixels */
    double *sb=0;	/* Image star image flux in counts */
    int *sp=0;		/* Image star peak fluxes in counts */
    int ns;		/* Number of image stars */
    int nbs;		/* Number of brightest image stars actually used */
    double cra, cdec;	/* Nominal center in degrees from RA/DEC FITS fields */
    double dra, ddec;	/* Image half-widths in degrees */
    double secpix;	/* Pixel size in arcseconds */
    int imw, imh;	/* Image size, pixels */
    int imsearch = 1;	/* Flag set if image should be searched for sources */
    double mag1,mag2;
    double dxys;
    int minstars;
    int ngmax;
    int nbin, nbytes;
    int iterate = iterate0;
    int recenter = recenter0;
    int ret = 0;
    int is, ig, igs;
    char rstr[32], dstr[32];
    double refeq, refep;
    int refsys;
    char refcoor[8];
    char title[80];
    char *imcatname;	/* file name for image star catalog, if used */
    struct WorldCoor *wcs=0;	/* WCS structure */
    double *sx1, *sy1, *gra1, *gdec1, *gnum1, *gm1;
    double imfrac = imfrac0;

    /* Set reference catalog coordinate system and epoch */
    refcat = RefCat (refcatname, title, &refsys, &refeq, &refep);
    wcscstr (refcoor, refsys, refeq, refep);

    /* get nominal position and scale */
getfield:
    wcs = GetFITSWCS (header,verbose,&cra,&cdec,&dra,&ddec,&secpix,&imw,&imh,
		      &refsys, &refeq);
    if (nowcs (wcs)) {
	ret = 0;
	goto out;
	}

    if (nofit) {
	SetFITSWCS (header, wcs);
	ret = 1;
	goto out;
	}
    wcs->prjcode = WCS_TAN;

    wcseqset (wcs, refeq);

    if (imfrac > 0.0) {
	dra = dra * imfrac;
	ddec = ddec * imfrac;
	}

    if (refmag1 == refmag2) {
	mag1 = 0.0;
	mag2 = 0.0;
	}
    else {
	mag1 = refmag1;
	mag2 = refmag2;
	}

    /* Allocate arrays for results of reference star search */
    ngmax = maxcat;
    if (imfrac > 1.0)
	ngmax = (int) ((double) ngmax * imfrac * imfrac);
    nbytes = ngmax * sizeof (double);
    gnum = (double *) malloc (nbytes);
    gra = (double *) malloc (nbytes);
    gdec = (double *) malloc (nbytes);
    gm = (double *) malloc (nbytes);
    gmb = (double *) malloc (nbytes);
    gc = (int *) malloc (nbytes);

    /* Find the nearby reference stars, in ra/dec */
    if (refcat == UAC || refcat == UA1 || refcat == UA2 ||
	refcat == USAC || refcat == USA1 || refcat == USA2)
	ng = uacread (refcatname,cra,cdec,dra,ddec,0.0,refsys,refeq,wcs->epoch,
		      mag1,mag2,uplate,ngmax,gnum,gra,gdec,gm,gmb,gc,verbose*1000);
    else if (refcat == UJC)
	ng = ujcread (cra,cdec,dra,ddec,0.0,refsys,refeq,wcs->epoch,
		      mag1,mag2,uplate,ngmax,gnum,gra,gdec,gm,gc,verbose);
    else if (refcat == GSC)
	ng = gscread (cra,cdec,dra,ddec,0.0,refsys,refeq,wcs->epoch,mag1,mag2,
		      classd,ngmax,gnum,gra,gdec,gm,gc,verbose*100);
    else if (refcat == SAO)
	ng = binread ("SAOra",cra,cdec,dra,ddec,0.0,refsys,refeq,wcs->epoch,
		      mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,NULL,verbose*100);
    else if (refcat == PPM)
	ng = binread ("PPMra",cra,cdec,dra,ddec,0.0,refsys,refeq,wcs->epoch,
		      mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,NULL,verbose*100);
    else if (refcat == IRAS)
	ng = binread ("IRAS",cra,cdec,dra,ddec,0.0,refsys,refeq,wcs->epoch,
		      mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,NULL,verbose*100);
    else if (refcat == TYCHO)
	ng = binread ("tychora",cra,cdec,dra,ddec,0.0,refsys,refeq,wcs->epoch,
		      mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,NULL,verbose*100);
    else if (refcat == BINCAT)
	ng = binread (refcatname,cra,cdec,dra,ddec,0.0,refsys,refeq,wcs->epoch,
		      mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,NULL,verbose*100);
    else if (refcat == TABCAT)
	ng = tabread (refcatname,cra,cdec,dra,ddec,0.0,refsys,refeq,wcs->epoch,
		      mag1,mag2,ngmax,gnum,gra,gdec,gm,gc,verbose);
    else if (refcat == TXTCAT)
	ng = catread (refcatname,cra,cdec,dra,ddec,0.0,refsys,refeq,wcs->epoch,
		      mag1,mag2,ngmax,gnum,gra,gdec,gm,verbose);
    else {
	fprintf (stderr,"No reference star catalog specified\n");
	ret = 0;
	goto out;
	}

    minstars = MINSTARS;

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
    MagSortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, NULL, ng);

    /* Note how reference stars were selected */
    nbg = ngmax;
    if (ng > nbg) {
	if (verbose)
	    printf ("using %d / %d reference stars brighter than %.1f\n",
		     nbg, ng, gm[nbg-1]);
	}
    else {
	nbg = ng;
	if (verbose) {
	    if (refmag1 > 0.0 && refmag2 > 0.0)
		printf ("using all %d reference stars from %.1f to %.1f mag.\n",
			ng, refmag1, refmag2);
	    else if (refmag2 > 0.0)
		printf ("using all %d reference stars brighter than %.1f\n",
			ng,refmag2);
	    else
		printf ("using all %d reference stars\n", ng);
	    }
	}

    if (verbose) {
	printf ("%s:\n",refcatname);
	for (ig = 0; ig < ng; ig++) {
	    ra2str (rstr, 32, gra[ig], 3);
	    dec2str (dstr, 32, gdec[ig], 2);
	    if (refcat == UAC || refcat == UA1 || refcat == UA2 ||
		refcat == USAC || refcat == USA1 || refcat == USA2)
		printf (" %13.8f",gnum[ig]);
	    else if (refcat == UJC)
		printf (" %12.7f",gnum[ig]);
	    else if (refcat == GSC)
		printf (" %9.4f",gnum[ig]);
	    else
		printf (" %9.4f",gnum[ig]);
	    printf (" %s %s %5.2f %6.1f %6.1f\n",
		    rstr,dstr,gm[ig],gx[ig],gy[ig]);
	    }
	}

    if (ng < minstars) {
	if (ng < 0)
	    fprintf (stderr, "Error getting reference stars: %d\n", ng);
	else if (ng == 0)
	    fprintf (stderr,"No reference stars found in image area\n");
	else if (fitwcs)
	    fprintf (stderr, "Found only %d out of %d reference stars needed\n",
			     ng, minstars);
	if (ng <= 0 || fitwcs) {
	    ret = 0;
	    goto out;
	    }
	}

    /* Discover star-like things in the image, in pixels */
    if (imsearch) {
	ns = FindStars (header, image, &sx, &sy, &sb, &sp, verbose);
	if (ns < minstars) {
	    if (ns < 0)
		fprintf (stderr, "Error getting image stars: %d\n", ns);
	    else if (ns == 0)
		fprintf (stderr,"No stars found in image\n");
	    else if (fitwcs)
		fprintf (stderr, "Need at least %d image stars but only found %d\n",
			 minstars, ns);
	    if (ns <= 0 || fitwcs) {
		ret = 0;
		iterate = 0;
		recenter = 0;
		goto out;
		}
	    }
	imsearch = 0;
	}

    /* Fit a world coordinate system if requested */
    if (fitwcs) {

	/* Sort star-like objects in image by brightness */
	FluxSortStars (sx, sy, sb, sp, ns);

	/* If matching a catalog field the same size as the image field,
	   use only as many star-like objects as reference stars.  If using
	   a larger catalog field (imfrac > 0), increase the number of stars
	   proportionately (imfrac^2 * the number of image stars).  To adjust
	   for different bandpasses between the image and the reference catalog
	   use frac * the number of reference stars if more than image stars
	   or frac * the number of image stars if more than reference stars.
	*/
	if (imfrac > 0.0) {
	    double imfrac2 = imfrac * imfrac;
	    if (ns > (nbg / imfrac2)) {
		nbs = nbg * frac / imfrac2;
		if (nbs > ns)
		    nbs = ns;
		}
	    else {
		nbs = ns;
		nbg = (int) ((double) nbs * imfrac2);
		if (nbg > ng)
		    nbg = ng;
		}
	    }
	else {
	    if (ns > nbg) {
		nbs = nbg * frac;
		if (nbs > ns)
		    nbs = ns;
		}
	    else {
		nbs = ns;
		nbg = nbs * frac;
		if (nbg > ng)
		    nbg = ng;
		}
	    }
	if (verbose) {
	    if (nbg == ng)
		printf ("Using all %d reference stars\n", ng);
	    else
		printf ("Using brightest %d / %d reference stars\n", nbg, ng);
	    if (nbs == ns)
		printf ("Using all %d image stars\n", ns);
	    else
		printf ("Using brightest %d / %d image stars\n", nbs,ns);
	    }

	if (verbose) {
	    char rastr[32], decstr[32];
	    double xmag, mdiff, ra, dec;
	    for (is = 0; is < nbs; is++) {
		pix2wcs (wcs, sx[is], sy[is], &ra, &dec);
		ra2str (rastr, 32, ra, 3);
		dec2str (decstr, 32, dec, 2);
		xmag = -2.5 * log10 (sb[is]);
		if (!is) mdiff = gm[0] - xmag;
		xmag = xmag + mdiff;
		printf ("%4d %s %s %6.2f %6.1f %6.1f %d\n",
			is+1, rastr, decstr, xmag, sx[is], sy[is], sp[is]);
		}
	    }

	/* Match offsets between all pairs of image stars and reference stars
	   and fit WCS to matches */
	nbin = StarMatch (nbs,sx,sy, nbg,gra,gdec,gx,gy, tolerance, wcs,
			  verbose);

	if (nbin < 0) {
	    fprintf (stderr, "Star registration failed.\n");
	    ret = 0;
	    goto out;
	    }
	else if (verbose)
	    printf ("%d bin hits\n", nbin);

	SetFITSWCS (header, wcs);
	}

    /* Match reference and image stars */
    if (verbose || !fitwcs) {
	double ra,dec;
	double dx, dy, dx2, dy2, dxy, tol2;
	int nmatch=0;

	if (wcs->ncoeff1 > 0)
	    printf ("# %d-term x, %d-term y polynomial fit\n",
		    wcs->ncoeff1, wcs->ncoeff2);
	else
	    printf ("# Arcsec/Pixel: %.6f %.6f  Rotation: %.6f degrees\n",
		    3600.0*wcs->xinc, 3600.0*wcs->yinc, wcs->rot);
	ra2str (rstr, 32, wcs->xref, 3);
	dec2str (dstr, 32, wcs->yref, 2);
	printf ("# Optical axis: %s  %s %s at (%.2f,%.2f)\n",
		rstr,dstr, refcoor, wcs->xrefpix, wcs->yrefpix);
	ra = wcs->xref;
	dec = wcs->yref;
	if (refsys == WCS_J2000) {
	    (void)fk524e (&ra, &dec, wcs->epoch);
	    ra2str (rstr, 32, ra, 3);
	    dec2str (dstr, 32, dec, 2);
	    printf ("# Optical axis: %s  %s B1950 at (%.2f,%.2f)\n",
		    rstr,dstr, wcs->xrefpix, wcs->yrefpix);
	    }
	else {
	    (void)fk425e (&ra, &dec, wcs->epoch);
	    ra2str (rstr, 32, ra, 3);
	    dec2str (dstr, 32, dec, 2);
	    printf ("# Optical axis: %s  %s J2000 at (%.2f,%.2f)\n",
		    rstr,dstr, wcs->xrefpix, wcs->yrefpix);
	    }

	/* Find star matches for this offset and print them */
	tol2 = tolerance * tolerance;

	/* Use the fit WCS info to find catalog star x/y on image */
	for (ig = 0; ig < ng; ig++) {
	    gx[ig] = 0.0;
	    gy[ig] = 0.0;
	    wcs2pix (wcs, gra[ig], gdec[ig], &gx[ig], &gy[ig], &goff[ig]);
	    }

	/* Find best catalog matches to stars in image */
	nmatch = 0;
	nbytes = ns * sizeof (double);
	gra1 = (double *) malloc (nbytes);
	gdec1 = (double *) malloc (nbytes);
	gm1 = (double *) malloc (nbytes);
	gnum1 = (double *) malloc (nbytes);
	sx1 = (double *) malloc (nbytes);
	sy1 = (double *) malloc (nbytes);
	for (is = 0; is < ns; is++) {
	    dxys = -1.0;
	    igs = -1;
	    for (ig = 0; ig < ng; ig++) {
		if (!goff[ig]) {
		    dx = gx[ig] - sx[is];
		    dy = gy[ig] - sy[is];
		    dx2 = dx * dx;
		    dy2 = dy * dy;
		    dxy = dx2 + dy2;
		    if (dxy < tol2) {
			if (dxys < 0.0) {
			    dxys = dxy;
			    igs = ig;
			    }
			else if (dxy < dxys) {
			    dxys = dxy;
			    igs = ig;
			    }
		        }
		    }
		}
	    if (igs > -1) {
		gnum1[nmatch] = gnum[igs];
		gm1[nmatch] = gm[igs];
		gra1[nmatch] = gra[igs];
		gdec1[nmatch] = gdec[igs];
		sx1[nmatch] = sx[is];
		sy1[nmatch] = sy[is];
		nmatch++;
		}
	    }

	imcatname = getimcat ();
	if (strlen (imcatname) == 0)
	    imcatname = filename;

	/* If there were any matches found, print them */
	if (nmatch > 0) {
	    printf ("# %d matches between %s and %s:\n",
		    nmatch, refcatname, imcatname);

	    PrintRes (wcs, nmatch, sx1, sy1, gra1, gdec1, gm1, gnum1, refcat);

	    /* Fit the matched catalog and image stars with a polynomial */
	    if (!iterate && !recenter && fitplate) {

		if (verbose)
		    printf ("Fitting matched stars with a polynomial\n");

		/* Fit residuals */
		if (FitPlate (wcs, sx1, sy1, gra1, gdec1, nmatch, fitplate,
			      verbose))
		    printf ("FitPlate cannot fit matches\n");

		/* Print the new residuals */
		else {
		    printf ("# %d matches between %s and %s:\n",
			    nmatch, refcatname, imcatname);
		    PrintRes (wcs,nmatch,sx1,sy1,gra1,gdec1,gm1,gnum1,refcat);
		    SetFITSPlate (header, wcs);
		    }
		}
	    }
	else
	    fprintf (stderr, "SetWCSFITS: No matches between %s and %s:\n",
		    refcatname, imcatname);
	if (gra1) free ((char *)gra1);
	if (gdec1) free ((char *)gdec1);
	if (gm1) free ((char *)gm1);
	if (gnum1) free ((char *)gnum1);
	}

    ret = 1;

    out:

    /* Free catalog source arrays */
    if (gra) free ((char *)gra);
    if (gdec) free ((char *)gdec);
    if (gm) free ((char *)gm);
    if (gmb) free ((char *)gmb);
    if (gnum) free ((char *)gnum);
    if (gx) free ((char *)gx);
    if (gy) free ((char *)gy);
    if (gc) free ((char *)gc);

    if (wcs) free (wcs);
    if (iterate) {
	setdcenter (wcs->xref, wcs->yref);
	setrefpix (wcs->xrefpix, wcs->yrefpix);
	setsecpix (3600.0 * wcs->xinc);
	setsecpix2 (3600.0 * wcs->yinc);
	setrot (wcs->rot);
	iterate = 0;
	imfrac = 0.0;
	goto getfield;
	}
    if (recenter) {
	double ra, dec, x, y;
	x = 0.5*wcs->nxpix;
	y = 0.5*wcs->nypix;
	pix2wcs (wcs, x, y, &ra, &dec);
	setdcenter (ra, dec);
	setrefpix (x, y);
	setsecpix (3600.0 * wcs->xinc);
	setsecpix2 (3600.0 * wcs->yinc);
	setrot (wcs->rot);
	recenter = 0;
	imfrac = 0.0;
	goto getfield;
	}

    /* Free image source arrays */
    if (sx) free ((char *)sx);
    if (sy) free ((char *)sy);
    if (sb) free ((char *)sb);
    if (sp) free ((char *)sp);

    return (ret);
}


static void
PrintRes (wcs, nmatch, sx1, sy1, gra1, gdec1, gm1, gnum1, refcat)

struct WorldCoor *wcs;
int	nmatch;		/* Number of image/catalog matches */
double	*sx1, *sy1;	/* Image star pixel coordinates */
double	*gra1, *gdec1;	/* Reference catalog sky coordinates */
double	*gm1;		/* Reference catalog magnitudes */
double	*gnum1;		/* Reference catalog numbers */
int	refcat;		/* Reference catalog code */

{
    int i, goff;
    double gx, gy, dx, dy, dx2, dy2, dxy;
    double sep, rsep, rsep2, dsep, dsep2;
    double dmatch, sra, sdec;
    double sepsum = 0.0;
    double rsepsum = 0.0;
    double rsep2sum = 0.0;
    double dsepsum = 0.0;
    double dsep2sum = 0.0;
    double dxsum = 0.0;
    double dysum = 0.0;
    double dx2sum = 0.0;
    double dy2sum = 0.0;
    double dxysum = 0.0;
    char rstr[32], dstr[32];

    for (i = 0; i < nmatch; i++) {
	wcs2pix (wcs, gra1[i], gdec1[i], &gx, &gy, &goff);
	dx = gx - sx1[i];
	dy = gy - sy1[i];
	dx2 = dx * dx;
	dy2 = dy * dy;
	dxy = dx2 + dy2;
	dxsum = dxsum + dx;
	dysum = dysum + dy;
	dx2sum = dx2sum + dx2;
	dy2sum = dy2sum + dy2;
	dxysum = dxysum + sqrt (dxy);
	pix2wcs (wcs, sx1[i], sy1[i], &sra, &sdec);
	sep = 3600.0 * wcsdist(gra1[i],gdec1[i],sra,sdec);
	rsep = 3600.0 * ((gra1[i]-sra) / cos(degrad(sdec)));
	rsep2 = rsep * rsep;
	dsep = 3600.0 * (gdec1[i] - sdec);
	dsep2 = dsep * dsep;
	sepsum = sepsum + sep;
	rsepsum = rsepsum + rsep;
	dsepsum = dsepsum + dsep;
	rsep2sum = rsep2sum + rsep2;
	dsep2sum = dsep2sum + dsep2;
	ra2str (rstr, 32, gra1[i], 3);
	dec2str (dstr, 32, gdec1[i], 2);
	if (refcat == UAC || refcat == UA1 || refcat == UA2 ||
	    refcat == USAC || refcat == USA1 || refcat == USA2)
	    printf (" %13.8f",gnum1[i]);
	else if (refcat == UJC)
	    printf (" %12.7f",gnum1[i]);
	else if (refcat == GSC)
	    printf (" %9.4f",gnum1[i]);
	else
	    printf (" %9.4f",gnum1[i]);
	printf (" %s %s %5.2f", rstr, dstr, gm1[i]);
	printf (" %6.1f %6.1f %6.2f %6.2f %6.2f\n",
		sx1[i], sy1[i], rsep, dsep, sep);
	}
    dmatch = (double) nmatch;
    dx = dxsum / dmatch;
    dy = dysum / dmatch;
    dx2 = sqrt (dx2sum / dmatch);
    dy2 = sqrt (dy2sum / dmatch);
    dxy = dxysum / dmatch;
    rsep = rsepsum / dmatch;
    dsep = dsepsum / dmatch;
    rsep2 = sqrt (rsep2sum / dmatch);
    dsep2 = sqrt (dsep2sum / dmatch);
    sep = sepsum / dmatch;
    printf ("# Mean  dx= %.4f/%.4f  dy= %.4f/%.4f  dxy= %.4f\n",
	    dx, dx2, dy, dy2, dxy);
    printf ("# Mean dra= %.4f/%.4f  ddec= %.4f/%.4f sep= %.4f\n",
	    rsep, rsep2, dsep, dsep2, sep);

    return;
}


/* Subroutines to initialize various parameters */
void
settolerance (tol)
double tol;
{ tolerance = tol; return; }

void
setmatch (cat)
char *cat;
{
    strcpy (matchcat, cat);
    return;
}

void
setrefcat (cat)
char *cat;
{  RefCat (cat);
    strcpy (refcatname, cat); return; }

void
setreflim (lim1, lim2)
double lim1, lim2;
{ refmag2 = lim2;
  if (lim1 > -2.0) refmag1 = lim1;
  return; }

void
setclass (class)
int class;
{ classd = class;
  return; }

void
setfitwcs (wfit)
int wfit;
{ fitwcs = wfit;
  return; }

void
setfitplate (nc)
int nc;
{ fitplate = nc;
  return; }

void
setplate (plate)
int plate;
{ uplate = plate;
  return; }

void
setnofit ()
{ nofit = 1;
  return; }

void
setfrac (frac0)
double frac0;
{ if (frac0 < 1.0) frac = 1.0 + frac0;
    else frac = frac0;
  return; }

/* Fraction by which to increase dimensions of area to be searched */
void
setimfrac (frac0)
double frac0;
{ if (frac0 > 0.0) imfrac0 = frac0;
  else imfrac0 = 1.0;
  return; }

void
setmaxcat (ncat)
int ncat;
{ if (ncat < 1) maxcat = 25;
  else if (ncat > 200) maxcat = 200;
  else maxcat = ncat;
  return; }

void
setiterate (iter)
int iter;
{ iterate0 = iter;
  return; }

void
setrecenter (recenter)
int recenter;
{ recenter0 = recenter;
  return; }

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
 * Aug 26 1996	Drop unused variable EQ in setfitswcs
 * Aug 28 1996	Improve output format for matched stars
 * Sep  1 1996	Set some defaults in lwcs.h
 * Sep  3 1996	Fix bug to set plate scale from command line
 * Sep  4 1996	Print reference catalog name on separate line from entries
 * Sep 17 1996	Clean up code
 * Oct 15 1996	Break off getfitswcs into separate file
 * Nov 18 1996	Add USNO A catalog searching
 * Nov 18 1996	Write same number into CROAT2 as CROTA1
 * Nov 19 1996	If EPOCH was equinox in original image or not set, set it
 * Dec 10 1996	Make equinox double in getfitswcs call
 *
 * Mar 17 1997	Print found reference stars even when there are not enough
 * Jul 14 1997	If nfit is negative return with header set for nominal WCS
 * Aug  4 1997	Reset nfit limit to 7 for reference pixel fit and fix nfit0 bug
 * Aug 20 1997	Make maximum number of reference stars settable on the command line
 * Aug 28 1997	Print star ID numbers in appropriate format for each catalog
 * Aug 28 1997	Add option to match image to reference stars without fitting WCS
 * Sep  3 1997	Add option to change image dimensions by a fraction
 * Sep  9 1997	Return with default WCS in header only if nfit < -7
 * Sep  9 1997	Print separate right ascension and declination residuals
 * Sep 11 1997	Print average magnitude as well as value of residuals
 * Oct 16 1997	Print same information for image stars as for reference stars
 * Oct 22 1997	Print result of chip rotation as well as optical axis rotation
 * Nov  6 1997	Move nfit entirely to matchstar
 * Nov  6 1997	Rearrange output for IMMATCH use, adding filename argument
 * Nov 14 1997	Increase, not multiply, dimensions by IMFRAC
 * Nov 17 1997	Initialize both magnitude limits
 * Dec  1 1997	Fix bug computing RA separation
 * Dec 16 1997	Fix bug printing no match error message
 *
 * Jan 26 1998	Remove chip rotation code
 * Jan 27 1998	Add switch to use either Calabretta or classic AIPS WCS
 * Jan 29 1998	Fix summary to keep only one reference star per image star
 * Feb  3 1998	Add option to improve WCS by fitting residuals
 * Feb 12 1998	Add USNO SA catalog to reference catalog options
 * Feb 12 1998	Match stars even if less than the minimum if no WCS fit
 * Feb 19 1998	Increase number of reference stars if IMFRAC > 1
 * Feb 22 1998	Fix residual fitting
 * Mar  3 1998	Repeat field search and match after first try
 * Mar  4 1998	Use imfrac on first pass, but not second
 * Mar  5 1998	Correct number of stars used if IMFRAC > 1
 * Mar  6 1998	Add option to use image center for second pass
 * Mar 25 1998	Change residual polynomial fit to full plate-style polynomial
 * Mar 25 1998	Move residual printing to a subroutine
 * Mar 26 1998	Do not fit polynomial until both recenter and iterate done
 * Mar 27 1998	Save plate fit coefficients to FITS header
 * Apr  8 1998	Reset equinox to that of reference catalog
 * Apr 30 1998	Handle prematched star/pixel file
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Sep 17 1998	Allow use of catalogs with other than J2000 coordinates
 * Sep 17 1998	Add coordinate system argument to GetFITSWCS()
 * Sep 28 1998	Add SAO binary format catalogs (SAO, PPM, IRAS, Tycho)
 * Sep 28 1998	Pass system, equinox, and epoch to all catalog search programs
 * Oct  2 1998	Fix arguments in call to GetFITSWCS
 * Oct  7 1998	Set projection to TAN before fitting
 * Oct 16 1998	Add option to read from TDC ASCII format catalog
 * Oct 26 1998	Use passed refcatname and new RefCat subroutine and wcscat.h
 * Oct 26 1998	Add TDC binary catalog option
 * Oct 28 1998	Only search for sources in image once
 * Nov 19 1998	Add catalog name to uacread() call
 * Dec  1 1998	Add version 2.0 of USNO A and SA catalogs
 */
