/* File libwcs/imsetwcs.c
 * November 6, 1997
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
#define UAC		3	/* refcat value for USNO A Star Catalog */
#define USAC		4	/* refcat value for USNO SA Star Catalog */

extern int FindStars ();
extern int TriMatch ();
extern int FocasMatch ();
extern int StarMatch ();
extern int findStars ();
extern int gscread();
extern int uacread();
extern int usaread();
extern int ujcread();
extern int tabread();
extern void MagSortStars ();
extern void FluxSortStars ();
extern struct WorldCoor *GetFITSWCS ();
extern char *getimcat();

extern void fk425(), fk425e(), fk524e();

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
static int refcat = GSC;		/* reference catalog switch */
static char refcatname[32]="GSC";	/* reference catalog name */
static int classd = -1;			/* Guide Star Catalog object classes */
static int uplate = 0;			/* UJ Catalog plate number to use */
static double frac = 1.0;		/* Additional catalog/image stars */
static char wcstype[8]="TAN";		/* WCS projection name */
static int nofit = 0;			/* if =1, do not fit WCS */
static int maxcat = MAXSTARS;		/* Maximum number of catalog stars to use */
static int fitwcs = 1;			/* If 1, fit WCS, else use current WCS */
static double imfrac = 1.0;		/* Multiply image dimensions by this for search */


/* set the C* WCS fields in the input image header based on the given limiting
 * reference mag.
 * Finding stars in the input image and in the reference catalog down to
 * reflim and compute the angle and offsets which result in the best fit.
 * verbose generates extra info on stdout.
 * return 0 if all ok, else -1
 */

int
SetWCSFITS (filename, header, image, verbose)

char	*filename;	/* image file name */
char	*header;	/* FITS header */
char	*image;		/* Image pixels */
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
    double mag1,mag2;
    int minstars;
    int ngmax;
    int nbin, nbytes;
    int ret = 0;
    int is, ig;
    char rstr[32], dstr[32];
    char *imcatname;	/* file name for imagestar catalog, if used */
    struct WorldCoor *wcs=0;	/* WCS structure */

    /* get nominal position and scale */
    wcs = GetFITSWCS (header,verbose,&cra,&cdec,&dra,&ddec,&secpix,&imw,&imh,2000.0);
    if (nowcs (wcs)) {
	ret = 0;
	goto out;
	}

    if (nofit) {
	SetFITSWCS (header, wcs);
	ret = 1;
	goto out;
	}

    dra = dra * imfrac;
    ddec = ddec * imfrac;

    if (reflim > 0.0) {
	mag1 = -2.0;
	mag2 = reflim;
	}
    else {
	mag1 = 0.0;
	mag2 = 0.0;
	}

    /* Allocate arrays for results of reference star search */
    ngmax = maxcat;
    nbytes = ngmax * sizeof (double);
    gnum = (double *) malloc (nbytes);
    gra = (double *) malloc (nbytes);
    gdec = (double *) malloc (nbytes);
    gm = (double *) malloc (nbytes);
    gmb = (double *) malloc (nbytes);
    gc = (int *) malloc (nbytes);

    /* Find the nearby reference stars, in ra/dec */
    if (refcat == UJC)
	ng = ujcread (cra,cdec,dra,ddec,0.0,mag1,mag2,uplate,ngmax,
		      gnum,gra,gdec,gm,gc,verbose);
    else if (refcat == UAC)
	ng = uacread (cra,cdec,dra,ddec,0.0,mag1,mag2,uplate,ngmax,
		      gnum,gra,gdec,gm,gmb,gc,verbose*100);
    else if (refcat == USAC)
	ng = usaread (cra,cdec,dra,ddec,0.0,mag1,mag2,uplate,ngmax,
		      gnum,gra,gdec,gm,gmb,gc,verbose*100);
    else if (refcat == GSC)
	ng = gscread (cra,cdec,dra,ddec,0.0,mag1,mag2,classd,ngmax,
		      gnum,gra,gdec,gm,gc,verbose*100);
    else if (refcatname[0] > 0)
	ng = tabread (refcatname,cra,cdec,dra,ddec,0.0,mag1,mag2,ngmax,
		      gnum,gra,gdec,gm,gc,verbose);
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
    MagSortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, ng);

    /* Use only the brightest reference stars */
    nbg = maxcat;
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

    if (verbose) {
	printf ("%s:\n",refcatname);
	for (ig = 0; ig < ng; ig++) {
	    ra2str (rstr, gra[ig], 3);
	    dec2str (dstr, gdec[ig], 2);
	    if (refcat == GSC)
		printf (" %9.4f",gnum[ig]);
	    else if (refcat == UAC || refcat == USAC)
		printf (" %13.8f",gnum[ig]);
	    else if (refcat == UJC)
		printf (" %12.7f",gnum[ig]);
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
	    fprintf (stderr,"Need >= %d reference stars; found none\n",
							minstars);
	else
	    fprintf (stderr, "Need >= %d reference stars; found only %d\n",
							minstars, ng);
	ret = 0;
	goto out;
	}

    /* Discover star-like things in the image, in pixels */
    ns = FindStars (header, image, &sx, &sy, &sb, &sp, verbose);
    if (ns < minstars) {
	fprintf (stderr, "Need at least %d image stars but only found %d\n",
							minstars, ns);
	ret = 0;
	goto out;
	}

    if (fitwcs) {

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
	    char rastr[32], decstr[32];
	    double xmag, mdiff, ra, dec;
	    for (is = 0; is < nbs; is++) {
		pix2wcs (wcs, sx[is], sy[is], &ra, &dec);
		ra2str (rastr, ra, 3);
		dec2str (decstr, dec, 2);
		xmag = -2.5 * log10 (sb[is]);
		if (!is) mdiff = gm[0] - xmag;
		xmag = xmag + mdiff;
		printf ("%4d %s %s %6.2f %6.1f %6.1f %d\n",
			is+1, rastr, decstr, xmag, sx[is], sy[is], sp[is]);
		}
	    }

	/* Check offsets between all pairs of image stars and reference stars */
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

    if (verbose || !fitwcs) {
	double ra,dec;
	double x,y,dx,dy,dx2,dy2,dxy, rsep,dsep,sep,xysep, tol2, rsep2, dsep2;
	double dxsum=0.0, dysum=0.0, dxysum=0.0, dx2sum=0.0, dy2sum=0.0;
	double sepsum=0.0, rsepsum=0.0, dsepsum=0.0, rsep2sum=0.0, dsep2sum=0.0;
	int offscl, nmatch=0;

	if (wcs->mrot == 0.0)
	    printf ("# Arcsec/Pixel: %.6f %.6f  Rotation: %.6f degrees\n",
		    3600.0*wcs->xinc, 3600.0*wcs->yinc, wcs->rot);
	else
	    printf ("# Arcsec/Pixel: %.5f %.5f  Rotation: %.5f degrees  Chip Rot: %.5f degrees\n",
		    3600.0*wcs->xinc, 3600.0*wcs->yinc, wcs->rot, wcs->mrot);
	ra2str (rstr, wcs->xref, 3);
	dec2str (dstr, wcs->yref, 2);
	printf ("# Optical axis: %s  %s J2000 at (%.2f,%.2f)\n",
		rstr,dstr, wcs->xrefpix, wcs->yrefpix);
	ra = wcs->xref;
	dec = wcs->yref;
	(void)fk524e (&ra, &dec, wcs->epoch);
	ra2str (rstr, ra, 3);
	dec2str (dstr, dec, 2);
	printf ("# Optical axis: %s  %s B1950 at (%.2f,%.2f)\n",
		rstr,dstr, wcs->xrefpix, wcs->yrefpix);

    /* Find star matches for this offset and print them */
	tol2 = tolerance * tolerance;
	for (ig = 0; ig < ng; ig++) {
	    wcs2pix (wcs, gra[ig], gdec[ig], &x, &y, &offscl);
	    if (!offscl) {
		for (is = 0; is < ns; is++) {
		    dx = x - sx[is];
		    dy = y - sy[is];
		    dx2 = dx * dx;
		    dy2 = dy * dy;
		    dxy = dx2 + dy2;
		    if (dxy < tol2)
			nmatch++;
		    }
		}
	    }
	imcatname = getimcat ();
	if (imcatname == NULL)
	    imcatname = filename;
	if (nmatch > 0) {
	    printf ("# %d matches between %s and %s:\n",
		    nmatch, refcatname, imcatname);
	    for (ig = 0; ig < ng; ig++) {
		wcs2pix (wcs, gra[ig], gdec[ig], &x, &y, &offscl);
		if (!offscl) {
		    for (is = 0; is < ns; is++) {
			dx = x - sx[is];
			dy = y - sy[is];
			dx2 = dx * dx;
			dy2 = dy * dy;
			dxy = dx2 + dy2;
			if (dxy < tol2) {
			    dxsum = dxsum + dx;
			    dysum = dysum + dy;
			    dx2sum = dx2sum + dx2;
			    dy2sum = dy2sum + dy2;
			    dxysum = dxysum + sqrt (dxy);
			    pix2wcs (wcs, sx[is], sy[is], &sra, &sdec);
			    xysep = sqrt (dxy);
			    sep = 3600.0 * wcsdist (gra[ig],gdec[ig],sra,sdec);
			    rsep = 3600.0 * ((gra[ig]-sra) / cos(degrad(sdec)));
			    rsep2 = rsep * rsep;
			    dsep = 3600.0 * (gdec[ig] - sdec);
			    dsep2 = dsep * dsep;
			    sepsum = sepsum + sep;
			    rsepsum = rsepsum + rsep;
			    dsepsum = dsepsum + dsep;
			    rsep2sum = rsep2sum + rsep2;
			    dsep2sum = dsep2sum + dsep2;
			    ra2str (rstr, gra[ig], 3);
			    dec2str (dstr, gdec[ig], 2);
			    if (refcat == GSC)
				printf (" %9.4f",gnum[ig]);
			    else if (refcat == UAC || refcat == USAC)
				printf (" %13.8f",gnum[ig]);
			    else if (refcat == UJC)
				printf (" %12.7f",gnum[ig]);
			    else
				printf (" %9.4f",gnum[ig]);
			    printf (" %s %s %5.2f", rstr, dstr, gm[ig]);
			    printf (" %6.1f %6.1f %6.2f %6.2f %6.2f\n",
			        sx[is], sy[is], rsep, dsep, sep);
			    }
		        }
		    }
	        }
	    dx = dxsum / (double)nmatch;
	    dy = dysum / (double)nmatch;
	    dx2 = sqrt (dx2sum / (double)nmatch);
	    dy2 = sqrt (dy2sum / (double)nmatch);
	    dxy = dxysum / (double)nmatch;
	    rsep = rsepsum / (double)nmatch;
	    dsep = dsepsum / (double)nmatch;
	    rsep2 = sqrt (rsep2sum / (double)nmatch);
	    dsep2 = sqrt (dsep2sum / (double)nmatch);
	    sep = sepsum / (double)nmatch;
	    printf ("# Mean  dx= %.2f/%.2f  dy= %.2f/%.2f  dxy= %.2f\n",
		    dx, dx2, dy, dy2, dxy);
	    printf ("# Mean dra= %.2f/%.2f  ddec= %.2f/%.2f sep= %.2f\n",
		    rsep, rsep2, dsep, dsep2, sep);
	    }
	else
	    fprintf (stderr, "SetWCSFITS: No matches between %s and %s:\n",
		    nmatch, refcatname, imcatname);
	}

    ret = 1;

    out:
    if (gra) free ((char *)gra);
    if (gdec) free ((char *)gdec);
    if (gm) free ((char *)gm);
    if (gmb) free ((char *)gmb);
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
    if (hgetr8 (header, "WEPOCH", &ep))
	hputr8 (header, "EPOCH", wcs->equinox);
    else if (!hgetr8 (header, "EPOCH", &ep))
	hputr8 (header, "EPOCH", wcs->equinox);
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
    hputnr8 (header, "CROTA2", 3, wcs->rot);
    return;
}


/* Subroutines to initialize various parameters */
void
settolerance (tol)
double tol;
{ tolerance = tol; return; }

void
setrefcat (cat)
char *cat;
{  if (strcmp(cat,"gsc")==0 || strcmp(cat,"GSC")==0)
	refcat = GSC;
   else if (strcmp(cat,"ujc")==0 || strcmp(cat,"UJC")==0)
	refcat = UJC;
   else if (strcmp(cat,"uac")==0 || strcmp(cat,"UAC")==0)
	refcat = UAC;
    else
	refcat = 0;
    strcpy (refcatname, cat); return; }

void
setwcstype (type)
char *type;
{ strcpy (wcstype, type); return; }

void
setreflim (lim)
double lim;
{ reflim = lim;
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

void
setimfrac (frac0)
double frac0;
{ if (frac0 > 0.0) imfrac = frac0;
  else imfrac = 1.0;
  return; }

void
setmaxcat (ncat)
int ncat;
{ if (ncat < 1) maxcat = 25;
  else if (ncat > 200) maxcat = 200;
  else maxcat = ncat;
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
 */
