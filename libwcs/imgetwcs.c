/* File libwcs/imgetwcs.c
 * December 8, 1997
 * By Doug Mink, remotely based on UIowa code
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fitshead.h"
#include "wcs.h"
#include "lwcs.h"

extern void fk425e(), fk524e(), fk425(), fk524(), fk5prec(), fk4prec();

/* Get the C* WCS fields in  a FITS header based on a reference catalog
 * do it by finding stars in the image and in the reference catalog and
 * finding the rotation and offsets which result in a best-fit.
 * verbose generates extra info on stderr.
 * try using deeper reference star catalog searches if there is trouble.
 * return 1 if all ok, else 0
 */

/* These parameters can be set on the command line */
static double secpix0 = PSCALE;		/* Set image scale--override header */
static double rot0 = 361.0;		/* Initial image rotation */
static int fk4 = 0;			/* Command line center is FK4 */
static double ra0 = -99.0;		/* Initial center RA in degrees */
static double dec0 = -99.0;		/* Initial center Dec in degrees */
static double rad0 = 10.0;		/* Search box radius in arcseconds */
static double xref0 = -99999.0;		/* Reference pixel X coordinate */
static double yref0 = -99999.0;		/* Reference pixel Y coordinate */


/* Set a nominal world coordinate system from image header info.
 * If the image center is not FK5 (J2000) equinox, convert it
 * Return a WCS structure if OK, else return NULL
 */

struct WorldCoor *
GetFITSWCS (header, verbose, cra, cdec, dra, ddec, secpix, wp, hp, eq2)

char	*header;	/* Image FITS header */
int	verbose;	/* Extra printing if =1 */
double	*cra;		/* Center right ascension in degrees (returned) */
double	*cdec;		/* Center declination in degrees (returned) */
double	*dra;		/* Right ascension half-width in degrees (returned) */
double	*ddec;		/* Declination half-width in degrees (returned) */
double	*secpix;	/* Arcseconds per pixel (returned) */
int	*wp;		/* Image width in pixels (returned) */
int	*hp;		/* Image height in pixels (returned) */
double	eq2;		/* Equinox to return (0=keep equinox) */
{
    int nax;
    int equinox, eqcoor;
    double eq1, epoch, xref, yref, degpix;
    struct WorldCoor *wcs;
    int eqref;
    char rstr[32],dstr[32];

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

    /* Set plate center from command line, if it is there */
    if (ra0 > -99.0 && dec0 > -99.0) {
	hputnr8 (header, "CRVAL1" ,8,ra0);
	hputnr8 (header, "CRVAL2" ,8,dec0);
	hputra (header, "RA", ra0);
	hputdec (header, "DEC", dec0);
	hputc (header, "CTYPE1", "RA---TAN");
	hputc (header, "CTYPE2", "DEC--TAN");
	if (fk4) {
	    hputi4 (header, "EPOCH", 1950);
	    hputi4 (header, "EQUINOX", 1950);
	    }
	else {
	    hputi4 (header, "EPOCH", 2000);
	    hputi4 (header, "EQUINOX", 2000);
	    }
	if (hgetr8 (header, "SECPIX", secpix)) {
	    degpix = *secpix / 3600.0;
	    hputnr8 (header, "CDELT1", 8, -degpix);
	    hputnr8 (header, "CDELT2", 8, degpix);
	    } 
	}

    /* Set reference pixel from command line, if it is there */
    if (xref0 > -99999.0 && yref0 > -99999.0) {
	hputr8 (header, "CRPIX1", xref0);
	hputr8 (header, "CRPIX2", yref0);
	}
    else if (hgetr8 (header, "CRPIX1", &xref) < 1) {
	xref = (double) *wp / 2.0;
	yref = (double) *hp / 2.0;
	hputnr8 (header, "CRPIX1", 3, xref);
	hputnr8 (header, "CRPIX2", 3, yref);
	}

    /* Set plate scale from command line, if it is there */
    if (secpix0 > 0.0) {
	*secpix = secpix0;
	hputnr8 (header, "SECPIX", 5, *secpix);
	degpix = *secpix / 3600.0;
	hputnr8 (header, "CDELT1", 8, -degpix);
	hputnr8 (header, "CDELT2", 8, degpix);
	}

    /* Set rotation angle from command line, if it is there */
    if (rot0 < 361.0) {
	hputnr8 (header, "CROTA1", 5, rot0);
	hputnr8 (header, "CROTA2", 5, rot0);
	}

    /* Initialize WCS structure from FITS header */
    wcs = wcsinit (header);

    /* If incomplete WCS in header, drop out */
    if (nowcs (wcs)) {
	if (verbose)
	    fprintf (stderr,"Insufficient information for initial WCS\n");
	return (NULL);
	}

    /* Set flag to get appropriate equinox for catalog search */
    equinox = (int) wcs->equinox;
    eq1 = wcs->equinox;
    if (eq2 == 0.0)
	eqref = equinox;
    else
	eqref = (int) eq2;
    if (eqref == 1950)
	wcsoutinit (wcs, "FK4");
    else
	wcsoutinit (wcs, "FK5");

    /* Get center and size for catalog searching */
    wcssize (wcs, cra, cdec, dra, ddec);

    /* Set reference pixel to center of image if it has not been set */
    if (wcs->xref == 0.0 && wcs->yref == 0.0) {
	wcs->xref = *cra;
	wcs->yref = *cdec;
	if (wcs->xrefpix == 0.0 && wcs->yrefpix == 0.0) {
	    wcs->xrefpix = (double) wcs->nxpix * 0.5;
	    wcs->yrefpix = (double) wcs->nypix * 0.5;
	    }
	wcs->xinc = *dra / wcs->xrefpix;
	wcs->yinc = *ddec / wcs->yrefpix;
	hchange (header,"PLTRAH","PLT0RAH");
	wcs->plate_fit = 0;
	}

    /* Compute plate scale to return if it was not set on the command line */
    if (secpix0 <= 0.0)
	*secpix = 3600.0 * 2.0 * *ddec / (double) *hp;

    /* Reset to reference pixel position from command line, if there is one */
    if (ra0 > -99.0 && dec0 > -99.0 && verbose) {
	ra2str (rstr, ra0, 3);
        dec2str (dstr, dec0, 2);
	if (fk4)
	    printf ("Reference pixel (%.2f,%.2f) %s %s B1950\n",
		    wcs->xrefpix, wcs->yrefpix, rstr, dstr);
	else
	    printf ("Reference pixel (%.2f,%.2f) %s %s J2000\n",
		    wcs->xrefpix, wcs->yrefpix, rstr, dstr);
	}

    /* Image size for catalog search */
    if (verbose) {
	char rstr[64], dstr[64];
	ra2str (rstr, *cra, 3);
	dec2str (dstr, *cdec, 2);
	if (eqref == 2000)
	    printf ("Search at %s %s J2000", rstr, dstr);
	else
	    printf ("Search at %s %s B1950", rstr, dstr);
	ra2str (rstr, *dra, 3);
	dec2str (dstr, *ddec, 2);
	printf (" +- %s %s\n", rstr, dstr);
	printf ("Image width=%d height=%d, %g arcsec/pixel\n",
				*wp, *hp, *secpix);
	}

    return (wcs);
}


/* Get a center and radius for a search area.  If the image center is not
 * given in the equinox of the reference catalog, convert it.
 * Return 0 if OK, else -1
 */

int
GetArea (verbose, eqref, cra, cdec, dra, ddec)

int	verbose;	/* Extra printing if =1 */
int	eqref;		/* Equinox of reference catalog */
double	*cra;		/* Center right ascension in degrees (returned) */
double	*cdec;		/* Center declination in degrees (returned) */
double	*dra;		/* Right ascension half-width in degrees (returned) */
double	*ddec;		/* Declination half-width in degrees (returned) */
{
    int eqcoor;
    char rstr[32], dstr[32];

    /* Set plate center from command line, if it is there */
    if (ra0 < 0.0 && dec0 < -90.0) {
	if (verbose)
	    fprintf (stderr, "GetArea: Illegal center, ra= %.5f, dec= %.5f\n",
		     ra0,dec0);
	return (-1);
	}
    else {
	*cra = ra0;
	*cdec = dec0;
	}

    /* If coordinate equinox not reference catalog equinox, convert */
    if (fk4)
	eqcoor = 1950;
    else
	eqcoor = 2000;
    if (eqcoor != eqref) {
	if (eqref > 1980)
	    fk425 (cra, cdec);
	else
	    fk524 (cra, cdec);
	}

    /* Set search box radius from command line, if it is there */
    if (rad0 > 0.0) {
	*ddec = rad0 / 3600.0;
	if (*cdec < 90.0 && *cdec > -90.0)
	    *dra = *ddec / cos (degrad (*cdec));
	else
	    *dra = 180.0;
	}
    else {
	if (verbose)
	    fprintf (stderr, "GetArea: Illegal radius, rad= %.5f\n",rad0);
	return (-1);
	}


    if (verbose) {
	if (eqcoor != eqref) {
	    ra2str (rstr, ra0, 3);
            dec2str (dstr, dec0, 2);
	    if (fk4)
		fprintf (stderr,"Center:  %s   %s (B1950)\n", rstr, dstr);
	    else
		fprintf (stderr,"Center:  %s   %s (J2000)\n", rstr, dstr);
	    }
	ra2str (rstr, *cra, 3);
        dec2str (dstr, *cdec, 2);
	if (eqref == 2000)
	    fprintf (stderr,"Center:  %s   %s (J2000)\n", rstr, dstr);
	else
	    fprintf (stderr,"Center:  %s   %s (B1950)\n", rstr, dstr);
	ra2str (rstr, *dra * 2.0, 2); 
	dec2str (dstr, *ddec * 2.0, 2); 
	fprintf (stderr,"Area:    %s x %s\n", rstr, dstr);
	}

    return (0);
}


void setrot (rot)
double rot;
{
    rot0 = rot;
    return;
}

void setsecpix (secpix)
double secpix;
{
    secpix0 = secpix;
    return;
}

void setfk4 ()
{
    fk4 = 1;
    return;
}

void setcenter (rastr, decstr)
char *rastr, *decstr;
{
    ra0 = str2ra (rastr);
    dec0 = str2dec (decstr);
    return;
}

void setrefpix (x, y)
double x, y;
{
    xref0 = x;
    yref0 = y;
    return;
}

void setradius (rad)
double rad;
{
    rad0 = rad;
    return;
}


/* Feb 29 1996	New program
 * May 23 1996	Use pre-existing WCS for center, if it is present
 * May 29 1996	Simplify program by always using WCS structure
 * Jun 12 1996	Be more careful with nominal WCS setting
 * Jul  3 1996	Set epoch from old equinox if not already set
 * Jul 19 1996	Set image center in WCS if DSS WCS
 * Aug  5 1996	Check for SECPIX1 as well as SECPIX
 * Aug  7 1996	Save specified number of decimal places in header parameters
 * Aug  7 1996	Rename old center parameters
 * Aug 26 1996	Decrease default pixel tolerance from 20 to 10
 * Sep  1 1996	Set plate scale default in lwcs.h
 * Sep  3 1996	Fix bug to set plate scale from command line
 * Oct 15 1996	Break off from imsetwcs.c
 * Oct 16 1996	Clean up center setting so eqref is used
 * Oct 17 1996	Do not print error messages unless verbose is set
 * Oct 30 1996	Keep equinox from image if EQREF is zero
 * Nov  1 1996	Declare undeclared subroutines; remove unused variables
 * Nov  4 1996	Add reference pixel and projection to wcsset() call
 * Nov 14 1996	Add GetLimits() to deal with search limits around the poles
 * Nov 15 1996	Drop GetLimits(); code moved to individual catalog routines
 * Dec 10 1996	Fix precession and make equinox double

 * Feb 19 1997	If eq2 is 0, use equinox of image
 * Feb 24 1997	Always convert center to output equinox (bug fix)
 * Mar 20 1997	Declare EQ2 double instead of int, fixing a bug
 * Jul 11 1997	Allow external (command line) setting of reference pixel coords
 * Sep 26 1997	Set both equinox and epoch if input center coordinates
 * Nov  3 1997	Separate WCS reference pixel from search center
 * Dec  8 1997	Set CDELTn using SECPIX if it is in the header
 *
 * Jan  6 1997	Do not print anything unless verbose is set
 */
