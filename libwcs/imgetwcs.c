/* File libwcs/imgetwcs.c
 * November 15, 1996
 * By Doug Mink, based on UIowa code
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "fitshead.h"
#include "wcs.h"
#include "lwcs.h"

extern void fk425e(), fk524e(), fk425(), fk524();

/* Get the C* WCS fields in  a FITS header based on a reference catalog
 * do it by finding stars in the image and in the reference catalog and
 * finding the rotation and offsets which result in a best-fit.
 * verbose generates extra info on stderr.
 * try using deeper reference star catalog searches if there is trouble.
 * return 1 if all ok, else 0
 */

/* These parameters can be set on the command line */
static double secpix0 = PSCALE;		/* Set image scale--override header */
static double rot0 = 0.0;		/* Initial image rotation */
static int fk4 = 0;			/* Command line center is FK4 */
static double ra0 = -99.0;		/* Initial center RA in degrees */
static double dec0 = -99.0;		/* Initial center Dec in degrees */
static double rad0 = 10.0;		/* Search box radius in arcseconds */


/* Set a nominal world coordinate system from image header info.
 * If the image center is not FK5 (J2000) equinox, convert it
 * Return a WCS structure if OK, else return NULL
 */

struct WorldCoor *
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
int	eqref;		/* Equinox of reference catalog (0=keep equinox) */
{
    int nax;
    int equinox, eqcoor;
    double epoch, xref, yref;
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

    /* Set plate scale from command line, if it is there */
    if (secpix0 > 0.0)
	hputnr8 (header, "SECPIX", 5, secpix0);

    /* Set plate center from command line, if it is there */
    if (ra0 > -99.0 && dec0 > -99.0) {
	hputra (header, "RA", ra0);
	hputdec (header, "DEC", ra0);
	if (fk4)
	    hputi4 (header, "EPOCH", 1950);
	else
	    hputi4 (header, "EPOCH", 2000);
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
	if (eqref != 0.0 && wcs->equinox != eqref) {
	    if (eqref > 1980) {
		fk425e (cra, cdec, wcs->epoch);
		wcsshift (wcs, *cra, *cdec, "FK5");
		}
	    else {
		fk524e (cra, cdec, wcs->epoch);
		wcsshift (wcs, *cra, *cdec, "FK4");
		}
	    }
	}

    /* Otherwise use nominal center from RA and DEC fields */
    else {
	*cra = 0.0;
	*cdec = 0.0;
	if (hgetra (header, "RA", cra) == 0) {
	    if (verbose)
		fprintf (stderr, "No RA field\n");
	    return (NULL);
	    }
	if (hgetdec (header, "DEC", cdec) == 0) {
	    if (verbose)
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
	    if (eqref > 1980)
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
		    if (hgetr8 (header, "PLTSCALE", secpix) == 0) {
			if (verbose)
			    fprintf (stderr, "Cannot find SECPIX in header\n");
			}
		    return (NULL);
		    }
		}
	    }
	}

    /* Set WCS structure if it has not already been set */
    if (nowcs (wcs)) {
	xref = (double)*wp * 0.5;
	yref = (double)*hp * 0.5;
	wcs = wcsset (*cra, *cdec, *secpix, xref, yref, *wp, *hp, rot0,
		      eqref, epoch ,"-TAN");
	}

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
	    if (eqref > 1980)
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
 */
