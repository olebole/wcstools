/* File platepos.c
 * May 5, 1998
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics

 * Module:	platepos.c (Plate solution WCS conversion
 * Purpose:	Compute WCS from plate fit
 * Subroutine:	platepos() converts from pixel location to RA,Dec 
 * Subroutine:	platepix() converts from RA,Dec to pixel location   
 */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "wcs.h"

int
platepos (xpix, ypix, wcs, xpos, ypos)

/* Routine to determine accurate position for pixel coordinates */
/* returns 0 if successful otherwise 1 = angle too large for projection; */
/* based on amdpos() from getimage */

/* Input: */
double	xpix;		/* x pixel number  (RA or long without rotation) */
double	ypix;		/* y pixel number  (dec or lat without rotation) */
struct WorldCoor *wcs;	/* WCS parameter structure */

/* Output: */
double	*xpos;		/* Right ascension or longitude in degrees */
double	*ypos;		/* Declination or latitude in degrees */

{
    double x, y, x2, y2, x3, y3, r2;
    double xi, xir, eta, etar, raoff, ra, dec, ra0, dec0;
    double twopi = 6.28318530717959;
    double ctan, ccos;
    int ncoeff1 = wcs->ncoeff1;
    int ncoeff2 = wcs->ncoeff2;

    /* Convert from image pixels to pixels from reference pixel */
    x = xpix - wcs->crpix[0];
    y = ypix - wcs->crpix[1];

    /* Convert from pixels to angle from reference */
    polyrev (x, y, wcs, xi, eta);

    /* Convert to radians */
    xir = degrad (xi);
    etar = degrad (eta);

    /* Convert to RA and Dec */
    ra0 = degrad (wcs->crval[0]);
    dec0 = degrad (wcs->crval[1]);
    ctan = tan (dec0);
    ccos = cos (dec0);
    raoff = atan2 (xir / ccos, 1.0 - etar * ctan);
    ra = raoff + ra0;
    if (ra < 0.0) ra = ra + twopi;
    *xpos = raddeg (ra);

    dec = atan (cos (raoff) / ((1.0 - (etar * ctan)) / (etar + ctan)));
    *ypos = raddeg (dec);
    return 0;
}


int
platepix (xpos, ypos, wcs, xpix, ypix)

/* Routine to determine pixel coordinates for sky position */
/* returns 0 if successful otherwise 1 = angle too large for projection; */
/* based on amdinv() from getimage */

/* Input: */
double	xpos;		/* Right ascension or longitude in degrees */
double	ypos;		/* Declination or latitude in degrees */
struct WorldCoor *wcs;	/* WCS parameter structure */

/* Output: */
double	*xpix;		/* x pixel number  (RA or long without rotation) */
double	*ypix;		/* y pixel number  (dec or lat without rotation) */

{
    double div,xi,eta,x,y,xy,x2,y2,x2y,y2x,x3,y3,r2,dx,dy;
    double tdec,ctan,ccos,traoff, craoff, etar, xir;
    double f,fx,fy,g,gx,gy;
    double ra0, dec0, ra, dec;
    double tolerance = 0.0000005;
    int    max_iterations = 50;
    int    i;
    int	ncoeff1 = wcs->ncoeff1;
    int	ncoeff2 = wcs->ncoeff2;
    double xr, yr; 	/* position in radians */

    /* Convert RA and Dec in radians to reference system of image */
    ra = degrad (xpos);
    dec = degrad (ypos);
    tdec = tan (dec);
    ra0 = degrad (wcs->crval[0]);
    dec0 = degrad (wcs->crval[1]);
    ctan = tan (dec0);
    ccos = cos (dec0);
    traoff = tan (ra - ra0);
    craoff = cos (ra - ra0);
    etar = (1.0 - ctan * craoff / tdec) / (ctan + (craoff / tdec));
    xir = traoff * ccos * (1.0 - (etar * ctan));
    xi = raddeg (xir);
    eta = raddeg (etar);

    polyfwd (xi, eta, wcs, &x, &y);

    /* Convert from image pixels from relative pixels */
    *xpix = x + wcs->crpix[0];
    *ypix = y + wcs->crpix[1];

    /* If position is off of the image, return offscale code */
    if (*xpix < 0.5 || *xpix > wcs->nxpix+0.5)
	return -1;
    if (*ypix < 0.5 || *ypix > wcs->nypix+0.5)
	return -1;

    return 0;
}


/* Set plate fit coefficients in structure from arguments */
int
SetPlate (wcs, ncoeff1, ncoeff2, coeff)

struct WorldCoor *wcs;  /* World coordinate system structure */
int	ncoeff1;		/* Number of coefficients for x */
int	ncoeff2;		/* Number of coefficients for y */
double	*coeff;		/* Plate fit coefficients */

{
    int i;

    if (nowcs (wcs) || (ncoeff1 < 1 && ncoeff2 < 1))
	return 1;

    wcs->ncoeff1 = ncoeff1;
    wcs->ncoeff2 = ncoeff2;
    wcs->prjcode = WCS_PLT;

    for (i = 0; i < 20; i++) {
	if (i < ncoeff1)
	    wcs->x_coeff[i] = coeff[i];
	else
	    wcs->x_coeff[i] = 0.0;
	}

    for (i = 0; i < 20; i++) {
	if (i < ncoeff2)
	    wcs->y_coeff[i] = coeff[ncoeff1+i];
	else
	    wcs->y_coeff[i] = 0.0;
	}
    return 0;
}


/* Return plate fit coefficients from structure in arguments */
int
GetPlate (wcs, ncoeff1, ncoeff2, coeff)

struct WorldCoor *wcs;  /* World coordinate system structure */
int	*ncoeff1;	/* Number of coefficients for x */
int	*ncoeff2;	/* Number of coefficients for y) */
double	*coeff;		/* Plate fit coefficients */

{
    int i;

    if (nowcs (wcs))
	return 1;

    *ncoeff1 = wcs->ncoeff1;
    *ncoeff2 = wcs->ncoeff2;

    for (i = 0; i < *ncoeff1; i++)
	coeff[i] = wcs->x_coeff[i];

    for (i = 0; i < *ncoeff2; i++)
	coeff[*ncoeff1+i] = wcs->y_coeff[i];

    return 0;
}


/* Set FITS header plate fit coefficients from structure */
void
SetFITSPlate (header, wcs)

char    *header;        /* Image FITS header */
struct WorldCoor *wcs;  /* WCS structure */

{
    char keyword[16];
    int i;

    for (i = 0; i < wcs->ncoeff1; i++) {
	sprintf (keyword,"CO1_%d",i+1);
	hputnr8 (header, keyword, -15, wcs->x_coeff[i]);
	}
    for (i = 0; i < wcs->ncoeff2; i++) {
	sprintf (keyword,"CO2_%d",i+1);
	hputnr8 (header, keyword, -15, wcs->y_coeff[i]);
	}
    return;
}

/* Mar 27 1998	New subroutines for direct image pixel <-> sky polynomials
 * Apr 10 1998	Make terms identical for both x and y polynomials
 * Apr 10 1998	Allow different numbers of coefficients for x and y
 * Apr 16 1998	Drom NCOEFF header parameter
 * Apr 28 1998  Change projection flags to WCS_*
 */
