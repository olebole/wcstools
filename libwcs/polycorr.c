/* File saoimage/wcslib/polycorr.c
 * May 4, 1998
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics

 * Module:	polycorr.c 
 * Purpose:	Polynomial correction for WCS
 * Subroutine:	polyfwd() converts to pixel location
 * Subroutine:	polyrev() converts from pixel location   

    These functions are based on the astrmcal.c portion of GETIMAGE by
    J. Doggett and the documentation distributed with the Digital Sky Survey.

*/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include "wcs.h"

/* Routine to determine accurate position for pixel coordinates */
/* based on amdpos() from getimage */

int
polyrev (x, y, wcs, xc, yc)

double	x, y;		/* Uncorrected coordinates */
struct WorldCoor *wcs;	/* WCS parameter structure */
double	*xc, *yc;	/* Corrected coordinates (output)*/

{
    double x, y, x2, y2, x3, y3, r2;
    int ncoeff1 = wcs->ncoeff1;
    int ncoeff2 = wcs->ncoeff2;
    double xi, xir, eta, etar, raoff, ra, dec, ra0, dec0;
    double ctan, ccos;

    x = xpix - wcs->crpix[0];
    y = ypix - wcs->crpix[1];

    x2 = x * x;
    y2 = y * y;
    x3 = x * x2;
    y3 = y * y2;
    r2 = x2 + y2;

    /*  Compute xc,yc coordinates in degrees from x,y and plate model */
    *xc = wcs->x_coeff[ 0]	+ wcs->x_coeff[ 1]*x +
	  wcs->x_coeff[ 2]*y	+ wcs->x_coeff[ 3]*x2 +
	  wcs->x_coeff[ 4]*y2	+ wcs->x_coeff[ 5]*x*y;

    if (ncoeff1 > 6)
	  *xc = *xc + wcs->x_coeff[ 6]*x3	+ wcs->x_coeff[ 7]*y3;

    if (ncoeff1 > 8) {
	*xc = *xc + wcs->x_coeff[ 8]*x2*y	+ wcs->x_coeff[ 9]*x*y2 +
		    wcs->x_coeff[10]*(r2)	+ wcs->x_coeff[11]*x*r2 +
		    wcs->x_coeff[12]*y*r2;
	}

    *yc = wcs->y_coeff[ 0]	+ wcs->y_coeff[ 1]*x +
	  wcs->y_coeff[ 2]*y	+ wcs->y_coeff[ 3]*x2 +
	  wcs->y_coeff[ 4]*y2	+ wcs->y_coeff[ 5]*x*y;

    if (ncoeff2 > 6)
	*yc = *yc + wcs->y_coeff[ 6]*x3	+ wcs->y_coeff[ 7]*y3;

    if (ncoeff2 > 8) {
	*yc = *yc + wcs->y_coeff[ 8]*x2*y + wcs->y_coeff[ 9]*y2*x +
		    wcs->y_coeff[10]*r2   + wcs->y_coeff[11]*x*r2 +
		    wcs->y_coeff[12]*y*r2;
	}

    return 0;
}


/* Routine to determine coordinates from corrected coordinates */
/* based on amdinv() from getimage */

int
polyfwd (xc, yc, wcs, x, y)

double	xc, yc;		/* Corrected coordinates */
struct WorldCoor *wcs;	/* WCS parameter structure */
double	*x, *y;		/* Uncorrected coordinates (output) */

{
    double div,xi,eta,x,y,xy,x2,y2,x2y,y2x,x3,y3,r2,dx,dy;
    double f,fx,fy,g,gx,gy;
    double tolerance = 0.0000005;
    int    max_iterations = 50;
    int    i;
    int	ncoeff1 = wcs->ncoeff1;
    int	ncoeff2 = wcs->ncoeff2;

    /* Set initial value for x,y */
    x = (xc - wcs->x_coeff[0]) / wcs->x_coeff[1];
    y = (yc - wcs->y_coeff[0]) / wcs->y_coeff[2];

    /* Iterate by Newton's method */
    for (i = 0; i < max_iterations; i++) {

	/* X plate model */
	xy = x * y;
	x2 = x * x;
	y2 = y * y;
	x3 = x2 * x;
	y3 = y2 * y;
	x2y = x2 * y;
	y2x = y2 * x;
	r2 = x2 + y2;

	f = wcs->x_coeff[0]	+ wcs->x_coeff[1]*x +
	    wcs->x_coeff[2]*y	+ wcs->x_coeff[3]*x2 +
	    wcs->x_coeff[4]*y2	+ wcs->x_coeff[5]*xy;

	/*  Derivative of X model wrt x */
	fx = wcs->x_coeff[1]	+ wcs->x_coeff[3]*2.0*x +
	     wcs->x_coeff[5]*y;

	/* Derivative of X model wrt y */
	fy = wcs->x_coeff[2]	+ wcs->x_coeff[4]*2.0*y +
	     wcs->x_coeff[5]*x;

	if (ncoeff1 > 6) {
	    f = f + wcs->x_coeff[6]*x3	+ wcs->x_coeff[7]*y3;
	    fx = fx + wcs->x_coeff[6]*3.0*x2;
	    fy = fy + wcs->x_coeff[7]*3.0*y2;
	    }

	if (ncoeff1 > 8) {
	    f = f +
		wcs->x_coeff[8]*x2y	+ wcs->x_coeff[9]*y2x +
		wcs->x_coeff[10]*r2 + wcs->x_coeff[11]*x*r2 +
		wcs->x_coeff[12]*y*r2;

	    fx = fx +	wcs->x_coeff[8]*2.0*xy + 
			wcs->x_coeff[9]*y2 +
	 		wcs->x_coeff[10]*2.0*x +
			wcs->x_coeff[11]*(3.0*x2+y2) +
			wcs->x_coeff[12]*2.0*xy;

	    fy = fy +	wcs->x_coeff[8]*x2 +
			wcs->x_coeff[9]*2.0*xy +
			wcs->x_coeff[10]*2.0*y +
			wcs->x_coeff[11]*2.0*xy +
			wcs->x_coeff[12]*(3.0*y2+x2);
	    }

	/* Y plate model */
	g = wcs->y_coeff[0]	+ wcs->y_coeff[1]*x +
	    wcs->y_coeff[2]*y	+ wcs->y_coeff[3]*x2 +
	    wcs->y_coeff[4]*y2	+ wcs->y_coeff[5]*xy;

	/* Derivative of Y model wrt x */
	gx = wcs->y_coeff[1]	+ wcs->y_coeff[3]*2.0*x +
	     wcs->y_coeff[5]*y;

	/* Derivative of Y model wrt y */
	gy = wcs->y_coeff[2]	+ wcs->y_coeff[4]*2.0*y +
	     wcs->y_coeff[5]*x;

	if (ncoeff2 > 6) {
	    g = g + wcs->y_coeff[6]*x3	+ wcs->y_coeff[7]*y3;
	    gx = gx + wcs->y_coeff[6]*3.0*x2;
	    gy = gy + wcs->y_coeff[7]*3.0*y2;
	    }

	if (ncoeff2 > 8) {
	    g = g +
		wcs->y_coeff[8]*x2y	+ wcs->y_coeff[9]*y2x +
		wcs->y_coeff[10]*r2	+ wcs->y_coeff[11]*x*r2 +
		wcs->y_coeff[12]*y*r2;

	    gx = gx +	wcs->y_coeff[8]*2.0*xy + 
			wcs->y_coeff[9]*y2 +
	 		wcs->y_coeff[10]*2.0*x +
			wcs->y_coeff[11]*(3.0*x2+y2) +
			wcs->y_coeff[12]*2.0*xy;

	    gy = gy +	wcs->y_coeff[8]*x2 +
			wcs->y_coeff[9]*2.0*xy +
			wcs->y_coeff[10]*2.0*y +
			wcs->y_coeff[11]*2.0*xy +
			wcs->y_coeff[12]*(3.0*y2+x2);
	    }

	f = f - xi;
	g = g - eta;
	dx = ((-f * gy) + (g * fy)) / ((fx * gy) - (fy * gx));
	dy = ((-g * fx) + (f * gx)) / ((fx * gy) - (fy * gx));
	x = x + dx;
	y = y + dy;
	if ((fabs(dx) < tolerance) && (fabs(dy) < tolerance)) break;
	}

    return 0;
}


/* Set polynomial fit coefficients in structure from arguments */

int
SetPoly (wcs, ncoeff1, ncoeff2, coeff)

struct WorldCoor *wcs;  /* World coordinate system structure */
int	ncoeff1;	/* Number of coefficients for x */
int	ncoeff2;	/* Number of coefficients for y */
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


/* Return polynomial fit coefficients from structure in arguments */

int
GetPoly (wcs, ncoeff1, ncoeff2, coeff)

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


/* Set FITS header polynomial fit coefficients from structure */

void
SetFITSPoly (header, wcs)

char    *header;        /* Image FITS header */
struct WorldCoor *wcs;  /* WCS structure */

{
    char keyword[16];
    int i;

    for (i = 0; i < wcs->ncoeff1; i++) {
	if (wcs->x_coeff[i] != 0.0) {
	    sprintf (keyword,"CO1_%d",i+1);
	    hputnr8 (header, keyword, -15, wcs->x_coeff[i]);
	    }
	}
    for (i = 0; i < wcs->ncoeff2; i++) {
	if (wcs->y_coeff[i] != 0.0) {
	    sprintf (keyword,"CO2_%d",i+1);
	    hputnr8 (header, keyword, -15, wcs->y_coeff[i]);
	    }
	}
    return;
}

/* May  4 1998	New subroutines for direct image pixel <-> sky polynomials
 */
