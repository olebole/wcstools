/* File wcstools/libwcs/wcspoly.c
 * February 24, 1999
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics

 * Module:	wcspoly.c
 * Purpose:	Polynomial correction to pixel coordinates
 * Subroutine:	polypos() converts from pixel location to corrected pixel
 * Subroutine:	polypix() converts from corrected pixel to pixel location   
 */

#include <string.h>
#include <stdio.h>
#include "wcs.h"
#include "fitshead.h"

#define WCS_NCOEFF	15

static double wcspolycomp();
static double wcspolydx();
static double wcspolydy();

int
wcspolypos (xpix, ypix, wcs, xpix1, ypix1)

/* Routine to determine accurate position for pixel coordinates */
/* returns 0 if successful otherwise 1 = angle too large for projection; */
/* based on amdpos() from getimage */

/* Input: */
double	xpix;		/* x pixel number  (RA or long without rotation) */
double	ypix;		/* y pixel number  (dec or lat without rotation) */
struct WorldCoor *wcs;	/* WCS parameter structure */

/* Output: */
double	*xpix1;		/* Corrected x pixel number */
double	*ypix1;		/* Corrected y pixel number */

{
    double x, y;

    /* Make pixels relative to reference pixels */
    x = xpix - wcs->crpix[0];
    y = ypix - wcs->crpix[1];

    /*  Compute new x and y from polynomial model */
    *xpix1 = wcspolycomp (x, y, wcs->x_coeff);
    *ypix1 = wcspolycomp (x,y, wcs->y_coeff);

    return (0);
}

int
wcspolypix (xpix1, ypix1, wcs, xpix, ypix)

/* Routine to determine pixel coordinates for sky position */
/* returns 0 if successful otherwise 1 = angle too large for projection; */
/* based on amdinv() from getimage */

/* Input: */
double	xpix1;		/* Right ascension or longitude in degrees */
double	ypix1;		/* Declination or latitude in degrees */
struct WorldCoor *wcs;	/* WCS parameter structure */

/* Output: */
double	*xpix;		/* x pixel number  (RA or long without rotation) */
double	*ypix;		/* y pixel number  (dec or lat without rotation) */

{
    double x, y, dx, dy;
    double f,fx,fy,g,gx,gy;
    double tolerance = 0.0000005;
    int    max_iterations = 50;
    int    i;

    /* Set initial value for x,y */
    if (wcs->x_coeff[1] == 0.0)
	x = xpix1 - wcs->x_coeff[0];
    else
	x = (xpix1 - wcs->x_coeff[0]) / wcs->x_coeff[1];
    if (wcs->y_coeff[2] == 0.0)
	y = ypix1 - wcs->y_coeff[0];
    else
	y = (ypix1 - wcs->y_coeff[0]) / wcs->y_coeff[2];

    /* Iterate by Newton's method */
    for (i = 0; i < max_iterations; i++) {

	f = wcspolycomp (x, y, wcs->x_coeff);

	/*  Derivative of X model wrt x */
	fx = wcspolydx (x, y, wcs->x_coeff);

	/* Derivative of X model wrt y */
	fy = wcspolydy (x, y, wcs->x_coeff);

	/* Y plate model */
	g = wcspolycomp (x, y, wcs->y_coeff);

	/* Derivative of Y model wrt x */
	gx = wcspolydx (x, y, wcs->y_coeff);

	/* Derivative of Y model wrt y */
	gy = wcspolydy (x, y, wcs->y_coeff);

	f = f - xpix1;
	g = g - ypix1;
	dx = ((-f * gy) + (g * fy)) / ((fx * gy) - (fy * gx));
	dy = ((-g * fx) + (f * gx)) / ((fx * gy) - (fy * gx));
	x = x + dx;
	y = y + dy;
	if ((fabs(dx) < tolerance) && (fabs(dy) < tolerance)) break;
	}

    /* Convert from plate pixels to image pixels */
    *xpix = x;
    *ypix = y;

    /* If position is off of the image, return offscale code */
    if (*xpix < 0.5 || *xpix > wcs->nxpix+0.5)
	return -1;
    if (*ypix < 0.5 || *ypix > wcs->nypix+0.5)
	return -1;

    return 0;
}


/* Set plate fit coefficients in structure from arguments */
int
wcspolyset (wcs, xcoeff, ycoeff, coeff)

struct WorldCoor *wcs;  /* World coordinate system structure */
int	*xcoeff;	/* 1 if corresponding x coefficient used */
int	*ycoeff;	/* 1 if corresponding y coefficient used */
double	*coeff;		/* Plate fit coefficients, only non-zero ones */

{
    int i, j;

    if (nowcs (wcs))
	return 1;

    j = 0;
    for (i = 0; i < WCS_NCOEFF; i++) {
	if (xcoeff[i])
	    wcs->x_coeff[i] = coeff[j++];
	else
	    wcs->x_coeff[i] = 0.0;
	}

    for (i = 0; i < WCS_NCOEFF; i++) {
	if (ycoeff[i])
	    wcs->y_coeff[i] = coeff[j++];
	else
	    wcs->y_coeff[i] = 0.0;
	}
    return 0;
}


/* Return plate fit coefficients from structure in argument */

int
wcspolyget (wcs, xcoeff, ycoeff, coeff)

struct WorldCoor *wcs;  /* World coordinate system structure */
int	*xcoeff;	/* 1 if corresponding x coefficient used */
int	*ycoeff;	/* 1 if corresponding y coefficient used */
double	*coeff;		/* Plate fit coefficients, only non-zero ones */

{
    int i, j;

    if (nowcs (wcs))
	return 1;

    j = 0;
    for (i = 0; i < WCS_NCOEFF; i++) {
	if (xcoeff[i])
	    coeff[i] = wcs->x_coeff[j++];
	}

    for (i = 0; i < WCS_NCOEFF; i++) {
	if (ycoeff[i])
	    coeff[i] = wcs->y_coeff[j++];
	}

    return 0;
}


/* Set FITS header plate fit coefficients from structure */
void
wcspolyfset (header, wcs)

char    *header;        /* Image FITS header */
struct WorldCoor *wcs;  /* WCS structure */

{
    char keyword[16];
    int i;

    for (i = 0; i < WCS_NCOEFF; i++) {
	if (wcs->x_coeff[i] != 0.0) {
	    sprintf (keyword,"PX1_%d",i+1);
	    hputnr8 (header, keyword, -15, wcs->x_coeff[i]);
	    }
	}
    for (i = 0; i < WCS_NCOEFF; i++) {
	if (wcs->y_coeff[i] != 0.0) {
	    sprintf (keyword,"PX2_%d",i+1);
	    hputnr8 (header, keyword, -15, wcs->y_coeff[i]);
	    }
	}
    return;
}


/* Read FITS header plate fit coefficients into structure */
void
wcspolyfget (header, wcs)

char    *header;        /* Image FITS header */
struct WorldCoor *wcs;  /* WCS structure */

{
    char keyword[16];
    int i;

    for (i = 0; i < WCS_NCOEFF; i++) {
	sprintf (keyword,"PX1_%d",i+1);
	if (!hgetr8 (header, keyword, wcs->x_coeff[i]))
	    wcs->x_coeff[i] = 0.0;
	}
    for (i = 0; i < WCS_NCOEFF; i++) {
	sprintf (keyword,"PX2_%d",i+1);
	if (!hgetr8 (header, keyword, wcs->y_coeff[i]))
	    wcs->y_coeff[i] = 0.0;
	}
    return;
}

/* WCSPOLYCOMP: Compute polynomial with cross terms */

static double wx0 = 0.0;
static double wy0 = 0.0;
static double x2 = 0.0;
static double y2 = 0.0;
static double xy = 0.0;
static double x3 = 0.0;
static double y3 = 0.0;
static double x2y = 0.0;
static double xy2 = 0.0;
static double r2 = 0.0;
static double xr2 = 0.0;
static double yr2 = 0.0;
static double r4 = 0.0;
static double xr4 = 0.0;
static double yr4 = 0.0;

static double
wcspolycomp (x, y, coeff)

double	x, y;		/* Coordinates to be used in transform */
double	*coeff;		/* Polynomial coefficients */
{
    double v;

    /* If new x or y, compute new terms */
    if (x != wx0 || y != wy0) {
	wx0 = x;
	wy0 = y;
	x2 = x * x;
	y2 = y * y;
	xy = x * y;
	x3 = x * x2;
	y3 = y * y2;
	x2y = x2 * y;
	xy2 = x * y2;
	r2 = x2 + y2;
	r4 = r2 * r2;
	xr4 = x * r4;
	yr4 = y * r4;
	}

    v = coeff[ 0]	+ coeff[ 1]*x	+ coeff[ 2]*y	+ coeff[ 3]*x2 +
	coeff[ 4]*y2	+ coeff[ 5]*xy	+ coeff[ 6]*x3	+ coeff[ 7]*y3 +
	coeff[ 8]*x2y	+ coeff[ 9]*xy2	+ coeff[10]*r2	+ coeff[11]*xr2 +
	coeff[12]*yr2	+ coeff[13]*r4	+ coeff[14]*xr4	+ coeff[15]*yr4;

    return (v);
}


static double dx0 = 0.0;
static double dy0 = 0.0;
static double tx = 0.0;
static double ty = 0.0;
static double tx2 = 0.0;
static double ty2 = 0.0;
static double txy = 0.0;
static double fr2 = 0.0;
static double fx3 = 0.0;
static double fy3 = 0.0;
static double x4 = 0.0;
static double y4 = 0.0;

/*  Compute derivative of one coordinate with respect to x */

static double
wcspolydx (x, y, coeff)

double	x, y;		/* Coordinates to be used in transform */
double	*coeff;		/* Polynomial coefficients */
{
    double dx;

    /* If new x or y, compute new terms */
    if (x != dx0 || y != dy0) {
	dx0 = x;
	dy0 = y;
	tx = 2.0 * x;
	tx2 = 3.0 * x2;
	ty = 2.0 * y;
	ty2 = 3.0 * y2;
	txy = 2.0 * xy;
	fr2 = 4.0 * r2;
	fx3 = 4.0 * x3;
	fy3 = 4.0 * y3;
	x4 = x2 * x2;
	y4 = y2 * y2;
	}

    dx = coeff[1]	+ coeff[3]*tx	+ coeff[5]*y	+ coeff[6]*tx2 +
	 coeff[8]*txy	+ coeff[9]*y2	+ coeff[10]*tx	+ coeff[11]*(tx2+y2) +
	 coeff[12]*txy	+ coeff[13]*fr2*wx0 + coeff[14]*(5.0*x4 + 6.0*x2*y2) +
	 coeff[15]*(fy3*wx0 + fx3);

    return (dx);
}


/*  Compute derivative of one coordinate with respect to y */
static double
wcspolydy (x, y, coeff)

double	x, y;		/* Coordinates to be used in transform */
double	*coeff;		/* Polynomial coefficients */
{
    double dy;

    /* If new x or y, compute new terms */
    if (x != dx0 || y != dy0) {
	dx0 = x;
	dy0 = y;
	tx = 2.0 * x;
	tx2 = 3.0 * x2;
	ty = 2.0 * y;
	ty2 = 3.0 * y2;
	txy = 2.0 * xy;
	fr2 = 4.0 * r2;
	fx3 = 4.0 * x3;
	fy3 = 4.0 * y3;
	x4 = x2 * x2;
	y4 = y2 * y2;
	}

    dy = coeff[2]	+ coeff[4]*ty	+ coeff[5]*x	+ coeff[7]*ty2 +
	 coeff[8]*x2	+ coeff[9]*txy	+ coeff[10]*ty	+ coeff[11]*txy +
	 coeff[12]*(ty2+x2) + coeff[13]*fr2*wy0 + coeff[14]*(fx3*wy0 + fy3) +
	 coeff[15]*(5.0*y4 + 6.0*x2*y2);

    return (dy);
}

/* Feb 24 1999	New subroutines for direct image pixel <-> corrected image pixel
 */
