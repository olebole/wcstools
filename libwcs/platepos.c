/* File saoimage/wcslib/platepos.c
 * February 25, 1996
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics

 * Module:	platepos.c (Plate solution WCS conversion
 * Purpose:	Compute WCS from Digital Sky Survey plate fit
 * Subroutine:	platepos() converts from pixel location to RA,Dec 
 * Subroutine:	platepix() converts from RA,Dec to pixel location   

    These functions are based on the astrmcal.c portion of GETIMAGE by
    J. Doggett and the documentation distributed with the Digital Sky Survey.

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
  double x, y, xmm, ymm, xmm2, ymm2, xmm3, ymm3, x2y2;
  double xi, xir, eta, etar, raoff, ra, dec;
  double cond2r = 1.745329252e-2;
  double cons2r = 206264.8062470964;
  double twopi = 6.28318530717959;
  double ctan, ccos;

/*  Ignore magnitude and color terms 
  double mag = 0.0;
  double color = 0.0; */

/* Convert from image pixels to plate pixels */
  x = xpix + wcs->x_pixel_offset - 1.0 + 0.5;
  y = ypix + wcs->y_pixel_offset - 1.0 + 0.5;

/* Convert from pixels to millimeters */
  xmm = (wcs->ppo_coeff[2] - x * wcs->x_pixel_size) / 1000.0;
  ymm = (y * wcs->y_pixel_size - wcs->ppo_coeff[5]) / 1000.0;
  xmm2 = xmm * xmm;
  ymm2 = ymm * ymm;
  xmm3 = xmm * xmm2;
  ymm3 = ymm * ymm2;
  x2y2 = xmm2 + ymm2;

/*  Compute coordinates from x,y and plate model */

  xi =  wcs->amd_x_coeff[ 0]*xmm	+ wcs->amd_x_coeff[ 1]*ymm +
	wcs->amd_x_coeff[ 2]		+ wcs->amd_x_coeff[ 3]*xmm2 +
	wcs->amd_x_coeff[ 4]*xmm*ymm	+ wcs->amd_x_coeff[ 5]*ymm2 +
	wcs->amd_x_coeff[ 6]*(x2y2)	+ wcs->amd_x_coeff[ 7]*xmm3 +
	wcs->amd_x_coeff[ 8]*xmm2*ymm	+ wcs->amd_x_coeff[ 9]*xmm*ymm2 +
	wcs->amd_x_coeff[10]*ymm3	+ wcs->amd_x_coeff[11]*xmm*(x2y2) +
	wcs->amd_x_coeff[12]*xmm*x2y2*x2y2;

/*  Ignore magnitude and color terms 
	+ wcs->amd_x_coeff[13]*mag	+ wcs->amd_x_coeff[14]*mag*mag +
	wcs->amd_x_coeff[15]*mag*mag*mag + wcs->amd_x_coeff[16]*mag*xmm +
	wcs->amd_x_coeff[17]*mag*x2y2	+ wcs->amd_x_coeff[18]*mag*xmm*x2y2 +
	wcs->amd_x_coeff[19]*color; */

  eta =	wcs->amd_y_coeff[ 0]*ymm	+ wcs->amd_y_coeff[ 1]*xmm +
	wcs->amd_y_coeff[ 2]		+ wcs->amd_y_coeff[ 3]*ymm2 +
	wcs->amd_y_coeff[ 4]*xmm*ymm	+ wcs->amd_y_coeff[ 5]*xmm2 +
	wcs->amd_y_coeff[ 6]*(x2y2)	+ wcs->amd_y_coeff[ 7]*ymm3 +
	wcs->amd_y_coeff[ 8]*ymm2*xmm	+ wcs->amd_y_coeff[ 9]*ymm*xmm2 +
	wcs->amd_y_coeff[10]*xmm3	+ wcs->amd_y_coeff[11]*ymm*(x2y2) +
	wcs->amd_y_coeff[12]*ymm*x2y2*x2y2;

/*  Ignore magnitude and color terms 
	+ wcs->amd_y_coeff[13]*mag	+ wcs->amd_y_coeff[14]*mag*mag +
	wcs->amd_y_coeff[15]*mag*mag*mag + wcs->amd_y_coeff[16]*mag*ymm +
	wcs->amd_y_coeff[17]*mag*x2y2)	+ wcs->amd_y_coeff[18]*mag*ymm*x2y2 +
	wcs->amd_y_coeff[19]*color; */

/* Convert to radians */

  xir = xi / cons2r;
  etar = eta / cons2r;

/* Convert to RA and Dec */

  ctan = tan (wcs->plate_dec);
  ccos = cos (wcs->plate_dec);
  raoff = atan2 (xir / ccos, 1.0 - etar * ctan);
  ra = raoff + wcs->plate_ra;
  if (ra < 0.0) ra = ra + twopi;
  *xpos = ra / cond2r;

  dec = atan (cos (raoff) / ((1.0 - (etar * ctan)) / (etar + ctan)));
  *ypos = dec / cond2r;
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
  double div,xi,eta,x,y,xy,x2,y2,x2y,y2x,x3,y3,x4,y4,x2y2,cjunk,dx,dy;
  double sypos,cypos,syplate,cyplate,sxdiff,cxdiff;
  double f,fx,fy,g,gx,gy, xmm, ymm;
  double conr2s = 206264.8062470964;
  double tolerance = 0.0000005;
  int    max_iterations = 50;
  int    i;
  double xr, yr; 	/* position in radians */
  double cond2r = 1.745329252e-2;

/* Convert RA and Dec in radians to standard coordinates on a plate */
  xr = xpos * cond2r;
  yr = ypos * cond2r;
  sypos = sin (yr);
  cypos = cos (yr);
  syplate = sin (wcs->plate_dec);
  cyplate = cos (wcs->plate_dec);
  sxdiff = sin (xr - wcs->plate_ra);
  cxdiff = cos (xr - wcs->plate_ra);
  div = (sypos * syplate) + (cypos * cyplate * cxdiff);
  xi = cypos * sxdiff * conr2s / div;
  eta = ((sypos * cyplate) - (cypos * syplate * cxdiff)) * conr2s / div;

/* Set initial value for x,y */
  xmm = xi / wcs->plate_scale;
  ymm = eta / wcs->plate_scale;

/* Iterate by Newton's method */
  for (i = 0; i < max_iterations; i++) {

    /* X plate model */
    xy = xmm * ymm;
    x2 = xmm * xmm;
    y2 = ymm * ymm;
    x2y = x2 * ymm;
    y2x = y2 * xmm;
    x2y2 = x2 + y2;
    cjunk = x2y2 * x2y2;
    x3 = x2 * xmm;
    y3 = y2 * ymm;
    x4 = x2 * x2;
    y4 = y2 * y2;
    f = wcs->amd_x_coeff[0]*xmm      + wcs->amd_x_coeff[1]*ymm +
        wcs->amd_x_coeff[2]          + wcs->amd_x_coeff[3]*x2 +
        wcs->amd_x_coeff[4]*xy       + wcs->amd_x_coeff[5]*y2 +
        wcs->amd_x_coeff[6]*x2y2     + wcs->amd_x_coeff[7]*x3 +
        wcs->amd_x_coeff[8]*x2y      + wcs->amd_x_coeff[9]*y2x +
        wcs->amd_x_coeff[10]*y3      + wcs->amd_x_coeff[11]*xmm*x2y2 +
        wcs->amd_x_coeff[12]*xmm*cjunk;
    /* magnitude and color terms ignored
      + wcs->amd_x_coeff[13]*mag +
        wcs->amd_x_coeff[14]*mag*mag   + wcs->amd_x_coeff[15]*mag*mag*mag +
        wcs->amd_x_coeff[16]*mag*xmm   + wcs->amd_x_coeff[17]*mag*(x2+y2) +
        wcs->amd_x_coeff[18]*mag*xmm*(x2+y2)  + wcs->amd_x_coeff[19]*color;
    */

    /*  Derivative of X model wrt x */
    fx = wcs->amd_x_coeff[0]           + wcs->amd_x_coeff[3]*2.0*xmm +
         wcs->amd_x_coeff[4]*ymm       + wcs->amd_x_coeff[6]*2.0*xmm +
         wcs->amd_x_coeff[7]*3.0*x2    + wcs->amd_x_coeff[8]*2.0*xy +
         wcs->amd_x_coeff[9]*y2        + wcs->amd_x_coeff[11]*(3.0*x2+y2) +
         wcs->amd_x_coeff[12]*(5.0*x4 +6.0*x2*y2+y4);
    /* magnitude and color terms ignored
         wcs->amd_x_coeff[16]*mag      + wcs->amd_x_coeff[17]*mag*2.0*xmm +
         wcs->amd_x_coeff[18]*mag*(3.0*x2+y2);
    */

    /* Derivative of X model wrt y */
    fy = wcs->amd_x_coeff[1]           + wcs->amd_x_coeff[4]*xmm +
         wcs->amd_x_coeff[5]*2.0*ymm   + wcs->amd_x_coeff[6]*2.0*ymm +
         wcs->amd_x_coeff[8]*x2        + wcs->amd_x_coeff[9]*2.0*xy +
         wcs->amd_x_coeff[10]*3.0*y2   + wcs->amd_x_coeff[11]*2.0*xy +
         wcs->amd_x_coeff[12]*4.0*xy*x2y2;
    /* magnitude and color terms ignored
         wcs->amd_x_coeff[17]*mag*2.0*ymm +
         wcs->amd_x_coeff[18]*mag*2.0*xy;
    */

    /* Y plate model */
    g = wcs->amd_y_coeff[0]*ymm       + wcs->amd_y_coeff[1]*xmm +
       wcs->amd_y_coeff[2]            + wcs->amd_y_coeff[3]*y2 +
       wcs->amd_y_coeff[4]*xy         + wcs->amd_y_coeff[5]*x2 +
       wcs->amd_y_coeff[6]*x2y2       + wcs->amd_y_coeff[7]*y3 +
       wcs->amd_y_coeff[8]*y2x        + wcs->amd_y_coeff[9]*x2y +
       wcs->amd_y_coeff[10]*x3        + wcs->amd_y_coeff[11]*ymm*x2y2 +
       wcs->amd_y_coeff[12]*ymm*cjunk;
    /* magnitude and color terms ignored
       wcs->amd_y_coeff[13]*mag        + wcs->amd_y_coeff[14]*mag*mag +
       wcs->amd_y_coeff[15]*mag*mag*mag + wcs->amd_y_coeff[16]*mag*ymm +
       wcs->amd_y_coeff[17]*mag*x2y2 +
       wcs->amd_y_coeff[18]*mag*ymm*x2y2 + wcs->amd_y_coeff[19]*color;
    */

    /* Derivative of Y model wrt x */
    gx = wcs->amd_y_coeff[1]           + wcs->amd_y_coeff[4]*ymm +
         wcs->amd_y_coeff[5]*2.0*xmm   + wcs->amd_y_coeff[6]*2.0*xmm +
         wcs->amd_y_coeff[8]*y2       + wcs->amd_y_coeff[9]*2.0*xy +
         wcs->amd_y_coeff[10]*3.0*x2  + wcs->amd_y_coeff[11]*2.0*xy +
         wcs->amd_y_coeff[12]*4.0*xy*x2y2;
    /* magnitude and color terms ignored
         wcs->amd_y_coeff[17]*mag*2.0*xmm +
         wcs->amd_y_coeff[18]*mag*ymm*2.0*xmm;
    */

    /* Derivative of Y model wrt y */
    gy = wcs->amd_y_coeff[0]            + wcs->amd_y_coeff[3]*2.0*ymm +
         wcs->amd_y_coeff[4]*xmm        + wcs->amd_y_coeff[6]*2.0*ymm +
         wcs->amd_y_coeff[7]*3.0*y2     + wcs->amd_y_coeff[8]*2.0*xy +
         wcs->amd_y_coeff[9]*x2         + wcs->amd_y_coeff[11]*(x2+3.0*y2) +
         wcs->amd_y_coeff[12]*(5.0*y4 + 6.0*x2*y2 + x4);
    /* magnitude and color terms ignored
         wcs->amd_y_coeff[16]*mag       + wcs->amd_y_coeff[17]*mag*2.0*ymm +
         wcs->amd_y_coeff[18]*mag*(x2+3.0*y2);
    */

    f = f - xi;
    g = g - eta;
    dx = ((-f * gy) + (g * fy)) / ((fx * gy) - (fy * gx));
    dy = ((-g * fx) + (f * gx)) / ((fx * gy) - (fy * gx));
    xmm = xmm + dx;
    ymm = ymm + dy;
    if ((fabs(dx) < tolerance) && (fabs(dy) < tolerance)) break;
    }

/* Convert mm from plate center to plate pixels */
  x = (wcs->ppo_coeff[2] - xmm*1000.0) / wcs->x_pixel_size;
  y = (wcs->ppo_coeff[5] + ymm*1000.0) / wcs->y_pixel_size;

/* Convert from plate pixels to image pixels */
  *xpix = x - wcs->x_pixel_offset + 1.0 - 0.5;
  *ypix = y - wcs->y_pixel_offset + 1.0 - 0.5;

/* If position is off of the image, return offscale code */
  if (*xpix < 0.5 || *xpix > wcs->nxpix+0.5)
    return -1;
  if (*ypix < 0.5 || *ypix > wcs->nypix+0.5)
    return -1;

  return 0;
}
/* Mar  6 1995	Original version of this code
   May  4 1995	Fix eta cross terms which were all in y
   Jun 21 1995	Add inverse routine
   Oct 17 1995	Fix inverse routine (degrees -> radians)
   Nov  7 1995	Add half pixel to image coordinates to get astrometric
                  plate coordinates
   Feb 26 1996	Fix plate to image pixel conversion error
 */
