/*  worldpos.c -- WCS Algorithms from Classic AIPS.
 *  October 8, 1997
 *  By Doug Mink <dmink@cfa.harvard.edu> after Bill Cotton and Don Wells

 * Module:	worldpos.c
 * Purpose:	Perform forward and reverse WCS computations for WCSLIB projections
 * Subroutine:	worldpos() converts from pixel location to RA,Dec 
 * Subroutine:	worldpix() converts from RA,Dec         to pixel location   

    These two ANSI C functions, worldpos() and worldpix(), perform forward and
    reverse WCS computations for 25 types of projective geometries:

	worldpos() converts from pixel location to RA,Dec 
	worldpix() converts from RA,Dec         to pixel location   

    where "(RA,Dec)" are more generically (long,lat). These functions
    are based on the draft WCS proposal by Greisen and Calabretta at:
    
	ftp://fits.cv.nrao.edu/fits/documents/wcs/wcs.all.ps.Z
    
    The original version of this code was Emailed to D.Wells on
    Friday, 23 September by Bill Cotton <bcotton@gorilla.cv.nrao.edu>,
    who described it as a "..more or less.. exact translation from the
    AIPSish..". Changes were made by Don Wells <dwells@nrao.edu>
    during the period October 11-13, 1994:
    1) added GNU license and header comments
    2) added testpos.c program to perform extensive circularity tests
    3) changed float-->double to get more than 7 significant figures
    4) testpos.c circularity test failed on MER and AIT. B.Cotton
       found that "..there were a couple of lines of code [in] the wrong
       place as a result of merging several Fortran routines." 
    5) testpos.c found 0h wraparound in worldpix() and worldpos().
    6) E.Greisen recommended removal of various redundant if-statements,
       and addition of a 360d difference test to MER case of worldpos(). 
    By October 19, 1994, Doug Mink changed the input to use his wcs data
    structure.  In December 1994, he implemented the CD rotation matrix,
    and over the next three years, various small changes were made
    In October 1997 replaced the projection code with calls to Mark
    Calabretta's WCSLIB subroutines.  Almost all of the code which is
    left was written by Doug Mink, so the old NRAO copyright was removed.
*/
#include <math.h>
#include <string.h>
#include <stdio.h>
#include "wcs.h"

int
worldpos (xpix, ypix, wcs, xpos, ypos)

/* Routine to determine accurate position for pixel coordinates */
/* returns 0 if successful otherwise 1 = angle too large for projection; */

/* Input: */
double	xpix;		/* x pixel number  (RA or long without rotation) */
double	ypix;		/* y pixel number  (dec or lat without rotation) */
struct WorldCoor *wcs;	/* WCS parameter structure */

/* Output: */
double	*xpos;		/* x (RA) coordinate (deg) */
double	*ypos;		/* y (dec) coordinate (deg) */

{
  double cosr, sinr, dx, dy, temp, x, y, phi, theta;

  /* Structure elements */
  double xref;		/* X reference coordinate value (deg) */
  double yref;		/* Y reference coordinate value (deg) */
  double xrefpix;	/* X reference pixel */
  double yrefpix;	/* Y reference pixel */
  double xinc;		/* X coordinate increment (deg) */
  double yinc;		/* Y coordinate increment (deg) */
  double rot;		/* Optical axis rotation (deg)  (N through E) */
  int itype = wcs->pcode;

/* Set local projection parameters */
  xref = wcs->xref;
  yref = wcs->yref;
  xrefpix = wcs->xrefpix;
  yrefpix = wcs->yrefpix;
  xinc = wcs->xinc;
  yinc = wcs->yinc;

/* Offset from ref pixel */
  dx = xpix - xrefpix;
  dy = ypix - yrefpix;

/* Scale and rotate using CD matrix */
  if (wcs->rotmat) {
    temp = dx * wcs->cd11 + dy * wcs->cd12;
    dy = dx * wcs->cd21 + dy * wcs->cd22;
    dx = temp;
    }
  else {

/* Check axis increments - bail out if either 0 */
    if ((xinc==0.0) || (yinc==0.0)) {
      *xpos=0.0;
      *ypos=0.0;
      return 2;
      }

/* Scale using CDELT */
    dx = dx * xinc;
    dy = dy * yinc;

/* Take out rotation from CROTA */
    rot = wcs->rot;
    cosr = wcs->crot;
    sinr = wcs->srot;
    if (rot != 0.0) {
      temp = dx * cosr - dy * sinr;
      dy = dy * cosr + dx * sinr;
      dx = temp;
      }
    }

/* Flip axes if required */
  if (wcs->coorflip) {
    temp = dx;
    dx = dy;
    dy = temp;
    }
  x = xref + dx;
  y = yref + dx;

/* Default, linear result for error or pixel return  */
  *xpos = x;
  *ypos = y;
  if (itype < 0)
    return 0;

/* PIXEL or LINEAR */
  if (itype < 1) {
    *xpos =  x;
    if (dx > 180.0)
      x = x - 360.0;
    if (dx < -180.0)
      x = x + 360.0;
    if (*xpos < 0.0)
      x += 360.0; /* added by DCW 10/12/94 */
    *xpos = x;
    *ypos = y;
    }

/* WCSLIB projections */ 
  else {
    if (celrev (wcs->ptype, x, y, wcs->prj, &phi, &theta,
		wcs->cel, xpos, ypos))
      return 1;
    }

  return 0;
}  /* End of worldpos */


int
worldpix (xpos, ypos, wcs, xpix, ypix)

/*-----------------------------------------------------------------------*/
/* routine to determine accurate pixel coordinates for an RA and Dec     */
/* returns 0 if successful otherwise:                                    */
/*  1 = angle too large for projection;                                  */
/*  2 = bad values                                                       */
/* does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections            */
/* anything else is linear                                               */

/* Input: */
double	xpos;		/* x (RA) coordinate (deg) */
double	ypos;		/* y (dec) coordinate (deg) */
struct WorldCoor *wcs;	/* WCS parameter structure */

/* Output: */
double	*xpix;		/* x pixel number  (RA or long without rotation) */
double	*ypix;		/* y pixel number  (dec or lat without rotation) */
{
  double dx, dy, ra0, dec0, ra, dec, dt, sint;
  double sinr, cosr, temp, phi, theta;

/* Structure elements */
  double xref;		/* x reference coordinate value (deg) */
  double yref;		/* y reference coordinate value (deg) */
  double xrefpix;	/* x reference pixel */
  double yrefpix;	/* y reference pixel */
  double xinc;		/* x coordinate increment (deg) */
  double yinc;		/* y coordinate increment (deg) */
  double rot;		/* rotation (deg)  (from N through E) */
  int itype;

/* Set local projection parameters */
  xref = wcs->xref;
  yref = wcs->yref;
  xrefpix = wcs->xrefpix;
  yrefpix = wcs->yrefpix;
  xinc = wcs->xinc;
  yinc = wcs->yinc;
  rot = wcs->rot;
  cosr = wcs->crot;
  sinr = wcs->srot;

/* Projection type */
  itype = wcs->pcode;

  /* 0h wrap-around tests added by D.Wells 10/12/94: */
  dt = (xpos - xref);
  if (itype >= 0) {
    if (dt > 180.0) xpos -= 360.0;
    if (dt < -180.0) xpos += 360.0;
    /* NOTE: changing input argument xpos is OK (call-by-value in C!) */
    }

/* Nonlinear position */
  if (itype > 0) {
    if (wcs->coorflip) {
      dec0 = degrad (xref);
      ra0 = degrad (yref);
      }
    else {
      ra0 = degrad (xref);
      dec0 = degrad (yref);
      }
    ra = degrad (xpos);
    dec = degrad (ypos);
    }

/* For linear or pixel projection */
  if (itype < 1) {
    dx = xpos - xref;
    dy = ypos - yref;
    }

/* Process WCSLIB projections */
  else {
    if (celfwd (wcs->ptype, xpos, ypos, wcs->prj, &phi, &theta,
	        wcs->cel, &dx, &dy))
      return 1;
    }

/* Scale and rotate using CD matrix */
  if (wcs->rotmat) {
    temp = dx * wcs->dc11 + dy * wcs->dc12;
    dy = dx * wcs->dc21 + dy * wcs->dc22;
    dx = temp;
    }
  else {

/* Correct for rotation */
    if (rot!=0.0) {
      temp = dx*cosr + dy*sinr;
      dy = dy*cosr - dx*sinr;
      dx = temp;
      }

/* Scale using CDELT */
    if (xinc != 0.)
      dx = dx / xinc;
    if (yinc != 0.)
      dy = dy / yinc;
    }

  if (wcs->coorflip) {
    temp = dx;
    dx = dy;
    dy = temp;
    }

/* Convert to pixels  */
  *xpix = dx + xrefpix;
  *ypix = dy + yrefpix;

  return 0;
}  /* end worldpix */

/* Oct 26 1995	Fix bug which interchanged RA and Dec twice when coorflip
 * Oct 31 1996	Fix CD matrix use in WORLDPIX
 * Nov  4 1996	Eliminate extra code for linear projection in WORLDPIX
 * Nov  5 1996	Add coordinate flip in WORLDPIX
 *
 * May 22 1997	Avoid angle wraparound when CTYPE is pixel
 * Jun  4 1997	Return without angle conversion from worldpos if type is PIXEL
 * Oct  8 1997	Use Mark Calabretta's WCSLIB for coordinate mapping
 */
