/*** File libwcs/distort.c
 *** April 2, 2003
 *** By Doug Mink, dmink@cfa.harvard.edu, 
 *** Based on code written by Jing Li, IPAC
 *** Harvard-Smithsonian Center for Astrophysics
 *** Copyright (C) 2003
 *** Smithsonian Astrophysical Observatory, Cambridge, MA, USA

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Correspondence concerning WCSTools should be addressed as follows:
           Internet email: dmink@cfa.harvard.edu
           Postal address: Doug Mink
                           Smithsonian Astrophysical Observatory
                           60 Garden St.
                           Cambridge, MA 02138 USA

 * Module:	distort.c (World Coordinate Systems)
 * Purpose:	Convert focal plane coordinates to pixels and vice versa:
 * Subroutine:	pix2foc (wcs, x, y, u, v) pixel coordinates -> focal plane coordinates
 * Subroutine:	foc2pix (wcs, u, v, x, y) focal plane coordinates -> pixel coordinates
 */

#include "wcs.h"

void
foc2pix (wcs, x, y, u, v)

struct WorldCoor *wcs;  /* World coordinate system structure */
double	x, y;		/* Focal plane coordinates */
double	*u, *v;		/* Image pixel coordinates (returned) */
{
    int m, n, i, j, k;
    double s[DISTMAX], sum;
    double temp_x, temp_y;

    /* SIRTF distortion */
    if (wcs->distcode == DISTORT_SIRTF) {
	m = wcs->distort.ap_order;
	n = wcs->distort.bp_order;

	temp_x = x - wcs->xrefpix;
	temp_y = y - wcs->yrefpix;

	/* compute u */
	for (j = 0; j <= m; j++) {
	    s[j] = wcs->distort.ap[m-j][j];
	    for (k = j-1; k >= 0; k--) {
	   	s[j] = (temp_y * s[j]) + wcs->distort.ap[m-j][k];
		}
	    }
  
	sum = s[0];
	for (i=m; i>=1; i--){
	    sum = (temp_x * sum) + s[m-i+1];
	    }
	*u = sum;

	/* compute v*/
	for (j = 0; j <= n; j++) {
	    s[j] = wcs->distort.bp[n-j][j];
	    for (k = j-1; k >= 0; k--) {
		s[j] = temp_y*s[j] + wcs->distort.bp[n-j][k];
		}
	    }
   
	sum = s[0];
	for (i = n; i >= 1; i--)
	    sum = temp_x * sum + s[n-i+1];

	*v = sum;

	*u = x + *u;
	*v = y + *v;
	}

    /* If no distortion, return pixel positions unchanged */
    else {
	*u = x;
	*v = y;
	}

    return;
}


void
pix2foc (wcs, u, v, x, y)

struct WorldCoor *wcs;  /* World coordinate system structure */
double u, v;		/* Image pixel coordinates */
double *x, *y;		/* Focal plane coordinates (returned) */
{
    int m, n, i, j, k;
    double s[DISTMAX], sum;
    double temp_u, temp_v;

    /* SIRTF distortion */
    if (wcs->distcode == DISTORT_SIRTF) {
	m = wcs->distort.a_order;
	n = wcs->distort.b_order;

	temp_u = u - wcs->xrefpix;
	temp_v = v - wcs->yrefpix;

	/* compute u */
	for (j = 0; j <= m; j++) {
	    s[j] = wcs->distort.a[m-j][j];
	    for (k = j-1; k >= 0; k--) {
		s[j] = (temp_v * s[j]) + wcs->distort.a[m-j][k];
		}
	    }
  
	sum = s[0];
	for (i=m; i>=1; i--){
	    sum = temp_u*sum + s[m-i+1];
	    }
	*x = sum;

	/* compute v*/
	for (j=0; j<=n; j++) {
	    s[j] = wcs->distort.b[n-j][j];
	    for (k=j-1; k>=0; k--) {
		s[j] =temp_v*s[j] + wcs->distort.b[n-j][k];
		}
	    }
   
	sum = s[0];
	for (i=n; i>=1; i--)
	    sum = temp_u*sum + s[n-i+1];

	*y = sum;
  
	*x = u + *x;
	*y = v + *y;

/*	*x = u + *x + coeff.crpix1; */
/*	*y = v + *y + coeff.crpix2; */
	}

    /* If no distortion, return pixel positions unchanged */
    else {
	*x = u;
	*y = v;
	}

    return;
}
/* Apr  2 2003	New subroutines
 */
