/* File libwcs/imrotate.c
 * February 23, 1998
 * By Doug Mink
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitshead.h"

/* Rotate an image by 90, 180, or 270 degrees, with an optional
 * reflection across the vertical axis.
 * verbose generates extra info on stdout.
 * return 0 if successful or -1.
 */

int
RotFITS (pathname, header, image0, rotate, mirror, bitpix2, verbose)

char	*pathname;	/* Name of file which is being changed */
char	*header;	/* FITS header */
char	**image0;	/* Image pixels */
int	rotate;		/* Angle to by which to rotate image (90, 180, 270) */
int	mirror;		/* 1 to reflect image around vertical axis */
int	bitpix2;	/* Number of bits per pixel in output image */
int	verbose;

{
    int bitpix1, ny, nx, nax;
    int x1, y1, x2, y2, nbytes;
    char *rotimage;
    char history[72];
    char *filename;
    char *image;		/* Image pixels */
    extern int DelWCSFITS();

    image = *image0;

    if (rotate == 1)
	rotate = 90;
    else if (rotate == 2)
	rotate = 180;
    else if (rotate == 3)
	rotate = 270;
    else if (rotate < 0)
	rotate = rotate + 360;

    filename = strrchr (pathname,'/');
    if (filename)
	filename = filename + 1;
    else
	filename = pathname;

    /* Get image size */
    nax = 0;
    if (hgeti4 (header,"NAXIS",&nax) < 1)
	return (-1);
    else {
	if (hgeti4 (header,"NAXIS1",&nx) < 1)
	    return (-1);
	else {
	    if (hgeti4 (header,"NAXIS2",&ny) < 1)
		return (-1);
	    }
	}
    bitpix1 = 16;
    hgeti4 (header,"BITPIX", &bitpix1);
    if (bitpix2 == 0)
	bitpix2 = bitpix1;

    /* Delete WCS fields in header */
    if (rotate != 0 || mirror)
	(void) DelWCSFITS (header, verbose);

    /* Allocate buffer for rotated image */
    switch (bitpix2) {
	case 16:
	    nbytes = nx * ny * 2;
	    break;
	case 32:
	    nbytes = nx * ny * 4;
	    break;
	case -16:
	    nbytes = nx * ny * 2;
	    break;
	case -32:
	    nbytes = nx * ny * 4;
	    break;
	case -64:
	    nbytes = nx * ny * 8;
	    break;
	default:
	    return (-1);
	}

    rotimage = (char *) malloc (nbytes);
    if (rotimage == NULL)
	return (-1);

    if (bitpix1 != bitpix2) {
	sprintf (history,"Copy of image %s bits per pixel %d -> %d",
		filename, bitpix1, bitpix2);
	hputc (header,"HISTORY",history);
	if (verbose)
	    printf ("%s\n",history);
	}

    /* Mirror image without rotation */
    if (rotate < 45.0 && rotate > -45.0) {
	if (mirror) {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = nx - x1 - 1;
		    y2 = y1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,nx,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s reflected",filename);
	    hputc (header,"HISTORY",history);
            hputc (header,"HISTORY",history);
	    }
	else {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = x1;
		    y2 = y1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,nx,x2,y2);
		    }
		}
	    }
	}

    /* Rotate by 90 degrees */
    else if (rotate >= 45 && rotate < 135) {
	if (mirror) {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = ny - y1 - 1;
		    y2 = nx - x1 - 1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,ny,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s reflected and rotated 90 degrees",filename);
            hputc (header,"HISTORY",history);
	    }
	else {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = ny - y1 - 1;
		    y2 = x1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,ny,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s rotated 90 degrees",filename);
            hputc (header,"HISTORY",history);
	    }
	hputi4 (header,"NAXIS1",ny);
	hputi4 (header,"NAXIS2",nx);
	}

    /* Rotate by 180 degrees */
    else if (rotate >= 135 && rotate < 225) {
	if (mirror) {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = x1;
		    y2 = ny - y1 - 1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,nx,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s reflected and rotated 180 degrees",filename);
            hputc (header,"HISTORY",history);
	    }
	else {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = nx - x1 - 1;
		    y2 = ny - y1 - 1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,nx,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s rotated 180 degrees",filename);
            hputc (header,"HISTORY",history);
	    }
	}

    /* Rotate by 270 degrees */
    else if (rotate >= 225 && rotate < 315) {
	if (mirror) {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = y1;
		    y2 = x1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,ny,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s reflected and rotated 270 degrees",filename);
            hputc (header,"HISTORY",history);
	    }
	else {
	    for (y1 = 0; y1 < ny; y1++) {
		for (x1 = 0; x1 < nx; x1++) {
		    x2 = y1;
		    y2 = nx - x1 - 1;
		    movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,ny,x2,y2);
		    }
		}
	    sprintf (history,"Copy of image %s rotated 270 degrees",filename);
            hputc (header,"HISTORY",history);
	    }
	hputi4 (header,"NAXIS1",ny);
	hputi4 (header,"NAXIS2",nx);
	}

    /* If rotating by more than 315 degrees, assume top-bottom reflection */
    else if (rotate >= 315 && mirror) {
	for (y1 = 0; y1 < ny; y1++) {
	    for (x1 = 0; x1 < nx; x1++) {
		x2 = y1;
		y2 = x1;
		movepix (image,bitpix1,nx,x1,y1,rotimage,bitpix2,ny,x2,y2);
		}
	    }
	sprintf (history,"Copy of image %s reflected top to bottom",filename);
        hputc (header,"HISTORY",history);
	}
    
    if (verbose)
	printf ("%s\n",history);

    free (*image0);
    *image0 = rotimage;
    return (0);
}
/* May 29 1996	Change name from rotFITS to RotFITS
 * Jun  4 1996	Fix bug when handling assymetrical images
 * Jun  5 1996	Print filename, not pathname, in history
 * Jun 10 1996	Remove unused variables after running lint
 * Jun 13 1996	Replace image with rotated image
 * Jun 18 1996	Fix formatting bug in history
 *
 * Jul 11 1997	If rotation is 360, flip top bottom if mirror flat is set
 *
 * Feb 23 1998	Do not delete WCS if image not rotated or mirrored
 */
