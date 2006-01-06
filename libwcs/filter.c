/* File libwcs/filter.c
 * October 25, 2005
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

/* Slide box across an image */

#include <string.h>             /* NULL, strlen, strstr, strcpy */
#include <stdio.h>
#include <stdlib.h>
#include "fitshead.h"

#define MEDIAN 1
#define MEAN 2
#define GAUSSIAN 3

char *medfilt();
char *meanfilt();

short medpixi2();
short meanpixi2();
int medpixi4();
int meanpixi4();
float medpixr4();
float meanpixr4();
double medpixr8();
double meanpixr8();
static double fwidth = 1.0;		/* Filter halfwidth, if relevant */

char *
FiltFITS (header, image, filter, xsize, ysize, nlog)

char	*header;	/* Image header */
char	*image;		/* Image bytes to be filtered */
int	filter;		/* Smoothing filter (median,mean,gaussian) */
int	xsize;		/* Number of pixels in x (odd, horizontal) */
int	ysize;		/* Number of pixels in y (odd, vertical) */
int	nlog;		/* Logging interval in lines */

{
    if (filter == MEDIAN)
	return (medfilt (image, header, xsize, ysize, nlog));
    else
	return (meanfilt (image, header, xsize, ysize, nlog));
}


/* Median filter an image */

static short *vi2;	/* Working vector to sort for median */
static int *vi4;	/* Working vector to sort for median */
static float *vr4;	/* Working vector to sort for median */
static double *vr8;	/* Working vector to sort for median */

char *
medfilt (buff, header, ndx, ndy, nlog)

char	*buff;	/* Image buffer */
char	*header; /* FITS image header */
int	ndx;	/* Number of columns over which to compute the median */
int	ndy;	/* Number of rows over which to compute the median */
int	nlog;	/* Logging interval in pixels */

{
char	*buffret;
int	nx,ny;	/* Number of columns and rows in image */
int	ix,iy;	/* Pixel around which to compute median */
int	npix;	/* Number of pixels in image */
int	bitpix;	/* Number of bits per pixel (<0=floating point) */
int	naxes;

    hgeti4 (header, "BITPIX", &bitpix);
    hgeti4 (header, "NAXIS", &naxes);
    hgeti4 (header, "NAXIS1", &nx);
    if (naxes > 1)
	hgeti4 (header, "NAXIS2", &ny);
    else
	ny = 1;
    npix = nx * ny;

    buffret = NULL;
    if (bitpix == 16) {
	short *b, *buffout;
	vi2 = NULL;
	buffret = (char *) calloc (npix, sizeof (short));
	buffout = (short *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b++ = medpixi2 (buff, ix, iy, nx, ny, ndx, ndy);
		}
	    if ((iy+1)%nlog == 0)
		fprintf (stderr,"MEDFILT: %d lines filtered\r", iy+1);
	    }
	fprintf (stderr,"\n");
	free (vi2);
	vi2 = NULL;
	}
    else if (bitpix == 32) {
	int *b, *buffout;
	vi4 = NULL;
	buffret = (char *) calloc (npix, sizeof (int));
	buffout = (int *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b++ = medpixi4 (buff, ix, iy, nx, ny, ndx, ndy);
		}
	    if ((iy+1)%nlog == 0)
		fprintf (stderr,"MEDFILT: %d lines filtered\r", iy+1);
	    }
	fprintf (stderr,"\n");
	free (vi4);
	vi4 = NULL;
	}
    else if (bitpix == -32) {
	float *b, *buffout;
	buffret = (char *) calloc (npix, sizeof (float));
	buffout = (float *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b++ = medpixr4 (buff, ix, iy, nx, ny, ndx, ndy);
		}
	    if ((iy+1)%nlog == 0)
		fprintf (stderr,"MEDFILT: %d lines filtered\r", iy+1);
	    }
	fprintf (stderr,"\n");
	free (vr4);
	vr4 = NULL;
	}
    else if (bitpix == -64) {
	double *b, *buffout;
	buffret = (char *) calloc (npix, sizeof (double));
	buffout = (double *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b++ = medpixr8 (buff, ix, iy, nx, ny, ndx, ndy);
		}
	    if ((iy+1)%nlog == 0)
		fprintf (stderr,"MEDFILT: %d lines filtered\r", iy+1);
	    }
	fprintf (stderr,"\n");
	free (vr8);
	vr8 = NULL;
	}
    return (buffret);
}


/* Compute median of rectangular group of pixels */

short
medpixi2 (x, ix, iy, nx, ny, ndx, ndy)

short	*x;	/* Image buffer */
int	ix,iy;	/* Pixel around which to compute median */
int	nx,ny;	/* Number of columns and rows in image */
int	ndx;	/* Number of columns over which to compute the median */
int	ndy;	/* Number of rows over which to compute the median */

{
    short xx, *vecj, *img;
    int n, l, ir, i, j, n2;
    int  nx2, ny2, npix;
    int jx, jx1, jx2, jy, jy1, jy2;
    int nomed;

    /* Allocate working buffer if it hasn't already been allocated */
    npix = ndx * ndy;
    if (vi2 == NULL) {
	vi2 = (short *) calloc (npix, sizeof (short));
	if (vi2 == NULL) {
	    fprintf (stderr, "MEDPIXI2: Could not allocate %d-pixel buffer\n");
	    return (0);
	    }
	}

    n = ndx * ndy;
    if (n <= 0)
	return (0.0);
    else if (n == 1)
	return (*(x + (iy * ny) + ix));

    /* Compute limits for this pixel */
    nx2 = ndx / 2;
    jx1 = ix - nx2;
    if (jx1 < 0)
	jx1 = 0;
    jx2 = ix + nx2 + 1;
    if (jx2 > nx)
	jx2 = nx;
    ny2 = ndy / 2;
    jy1 = iy - ny2;
    if (jy1 < 0)
	jy1 = 0;
    jy2 = iy + ny2 + 1;
    if (jy2 > ny)
	jy2 = ny;

    /* Compute actual number of pixels used for this pixel */
    n = (jx2 - jx1 + 1) * (jy2 - jy1 + 1);

    /* Set up working vector for this pixel */
    vecj = vi2;
    for (jy = jy1; jy < jy2; jy++) {
	img = x + (jy * nx) + jx1;
	for (jx = jx1; jx < jx2; jx++) {
	    *vecj++ = *img++;
	    }
	}

    /* Sort numbers in working vector */
    for (j = 2; j <= n; j++) {
	xx = vi2[j];
	i = j - 1;
	while (i > 0 && vi2[i] > xx) {
	    vi2[i+1] = vi2[i];
	    i--;
	    }
	vi2[i+1] = xx;
	}

    /* Middle number is the median */
    return (vi2[n/2]);
}


/* Compute median of rectangular group of pixels */

int
medpixi4 (x, ix, iy, nx, ny, ndx, ndy)

int	*x;	/* Image buffer */
int	ix,iy;	/* Pixel around which to compute median */
int	nx,ny;	/* Number of columns and rows in image */
int	ndx;	/* Number of columns over which to compute the median */
int	ndy;	/* Number of rows over which to compute the median */

{
    int xx, *vecj, *img;
    int n, l, ir, i, j, n2;
    int  nx2, ny2, npix;
    int jx, jx1, jx2, jy, jy1, jy2;
    int nomed;

    /* Allocate working buffer if it hasn't already been allocated */
    npix = ndx * ndy;
    if (vi4 == NULL) {
	vi4 = (int *) calloc (npix, sizeof (int));
	if (vi4 == NULL) {
	    fprintf (stderr, "MEDIANI4: Could not allocate %d-pixel buffer\n");
	    return (0);
	    }
	}

    n = ndx * ndy;
    if (n <= 0)
	return (0.0);
    else if (n == 1)
	return (*(x + (iy * ny) + ix));

    /* Compute limits for this pixel */
    nx2 = ndx / 2;
    jx1 = ix - nx2;
    if (jx1 < 0)
	jx1 = 0;
    jx2 = ix + nx2 + 1;
    if (jx2 > nx)
	jx2 = nx;
    ny2 = ndy / 2;
    jy1 = iy - ny2;
    if (jy1 < 0)
	jy1 = 0;
    jy2 = iy + ny2 + 1;
    if (jy2 > ny)
	jy2 = ny;

    /* Compute actual number of pixels used for this pixel */
    n = (jx2 - jx1 + 1) * (jy2 - jy1 + 1);

    /* Set up working vector for this pixel */
    vecj = vi4;
    for (jy = jy1; jy < jy2; jy++) {
	img = x + (jy * nx) + jx1;
	for (jx = jx1; jx < jx2; jx++) {
	    *vecj++ = *img++;
	    }
	}

    /* Sort numbers in working vector */
    for (j = 2; j <= n; j++) {
	xx = vi4[j];
	i = j - 1;
	while (i > 0 && vi4[i] > xx) {
	    vi4[i+1] = vi4[i];
	    i--;
	    }
	vi4[i+1] = xx;
	}

    /* Middle number is the median */
    return (vi4[n/2]);
}


float
medpixr4 (x, ix, iy, nx, ny, ndx, ndy)

float	*x;	/* Image buffer */
int	ix,iy;	/* Pixel around which to compute median */
int	nx,ny;	/* Number of columns and rows in image */
int	ndx;	/* Number of columns over which to compute the median */
int	ndy;	/* Number of rows over which to compute the median */

{
    float xx, *vecj, *img;
    int n, l, ir, i, j, n2;
    int  nx2, ny2, npix;
    int jx, jx1, jx2, jy, jy1, jy2;
    int nomed;

    /* Allocate working buffer if it hasn't already been allocated */
    npix = ndx * ndy;
    if (vr4 == NULL) {
	vr4 = (float *) calloc (npix, sizeof (float));
	if (vr4 == NULL) {
	    fprintf (stderr, "MEDIANR4: Could not allocate %d-pixel buffer\n");
	    return (0);
	    }
	}

    n = ndx * ndy;
    if (n <= 0)
	return (0.0);
    else if (n == 1)
	return (*(x + (iy * ny) + ix));

    /* Compute limits for this pixel */
    nx2 = ndx / 2;
    jx1 = ix - nx2;
    if (jx1 < 0)
	jx1 = 0;
    jx2 = ix + nx2 + 1;
    if (jx2 > nx)
	jx2 = nx;
    ny2 = ndy / 2;
    jy1 = iy - ny2;
    if (jy1 < 0)
	jy1 = 0;
    jy2 = iy + ny2 + 1;
    if (jy2 > ny)
	jy2 = ny;

    /* Compute actual number of pixels used for this pixel */
    n = (jx2 - jx1) * (jy2 - jy1);

    /* Set up working vector for this pixel */
    vecj = vr4;
    for (jy = jy1; jy < jy2; jy++) {
	img = x + (jy * nx) + jx1;
	for (jx = jx1; jx < jx2; jx++) {
	    *vecj++ = *img++;
	    }
	}

    /* Sort numbers in working vector */
    for (j = 2; j <= n; j++) {
	xx = vr4[j];
	i = j - 1;
	while (i > 0 && vr4[i] > xx) {
	    vr4[i+1] = vr4[i];
	    i--;
	    }
	vr4[i+1] = xx;
	}

    /* Middle number is the median */
    return (vr4[n/2]);
}


double
medpixr8 (x, ix, iy, nx, ny, ndx, ndy)

double	*x;	/* Image buffer */
int	ix,iy;	/* Pixel around which to compute median */
int	nx,ny;	/* Number of columns and rows in image */
int	ndx;	/* Number of columns over which to compute the median */
int	ndy;	/* Number of rows over which to compute the median */

{
    double xx, *vecj, *img;
    int n, l, ir, i, j, n2;
    int  nx2, ny2, npix;
    int jx, jx1, jx2, jy, jy1, jy2;
    int nomed;

    /* Allocate working buffer if it hasn't already been allocated */
    npix = ndx * ndy;
    if (vr8 == NULL) {
	vr8 = (double *) calloc (npix, sizeof (double));
	if (vr8 == NULL) {
	    fprintf (stderr, "MEDIANR8: Could not allocate %d-pixel buffer\n");
	    return (0);
	    }
	}

    n = ndx * ndy;
    if (n <= 0)
	return (0.0);
    else if (n == 1)
	return (*(x + (iy * ny) + ix));

    /* Compute limits for this pixel */
    nx2 = ndx / 2;
    jx1 = ix - nx2;
    if (jx1 < 0)
	jx1 = 0;
    jx2 = ix + nx2 + 1;
    if (jx2 > nx)
	jx2 = nx;
    ny2 = ndy / 2;
    jy1 = iy - ny2;
    if (jy1 < 0)
	jy1 = 0;
    jy2 = iy + ny2 + 1;
    if (jy2 > ny)
	jy2 = ny;

    /* Compute actual number of pixels used for this pixel */
    n = (jx2 - jx1 + 1) * (jy2 - jy1 + 1);

    /* Set up working vector for this pixel */
    vecj = vr8;
    for (jy = jy1; jy < jy2; jy++) {
	img = x + (jy * nx) + jx1;
	for (jx = jx1; jx < jx2; jx++) {
	    *vecj++ = *img++;
	    }
	}

    /* Sort numbers in working vector */
    for (j = 2; j <= n; j++) {
	xx = vr8[j];
	i = j - 1;
	while (i > 0 && vr8[i] > xx) {
	    vr8[i+1] = vr8[i];
	    i--;
	    }
	vr8[i+1] = xx;
	}

    /* Middle number is the median */
    return (vr8[n/2]);
}


/* Mean filter an image */

char *
meanfilt (buff, header, ndx, ndy, nlog)

char	*buff;	/* Image buffer */
char	*header; /* FITS image header */
int	ndx;	/* Number of columns over which to compute the mean */
int	ndy;	/* Number of rows over which to compute the mean */
int	nlog;	/* Logging interval in pixels */

{
char	*buffret;
int	nx,ny;	/* Number of columns and rows in image */
int	ix,iy;	/* Pixel around which to compute mean */
int	npix;	/* Number of pixels in image */
int	bitpix;	/* Number of bits per pixel (<0=floating point) */
int	naxes;

    hgeti4 (header, "BITPIX", &bitpix);
    hgeti4 (header, "NAXIS", &naxes);
    hgeti4 (header, "NAXIS1", &nx);
    if (naxes > 1)
	hgeti4 (header, "NAXIS2", &ny);
    else
	ny = 1;
    npix = nx * ny;

    buffret = NULL;
    if (bitpix == 16) {
	short *b, *buffout;
	vi2 = NULL;
	buffret = (char *) calloc (npix, sizeof (short));
	buffout = (short *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b++ = meanpixi2 (buff, ix, iy, nx, ny, ndx, ndy);
		}
	    if ((iy+1)%nlog == 0)
		fprintf (stderr,"MEANFILT: %d lines filtered\r", iy+1);
	    }
	fprintf (stderr,"\n");
	free (vi2);
	vi2 = NULL;
	}
    else if (bitpix == 32) {
	int *b, *buffout;
	vi4 = NULL;
	buffret = (char *) calloc (npix, sizeof (int));
	buffout = (int *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b++ = meanpixi4 (buff, ix, iy, nx, ny, ndx, ndy);
		}
	    if ((iy+1)%nlog == 0)
		fprintf (stderr,"MEANFILT: %d lines filtered\r", iy+1);
	    }
	fprintf (stderr,"\n");
	free (vi4);
	vi4 = NULL;
	}
    else if (bitpix == -32) {
	float *b, *buffout;
	buffret = (char *) calloc (npix, sizeof (float));
	buffout = (float *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b++ = meanpixr4 (buff, ix, iy, nx, ny, ndx, ndy);
		}
	    if ((iy+1)%nlog == 0)
		fprintf (stderr,"MEANFILT: %d lines filtered\r", iy+1);
	    }
	fprintf (stderr,"\n");
	free (vr4);
	vr4 = NULL;
	}
    else if (bitpix == -64) {
	double *b, *buffout;
	buffret = (char *) calloc (npix, sizeof (double));
	buffout = (double *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b++ = meanpixr8 (buff, ix, iy, nx, ny, ndx, ndy);
		}
	    if ((iy+1)%nlog == 0)
		fprintf (stderr,"MEANFILT: %d lines filtered\r", iy+1);
	    }
	fprintf (stderr,"\n");
	free (vr8);
	vr8 = NULL;
	}
    return (buffret);
}


/* Compute mean of rectangular group of pixels */

short
meanpixi2 (x, ix, iy, nx, ny, ndx, ndy)

short	*x;	/* Image buffer */
int	ix,iy;	/* Pixel around which to compute median */
int	nx,ny;	/* Number of columns and rows in image */
int	ndx;	/* Number of columns over which to compute the median */
int	ndy;	/* Number of rows over which to compute the median */

{
    double sum;
    short *img;
    int n, n2, nx2, ny2;
    int jx, jx1, jx2, jy, jy1, jy2;

    n = ndx * ndy;
    if (n <= 0)
	return (0.0);
    else if (n == 1)
	return (*(x + (iy * ny) + ix));

    /* Compute limits for this pixel */
    nx2 = ndx / 2;
    jx1 = ix - nx2;
    if (jx1 < 0)
	jx1 = 0;
    jx2 = ix + nx2 + 1;
    if (jx2 > nx)
	jx2 = nx;
    ny2 = ndy / 2;
    jy1 = iy - ny2;
    if (jy1 < 0)
	jy1 = 0;
    jy2 = iy + ny2 + 1;
    if (jy2 > ny)
	jy2 = ny;

    /* Compute actual number of pixels used for this pixel */
    n = (jx2 - jx1 + 1) * (jy2 - jy1 + 1);

    /* Compute total counts around this pixel */
    sum = 0.0;
    for (jy = jy1; jy < jy2; jy++) {
	img = x + (jy * nx) + jx1;
	for (jx = jx1; jx < jx2; jx++) {
	    sum = sum + (double) *img++;
	    }
	}

    return ((short) (sum / (double) n));
}


int
meanpixi4 (x, ix, iy, nx, ny, ndx, ndy)

int	*x;	/* Image buffer */
int	ix,iy;	/* Pixel around which to compute median */
int	nx,ny;	/* Number of columns and rows in image */
int	ndx;	/* Number of columns over which to compute the median */
int	ndy;	/* Number of rows over which to compute the median */

{
    double sum;
    int *img;
    int n, n2, nx2, ny2;
    int jx, jx1, jx2, jy, jy1, jy2;

    n = ndx * ndy;
    if (n <= 0)
	return (0.0);
    else if (n == 1)
	return (*(x + (iy * ny) + ix));

    /* Compute limits for this pixel */
    nx2 = ndx / 2;
    jx1 = ix - nx2;
    if (jx1 < 0)
	jx1 = 0;
    jx2 = ix + nx2 + 1;
    if (jx2 > nx)
	jx2 = nx;
    ny2 = ndy / 2;
    jy1 = iy - ny2;
    if (jy1 < 0)
	jy1 = 0;
    jy2 = iy + ny2 + 1;
    if (jy2 > ny)
	jy2 = ny;

    /* Compute actual number of pixels used for this pixel */
    n = (jx2 - jx1 + 1) * (jy2 - jy1 + 1);

    /* Set up working vector for this pixel */
    sum = 0.0;
    for (jy = jy1; jy < jy2; jy++) {
	img = x + (jy * nx) + jx1;
	for (jx = jx1; jx < jx2; jx++) {
	    sum = sum + (double) *img++;
	    }
	}

    return ((int) (sum / (double) n));
}


float
meanpixr4 (x, ix, iy, nx, ny, ndx, ndy)

float	*x;	/* Image buffer */
int	ix,iy;	/* Pixel around which to compute median */
int	nx,ny;	/* Number of columns and rows in image */
int	ndx;	/* Number of columns over which to compute the median */
int	ndy;	/* Number of rows over which to compute the median */

{
    double sum;
    float *img;
    int n, n2, nx2, ny2;
    int jx, jx1, jx2, jy, jy1, jy2;

    n = ndx * ndy;
    if (n <= 0)
	return (0.0);
    else if (n == 1)
	return (*(x + (iy * ny) + ix));

    /* Compute limits for this pixel */
    nx2 = ndx / 2;
    jx1 = ix - nx2;
    if (jx1 < 0)
	jx1 = 0;
    jx2 = ix + nx2 + 1;
    if (jx2 > nx)
	jx2 = nx;
    ny2 = ndy / 2;
    jy1 = iy - ny2;
    if (jy1 < 0)
	jy1 = 0;
    jy2 = iy + ny2 + 1;
    if (jy2 > ny)
	jy2 = ny;

    /* Compute actual number of pixels used for this pixel */
    n = (jx2 - jx1) * (jy2 - jy1);

    /* Set up working vector for this pixel */
    sum = 0.0;
    for (jy = jy1; jy < jy2; jy++) {
	img = x + (jy * nx) + jx1;
	for (jx = jx1; jx < jx2; jx++) {
	    sum = sum + (double) *img++;
	    }
	}

    return ((float) (sum / (double) n));
}


double
meanpixr8 (x, ix, iy, nx, ny, ndx, ndy)

double	*x;	/* Image buffer */
int	ix,iy;	/* Pixel around which to compute median */
int	nx,ny;	/* Number of columns and rows in image */
int	ndx;	/* Number of columns over which to compute the median */
int	ndy;	/* Number of rows over which to compute the median */

{
    double *img;
    double sum;
    int n, n2, nx2, ny2;
    int jx, jx1, jx2, jy, jy1, jy2;

    n = ndx * ndy;
    if (n <= 0)
	return (0.0);
    else if (n == 1)
	return (*(x + (iy * ny) + ix));

    /* Compute limits for this pixel */
    nx2 = ndx / 2;
    jx1 = ix - nx2;
    if (jx1 < 0)
	jx1 = 0;
    jx2 = ix + nx2 + 1;
    if (jx2 > nx)
	jx2 = nx;
    ny2 = ndy / 2;
    jy1 = iy - ny2;
    if (jy1 < 0)
	jy1 = 0;
    jy2 = iy + ny2 + 1;
    if (jy2 > ny)
	jy2 = ny;

    /* Compute actual number of pixels used for this pixel */
    n = (jx2 - jx1 + 1) * (jy2 - jy1 + 1);

    /* Set up working vector for this pixel */
    sum = 0.0;
    for (jy = jy1; jy < jy2; jy++) {
	img = x + (jy * nx) + jx1;
	for (jx = jx1; jx < jx2; jx++) {
	    sum = sum + *img++;
	    }
	}

    return (sum / (double) n);
}

/* Oct 25 2005	New subroutine translated from Fortran imlib/smooth.f
 */
