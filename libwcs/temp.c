

/* Mean filter an image */

char *
meanfilt (buff, header, ndx, ndy, nlog)

char	*buff;	/* Image buffer */
char	*header; /* FITS image header */
int	ndx;	/* Number of columns over which to compute the mean */
int	ndy;	/* Number of rows over which to compute the mean */
int	nlog;	/* Logging interval in pixels */

{
char	*bufret;
int	nx,ny;	/* Number of columns and rows in image */
int	ix,iy;	/* Pixel around which to compute mean */
int	npix;	/* Number of pixels in image */
int	bitpix;	/* Number of bits per pixel (<0=floating point) */

    hgeti4 (header, "BITPIX", &bitpix);
    hgeti4 (header, "NAXIS", &naxes);
    hgeti4 (header, "NAXIS1", &nx);
    if (naxes > 1)
	hgeti4 (header, "NAXIS2", &ny);
    else
	ny = 1;
    npix = nx * ny;

    bufret = NULL;
    if (bitpix == 16) {
	short *b, *buffout;
	vi2 = NULL;
	buffret = (char *) calloc (npix, sizeof (short));
	buffout = (short *) buffret;
	b = buffout;
	for (iy = 0; iy < ny; iy++) {
	    for (ix = 0; ix < nx; ix++) {
		*b = meanpixi2 (buff, ix, iy, nx, ny, ndx, ndy);
		}
	    }
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
		*b = meanpixi4 (buff, ix, iy, nx, ny, ndx, ndy);
		}
	    }
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
		*b = meanpixr4 (buff, ix, iy, nx, ny, ndx, ndy);
		}
	    }
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
		*b = meanpixr8 (buff, ix, iy, nx, ny, ndx, ndy);
		}
	    }
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
    if (jx2 > ny)
	jx2 = ny;
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
	img = x + (jy * ny) + nx;
	for (jx = jx1; jx < jx2; jx++) {
	    sum = sum + (double) *img++;
	    }
	}

    return ((short) (sum / (double) n));
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
    if (jx2 > ny)
	jx2 = ny;
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
	img = x + (jy * ny) + nx;
	for (jx = jx1; jx < jx2; jx++) {
	    sum = sum + (double) *img++;
	    }
	}

    return ((int) (sum / (double) n));
}


int
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
    if (jx2 > ny)
	jx2 = ny;
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
	img = x + (jy * ny) + nx;
	for (jx = jx1; jx < jx2; jx++) {
	    sum = sum + (double) *img++;
	    }
	}

    return ((float) (sum / (double) n));
}


int
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
    if (jx2 > ny)
	jx2 = ny;
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
	img = x + (jy * ny) + nx;
	for (jx = jx1; jx < jx2; jx++) {
	    sum = sum + *img++;
	    }
	}

    return (sum / (double) n);
}
