/*** File findstar.c
 *** February 21. 1996
 *** By Elwood Downey, revised by Doug Mink
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitshead.h"

#define	FSBORDER	16	/* findStars() ignore this much of the edges */
#define	NSTATPIX	20	/* findStars() stats row buff len */
#define	MAXRAD		20	/* Maximum radius we accept as a star */
#define	MAXWALK		20	/* Farthest distance to walk from seed */
#define	SIGMAS		5	/* Stars are this many sigmas above mean */
#define	BURNEDOUT	65535	/* Clamp pixels brighter than this */
#define NITERATE	3	/* Number of iterations for sigma clipping */

static int hotPixel();
static int starRadius();
static void starCentroid();
static int brightWalk ();
double FITSnoise();
double getpix();
void putpix();
void getvec();
void putvec();

/* Find the location and brightest pixel of stars in the given image.
 * Return malloced arrays of x and y and b.
 * N.B. caller must free *xa and *ya and *ba even if we return 0 (but not -1).
 * N.B. we ignore pixels outside FSBORDER.
 * N.B. we ignore isolated hot pixels.
 * return number of stars (might well be 0 :-), or -1 if trouble.
 */
int
findStars (header, image, xa, ya, ba)

char	*header;	/* FITS header */
char	*image;		/* image pixels */
double	**xa, **ya;	/* we set *xa and *ya to memory we malloc */
double	**ba;		/* we set *ba to memory we malloc */

{
    double noise;
    int nstars;
    double minll;
    int bitpix;
    int w, h;
    int x, y;
    int wrap;
    int n;
    int p, np;
    double sp;
    double lmean, lsigma;	/* left and right stats */
    double rmean, rsigma;
    double lll, rll;		/* left and right lower limits*/
    double sumx, sumxx;
    double *svec, *sv;
    double nspix, nspix1;
    int i;

    hgeti4 (header,"NAXIS1", &w);
    hgeti4 (header,"NAXIS2", &h);
    hgeti4 (header,"BITPIX", &bitpix);

    /* Init the malloced arrays.
     * it's ok to do now because we claim caller should always free these.
     */
    *xa = (double *) malloc (sizeof(double));
    *ya = (double *) malloc (sizeof(double));
    *ba = (double *) malloc (sizeof(double));
    svec =  (double *) malloc (w * sizeof (double));

    /* Compute image noise from a central swath */
    noise = FITSnoise (image, bitpix, w, FSBORDER, h/2, w-2*FSBORDER, 10);

    /* Scan for stars based on surrounding local noise figure */
    nstars = 0;
    nspix = NSTATPIX;
    nspix1 = nspix - 1.0;
    for (y = FSBORDER; y < h-FSBORDER; y++) {
        int ipix = 0;

	/* Get one line of the image and fill borders with noise */
	getvec (image, bitpix, w * (y-1), w, svec);
	for (x = 0; x < FSBORDER; x++)
	    svec[x] = noise;
	for (x = w-1; x > w - FSBORDER; x--)
	    svec[x] = noise;

	/* Search row for bright pixels */
	for (x = FSBORDER; x < w-FSBORDER; x++, p++) {

	    /* Redo stats once for every four pixels
	     * N.B. this assumes FSBORDER is a multiple of 4!!
	     */
	    if (ipix++ % 4 == 0) {

		/* Find stats to the left */
		sumx = sumxx = 0;
		sv = svec + x - NSTATPIX;
		for (i = 0; i < NSTATPIX; i++) {
		    sp = *sv++;
		    sumx += sp;
		    sumxx += sp * sp;
		    }
		lsigma = sqrt ((sumxx - sumx*sumx / nspix) / nspix1);
		lmean = sumx / nspix;
		lll = lmean + (SIGMAS * lsigma);
		sv++;

		/* Find stats to the right */
		sumx = sumxx = 0;
		for (i = 0; i < NSTATPIX; i++) {
		    sp = *sv++;
		    sumx += sp;
		    sumxx += sp * sp;
		    }
		rsigma = sqrt ((sumxx - sumx*sumx / nspix) / nspix1);
		rmean = sumx / nspix;
		rll = rmean + SIGMAS * rsigma;

		/* pick lower as noise level */
		minll = lll < rll ? lll : rll;
		}

	    /* Pixel is a candidate if above the noise */
	    if (svec[x] > minll) {
		int sx, sy, r;
		double b;
		int maxr = MAXRAD;
		int maxw = MAXWALK;
		int i, xai, yai;

		/* Ignore hot pixels */
		if (hotPixel (image,bitpix,w, x, y, minll) == 0)
		    continue;

		/* Walkabout a little to find brightest in neighborhood.  */
		if (brightWalk (image, bitpix, w, x, y, maxw, &sx, &sy, &b) < 0)
		    continue;

		/* Ignore really bright stars */
		if (b >= BURNEDOUT)
		    continue;

		/* Do not do the same one again */
		for (i = 0; i < nstars; i++) {
		    xai = (int) ((*xa)[i] + 0.5);
		    yai = (int) ((*ya)[i] + 0.5);
		    if (abs (xai-sx) <= maxr && abs (yai-sy) <= maxr)
			break;
		    }
		if (i < nstars)
		    continue;

		/* Keep it if it is small enough to be a star */
		r = starRadius(image, bitpix, w, sx, sy, b);
		if (r <= maxr) {
		    /* finally! */
		    nstars++;
		    *xa= (double *) realloc(*xa, nstars*sizeof(double));
		    *ya= (double *) realloc(*ya, nstars*sizeof(double));
		    *ba= (double *) realloc(*ba, nstars*sizeof(double));
		    starCentroid (image, bitpix, w, sx, sy, 
				  &(*xa)[nstars-1], &(*ya)[nstars-1]);
		    (*ba)[nstars-1] = b;
		    }
		}
	    }
	}

    free ((char *)image);

    return (nstars);
}


/* Check pixel at x/y for being "hot", ie, a pixel surrounded by noise.
 * If any are greater than pixel at x/y then return -1.
 * Else set the pixel at x/y to llimit and return 0.
 */

static int
hotPixel (image, bitpix, w, x, y, llimit)

char	*image;
int	bitpix;
int	w;
int	x, y;
double	llimit;

{
    double *left;
    double pix1, pix2, pix3;

    pix1 = getpix (image,bitpix,w,x-1,y-1);
    pix2 = getpix (image,bitpix,w,x,y-1);
    pix3 = getpix (image,bitpix,w,x+1,y-1);
    if (pix1 > llimit || pix2 > llimit || pix3 > llimit)
	return (-1);
    pix1 = getpix (image,bitpix,w,x-1,y);
    pix3 = getpix (image,bitpix,w,x+1,y);
    if (pix1 > llimit || pix3 > llimit)
	return (-1);
    pix1 = getpix (image,bitpix,w,x-1,y+1);
    pix2 = getpix (image,bitpix,w,x,y+1);
    pix3 = getpix (image,bitpix,w,x+1,y+1);
    if (pix1 > llimit || pix2 > llimit || pix3 > llimit)
	return (-1);

    putpix (image, bitpix, w, x, y, llimit);
    return (0);
}


/* Compute and return the radius of the star centered at x0, y0.
 * N.B. we assume there is a guard band on the image.
 * N.B. it is expected our caller will reject this thing if r > MAXRAD.
 */

static int
starRadius (imp, bitpix, w, x0, y0, b)

char	*imp;
int	bitpix;
int	w;
int	x0, y0;
double	b;

{
    int r;

    /* compute star's radius.
     * scan in ever-greater circles until find one such that the peak is
     * SIGMAS above the mean at that radius.
     */
    for (r = 2; r <= MAXRAD; r++) {
	double sigma, mean;
	int inrr = r*r;
	int outrr = (r+1)*(r+1);
	int np = 0;
	double sum = 0.0, sum2 = 0.0;
	int x, y;

	for (y = -r; y <= r; y++) { 
	int yrr = y*y;
	for (x = -r; x <= r; x++) {
	    int xyrr = x*x + yrr;
	    if (xyrr >= inrr && xyrr < outrr) {
		double dp;
		dp = getpix (imp,bitpix,w,x0+x,y0+y);
		sum += dp;
		sum2 += dp*dp;
		np++;
	    }
	}
	}

	sigma = sqrt((sum2 - sum*sum/np)/(np-1));
	mean = sum/np;
	if (b > mean + SIGMAS*sigma)
	break;
    }

    return (r);
}

/* Compute the fine location of the star located at [x0,y0].  */

static void
starCentroid (imp, bitpix, w, x0, y0, xp, yp)

char	*imp;
int	bitpix;
int	w;
int	x0, y0;
double	*xp, *yp;

{
    double p1, p2, p22, p3, d;

    /* find maximum of best-fit parabola in each direction.
     * see Bevington, page 210
     */

    p1 = getpix (imp,bitpix,w,x0-1,y0);
    p2 = getpix (imp,bitpix,w,x0,y0);
    p22 = 2*p2;
    p3 = getpix (imp,bitpix,w,x0+1,y0);
    d = p3 - p22 + p1;
    *xp = (d == 0) ? x0 : x0 + 0.5 - (p3 - p2)/d;
    *xp = *xp + 1.0;

    p1 = getpix (imp,bitpix,w,x0,y0-1);
    p3 = getpix (imp,bitpix,w,x0,y0+1);
    d = p3 - p22 + p1;
    *yp = (d == 0) ? y0 : y0 + 0.5 - (p3 - p2)/d;
    *yp = *yp + 1.0;
}


/* given an image and a starting point, walk the gradient to the brightest
 * pixel and return its location. we never take more maxr away.
 * return 0 if find brightest pixel within maxsteps else -1.
 */

static int dx[8]={1,0,-1,1,-1,1,0,-1};
static int dy[8]={1,1,1,0,0,-1,-1,-1};

static int
brightWalk (image, bitpix, w, x0, y0, maxr, xp, yp, bp)

char	*image;
int	bitpix;
int	w;
int	x0;
int	y0;
int	maxr;
int	*xp;
int	*yp;
double	*bp;

{

    double b, tmpb, newb;
    int x, y, x1, y1, i;

    /* start by assuming seed point is brightest */
    b = getpix (image,bitpix,w, x0,y0);
    x = x0;
    y = y0;

    /* walk towards any brighter pixel */
    for (;;) {
	int newx, newy;

    /* Find brightest pixel in 3x3 region */
	newb = b;
	for (i = 0; i < 8; i++) {
	x1 = x + dx[i];
	y1 = y + dy[i];
	tmpb = getpix (image,bitpix,w, x1, y1);
	if (tmpb > newb) {
	    newx = x1;
	    newy = y1;
	    newb = tmpb;
	    }
	}
 
    /* If brightest pixel is one in center of region, quit */
	if (newb == b)
	break;

	x = newx;
	y = newy;
	b = newb;
	if (abs(x-x0) > maxr || abs(y-y0) > maxr)
	return (-1);
    }

    *xp = x;
    *yp = y;
    *bp = b;
    return (0);
}

/* compute stats in the give region of the image of width w pixels.
 * N.B. we do not check bounds.
 */

double
FITSnoise (image, bitpix, w, x0, y0, nx, ny)

char	*image;
int	bitpix;
int	w;
int	x0;
int	y0;
int	nx;
int	ny;

{
    double pmin, pmax, pmean;
    double sd, sd2;
    int x, y;
    int i, n;

    pmin = -1.0e20;
    pmax = 1.0e20;

    for (i = 0; i < NITERATE; i++ ) {
	double sum = 0.0;
	double sum2 = 0.0;
	double dnpix;
	int npix = 0;
	for (y = 0; y < ny; y++) {
	for (x = 0; x < nx; x++) {
	    double p = getpix (image,bitpix,w, x, y);
	    if (p > pmin && p < pmax) {
		sum += p;
		sum2 += p * p;
		npix++;
		}
	    }
	}
	dnpix = (double) npix;
	pmean = sum / dnpix;
	sd2 = (sum2 - (pmean * pmean)) / (dnpix - 1.0);
	sd = sd2 <= 0.0 ? 0.0 : sqrt(sd2);
	pmin = pmean - sd * SIGMAS;
	pmax = pmean + sd * SIGMAS;
	}

    return (pmean);
}


double
getpix (image, bitpix, w, x, y)

char	*image;
int	bitpix;
int	w;
int	x;
int	y;

{
    short *im2;
    int *im4;
    unsigned int *imu;
    float *imr;
    double *imd;

    switch (bitpix) {

	case 16:
	im2 = (short *)image;
	return ((double) im2[(y*w) + x]);
	break;

	case 32:
	im4 = (int *)image;
	return ((double) im4[(y*w) + x]);
	break;

	case -16:
	imu = (unsigned int *)image;
	return ((double) imu[(y*w) + x]);
	break;

	case -32:
	imr = (float *)image;
	return ((double) imr[(y*w) + x]);
	break;

	case -64:
	imd = (double *)image;
	return (imd[(y*w) + x]);
	break;

	default:
	return (0.0);
	}
}


void
putpix (image, bitpix, w, x, y, dpix)

char	*image;
int	bitpix;
int	w;
int	x;
int	y;
double	dpix;

{
    short *im2;
    int *im4;
    unsigned int *imu;
    float *imr;
    double *imd;

    switch (bitpix) {

	case 16:
	    im2 = (short *)image;
	    im2[(y*w) + x] = (short) dpix;
	    break;

	case 32:
	    im4 = (int *)image;
	    im4[(y*w) + x] = (int) dpix;
	    break;

	case -16:
	    imu = (unsigned int *)image;
	    imu[(y*w) + x] = (unsigned int) dpix;
	    break;

	case -32:
	    imr = (float *)image;
	    imr[(y*w) + x] = (float) dpix;
	    break;

	case -64:
	    imd = (double *)image;
	    imd[(y*w) + x] = dpix;
	    break;

	}
    return;
}


void
getvec (image, bitpix, pix1, pixoff, dpix)

char	*image;
int	bitpix;
int	pix1;
int	pixoff;
double	*dpix;

{
    short *im2;
    int *im4;
    unsigned int *imu;
    float *imr;
    double *imd;
    int ipix, pix2;
    double *dp;

    pix2 = pix1 + pixoff;
    dp = dpix;

    switch (bitpix) {

	case 16:
	    im2 = (short *)image;
	    for (ipix = pix1; ipix < pix2; ipix++)
		*dpix++ = (double) *(im2+ipix);
	    break;

	case 32:
	    im4 = (int *)image;
	    for (ipix = pix1; ipix < pix2; ipix++)
		*dpix++ = (double) *(im4+ipix);
	    break;

	case -16:
	    imu = (unsigned int *)image;
	    for (ipix = pix1; ipix < pix2; ipix++)
		*dpix++ = (double) *(imu+ipix);
	    break;

	case -32:
	    imr = (float *)image;
	    for (ipix = pix1; ipix < pix2; ipix++)
		*dpix++ = (double) *(imr+ipix);
	    break;

	case -64:
	    imd = (double *)image;
	    for (ipix = pix1; ipix < pix2; ipix++)
		*dpix++ = *(imd+ipix);
	    break;

	}
    return;
}


void
putvec (image, bitpix, pix1, pixoff, dpix)

char	*image;
int	bitpix;
int	pix1;
int	pixoff;
double	dpix;

{
    short *im2;
    int *im4;
    unsigned int *imu;
    float *imr;
    double *imd;
    int ipix, pix2;

    pix2 = pix1 + pixoff;

    switch (bitpix) {

	case 16:
	    im2 = (short *)image;
	    for (ipix = pix1; ipix < pix2; ipix++)
		*(im2+ipix) = (short) dpix;
	    break;

	case 32:
	    im4 = (int *)image;
	    for (ipix = pix1; ipix < pix2; ipix++)
		*(im4+ipix) = (int) dpix;
	    break;

	case -16:
	    imu = (unsigned int *)image;
	    for (ipix = pix1; ipix < pix2; ipix++)
		*(imu+ipix) = (unsigned int) dpix;
	    break;

	case -32:
	    imr = (float *)image;
	    for (ipix = pix1; ipix < pix2; ipix++)
		*(imr+ipix) = (float) dpix;
	    break;

	case -64:
	    imd = (double *)image;
	    for (ipix = pix1; ipix < pix2; ipix++)
		*(imd+ipix) = (double) dpix;
	    break;

	}
    return;
}
