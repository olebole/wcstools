/*** File libwcs/findstar.c
 *** August 6, 1996
 *** By Elwood Downey, revised by Doug Mink
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitshead.h"

#define	FSBORDER	16	/* findStars() ignore this much of the edges */
#define	NSTATPIX	20	/* findStars() stats row buff len */
#define	MAXWALK		20	/* Farthest distance to walk from seed */
#define	BURNEDOUT	65535	/* Clamp pixels brighter than this */
#define NITERATE	3	/* Number of iterations for sigma clipping */
#define ABS(a) ((a) < 0 ? (-(a)) : (a))

static int hotPixel();
static int starRadius();
static void starCentroid();
static int brightWalk ();
static double FindFlux ();
static void FITSnoise();

static double starsig = 5.0;	/* Stars are this many sigmas above mean */
void setstarsig (sig)
double sig;
{ starsig = sig; return; }

static int maxrad = 20;		/* Maximum radius we accept as a star */
void setmaxrad (max)
int max;
{ maxrad = max; return; }

static double bmin = 10;	/* Minimum peak for star */
void setbmin (min)
double min;
{ bmin = min; return; }

/* Find the location and brightest pixel of stars in the given image.
 * Return malloced arrays of x and y and b.
 * N.B. caller must free *xa and *ya and *ba even if we return 0 (but not -1).
 * N.B. we ignore pixels outside FSBORDER.
 * N.B. we ignore isolated hot pixels.
 * return number of stars (might well be 0 :-), or -1 if trouble.
 */
int
FindStars (header, image, xa, ya, ba, pa, verbose)

char	*header;	/* FITS header */
char	*image;		/* image pixels */
double	**xa, **ya;	/* X and Y coordinates of stars, array returned */
double	**ba;		/* Fluxes of stars in counts, array returned */
int	**pa;		/* Peak counts of stars in counts, array returned */
int	verbose;	/* 1 to print each star's position */

{
    double noise, nsigma;
    int nstars;
    double minll;
    int bitpix;
    int w, h;
    int x, y;
    double xai, yai, bai;
    double lmean, lsigma;	/* left and right stats */
    double rmean, rsigma;
    double lll, rll;		/* left and right lower limits*/
    double sumx, sumxx, dsumxx;
    double *svec, *sv, *svb, *sv1, *sv2, *svx, *svlim;
    double background;
    double rmax, nsig;
    int lwidth, npix;
    int nspix = NSTATPIX;
    int border = FSBORDER;
    int nextline;

    hgeti4 (header,"NAXIS1", &w);
    hgeti4 (header,"NAXIS2", &h);
    hgeti4 (header,"BITPIX", &bitpix);

    /* Allocate the position, flux, and peak intensity arrays
     * it's ok to do now because we claim caller should always free these.
     */
    *xa = (double *) malloc (sizeof(double));
    *ya = (double *) malloc (sizeof(double));
    *ba = (double *) malloc (sizeof(double));
    *pa = (int *) malloc (sizeof(int));

    /* Allocate a buffer to hold one image line */
    svec =  (double *) malloc (w * sizeof (double));

    /* Compute image noise from a central swath */
    FITSnoise (image, bitpix, w, h, w-2*FSBORDER, 10, &noise, &nsigma);

    /* Fill in borders of the image line buffer with noise */
    svlim = svec + w;
    svb = svec + border;
    for (sv = svec; sv < svb; sv++)
	*sv = noise;
    for (sv = svlim - border; sv < svlim; sv++)
	*sv = noise;

    /* Scan for stars based on surrounding local noise figure */
    nstars = 0;
    lwidth = w - (2 * border);
    for (y = border; y < h-border; y++) {
        int ipix = 0;

	/* Get one line of the image minus the noise-filled borders */
	nextline = (w * (y-1)) + border - 1;
	getvec (image, bitpix, nextline, lwidth, svb);

	/* Search row for bright pixels */
	for (x = border; x < w-border; x++) {

	    /* Redo stats once for every four pixels */
	    if (ipix++ % 4 == 0) {

		/* Find stats to the left */
		sumx = sumxx = 0;
		svx = svec + x;
		sv1 = svec + x - nspix;
		if (sv1 < svec) sv1 = svec;
		npix = 0;
		for (sv = sv1; sv < svx; sv++) {
		    sumx += *sv;
		    sumxx += *sv * *sv;
		    npix++;
		    }
		if (npix > 0) {
		    dsumxx = (sumx * sumx) / (double) npix;
		    if (dsumxx > sumxx || npix < 2)
			lsigma = nsigma;
		    else
			lsigma = sqrt ((sumxx - dsumxx) / (double) (npix - 1));
		    if (lsigma < nsigma) lsigma = nsigma;
		    lmean = sumx / (double) npix;
		    }
		else {
		    lmean = noise;
		    lsigma = nsigma;
		    }
		lll = lmean + (starsig * lsigma);
		sv++;

		/* Find stats to the right */
		sumx = sumxx = 0;
		sv2 = svec + x + nspix;
		if (sv2 > svlim) sv2 = svlim;
		npix = 0;
		for (sv = svx + 1; sv < sv2; sv++) {
		    sumx += *sv;
		    sumxx += *sv * *sv;
		    npix++;
		    }
		if (npix > 0) {
		    dsumxx = (sumx * sumx) / (double) npix;
		    if (dsumxx > sumxx || npix < 2)
			rsigma = nsigma;
		    else
			rsigma = sqrt ((sumxx - dsumxx) / (double) (npix - 1));
		    if (rsigma < nsigma) rsigma = nsigma;
		    rmean = sumx / (double) npix;
		    }
		else {
		    lmean = noise;
		    lsigma = nsigma;
		    }
		rll = rmean + starsig * rsigma;

		/* pick lower as noise level */
		minll = lll < rll ? lll : rll;
		}

	    /* Pixel is a candidate if above the noise */
	    if (svec[x] > minll) {
		int sx, sy, r, rf;
		double b;
		int maxr = maxrad;
		int maxw = MAXWALK;
		int i, ix, iy;

		/* Ignore faint stars */
		if (svec[x] < bmin)
		    continue;

		/* Ignore hot pixels */
		if (hotPixel (image,bitpix,w,h, x, y, minll) == 0)
		    continue;

		/* Walkabout a little to find brightest in neighborhood.  */
		if (brightWalk (image, bitpix, w, h, x, y, maxw, &sx, &sy, &b) < 0)
		    continue;

		/* Ignore really bright stars */
		if (b >= BURNEDOUT)
		    continue;

		/* Do not do the same one again */
		for (i = 0; i < nstars; i++) {
		    ix = (int) ((*xa)[i] + 0.5);
		    iy = (int) ((*ya)[i] + 0.5);
		    if (abs (ix-sx) <= maxr && abs (iy-sy) <= maxr)
			break;
		    }
		if (i < nstars)
		    continue;

		/* Keep it if it is small enough to be a star */
		rmax = maxrad;
		nsig = starsig;
		r = starRadius (image, bitpix, w, h, sx, sy, b, rmax, starsig,
			       &background);
		if (r <= maxr) {

		/* Find center of star */
		    nstars++;
		    *xa= (double *) realloc(*xa, nstars*sizeof(double));
		    *ya= (double *) realloc(*ya, nstars*sizeof(double));
		    *ba= (double *) realloc(*ba, nstars*sizeof(double));
		    *pa= (int *) realloc(*pa, nstars*sizeof(int));
		    starCentroid (image, bitpix, w, h, sx, sy, &xai, &yai); 
		    (*xa)[nstars-1] = xai;
		    (*ya)[nstars-1] = yai;
		    (*pa)[nstars-1] = (int) b;

		/* Find size of star for photometry */
		/* (points more than three noise sigma above surroundings) */
		    sx = (int) (xai + 0.5);
		    sy = (int) (yai + 0.5);
		    rmax = 2.0 * (double) maxrad;
		    nsig = starsig * 2.0;
		    rf = starRadius (image, bitpix, w, h, sx, sy, b, rmax, nsig,
				    &background);

		/* Find flux from star */
		    bai = FindFlux (image, bitpix, w, h, sx, sy, rf, background);
		    (*ba)[nstars-1] = bai;
		    if (verbose) {
			fprintf (stderr," %d: (%d %d) -> (%7.3f %7.3f)",
				 nstars, sx, sy, xai, yai);
			fprintf (stderr," %8.1f -> %10.1f  %d -> %d    ",
				 b, bai, r, rf);
			(void)putc (13,stderr);
			}
		    }
		/* else {
		    fprintf (stderr," %d: (%d %d) %d > %d\n",
			     nstars, sx, sy, r, maxrad);
		    } */
		}
	    }
	}

    return (nstars);
}


/* Check pixel at x/y for being "hot", ie, a pixel surrounded by noise.
 * If any are greater than pixel at x/y then return -1.
 * Else set the pixel at x/y to llimit and return 0.
 */

static int
hotPixel (image, bitpix, w, h, x, y, llimit)

char	*image;
int	bitpix;
int	w;
int	h;
int	x, y;
double	llimit;

{
    double pix1, pix2, pix3;

    pix1 = getpix (image,bitpix,w,h,x-1,y-1);
    pix2 = getpix (image,bitpix,w,h,x,y-1);
    pix3 = getpix (image,bitpix,w,h,x+1,y-1);
    if (pix1 > llimit || pix2 > llimit || pix3 > llimit)
	return (-1);
    pix1 = getpix (image,bitpix,w,h,x-1,y);
    pix3 = getpix (image,bitpix,w,h,x+1,y);
    if (pix1 > llimit || pix3 > llimit)
	return (-1);
    pix1 = getpix (image,bitpix,w,h,x-1,y+1);
    pix2 = getpix (image,bitpix,w,h,x,y+1);
    pix3 = getpix (image,bitpix,w,h,x+1,y+1);
    if (pix1 > llimit || pix2 > llimit || pix3 > llimit)
	return (-1);

    putpix (image, bitpix, w, h, x, y, llimit);
    return (0);
}


/* Compute and return the radius of the star centered at x0, y0.
 * A guard band is assumed to exist on the image.
 * Calling program is assumed to reject object if r > rmax.
 */

static int
starRadius (imp, bitpix, w, h, x0, y0, b, rmax, nsig, mean)

char	*imp;
int	bitpix;
int	w;
int	h;
int	x0, y0;
double	b;
double	rmax;
double	nsig;
double	*mean;

{
    double r, sigma, sum2p;

    /* Compute star's radius.
     * Scan in ever-greater circles until find one such that the peak is
     * nsig above the mean at that radius.
     */
    for (r = 2; r <= rmax; r++) {
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
		dp = getpix (imp,bitpix,w,h,x0+x,y0+y);
		sum += dp;
		sum2 += dp*dp;
		np++;
		}
	    }
	}

	*mean = sum / np;
	sum2p = sum * sum / np;
	if (sum2p > sum2 || np < 2)
	    sum2 = sum2p;
	sigma = sqrt ((sum2 - sum2p) / (double)(np - 1));
	if (b > *mean + nsig*sigma)
	    break;
    }

    return (r);
}

/* Compute the fine location of the star located at [x0,y0].  */

static void
starCentroid (imp, bitpix, w, h, x0, y0, xp, yp)

char	*imp;
int	bitpix;
int	w;
int	x0, y0;
double	*xp, *yp;

{
    double p1, p2, p22, p3, d;

    /* Find maximum of best-fit parabola in each direction.
     * see Bevington, page 210
     */

    p1 = getpix (imp,bitpix,w,h,x0-1,y0);
    p2 = getpix (imp,bitpix,w,h,x0,y0);
    p22 = 2*p2;
    p3 = getpix (imp,bitpix,w,h,x0+1,y0);
    d = p3 - p22 + p1;
    *xp = (d == 0) ? x0 : x0 + 0.5 - (p3 - p2)/d;
    *xp = *xp + 1.0;

    p1 = getpix (imp,bitpix,w,h,x0,y0-1);
    p3 = getpix (imp,bitpix,w,h,x0,y0+1);
    d = p3 - p22 + p1;
    *yp = (d == 0) ? y0 : y0 + 0.5 - (p3 - p2)/d;
    *yp = *yp + 1.0;
}


/* Given an image and a starting point, walk the gradient to the brightest
 * pixel and return its location. we never take more maxr away.
 * Return 0 if find brightest pixel within maxsteps else -1.
 */

static int dx[8]={1,0,-1,1,-1,1,0,-1};
static int dy[8]={1,1,1,0,0,-1,-1,-1};

static int
brightWalk (image, bitpix, w, h, x0, y0, maxr, xp, yp, bp)

char	*image;
int	bitpix;
int	w;
int	h;
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
    b = getpix (image,bitpix,w,h, x0,y0);
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
	tmpb = getpix (image,bitpix,w,h, x1, y1);
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

/* Compute stats in the give region of the image of width w pixels.
 * Bounds are not checked.
 */

static void
FITSnoise (image, bitpix, w, h, nx, ny, mean, sigma)

char	*image;
int	bitpix;
int	w;
int	h;
int	nx;
int	ny;
double	*mean;
double	*sigma;

{
    double pmin, pmax, pmean;
    double sd, sd2;
    int x, y;
    int i;

    pmin = -1.0e20;
    pmax = 1.0e20;

    for (i = 0; i < NITERATE; i++ ) {
	double sum = 0.0;
	double sum2 = 0.0;
	double dnpix;
	int npix = 0;
	for (y = 0; y < ny; y++) {
	for (x = 0; x < nx; x++) {
	    double p = getpix (image,bitpix,w,h, x, y);
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
	pmin = pmean - sd * starsig;
	pmax = pmean + sd * starsig;
	}

    *mean = pmean;
    *sigma = sd;
    return;
}


/* Find total flux within a circular region minus a mean background level */

static double
FindFlux (image, bitpix, w, h, x0, y0, r, background)

char	*image;
int	bitpix;
int	w;
int	x0;
int	y0;
int	r;
double	background;
{
    double sum = 0.0;
    int x, y, x1, x2, y1, y2, yy, xxyy;
    int rr = r * r;
    double dp;

/* Keep X within image */
    x1 = -r;
    if (x0-r < 0)
	x1 = 0;
    x2 = r;
    if (x0+r > 0)
	x2 = w;

/* Keep Y within image */
    y1 = -r;
    if (y0-r < 0)
	y1 = 0;
    y2 = r;
    if (y0+r > 0)
	y2 = h;

/* Integrate circular region around a star */
    for (y = y1; y <= y2; y++) { 
	yy = y*y;
	for (x = x1; x <= x2; x++) {
	    xxyy = x*x + yy;
	    if (xxyy <= rr) {
		dp = getpix (image, bitpix, w,h, x0+x, y0+y);
		if (dp > background) {
		    sum += dp - background;
		    }
		}
	    }
	}

    return (sum);
}
/* May 21 1996	Return peak flux in counts
 * May 22 1996	Add arguments so GETPIX and PUTPIX can check coordinates
 * Jun  6 1996	Change name from findStars to FindStars
 * Jun 12 1996	Remove unused variables after using lint
 * Jun 13 1996	Removed leftover free of image
 * Aug  6 1996	Fixed small defects after lint
 */
