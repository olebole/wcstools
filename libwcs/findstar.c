/*** File libwcs/findstar.c
 *** May 27, 1998
 *** By Elwood Downey, revised by Doug Mink
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fitshead.h"
#include "imio.h"
#include "lwcs.h"

#define ABS(a) ((a) < 0 ? (-(a)) : (a))

extern int daoread();

static int HotPixel();
static int starRadius();
static void starCentroid();
static int BrightWalk ();
static double FindFlux ();
static void mean2d();
static void mean1d();

/* Stars must be at least this many standard deviations above the mean */
static double starsig = STARSIGMA;
void setstarsig (sig)
double sig;
{ starsig = sig; return; }

/* Ignore this much of the edge */
static int fsborder = BORDER;
void setborder (brd)
double brd;
{ fsborder = brd; return; }

/* Set input catalog for image stars */
static char imcatname[32] = "";
void setimcat (cat)
char *cat;
{strcpy (imcatname, cat); return; }

/* Get input catalog for image stars */
char *getimcat ()
{return (imcatname); }

static int minsep = MINSEP;	/* Minimum separation for stars */
static int minrad = MINRAD;	/* Minimum radius for a star */
static int maxrad = MAXRAD;	/* Maximum radius for a star */
void setmaxrad (max)
int max;
{ maxrad = max; return; }

static double bmin = MINPEAK;	/* Minimum peak for a star */
void setbmin (min)
double min;
{ bmin = min; return; }

/* Find the location and brightest pixel of stars in the given image.
 * Return malloced arrays of x and y and b.
 * N.B. Caller must free *xa and *ya and *ba even if 0 stars are returned.
 * N.B. Pixels outside fsborder are ignored.
 * N.B. Isolated hot pixels are ignored.
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
    int x, y, x1, x2, y1, y2;
    double xai, yai, bai;
    double lmean, lsigma;	/* left and right stats */
    double rmean, rsigma;
    double lll, rll;		/* left and right lower limits*/
    double minsig;
    double *svec, *sv, *svb, *sv1, *sv2, *svlim;
    double background;
    double rmax;
    int lwidth;
    int nspix = NSTATPIX;
    int ispix = ISTATPIX;
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

    /* Read star list from file */
    if (imcatname[0] != 0) {
	int nlog = 0;
	if (verbose) nlog = 10;
	nstars = daoread (imcatname, xa, ya, ba, pa, nlog);
	return (nstars);
	}

    /* Allocate a buffer to hold one image line */
    svec =  (double *) malloc (w * sizeof (double));

    /* Compute image noise from a central swath */
    x1 = fsborder;
    x2 = w - fsborder;
    if (x1 > x2) {
	x1 = 1;
	x2 = w;
	}
    y1 = (h / 2) - 25;
    if (y1 < 1)
	y1 = 0;
    y2 = (h / 2) + 25;
    if (y2 > h)
	y2 = h;
    mean2d (image, bitpix, w, h, x1, x2, y1, y2, &noise, &nsigma);

    /* Fill in borders of the image line buffer with noise */
    svlim = svec + w;
    svb = svec + fsborder;
    for (sv = svec; sv < svb; sv++)
	*sv = noise;
    for (sv = svlim - fsborder; sv < svlim; sv++)
	*sv = noise;

    /* Scan for stars based on surrounding local noise figure */
    nstars = 0;
    lwidth = w - (2 * fsborder);
    for (y = fsborder; y < h-fsborder; y++) {
        int ipix = 0;

	/* Get one line of the image minus the noise-filled borders */
	nextline = (w * (y-1)) + fsborder - 1;
	getvec (image, bitpix, nextline, lwidth, svb);

	/* Search row for bright pixels */
	for (x = fsborder; x < w-fsborder; x++) {

	    /* Redo stats once for every several pixels */
	    if (ipix++ % ispix == 0) {

		/* Find stats to the left */
		sv1 = svec + x - nspix;
		if (sv1 < svec) sv1 = svec;
		sv2 = svec + x;
		if (sv2 > sv1+1)
		    mean1d (sv1, sv2, &lmean, &lsigma);
		else {
		    lmean = noise;
		    lsigma = nsigma;
		    }
		lll = lmean + (starsig * lsigma);
		sv++;

		/* Find stats to the right */
		sv1 = svec + x;
		sv2 = svec + x + nspix;
		if (sv2 > svlim) sv2 = svlim;
		if (sv2 > sv1+2)
		    mean1d (sv1, sv2, &rmean, &rsigma);
		else {
		    rmean = noise;
		    rsigma = nsigma;
		    }
		rll = rmean + starsig * rsigma;

		/* pick lower as noise level */
		minll = lll < rll ? lll : rll;
		minsig = lsigma < rsigma ? lsigma : rsigma;
		}

	    /* Pixel is a candidate if above the noise */
	    if (svec[x] > minll) {
		int sx, sy, r, rf;
		double b;
		int maxw = MAXWALK;
		int i, ix, iy;

		/* Ignore faint stars */
		if (svec[x] < bmin)
		    continue;

		/* Ignore hot pixels */
		if (!HotPixel (image,bitpix,w,h, x, y, minll))
		    continue;

		/* Walkabout a little to find brightest in neighborhood.  */
		if (BrightWalk (image, bitpix, w, h, x, y, maxw, &sx, &sy, &b) < 0)
		    continue;

		/* Ignore really bright stars */
		if (b >= BURNEDOUT)
		    continue;

		/* Do not do the same one again */
		for (i = 0; i < nstars; i++) {
		    ix = (int) ((*xa)[i] + 0.5);
		    iy = (int) ((*ya)[i] + 0.5);
		    if (abs (ix-sx) <= minsep && abs (iy-sy) <= minsep)
			break;
		    }
		if (i < nstars)
		    continue;

		/* Keep it if it is within the size range for stars */
		rmax = maxrad;
		r = starRadius (image, bitpix, w, h, sx, sy, b, rmax, minsig,
			       &background);
		if (r > minrad && r <= maxrad) {

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
		    rf = starRadius (image, bitpix, w, h, sx, sy, b, rmax,
				    minsig, &background);

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

    free ((char *)svec);
    return (nstars);
}


/* Check pixel at x/y for being "hot", ie, a pixel surrounded by noise.
 * If any are greater than pixel at x/y then return -1.
 * Else set the pixel at x/y to llimit and return 0.
 */

static int
HotPixel (image, bitpix, w, h, x, y, llimit)

char	*image;
int	bitpix;
int	w;
int	h;
int	x, y;
double	llimit;

{
    double pix1, pix2, pix3;

    /* Check for hot row */
    pix1 = getpix (image,bitpix,w,h,x-1,y-1);
    pix2 = getpix (image,bitpix,w,h,x,y-1);
    pix3 = getpix (image,bitpix,w,h,x+1,y-1);
    if (pix1 > llimit || pix2 > llimit || pix3 > llimit)
	return (-1);
    pix1 = getpix (image,bitpix,w,h,x-1,y+1);
    pix2 = getpix (image,bitpix,w,h,x,y+1);
    pix3 = getpix (image,bitpix,w,h,x+1,y+1);
    if (pix1 > llimit || pix2 > llimit || pix3 > llimit)
	return (-1);

    /* Check for hot column */
    pix1 = getpix (image,bitpix,w,h,x-1,y-1);
    pix2 = getpix (image,bitpix,w,h,x-1,y);
    pix3 = getpix (image,bitpix,w,h,x-1,y+1);
    if (pix1 > llimit || pix2 > llimit || pix3 > llimit)
	return (-1);
    pix1 = getpix (image,bitpix,w,h,x+1,y-1);
    pix2 = getpix (image,bitpix,w,h,x+1,y);
    pix3 = getpix (image,bitpix,w,h,x+1,y+1);
    if (pix1 > llimit || pix2 > llimit || pix3 > llimit)
	return (-1);

    /* Check for hot pixel */
    pix1 = getpix (image,bitpix,w,h,x-1,y);
    pix3 = getpix (image,bitpix,w,h,x+1,y);
    if (pix1 > llimit || pix3 > llimit)
	return (-1);
    pix1 = getpix (image,bitpix,w,h,x,y-1);
    pix3 = getpix (image,bitpix,w,h,x,y+1);
    if (pix1 > llimit || pix3 > llimit)
	return (-1);

    putpix (image, bitpix, w, h, x, y, llimit);
    return (0);
}


/* Compute and return the radius of the star centered at x0, y0.
 * A guard band is assumed to exist on the image.
 * Calling program is assumed to reject object if r > rmax.
 */

static int
starRadius (imp, bitpix, w, h, x0, y0, b, rmax, minsig, mean)

char	*imp;
int	bitpix;
int	w;
int	h;
int	x0, y0;
double	b;
double	rmax;
double	minsig;
double	*mean;

{
    int r, irmax;
    irmax = (int) rmax;

    /* Compute star's radius.
     * Scan in ever-greater circles until find one such that the peak is
     * n sigma above the mean at that radius.
     */
    for (r = 2; r <= irmax; r++) {
	int inrr = r*r;
	int outrr = (r+1)*(r+1);
	int np = 0;
	double sum = 0.0;
	int x, y;

	for (y = -r; y <= r; y++) { 
	    int yrr = y*y;
	    for (x = -r; x <= r; x++) {
		int xyrr = x*x + yrr;
		if (xyrr >= inrr && xyrr < outrr) {
		    double dp;
		    dp = getpix (imp,bitpix,w,h,x0+x,y0+y);
		    sum += dp;
		    np++;
		    }
		}
	    }

	*mean = sum / np;
	if (b > *mean + minsig)
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
 * pixel and return its location. we never take more maxrad away.
 * Return 0 if find brightest pixel within maxsteps else -1.
 */

static int dx[8]={1,0,-1,1,-1,1,0,-1};
static int dy[8]={1,1,1,0,0,-1,-1,-1};

static int
BrightWalk (image, bitpix, w, h, x0, y0, maxr, xp, yp, bp)

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
mean2d (image, bitpix, w, h, x1, x2, y1, y2, mean, sigma)

char	*image;
int	bitpix;
int	w;
int	h;
int	x1,x2;
int	y1, y2;
double	*mean;
double	*sigma;

{
    double p, pmin, pmax, pmean;
    double sd;
    int x, y;
    int i;

    pmin = -1.0e20;
    pmax = 1.0e20;

    for (i = 0; i < NITERATE; i++ ) {
	double sum = 0.0;
	double dnpix;
	int npix = 0;

    /* Compute mean */
	for (y = y1; y < y2; y++) {
	    for (x = x1; x < x2; x++) {
		p = getpix (image,bitpix,w,h, x, y);
		if (p > pmin && p < pmax) {
		    sum += p;
		    npix++;
		    }
		}
	    }
	dnpix = (double) npix;
	pmean = sum / dnpix;

    /* Compute average deviation */
	npix = 0;
	sum = 0.0;
	for (y = y1; y < y2; y++) {
	    for (x = x1; x < x2; x++) {
		p = getpix (image,bitpix,w,h, x, y);
		if (p > pmin && p < pmax) {
		    sum += fabs (p - pmean);
		    npix++;
		    }
		}
	    }

	if (npix > 0)
	    sd = sum / dnpix;
	else
	    sd = 0.0;
	pmin = pmean - sd * starsig;
	pmax = pmean + sd * starsig;
	}

    *mean = pmean;
    *sigma = sd;
    return;
}


static void
mean1d (sv1, sv2, mean, sigma)

double *sv1, *sv2;	/* starting and ending pixels for statistics */
double *mean;		/* Mean value of pixels (returned) */
double *sigma;		/* Average deviation of pixels (returned) */
{
    double *sv;
    double sumx  = 0.0;
    double sumxx = 0.0;
    int npix = 0;

    /* Compute mean */
    for (sv = sv1; sv < sv2; sv++) {
	sumx += *sv;
	npix++;
	}
    *mean = sumx / (double) npix;

    /* Compute average deviation */
    sumxx = 0.0;
    for (sv = sv1; sv < sv2; sv++) {
	sumxx += fabs (*sv - *mean);
	}
    *sigma = sumxx / (double) npix;
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
 * Aug 26 1996	Drop unused variables NH and NW
 * Aug 30 1996	Modify sigma computation; allow border to be set
 * Sep  1 1996	Set constants in lwcs.h
 * Oct 15 1996	Drop unused variables
 * Dec 10 1996	Check for hot columns as well as hot rows
 * Dec 10 1996	Add option to read image stars from DAOFIND file
 *
 * Mar 20 1997	Declare external subroutine DAOREAD
 * Nov  6 1997	Add subroutine to return image catalog filename
 * Dec 15 1997	Change calls to ABS to FABS when doubles are involved
 *
 * May 27 1998	Include imio.h
 */
