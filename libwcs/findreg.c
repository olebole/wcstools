/* this function is not for public use -- it is just used by setWCSFITS */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "wcs.h"

#define	COARSETOL	10	/* +/- this many pixels is a hit */
#define	FTOL		0.001	/* fractional change of chisqr() we call done */
#define	NPEAKS		10	/* binning peak history */

/* these statics are used by the chisqr evaluator */
static double *sx_p;
static double *sy_p;
static double *gx_p;
static double *gy_p;
static double w_p, h_p;
static int nbin_p;

static double chisqr ();
static void compute_model ();
static void call_amoeba ();

/* from Numerical Recipes */
static void amoeba();
static double amotry();
static void nrerror();
static double *vector();
static void free_vector();

/* find shift and center rotation of s's to best-match g's
 * N.B. we assume rotation will be "small enough" so that initial guesses can
 *   be done using just shifts.
 * return count of total coincidences found, else 0 if none or -1 if trouble.
 */

int
findRegistration (sx, sy, ns, gx, gy, ng, w, h, minbin, tol, rotp, dxp, dyp)

double *sx, *sy;	/* test set of x/y pairs, pixels */
int ns;			/* number of entries in sx[] and sy[] */
double *gx, *gy;	/* reference set of x/y pairs, pixels */
int ng;			/* number of entries in gx[] and gy[] */
int w, h;		/* image size, pixels */
int minbin;		/* minimum number of coincidence hits we require */
int tol;		/* +/- this many pixels is a hit */
double *rotp;		/* best-fit rotation, rads */
double *dxp, *dyp;	/* best-fit displacements */

{
	int dx, bestdx;
	int dy, bestdy;
	int bestnbin;
	int s, g;
	int nbin;
	double *sbx, *sby;	/* malloced array of s stars in best bin */
	double *gbx, *gby;	/* malloced array of g stars in best bin */
	int peaks[NPEAKS];	/* history of bin counts */
	int dxpeaks[NPEAKS], dypeaks[NPEAKS]; /* history of dx/dy at peaks */
	int npeaks;		/* entrires in use in peaks[] */

	/* do coarse alignment assuming no rotation required.
	 * this will allow us to collect a set of stars that correspond and
	 * establish an initial guess of the solution.
	 */
	npeaks = 0;
	bestnbin = 0;
	for (dx = -w/2; dx <= w/2; dx += tol-1) {
	    for (dy = -h/2; dy <= h/2; dy += tol-1) {
		nbin = 0;
		for (s = 0; s < ns; s++) {
		    int sxi = sx[s] + dx + 0.5;
		    int syi = sy[s] + dy + 0.5;
		    for (g = 0; g < ng; g++) {
			int gxi = gx[g] + 0.5;
			int gyi = gy[g] + 0.5;

			if (abs(sxi-gxi)<=tol && abs(syi-gyi)<=tol)
			    nbin++;
		    }
		}

		if (nbin >= bestnbin) {
		    int i;

		    bestnbin = nbin;
		    bestdx = dx;
		    bestdy = dy;

		    /* keep last NPEAKS bestnbins, dx and dy;
		     * put newest first in arrays
		     */
		    for (i = npeaks-1; i > 0; i--) {
			peaks[i] = peaks[i-1];
			dxpeaks[i] = dxpeaks[i-1];
			dypeaks[i] = dypeaks[i-1];
		    }
		    peaks[0] = bestnbin;
		    dxpeaks[0] = bestdx;
		    dypeaks[0] = bestdy;
		    if (npeaks < NPEAKS)
			npeaks++;
		}
	    }
	}

#undef TRACE
#ifdef TRACE
	{  int i;
	    printf ("Bin history (ns=%d ng=%d tol=%d minbin=%d):\n", ns, ng,
							    tol, minbin);
	    for (i = 0; i < npeaks; i++)
		printf (" %2d bins at dx=%3d dy=%3d\n", peaks[i], dxpeaks[i],
								dypeaks[i]);
	}
#endif /* TRACE */

	if (npeaks < 2 || peaks[1] == peaks[0]) {
	    /* peak is too broad */
#ifdef TRACE
	    printf ("  Broad peak of %d bins at dx=%d dy=%d\n", peaks[0],
								bestdx, bestdy);
#endif /* TRACE */
	}

	if (bestnbin < minbin) {
	    /* too few hits */
	    return (bestnbin);
	}

	/* one more pass to find the stars that made up the best binning */
	sbx = (double *) malloc (bestnbin * sizeof(double));
	sby = (double *) malloc (bestnbin * sizeof(double));
	gbx = (double *) malloc (bestnbin * sizeof(double));
	gby = (double *) malloc (bestnbin * sizeof(double));
	nbin = 0;
	for (s = 0; s < ns; s++) {
	    int sxi = sx[s] + bestdx + 0.5;
	    int syi = sy[s] + bestdy + 0.5;
	    for (g = 0; g < ng; g++) {
		int gxi = gx[g] + 0.5;
		int gyi = gy[g] + 0.5;

		if (abs(sxi-gxi)<=tol && abs(syi-gyi)<=tol) {
		    sbx[nbin] = sx[s];
		    sby[nbin] = sy[s];
		    gbx[nbin] = gx[g];
		    gby[nbin] = gy[g];
		    nbin++;
		}
	    }
	}
	if (nbin != bestnbin) {
	    /* this can't happen :-) */
	    fprintf (stderr, "findReg: nbin mismatch: nbin=%d bestnbin=%d!\n",
								nbin, bestnbin);
	    exit (1);
	}

	/* provide non-parametric access to the star lists */
	sx_p = sbx;
	sy_p = sby;
	gx_p = gbx;
	gy_p = gby;
	nbin_p = nbin;
	w_p = w/2;
	h_p = h/2;

	/* solve for best rotation and offsets. */
	call_amoeba (0.0, (double)bestdx, (double)bestdy, tol, rotp, dxp, dyp);

#ifdef TRACE
	printf ("Amoeba:\n");
	printf ("  initial guess: dx=%5.2f dy=%5.2f rot=%5.2f\n", 
					(double)bestdx, (double)bestdy, 0.0);
	printf ("  first soltion: dx=%5.2f dy=%5.2f rot=%5.2f\n",
						    *dxp, *dyp, raddeg(*rotp));
#endif /* TRACE */

#define RESID_REFINE
#ifdef RESID_REFINE
	/* if we have extra bins, repeat with the minbin best ones */
	if (bestnbin > minbin) {
	    double *resid = (double *) malloc (bestnbin * sizeof(double));
	    int i, j;

	    /* compute residuals at each star location */
	    for (i = 0; i < bestnbin; i++) {
		double mx, my, xe, ye;;

		compute_model (sbx[i], sby[i], *rotp, *dxp, *dyp, &mx, &my);
		xe = mx - gbx[i];
		ye = my - gby[i];
		resid[i] = xe*xe + ye*ye;
	    }

	    /* sort by increasing total residual */
	    for (i = 0; i < bestnbin-1; i++) {
		for (j = i+1; j < bestnbin; j++) {
		    if (resid[j] < resid[i]) {
			double tmp;

			tmp = sbx[i]; sbx[i] = sbx[j]; sbx[j] = tmp;
			tmp = sby[i]; sby[i] = sby[j]; sby[j] = tmp;
			tmp = gbx[i]; gbx[i] = gbx[j]; gbx[j] = tmp;
			tmp = gby[i]; gby[i] = gby[j]; gby[j] = tmp;
			tmp = resid[i]; resid[i] = resid[j]; resid[j] = tmp;
		    }
		}
	    }

	    nbin_p = minbin;
	    call_amoeba (*rotp, *dxp, *dyp, tol, rotp, dxp, dyp);

#ifdef TRACE
	printf ("  resid soltion: dx %5.2f dy %5.2f rot %5.2f\n",
						*dxp, *dyp, raddeg(*rotp));
#endif /* TRACE */

	    free ((char *)resid);
	}
#endif /* RESID_REFINE */

	free ((char *)sbx);
	free ((char *)sby);
	free ((char *)gbx);
	free ((char *)gby);

	return (bestnbin);
}

/* compute the chisqr of the vector v, where v[1]=x0, v[2]=y0, v[3]=theta.
 * (v[0] is unused)
 */

static double
chisqr (v)

double v[4];

{
	double x0 = v[1];
	double y0 = v[2];
	double th = v[3];
	double chisqr;
	int i;

	chisqr = 0.0;
	for (i = 0; i < nbin_p; i++) {
	    double xmp, ymp, dx, dy;

	    compute_model (sx_p[i], sy_p[i], th, x0, y0, &xmp, &ymp);
	    dx = xmp - gx_p[i];
	    dy = ymp - gy_p[i];
	    chisqr += dx*dx + dy*dy;
	}

#ifdef TRACE_CHSQR
	printf ("  %5.2f %5.2f %5.2f -> %5.2f\n", x0, y0, raddeg(th), chisqr);
#endif
	return (chisqr);
}


/* given a star loc at [x,y], the current model params x0/y0/th, compute the
 * resulting location.
 */

static void
compute_model (x, y, th, x0, y0, xp, yp)

double	x, y;
double	th;
double	x0, y0;
double	*xp, *yp;

{
	*xp =  x + (y - h_p)*th + x0;
	*yp = -(x - w_p)*th + y + y0;
}


/* Set up the necessary temp arrays and call the amoeba() multivariat solver */

static void
call_amoeba (rotguess, dxguess, dyguess, tol, rotp, dxp, dyp)

double rotguess;
double dxguess;
double dyguess;
int tol;
double *rotp;
double *dxp;
double *dyp;

{
    double *p[5];				  /* used as p[1..4][1..3] */
    double p0[4], p1[4], p2[4], p3[4], p4[4]; /* used as px[1..3] */
    double y[5];				  /* used as y[1..4] */
    double sum;
    int iter;
    int i;

/* set up matrix of 4 initial guesses.
 * we use the supplied guess of course plus 3 more we make up.
 */
    p[0] = p0;
    p[1] = p1;
	p1[1] = dxguess;
	p1[2] = dyguess;
	p1[3] = rotguess;
	 y[1] = chisqr (p1);
    p[2] = p2;
	p2[1] = dxguess + tol;
	p2[2] = dyguess;
	p2[3] = rotguess;
	 y[2] = chisqr (p2);
    p[3] = p3;
	p3[1] = dxguess;
	p3[2] = dyguess + tol;
	p3[3] = rotguess;
	 y[3] = chisqr (p3);
    p[4] = p4;
	p4[1] = dxguess;
	p4[2] = dyguess;
	p4[3] = rotguess + degrad(1);
	 y[4] = chisqr (p4);

    amoeba (p, y, 3, FTOL, chisqr, &iter);

#undef	PDUMP
#ifdef	PDUMP
    {
	int i;
	for (i = 1; i <= 4; i++)
	    printf ("%d: x0=%5.1f y0=%5.1f th=%5.3f y=%g\n", i, p[i][1],
					p[i][2], raddeg(p[i][3]), y[i]);
    }
#endif

    /* on return, all entries in p[1..4] are within FTOL;
     * average them?? pick first one?
     */
    for (sum = 0, i = 1; i <= 4; i++)
	sum += p[i][1];
    *dxp = sum/4;
    for (sum = 0, i = 1; i <= 4; i++)
	sum += p[i][2];
    *dyp = sum/4;
    for (sum = 0, i = 1; i <= 4; i++)
	sum += p[i][3];
    *rotp = sum/4;

#undef RESIDDUMP
#ifdef RESIDDUMP
    printf ("iter=%d rot=%g dx=%g dy=%g\n", iter, raddeg(*rotp), *dxp,*dyp);
    {
	int i;
	for (i = 0; i < nbin_p; i++) {
	    double mx, my, ex, ey;

	    compute_model (sx_p[i], sy_p[i], *rotp, *dxp, *dyp, &mx, &my);
	    ex = mx-gx_p[i];
	    ey = my-gy_p[i];

	    printf ("%2d:", i);
	    printf ("X:s=%5.1f g=%5.1f s'=%5.1f e=%4.1f ",
		    sx_p[i], gx_p[i], mx, ex);
	    printf ("Y:s=%5.1f g=%5.1f s'=%5.1f e=%4.1f ",
		    sy_p[i], gy_p[i], my, ey);
	    printf ("r=%4.1f", sqrt(ex*ex + ey*ey));
	    putchar ('\n');
	    }
    }
#endif
}

/* The following subroutines are from Numerical Recipes in C */

/* amoeba.c */

#define NMAX 5000
#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

#define GET_PSUM for (j=1;j<=ndim;j++) { for (i=1,sum=0.0;i<=mpts;i++)\
						sum += p[i][j]; psum[j]=sum;}

static void
amoeba (p, y, ndim, ftol, funk, nfunk)

double	**p;
double	y[];
double	ftol;
double	(*funk)();
int	ndim;
int	*nfunk;

{
int i,j,ilo,ihi,inhi,mpts=ndim+1;
double ytry,ysave,sum,rtol,amotry(),*psum,*vector();
void nrerror(),free_vector();

    psum=vector(1,ndim);
    *nfunk=0;
    GET_PSUM
    for (;;) {
	ilo=1;
	ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
	for (i = 1; i <= mpts; i++) {
	    if (y[i] < y[ilo])
		ilo=i;
	    if (y[i] > y[ihi]) {
		inhi=ihi;
		ihi=i;
		}
	    else if (y[i] > y[inhi])
		if (i != ihi)
		    inhi=i;
	    }
	rtol = 2.0 * fabs(y[ihi]-y[ilo]) / (fabs(y[ihi]) + fabs(y[ilo]));
	if (rtol < ftol)
	    break;
	if (*nfunk >= NMAX)
	    nrerror("Too many iterations in AMOEBA");
	ytry = amotry (p, y, psum, ndim, funk, ihi, nfunk, -ALPHA);
	if (ytry <= y[ilo])
	    ytry = amotry (p, y, psum, ndim, funk, ihi, nfunk, GAMMA);
	else if (ytry >= y[inhi]) {
	    ysave = y[ihi];
	    ytry = amotry (p,y,psum,ndim,funk,ihi,nfunk,BETA);
	    if (ytry >= ysave) {
		for (i = 1; i <= mpts; i++) {
		    if (i != ilo) {
			for (j = 1; j <= ndim; j++) {
			    psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
			    p[i][j] = psum[j];
			    }
			y[i]=(*funk)(psum);
			}
		    }
		*nfunk += ndim;
		GET_PSUM
		}
	    }
	}
    free_vector(psum,1,ndim);
}


static double
amotry (p, y, psum, ndim, funk, ihi, nfunk, fac)

double	**p;
double	*y;
double	*psum;
double	(*funk)();
double	fac;
int	ndim;
int	ihi;
int	*nfunk;

{
    int j;
    double fac1,fac2,ytry,*ptry,*vector();
    void nrerror(),free_vector();

    ptry = vector (1,ndim);
    fac1 = (1.0 - fac) / ndim;
    fac2 = fac1 - fac;
    for (j = 1; j <= ndim; j++)
	ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
    ytry = (*funk)(ptry);
    ++(*nfunk);
    if (ytry < y[ihi]) {
	y[ihi] = ytry;
	for (j = 1; j <= ndim; j++) {
    	    psum[j] +=  ptry[j] - p[ihi][j];
    	    p[ihi][j] = ptry[j];
	    }
	}
    free_vector (ptry,1,ndim);
    return ytry;
}

#undef ALPHA
#undef BETA
#undef GAMMA
#undef NMAX


/* nrutil.c */

static void
nrerror (error_text)

char	error_text[];

{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}



static double *
vector (nl,nh)

int	nl;
int	nh;

{
    double *v;

    v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
    if (!v) nrerror("allocation failure in vector()");
    return v-nl;
}


static void
free_vector (v, nl, nh)

double	*v;
int	nl;
int	nh;
{
    free((char*) (v+nl));
    return;
}
