/* File libwcs/matchstar.c
 * November 19, 1996
 * By Doug Mink, Smithsonian Astrophyscial Observatory
 */

/* StarMatch (ns, sx, sy, ng, gra, gdec, gx, gy, tol, wcs, nfit, debug)
 *  Find shift, scale, and rotation of image stars to best-match reference stars
 * call_amoeba (wcs0) Set up temp arrays and call multivariate solver
 * chisqr (v) Compute the chisqr of the vector v
 * amoeba (p, y, ndim, ftol, itmax, funk, nfunk)
 *    Multivariate solver from Numerical Recipes
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "fitshead.h"
#include "wcs.h"
#include "lwcs.h"

#define ABS(a) ((a) < 0 ? (-(a)) : (a))

#undef RESID_REFINE

static void call_amoeba ();

/* Statics used by the chisqr evaluator */
static double	*sx_p;
static double	*sy_p;
static double	*gx_p;
static double	*gy_p;
static double	xref_p, yref_p;
static int	nbin_p;
static int	nfit;	/* Number of parameters to fit */

/* Find shift, scale, and rotation of image stars to best-match reference stars
 * Get best match by finding which offsets between pairs of s's and g's
 * work for the most other pairs of s's and g's
 * N.B. we assume rotation will be "small enough" so that initial guesses can
 *   be done using just shifts.
 * Return count of total coincidences found, else 0 if none or -1 if trouble.
 */

int
StarMatch (ns, sx, sy, ng, gra, gdec, gx, gy, tol, wcs, nfit0, debug)

int	ns;		/* Number of image stars */
double	*sx;		/* Image star X coordinates in pixels */
double	*sy;		/* Image star Y coordinates in pixels */
int	ng;		/* Number of reference stars */
double	*gra;		/* Reference star right ascensions in degrees */
double	*gdec;		/* Reference star right ascensions in degrees */
double	*gx;		/* Reference star X coordinates in pixels */
double	*gy;		/* Reference star Y coordinates in pixels */
double	tol;		/* +/- this many pixels is a hit */
struct WorldCoor *wcs;	/* World coordinate structure (fit returned) */
int	nfit0;		/* Number of parameters to fit (0=set from matches) */
int	debug;

{
    double dx, bestdx, dxi;
    double dy, bestdy, dyi;
    int nmatch;
    int s, g, si, gi;
    int nbin;
    double *sbx, *sby;	/* malloced array of s stars in best bin */
    double *gbra, *gbdec;	/* malloced array of g stars in best bin */
    double vguess[5];	/* Initial values for fit */
    int peaks[NPEAKS];	/* history of bin counts */
    int dxpeaks[NPEAKS], dypeaks[NPEAKS]; /* history of dx/dy at peaks */
    int npeaks;		/* entries in use in peaks[] */
    int maxnbin, i;
    int *is, *ig, *ibs, *ibg;
    char rastr[16], decstr[16];
    int minbin;		/* Minimum number of coincidence hits needed */

#ifdef RESID_REFINE
    int bestbin;	/* Number of coincidences for refit */
#endif /* RESID_REFINE */

    /* Do coarse alignment assuming no rotation required.
     * This will allow us to collect a set of stars that correspond and
     * establish an initial guess of the solution.
     */
    npeaks = 0;
    nmatch = 0;
    maxnbin = ns;
    minbin = 2;
    for (i = 0; i < NPEAKS; i++) {
	peaks[i] = 0;
	dxpeaks[i] = 0;
	dypeaks[i] = 0;
	}
    if (ng > ns)
	maxnbin = ng;
    is = (int *) malloc (maxnbin * sizeof(int));
    ig = (int *) malloc (maxnbin * sizeof(int));
    ibs = (int *) malloc (maxnbin * sizeof(int));
    ibg = (int *) malloc (maxnbin * sizeof(int));
    for (g = 0; g < ng; g++) {
	for (s = 0; s < ns; s++) {
	    dx = gx[g] - sx[s];
	    dy = gy[g] - sy[s];
	    nbin = 1;
	    is[0] = s;
	    ig[0] = g;
	    for (gi = 0; gi < ng; gi++) {
		for (si = 0; si < ns; si++) {
		    if (si != s && gi != g) {
			dxi = ABS (gx[gi] - sx[si] - dx);
			dyi = ABS (gy[gi] - sy[si] - dy);
			if (dxi <= tol && dyi <= tol) {
			    /* if (debug)
				printf ("%d %d %d %d %5.1f %5.1f %5.1f %5.1f\n",
					g,s,gi,si,dx,dy,dxi,dyi); */
			    is[nbin] = si;
			    ig[nbin] = gi;
			    nbin++;
			    }
			}
		    }
		}
	    /* if (debug)
		printf ("%d %d %d %d %d\n", g,s,gi,si,nbin); */
	    if (nbin > 1 && nbin >= nmatch) {
		int i;
		nmatch = nbin;
		bestdx = (double) dx;
		bestdy = (double) dy;
		for (i = 0; i < nbin; i++) {
		    ibs[i] = is[i];
		    ibg[i] = ig[i];
		    }
	
		/* keep last NPEAKS nmatchs, dx and dy;
		 * put newest first in arrays */
		if (npeaks > 0) {
		    for (i = npeaks; i > 0; i--) {
			peaks[i] = peaks[i-1];
			dxpeaks[i] = dxpeaks[i-1];
			dypeaks[i] = dypeaks[i-1];
			}
		    }
		peaks[0] = nmatch;
		dxpeaks[0] = (int) (bestdx + 0.5);
		dypeaks[0] = (int) (bestdy + 0.5);
		if (npeaks < NPEAKS)
		    npeaks++;
		}
	    }
	}

    if (debug) {
	int i;
	fprintf (stderr,"Bin history (ns=%d ng=%d tol=%3.0f minbin=%d):\n",
		ns, ng, tol, minbin);
	for (i = 0; i < npeaks; i++)
	    fprintf (stderr," %d bins at dx=%d dy=%d\n",
		    peaks[i], dxpeaks[i], dypeaks[i]);
	}

    /* peak is broad */
    if (npeaks < 2 || peaks[1] == peaks[0]) {
	if (debug)
	    fprintf (stderr,"  Broad peak of %d bins at dx=%.0f dy=%.0f\n",
		     peaks[0], bestdx, bestdy);
	}

    /* too few hits */
    if (nmatch < minbin)
	return (nmatch);

    /* Get X and Y coordinates of matches from best binning */
    sbx = (double *) malloc (nmatch * sizeof(double));
    sby = (double *) malloc (nmatch * sizeof(double));
    gbra = (double *) malloc (nmatch * sizeof(double));
    gbdec = (double *) malloc (nmatch * sizeof(double));
    for (nbin = 0; nbin < nmatch; nbin++) {
	sbx[nbin] = sx[ibs[nbin]];
	sby[nbin] = sy[ibs[nbin]];
	gbra[nbin] = gra[ibg[nbin]];
	gbdec[nbin] = gdec[ibg[nbin]];
	}

    /* Reset image center based on star matching */
    wcs->xref = wcs->xref + (dx * wcs->xinc);
    wcs->yref = wcs->yref + (dy * wcs->yinc);

    /* Provide non-parametric access to the star lists */
    sx_p = sbx;
    sy_p = sby;
    gx_p = gbra;
    gy_p = gbdec;
    xref_p = wcs->xref;
    yref_p = wcs->yref;
    nbin_p = nbin;

    vguess[0] = wcs->xref;
    vguess[1] = wcs->yref;
    vguess[2] = wcs->xinc;
    vguess[3] = wcs->rot;
    vguess[4] = wcs->yinc;

    if (nfit0 > 0)
	nfit = nfit0;
    else if (nbin < 4)
	nfit = 2;
    else if (nbin > 5)
	nfit = 4;
    else
	nfit = 3;

    /* Fit image star coordinates to reference star positions */
    call_amoeba (wcs);

    if (debug) {
	fprintf (stderr,"Amoeba:\n");
	ra2str (rastr, vguess[0], 3);
	dec2str (decstr, vguess[1], 2);
	fprintf (stderr,"   initial guess:\n");
	fprintf (stderr," cra= %s cdec= %s rot=%7.4f del=%7.4f %7.4f\n", 
		rastr, decstr,vguess[3],vguess[2]*3600.0, vguess[4]*3600.0);
	ra2str (rastr, wcs->xref, 3);
	dec2str (decstr, wcs->yref, 2);
	fprintf (stderr,"first solution:\n");
	fprintf (stderr," cra= %s cdec= %s rot=%7.4f del=%7.4f %7.4f\n", 
		 rastr, decstr, wcs->rot, 3600.0*wcs->xinc, 3600.0*wcs->yinc);
	}

#ifdef RESID_REFINE
    /* If we have extra bins, repeat with the best ones */
    bestbin = nfit + 1;
    if (nmatch > bestbin) {
	double *resid = (double *) malloc (nmatch * sizeof(double));
	int i, j;

	/* Compute residuals at each star location */
	for (i = 0; i < nmatch; i++) {
	    double mx, my, xe, ye;;

	    pix2wcs (wcs, sbx[i], sby[i], &mx, &my);
	    xe = mx - gbra[i];
	    ye = my - gbdec[i];
	    resid[i] = xe*xe + ye*ye;
	    }

	/* sort by increasing total residual */
	for (i = 0; i < nmatch-1; i++) {
	    for (j = i+1; j < nmatch; j++) {
		if (resid[j] < resid[i]) {
		    double tmp;

		    tmp = sbx[i]; sbx[i] = sbx[j]; sbx[j] = tmp;
		    tmp = sby[i]; sby[i] = sby[j]; sby[j] = tmp;
		    tmp = gbra[i]; gbra[i] = gbra[j]; gbra[j] = tmp;
		    tmp = gbdec[i]; gbdec[i] = gbdec[j]; gbdec[j] = tmp;
		    tmp = resid[i]; resid[i] = resid[j]; resid[j] = tmp;
		    }
		}
	    }

	xref_p = wcs->xref;
	yref_p = wcs->yref;
	nbin_p = bestbin;
	call_amoeba (wcs);

	if (debug)
	    ra2str (rastr, wcs->xref, 3);
	    dec2str (decstr, wcs->yref, 2);
	    fprintf (stderr,"resid solution:\n");
	    fprintf (stderr," cra=%s cdec=%s rot=%7.4f secpix=%7.4f %7.4f\n", 
		 rastr, decstr, wcs->rot, 3600.0*wcs->xinc, 3600.0*wcs->yinc);

	free ((char *)resid);
	}
#endif /* RESID_REFINE */

    free ((char *)sbx);
    free ((char *)sby);
    free ((char *)gbra);
    free ((char *)gbdec);
    free ((char *)is);
    free ((char *)ig);
    free ((char *)ibs);
    free ((char *)ibg);

    return (nmatch);
}
struct WorldCoor *wcsf;

static double chisqr ();

/* From Numerical Recipes */
static void amoeba();
static double amotry();

#define NPAR 5
#define NPAR1 6

/* Set up the necessary temp arrays and call the amoeba() multivariate solver */

static void
call_amoeba (wcs0)

struct WorldCoor *wcs0;

{
    double vguess[NPAR], vp[NPAR];
    double *p[NPAR+1];				  /* used as p[NPAR1][NPAR] */
    double p0[NPAR], p1[NPAR], p2[NPAR], p3[NPAR], p4[NPAR],
	   p5[NPAR]; /* used as px[0..NPAR-1] */
    double y[NPAR1];				  /* used as y[1..NPAR] */
    int iter;
    int i, j;
    int nfit1;
    char rastr[16],decstr[16];
    int nitmax;

    nitmax = NMAX;
    nfit1 = nfit + 1;
    wcsf = wcs0;
    vguess[0] = 0.0;
    vguess[1] = 0.0;
    vguess[2] = wcsf->xinc;
    vguess[3] = wcsf->rot;
    if (nfit > 2 && nfit < 5) {
	if (wcsf->xinc < 0)
	    vguess[4] = -wcsf->xinc;
	else
	    vguess[4] = wcsf->xinc;
	}
    else
	vguess[4] = wcsf->yinc;

/* Set up matrix of 6 initial guesses.
 * The supplied guess, plus one for each parameter alltered by a small amount
 */
    p[0] = p0;
	p0[0] = vguess[0];
	p0[1] = vguess[1];
	p0[2] = vguess[2];
	p0[3] = vguess[3];
	p0[4] = vguess[4];
	 y[0] = chisqr (p0, 0);
    p[1] = p1;
	p1[0] = vguess[0] + (5 * wcsf->xinc);
	p1[1] = vguess[1];
	p1[2] = vguess[2];
	p1[3] = vguess[3];
	p1[4] = vguess[4];
	 y[1] = chisqr (p1, 0);
    p[2] = p2;
	p2[0] = vguess[0];
	p2[1] = vguess[1] + (5 * wcsf->yinc);
	p2[2] = vguess[2];
	p2[3] = vguess[3];
	p2[4] = vguess[4];
	 y[2] = chisqr (p2, 0);
    p[3] = p3;
	p3[0] = vguess[0];
	p3[1] = vguess[1];
	p3[2] = vguess[2] * 1.01;
	p3[3] = vguess[3];
	p3[4] = vguess[4];
	 y[3] = chisqr (p3, 0);
    p[4] = p4;
	p4[0] = vguess[0];
	p4[1] = vguess[1];
	p4[2] = vguess[2];
	p4[3] = vguess[3] + 1.0;
	p4[4] = vguess[4];
	 y[4] = chisqr (p4, 0);
    p[5] = p5;
	p5[0] = vguess[0];
	p5[1] = vguess[1];
	p5[2] = vguess[2];
	p5[3] = vguess[3];
	p5[4] = vguess[4] * 1.01;
	 y[5] = chisqr (p5, 0);

    amoeba (p, y, nfit, FTOL, nitmax, chisqr, &iter);

#define	PDUMP
#ifdef	PDUMP
    for (i = 0; i < nfit1; i++) {
	ra2str (rastr,p[i][0]+xref_p,3);
	dec2str (decstr,p[i][1]+yref_p,2);
	if (nfit > 2 && nfit < 5)
	    p[i][4] = p[i][2];
	printf ("%d: %s %s rot=%5.3f del=%6.4f %6.4f y=%g\n",
		i,rastr,decstr,p[i][3], 3600.0*p[i][2], 3600.0*p[i][4], y[i]);
	}
#endif

    /* on return, all entries in p[1..5] are within FTOL;
     * average them?? pick first one?
     */
    for (j = 0; j < NPAR; j++) {
	double sum = 0.0;
        for (i = 0; i < nfit1; i++)
	    sum += p[i][j];
	vp[j] = sum / (double)nfit1;
	}
    wcsf->xref = xref_p + vp[0];
    wcsf->yref = yref_p + vp[1];
    if (nfit > 2)
	wcsf->xinc = vp[2];
    if (nfit > 3) {
	wcsf->rot = vp[3];
	wcsf->srot = sin (degrad (vp[3]));
	wcsf->crot = cos (degrad (vp[3]));
	}
    if (nfit > 4)
	wcsf->yinc = vp[4];
    else if (nfit > 2) {
	if (wcsf->xinc < 0)
	    wcsf->yinc = -vp[2];
	else
	    wcsf->yinc = vp[2];
	}

#define RESIDDUMP
#ifdef RESIDDUMP
    ra2str (rastr,wcsf->xref,3);
    dec2str (decstr,wcsf->yref,2);

    printf ("iter=%d\n cra= %s cdec= %s rot=%7.4f del=%7.4f %7.4f\n",
	    iter, rastr, decstr, vp[3], vp[2]*3600.0, vp[4]*3600.0);
    for (i = 0; i < nbin_p; i++) {
	double mx, my, ex, ey;
	char gstr[16], sstr[16];

	pix2wcs (wcsf, sx_p[i], sy_p[i], &mx, &my);
	ex = mx - gx_p[i];
	ey = my - gy_p[i];

	ra2str (gstr, gx_p[i], 3);
	ra2str (sstr, mx, 3);
	printf ("%2d: X: g %s s' %s e: %6.4f\n",
		i+1, gstr, sstr, 3600.0*ex);
	dec2str (gstr, gy_p[i], 2);
	dec2str (sstr, my, 2);
	printf ("    Y: g: %s s': %s e: %6.2f ",
		gstr, sstr, 3600.0*ey);
	printf ("r=%6.2f", 3600.0*sqrt(ex*ex + ey*ey));
	putchar ('\n');
	}
#endif
}


/* Compute the chisqr of the vector v, where
 * v[0]=cra, v[1]=cdec, v[3]=theta, v[2]=ra deg/pix, v[4]=dec deg/pix
 * chisqr is in arcsec^2
 */

static double
chisqr (v, iter)

double	*v;	/* Vector of parameter values */
int	iter;	/* Number of iterations */

{
    double chsq;
    char rastr[16],decstr[16];
    double xmp, ymp, dx, dy;
    int i, offscale;

    /* Set WCS parameters from fit parameter vector */
    wcsf->xref = xref_p + v[0];
    wcsf->yref = yref_p + v[1];
    if (nfit > 2)
	wcsf->xinc = v[2];
    if (nfit > 3) {
	wcsf->rot = v[3];
	wcsf->crot = cos (degrad(v[3]));
	wcsf->srot = sin (degrad(v[3]));
	}
    if (nfit > 4)
	wcsf->yinc = v[4];
    else if (nfit > 2) {
	if (wcsf->xinc < 0)
	    wcsf->yinc = -v[2];
	else
	    wcsf->yinc = v[2];
	}

    /* Compute sum of squared residuals for these parameters */
    chsq = 0.0;
    for (i = 0; i < nbin_p; i++) {
	wcs2pix (wcsf, gx_p[i], gy_p[i], &xmp, &ymp, &offscale);
	/* if (!offscale) { */
	    dx = xmp - sx_p[i];
	    dy = ymp - sy_p[i];
	    chsq += dx*dx + dy*dy;
	    /* } */
	}

#define TRACE_CHSQR
#ifdef TRACE_CHSQR
    ra2str (rastr,wcsf->xref,3);
    dec2str (decstr,wcsf->yref,2);
    fprintf (stderr,"%4d: %s %s %8.5f %9.7f %9.7f -> %f",
	    iter, rastr, decstr, wcsf->rot,
	    wcsf->xinc*3600.0, wcsf->yinc*3600.0, chsq);
    (void)putc (13,stderr);
#endif
    return (chsq);
}

/* The following subroutines are based on those in Numerical Recipes in C */

/* amoeba.c */

#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

#define GET_PSUM for (j=0; j<ndim; j++) { for (i=0,sum=0.0; i<ndim1; i++)\
						sum += p[i][j]; psum[j]=sum;}

static void
amoeba (p, y, ndim, ftol, itmax, funk, nfunk)

double	**p;
double	y[];
double	ftol;
int	itmax;
double	(*funk)();
int	ndim;
int	*nfunk;

{
int i,j,ilo,ihi,inhi,ndim1=ndim+1;
double ytry,ysave,sum,rtol,*psum;

    psum = (double *) malloc ((unsigned)ndim * sizeof(double));
    *nfunk = 0;
    GET_PSUM
    for (;;) {
	ilo=1;
	if (y[0] > y[1]) {
	    inhi = 1;
	    ihi = 0;
	    }
	else {
	    inhi = 0;
	    ihi = 1;
	    }
	for (i = 0; i < ndim1; i++) {
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
	if (*nfunk >= itmax) {
	    fprintf (stderr,"Numerical Recipes run-time error...\n");
	    fprintf (stderr,"%s\n","Too many iterations in AMOEBA");
	    fprintf (stderr,"...now exiting to system...\n");
	    exit (1);
	    }
	ytry = amotry (p, y, psum, ndim, funk, ihi, nfunk, -ALPHA);
	if (ytry <= y[ilo])
	    ytry = amotry (p, y, psum, ndim, funk, ihi, nfunk, GAMMA);
	else if (ytry >= y[inhi]) {
	    ysave = y[ihi];
	    ytry = amotry (p,y,psum,ndim,funk,ihi,nfunk,BETA);
	    if (ytry >= ysave) {
		for (i = 0; i < ndim1; i++) {
		    if (i != ilo) {
			for (j = 0; j < ndim; j++) {
			    psum[j] = 0.5 * (p[i][j] + p[ilo][j]);
			    p[i][j] = psum[j];
			    }
			y[i]=(*funk)(psum, *nfunk);
			}
		    }
		*nfunk += ndim;
		GET_PSUM
		}
	    }
	}
    free (psum);
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
    double fac1,fac2,ytry,*ptry;

    ptry = (double *) malloc ((unsigned) ndim * sizeof(double));
    fac1 = (1.0 - fac) / ndim;
    fac2 = fac1 - fac;
    for (j = 0; j < ndim; j++)
	ptry[j] = psum[j] * fac1 - p[ihi][j] * fac2;
    ytry = (*funk)(ptry, *nfunk);
    ++(*nfunk);
    if (ytry < y[ihi]) {
	y[ihi] = ytry;
	for (j = 0; j < ndim; j++) {
    	    psum[j] +=  ptry[j] - p[ihi][j];
    	    p[ihi][j] = ptry[j];
	    }
	}
    free (ptry);
    return ytry;
}
/* Aug  6 1996	New subroutine
 * Sep  1 1996	Move constants to lwcs.h
 * Sep  3 1996	Use offscale pixels for chi^2 computation
 * Sep  3 1996	Overprint chi^2 in verbose mode
 * Oct 15 1996	Fix am* subroutine declarations
 * Nov 19 1996	Fix bug regarding rotation
 */ 
