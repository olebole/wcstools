/* File libwcs/matchstar.c
 * October 20, 1997
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

#define NPAR 7
#define NPAR1 8

#define ABS(a) ((a) < 0 ? (-(a)) : (a))

static void call_amoeba ();

/* Statics used by the chisqr evaluator */
static double	*sx_p;
static double	*sy_p;
static double	*gx_p;
static double	*gy_p;
static double	xref_p, yref_p;
static double	xrefpix, yrefpix;
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
    double vguess[NPAR];	/* Initial values for fit */
    int peaks[NPEAKS+1];	/* history of bin counts */
    int dxpeaks[NPEAKS+1], dypeaks[NPEAKS+1]; /* history of dx/dy at peaks */
    int npeaks;		/* entries in use in peaks[] */
    int maxnbin, i, nmatchd;
    int *is, *ig, *ibs, *ibg;
    char rastr[16], decstr[16];
    int minbin;		/* Minimum number of coincidence hits needed */
    int resid_refine = 0;
    int bestbin;	/* Number of coincidences for refit */
    double xinc1, yinc1;
    double mrot;

    /* Do coarse alignment assuming no rotation required.
     * This will allow us to collect a set of stars that correspond and
     * establish an initial guess of the solution.
     */
    npeaks = 0;
    nmatch = 0;
    minbin = 2;
    for (i = 0; i < NPEAKS; i++) {
	peaks[i] = 0;
	dxpeaks[i] = 0;
	dypeaks[i] = 0;
	}
    maxnbin = ns;
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
    nmatchd = nmatch * sizeof (double);
    if (!(sbx = (double *) malloc (nmatchd)))
	fprintf (stderr," Could not allocate %d bytes for SBX\n", nmatchd);
    if (!(sby = (double *) malloc (nmatchd)))
	fprintf (stderr," Could not allocate %d bytes for SBY\n", nmatchd);
    if (!(gbra = (double *) malloc (nmatchd)))
	fprintf (stderr," Could not allocate %d bytes for GBRA\n", nmatchd);
    if (!(gbdec = (double *) malloc (nmatchd)))
	fprintf (stderr," Could not allocate %d bytes for GBDEC\n", nmatchd);
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
    xrefpix = wcs->xrefpix;
    yrefpix = wcs->yrefpix;
    nbin_p = nbin;

    vguess[0] = wcs->xref;
    vguess[1] = wcs->yref;
    vguess[2] = wcs->xinc;
    vguess[3] = wcs->rot;
    vguess[4] = wcs->yinc;
    if (nfit == 6)
	vguess[5] = wcs->mrot;
    else if (nfit > 6) {
	vguess[5] = wcs->xrefpix;
	vguess[6] = wcs->yrefpix;
	}

    /* Number of parameters to fit from command line or number of matches */
    if (nfit0 > -8) {
	if (nfit0 < 0) {
	    nfit = -nfit0;
	    resid_refine = 1;
	    }
	else
	    nfit = nfit0;
	}
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
	if (nfit == 6)
	    mrot = vguess[5];
	else
	    mrot = wcs->mrot;
	fprintf (stderr,"   initial guess:\n");
	fprintf (stderr," cra= %s cdec= %s rot=%7.4f del=%7.4f,%7.4f mrot=%7.4f (%8.2f,%8.2f\n",
		rastr, decstr, vguess[3], vguess[2]*3600.0, vguess[4]*3600.0,
		mrot, xrefpix, yrefpix);
	ra2str (rastr, wcs->xref, 3);
	dec2str (decstr, wcs->yref, 2);
	fprintf (stderr,"first solution:\n");
	fprintf (stderr," cra= %s cdec= %s rot=%7.4f del=%7.4f,%7.4f mrot=%7.4f (%8.2f,%8.2f)\n", 
		 rastr, decstr, wcs->rot, 3600.0*wcs->xinc, 3600.0*wcs->yinc,
		 wcs->mrot, wcs->xrefpix, wcs->yrefpix);
	}

    /* If we have extra bins, repeat with the best ones */
    bestbin = nfit + 1;
    if (resid_refine && nmatch > bestbin) {
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
	xrefpix = wcs->xrefpix;
	yrefpix = wcs->yrefpix;
	nbin_p = bestbin;
	call_amoeba (wcs);

	if (debug) {
	    ra2str (rastr, wcs->xref, 3);
	    dec2str (decstr, wcs->yref, 2);
	    fprintf (stderr,"resid solution:\n");
	    fprintf (stderr," cra=%s cdec=%s rot=%7.4f secpix=%7.4f,%7.4f mrot=%7.4f (%8.2f,%8.2f)\n", 
	      rastr, decstr, wcs->rot, 3600.0*wcs->xinc, 3600.0*wcs->yinc,
		wcs->mrot,wcs->xrefpix, wcs->yrefpix);
	    }
	free (resid);
	}

    free (sbx);
    free (sby);
    free (gbra);
    free (gbdec);
    free (is);
    free (ig);
    free (ibs);
    free (ibg);

    return (nmatch);
}
struct WorldCoor *wcsf;

static double chisqr ();

/* From Numerical Recipes */
static void amoeba();
static double amotry();


/* Set up the necessary temp arrays and call the amoeba() multivariate solver */

static void
call_amoeba (wcs0)

struct WorldCoor *wcs0;

{
    double *p[NPAR1];				  /* used as p[NPAR1][NPAR] */
    double vguess[NPAR], vp[NPAR], vdiff[NPAR];
    double p0[NPAR], p1[NPAR], p2[NPAR], p3[NPAR], p4[NPAR],
	   p5[NPAR], p6[NPAR], p7[NPAR], p8[NPAR]; /* used as px[0..NPAR-1] */
    double y[NPAR1];				  /* used as y[1..NPAR] */
    double xinc1, yinc1, xrefpix1, yrefpix1, rot, mrot;
    int iter;
    int i, j;
    int nfit1;
    char rastr[16],decstr[16];
    int nitmax;

    nitmax = NMAX;
    if (nfit > NPAR)
	nfit = NPAR;
    nfit1 = nfit + 1;
    wcsf = wcs0;

/* Optical axis center (RA and Dec degrees)*/
    vguess[0] = 0.0;
    vguess[1] = 0.0;
    vdiff[0] = 5.0 * wcsf->xinc;
    vdiff[1] = 5.0 * wcsf->yinc;

/* Plate scale at optical axis right ascension or both (degrees/pixel) */
    vguess[2] = wcsf->xinc;
    vdiff[2] = wcsf->xinc * 0.03;

/* Rotation about optical axis in degrees */
    vguess[3] = wcsf->rot;
    vdiff[3] = 0.5;

/* Plate scale in declination at optical axis (degrees/pixel) */
    vguess[4] = wcsf->yinc;
    vdiff[4] = wcsf->yinc * 0.03;

/* Rotation about chip center (degrees) */
    if (nfit == 6) {
	vguess[5] = wcsf->mrot;
	vdiff[5] = 0.5;
	vguess[6] = 0.0;
	vdiff[6] = 0.0;
	}

/* Reference pixel (optical axis) */
    else if (nfit > 6) {
	vguess[5] = wcsf->xrefpix;
	vguess[6] = wcsf->yrefpix;
	vdiff[5] = 10.0;
	vdiff[6] = 10.0;
	}

/* Set up matrix of nfit+1 initial guesses.
 * The supplied guess, plus one for each parameter altered by a small amount
 */
    p[0] = p0;
	p0[0] = vguess[0];
	p0[1] = vguess[1];
	p0[2] = vguess[2];
	p0[3] = vguess[3];
	p0[4] = vguess[4];
	p0[5] = vguess[5];
	p0[6] = vguess[6];
	 y[0] = chisqr (p0, 0);

    /* Change optical axis right ascension */
    p[1] = p1;
	p1[0] = vguess[0] + vdiff[0];
	p1[1] = vguess[1];
	p1[2] = vguess[2];
	p1[3] = vguess[3];
	p1[4] = vguess[4];
	p1[5] = vguess[5];
	p1[6] = vguess[6];
	 y[1] = chisqr (p1, -1);

    /* Change optical axis declination */
    p[2] = p2;
	p2[0] = vguess[0];
	p2[1] = vguess[1] + vdiff[1];
	p2[2] = vguess[2];
	p2[3] = vguess[3];
	p2[4] = vguess[4];
	p2[5] = vguess[5];
	p2[6] = vguess[6];
	 y[2] = chisqr (p2, -2);

    /* Change right ascension plate scale */
    p[3] = p3;
	p3[0] = vguess[0];
	p3[1] = vguess[1];
	p3[2] = vguess[2] + vdiff[2];
	p3[3] = vguess[3];
	p3[4] = vguess[4];
	p3[5] = vguess[5];
	p3[6] = vguess[6];
	 y[3] = chisqr (p3, -3);

    /* Change rotation about optical axis */
    p[4] = p4;
	p4[0] = vguess[0];
	p4[1] = vguess[1];
	p4[2] = vguess[2];
	p4[3] = vguess[3] + vdiff[3];
	p4[4] = vguess[4];
	p4[5] = vguess[5];
	p4[6] = vguess[6];
	 y[4] = chisqr (p4, -4);

    /* Change declination plate scale */
    p[5] = p5;
	p5[0] = vguess[0];
	p5[1] = vguess[1];
	p5[2] = vguess[2];
	p5[3] = vguess[3];
	p5[4] = vguess[4] + vdiff[4];
	p5[5] = vguess[5];
	p5[6] = vguess[6];
	 y[5] = chisqr (p5, -5);

    /* Change rotation about chip center */
    if (nfit == 6) {
	p[6] = p6;
	    p6[0] = vguess[0];
	    p6[1] = vguess[1];
	    p6[2] = vguess[2];
	    p6[3] = vguess[3];
	    p6[4] = vguess[4];
	    p6[5] = vguess[5] + vdiff[5];
	    p6[6] = vguess[6];
	     y[6] = chisqr (p6, -6);
	}
    else if (nfit > 6) {
	p[6] = p6;
	    p6[0] = vguess[0];
	    p6[1] = vguess[1];
	    p6[2] = vguess[2];
	    p6[3] = vguess[3];
	    p6[4] = vguess[4];
	    p6[5] = vguess[5] + vdiff[5];
	    p6[6] = vguess[6];
	     y[6] = chisqr (p6, -6);
	p[7] = p7;
	    p7[0] = vguess[0];
	    p7[1] = vguess[1];
	    p7[2] = vguess[2];
	    p7[3] = vguess[3];
	    p7[4] = vguess[4];
	    p7[5] = vguess[5];
	    p7[6] = vguess[6] + vdiff[6];
	     y[7] = chisqr (p7, -7);
	}

#define	PDUMP
#ifdef	PDUMP
    printf ("Before:\n");
    for (i = 0; i < nfit1; i++) {
	ra2str (rastr,p[i][0]+xref_p,3);
	dec2str (decstr,p[i][1]+yref_p,2);
	xinc1 = p[i][2];
	rot = p[i][3];
	if (nfit > 4)
	    yinc1 = p[i][4];
	else
	    yinc1 = xinc1;
	if (nfit == 6)
	    mrot = p[i][5];
	else
	    mrot = wcsf->mrot;
	if (nfit > 6) {
	    xrefpix1 = xrefpix + p[i][5];
	    yrefpix1 = xrefpix + p[i][6];
	    }
	else {
	    xrefpix1 = wcsf->xrefpix;
	    yrefpix1 = wcsf->yrefpix;
	    }
	printf ("%d: %s %s rot=%5.3f del=%6.4f,%6.4f mrot=%5.3f (%8.2f,%8.2f) y=%g\n",
		i,rastr,decstr, rot, 3600.0*xinc1, 3600.0*yinc1, mrot,
		xrefpix1, yrefpix1, y[i]);
	}
#endif

    amoeba (p, y, nfit, FTOL, nitmax, chisqr, &iter);

#define	PDUMP
#ifdef	PDUMP
    printf ("\nAfter:\n");
    for (i = 0; i < nfit1; i++) {
	ra2str (rastr,p[i][0]+xref_p,3);
	dec2str (decstr,p[i][1]+yref_p,2);
	xinc1 = p[i][2];
	if (nfit > 3)
	   rot = p[i][3];
	if (nfit > 4)
	    yinc1 = p[i][4];
	else if (nfit > 2) {
	    if (xinc1 < 0)
		yinc1 = -xinc1;
	    else
		yinc1 = xinc1;
	    }
	if (nfit == 6)
	    mrot = p[i][5];
	else
	    mrot = 0.0;
	if (nfit > 6) {
	    xrefpix1 = p[i][5];
	    yrefpix1 = p[i][6];
	    }
	else {
	    xrefpix1 = wcsf->xrefpix;
	    yrefpix1 = wcsf->yrefpix;
	    }
	printf ("%d: %s %s rot=%5.3f del=%6.4f,%6.4f crot=%5.3f (%8.2f,%8.2f) y=%g\n",
		i,rastr,decstr, rot, 3600.0*xinc1, 3600.0*yinc1, mrot,
		xrefpix1, yrefpix1, y[i]);
	}
#endif

    /* on return, all entries in p[1..NPAR] are within FTOL;
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
    if (nfit == 6)
	wcsf->mrot = vp[5];
    if (nfit > 6) {
	wcsf->xrefpix = xrefpix + vp[5];
	wcsf->yrefpix = yrefpix + vp[6];
	}

#define RESIDDUMP
#ifdef RESIDDUMP
    ra2str (rastr,wcsf->xref,3);
    dec2str (decstr,wcsf->yref,2);

    printf ("iter=%d\n cra= %s cdec= %s rot=%7.4f del=%7.4f,%7.4f mrot=%7.4f (%8.2f,%8.2f)\n",
	    iter, rastr, decstr, wcsf->rot, wcsf->xinc*3600.0, wcsf->yinc*3600.0,
	    wcsf->mrot, wcsf->xrefpix, wcsf->yrefpix);
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
 * v[0]=cra, v[1]=cdec, v[2]=ra deg/pix, v[3]=theta,
 * v[4]=dec deg/pix, and v[5]=chip rotation
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
	    wcsf->yinc = -wcsf->xinc;
	else
	    wcsf->yinc = wcsf->xinc;
	}
    if (nfit == 6) {
	wcsf->mrot = v[5];
	wcsf->cmrot = cos (degrad(v[5]));
	wcsf->smrot = sin (degrad(v[5]));
	}
    if (nfit > 6) {
	wcsf->xrefpix = xrefpix + v[5];
	wcsf->yrefpix = yrefpix + v[6];
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
    fprintf (stderr,"%4d: %s %s %8.5f %9.7f,%9.7f %8.5f (%8.2f,%8.2f) -> %f\r",
	    iter, rastr, decstr, wcsf->rot, wcsf->xinc*3600.0, wcsf->yinc*3600.0,
	    wcsf->mrot, wcsf->xrefpix, wcsf->yrefpix, chsq);
#endif
    return (chsq);
}

/* The following subroutines are based on those in Numerical Recipes in C */

/* amoeba.c */

#define ALPHA 1.0
#define BETA 0.5
#define GAMMA 2.0

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
    for (j=0; j<ndim; j++) {
	for (i=0,sum=0.0; i<ndim1; i++)
	    sum += p[i][j]; psum[j]=sum;
	}
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
	    fprintf (stderr,"Too many iterations in AMOEBA %d > %d",*nfunk,itmax);
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
		for (j=0; j<ndim; j++) {
		    for (i=0,sum=0.0; i<ndim1; i++)
			sum += p[i][j]; psum[j]=sum;
		    }
		}
	    }
	}
    free (psum);
    return;
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
 *
 * Jul 21 1997	Add reference pixel position fitting
 * Aug  4 1997	Increase maximum iterations from 750 to 1000 in lwcs.h
 * Aug 28 1997	Fix VGUESS dimension bug
 * Sep  9 1997	Print RA and Dec offsets in residual listing
 * Sep  9 1997	Turn on resid_refinement if number of parameters to fit negated
 * Sep  9 1997	Fit separate horizontal and vertical plate scales if nfit=5
 * Sep  9 1997	Fix bugs associated with fitting optical axis
 * Sep 12 1997	Add chip rotation instead of second plate scale
 * Oct  2 1997	Keep second plate scale AND chip rotation
 * Oct 16 1997	Try to deal with reference pixel position correctly
 */ 
