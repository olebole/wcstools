/*** File libwcs/tnxfit.c
 *** January 11, 2001
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 */

/*  Nonlinear least squares fitting program using data array 
 *  (x, y) to fit array (x1, y1)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wcs.h"
#include "lwcs.h"

static void tnx_amoeba();
static double tnx_chisqr();
static int ncoeff=0;
static double   *sx_p;
static double   *sy_p;
static double   *gx_p;
static double   *gy_p;
static int	nbin_p;

#define MAXPAR 26
#define MAXPAR1 27
#define NITMAX 2500

int
TNXFit (wcs, x, y, x1, y1, np, ncoeff0, debug)

struct WorldCoor *wcs;	/* World coordinate system structure */
double	*x, *y;		/* Image WCS coordinates */
double	*x1, *y1;	/* Image pixel coordinates */
int	np;		/* Number of points to fit */
int	ncoeff0;	/* Order of polynomial terms in x and y */
int	debug;

{
    sx_p = x;
    sy_p = y;
    gx_p = x1;
    gy_p = y1;
    nbin_p = np;
    ncoeff = ncoeff0;

    /* Fit polynomials */
    tnx_amoeba (wcs);

    return (0);
}

static struct WorldCoor *wcsp;

/* Set up the necessary temp arrays and call the amoeba() multivariate solver */

static void
tnx_amoeba (wcs0)

struct WorldCoor *wcs0;

{
    double *p[MAXPAR1];				  /* used as p[NPAR1][NPAR] */
    double vguess[MAXPAR], vp[MAXPAR], vdiff[MAXPAR];
    double y[MAXPAR1];				  /* used as y[1..NPAR] */
    double sumx, sumy, sumr;
    int nbytes;
    int iter;
    int i, j;
    int nfit, nfit1;
    int nitmax;
    extern void amoeba();

    /* Allocate memory for fit */
    nfit = ncoeff * 2;
    nfit1 = nfit + 1;
    nbytes = nfit * sizeof (double);
    for (i = 0; i < nfit1; i++)
	p[i] = (double *) malloc (nbytes);

    nitmax = NITMAX;
    wcsp = wcs0;

/* Zero guess and difference vectors */
    for (i = 0; i < MAXPAR; i++) {
	vguess[i] = 0.0;
	vdiff[i] = 0.0;
	}

    if (nfit > 0) {
	vguess[0] = 1.0;
	vdiff[0] = 0.01;
	vguess[1] = 0.0;
	vdiff[1] = 0.01;
	vguess[2] = 0.0;
	vdiff[2] = 0.01;
	vguess[3] = 0.0;
	vdiff[3] = 0.01;
	vguess[4] = 0.0;
	vdiff[4] = 0.01;
	vguess[5] = 0.0;
	vdiff[5] = 0.01;
	vguess[6] = 1.0;
	vdiff[6] = 0.01;
	vguess[7] = 0.0;
	vdiff[7] = 0.01;
	vguess[8] = 0.0;
	vdiff[8] = 0.01;
	vguess[9] = 0.0;
	vdiff[9] = 0.01;
	vguess[10] = 0.0;
	vdiff[10] = 0.01;
	vguess[11] = 0.0;
	vdiff[11] = 0.01;
	}

    /* Set up matrix of nfit+1 initial guesses.
     * The supplied guess, plus one for each parameter altered by a small amount
     */
    for (i = 0; i < nfit1; i++) {
	for (j = 0; j < nfit; j++)
	    p[i][j] = vguess[j];
	if (i > 0)
	    p[i][i-1] = vguess[i-1] + vdiff[i-1];
	y[i] = tnx_chisqr (p[i], -i);
	}

#define	PDUMP
#ifdef	PDUMP
    fprintf (stderr,"Before:\n");
    for (i = 0; i < nfit1; i++) {
	fprintf (stderr,"%3d: ", i);
	for (j = 0; j < ncoeff; j++)
	    fprintf (stderr," %9.7f",p[i][j]);
	fprintf (stderr,"\n     ");
	for (j = 0; j < ncoeff; j++)
	    fprintf (stderr," %9.7f",p[i][ncoeff+j]);
	fprintf (stderr,"\n");
	}
#endif

    amoeba (p, y, nfit, FTOL, nitmax, tnx_chisqr, &iter);

#define	PDUMP
#ifdef	PDUMP
    fprintf (stderr,"\nAfter:\n");
    for (i = 0; i < nfit1; i++) {
	fprintf (stderr,"%3d: ", i);
	for (j = 0; j < ncoeff; j++)
	    fprintf (stderr," %9.7f",p[i][j]);
	fprintf (stderr,"\n     ");
	for (j = 0; j < ncoeff; j++)
	    fprintf (stderr," %9.7f",p[i][ncoeff+j]);
	fprintf (stderr,"\n");
	}
#endif

    /* on return, all entries in p[1..NPAR] are within FTOL; average them */
    for (j = 0; j < MAXPAR; j++) {
	double sum = 0.0;
        for (i = 0; i < nfit1; i++)
	    sum += p[i][j];
	vp[j] = sum / (double)nfit1;
	}
    tnxset (wcsp, nfit, vp);

#define RESIDDUMP
#ifdef RESIDDUMP
    fprintf (stderr,"iter=%4d\n  ", iter);
    for (j = 0; j < ncoeff; j++)
	    fprintf (stderr," %9.7f",vp[j]);
    fprintf (stderr,"\n    ");
    for (j = 0; j < ncoeff; j++)
	    fprintf (stderr," %9.7f",vp[j+6]);
    fprintf (stderr,"\n");

    sumx = 0.0;
    sumy = 0.0;
    sumr = 0.0;
    for (i = 0; i < nbin_p; i++) {
	double mx, my, ex, ey, er;
	char rastr[16], decstr[16];

	pix2wcs (wcsp, sx_p[i], sy_p[i], &mx, &my);
	ex = 3600.0 * (mx - gx_p[i]);
	ey = 3600.0 * (my - gy_p[i]);
	er = sqrt (ex * ex + ey * ey);
	sumx = sumx + ex;
	sumy = sumy + ey;
	sumr = sumr + er;

	ra2str (rastr, 16, gx_p[i], 3);
	dec2str (decstr, 16, gy_p[i], 2);
	ra2str (rastr, 16, mx, 3);
	dec2str (decstr, 16, my, 2);
	fprintf (stderr,"%2d: c: %s %s ", i+1, rastr, decstr);
	fprintf (stderr,"i: %s %s %6.3f %6.3f %6.3f\n",
		rastr, decstr, 3600.0*ex, 3600.0*ey,
		3600.0*sqrt(ex*ex + ey*ey));
	}
    sumx = sumx / (double)nbin_p;
    sumy = sumy / (double)nbin_p;
    sumr = sumr / (double)nbin_p;
    fprintf (stderr,"mean dra: %6.3f, ddec: %6.3f, dr = %6.3f\n", sumx, sumy, sumr);
#endif

    for (i = 0; i < nfit1; i++)
	free (p[i]);
    return;
}


/* Compute the chisqr of the vector v, where v[i]=plate fit coeffients
 * chisqr is in arcsec^2
 */

static double
tnx_chisqr (v, iter)

double	*v;	/* Vector of parameter values */
int	iter;	/* Number of iterations */

{
    double chsq;
    double xmp, ymp, dx, dy;
    int i, j, offscale;
    extern int tnxpset();

    /* Set plate constants from fit parameter vector */
    if (tnxpset (wcsp, 3, 3, 2, v)) {
	fprintf (stderr,"CHISQR: Cannot reset WCS!\n");
	return (0.0);
	}

    /* Compute sum of squared residuals for these parameters */
    chsq = 0.0;
    for (i = 0; i < nbin_p; i++) {
	wcs2pix (wcsp, gx_p[i], gy_p[i], &xmp, &ymp, &offscale);
	/* if (!offscale) { */
	    dx = xmp - sx_p[i];
	    dy = ymp - sy_p[i];
	    chsq += dx*dx + dy*dy;
	    /* } */
	}

#define TRACE_CHSQR
#ifdef TRACE_CHSQR
    fprintf (stderr,"%4d:", iter);
    for (j = 0; j < ncoeff; j++)
	fprintf (stderr," %9.7f",v[j]);
    for (j = 0; j < ncoeff; j++)
	fprintf (stderr," %9.7f",v[ncoeff+j]);
    fprintf (stderr," -> %f\n", chsq);
#endif
    return (chsq);
}

/* Mar 26 1998	New subroutines
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Oct 15 1999	Include stdlib.h for malloc() declaration
 *
 * Jan 11 2001	Print all messages to stderr
 */
