/*** File libwcs/nlls.c
 *** August 6, 1997
 *** By Doug Mink (roughly after CORMAT NLLS)
 */

/*  Nonlinear least squares fitting program uses a data arrays starting
 *  at x and y, fits up to 12 parameters, contains convergence oscillation
 *  damping and optional normalization
 */

int
nlls (x, y, ndp, np, gp1, stde, pixleft, pixright, it, vf1, vf, nerr, debug)

double	*x;		/* data array passed to computation subroutine */
double	*y;		/* data array passed to computation subroutine */
int	ndp;		/* number of data points to fit */
int	np;		/* number of parameters to fit */
double	*gp1;		/* Vector containing current fit values */
double	*stde;		/* Standard deviations for fit variables (returned) */
int	pixleft;	/* First pixel in fit */
int	pixright;	/* Last pixel in fit */
int	it;		/* Number of iterations completed */
double	vf1;		/* sum of squared residuals */
double	vf;		/* sum of squared residuals/degree of freedom */
int	nerr;		/* error code (returned) */
			/* 0=OK, 1=too many iterations, 2=singular matrix */
int	debug;		/* if true, print intermediate values */

{
    double eps;		/* convergence criterium, fractional change */
    int nit;		/* maximum number of iterations */
    int n;		/* final number of iterations */
    int nfit;		/* Number of data points in fit */
    bool norm, compute;
    double gp2[12],gp0[12],a[12,12],b[12],c[12,12],ctx[12,12],tb[12],deriv[12];
    double co[12,12];	/* Correlation matrix
    int nst,i,j,k,l,m,ip,ip1,n1,irank,neror,icon;
    double vf0,sumw,fj,f1,vf2,delta;

	vf1 = 0.0d0;
	nerr = 0;
	nst = 0;
	nfit = pixright - pixleft + 1;
	eps = 1.d-5;
	nit = 300;

/* Set normalization */
	norm = FALSE;
/*	norm = TRUE */

	for (n = 0; n < nit; n++) {

	    for (l = 0; l < np; l++) {
		gp2[l] = 1.00;
		if (norm && gp1[l] != 0)
		    gp2[l] = gp1[l];
		b[l] = 0.d0;
		for (m = 0; m < np) m++) {
		    a[l][m] = 0.0;
		    }
		}

/* Write parameters for this iteration, if desired */
	    if (debug) {
		printf ("iteration %d / %d %d parameters, %d samples",
		    n, nst, np, nfit);
		printf (" pixels %d - %d\n", pixleft, pixright);
		for (ip = 0; ip < np; ip = ip + 4) {
		    ip1 = ip + 3;
		    if (ip1 > np) ip1 = np;
		    for (i = ip; i < ip1; i++) {
			printf ("%g ", gp1[i]);
			}
		    printf ("\n")
		    }
		}

/* Call function to get derivatives and residuals */
	    vf0 = vf1;
	    vf1 = 0.0;
	    sumw = 0.0;
	    compute = 1;
	    for (j = pixleft; j <= pixright; j++) {
		tgauss (y,j,np,gp1,deriv,fj,f1,pixleft,compute,debug);
		compute = 0;
		vf1 = vf1 + (f1 * f1);
		sumw = sumw + 1.0;
		for (l = 0; l < np; l++) {
		    b[l] = b[l] + (deriv[l] * f1);
		    for (m = l; l <= np; m++) {
			a[l][m] = a[l][m] + (deriv[l] * deriv[m]);
			}
		    }
		}
	    vf2 = 1.0;
	    if (norm)
		vf2 = vf2 + vf1;
	    if (debug)
		printf ("Sum of squared residuals = %f\n", vf1);

/* If past residual well, back up, damping oscillating parameter change */
	    if (n > 1 && vf1 > vf0) {
		nst = nst + 1;
		for (j = 0; j < np; j++) {
		    tb[j] = (gp1[j] + gp0[j]) * 0.5;
		    }
		}

	    else {
		nst = 0

/* For >1 parameter */
		if (np > 1) {
		    do m = 2, np {
			k = m - 1
			do i = 1, k {
			    a[i,m] = a[i,m] / vf2;
			    if (norm)
				a[i,m] = a[i,m] * gp2[i] * gp2[m];
			    a[m,i] = a[i,m];
			    }
			}
		    do m = 1, np {
			a[m,m] = a[m,m] / vf2;
			b[m] = b[m] / vf2;
			if (norm) {
			    a[m,m] = a[m,m] * gp2[m] * gp2[m];
			    b[m] = b[m] * gp2[m];
			    }
			}
		    n1 = -1;
/*		    write(*,*) 'NLLS:  calling SOLVE2' */

		    if (solve2 (a,c,ctx,np,n1,12,delta,irank,1)) {
			nerr = 2;
			it = n;
			return (nerr);
			}
/*		    write(*,*) 'NLLS:  back from SOLVE2' */
		    }

/*  for one parameter fit */
		else {
		    if (a[0][0] <= 1.0e-70) {
			nerr = 2;
			it = n;
			return (nerr);
			}
		    c[0][0] = 1.0 / a[0][0];
		    }

/*  Calculate parameters for next iteration */
		for (i = 0; i < np; i++) {
		    tb[i] = 0.0;
		    for (j = 0; j < np; j++) {
			tb[i] = tb[i] + (c[i][j] * b[j]);
			}
		    tb[i] = gp1[i] + (gp2[i] * tb[i]);
		    }
		}

/* Test for convergence of fit */
	    icon = 0;
	    for (i = 0; i < np; i++) {
		if (abs ((gp1[i] / tb[i]) - 1.0) > eps) icon = 1;
		}

/* Set parameters for next iteration */
	    for (i = 0; i < np; i++) {
		gp0[i] = gp1[i];
		gp1[i] = tb[i];
		}

/* Drop out of loop if all parameters converged */
	    if (icon == 0) break;

	    }

/* No convergence */
	if (icon != 0) {
	    nerr = 1;
	    it = nit;
	    }

/* Calculate standard deviations and correlation matrix */
	vf = vf1 / (double) (nfit - np);
	for (i = 0; i < np; i++) {
	    stde[i] = gp2[i] * sqrt (c[i][i] / vf2 * vf);
	    for (j = 0; j < np; j++) {
		co[i][j] = c[i][j] / sqrt (c[i,i] * c[j,j]);
		}
	    }
	return (nerr);

end

double t[12];

solve2 (a,x,b,nv,nrhs,jsize,det,irank,maxit)

/* The equation ax=b is to be solved for x.
 * 	dimension assumptions:
 */
int	jsize		/* left dimension of a, x, and b.
int	nv		/* number of variables.
double	a[jsize][nv]	/* a(jsize,k)   k ge nv;
double	x[jsize][nv]	/* x(jsize,l)   l ge nrhs;
double	b[jsize][nv]	/* b(jsize,l)   l ge nrhs.
int	nrhs;		/* Number of right hand sides.
double	det;		/* Determinant of a.
int	irank;		/* Rank of a.
int	maxit;		/* Maximum number of iterations allowed.
{
    int ir[12];		/* Contains the pivot information.
    double s[12][12];	/* Used for the l-u decomposition of a.

/*     the largest system which may be solved is a 12 x 12.
 *     this routine does no iterative improvement.
 *     it does not destroy either a or b
 */

    int i,j,jdim,it,in, solved;
    double zero, c, xn0;

	jdim = 12;

	zero = 0.0;
	solved = 1;
	if (jsize < nv) {
	    printf ("SOLVE2: %d data points < %d variables\n",jsize,nv);
	    if (in != 0) nrhs = -1;
	    return (1);
	    }
	if (jdim < nv) {
	    printf ("SOLVE2: Too many variables to fit: %d > %d\n",nv,jdim);
	    if (in != 0) nrhs = -1;
	    return (1);
	    }
	xn0 = 0.0;
	in = 0;

/* Calculate zero and store a in s */
	for (i = 0; i< nv; i++) {
	    for (j = 0; j < nv; j++) {
		zero = zero + (a[i][j] * a[i][j]);
		s[i][j] = a[i][j];
		}
	    }
	c = nv * nv;
	zero = sqrt (zero) * 1.0e-16 / c;

	if (nrhs == 0) {
	    solved = 0;
	    zero = 0.0;
	    lineq (nv,jdim,s,x,jsize,nrhs,det,ir,zero,irank,1);
	    return;
	    }

/* Calculate inverse */

/* Set up b */
	if (nrhs < 0) {
	    for (i = 0; i < nv; i++) {
		for (j = 0; j < nv; j++) {
		    b[i][j] = 0.0;
		    }
		b[i][i] = 1.0;
		}
	    nrhs = nv;
	    in = 1;
	    }

/*  Solve the system */
	it = 0;
	for (i = 0; i < nrhs; i++) {
	    for (j = 0; j < nv; j++) {
		x[i][j] = b[i][j];
		}
	    }
	lineq (nv,jdim,s,x,jsize,nrhs,det,ir,zero,irank,1);
	if (nv > irank) 
	    return;
	solved = 1;

/* Unshuffle solution vector x */
	for (j = 0; j < nrhs; j++) {
	    for (i = 0; i < nv; i++) {
		t[i] = x[ir[i]][j];
		}
	    for (i = 0; i < nv; i++) {
		x[i][j] = t[i];
		}
	    }

	if (in != 0)
	    nrhs = -1;

	return (0);
}


void
lineq (nv,jdim,xmat,soln,jsize,nrhs,determ,ir,tol,irank,igo)

int	jdim,nv;
double	xmat[jdim][nv];		/* Matrix of left dimension jdim */
int	jsize,nrhs;
double	soln[jsize][nrhs];	/* Right hand side to contain the answer at the end */
double	determ;
int	ir[nv];
double	tol;
int	irank;
int	igo;		/* controls the program flow as follows: */
			/* =1  form the lu decomposition of xmat and get first sol */
			/* =2  use lu information and get solution */

int i,j,k,l,idiag,irp;
double hold,sum;

define	maxpiv_	10
define	elim_	20
define	entigo2_ 30

{

/*  branch on igo */
    idiag = 1;
    if (igo > 1) {
	entigo2 = 1;
	maxpiv = 0;
	}
    else {
	entigo = 0;
	maxpiv = 1;
	}

/* Initialize for lu decomposition */
	irank = nv
	for (i = 0; i < nv; i++) {
	    ir[i] = i;
	    }
	determ = 1.0;

/* Find max pivot */
maxpiv_
    if (maxpiv) {
	hold = 0.0;
	irp = idiag;
	do i = idiag, nv {
	    if (hold < dabs (xmat[ir[i],idiag])) {
		irp = i;
		hold = dabs (xmat[ir[irp],idiag]);
		}
	    }

/* Update determinant and check for pivot too small */
	determ = determ * xmat[ir[irp],idiag];
	if (hold <= tol) {
	    irank = idiag - 1;
	    return;
	    }

/* Do row and column switches if necessary */
	if (ir[irp] != ir[idiag]) {
	    determ = -determ;
	    i = ir[irp];
	    ir[irp] = ir[idiag];
	    ir[idiag] = i;
	    }
	}

/* Actual elimination and entry for igo = 2 */
elim_

/* Back substitution
	if (idiag == nv)  {
	    if (nrhs == 0)
		return;

	    hold = xmat[ir[nv]][nv];
	    for (i = 0; i < nrhs; i++) {
		soln[ir[nv]][i] = soln[ir[nv]][i] / hold;
		}
	    i = nv-1;
	    k = nv - 2;
	    do i = nv, 2, -1 {
		hold = xmat[ir[k]][k];
		for (j = 0; j < nrhs; j++) {
		    sum = 0.0;
		    for (l = i; l < nv; l++) {
			sum = sum + (xmat[ir[k]][l] * soln[ir[l]][j]);
			}
		    soln[ir[k]][j] = (soln[ir[k]][j] - sum) / hold;
		    }
		k = k - 1;
		}
	    return;
	    }

entigo2_
	k = idiag + 1
	do i = k,nv {
	    hold = xmat(ir[i],idiag) / xmat(ir(idiag),idiag)
	    if (igo == 1) {
		do j = k, nv {
		    xmat(ir[i],j) = xmat(ir[i],j) - (xmat[ir[idiag],j] * hold)
		    }
		}
	    if (igo == 2 || nrhs != 0) {
		do j = 1, nrhs {
		    soln(ir[i],j) = soln(ir[i],j) - (soln[ir[idiag],j] * hold)
		    }
		}
	    }
	idiag = idiag + 1
	if (igo == 1)
	    go to maxpiv_
	go to elim_

end


/*--- Compute triple gaussian and derivatives for nlls fit */

tgauss (sbuff, j, np, gp1, deriv, dyj, diff, j0, compute, debug)

double	*xbuff;		/* Data points */
double	*ybuff;		/* Data points */
int	j;		/* Index into data array */
int	np;		/* Number of parameters being fit */
double	*gp1;		/* Values for this round of parameters */
double	*deriv;		/* Computed derivatives (returned) */
double	dyj;		/* Computed value for this point */
double	diff;		/* Observed - computed data for this point */
int	j0;		/* Initial pixel for background and slope */
int	compute;	/* If nonzero, recompute constants */
int	debug;		/* if nonzero, print intermediate values */

double	sj;
double	db,dt1,dt1s,dg1,dt2,dt2s,dg2,dt3,dt3s,dg3,df;
int	i;

common/elfit/ back,bslope,t1,a1,s1,t2,a2,s2,t3,a3,s3,sj0

double	t1;		/* pixel of first peak */
double	t2;		/* pixel of second peak */
double	t3;		/* pixel of third peak */
double	a1;		/* height of first gaussian peak */
double	a2;		/* height of second gaussian peak */
double	a3;		/* height of third gaussian peak */
double	s1;		/* half-width of first gaussian */
double	s2;		/* half-width of second gaussian */
double	s3;		/* half-width of third gaussian */
double	back;		/* background level at pixel j0 */
double	bslope;		/* background slope from pixel j0 */
double	sj0;		/* pixel for background */

common/prsub/ s1s,s1c,ts1s,s2s,s2c,ts2s,s3s,s3c,ts3s
double	s1s,s2s,s3s	/* sigma squared for three line profiles */
double	s1c,s2c,s3c	/* sigma cubed for three line profiles */
double	ts1s,ts2s,ts3s	/* 2 * sigma squared for three line profiles */

begin

/* Zero derivatives */
	if (np > 0) {
	    if (j > 0) {
		do i = 1, 12 {
		    deriv[i] = 0.0d0
		    }
		}

/* Set values of parameters currently being fit */
	    if (compute) {
		sj0 = double (j0)
		back = gp1[1]
		bslope = gp1[2]
		if (np > 2) {
		    t1 = gp1[3]
		    a1 = gp1[4]
		    s1 = gp1[5]
		    }
		else {
		    t1 = 0.d0
		    a1 = 0.d0
		    s1 = 0.d0
		    }
		if (np > 5) {
		    t2 = gp1[6]
		    a2 = gp1[7]
		    s2 = gp1[8]
		    }
		else {
		    t2 = 0.d0
		    a2 = 0.d0
		    s2 = 0.d0
		    }
		if (np > 8) {
		    t3 = gp1[9]
		    a3 = gp1[10]
		    s3 = gp1[11]
		    }
		else {
		    t3 = 0.d0
		    a3 = 0.d0
		    s3 = 0.d0
		    }
		}
	    }

/* Set constants to speed computations in loop */
	if (compute) {
	    s1s = s1 * s1
	    s1c = s1s * s1
	    ts1s = 2.d0 * s1s

	    s2s = s2 * s2
	    s2c = s2s * s2
	    ts2s = 2.d0 * s2s

	    s3s = s3 * s3
	    s3c = s3s * s3
	    ts3s = 2.d0 * s3s
	    }

	sj = double (j)

/* Compute background */
	db = back + (bslope * (sj - sj0))

/* Compute Gaussian profiles for each part of n-line */
	if ((a1 != 0) && (s1 > 0)) {
	    dt1 = sj - t1
	    dt1s = dt1 * dt1
	    dg1 = a1 * dexp (-dt1s / ts1s)
	    }
	else
	    dg1 = 0.d0

	if ((a2 != 0) && (s2 > 0)) {
	    dt2 = sj - t2
	    dt2s = dt2 * dt2
	    dg2 = a2 * dexp (-dt2s / ts2s)
	    }
	else
	    dg2 = 0.d0

	if ((a3 != 0) && (s3 > 0)) {
	    dt3 = sj - t3
	    dt3s = dt3 * dt3
	    dg3 = a3 * dexp (-dt3s / ts3s)
	    }
	else
	    dg3 = 0.d0

/* Add background and gaussians for all 3 lines */
	dyj = db + dg1 + dg2 + dg3

/* Derivatives for background */
	deriv[1] = 1.d0
	deriv[2] = sj - sj0

/* Derivatives  of center, amplitude, and half-width of first line */
	if (np > 2) {
	    deriv[3] = dg1 * dt1 / s1s
	    deriv[4] = dg1 / a1
	    deriv[5] = dg1 * dt1s / s1c
	    }

/* Derivatives  of center, amplitude, and half-width of second line */
	if (np > 5) {
	    deriv[6] = dg2 * dt2 / s2s
	    deriv[7] = dg2 / a2
	    deriv[8] = dg2 * dt2s / s2c
	    }

/* Derivatives  of center, amplitude, and half-width of third line */
	if (np > 8) {
	    deriv[9] = dg3 * dt3 / s3s;
	    deriv[10] = dg3 / a3;
	    deriv[11] = dg3 * dt3s / s3c;
	    }

/* Compute difference  */
	if (np > 0) {
	    diff = (double)sbuff[j] - dyj;
	    df = abs (diff);
	    }
/*	if (debug) {
 *	    printf ("%d: %.4f %.4f %.4f\n", j, sbuff[j], dyj, diff);
 *	    printf ("  %.4f %.4f %.4f %.4f\n", db, dg1, dg2, dg3);
 *	    for (i = 1; i < np; i++) {
 *		printf (" %g", deriv[i]);
 *		}
 *	    printf ("\n");
 *	    }
 */

	return;
end

/* Aug 29 1997	New subroutines, based on RVSAO.EMSAO SPP fitting code
 */
