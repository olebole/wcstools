/*** File fitlib/polfit.c
 *** February 23, 1998
 *** By Doug Mink, after Bevington, page 141

 *--- Polynomial least squares fitting program, almost identical to the
 *    one in Bevington, "Data Reduction and Error Analysis for the
 *    Physical Sciences," page 141.  I changed the argument list and
 *    removed the weighting.
 *      y = a(1) + a(2)*(x-x0) + a(3)*(x-x0)**2 + a(3)*(x-x0)**3 + . . .
 */

static double determ();

void
polfit (x, y, x0, npts, nterms, a, chisqr)

double *x;	/* Array of independent variable points */
double *y;	/* Array of dependent variable points */
double x0;	/* Offset to independent variable */
int npts;	/* Number of data points to fit */
int nterms;	/* Number of parameters to fit */
double *a;	/* Vector containing current fit values */
double *chisqr;
{
    double chisq;
    double xterm,yterm,xi,yi,free;
    double sumx[19],sumy[10],array[100];
    int i,j,k,l,n,nmax;
    double delta;
    int adim = 10;

    /* accumulate weighted sums */
    nmax = 2 * nterms - 1;
    for (n = 0; n < nmax; n++)
	sumx[n] = 0.0;
    for (j = 0; j < nterms; j++)
	sumy[j] = 0.0;
    chisq = 0.0;
    for (i = 0; i < npts; i++) {
	xi = x[i] - x0;
	yi = y[i];
	xterm = 1.0;
	for (n = 0; n < nmax; n++) {
	    sumx[n] = sumx[n] + xterm;
	    xterm = xterm * xi;
	    }
	yterm = yi;
	for (n = 0; n < nterms; n++) {
	    sumy[n] = sumy[n] + yterm;
	    yterm = yterm * xi;
	    }
	chisq = chisq + yi*yi;
	}

    /* Construct matrices and calculate coeffients */
    for (j = 0; j < nterms; j++) {
	for (k = 1; k < nterms; k++) {
	    n = j + k - 1;
	    array[j+k*adim] = sumx[n];
	    }
	}
    delta = determ (array, nterms, adim);
    if (delta == 0.0) {
	*chisqr = 0.;
	for (j = 0; j < nterms; j++)
	    a[j] = 0. ;
	return;
	}

    for (l = 0; l < nterms; l++) {
	for (j = 0; j < nterms; j++) {
	    for (k = 0; k < nterms; k++) {
		n = j + k;
		array[j+k*adim] = sumx[n];
		}
	    array[j+l*adim] = sumy[j];
	    }
	a[l] = determ (array, nterms, adim) / delta;
	}

    /* Calculate chi square */
    for (j = 0; j < nterms; j++) {
	chisq = chisq - (2.0 * a[j] * sumy[j]);
	for (k = 0; k < nterms; k++) {
	    n = j + k - 1;
	    chisq = chisq + (a[j] * a[k] * sumx[n]);
	    }
	}
    free = npts - nterms;
    *chisqr = chisq / free;

    return;
}


/*--- Calculate the determinant of a square matrix
 *    This subprogram destroys the input matrix array
 *    From Bevington, page 294.
 */

static double
determ (array, norder, adim)

double	*array;		/* Input matrix array */
int	norder;		/* Order of determinant (degree of matrix) */
int	adim;		/* Dimension of 2-D array */

{
    double save, det;
    int i,j,k,k1, zero;

    det = 1.0;
    for (k = 0; k < norder; k++) {

	/* Interchange columns if diagonal element is zero */
	if (array[k+k*adim] == 0) {
	    zero = 1;
	    for (j = k; j < norder; j++) {
		if (array[k+j*adim] != 0.0)
		    zero = 0;
		}
	    if (zero)
		return (0.0);

	    for (i = k; i < norder; i++) {
		save = array[i+j*adim]; 
		array[i+j*adim] = array[i+k*adim];
		array[i+k*adim] = save ;
		}
	    det = -det;
	    }

	/* Subtract row k from lower rows to get diagonal matrix */
	det = det * array[k+k*adim];
	if (k < norder - 1) {
	    k1 = k + 1;
	    for (i = k1; i < norder; i++) {
		for (j = k1; j < norder; j++) {
		    array[i+j*adim] = array[i+j*adim] -
				      (array[i+k*adim] * array[k+j*adim] /
				      array[k+k*adim]);
		    }
		}
	    }
	}
	return (det);
}
/* Sep 10 1987	Program written
 *
 * Mar 17 1993	Add x offset
 *
 * Feb 23 1998	Translate to C
 */
