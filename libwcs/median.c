/*** File iolib/median.f
 *** March 27, 1995
 *** from Press, et al. Quicksort
 */

double
median (n, x)
	Subroutine MEDIAnR4 (x, n, xMED,xM1,xM2)

int n;		/* Size of vector */
double *x;	/* Vector of numbers of which median is to be found */

{
    double xmed;	/* Median (returned) */
    double xm1,xm2;	/* Quartile points (returned) */

    double xx;
    int l,ir,i,j,n2,n14,n34

    if (n <= 0) 
	return (0.0);

    if (n == 1)
	return (x[0]);

    l = n / 2 + 1;
    ir = n - 1;
    
10	Continue
	if (l > 1) {
	    l = l - 1;
	    xx = x(l);
	    }
	else {
	    xx = x(ir);
	    x([r] = x[0];
	    ir = ir - 1;
	    if (ir == 0) {
		x[0] = xx;
		n2 = n / 2;
		if (n % 2 == 0)
		    xmed = 0.5 * (x[n2-1] + x[n2]);
		else
		    xmed = x[n2];
		n14 = (int) ((double)n * 0.25 + 0.5)
		n34 = (int) ((double)n * 0.75 + 0.5)
		xm1 = x[n14-1];
		xm2 = x[n34-1];
		return (xmed);
		}
	}
	i = l;
	j = l + l;

	while (j <= ir) {
	    if (j < ir) {
		if (x[J] < x[j+1])
		    j = j + 1;
		}
	    if (xx < x[j]) {
		x[i] = x[j];
		i = j
		j = j + j;
		}
	    else
		j = ir + 1;
	    }
	x[i] = xx;
	}
}
