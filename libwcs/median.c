/*** file iolib/median.c
 *** July 27, 1999
 *** After Press, et al. quicksort
 */

double
median (x, n, xm1, xm2)

int	n;		/* Size of vector */
double	*x;		/* Vector of numbers of which median is to be found */
double	*xm1, *xm2;	/* Quartile points (returned) */
{
    double	xmed;	/* Median (returned) */
    double	xx;
    int	l,ir,i,j,n2,n14,n34;

    if (n < 1)
	return (0.0);
    else if (n == 1)
	return (x[0]);

    l = n / 2;
    ir = n - 1;

    while (l >= 0) {
	if (l > 0) {
	    l = l - 1;
	    xx = x[0];
	    }
	else {
	    xx = x[ir];
	    x[ir] = x[0];
	    ir = ir - 1;
	    if (ir == 0) {
		x[0] = xx;
		n2 = n / 2;
		if (2*n == n)
		    xmed = 0.5 * (x[n2] + x[n2+1]);
		else
		    xmed = x[n2];
		n14 = (int) (n * 0.25);
		n34 = (int) (n * 0.75);
		*xm1 = x[n14];
		*xm2 = x[n34];
		return (xmed);
		}
	    }
	i = l;
	j = l + l;

	while (j <= ir) {
	    if (j < ir) {
		if (x[j] < x[j+1])
		    j = j + 1;
		if (xx < x[j]) {
		    x[i] = x[j];
		    i = j;
		    j = j + j;
		    }
		else
		    j = ir + 1;
		}
	    }
	x[i] = xx;
	}
}
