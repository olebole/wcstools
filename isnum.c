/* File isnum.c
 * November 7, 2001
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 *
 * Return 1 if argument is an integer, 2 if it is floating point, else 0
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libwcs/fitshead.h"

static void usage();

main (ac, av)
int ac;
char **av;
{
    char *fn;
    char *str;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help")) {
	fprintf (stderr,"ISNUM: Return 1 if argument is an integer, ");
	fprintf (stderr,"2 if it is floating point, else 0\n");
	exit (1);
	}

    /* check to see if this is a number */
    else
	printf ("%d\n", isnum (str));

    exit (0);
}
