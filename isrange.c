/* File isrange.c
 * April 11, 2005
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 *
 * Return 1 if argument is a range
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static int isrange();

main (ac, av)
int ac;
char **av;
{
    char *fn;
    char *str;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help")) {
	fprintf (stderr,"Usage:  Return 1 if argument is a range of numbers: n1[-n2[xs]],...\n");
	fprintf (stderr,"        where n1=first number, n2=last number, and s=step size.\n");
	exit (1);
	}

    else if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	exit (1);
	}

    /* check to see if this is a range */
    else
	printf ("%d\n", isrange (str));

    exit (0);
}


/* ISRANGE -- Return 1 if string is a range, else 0 */

static int
isrange (string)

char *string;		/* String which might be a range of numbers */

{
    int i, lstr;

    /* If range separators present, check to make sure string is range */
    if (strchr (string+1, '-') || strchr (string+1, ',')) {
	lstr = strlen (string);
	for (i = 0; i < lstr; i++) {
	    if (strchr ("0123456789-,.x", (int)string[i]) == NULL)
		return (0);
	    }
	return (1);
	}
    else
	return (0);
}
/* Dec 14 2001	New program
 *
 * Apr 11 2005	Print version; improve online documentation
 */
