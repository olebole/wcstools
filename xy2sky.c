/* File xy2sky.c
 * February 23, 1996
 * By Elwood Downey, revised by Doug Mink
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include "wcs.h"

static void usage();
extern struct WorldCoor *getWCS ();	/* Read WCS from FITS or IRAF header */
extern int pix2wcst();

static int verbose = 0;		/* verbose/debugging flag */
static char *RevMsg = "XY2SKY version 1.0, 23 February 1996";


main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    char wcstring[32];
    char lstr = 32;
    double x, y;
    FILE *fd;
    char *ln, *listname;
    char line[80];
    char *fn;
    struct WorldCoor *wcs;
    char xstr[32], ystr[32];

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
    	switch (c) {
    	case 'v':	/* more verbosity */
    	    verbose++;
    	    break;
    	default:
    	    usage(progname);
    	    break;
    	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage (progname);

    fn = *av++;
    if (verbose)
	printf ("%s:\n", fn);
    wcs = getWCS (fn);
    while (ac-- > 1) {
	char c;
	listname = *av;
	if (listname[0] == '@') {
	    ln = listname;
	    while (*ln++)
		*(ln-1) = *ln;
	    if (fd = fopen (listname, "r")) {
		while (fgets (line, 80, fd)) {
		    sscanf (line,"%s %s", xstr, ystr);
		    x = atof (xstr);
		    y = atof (ystr);
		    if (pix2wcst (wcs, x, y, wcstring, lstr))
			printf ("%.3f %.3f -> %s\n",x, y, wcstring);
		    }
		}
	    else
		fprintf (stderr, "Cannot read file %s\n", listname);
	    av++;
	    }
	else if (ac > 1) {
	    x = atof (*av);
	    ac--;
	    y = atof (*++av);
	    if (pix2wcst (wcs, x, y, wcstring, lstr))
		printf ("%.3f %.3f -> %s\n",x, y, wcstring);
	    av++;
	    }
	}

    if (wcs)
	free ((char *)wcs);
    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"%s\n",RevMsg);
    fprintf (stderr,"Compute RA Dec from X Y using WCS in FITS and IRAF image files\n");
    fprintf (stderr, "By D. Mink, SAO\n");
    fprintf(stderr,"%s: usage: [-v] file.fts x1 y1 ... xn yn\n", progname);
    fprintf(stderr,"%s: usage: [-v] file.fts @listfile\n", progname);
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}
