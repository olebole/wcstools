/* File sky2xy.c
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
static char *RevMsg = "SKY2XY version 1.0, 23 February 1996";


main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    char wcstring[32];
    char lstr = 32;
    double x, y, ra, dec;
    FILE *fd;
    char *ln, *listname;
    char line[80];
    char *fn;
    struct WorldCoor *wcs;
    char rastr[32], decstr[32];
    int offscale;

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
		    sscanf (line,"%s %s", rastr, decstr);
		    str2ra (rastr, &ra);
		    str2dec (decstr, &dec);
		    wcs2pix (wcs, ra, dec, &x, &y, &offscale);
		    if (offscale)
			printf ("%s %s -> offscale\n",rastr, decstr);
		    else
			if (verbose)
			    printf ("%s %s -> %.5f %.5f -> %.3f %.3f\n",
				    rastr, decstr, ra, dec, x, y);
			else
			    printf ("%s %s -> %.3f %.3f\n",rastr, decstr, x, y);
		    }
		}
	    else
		fprintf (stderr, "Cannot read file %s\n", listname);
	    av++;
	    }
	else if (ac > 1) {
	    strcpy (rastr, *av);
	    str2ra (rastr, &ra);
	    ac--;
	    strcpy (decstr, *++av);
	    str2dec (decstr, &dec);
	    wcs2pix (wcs, ra, dec, &x, &y, &offscale);
	    if (offscale)
		printf ("%s %s -> offscale\n",rastr, decstr);
	    else
		printf ("%s %s -> %.3f %.3f\n",rastr, decstr, x, y);
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
    fprintf (stderr,"Compute X Y from RA Dec using WCS in FITS and IRAF image files\n");
    fprintf (stderr, "By D. Mink, SAO\n");
    fprintf(stderr,"%s: usage: [-v] file.fts ra1 dec1 ... ran decn\n", progname);
    fprintf(stderr,"%s: usage: [-v] file.fts @listfile\n", progname);
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}
