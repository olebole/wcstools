/* File xy2sky.c
 * November 1, 1996
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include "libwcs/wcs.h"

static void usage();
extern struct WorldCoor *GetWCSFITS ();	/* Read WCS from FITS or IRAF header */
extern int pix2wcst();

static int verbose = 0;		/* verbose/debugging flag */
static char coorsys[16];

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    char wcstring[40];
    char lstr = 40;
    int ndecset = 0;
    int degout = 0;
    int ndec = 3;
    double x, y;
    FILE *fd;
    char *ln, *listname;
    char line[200];
    char *fn;
    struct WorldCoor *wcs;
    char xstr[32], ystr[32];
    *coorsys = 0;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
    	switch (c) {

    	case 'v':	/* more verbosity */
    	    verbose++;
    	    break;

	case 'b':
	    strcpy (coorsys,"b1950");
    	    break;

       case 'd':
            degout++;
            if (!ndecset) {
                ndec = 5;
		ndecset++;
		}
            break;

	case 'g':
	    strcpy (coorsys,"galactic");
    	    break;

	case 'j':
	    strcpy (coorsys,"j2000");
    	    break;

	case 'n':
	    if (ac < 2)
		usage();
	    ndec = atoi (*++av);
	    ndecset++;
	    ac--;
	    break;

    	default:
    	    usage();
    	    break;
    	}
    }

    /* There are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage ();

    fn = *av++;
    if (verbose)
	printf ("%s:\n", fn);
    wcs = GetWCSFITS (fn);
    if (nowcs (wcs))
	exit (1);
    if (*coorsys)
	wcsoutinit (wcs, coorsys);
    if (ndecset)
	wcs->ndec = ndec;
    if (degout)
	wcs->degout = degout;
    while (ac-- > 1) {
	listname = *av;
	if (listname[0] == '@') {
	    ln = listname;
	    while (*ln++)
		*(ln-1) = *ln;
	    if (fd = fopen (listname, "r")) {
		while (fgets (line, 200, fd)) {
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
usage ()
{
    fprintf (stderr,"Compute RA Dec from X Y using WCS in FITS and IRAF image files\n");
    fprintf (stderr,"Usage: [-vbjg] [-n ndec] file.fts x1 y1 ... xn yn\n");
    fprintf (stderr,"Usage: [-bdjgv] [-n ndec] file.fts @listfile\n");
    fprintf (stderr,"  -b: B1950 (FK4) output\n");
    fprintf (stderr,"  -d: RA and Dec output in degrees\n");
    fprintf (stderr,"  -j: J2000 (FK5) output\n");
    fprintf (stderr,"  -g: galactic longitude and latitude output\n");
    fprintf (stderr,"  -n: number of decimal places in output RA seconds\n");
    fprintf (stderr,"  -v: verbose\n");
    exit (1);
}
/*
 * Feb 23 1996	New program
 * Apr 24 1996	Version 1.1: Add B1950, J2000, or galactic coordinate output options
 * Jun 10 1996	Change name of subroutine which reads WCS
 * Aug 28 1996	Remove unused variables after lint
 * Nov  1 1996	Add options to set number of decimal places and output degrees
 */
