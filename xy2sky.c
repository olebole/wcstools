/* File xy2sky.c
 * April 24, 1996
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
static char *RevMsg = "XY2SKY 1.1, 8 August 1996, Doug Mink, SAO";
static char coorsys[16];

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
	case 'j':
	    strcpy (coorsys,"j2000");
    	    break;
	case 'g':
	    strcpy (coorsys,"galactic");
    	    break;
    	default:
    	    usage(progname);
    	    break;
    	}
    }

    /* There are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage (progname);

    fn = *av++;
    if (verbose)
	printf ("%s:\n", fn);
    wcs = GetWCSFITS (fn);
    if (*coorsys)
	wcsoutinit (wcs, coorsys);
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
    fprintf(stderr,"%s: usage: [-vbjg] file.fts x1 y1 ... xn yn\n", progname);
    fprintf(stderr,"%s: usage: [-vbjg] file.fts @listfile\n", progname);
    fprintf(stderr,"  -b: B1950 (FK4) output\n");
    fprintf(stderr,"  -j: J2000 (FK5) output\n");
    fprintf(stderr,"  -g: galactic longitude and latitude output\n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}
/*
 * Feb 23 1996	New program
 * Apr 24 1996	Version 1.1: Add B1950, J2000, or galactic coordinate output options
 * Jun 10 1996	Change name of subroutine which reads WCS
 */
