/* File sky2xy.c
 * November 4, 1997
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
#include "libwcs/fitshead.h"

static void usage();
extern struct WorldCoor *GetWCSFITS ();	/* Read WCS from FITS or IRAF header */

static int verbose = 0;		/* verbose/debugging flag */
static char coorsys[16];

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    double x, y, ra, dec, ra0, dec0;
    FILE *fd;
    char *ln, *listname;
    char line[80];
    char *fn;
    char csys[16];
    struct WorldCoor *wcs;
    char rastr[32], decstr[32];
    int offscale;

    *coorsys = 0;

    /* Decode arguments */
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
    if (nowcs (wcs)) {
	fprintf (stderr, "No WCS in image file %s\n", fn);
	exit (1);
	}
    if (*coorsys)
	wcsoutinit (wcs, coorsys);

    while (ac-- > 1) {
	listname = *av;
	if (listname[0] == '@') {
	    ln = listname;
	    while (*ln++)
		*(ln-1) = *ln;
	    if (fd = fopen (listname, "r")) {
		while (fgets (line, 80, fd)) {
		    sscanf (line,"%s %s", rastr, decstr);
		    ra = str2ra (rastr);
		    dec = str2dec (decstr);
		    wcs2pix (wcs, ra, dec, &x, &y, &offscale);
		    if (verbose)
			printf ("%s %s -> %.5f %.5f -> %.3f %.3f",
				 rastr, decstr, ra, dec, x, y);
		    else
			printf ("%s %s -> %.3f %.3f",rastr, decstr, x, y);
		    if (offscale)
			printf (" (offscale)\n");
		    else
			printf ("\n");
		    }
		}
	    else {
		fprintf (stderr, "Cannot read file %s\n", listname);
		exit (1);
		}
	    av++;
	    }
	else if (ac > 1) {
	    strcpy (rastr, *av);
	    ac--;
	    av++;
	    strcpy (decstr, *av);
	    ra0 = str2ra (rastr);
	    dec0 = str2dec (decstr);
	    ra = ra0;
	    dec = dec0;
	    av++;

	/* Convert coordinates system to that of image */
	    if (ac > 1) {
		strcpy (csys, *av);
		if (csys[0] == 'B' || csys[0] == 'b') {
		    if (wcs->equinox == 2000.0)
			fk425e (&ra, &dec, wcs->epoch);
		    ac--;
		    av++;
		    }
		else if (csys[0] == 'J' || csys[0] == 'j') {
		    if (wcs->equinox == 1950.0)
			fk524e (&ra, &dec, wcs->epoch);
		    ac--;
		    av++;
		    }
		else if ((!strcmp (csys,"FK4") || !strcmp (csys,"fk4")) &&
			wcs->equinox == 2000.0) {
			fk425e (&ra, &dec, wcs->epoch);
		    ac--;
		    av++;
		    }
		else if ((!strcmp (csys,"FK5") || !strcmp (csys,"fk5")) &&
			wcs->equinox == 1950.0) {
			fk524e (&ra, &dec, wcs->epoch);
		    ac--;
		    av++;
		    }
		else {
		    if (wcs->equinox == 1950.0)
			strcpy (csys, "B1950");
		    else
			strcpy (csys, "J2000");
		    }
		}
	    else {
		if (wcs->pcode < 0)
		    strcpy (csys, "PIXEL");
		else if (wcs->equinox == 1950.0)
		    strcpy (csys, "B1950");
		else
		    strcpy (csys, "J2000");
		}

	    if (ra != ra0 || verbose) {
		printf ("%s %s %s -> ", rastr, decstr, csys);
		ra2str (rastr, ra, 3);
		dec2str (decstr, dec, 2);
		if (wcs->equinox == 1950.0)
		    strcpy (csys, "B1950");
		else if (wcs->equinox == 2000.0)
		    strcpy (csys, "J2000");
		printf ("%s %s %s\n", rastr, decstr, csys);
		}
	    wcs2pix (wcs, ra, dec, &x, &y, &offscale);
	    printf ("%s %s %s -> %.3f %.3f",rastr, decstr, csys, x, y);
	    if (offscale)
		printf (" (offscale)\n");
	    else
		printf ("\n");
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
    fprintf (stderr,"Compute X Y from RA Dec using WCS in FITS and IRAF image files\n");
    fprintf(stderr,"%s: usage: [-vbjg] file.fts ra1 dec1 ... ran decn\n", progname);
    fprintf(stderr,"%s: usage: [-vbjg] file.fts @listfile\n", progname);
    fprintf(stderr,"  -b: B1950 (FK4) input\n");
    fprintf(stderr,"  -j: J2000 (FK5) input\n");
    fprintf(stderr,"  -g: galactic longitude and latitude input\n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}
/* Feb 23 1996	New program
 * Apr 24 1996	Version 1.1: Add B1950, J2000, or galactic coordinate input options
 * Jun 10 1996	Change name of WCS subroutine
 * Aug 27 1996	Clean up code after lint
 * Oct 29 1996	Allow alternate coordinate systems for input coordinates
 * Oct 30 1996	Exit if image file is not found
 * Nov  1 1996	Fix bug so systemless coordinates do not cause crash
 * Nov  5 1996	Fix multiple sets of coordinates on command line
 *
 * Jun  4 1997	Add PIXEL wcs for linear non-sky projections
 * Nov  4 1997	If verbose mode, always print converted input string
 */
