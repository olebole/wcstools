/* File sky2xy.c
 * August 6, 1998
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
#include "wcs.h"
#include "fitsfile.h"

static void usage();
extern struct WorldCoor *GetWCSFITS ();	/* Read WCS from FITS or IRAF header */

static int verbose = 0;		/* verbose/debugging flag */
static char coorsys[16];
static double eqin = 0.0;
static double eqout = 0.0;

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
    int sysin;
    struct WorldCoor *wcs;
    char rastr[32], decstr[32];
    int offscale, n;

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
		strcpy (coorsys,"B1950");
    		break;

	    case 'e':
		strcpy (coorsys,"ecliptic");
    		break;

	    case 'g':
		strcpy (coorsys,"galactic");
    		break;

	    case 'j':
		strcpy (coorsys,"J2000");
    		break;

	    case 'z':       /* Use AIPS classic WCS */
		setdefwcs (1);
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

    while (ac-- > 1) {
	listname = *av;
	if (listname[0] == '@') {
	    if (*coorsys)
		wcsininit (wcs, coorsys);
	    ln = listname;
	    while (*ln++)
		*(ln-1) = *ln;
	    if (fd = fopen (listname, "r")) {
		while (fgets (line, 80, fd)) {
		    n = sscanf (line,"%s %s %s", rastr, decstr, csys);
		    ra = str2ra (rastr);
		    dec = str2dec (decstr);
		    if (n > 2)
			sysin = wcscsys (csys);
		    else
			sysin = -1;
		    if (sysin > -1)
			wcsc2pix (wcs, ra, dec, csys, &x, &y, &offscale);
		    else {
			wcs2pix (wcs, ra, dec, &x, &y, &offscale);
			strcpy (csys, coorsys);
			}
		    if (verbose)
			printf ("%s %s %s -> %.5f %.5f -> %.3f %.3f",
				 rastr, decstr, csys, ra, dec, x, y);
		    else
			printf ("%s %s %s -> %.3f %.3f",
				rastr, decstr, csys, x, y);
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
		if (csys[0] == 'B' || csys[0] == 'b' || csys[2] == '4' ||
		    csys[0] == 'J' || csys[0] == 'j' || csys[2] == '5' ||
		    csys[0] == 'G' || csys[0] == 'g' ||
		    csys[0] == 'E' || csys[0] == 'e' ||
		    csys[0] == 'L' || csys[0] == 'l') {
		    ac--;
		    av++;
		    }
		else {
		    strcpy (csys, wcs->radecsys);
		    }
		}
	    else {
		if (wcs->prjcode < 0)
		    strcpy (csys, "PIXEL");
		else {
		    strcpy (csys, wcs->radecsys);
		    }
		}

	    sysin = wcscsys (csys);
	    eqin = wcsceq (csys);
	    wcscon (sysin, wcs->syswcs, eqin, eqout, &ra, &dec, wcs->epoch);
	    if (sysin != wcs->syswcs && verbose) {
		printf ("%s %s %s -> ", rastr, decstr, csys);
		ra2str (rastr, 32, ra, 3);
		dec2str (decstr, 32, dec, 2);
		printf ("%s %s %s\n", rastr, decstr, wcs->radecsys);
		}
	    wcsc2pix (wcs, ra0, dec0, csys, &x, &y, &offscale);
	    printf ("%s %s %s -> %.3f %.3f",rastr, decstr, csys, x, y);
	    if (wcs->wcsl.cubeface)
		printf (" %d", wcszout (wcs));
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
    fprintf(stderr,"%s: usage: [-vbjg] file.fts ra1 dec1 sys1 ... ran decn sysn\n", progname);
    fprintf (stderr,"%s: usage: [-vbjg] file.fts @listfile\n", progname);
    fprintf (stderr,"  -v: verbose\n");
    fprintf (stderr,"  -z: use AIPS classic projections instead of WCSLIB\n");
    fprintf (stderr,"These flags are best used for files of coordinates in the same system:\n");
    fprintf (stderr,"  -b: B1950 (FK4) input\n");
    fprintf (stderr,"  -e: ecliptic longitude and latitude input\n");
    fprintf (stderr,"  -j: J2000 (FK5) input\n");
    fprintf (stderr,"  -g: galactic longitude and latitude input\n");
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
 * Dec 15 1997	Handle new IRAF 2.11 image header format
 *
 * Jan 28 1998  Implement Mark Calabretta's WCSLIB
 * Jan 29 1998  Add -z for AIPS classic WCS projections
 * Feb 17 1998	Add support for galactic coordinates as input
 * Feb 18 1998	Version 2.0: Full Calabretta implementation
 * Mar 27 1998	Version 2.2: Add TNX and polynomial plate fit
 * Apr 13 1998	Compute pixel from galactic coordinates correctly
 * Apr 14 1998	Add ecliptic coordinates
 * Apr 24 1998	Handle linear coodinates
 * Apr 28 1998	Implement separate coordinate system for input
 * May 13 1998	Implement arbitrary equinox for input
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Jun 25 1998	Set WCS subroutine choice with SETDEFWCS()
 * Jul 16 1998	Print face if cube face is returned
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 */
