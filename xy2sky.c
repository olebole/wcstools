/* File xy2sky.c
 * May 13, 1998
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
static int append = 0;		/* append input line flag */
static int tabtable = 0;	/* tab table output flag */
static int oldwcs = 0;		/* AIPS classic WCS flag */
static char coorsys[16];
static int linmode = -1;

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    char wcstring[64];
    char lstr = 64;
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
    char temp[32];
    *coorsys = 0;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
    	switch (c) {

    	case 'v':	/* more verbosity */
    	    verbose++;
    	    break;

	case 'a':	/* Append input line to sky position */
	    append++;
    	    break;

	case 'b':	/* Output B1950 coordinates */
	    strcpy (coorsys,"b1950");
    	    break;

       case 'd':	/* Output degrees instead of hh:mm:ss dd:mm:ss */
            degout++;
            break;

	case 'e':	/* Output galactic coordinates */
	    strcpy (coorsys,"ecliptic");
            degout++;
    	    break;

	case 'g':	/* Output galactic coordinates */
	    strcpy (coorsys,"galactic");
            degout++;
    	    break;

	case 'j':	/* Output J2000 coordinates */
	    strcpy (coorsys,"j2000");
    	    break;

	case 'm':	/* Mode for output of linear coordinates */
	    if (ac < 2)
		usage();
	    linmode = atoi (*++av);
	    ac--;
	    break;

	case 'n':	/* Number of decimal places in output sec or deg */
	    if (ac < 2)
		usage();
	    ndec = atoi (*++av);
	    ndecset++;
	    ac--;
	    break;

	case 'q':	/* Equinox for output */
	    if (ac < 2)
		usage();
	    strcpy (coorsys, *++av);
	    ac--;
	    break;

	case 't':	/* Output tab table */
	    tabtable++;
    	    break;

    	case 'z':	/* Use AIPS classic WCS */
    	    oldwcs++;
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
    if (nowcs (wcs)) {
	printf ("%s: No WCS for file, cannot compute image size\n", fn);
	exit(1);
	}
    wcs->oldwcs = oldwcs;
    if (linmode > -1)
	setlinmode (wcs, linmode);
    if (*coorsys)
	wcsoutinit (wcs, coorsys);
    if (tabtable) {
	wcs->tabsys = 1;
	if (!append) {
	    if (wcs->sysout == WCS_B1950 || wcs->sysout == WCS_J2000)
		printf ("X    	Y    	RA      	Dec     	Equinox\n");
	    else if (wcs->sysout == WCS_GALACTIC)
		printf ("X    	Y    	Gal Long 	Gal Lat 	Equinox\n");
	    else if (wcs->sysout == WCS_ECLIPTIC)
		printf ("X    	Y    	Ecl Long 	Ecl Lat 	Equinox\n");
	    printf ("-----	-----	--------	--------	-------\n");
	    }
	}
    if (degout) {
	wcs->degout = degout;
        if (!ndecset) {
            ndec = 5;
	    ndecset++;
	    }
	}
    if (ndecset)
	wcs->ndec = ndec;
    while (ac-- > 1) {
	listname = *av;
	if (listname[0] == '@') {
	    ln = listname;
	    while (*ln++)
		*(ln-1) = *ln;
	    if (strcmp (listname,"STDIN")==0 || strcmp (listname,"stdin")==0)
		fd = stdin;
	    else
		fd = fopen (listname, "r");
	    if (fd != NULL) {
		while (fgets (line, 200, fd)) {
		    sscanf (line,"%s %s", xstr, ystr);
		    x = atof (xstr);
		    y = atof (ystr);
		    if (pix2wcst (wcs, x, y, wcstring, lstr)) {
			if (wcs->sysout == WCS_ECLIPTIC) {
			    sprintf(temp,"%.5f",wcs->epoch);
			    strcat (wcstring, " ");
			    strcat (wcstring, temp);
			    }
			if (append) {
			    if (tabtable)
				printf ("%s	%s", wcstring, line);
			    else
				printf ("%s %s", wcstring, line);
			    }
			else if (degout) {
			    if (wcs->nxpix > 9999 || wcs->nypix > 9999) {
				if (tabtable)
				    printf ("%9.3f	%9.3f	%s\n",x, y, wcstring);
				else
				    printf ("%9.3f %9.3f %s\n",x, y, wcstring);
				}
			    else {
				if (tabtable)
				    printf ("%8.3f	%8.3f	%s\n",x, y, wcstring);
				else
				    printf ("%8.3f %8.3f %s\n",x, y, wcstring);
				}
			    }
			else if (tabtable)
			    printf ("%.3f	%.3f	%s\n",x, y, wcstring);
			else
			    printf ("%.3f %.3f -> %s\n",x, y, wcstring);
			}
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
	    if (pix2wcst (wcs, x, y, wcstring, lstr)) {
		if (wcs->sysout == WCS_ECLIPTIC) {
		    sprintf(temp,"%.5f",wcs->epoch);
		    strcat (wcstring, " ");
		    strcat (wcstring, temp);
		    }
		if (tabtable)
		    printf ("%.3f	%.3f	%s\n",x, y, wcstring);
		else
		    printf ("%.3f %.3f -> %s\n",x, y, wcstring);
		}
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
    fprintf (stderr,"Usage: [-abdjgv] [-n ndec] file.fits x1 y1 ... xn yn\n");
    fprintf (stderr,"Usage: [-abdjgv] [-n ndec] file.fits @listfile\n");
    fprintf (stderr,"  -a: append input line after output position\n");
    fprintf (stderr,"  -b: B1950 (FK4) output\n");
    fprintf (stderr,"  -d: RA and Dec output in degrees\n");
    fprintf (stderr,"  -e: ecliptic longitude and latitude output\n");
    fprintf (stderr,"  -g: galactic longitude and latitude output\n");
    fprintf (stderr,"  -j: J2000 (FK5) output\n");
    fprintf (stderr,"  -m: mode for output of LINEAR WCS coordinates\n");
    fprintf (stderr,"  -n: number of decimal places in output RA seconds\n");
    fprintf (stderr,"  -q: output equinox if not 2000 or 1950\n");
    fprintf (stderr,"  -t: tab table output\n");
    fprintf (stderr,"  -v: verbose\n");
    fprintf (stderr,"  -z: use AIPS classic projections instead of WCSLIB\n");
    exit (1);
}
/*
 * Feb 23 1996	New program
 * Apr 24 1996	Version 1.1: Add B1950, J2000, or galactic coordinate output options
 * Jun 10 1996	Change name of subroutine which reads WCS
 * Aug 28 1996	Remove unused variables after lint
 * Nov  1 1996	Add options to set number of decimal places and output degrees
 *
 * Dec 15 1997	Print message if no WCS; read IRAF 2.11 header format
 * Dec 15 1997	Drop -> if output sky coordinates are in degrees
 * Dec 31 1997	Allow entire input line to be appended to sky position
 *
 * Jan  7 1998	Apply WFPC and WFPC2 pixel corrections if requested
 * Jan  7 1998	Add tab table output using -t
 * Jan 26 1998	Implement Mark Calabretta's WCSLIB
 * Jan 29 1998	Add -z for AIPS classic WCS projections
 * Feb 18 1998	Version 2.0: Full Calabretta implementation
 * Mar 12 1998	Version 2.1: IRAF TNX projection added
 * Mar 27 1998	Version 2.2: Polynomial plate fit added
 * Apr 24 1998	Increase size of WCS string from 40 to 64
 * Apr 28 1998	Change coordinate system flag to WCS_*
 * Apr 28 1998	Add output mode for linear coordinates
 * Apr 28 1998	Add ecliptic coordinate system output
 * May 13 1998	Allow arbitrary equinox for output coordinates
 */
