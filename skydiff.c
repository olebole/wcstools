/* File skydiff.c
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

static void usage();
static void skycons();

static int verbose = 0;		/* verbose/debugging flag */
static double epoch = 0.0;
static double eqout = 0.0;
static double eqin = 0.0;

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    FILE *fd;
    char *ln, *listname;
    char line[80];
    char rastr0[32], decstr0[32];
    char rastr1[32], decstr1[32];
    char csys[32];
    int sys0;
    int sys1 = -1;
    int sys2 = -1;
    double ra, dec;
    int degout = 0;
    int ndec = 3;		/* Number of decimal places in RA seconds */
    int ndecset = 0;
    char coorsys[4][12];
    char coorout[16];
    strcpy (coorsys[0], "J2000");
    strcpy (coorsys[1], "B1950");
    strcpy (coorsys[2], "galactic");
    strcpy (coorsys[3], "ecliptic");

    listname = NULL;

    /* There are ac arguments starting at av[0] */
    if (ac == 1)
	usage (progname);

    /* Decode arguments */
    for (av++; --ac > 0 && (*(str = *av) == '-'); av++) {
	char c;
	while (c = *++str)
    	switch (c) {

    	case 'v':	/* more verbosity */
    	    verbose++;
    	    break;

	case 'b':
	    sys2 = WCS_B1950;
    	    break;

	case 'd':
	    degout++;
	    if (!ndecset)
		ndec = 5;
    	    break;

	case 'e':
	    sys2 = WCS_ECLIPTIC;
	    if (!ndecset)
		ndec = 5;
    	    break;

	case 'g':
	    sys2 = WCS_GALACTIC;
	    if (!ndecset)
		ndec = 5;
    	    break;

	case 'j':
	    sys2 = WCS_J2000;
    	    break;

	case 'q':
	    if (ac < 2)
		usage();
	    strcpy (coorout, *++av);
	    ac--;
    	    break;

	case 'y':
	    if (ac < 2)
		usage();
	    epoch = atof (*++av);
	    ac--;
    	    break;

	case 'n':
	    if (ac < 2)
		usage();
	    ndec = atoi (*++av);
	    ndecset++;
	    ac--;
    	    break;

    	default:
	    if (notnum (str))
    		usage(progname);
    	    break;
    	}
    }

    if (verbose)

    while (ac-- > 0) {
	listname = *av;
	if (listname[0] == '@') {
	    ln = listname;
	    while (*ln++)
		*(ln-1) = *ln;
	    if (fd = fopen (listname, "r")) {
		while (fgets (line, 80, fd)) {
		   csys[0] = 0;
		    sscanf (line,"%s %s %s", rastr0, decstr0, csys);
	    	    if (csys[0] == 0) {
			if (sys2 = WCS_J2000)
			    sys0 = WCS_B1950;
			else
			    sys0 = WCS_J2000;
			}
		    else
			sys0 = wcscsys (csys);
		    if (sys2 < 0) {
			if (sys0 = WCS_J2000)
			    sys2 = WCS_B1950;
			else
			    sys2 = WCS_J2000;
			}
		    skycons (rastr0,decstr0,sys0,rastr1,decstr1,sys1,32,ndec);
		    if (verbose)
			printf ("%s %s %s -> %s %s %s\n",
			    rastr0, decstr0, coorsys[sys0],
			    rastr1, decstr1, coorsys[sys1]);
		    else
			printf ("%s %s %s\n", rastr1, decstr1, coorsys[sys1]);
		    }
		}
	    else
		fprintf (stderr, "Cannot read file %s\n", listname);
	    av++;
	    }
	else if (ac > 0) {
	    strcpy (rastr0, *av);
	    ac--;
	    av++;
	    strcpy (decstr0, *av);
	    av++;

	/* Convert coordinates system to that of image */
	    if (ac > 0) {
		ac--;
		strcpy (csys, *av);
		sys0 = wcscsys (csys);
		if (sys0 >= 0) {
		    av++;
		    eqin = wcsceq (csys);
		    }
		else if (sys2 == WCS_J2000)
		    sys0 = WCS_B1950;
		else
		    sys0 = WCS_J2000;
		}
	    else if (sys2 == WCS_B1950)
		sys0 = WCS_B1950;
	    else
		sys0 = WCS_J2000;
	    if (sys2 < 0) {
		if (sys0 == WCS_J2000)
		    sys2 = WCS_B1950;
		else
		    sys2 = WCS_J2000;
		}
	    if (strlen (coorout) == 0)
		strcpy (coorout, coorsys[sys1]);
	    eqout = wcsceq (coorout);

	    /* Read second coordinates */
	    if (ac > 0) {
		strcpy (rastr1, *av);
		ac--;
		av++;
		strcpy (decstr1, *av);
		av++;

		/* Convert coordinates system to that of image */
		if (ac > 0) {
		    ac--;
		    strcpy (csys, *av);
		    sys1 = wcscsys (csys);
		    if (sys1 >= 0) {
			av++;
			eqin = wcsceq (csys);
			}
		    else if (sys0 == WCS_J2000)
			sys1 = WCS_B1950;
		    else
			sys1 = WCS_J2000;
		    }
		else if (sys0 == WCS_J2000)
		    sys1 = WCS_B1950;
		else
		    sys1 = WCS_J2000;
		if (strlen (coor1) == 0)
		    strcpy (coor1, coorsys[sys1]);
		eq1 = wcsceq (coor1);

	    skydiff (rastr0, decstr0, sys0, rastr1, decstr1, sys1, 32, ndec);
	    if (degout) {
		ra = str2ra (rastr1);
		dec = str2dec (decstr1);
		deg2str (rastr1, 32, ra, 5);
		deg2str (decstr1, 32, dec, 5);
		}
	    if (verbose)
		printf ("%s %s %s -> %s %s %s\n",
		        rastr0, decstr0, coorsys[sys0],
			rastr1, decstr1, coorout);
	    else
		printf ("%s %s %s\n", rastr1, decstr1, coorout);
	    }
	}

    return (0);
}


static void
skycons (rastr0, decstr0, sys0, rastr1, decstr1, sys1, lstr, ndec)

char	*rastr0;	/* Input right ascension */
char	*decstr0;	/* Input declination */
int	sys0;		/* Input coordinate system */
char	*rastr1;	/* Output right ascension (returned) */
char	*decstr1;	/* Output declination (returned) */
int	sys1;		/* Output coordinate system */
int	lstr;		/* Length of output strings */
int	ndec;		/* Number of decimal places in output RA seconds */
{
    double ra, dec;
    ra = str2ra (rastr0);
    dec = str2dec (decstr0);

    wcscon (sys0, sys1, eqin, eqout, &ra, &dec, epoch);

    /* Convert to B1950 FK4 */
    if (sys1 == WCS_B1950) {
	ra2str (rastr1, lstr, ra, ndec);
	dec2str (decstr1, lstr, dec, ndec-1);
	}

    /* Convert to J2000 FK5 */
    else if (sys1 == WCS_J2000) {
	ra2str (rastr1, lstr, ra, ndec);
	dec2str (decstr1, lstr, dec, ndec-1);
	}

    /* Convert to galactic coordinates */
    else if (sys1 == WCS_GALACTIC) {
	deg2str (rastr1, lstr, ra, ndec);
	deg2str (decstr1, lstr, dec, ndec);
	}

    /* Convert to ecliptic coordinates */
    else if (sys1 == WCS_ECLIPTIC) {
	deg2str (rastr1, lstr, ra, ndec);
	deg2str (decstr1, lstr, dec, ndec);
	}
    return;
}


static void
usage ()
{
    fprintf (stderr,"Find distance between coordinates\n");
    fprintf (stderr,"Usage [-bdegjv] [-y epoch] [-q system] [-n ndec] ra1 dec1 [sys1] ra2 dec2 [sys2]\n");
    fprintf (stderr,"Usage: [-vbejg] [-y epoch] [-q system] [-n ndec] @listfile\n");
    fprintf (stderr,"  -b: B1950 (FK4) output\n");
    fprintf (stderr,"  -d: Difference output in degrees\n");
    fprintf (stderr,"  -e: Ecliptic longitude and latitude output\n");
    fprintf (stderr,"  -g: Galactic longitude and latitude output\n");
    fprintf (stderr,"  -j: J2000 (FK5) output\n");
    fprintf (stderr,"  -n: Number of decimal places in output arcseconds\n");
    fprintf (stderr,"  -q: Output system, including equinox\n");
    fprintf (stderr,"  -s: Difference output in arcseconds instead of d:m:s\n");
    fprintf (stderr,"  -v: Verbose\n");
    fprintf (stderr,"  -y: Epoch of coordinates in years\n");
    exit (1);
}
/* Aug  6 1996	New program
 */
