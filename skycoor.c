/* File skycoor.c
 * February 1, 2000
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
#include "libwcs/wcscat.h"
#include "libwcs/fitsfile.h"

static void usage();
static void skycons();
extern void s2v3();
extern void v2s3();

static int verbose = 0;		/* verbose/debugging flag */
static double epoch = 0.0;
static double eqout = 0.0;
static double eqin = 0.0;
static int version = 0;		/* If 1, print only program name and version */


main (ac, av)
int ac;
char **av;
{
    char *str;
    FILE *fd;
    char *ln, *listname;
    char line[80];
    char rastr0[32], decstr0[32];
    char rastr1[32], decstr1[32];
    char csys0[32], csys1[32];
    char csys[32];
    int sys0;
    int sys1 = -1;
    double ra, dec, r, ra1, dec1;
    int lstr = 32;
    int degout = 0;
    int ndec = 3;		/* Number of decimal places in RA seconds */
    int ndecset = 0;
    char coorout[16];
    double pos[3];

    listname = NULL;
    coorout[0] = (char) 0;

    /* There are ac arguments starting at av[0] */
    if (ac == 1)
	usage();

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    /* Decode arguments */
    for (av++; --ac > 0 && (*(str = *av) == '-'); av++) {
	char c;
	while (c = *++str)
    	switch (c) {

    	case 'v':	/* more verbosity */
    	    verbose++;
    	    break;

	case 'b':
	    sys1 = WCS_B1950;
    	    break;

	case 'd':
	    degout++;
	    if (!ndecset)
		ndec = 5;
    	    break;

	case 'e':
	    sys1 = WCS_ECLIPTIC;
	    if (!ndecset)
		ndec = 5;
    	    break;

	case 'g':
	    sys1 = WCS_GALACTIC;
	    if (!ndecset)
		ndec = 5;
    	    break;

	case 'j':
	    sys1 = WCS_J2000;
    	    break;

	case 'q':
	    if (ac < 2)
		usage();
	    strcpy (coorout, *++av);
	    ac--;
    	    break;

	case 'r':
	    if (ac < 5)
		usage();
	    ra = str2ra (*++av);
	    ac--;
	    dec = str2dec (*++av);
	    ac--;
	    ra1 = str2ra (*++av);
	    ac--;
	    dec1 = str2dec (*++av);
	    ac--;
	    ra2str (rastr0, lstr, ra, 3);
	    dec2str (decstr0, lstr, dec, 2);
	    printf ("ra1, dec1: %s %s\n",
		    rastr0, decstr0);
	    ra2str (rastr0, lstr, ra1, 3);
	    dec2str (decstr0, lstr, dec1, 2);
	    printf ("ra2, dec2: %s %s\n",
		    rastr0, decstr0);
	    r = wcsdist (ra, dec, ra1, dec1);
	    if (r < 1.0) {
		r = r * 3600.0;
		printf ("Distance is %.3f arcsec\n", r);
		}
	    else {
		printf ("Distance is %.5f degrees\n", r);
		}
	    break;

	case 'w':
	    if (ac < 3)
		usage();
	    ra = str2ra (*++av);
	    ac--;
	    dec = str2dec (*++av);
	    ac--;
	    r = 1.0;
	    s2v3 (degrad(ra), degrad(dec), r, pos);
	    printf (" x,y,z:  %.6f %.6f %.6f\n",
		    pos[0],pos[1],pos[2]);
	    break;

	case 'x':
	    if (ac < 4)
		usage();
	    pos[0] = atof (*++av);
	    ac--;
	    pos[1] = atof (*++av);
	    ac--;
	    pos[2] = atof (*++av);
	    ac--;
	    v2s3 (pos, &ra, &dec, &r);
	    ra2str (rastr0, lstr, raddeg(ra), 3);
	    dec2str (decstr0, lstr, raddeg(dec), 2);
	    printf ("ra, dec, r: %s %s %.6f\n",
		    rastr0, decstr0, r);
	    break;

	case 'y':	/* Epoch of coordinates */
	    if (ac < 2)
		usage();
	    epoch = fd2ep (*++av);
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
    		usage();
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
			if (sys1 == WCS_J2000)
			    sys0 = WCS_B1950;
			else
			    sys0 = WCS_J2000;
			}
		    else
			sys0 = wcscsys (csys);
		    if (sys1 < 0) {
			if (degout)
			    sys1 = sys0;
			else if (sys0 == WCS_J2000)
			    sys1 = WCS_B1950;
			else
			    sys1 = WCS_J2000;
			}
		    skycons (rastr0,decstr0,sys0,rastr1,decstr1,sys1,lstr,ndec);
		    wcscstr (csys0, sys0, 0.0, 0.0);
		    wcscstr (csys1, sys1, 0.0, 0.0);
		    if (degout) {
			ra = str2ra (rastr1);
			dec = str2dec (decstr1);
			deg2str (rastr1, 32, ra, 5);
			deg2str (decstr1, 32, dec, 5);
			}
		    if (verbose)
			printf ("%s %s %s -> %s %s %s\n",
			    rastr0, decstr0, csys0,
			    rastr1, decstr1, csys1);
		    else
			printf ("%s %s %s\n", rastr1, decstr1, csys1);
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
		av++;
		sys0 = wcscsys (csys);
		eqin = wcsceq (csys);
		}
	    else if (sys1 == WCS_J2000)
		sys0 = WCS_B1950;
	    else
		sys0 = WCS_J2000;
	    if (sys1 < 0) {
		if (degout)
		    sys1 = sys0;
		else if (sys0 == WCS_J2000)
		    sys1 = WCS_B1950;
		else
		    sys1 = WCS_J2000;
		}
	    if (strlen (coorout) == 0)
		wcscstr (coorout, sys1, 0.0, 0.0);
	    eqout = wcsceq (coorout);
	    skycons (rastr0, decstr0, sys0, rastr1, decstr1, sys1, lstr, ndec);
	    if (degout) {
		ra = str2ra (rastr1);
		dec = str2dec (decstr1);
		deg2str (rastr1, 32, ra, 5);
		deg2str (decstr1, 32, dec, 5);
		}
	    if (verbose) {
		wcscstr (csys0, sys0, 0.0, 0.0);
		printf ("%s %s %s -> %s %s %s\n",
		        rastr0, decstr0, csys0,
			rastr1, decstr1, coorout);
		}
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
    if (version)
	exit (-1);
    fprintf (stderr,"Convert coordinates\n");
    fprintf (stderr,"Usage [-bdegjv] [-y epoch] [-q system] [-n ndec] [-r RA Dec RA Dec]\n");
    fprintf (stderr,"      [-w RA Dec] [-x x y z] ra1 dec1 sys1 ... ran decn sysn\n");
    fprintf (stderr,"      [ra1 dec1 sys1 ... ran decn sysn] or [@listfile]\n");
    fprintf (stderr,"  -b: B1950 (FK4) output\n");
    fprintf (stderr,"  -d: RA and Dec output in degrees\n");
    fprintf (stderr,"  -e: Ecliptic longitude and latitude output\n");
    fprintf (stderr,"  -g: Galactic longitude and latitude output\n");
    fprintf (stderr,"  -j: J2000 (FK5) output\n");
    fprintf (stderr,"  -n: Number of decimal places in output RA seconds\n");
    fprintf (stderr,"  -q: Output system, including equinox\n");
    fprintf (stderr,"  -r: Angular separation between two RA, Dec pairs\n");
    fprintf (stderr,"  -v: verbose\n");
    fprintf (stderr,"  -w: Convert RA, Dec equatorial coordinates to x,y,z\n");
    fprintf (stderr,"  -x: Convert x,y,z equatorial coordinates to RA, Dec\n");
    fprintf (stderr,"  -y: Epoch of coordinates in years\n");
    exit (1);
}
/* Oct 30 1996	New program
 * Nov  1 1996	Use DEG2STR to get rounded degrees to output appropriately
 *
 * Feb 21 1997	Add option to input epoch of coordinates
 *
 * Apr 14 1998	Add ecliptic longitude, latitude output
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * May  1 1998	Increase coordinate system name array from 8 to 12 characters
 * May 13 1998	Add q command for output equinox; allow input equinox, too
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Aug  6 1998	Do not include fitshead.h; it is in wcs.h
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Apr 16 1999	Add xyz <-> RA/DEC conversions using w and x arguments
 * Jul  1 1999	Allow any legal FITS date format for epoch
 * Jul  8 1999	Fix bug in computing difference in arcseconds
 * Oct 22 1999	Drop unused variables after lint; fix 2 == bugs
 * Nov 29 1999	Include fitsfile.h for date conversion
 *
 * Feb  1 2000	If degrees output, assume same system unless told otherwise
 */
