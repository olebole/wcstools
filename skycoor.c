/* File skycoor.c
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
#include "libwcs/fitshead.h"

static void usage();
static void skycons();

static int verbose = 0;		/* verbose/debugging flag */

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
    char sys0[16], sys1[16];
    double ra, dec;
    int degout = 0;
    int ndec = 3;		/* Number of decimal places in RA seconds */
    int ndecset = 0;

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
	    strcpy (sys1,"B1950");
    	    break;
	case 'd':
	    degout++;
	    if (!ndecset)
		ndec = 5;
    	    break;
	case 'g':
	    strcpy (sys1,"galactic");
	    if (!ndecset)
		ndec = 5;
    	    break;
	case 'j':
	    strcpy (sys1,"J2000");
    	    break;
	case 'n':
	    if (ac < 2)
		usage();
	    ndec = atoi (*++av);
	    ndecset++;
	    ac--;
    	    break;
    	default:
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
		    sys0[0] = 0;
		    sscanf (line,"%s %s %s", rastr0, decstr0, sys0);
	    	    if (sys0[0] == 0) {
			if (!strcmp (sys1,"J2000"))
			    strcpy (sys0, "B1950");
			else
			    strcpy (sys0, "J2000");
			}
		    if (sys0[0] == 0)
			strcpy (sys0, "B1950");
		    skycons (rastr0, decstr0, sys0, rastr1, decstr1, sys1,ndec);
		    if (verbose)
			printf ("%s %s %s -> %s %s %s\n",
			    rastr0, decstr0, sys0, rastr1, decstr1, sys1);
		    else
			printf ("%s %s %s\n", rastr1, decstr1, sys1);
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
		strcpy (sys0, *av);
		av++;
		}
	    else if (!strcmp (sys1,"J2000"))
		strcpy (sys0, "B1950");
	    else
		strcpy (sys0, "J2000");
	    skycons (rastr0, decstr0, sys0, rastr1, decstr1, sys1, ndec);
	    if (degout) {
		ra = str2ra (rastr1);
		dec = str2dec (decstr1);
		deg2str (rastr1, ra, 5);
		deg2str (decstr1, dec, 5);
		}
	    if (verbose)
		printf ("%s %s %s -> %s %s %s\n",
		        rastr0, decstr0, sys0, rastr1, decstr1, sys1);
	    else
		printf ("%s %s %s\n", rastr1, decstr1, sys1);
	    }
	}

    return (0);
}


static void
skycons (rastr0, decstr0, sys0, rastr1, decstr1, sys1, ndec)

char *rastr0;	/* Input right ascension */
char *decstr0;	/* Input declination */
char *sys0;	/* Input coordinate system */
char *rastr1;	/* Output right ascension (returned) */
char *decstr1;	/* Output declination (returned) */
char *sys1;	/* Output coordinate system */
int ndec;	/* Number of decimal places in output RA seconds */
{
    double ra, dec;
    ra = str2ra (rastr0);
    dec = str2dec (decstr0);

    /* Convert to B1950 FK4 */
    if (sys1[0] == 'B' || sys1[0] == 'b' ||
	!strcmp (sys1,"FK4") || !strcmp (sys1, "fk4")) {
	if (sys0[0] == 'J' || sys0[0] == 'j' ||
	    !strcmp (sys0,"FK5") || !strcmp (sys0, "fk5"))
	    fk524 (&ra, &dec);
	else if (sys0[0] == 'G' || sys0[0] == 'g')
	    gal2fk4 (&ra, &dec);
	ra2str (rastr1,ra, ndec);
	dec2str (decstr1, dec, ndec-1);
	}

    /* Convert to J2000 FK5 */
    else if (sys1[0] == 'J' || sys1[0] == 'j' ||
	!strcmp (sys1,"FK5") || !strcmp (sys1, "fk5")) {
	if (sys0[0] == 'B' || sys0[0] == 'b' ||
	    !strcmp (sys0,"FK4") || !strcmp (sys0, "fk4"))
	    fk425 (&ra, &dec);
	else if (sys0[0] == 'G' || sys0[0] == 'g')
	    gal2fk5 (&ra, &dec);
	ra2str (rastr1,ra, ndec);
	dec2str (decstr1, dec, ndec-1);
	}

    /* Convert to galactic coordinates */
    else if (sys1[0] == 'G' || sys1[0] == 'g') {
	if (sys0[0] == 'B' || sys0[0] == 'b' ||
	    !strcmp (sys0,"FK4") || !strcmp (sys0, "fk4"))
	    fk42gal (&ra, &dec);
	else if (sys0[0] == 'J' || sys0[0] == 'j' ||
	    !strcmp (sys0,"FK5") || !strcmp (sys0, "fk5"))
	    fk52gal (&ra, &dec);
	deg2str (rastr1, ra, ndec);
	deg2str (decstr1, dec, ndec);
	}
    return;
}


static void
usage ()
{
    fprintf (stderr,"Convert coordinates\n");
    fprintf (stderr,"Usage [-bdgjv] [-n ndec] ra1 dec1 sys1 ... ran decn sysn\n");
    fprintf (stderr,"Usage: [-vbjg] [-n ndec] @listfile\n");
    fprintf (stderr,"  -b: B1950 (FK4) output\n");
    fprintf (stderr,"  -d: RA and Dec output in degrees\n");
    fprintf (stderr,"  -g: galactic longitude and latitude output\n");
    fprintf (stderr,"  -j: J2000 (FK5) output\n");
    fprintf (stderr,"  -n: number of decimal places in output RA seconds\n");
    fprintf (stderr,"  -v: verbose\n");
    exit (1);
}
/* Oct 30 1996	New program
 * Nov  1 1996	Use DEG2STR to get rounded degrees to output appropriately
 */
