/* File getcol.c
 * April 13, 1999
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
#include "libwcs/lwcs.h"
#include "libwcs/wcscat.h"

static void usage();

static int verbose = 0;		/* Verbose/debugging flag */
static int wfile = 0;		/* True to print output file */
static int classd = -1;		/* Guide Star Catalog object classes */
static int tabout = 0;		/* 1 for tab table to standard output */
static int debug = 0;		/* True for extra information */
static char *objname = NULL;	/* Object name for output */
static char *keyword = NULL;	/* Column to add to tab table output */
static char progpath[128];
static int ncat = 0;
static int version = 0;		/* If 1, print only program name and version */
static int match = 0;		/* If 1, match num exactly in BIN or ASC cats*/

main (ac, av)
int ac;
char **av;
{
    char *str;
    char rastr[16];
    char decstr[16];
    int i, lcat;
    char *refcatname[5];	/* reference catalog names */
    char *refcatn;
    char *temp;
    int systemp = 0;		/* Input search coordinate system */
    char *ranges = NULL;

    if (ac == 1)
        usage ();

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage ();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage ();
	}

    /* crack arguments */
    for (av++; --ac > 0; av++) {

	/* Set range and make a list of column numbers from it */
	if (strchr (*av + 1, '-') || strchr (*av + 1, ',')) {
	    if (ranges) {
		temp = ranges;
		ranges = (char *) calloc (strlen(ranges) + strlen(*av) + 2, 1);
		strcpy (ranges, temp);
		strcat (ranges, ",");
		strcat (ranges, *av);
		free (temp);
		}
	    else {
		ranges = (char *) calloc (strlen(*av) + 1, 1);
		if (strchr (*av,'.'))
		    match = 1;
		strcpy (ranges, *av);
		}
	    }

	/* Set column number */
	else if (isnum (*av)) {
	    if (strchr (*av,'.'))
		match = 1;
	    if (ranges) {
		temp = ranges;
		ranges = (char *)calloc (strlen(ranges)+strlen(*av)+2, 1);
		strcpy (ranges, temp);
		strcat (ranges, ",");
		strcat (ranges, *av);
		free (temp);
		}
	    else {
		 ranges = (char *) calloc (strlen(*av) + 1, 1);
		strcpy (ranges, *av);
		}
	    }

	/* Otherwise, read command */
	else if (*(str = *av) == '-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

	    case 'v':	/* more verbosity */
		verbose++;
		break;

	    case 'h':	/* ouput descriptive header */
		printhead++;
		break;

	    case 'n':	/* Number of lines to read */
		if (ac < 2)
		    usage ();
		nlines = atoi (*++av);
		ac--;
		break;

	    case 't':	/* tab table to stdout */
		tabout = 1;
		break;

	    default:
		usage ();
		break;
	    }
	    }
	else
	    filename = *av;
	}
    ListFile (filename, ranges);

    return (0);
}

static void
usage ()

{
    if (version)
	exit (-1);
    fprintf (stderr,"Extract specified columns from a file\n");
    fprintf (stderr,"Usage: [-tv] filenme\n");
    fprintf(stderr,"  -t: Tab table to standard output as well as file\n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}

static int
ListFile (filename, ranges)

char	*filename;	/* File name */
char	*ranges;	/* String with range of column numbers to list */

{
    int i, ngmax, nbytes;
    char line[256];
    FILE *fd;
    int nlog, closest;
    char temp[80];
    char isp[4];
    int refcat;		/* reference catalog switch */
    int icat, nndec;
    struct StarCat *starcat;

    if (verbose || printhead)

    if (debug)
	nlog = 1;
    else if (verbose)
	nlog = 100;
    else
	nlog = 0;

    if (!(refcat = RefCat (refcatname[icat], title, &sysref, &eqref, &epref))) {
	fprintf (stderr,"ListCat: Catalog '%s' is missing\n", refcatname[icat]);
	return (0);
	}

    /* Find columns specified by number */
    if (ranges != NULL) {
	int nfdef = 9;
	range = RangeInit (ranges, nfdef);
	nfind = rgetn (range);
	nbytes = nfind * sizeof (int);
	if (!(inum = (int *) calloc (nfind, sizeof(int))) {
	    fprintf (stderr, "Could not calloc %d bytes for unum\n", nbytes);
	    return;
	    }
	else {
	    for (i = 0; i < nfind; i++)
		inum[i] = rgeti4 (range);
	    }
	wfile = 0;

	if (!(fd = fopen (filename, "r"))) {
            return (0);
	    }
        }

    /* Free memory used for search results */
    if (inum) free ((char *)inum);

    return (nfind);
}

/* Feb  3 1999	New program
 * Apr 13 1999	Drop unnecessary use of progname variable
 */
