/* File getcol.c
 * November 3, 1999
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
#include "libwcs/wcscat.h"

static void usage();

static int verbose = 0;		/* Verbose/debugging flag */
static int wfile = 0;		/* True to print output file */
static int debug = 0;		/* True for extra information */
static char *objname = NULL;	/* Object name for output */
static char *keyword = NULL;	/* Column to add to tab table output */
static char progpath[128];
static int printhead = 0;	/* If 1, print program heading */
static int ncat = 0;
static int version = 0;		/* If 1, print only program name and version */
static int nread = 0;		/* Number of lines to read (0=all) */
static int nskip = 0;		/* Number of lines to skip */

main (ac, av)
int ac;
char **av;
{
    char *str;
    char *temp;
    char *filename;
    int systemp = 0;		/* Input search coordinate system */
    char *ranges = NULL;
    int match;

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

	/* Set column number */
	if (isnum (*av)) {
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

	/* Set range and make a list of column numbers from it */
	else if (isrange (*av)) {
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

	    case 'r':	/* Number of lines to read */
		if (ac < 2)
		    usage ();
		nread = atoi (*++av);
		ac--;
		break;

	    case 's':	/* Number of lines to skip */
		if (ac < 2)
		    usage ();
		nskip = atoi (*++av);
		ac--;
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

    free (ranges);
    return (0);
}

static void
usage ()

{
    if (version)
	exit (-1);
    fprintf (stderr,"Extract specified columns from an ASCII table file\n");
    fprintf (stderr,"Usage: [-rsv] filename [column number range]\n");
    fprintf(stderr,"  -r: Number of lines to read, if not all\n");
    fprintf(stderr,"  -s: Number of lines to skip\n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}

static int
ListFile (filename, ranges)

char	*filename;	/* File name */
char	*ranges;	/* String with range of column numbers to list */

{
    int i, il, nbytes;
    char line[1024];
    FILE *fd;
    int nlog;
    struct Tokens tokens;  /* Token structure */
    struct Range *range;
    int lline;
    int nfind, ntok, nt;
    int *inum;
    char *cwhite;
    char token[80];

    cwhite = NULL;

    if (verbose || printhead)

    if (debug)
	nlog = 1;
    else if (verbose)
	nlog = 100;
    else
	nlog = 0;
    if (nread < 1)
	nread = 100000;

    /* Open input file */
    if (!strcmp (filename, "stdin"))
	fd = stdin;
    else if (!(fd = fopen (filename, "r"))) {
        return (0);
	}

    /* Skip lines into input file */
    if (nskip > 0) {
	for (i = 0; i < nskip; i++) {
	    if (fgets (line, 80, fd) == NULL)
		break;
	    }
	}

    /* Print entire selected lines */
    if (ranges == NULL) {
	for (il = 0; il < nread; il++) {
	    if (fgets (line, 1024, fd) == NULL)
		break;
	    lline = strlen (line);
	    if (line[lline-1] < 32)
		line[lline-1] = 0;
	    printf ("%s\n", line);
	    }
	}

    /* Find columns specified by number */
    else {
	int nfdef = 9;
	range = RangeInit (ranges, nfdef);
	nfind = rgetn (range);
	nbytes = nfind * sizeof (int);
	if (!(inum = (int *) calloc (nfind, sizeof(int))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for unum\n", nbytes);
	    return;
	    }
	else {
	    for (i = 0; i < nfind; i++)
		inum[i] = rgeti4 (range);
	    }
	wfile = 0;
	for (il = 0; il < nread; il++) {
	    if (fgets (line, 1024, fd) == NULL)
		break;

	    ntok = setoken (&tokens, line, cwhite);
	    nt = 0;
	    for (i = 0; i < nfind; i++) {
		if (getoken (tokens, inum[i], token)) {
		    if (i < nfind-1)
			printf ("%s ", token);
		    else
			printf ("%s", token);
		    nt++;
		    }
		}
	    if (nt > 0)
		printf ("\n");
	    }
        }

    /* Free memory used for search results */
    if (inum) free ((char *)inum);

    return (nfind);
}

/* Nov  2 1999	New program
 * Nov  3 1999	Add option to read from stdin as input filename
 */
