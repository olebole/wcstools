/* File delwcs.c
 * February 23, 1996
 * By Elwood Downey, revised by Doug Mink
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>

#include "fitshead.h"

#define MAXHEADLEN 14400

static void usage();
static void delPos ();

static int verbose = 0;		/* verbose/debugging flag */
static char *RevMsg = "DELWCS version 1.0, 23 February 1996";

extern int setWCSFITS ();

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    double rot, scale, gsclim, frac;
    int tolerance;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {
	case 'v':	/* more verbosity */
	    verbose++;
	    break;
	default:
	    usage(progname);
	    break;
	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage (progname);

    while (ac-- > 0) {
	char *fn = *av++;
	delPos (fn);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"%s\n",RevMsg);
    fprintf (stderr,"Delete WCS in FITS and IRAF image files\n");
    fprintf (stderr, "By E. Downey, UIowa and D. Mink, SAO\n");
    fprintf(stderr,"%s: usage: [-v] file.fts ...\n", progname);
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
delPos (name)

char *name;

{
    char *header;
    char *image;
    int lhead;
    int iraffile;		/* 1 if IRAF image */
    int *irafheader;		/* IRAF image header */

    /* Allocate FITS header */
    header = malloc (MAXHEADLEN);
    lhead = MAXHEADLEN;

    /* Open IRAF image if .imh extension is present */
    if (strsrch (name,".imh") != NULL) {
	iraffile = 1;
	irafheader = irafrhead (name, lhead, header);
	if (irafheader == NULL) {
	    fprintf (stderr, "Cannot read header file for %s\n", name);
	    free (header);
	    return;
	    }
	image = irafrimage (name, irafheader, header);
	if (image == NULL) {
	    fprintf (stderr, "Cannot read pixel file for %s\n", name);
	    free (header);
	    free (irafheader);
	    return;
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	image = fitsrimage (name, lhead, header);
	if (image == NULL) {
	    fprintf (stderr, "Cannot read file %s\n", name);
	    free (header);
	    return;
	    }
	}
    if (verbose) {
	fprintf (stderr,"%s\n",RevMsg);
	fprintf (stderr,"Remove World Coordinate System from FITS or IRAF image files\n");
	fprintf (stderr, "By E. Downey, UIowa and D. Mink, SAO\n");
	}

    if (delWCSFITS (name, header, verbose) < 1) {
	if (verbose)
	    printf ("%s: no WCS fields found -- file unchanged\n", name);
	}
    else  {
	if (iraffile) {
	    if (irafwhead (name, irafheader, header) < 1)
		fprintf (stderr, "%s: Could not write FITS file\n", name);
	    else if (verbose)
		printf ("%s: rewritten successfully without WCS.\n", name);
	    }
	else {
	    if (fitswimage (name, header, image) < 1)
		fprintf (stderr, "%s: Could not write FITS file\n", name);
	    else if (verbose)
		printf ("%s: rewritten successfully without WCS.\n", name);
	    }
	}

    free (header);
    free (image);
    return;
}


/* delete all the C* fields.
 * return 0 if at least one such field is found, else -1.  */

int
delWCSFITS (name, header, verbose)

char *name;
char *header;
int verbose;

{
    static char *flds[] = {
	"CTYPE1", "CRVAL1", "CDELT1", "CRPIX1", "CROTA1",
	"CTYPE2", "CRVAL2", "CDELT2", "CRPIX2", "CROTA2" };
    int i;
    int n;

    n = 0;

    for (i = 0; i < sizeof(flds)/sizeof(flds[0]); i++) {
	if (hdel (header, flds[i])) {
	    n++;
	    if (verbose)
		printf ("%s: deleted\n", flds[i]);
	    }
	else if (verbose)
	    printf ("%s: not found\n", flds[i]);
	}

    if (ksearch (header, "EPOCH")) {
	if (hdel (header,"EQUINOX")) {
	    if (verbose)
		printf ("EQUINOX: deleted\n");
	    n++;
	    }
	else
	    printf ("%s: EPOCH, but not EQUINOX found\n");
	}

    return (n);
}
