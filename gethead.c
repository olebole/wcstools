/* File gethead.c
 * October 8, 1996
 * By Doug Mink Harvard-Smithsonian Center for Astrophysics)
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

#define MAXKWD 50

static void usage();
static void PrintValues ();

static int verbose = 0;		/* verbose/debugging flag */

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    char *kwd[MAXKWD];
    int nkwd;
    char *fn;

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

    if (ac-- > 0) {
	fn = *av++;
	}
    nkwd = 0;
    while (ac-- > 0 && nkwd < MAXKWD) {
	kwd[nkwd] = *av++;
	nkwd++;
	}
    if (nkwd > 0)
	PrintValues (fn, nkwd, kwd);

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"Print FITS or IRAF header keyword values\n");
    fprintf(stderr,"%s: usage: [-v] file.fit kw1 kw2 ... kwn\n", progname);
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
PrintValues (name, nkwd, kwd)

char	*name;	/* Name of FITS or IRAF image file */
int	nkwd;	/* Number of keywords for which to print values */
char	*kwd[];	/* Names of keywords for which to print values */

{
    char *header;	/* FITS image header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    int *irafheader;	/* IRAF image header */
    int iraffile;
    char string[80];
    char *kw, *kwe;
    int ikwd, lkwd;

    /* Open IRAF image if .imh extension is present */
    if (strsrch (name,".imh")) {
	iraffile = 1;
	if ((irafheader = irafrhead (name, &lhead))) {
	    header = iraf2fits (name, irafheader, lhead, &nbhead);
	    free (irafheader);
	    if (!header) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",name);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF file %s\n", name);
	    return;
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	if (!(header = fitsrhead (name, &lhead, &nbhead))) {
	    fprintf (stderr, "Cannot read FITS file %s\n", name);
	    return;
	    }
	}
    if (verbose) {
	fprintf (stderr,"Print Header Parameter Values from ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}

    /* Get keyword values one at a time */
    for (ikwd = 0; ikwd < nkwd; ikwd++) {
	lkwd = strlen (kwd[ikwd]);
	kwe = kwd[ikwd] + lkwd;
	for (kw = kwd[ikwd]; kw < kwe; kw++) {
	    if (*kw > 96 && *kw < 123)
		*kw = *kw - 32;
	    }
	if (hgets (header, kwd[ikwd], 80, string)) {
	    if (verbose)
		printf ("%s = %s", kwd[ikwd], string);
	    else
		printf ("%s",string);
	    }
	else if (verbose)
	    printf ("%s not found", kwd[ikwd]);
	else
	    printf ("___");

	if (verbose || ikwd == nkwd - 1)
	    printf ("\n");
	else
	    printf (" ");
	}

    free (header);
    return;
}

/* Sep  4 1996	New program
 * Oct  8 1996	Add newline after file name in verbose mode
 */
