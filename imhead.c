/* File imhead.c
 * November 19, 1996
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

static void usage();
static int PrintFITSHead ();
static void PrintHead ();

static int verbose = 0;		/* verbose/debugging flag */

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;

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
	PrintHead (fn);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"Print FITS or IRAF image header\n");
    fprintf(stderr,"%s: usage: [-v] file.fit ...\n", progname);
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
PrintHead (name)

char *name;

{
    char *header;	/* FITS image header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    int *irafheader;	/* IRAF image header */
    int iraffile;

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
	fprintf (stderr,"Print World Coordinate System from ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}

    if (!PrintFITSHead (header) && verbose)
	printf ("%s: no WCS fields found\n", name);

    free (header);
    return;
}


static int
PrintFITSHead (header)

char	*header;	/* Image FITS header */
{
    char line[80], *iline, *endhead;
    int i;

    endhead = ksearch (header, "END") + 80;

    for (iline = header; iline < endhead; iline = iline + 80) {
	strncpy (line, iline, 80);
	i = 79;
	while (line[i] <= 32)
	    line[i--] = 0;
	printf ("%s\n",line);
	}

    return (1);
}
/* Jul 10 1996	New program
 * Jul 16 1996	Update header I/O
 * Aug 15 1996	Drop unnecessary reading of FITS image; clean up code
 * Aug 27 1996	Drop unused variables after lint
 * Nov 19 1996	Add linefeeds after filename in verbose mode
 */
