/* File delwcs.c
 * August 6, 1998
 * By Doug Mink, after University of Iowa code
 * (Harvard-Smithsonian Center for Astrophysics)
 * Send bug reports to dmink@cfa.harvard.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>

#include "fitsfile.h"

static void usage();
static void DelWCS ();
extern int DelWCSFITS ();

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
	DelWCS (fn);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"Delete WCS in FITS and IRAF image files\n");
    fprintf (stderr, "By D. Mink, SAO, after E. Downey, UIowa\n");
    fprintf(stderr,"%s: usage: [-v] file.fts ...\n", progname);
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
DelWCS (name)

char *name;

{
    char *header;	/* FITS image header */
    char *image;	/* Image pixels */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    int iraffile;	/* 1 if IRAF image */
    char *irafheader;	/* IRAF image header */
    char pixname[128];	/* IRAF pixel file name */

    /* Open IRAF image if .imh extension is present */
    if (strsrch (name,".imh") != NULL) {
	iraffile = 1;
	if ((irafheader = irafrhead (name, &lhead)) != NULL) {
	    if ((header = iraf2fits (name, irafheader, lhead, &nbhead))==NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",name);
		free (irafheader);
		return;
		}
	    if ((image = irafrimage (header)) == NULL) {
		hgets (header,"PIXFILE", 64, pixname);
		fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
		free (irafheader);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", name);
	    return;
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	if ((header = fitsrhead (name, &lhead, &nbhead)) != NULL) {
	    if ((image = fitsrimage (name, nbhead, header)) == NULL) {
		fprintf (stderr, "Cannot read FITS image %s\n", name);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", name);
	    return;
	    }
	}
    if (verbose) {
	fprintf (stderr,"Remove World Coordinate System from ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}

    if (DelWCSFITS (header, verbose) < 1) {
	if (verbose)
	    printf ("%s: no WCS fields found -- file unchanged\n", name);
	}
    else  {
	if (iraffile) {
	    if (irafwhead (name, lhead, irafheader, header) < 1)
		fprintf (stderr, "%s: Could not write FITS file\n", name);
	    else {
		if (verbose)
		    printf ("%s: rewritten successfully without WCS.\n", name);
		}
	    }
	else {
	    if (fitswimage (name, header, image) < 1)
		fprintf (stderr, "%s: Could not write FITS file\n", name);
	    else {
		if (verbose)
		    printf ("%s: rewritten successfully without WCS.\n", name);
		}
	    }
	}

    free (header);
    free (image);
    return;
}
/*
 * Feb 23 1996	New program split off from SETWCS
 * Apr 15 1996	Move delWCSFITS subroutine to libwcs imdelwcs.c
 * Apr 15 1996	Drop name as argument to delWCSFITS
 * May 31 1996	Rename delPos to DelWCS
 * Aug 26 1996	Change HGETC call to HGETS
 * Aug 27 1996	Fix IRAFRHEAD arguments after lint
 * Oct 16 1996	Add newlines to heading
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 *
 * Apr 14 1998	Version 2.2: deletes more parameters
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jul 24 1998	Make irafheader char instead of int
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 */
