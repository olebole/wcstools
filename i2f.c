/* File i2f.c
 * October 17, 1996
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
static void IRAFtoFITS ();

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
	while (c = *++str) {
	    switch (c) {
		case 'v':	/* more verbosity */
		    verbose++;
		    break;
	        default:
		    usage(progname);
		    break;
		}
    	    }
	}

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage (progname);

    else {
	while (ac-- > 0) {
	    char *fn = *av++;
	    if (verbose)
		printf ("%s:\n", fn);
	    IRAFtoFITS (fn);
	    if (verbose)
		printf ("\n");
	    }
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"Write FITS files from IRAF image files\n");
    fprintf(stderr,"%s: usage: [-v] file.imh ...\n", progname);
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}

static void
IRAFtoFITS (name)
char *name;
{
    char *image;	/* FITS image */
    char *header;	/* FITS header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    int *irafheader;	/* IRAF image header */
    char irafname[64];	/* Name of IRAF file */
    char fitsname[64];	/* Name of FITS file */
    char *ext;		/* Pointer to start of extension */
    char pixname[128];	/* Pixel file name */

    /* Open IRAF image if .imh extension is present */
    if (strsrch (name,".imh") != NULL) {
	irafheader = irafrhead (name, &lhead);
	if (irafheader) {
	    header = iraf2fits (name, irafheader, lhead, &nbhead);
	    free (irafheader);
	    if (header == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",name);
		return;
		}
	    image = irafrimage (header);
	    if (image == NULL) {
		hgets (header,"PIXFILE", 64, pixname);
		fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
		free (irafheader);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr,"Cannot read IRAF header file %s\n", name);
	    return;
	    }
	strcpy (fitsname, name);
	ext = strsrch (fitsname,".imh");
	strcpy (ext,".fit");
	if (verbose) {
	    fprintf (stderr,"Write FITS files from IRAF image file %s\n", name);
	    }
	}

    /* Add .imh extension to make IRAF header file name if not present */
    else {
	strcpy (irafname, name);
	strcat (irafname,".imh");
	irafheader = irafrhead (irafname, &lhead);
	if (irafheader) {
	    header = iraf2fits (irafname, irafheader, lhead, &nbhead);
	    free (irafheader);
	    if (header == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",irafname);
		return;
		}
	    image = irafrimage (header);
	    if (image == NULL) {
		hgets (header,"PIXFILE",64, pixname);
		fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
		free (irafheader);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr,"Cannot read IRAF header file %s\n", irafname);
	    return;
	    }
	strcpy (fitsname, name);
	strcat (fitsname,".fit");
	if (verbose) {
	    fprintf (stderr,"Write FITS files from IRAF image file %s\n", irafname);
	    }
	}

    /* Write FITS image */
    if (fitswimage (fitsname, header, image) > 0 && verbose)
	printf ("%s: rewritten successfully.\n", fitsname);

    else if (verbose)
	printf ("%s: not written.\n", fitsname);

    free (header);
    free (image);
    return;
}
/* Jun  6 1996	New program
 * Jul 16 1996	Update header input
 * Aug 16 1996	Clean up code
 * Aug 26 1996	Change HGETC call to HGETS
 * Aug 27 1996	Drop unused variables after lint
 * Oct 17 1996	Clean up after lint
 */
