/* File getpix.c
 * August 6, 1998
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
#include "fitsfile.h"

static void usage();
static int PrintFITSHead ();
static void PrintPix ();

static int verbose = 0;		/* verbose/debugging flag */

main (ac, av)
int ac;
char **av;
{
    char *str;
    char *fn;
    int i, x[100], y[100];

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {
	case 'v':	/* more verbosity */
	    verbose++;
	    break;
	default:
	    usage();
	    break;
	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage ();

    fn = *av++;

    i = 0;
    while (--ac > 0 && i < 100) {
	x[i] = atoi (*av++);
	ac--;
	y[i] = atoi (*av++);
	i++;
	}
    if (i > 0)
	PrintPix (fn, i, x, y);

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"Print FITS or IRAF pixel values\n");
    fprintf(stderr,"Usage: getpix [-v] file.fit x y ...\n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
PrintPix (name, n, x, y)

char *name;
int n, *x, *y;

{
    char *header;	/* FITS image header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    char *irafheader;	/* IRAF image header */
    char *image;	/* FITS or IRAF image */
    int iraffile;
    double dpix;
    int bitpix,xdim,ydim, ipix, i;
    char pixname[128];

    /* Open IRAF image if .imh extension is present */
    if (strsrch (name,".imh") != NULL) {
	iraffile = 1;
	if ((irafheader = irafrhead (name, &lhead)) != NULL) {
	    header = iraf2fits (name, irafheader, lhead, &nbhead);
	    free (irafheader);
	    if (header == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",name);
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
	    fprintf (stderr, "Cannot read IRAF file %s\n", name);
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
	fprintf (stderr,"Print pixel values from ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}

/* Get value of specified pixel */
    hgeti4 (header,"BITPIX",&bitpix);
    hgeti4 (header,"NAXIS1",&xdim);
    hgeti4 (header,"NAXIS2",&ydim);

    for (i = 0; i < n; i++) {
        dpix = getpix (image, bitpix, xdim, ydim, x[i]-1, y[i]-1);
        if (bitpix > 0) {
	    if (dpix > 0)
        	ipix = (int) (dpix + 0.5);
	    else if (dpix < 0)
        	ipix = (int) (dpix - 0.5);
	    else
		ipix = 0;
	    printf ("%s[%d,%d] = %d\n",name,x[i],y[i],ipix);
	    }
	else
	    printf ("%s[%d,%d] = %.2f\n",name,x[i],y[i],dpix);
	}

    free (header);
    free (image);
    return;
}
/* Dec  6 1996	New program
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 * Dec 15 1997	Add capability of reading and writing IRAF 2.11 images
 *
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jul 24 1998	Make irafheader char instead of int
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 */
