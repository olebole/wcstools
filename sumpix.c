/* File sumpix.c
 * July 6, 1999
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
#include "wcscat.h"

static void usage();
static int PrintFITSHead ();
static void SumPix ();

static int verbose = 0;		/* verbose/debugging flag */
static int version = 0;	/* If 1, print only program name and version */

main (ac, av)
int ac;
char **av;
{
    char *str;
    char *fn;
    char *rrange;	/* Row range string */
    char *crange;	/* Column range string */

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    crange = NULL;
    rrange = NULL;
    fn = NULL;

    /* crack arguments */
    for (av++; --ac > 0; av++) {
	str = *av;
	if (str[0] == '-') {
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

	/* Read range in x or y */
        else if (isnum (*av) || isrange (*av)) {
	    if (crange == NULL)
		crange = *av;
	    else
		rrange = *av;
	    }

	/* read filename */
	else
	    fn = *av;

	if (fn && crange && rrange)
            SumPix (fn, crange, rrange);
        }

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Sum row, column, or region of a FITS or IRAF image\n");
    fprintf(stderr,"Usage: sumpix [-v] x_range  y_range file.fit ...\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  a range of 0 implies the full dimension\n");
    exit (1);
}


static void
SumPix (name, crange, rrange)

char *name;	/* FITS or IRAF .imh file name */
char *crange;	/* Column range string */
char *rrange;	/* Row range string */

{
    char *header;	/* FITS image header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    char *irafheader;	/* IRAF image header */
    char *image;	/* FITS or IRAF image */
    double bzero;	/* Zero point for pixel scaling */
    double bscale;	/* Scale factor for pixel scaling */
    int iraffile;
    double dpix, sum;
    int bitpix,xdim,ydim, ipix, i;
    int nx, ny, ix, iy, x, y;
    char pixname[128];
    struct Range *xrange;    /* X range structure */
    struct Range *yrange;    /* Y range structure */

    /* Open IRAF image if .imh extension is present */
    if (isiraf (name)) {
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
	if (!strcmp (crange, "0"))
	    fprintf (stderr,"Sum rows %s in ", rrange);
	else if (!strcmp (rrange, "0"))
	    fprintf (stderr,"Sum columns %s in ", crange);
	else
	    fprintf (stderr,"Sum rows %s, columns %s in ", rrange, crange);
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}

/* Get information about the image */
    hgeti4 (header,"BITPIX",&bitpix);
    xdim = 1;
    hgeti4 (header,"NAXIS1",&xdim);
    ydim = 1;
    hgeti4 (header,"NAXIS2",&ydim);
    bzero = 0.0;
    hgetr8 (header,"BZERO",&bzero);
    bscale = 1.0;
    hgetr8 (header,"BSCALE",&bscale);

    /* Sum entire columns */
    if (!strcmp (rrange, "0")) {
	xrange = RangeInit (crange, xdim);
	nx = rgetn (xrange);
	ny = ydim;
	for (ix = 0; ix < nx; ix++) {
	    x = rgeti4 (xrange) - 1;
	    sum = 0.0;
	    for (y = 0; y < ny; y++) {
        	dpix = getpix (image,bitpix,xdim,ydim,bzero,bscale,x,y);
		sum = sum + dpix;
		}
	    printf ("%f\n", sum);
	    }
	free (xrange);
	}

    /* Sum entire rows */
    else if (!strcmp (crange, "0")) {
	yrange = RangeInit (rrange, xdim);
	ny = rgetn (yrange);
	nx = xdim;
	for (iy = 0; iy < ny; iy++) {
	    y = rgeti4 (yrange) - 1;
	    sum = 0.0;
	    for (x = 0; x < nx; x++) {
        	dpix = getpix (image,bitpix,xdim,ydim,bzero,bscale,x,y);
		sum = sum + dpix;
		}
	    printf ("%f\n", sum);
	    }
	free (yrange);
	}

    /* Sum a region of a two-dimensional image */
    else {
	yrange = RangeInit (rrange, ydim);
	ny = rgetn (yrange);
	sum = 0.0;
	for (iy = 0; iy < ny; iy++) {
	    y = rgeti4 (yrange) - 1;
	    xrange = RangeInit (crange, xdim);
	    nx = rgetn (xrange);
	    for (ix = 0; ix < nx; ix++) {
		x = rgeti4 (xrange) - 1;
        	dpix = getpix (image,bitpix,xdim,ydim,bzero,bscale,x,y);
		sum = sum + dpix;
		if (verbose)
		    printf ("%s[%d,%d] = %f\n",name,x,y,dpix);
		}
	    }
	printf ("%f\n", sum);
	free (xrange);
	free (yrange);
	}

    free (header);
    free (image);
    return;
}
/* Jul  2 1999	New program
 * Jul  6 1999	Fix bug with x computation in patch adding section
 */
