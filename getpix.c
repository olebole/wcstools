/* File getpix.c
 * July 2, 1999
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
static void PrintPix ();

static int verbose = 0;		/* verbose/debugging flag */
static int version = 0;		/* If 1, print only program name and version */
static int nline = 10;		/* Number of pixels printer per line */
static char *pform = NULL;	/* Format in which to print pixels */
static int pixlabel = 0;	/* If 1, label pixels in output */

main (ac, av)
int ac;
char **av;
{
    char *str;
    char *fn;
    char *rrange;       /* Row range string */
    char *crange;       /* Column range string */

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

	/* output format */
	if (str[0] == '%') {
	    pform = *av;
	    }

	/* other command */
	else if (str[0] == '-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

		case 'v':	/* more verbosity */
		    verbose++;
		    break;

		case 'l':	/* label pixels */
		    pixlabel++;
		    break;

		case 'n':	/* Number of pixels per line */
		    if (ac < 2)
			usage();
		    nline = atoi (*++av);
		    ac--;
		    break;

		default:
		    usage();
		    break;
		}
	    }

	/* range of pixels to print */
        else if (isnum (str) || isrange (str)) {
	    if (crange == NULL)
		crange = str;
	    else
		rrange = str;
	    }

	/* file name */
	else
	    fn = str;

	if (fn && crange && rrange)
            PrintPix (fn, crange, rrange);
	}

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Print FITS or IRAF pixel values\n");
    fprintf(stderr,"Usage: getpix [-v][-n num][format] file.fit x_range y_range ...\n");
    fprintf(stderr,"  -n: number of pixel values printed per page\n");
    fprintf(stderr,"  -v: label pixels\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"   %: C format for each pixel value\n");
    exit (1);
}


static void
PrintPix (name, crange, rrange)

char *name;
char *crange;   /* Column range string */
char *rrange;   /* Row range string */

{
    char *header;	/* FITS image header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    char *irafheader;	/* IRAF image header */
    char *image;	/* FITS or IRAF image */
    double bzero;	/* Zero point for pixel scaling */
    double bscale;	/* Scale factor for pixel scaling */
    int iraffile;
    double dpix;
    char *c;
    int *yi;
    int bitpix,xdim,ydim, ipix, i, nx, ny, ix, iy, x, y, n, iline;
    char pixname[128];
    char nform[8];
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
	    fprintf (stderr,"Print rows %s in ", rrange);
	else if (!strcmp (rrange, "0"))
	    fprintf (stderr,"Print columns %s in ", crange);
	else
	    fprintf (stderr,"Print rows %s, columns %s in ", rrange, crange);
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}

/* Get value of specified pixel */
    hgeti4 (header,"BITPIX",&bitpix);
    xdim = 1;
    hgeti4 (header,"NAXIS1",&xdim);
    ydim = 1;
    hgeti4 (header,"NAXIS2",&ydim);
    bzero = 0.0;
    hgetr8 (header,"BZERO",&bzero);
    bscale = 1.0;
    hgetr8 (header,"BSCALE",&bscale);

/* Set format if not already set */
    if (pform == NULL) {
	pform = (char *) calloc (8,1);
	if (bitpix > 0)
	    strcpy (pform, "%d");
	else
	    strcpy (pform, "%.2f");
	}

/* Print entire columns */
    if (!strcmp (rrange, "0")) {
	xrange = RangeInit (crange, xdim);
	nx = rgetn (xrange);
	ny = ydim;
	for (ix = 0; ix < nx; ix++) {
	    rstart (xrange);
	    x = rgeti4 (xrange) - 1;
	    for (y = 0; y < ny; y++) {
        	dpix = getpix (image,bitpix,xdim,ydim,bzero,bscale,x,y);
	        if (bitpix > 0) {
		    if (dpix > 0)
	 		ipix = (int) (dpix + 0.5);
		    else if (dpix < 0)
		 	ipix = (int) (dpix - 0.5);
		    else
			ipix = 0;
		    }
		if (verbose) {
		    printf ("%s[%d,%d] = ",name,x,y);
		    if (bitpix > 0)
			printf (pform, ipix);
		    else
			printf (pform, dpix);
		    printf ("\n");
		    }
		else {
		    if (bitpix > 0)
			printf (pform, ipix);
		    else
			printf (pform, dpix);
		    if ((y+1) % nline == 0)
			printf ("\n");
		    else
			printf (" ");
		    }
		}
	    if (y % nline != 0)
		printf ("\n");
	    if (nx > 1)
		printf ("\n");
	    }
	free (xrange);
	}

/* Print entire rows */
    else if (!strcmp (crange, "0")) {
	yrange = RangeInit (rrange, xdim);
	ny = rgetn (yrange);
	nx = xdim;
	for (iy = 0; iy < ny; iy++) {
	    y = rgeti4 (yrange) - 1;
	    for (x = 0; x < nx; x++) {
        	dpix = getpix (image,bitpix,xdim,ydim,bzero,bscale,x,y);
	        if (bitpix > 0) {
		    if (dpix > 0)
	 		ipix = (int) (dpix + 0.5);
		    else if (dpix < 0)
		 	ipix = (int) (dpix - 0.5);
		    else
			ipix = 0;
		    }
		if (verbose) {
		    printf ("%s[%d,%d] = ",name,x,y);
		    if (bitpix > 0)
			printf (pform, ipix);
		    else
			printf (pform, dpix);
		    printf ("\n");
		    }
		else {
		    if (bitpix > 0)
			printf (pform, ipix);
		    else
			printf (pform, dpix);
		    if ((x+1) % nline == 0)
			printf ("\n");
		    else
			printf (" ");
		    }
		}
	    if (x % nline != 0)
		printf ("\n");
	    if (ny > 1)
		printf ("\n");
	    }
	free (yrange);
	}

/* Print a region of a two-dimensional image */
    else {
	xrange = RangeInit (crange, xdim);
	nx = rgetn (xrange);

	/* Make list of y coordinates */
	yrange = RangeInit (rrange, ydim);
	ny = rgetn (yrange);
	yi = (int *) calloc (ny, sizeof (int));
	for (i = 0; i < ny; i++) {
	    yi[i] = rgeti4 (yrange) - 1;
	    }

	/* Label horizontal pixels */
	if (!verbose && pixlabel) {
	    printf ("     ");
	    rstart (xrange);
	    strcpy (nform, pform);
	    if ((c = strchr (nform,'.')) != NULL) {
		*c = 'd';
		c[1] = (char) 0;
		}
	    else if ((c = strchr (nform,'f')) != NULL) {
		*c = 'd';
		}
	    for (ix = 0; ix < nx; ix++) {
		x = rgeti4 (xrange);
		printf (" ");
		printf (nform, x);
		}
	    printf ("\n");
	    }
	if (verbose)
	    iy = -1;
	else
	    iy = ny;

	/* Loop through rows starting with the last one */
	for (i = 0; i < ny; i++) {
	    if (verbose)
		iy++;
	    else
		iy--;
	    rstart (xrange);
	    if (!verbose && pixlabel)
		printf ("%4d: ",yi[iy]+1);

	    /* Loop through columns */
	    for (ix = 0; ix < nx; ix++) {
		x = rgeti4 (xrange) - 1;
        	dpix = getpix (image,bitpix,xdim,ydim,bzero,bscale,x,yi[iy]);
	        if (bitpix > 0) {
		    if ((c = strchr (pform,'f')) != NULL)
			*c = 'd';
		    if (dpix > 0)
	 		ipix = (int) (dpix + 0.5);
		    else if (dpix < 0)
		 	ipix = (int) (dpix - 0.5);
		    else
			ipix = 0;
		    }
		else {
		    if ((c = strchr (pform,'d')) != NULL)
			*c = 'f';
		    }
		if (verbose) {
		    printf ("%s[%d,%d] = ",name,x+1,yi[i]+1);
		    if (bitpix > 0)
			printf (pform, ipix);
		    else
			printf (pform, dpix);
		    printf ("\n");
		    }
		else {
		    if (bitpix > 0)
			printf (pform, ipix);
		    else
			printf (pform, dpix);
		    printf (" ");
		    }
		}
	    if (!verbose)
		printf ("\n");
	    }
	free (xrange);
	free (yrange);
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
 * Oct 14 1998	Use isiraf() to determine file type
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Feb 12 1999	Initialize dxisn to 1 so it works for 1-D images
 * Apr 29 1999	Add BZERO and BSCALE
 * Jun 29 1999	Fix typo in BSCALE setting
 * Jul  2 1999	Use ranges instead of individual pixels
 */
