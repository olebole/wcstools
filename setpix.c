/* File setpix.c
 * September 14, 1999
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
#include "fitsfile.h"
#include "wcs.h"
#include "wcscat.h"

static void usage();
static void SetPix();

static int newimage = 0;
static int verbose = 0;		/* verbose flag */
static int eachpix = 0;		/* If 1, print each pixel change */
static int version = 0;		/* If 1, print only program name and version */
static char *pform = NULL;	/* Format in which to print pixels */


main (ac, av)
int ac;
char **av;
{
    char *str;
    char *fn[100];
    char *crange[100], *rrange[100];
    char *value[100];
    char *listfile = NULL;
    char nextline[128];
    FILE *flist;
    int i = 0;
    int iv = 0;

    fn[0] = NULL;
    value[0] = NULL;
    rrange[0] = NULL;
    crange[0] = NULL;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    /* crack arguments */
    for (av++; --ac > 0; av++) {
	str = *av;

	if (str[0] == '-') {
	    char c;
	    while (c = *++str) {
		switch (c) {
		    case 'i':	/* print each change */
			eachpix++;
			break;
		    case 'n':	/* ouput new file */
			newimage++;
			break;
		    case 'v':	/* some verbosity */
			verbose++;
			break;
		    default:
			usage ();
			break;
		    }
		}
	    }

	/* Set up pixel changes from a file */
	else if (str[0] == '@') {
	    listfile = str+1;
	    if ((flist = fopen (listfile, "r")) == NULL) {
		fprintf (stderr,"CONPIX: List file %s cannot be read\n",
			 listfile);
		usage ();
		}
	    while (fgets (nextline, 64, flist) != NULL) { 
		crange[iv] = (char *) malloc (16);
		rrange[iv] = (char *) malloc (16);
		value[iv] = (char *) malloc (16);
		sscanf (nextline, "%s %s %s",
			crange[iv], rrange[iv], value[iv]);
		iv++;
		rrange[iv] = NULL;
		crange[iv] = NULL;
		value[iv] = NULL;
		}
	    }

	/* Set x or y range or new pixel value */
	else if (isnum (str) || isrange (str)) {
	    if (crange[iv] == NULL)
		crange[iv] = str;
	    else if (rrange[iv] == NULL)
		rrange[iv] = str;
	    else if (isnum (str)) {
		value[iv] = str;
		iv++;
		rrange[iv] = NULL;
		crange[iv] = NULL;
		value[iv] = NULL;
		}
	    }
	else {
	    fn[i] = str;
	    i++;
	    fn[i] = NULL;
	    }
	}

    if (i == 0 || iv == 0)
	usage();

    /* Loop through pixel changes for each image */
    i = 0;
    while (fn[i] != NULL) {
	SetPix (fn[i], crange, rrange, value);
	i++;
	}

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Edit pixels of FITS or IRAF image file\n");
    fprintf(stderr,"Usage: setpix [-vn] file.fts x_range y_range value ...\n");
    fprintf(stderr,"Usage: setpix [-vn] file.fts @valuefile ...\n");
    fprintf(stderr,"  -n: write new file, else overwrite \n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
SetPix (filename, crange, rrange, value)

char	*filename;	/* FITS or IRAF file filename */
char	**crange;	/* Column range string */
char	**rrange;	/* Row range string */
char	**value;	/* value to insert into pixel */

{
    char *image;		/* FITS image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    int iraffile;		/* 1 if IRAF image */
    char *irafheader;		/* IRAF image header */
    int i, nbytes, nhb, nhblk, lname, lext, lroot;
    int nx, ny, ix, iy, x, y, ipix;
    char *head, *imext, *imext1;
    double bzero;		/* Zero point for pixel scaling */
    double bscale;		/* Scale factor for pixel scaling */
    char newname[128];
    char pixname[128];
    char history[64];
    FILE *fd;
    char *ext, *fname;
    char *editcom;
    char newline[1];
    char echar;
    double dpix;
    char *c;
    int bitpix,xdim,ydim;
    struct Range *xrange;    /* X range structure */
    struct Range *yrange;    /* Y range structure */

    newline[0] = 10;

    /* Open IRAF image and header */
    if (isiraf (filename)) {
	iraffile = 1;
	if ((irafheader = irafrhead (filename, &lhead)) != NULL) {
	    if ((header = iraf2fits (filename, irafheader, lhead, &nbhead)) == NULL) {
		free (irafheader);
                fprintf (stderr, "Cannot translate IRAF header %s/n", filename);
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
	    fprintf (stderr, "Cannot read IRAF header file %s\n", filename);
	    free (header);
	    return;
	    }
	}

    /* Read FITS image and header */
    else {
	iraffile = 0;
	if ((header = fitsrhead (filename, &lhead, &nbhead)) != NULL) {
	    if ((image = fitsrimage (filename, nbhead, header)) == NULL) {
		fprintf (stderr, "Cannot read FITS image %s\n", filename);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    return;
	    }
	}
    if (verbose)

    /* Change value of specified pixel */
    hgeti4 (header,"BITPIX",&bitpix);
    xdim = 1;
    hgeti4 (header,"NAXIS1",&xdim);
    ydim = 1;
    hgeti4 (header,"NAXIS2",&ydim);
    bzero = 0.0;
    hgetr8 (header,"BZERO",&bzero);
    bscale = 1.0;
    hgetr8 (header,"BSCALE",&bscale);

    i = 0;
    while (value[i] != NULL) {

	/* Extract new value from command line string */
	if (strchr (value[i], (int)'.'))
	    dpix = (double) atoi (value[i]);
	else
	    dpix = atof (value[i]);

	/* Set format if not already set */
	if (pform == NULL) {
	    pform = (char *) calloc (8,1);
	    if (bitpix > 0)
		strcpy (pform, "%d");
	    else
		strcpy (pform, "%.2f");
	    }

	/* Set entire columns */
	if (!strcmp (rrange[i], "0")) {
	    xrange = RangeInit (crange[i], xdim);
	    nx = rgetn (xrange);
	    ny = ydim;
	    for (ix = 0; ix < nx; ix++) {
		rstart (xrange);
		x = rgeti4 (xrange) - 1;
		for (y = 0; y < ny; y++) {
		    putpix (image,bitpix,xdim,ydim,bzero,bscale,x,y,dpix);
	            if (bitpix > 0) {
			if (dpix > 0)
	 		    ipix = (int) (dpix + 0.5);
			else if (dpix < 0)
		 	    ipix = (int) (dpix - 0.5);
			else
			    ipix = 0;
			}
		    if (eachpix) {
			printf ("%s[%d,%d] = ", filename,x+1,y+1);
		        if (bitpix > 0)
			    printf (pform, ipix);
			else
			    printf (pform, dpix);
			printf ("\n");
			}
		    }
		}
	    if (isnum (crange[i]))
		sprintf (history, "SETPIX: pixels in column %s set to %s",
		     crange[i],value[i]);
	    else
		sprintf (history, "SETPIX: pixels in columns %s set to %s",
		     crange[i],value[i]);
	    free (xrange);
	    }

	/* Set entire rows */
	else if (!strcmp (crange[i], "0")) {
	    yrange = RangeInit (rrange[i], xdim);
	    ny = rgetn (yrange);
	    nx = xdim;
	    for (iy = 0; iy < ny; iy++) {
		y = rgeti4 (yrange) - 1;
		for (x = 0; x < nx; x++) {
		    putpix (image,bitpix,xdim,ydim,bzero,bscale,x,y,dpix);
		    if (eachpix) {
			if (bitpix > 0) {
			    if (dpix > 0)
	 			ipix = (int) (dpix + 0.5);
			    else if (dpix < 0)
		 		ipix = (int) (dpix - 0.5);
			    else
				ipix = 0;
			    }
			printf ("%s[%d,%d] = ", filename,x+1,y+1);
			if (bitpix > 0)
			    printf (pform, ipix);
			else
			    printf (pform, dpix);
			printf ("\n");
			}
		    }
		}
	    if (isnum (rrange[i]))
		sprintf (history, "SETPIX: pixels in row %s set to %s",
		     rrange[i],value[i]);
	    else
		sprintf (history, "SETPIX: pixels in rows %s set to %s",
		     rrange[i],value[i]);
	    free (yrange);
	    }

	/* Set a region of a two-dimensional image */
	else {
	    xrange = RangeInit (crange[i], xdim);
	    nx = rgetn (xrange);

	    /* Make list of y coordinates */
	    yrange = RangeInit (rrange[i], ydim);
	    ny = rgetn (yrange);

	    /* Loop through rows starting with the last one */
	    for (i = 0; i < ny; i++) {
		if (verbose)
		    iy++;
		else
		    iy--;
		rstart (xrange);

		/* Loop through columns */
		for (ix = 0; ix < nx; ix++) {
		    x = rgeti4 (xrange) - 1;
		    putpix (image,bitpix,xdim,ydim,bzero,bscale,x,y,dpix);
        	
		    if (eachpix) {
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
			printf ("%s[%d,%d] = ", filename,x+1,y+1);
			if (bitpix > 0)
			    printf (pform, ipix);
			else
			    printf (pform, dpix);
			printf ("\n");
			}
		    }
		}

	    /* Note addition as history line in header */
	    if (isnum (crange[i]) && isnum (rrange[i]))
		sprintf (history, "SETPIX: pixel at row %s, column %s set to %s",
		     rrange[i], crange[i], value[i]);
	    else if (isnum (rrange[i]))
		sprintf (history, "SETPIX: pixels in row %s, columns %s set to %s",
		     rrange[i], crange[i], value[i]);
	    else if (isnum (crange[i]))
		sprintf (history, "SETPIX: pixels in column %s, rows %s set to %s",
		     crange[i], rrange[i], value[i]);
	    else
		sprintf (history, "SETPIX: pixels in rows %s, columns %s set to %s",
		     rrange[i], crange[i], value[i]);
	    }
	hputc (header,"HISTORY",history);
	if (verbose)
	    printf ("%s\n", history);
	free (xrange);
	free (yrange);
	i++;
	}

    /* Make up name for new FITS or IRAF output file */
    if (newimage) {

    /* Remove directory path and extension from file name */
	fname = strrchr (filename, '/');
	if (fname)
	    fname = fname + 1;
	else
	    fname = filename;
	ext = strrchr (fname, '.');
	if (ext != NULL) {
	    lext = (fname + strlen (fname)) - ext;
	    lroot = ext - fname;
	    strncpy (newname, fname, lroot);
	    *(newname + lroot) = 0;
	    }
	else {
	    lext = 0;
	    lroot = strlen (fname);
	    strcpy (newname, fname);
	    }
	imext = strchr (fname, ',');
	imext1 = NULL;
	if (imext == NULL) {
	    imext = strchr (fname, '[');
	    if (imext != NULL) {
		imext1 = strchr (fname, ']');
		*imext1 = (char) 0;
		}
	    }
	if (imext != NULL) {
	    strcat (newname, "_");
	    strcat (newname, imext+1);
	    }
	if (fname)
	    fname = fname + 1;
	else
	    fname = filename;
	strcat (newname, "e");
	if (lext > 0) {
	    if (imext != NULL) {
		echar = *imext;
		*imext = (char) 0;
		strcat (newname, ext);
		*imext = echar;
		if (imext1 != NULL)
		    *imext1 = ']';
		}
	    else
		strcat (newname, ext);
	    }
	}
    else
	strcpy (newname, filename);

    /* Write fixed header to output file */
    if (iraffile) {
	if (irafwimage (newname,lhead,irafheader,header,image) > 0 && verbose)
	    printf ("%s rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
	free (irafheader);
	}
    else {
	if (fitswimage (newname, header, image) > 0 && verbose)
	    printf ("%s: rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
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
 * May 28 1998	Include fitsio.h instead of fitshead.h
 * Jul 24 1998	Make irafheader char instead of int
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Aug 14 1998	Preserve extension when creating new file name
 * Oct 14 1998	Use isiraf() to determine file type
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Feb 12 1999	Initialize dxisn to 1 so it works for 1-D images
 * Apr 29 1999	Add BZERO and BSCALE
 * Jun 29 1999	Fix typo in BSCALE setting
 * Jul 12 1999	Add ranges
 * Sep 14 1999	Add file of values to usage
 */
