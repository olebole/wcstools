/* File subpix.c
 * August 14, 1998
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

static void usage();
static int newimage = 0;
static int verbose = 0;		/* verbose flag */
static void SubPix();

main (ac, av)
int ac;
char **av;
{
    char *str;
    char *fn;
    char *value[100];
    int i, x[100], y[100];

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {
	case 'v':	/* more verbosity */
	    verbose++;
	    break;
	case 'n':	/* ouput new file */
	    newimage++;
	    break;
	default:
	    usage ();
	    break;
	}
    }

    /* now there are ac remaining arguments starting at av[0] */
    if (ac == 0)
	usage ();

    fn = *av++;

    i = 0;
    while (--ac > 2) {
	x[i] = atoi (*av++);
	ac--;
	y[i] = atoi (*av++);
	ac--;
	value[i] = *av++;
	i++;
	}
    if (i > 0)
        SubPix (fn,i,x,y,value);

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"Subtract from pixel of FITS or IRAF image file\n");
    fprintf(stderr,"Usage: subpix [-vn] file.fts x y value...\n");
    fprintf(stderr,"  -n: write new file, else overwrite \n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
SubPix (filename, n, x, y, value)

char	*filename;	/* FITS or IRAF file filename */
int	n;		/* number of pixels to change */
int	*x, *y;		/* Horizontal and vertical coordinates of pixel */
			/* (1-based) */
char	**value;	/* value to insert into pixel */

{
    char *image;		/* FITS image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    int iraffile;		/* 1 if IRAF image */
    char *irafheader;		/* IRAF image header */
    int i, nbytes, nhb, nhblk, lname, lext, lroot;
    char *head, *headend, *hlast, *imext, *imext1;
    char headline[160];
    char newname[128];
    char pixname[128];
    char tempname[128];
    char history[64];
    FILE *fd;
    char *ext, *fname;
    char *editcom;
    char echar;
    char newline[1];
    double dpix, dpix0, dpix1;
    int bitpix,xdim,ydim;

    newline[0] = 10;
    strcpy (tempname, "fitshead.temp");

    /* Open IRAF image and header if .imh extension is present */
    if (strsrch (filename,".imh") != NULL) {
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

    /* Read FITS image and header if .imh extension is not present */
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

    /* Add specified value to specified pixel */
    hgeti4 (header,"BITPIX",&bitpix);
    hgeti4 (header,"NAXIS1",&xdim);
    hgeti4 (header,"NAXIS2",&ydim);

    for (i = 0; i < n; i++) {
	if (strchr (value[i],(int)'.'))
	    dpix = (double) atoi (value[i]);
	else
	    dpix = atof (value[i]);
	dpix0 = getpix (image, bitpix, xdim, ydim, x[i]-1, y[i]-1);
	dpix1 = dpix0 - dpix;
	putpix (image, bitpix, xdim, ydim, x[i]-1, y[i]-1, dpix1);

	/* Note addition as history line in header */
	if (bitpix > 0) {
	    int ipix = (int)dpix;
	    sprintf (history, "SUBPIX: %d subtracted from pixel at row %d, column %d",
		     ipix,x[i],y[i]);
	    }
	else if (dpix < 1.0 && dpix > -1.0)
	    sprintf (history, "SUBPIX: %f subtracted from pixel at row %d, column %d",
		     dpix,x[i],y[i]);
	else
	    sprintf (history,"SUBPIX: %.2f subtracted from pixel at row %d, column %d",
		     dpix,x[i],y[i]);
	hputc (header,"HISTORY",history);
	if (verbose)
	    printf ("%s\n", history);
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

    free (image);
    free (header);
    return;
}

/* Dec  6 1996	New program
 *
 * Jan 15 1997	Print subtracted value rather than result in verbose mode
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 * Dec 15 1997	Add capability of reading and writing IRAF 2.11 images
 *
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jul 24 1998	Make irafheader char instead of int
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Aug 14 1998	Preserve extension when creating new file name
 */
