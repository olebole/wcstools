/* File setpix.c
 * April 29, 1999
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
static void SetPix();

static int newimage = 0;
static int verbose = 0;		/* verbose flag */
static int version = 0;		/* If 1, print only program name and version */


main (ac, av)
int ac;
char **av;
{
    char *str;
    char *fn;
    char *value[100];
    int i,x[100],y[100];

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

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
    while (--ac > 2 && i < 100) {
	x[i] = atoi (*av++);
	ac--;
	y[i] = atoi (*av++);
	ac--;
	value[i] = *av++;
	i++;
	}
    if (i > 0)
	SetPix (fn,i,x,y,value);

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Edit pixels of FITS or IRAF image file\n");
    fprintf(stderr,"Usage: setpix [-vn] file.fts x y value ...\n");
    fprintf(stderr,"  -n: write new file, else overwrite \n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
SetPix (filename,n,x,y,value)

char	*filename;	/* FITS or IRAF file filename */
int	n;		/* Number of pixels to set */
int	*x,*y;		/* Horizontal and vertical coordinates of pixel */
			/* (1-based) */
char	**value;		/* value to insert into pixel */

{
    char *image;		/* FITS image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    int iraffile;		/* 1 if IRAF image */
    char *irafheader;		/* IRAF image header */
    int i, nbytes, nhb, nhblk, lname, lext, lroot;
    char *head, *headend, *hlast, *imext, *imext1;
    double bzero;		/* Zero point for pixel scaling */
    double bscale;		/* Scale factor for pixel scaling */
    char headline[160];
    char newname[128];
    char pixname[128];
    char tempname[128];
    char history[64];
    FILE *fd;
    char *ext, *fname;
    char *editcom;
    char newline[1];
    char echar;
    double dpix;
    int bitpix,xdim,ydim;

    newline[0] = 10;
    strcpy (tempname, "fitshead.temp");

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
    hgetr8 (header,"BZERO",&bscale);

    for (i = 0; i < n; i++) {
	if (strchr (value[i],(int)'.'))
	    dpix = (double) atoi (value[i]);
	else
	    dpix = atof (value[i]);
	putpix (image,bitpix,xdim,ydim,bzero,bscale,x[i]-1,y[i]-1,dpix);

	/* Note addition as history line in header */
	if (bitpix > 0) {
	    int ipix = (int)dpix;
	    sprintf (history, "SETPIX: pixel at row %d, column %d set to %d",
		     x[i],y[i],ipix);
	    }
	else if (dpix < 1.0 && dpix > -1.0)
	    sprintf (history, "SETPIX: pixel at row %d, column %d set to %f",
		     x[i],y[i],dpix);
	else
	    sprintf (history, "SETPIX: pixel at row %d, column %d set to %.2f",
		     x[i],y[i],dpix);
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
 */
