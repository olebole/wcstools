/* File setpix.c
 * December 6, 1996
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
#include "libwcs/wcs.h"

static void usage();
static int newimage = 0;
static int verbose = 0;		/* verbose flag */
static void SetPix();

main (ac, av)
int ac;
char **av;
{
    char *str;
    char *fn;
    char *value[100];
    int i,x[100],y[100];

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
    int *irafheader;		/* IRAF image header */
    int i, nbytes, nhb, nhblk, lname, lext;
    char *head, *headend, *hlast;
    char headline[160];
    char newname[128];
    char pixname[128];
    char tempname[128];
    char history[64];
    FILE *fd;
    char *ext, *fname;
    char *editcom;
    char newline[1];
    double dpix;
    int bitpix,xdim,ydim;

    newline[0] = 10;
    strcpy (tempname, "fitshead.temp");

    /* Open IRAF image and header if .imh extension is present */
    if (strsrch (filename,".imh") != NULL) {
	iraffile = 1;
	if ((irafheader = irafrhead (filename, &lhead))) {
	    header = iraf2fits (filename, irafheader, lhead, &nbhead);
            if (header == NULL) {
		free (irafheader);
                fprintf (stderr, "Cannot translate IRAF header %s/n", filename);
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
	    fprintf (stderr, "Cannot read IRAF header file %s\n", filename);
	    free (header);
	    return;
	    }
	}

    /* Read FITS image and header if .imh extension is not present */
    else {
	iraffile = 0;
	if ((header = fitsrhead (filename, &lhead, &nbhead))) {
	    if (!(image = fitsrimage (filename, nbhead, header))) {
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
    hgeti4 (header,"NAXIS1",&xdim);
    hgeti4 (header,"NAXIS2",&ydim);
    for (i = 0; i < n; i++) {
	if (strchr (value[i],(int)'.'))
	    dpix = (double) atoi (value[i]);
	else
	    dpix = atof (value[i]);
	putpix (image, bitpix, xdim, ydim, x[i]-1, y[i]-1, dpix);

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
	ext = strrchr (filename, '.');
	fname = strrchr (filename, '/');
	if (fname)
	    fname = fname + 1;
	else
	    fname = filename;
	lname = strlen (fname);
	if (ext) {
	    lext = strlen (ext);
	    strncpy (newname, fname, lname - lext);
	    *(newname + lname - lext) = 0;
	    }
	else
	    strcpy (newname, fname);

    /* Add file extension preceded by a e */
	if (iraffile)
	    strcat (newname, "e.imh");
	else
	    strcat (newname, "e.fit");
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
 */
