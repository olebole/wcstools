/* File conpix.c
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

#define PIX_ADD	1
#define PIX_SUB	2
#define PIX_MUL	3
#define PIX_DIV	4

static void usage();
static void OpPix();

static int newimage = 0;
static int verbose = 0;		/* verbose flag */
static int version = 0;		/* If 1, print only program name and version */


main (ac, av)
int ac;
char **av;
{
    char *str;
    char *fn;
    int readlist = 0;
    char *lastchar;
    char filename[128];
    FILE *flist;
    char *listfile;
    int nop, op[10];
    double opcon[10];

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    nop = 0;

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

	case 'a':	/* add constant to all pixels */
	    op[nop] = PIX_ADD;
            if (ac < 2)
                usage ();
            opcon[nop++] = atof (*++av);
	    ac--;
	    break;

	case 's':	/* subtract constant from all pixels */
	    op[nop] = PIX_SUB;
            if (ac < 2)
                usage ();
            opcon[nop++] = atof (*++av);
	    ac--;
	    break;

	case 'm':	/* multiply all pixels by constant */
	    op[nop] = PIX_MUL;
            if (ac < 2)
                usage ();
            opcon[nop++] = atof (*++av);
	    ac--;
	    break;

	case 'd':	/* divide all pixels by constant */
	    op[nop] = PIX_DIV;
            if (ac < 2)
                usage ();
            opcon[nop++] = atof (*++av);
	    ac--;
	    break;

	case '@':       /* List of files to be read */
	    readlist++;
	    listfile = ++str;
	    str = str + strlen (str) - 1;
	    av++;
	    ac--;

	default:
	    usage ();
	    break;
	}
    }

    /* Process files in file of filenames */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"CONPIX: List file %s cannot be read\n",
		     listfile);
	    usage ();
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
	    OpPix (filename, nop, op, opcon);
	    if (verbose)
		printf ("\n");
	    }
	fclose (flist);
	}

    /* If no arguments left, print usage */
    if (ac == 0)
	usage ();

    /* Process files on command line */
    else {
	while (ac-- > 0) {
    	    char *fn = *av++;
    	    if (verbose)
    		printf ("%s:\n", fn);
	    OpPix (fn, nop, op, opcon);
    	    if (verbose)
    		printf ("\n");
	    }
	}

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Operate on all pixels of a FITS or IRAF image file\n");
    fprintf(stderr,"Usage: conpix [-vn]{-asmd constant] file.fits ...\n");
    fprintf(stderr,"       conpix [-vn]{-asmd constant] @filelist\n");
    fprintf(stderr,"  -a: add constant to all pixels\n");
    fprintf(stderr,"  -d: divide all pixels by constant\n");
    fprintf(stderr,"  -m: multiply all pixels by constant\n");
    fprintf(stderr,"  -n: write new file, else overwrite\n");
    fprintf(stderr,"  -s: subtract constant from all pixels\n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
OpPix (filename, nop, op, opcon)

char	*filename;	/* FITS or IRAF file filename */
int	nop;		/* Number of pixels to change */
int	*op;		/* List of operations to perform */
double	*opcon;		/* Constants for operations */

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
    char echar;
    double *imvec, *dvec, *endvec;
    int bitpix, xdim, ydim, x, y, pixoff, iop;
    double bzero;		/* Zero point for pixel scaling */
    double bscale;		/* Scale factor for pixel scaling */

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

    /* Add specified value to specified pixel */
    hgeti4 (header,"BITPIX",&bitpix);
    xdim = 1;
    hgeti4 (header,"NAXIS1",&xdim);
    ydim = 1;
    hgeti4 (header,"NAXIS2",&ydim);
    bzero = 0.0;
    hgetr8 (header,"BZERO",&bzero);
    bscale = 1.0;
    hgetr8 (header,"BZERO",&bscale);

    if (!(imvec = (double *) calloc (xdim, sizeof (double))))
	return;
    endvec = imvec + xdim;

    pixoff = 0;
    for (y = 0; y < ydim; y++) {
	getvec (image, bitpix, bzero, bscale, pixoff, xdim, imvec);
	for (iop = 0; iop < nop; iop++) {
	    double dpix = opcon[iop];
	    switch (op[iop]) {
		case PIX_ADD:
		    for (dvec = imvec; dvec < endvec; dvec++)
			*dvec = *dvec + dpix;
		    break;
		case PIX_SUB:
		    for (dvec = imvec; dvec < endvec; dvec++)
			*dvec = *dvec - dpix;
		    break;
		case PIX_MUL:
		    for (dvec = imvec; dvec < endvec; dvec++)
			*dvec = *dvec * dpix;
		    break;
		case PIX_DIV:
		    for (dvec = imvec; dvec < endvec; dvec++)
			*dvec = *dvec / dpix;
		    break;
		default:
		    break;
		}
	    }
	putvec (image, bitpix, pixoff, xdim, imvec);
	pixoff = pixoff + xdim;
	if (verbose) {
	    fprintf (stderr, "Row %4d operations complete", y);
	    (void) putc (13,stderr);
	    }
	}
    fprintf (stderr,"\n");

    /* Note addition as history line in header */
    for (iop = 0; iop < nop; iop++) {
	double dpix = opcon[iop];
	if (bitpix > 0) {
	    int ipix = (int)dpix;
	    switch (op[iop]) {
		case PIX_ADD:
		    sprintf (history, "CONPIX: %d added to all pixels", ipix);
		    break;
		case PIX_SUB:
		    sprintf (history,
			    "CONPIX: %d subtracted from all pixels", ipix);
		    break;
		case PIX_MUL:
		    sprintf (history,
			     "CONPIX: all pixels multiplied by %d", ipix);
		    break;
		case PIX_DIV:
		    sprintf (history,
			     "CONPIX: all pixels divided by %d", ipix);
		    break;
		}
	    }
	else if (dpix < 1.0 && dpix > -1.0) {
	    switch (op[iop]) {
		case PIX_ADD:
		    sprintf (history, "CONPIX: %f added to all pixels", dpix);
		    break;
		case PIX_SUB:
		    sprintf (history,
			    "CONPIX: %f subtracted from all pixels", dpix);
		    break;
		case PIX_MUL:
		    sprintf (history,
			     "CONPIX: all pixels multiplied by %f", dpix);
		    break;
		case PIX_DIV:
		    sprintf (history,
			     "CONPIX: all pixels divided by %f", dpix);
		    break;
		}
	    }
	else {
	    switch (op[iop]) {
		case PIX_ADD:
		    sprintf (history, "CONPIX: %.2f added to all pixels", dpix);
		    break;
		case PIX_SUB:
		    sprintf (history,
			    "CONPIX: %.2f subtracted from all pixels", dpix);
		    break;
		case PIX_MUL:
		    sprintf (history,
			     "CONPIX: all pixels multiplied by %.2f", dpix);
		    break;
		case PIX_DIV:
		    sprintf (history,
			     "CONPIX: all pixels divided by %.2f", dpix);
		    break;
		}
	    }
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

/* Dec  2 1998	New program
 *
 * Feb 12 1999	Initialize dimensions to one so it works with 1-D images
 * Apr 29 1999	Add BZERO and BSCALE
 */
