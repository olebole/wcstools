/* File fixpix.c
 * July 12, 1997
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include "libwcs/fitshead.h"
#include "libwcs/wcs.h"

#define MAXFIX 10
#define MAXFILES 50

static void FixPix();
static void FixReg();
static void usage();
static int newimage = 0;
static int verbose = 0;		/* verbose flag */
static int nfix = 0;		/* Number of regions to fix
				   If < 0 read regions from a file */
static int x1[MAXFIX],y1[MAXFIX]; /* Lower left corners of regions (1 based) */
static int x2[MAXFIX],y2[MAXFIX]; /* Upper right corners of regions (1 based) */

static void SetPix();

main (ac, av)
int ac;
char **av;
{
    char *str;
    char *fn[MAXFILES];
    int readlist = 0;
    int lfn, nfile, maxlfn;
    char *lastchar;
    char filename[128];
    int ifile;
    FILE *flist;
    char *listfile;
    char *regionlist;

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

	case '@':	/* List of files to be read */
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

    /* If there are no remaining arguments, print usage */
    if (ac == 0)
	usage ();

    /* Crack arguments */
    nfix = 0;
    nfile = 0;
    while (ac-- > 0  && nfile < MAXFILES && nfix < MAXFIX) {
	if (strsrch (*av,".fit") != NULL ||
	    strsrch (*av,".fts") != NULL ||
	    strsrch (*av,".FIT") != NULL ||
	    strsrch (*av,".FTS") != NULL ||
	    strsrch (*av,".imh") != NULL) {
	    fn[nfile] = *av;
	    lfn = strlen (fn[nfile]);
	    if (lfn > maxlfn)
		maxlfn = lfn;
	    nfile++;
	    av++;
	    }
	else if (*av[0] == '@') {
	    nfix = -1;
	    av++;
	    }
	else if (ac > 2) {
	    x1[nfix] = atoi (*av++);
	    ac--;
	    y1[nfix] = atoi (*av++);
	    ac--;
	    x2[nfix] = atoi (*av++);
	    ac--;
	    y2[nfix] = atoi (*av++);
	    nfix++;
	    }
	}

    /* Process only if a list of regions to fix has been found */
    if (nfix != 0) {

	/* Read through headers of images in listfile */
	if (readlist) {
	    if ((flist = fopen (listfile, "r")) == NULL) {
		fprintf (stderr,"FIXPIX: List file %s cannot be read\n",
		     listfile);
		usage ();
		}
	    while (fgets (filename, 128, flist) != NULL) {
		lastchar = filename + strlen (filename) - 1;
		if (*lastchar < 32) *lastchar = 0;
		FixPix (filename, regionlist);
		if (verbose)
		    printf ("\n");
		}
	    fclose (flist);
	    }

	/* Read image headers from command line list */
	else {
	    for (ifile = 0; ifile < nfile; ifile++)
		FixPix (fn[ifile], regionlist);
	    }
	}

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"Fix pixel regions of FITS or IRAF image file\n");
    fprintf(stderr,"Usage: fixpix [-vn] file.fits x1 y1 x2 y2...\n");
    fprintf(stderr,"Usage: fixpix [-vn] file.fits @regionlist\n");
    fprintf(stderr,"Usage: fixpix [-vn] @filelist x1 y1 x2 y2...\n");
    fprintf(stderr,"Usage: fixpix [-vn] @filelist @regionlist\n");
    fprintf(stderr,"  -n: write new file, else overwrite \n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
FixPix (filename, regionlist)

char	*filename;	/* FITS or IRAF file filename */
char	*regionlist;	/* Name of file of regions to fix, if nfix < 0 */

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
    char line[128];
    char history[64];
    FILE *fd, *freg;
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

    /* Fix values of specified area */
    hgeti4 (header,"BITPIX",&bitpix);
    hgeti4 (header,"NAXIS1",&xdim);
    hgeti4 (header,"NAXIS2",&ydim);

    /* Fix pixels over regions from a command line coordinate list */
    if (nfix > 0) {
	for (i = 0; i < nfix; i++) {
	    FixReg (image, bitpix, xdim, ydim, x1[i], y1[i], x2[i], y2[i]);

	    /* Note addition as history line in header */
	    sprintf (history, "FIXPIX: region x: %d-%d, y: %d-%d replaced",
		     x1[i],x2[i],y1[i],y2[i]);
	    hputc (header,"HISTORY",history);
	    if (verbose)
		printf ("%s\n", history);
	    }
	}

    /* Fix pixels over regions from a file */
    else {
	if ((freg = fopen (regionlist, "r")) == NULL) {
		fprintf (stderr,"FIXPIX: Region file %s cannot be read\n",
		     regionlist);
		usage ();
		}
	while (fgets (line, 128, freg) != NULL) {
	    sscanf (line,"%d %d %d %d", x1[1], y1[1], x1[1], y2[1]);
	    FixReg (image, bitpix, xdim, ydim, x1[i], y1[i], x2[i], y2[i]);

	    /* Note addition as history line in header */
	    sprintf (history, "FIXPIX: region x: %d-%d, y: %d-%d replaced",
			 x1[i],x2[i],y1[i],y2[i]);
	    hputc (header,"HISTORY",history);
	    if (verbose)
		printf ("%s\n", history);
	    }
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

    /* Write fixed image to output file */
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

static void
FixReg (image, bitpix, xdim, ydim, ix1, iy1, ix2, iy2)

char *image;		/* FITS image */
int bitpix;	/* Number of bits in each pixel */
int xdim;	/* Number of pixels in image horizontally */
int ydim;	/* Number of pixels in image vertically */
int ix1, iy1;	/* Lower left corner of region (1 based) */
int ix2, iy2;	/* Upper right corner of region (1 based) */

{
    int xdiff, ydiff, it, ix, iy;
    double pix1, pix2, dpix;

    /* Find dimensions of region to fix */
    if (ix1 > ix2) {
	it = ix2;
	ix2 = ix1;
	ix1 = it;
	}
    xdiff = ix2 - ix1 + 1;
    if (iy1 > iy2) {
	it = iy2;
	iy2 = iy1;
	iy1 = it;
	}
    ydiff = iy2 - iy1 + 1;

    /* Return if region contains no points */
    if (xdiff < 1 || ydiff < 1)
	return;

    /* If more horizontal than vertical, interpolate vertically */
    if (xdiff > ydiff) {
	if (iy1 - 1 < 0 || iy2 + 1 > ydim - 1)
	    return;
	for (ix = ix1; ix <= ix2; ix++) {
	    pix1 = getpix (image, bitpix, xdim, ydim, ix, iy1-1);
	    pix2 = getpix (image, bitpix, xdim, ydim, ix, iy2+1);
	    dpix = (pix2 - pix1) / (double)(ydiff + 1);
	    for (iy = iy1; iy <= iy2; iy++) {
		pix1 = pix1 + dpix;
		putpix (image, bitpix, xdim, ydim, ix, iy, pix1);
		}
	    }
	}

    /* If more vertical than horizontal, interpolate horizontally */
    else {
	if (ix1 - 1 < 0 || ix2 + 1 > xdim - 1)
	    return;
	for (iy = iy1; iy <= iy2; iy++) {
	    pix1 = getpix (image, bitpix, xdim, ydim, ix1-1, iy);
	    pix2 = getpix (image, bitpix, xdim, ydim, ix2+1, iy);
	    dpix = (pix2 - pix1) / (double)(ydiff + 1);
	    for (ix = ix1; ix <= ix2; ix++) {
		pix1 = pix1 + dpix;
		putpix (image, bitpix, xdim, ydim, ix, iy, pix1);
		}
	    }
	}

    return;
}

/* Jul 12 1997	New program
 */
