/* File getwcs.c
 * February 23, 1996
 * By Doug Mink
 */

#include <stdlib.h>
#include <stdio.h>
#include "fitshead.h"
#include "wcs.h"

#define MAXHEADLEN 14400

struct WorldCoor *wcsinit();	

struct WorldCoor *
getWCS (filename)

char *filename;	/* FITS or IRAF file filename */

{
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int iraffile;		/* 1 if IRAF image */
    int *irafheader;		/* IRAF image header */
    struct WorldCoor *wcs;	/* World coordinate system structure */

    /* Allocate FITS header */
    header = malloc (MAXHEADLEN);
    lhead = MAXHEADLEN;

    /* Open IRAF image if .imh extension is present */
    if (strsrch (filename,".imh") != NULL) {
	iraffile = 1;
	irafheader = irafrhead (filename, lhead, header);
	if (irafheader == NULL) {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", filename);
	    free (header);
	    return (NULL);
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	if (fitsrhead (filename, lhead, header) < 1) {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    free (header);
	    return (NULL);
	    }
	}

    /* Set the world coordinate system from the image header */
    wcs = wcsinit (header);

    free (header);
    if (iraffile)
	free (irafheader);
    return (wcs);
}
