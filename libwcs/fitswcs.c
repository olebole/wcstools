/* File libwcs/fitswcs.c
 * February 21, 1997
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics

 * Module:      fitswcs.c (FITS file WCS reading and deleting)
 * Purpose:     Read and delete FITS image world coordinate system keywords
 * Subroutine:  GetWCSFITS (filename)
 *		Open a FITS or IRAF image file and returns its WCS structure
 * Subroutine:  GetFITShead (filename)
 *		Open a FITS or IRAF image file and returns a FITS header
 * Subroutine:  DelWCSFITS (filename, verbose)
 *		Delete all standard WCS keywords from a FITS header
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fitshead.h"
#include "wcs.h"


struct WorldCoor *
GetWCSFITS (filename)

char *filename;	/* FITS or IRAF file filename */

{
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int *irafheader;		/* IRAF image header */
    struct WorldCoor *wcs;	/* World coordinate system structure */
    int nbiraf, nbfits;
    char *GetFITShead();

    /* Read the FITS or IRAF image file header */
    header = GetFITShead (filename);

    /* Set the world coordinate system from the image header */
    wcs = wcsinit (header);
    free (header);

    return (wcs);
}

char *
GetFITShead (filename)

char *filename;	/* FITS or IRAF file filename */

{
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int *irafheader;		/* IRAF image header */
    int nbiraf, nbfits;

    /* Open IRAF image if .imh extension is present */
    if (strsrch (filename,".imh") != NULL) {
	if ((irafheader = irafrhead (filename, &nbiraf)) != NULL) {
	    if ((header = iraf2fits (filename, irafheader, nbiraf, &lhead)) == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",filename);
		free (irafheader);
		return;
		}
	    free (irafheader);
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", filename);
	    return (NULL);
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	if ((header = fitsrhead (filename, &lhead, &nbfits)) == NULL) {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    return (NULL);
	    }
	}

    return (header);
}


/* delete all the C* fields.
 * return 0 if at least one such field is found, else -1.  */

int
DelWCSFITS (header, verbose)

char *header;
int verbose;

{
    static char *flds[] = {
	"CTYPE1", "CRVAL1", "CDELT1", "CRPIX1", "CROTA1",
	"CTYPE2", "CRVAL2", "CDELT2", "CRPIX2", "CROTA2", "IMWCS" };
    int i;
    int n;
    double eq;
    char rastr[16],decstr[16];

    n = 0;

    for (i = 0; i < sizeof(flds)/sizeof(flds[0]); i++) {
	if (hdel (header, flds[i])) {
	    n++;
	    if (verbose)
		printf ("%s: deleted\n", flds[i]);
	    }
	}
    if (verbose && n == 0)
	printf ("DelWCSFITS: No WCS in header\n");

    if (ksearch (header,"WRA")) {
	hdel (header, "RA");
	hchange (header, "WRA","RA");
	if (ksearch (header,"WDEC")) {
	    hdel (header, "DEC");
	    hchange (header, "WDEC", "DEC");
	    }
	if (ksearch (header,"WEPOCH")) {
	    hdel (header, "EPOCH");
	    hchange (header, "WEPOCH", "EPOCH");
	    }
	if (ksearch (header,"WEQUINOX")) {
	    hdel (header, "EQUINOX");
	    hchange (header, "WEQUINOX", "EQUINOX");
	    }
	if (verbose) {
	    hgets (header,"RA", 16, rastr);
	    hgets (header,"DEC", 16, decstr);
	    eq = 0.0;
	    hgetr8 (header,"EPOCH",&eq);
	    if (eq == 0.0)
		hgetr8 (header,"EQUINOX",&eq);
	    printf ("DelWCS: Center reset to %s %s %.1f\n", rastr,decstr, eq);
	    }
	}
    else if (ksearch (header, "EPOCH") && !ksearch (header, "PLTRAH")) {
	if (hdel (header,"EQUINOX")) {
	    if (verbose)
		printf ("EQUINOX: deleted\n");
	    n++;
	    }
	else if (verbose)
	    printf ("DelWCS: EPOCH, but not EQUINOX found\n");
	}

    return (n);
}
/* May 29 1996	Change name from delWCSFITS to DelWCSFITS
 * May 31 1996	Print single message if no WCS is found in header
 * May 31 1996	Use stream I/O instead of standard I/O
 * Jun 10 1996	Combine imgetwcs.c and imdelwcs.c into fitswcs.c
 * Jun 17 1996	Delete IMWCS record, too
 * Jul 16 1996	Update arguments for header-reading subroutines
 * Aug  6 1996  Fixed small defects after lint
 * Aug  8 1996  Restore old image center after deleting WCS
 * Aug 26 1996	Fix subroutine arguments after lint
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 * Feb 21 1997  Add GetFITShead subroutine and use it
 */
