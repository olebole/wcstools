/* File i2f.c
 * November 30, 1998
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

static void usage();
static void IRAFtoFITS ();

static int verbose = 0;		/* verbose/debugging flag */
static int version = 0;		/* If 1, print only program name and version */

main (ac, av)
int ac;
char **av;
{
    char *str;

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
	while (c = *++str) {
	    switch (c) {
		case 'v':	/* more verbosity */
		    verbose++;
		    break;
	        default:
		    usage();
		    break;
		}
    	    }
	}

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage ();

    else {
	while (ac-- > 0) {
	    char *fn = *av++;
	    if (verbose)
		printf ("%s:\n", fn);
	    IRAFtoFITS (fn);
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
    fprintf (stderr,"Write FITS files from IRAF image files\n");
    fprintf(stderr,"usage: i2f [-v] file.imh ...\n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}

static void
IRAFtoFITS (name)
char *name;
{
    char *image;	/* FITS image */
    char *header;	/* FITS header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    char *irafheader;	/* IRAF image header */
    char pixname[128];	/* Pixel file name */
    char history[128];	/* for HISTORY line */
    char *filename;	/* Pointer to start of file name */
    char irafname[128];	/* Name of IRAF file */
    char fitsname[128];	/* Name of FITS file */
    char *ext;		/* Pointer to start of extension */
    char *endchar;
    char *ltime;

    /* Open IRAF image if .imh extension is present */
    if (strsrch (name,".imh") != NULL) {
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
	    fprintf (stderr,"Cannot read IRAF header file %s\n", name);
	    return;
	    }
	strcpy (fitsname, name);
	ext = strsrch (fitsname,".imh");
	strcpy (ext,".fits");
	if (verbose) {
	    fprintf (stderr,"Write FITS files from IRAF image file %s\n", name);
	    }
	}

    /* Add .imh extension to make IRAF header file name if not present */
    else {
	strcpy (irafname, name);
	strcat (irafname,".imh");
	if ((irafheader = irafrhead (irafname, &lhead)) != NULL) {
	    header = iraf2fits (irafname, irafheader, lhead, &nbhead);
	    free (irafheader);
	    if (header == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",irafname);
		return;
		}
	    if ((image = irafrimage (header)) == NULL) {
		hgets (header,"PIXFILE",64, pixname);
		fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
		free (irafheader);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr,"Cannot read IRAF header file %s\n", irafname);
	    return;
	    }
	strcpy (fitsname, name);
	strcat (fitsname,".fits");
	if (verbose) {
	    fprintf (stderr,"Write FITS files from IRAF image file %s\n", irafname);
	    }
	}

    /* Add HISTORY notice of this conversion */
    filename = strrchr (name,'/');
    if (filename)
	filename = filename + 1;
    else
	filename = name;
    endchar = strchr (history, ',');
    *endchar = (char) 0;
    strcat (history, " ");
    ltime = getltime ();
    strcat (history, ltime);
    endchar = strrchr (history,':');
    *endchar = (char) 0;
    strcat (history, " FITS from ");
    strcat (history, filename);
    if (strlen (history) > 72)
	history[72] = 0;
    hputc (header, "HISTORY", history);

    /* Write FITS image */
    if (fitswimage (fitsname, header, image) > 0 && verbose)
	printf ("%s: rewritten successfully.\n", fitsname);

    else if (verbose)
	printf ("%s: not written.\n", fitsname);

    free (header);
    free (image);
    return;
}
/* Jun  6 1996	New program
 * Jul 16 1996	Update header input
 * Aug 16 1996	Clean up code
 * Aug 26 1996	Change HGETC call to HGETS
 * Aug 27 1996	Drop unused variables after lint
 * Oct 17 1996	Clean up after lint
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 *
 * Jan 14 1998	Version 1.3 to handle IRAF 2.11 .imh files
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jun  2 1998  Fix bug in hput()
 * Jul 24 1998	Make irafheader char instead of int
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Aug 14 1998	Write file.fits instead of file.fit
 * Aug 17 1998	Add HISTORY to header
 * Nov 30 1998	Add version and help commands for consistency
 */
