/* File wcshead.c
 * February 18, 1998
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
#include "libwcs/fitshead.h"
#include "libwcs/wcs.h"

static void usage();
static void PrintHead();

static int nskip = 0;		/* Number of bytes to skip */
static int nfiles = 0;		/* Nuber of files for headers */
static int verbose = 0;		/* verbose/debugging flag */
static int tabout = 0;		/* tab table output flag */
static int ndec = 3;		/* maximum number of decimal places for output*/
static int nchar = 16;		/* maximum number of characters for filename */
static int hms = 0;		/* 1 for output in hh:mm:ss dd:mm:ss */
static int nf = 0;

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    int readlist = 0;
    char *lastchar;
    char filename[128];
    FILE *flist;
    char *listfile;

    /* crack arguments */
    for (av++; --ac > 0 && (*(str = *av)=='-' || *str == '@'); av++) {
	char c;
	if (*str == '@')
	    str = str - 1;
	while (c = *++str)
	switch (c) {

	case 'h':	/* hh:mm:ss output for crval, cdelt in arcsec/pix */
	    hms++;
	    break;

	case 'n':	/* hh:mm:ss output */
	    tabout++;
	    break;

	case 't':	/* tab table output */
	    tabout++;
	    break;

	case 'v':	/* more verbosity */
	    verbose++;
	    break;

	case '@':	/* List of files to be read */
	    readlist++;
	    listfile = ++str;
	    str = str + strlen (str) - 1;
	    av++;
	    ac--;

	default:
	    usage(progname);
	    break;
	}
    }

    /* Find number of images to search and leave listfile open for reading */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMHEAD: List file %s cannot be read\n",
		     listfile);
	    usage (progname);
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
	    PrintHead (filename);
	    if (verbose)
		printf ("\n");
	    }
	fclose (flist);
	}

    /* If no arguments left, print usage */
    if (ac == 0)
	usage (progname);

    if (verbose)

    nfiles = ac;
    nf = 0;
    while (ac-- > 0) {
	char *fn = *av++;
	nf++;
	PrintHead (fn);
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"Print WCS part of FITS or IRAF image header\n");
    fprintf (stderr,"%s: usage: [-htv] file.fit ...\n", progname);
    fprintf (stderr,"  -h: print CRVALs as hh:mm:ss dd:mm:ss\n");
    fprintf (stderr,"  -t: print tab table output\n");
    fprintf (stderr,"  -v: verbose\n");
    exit (1);
}


static void
PrintHead (filename)

char	*filename;	/* FITS or IRAF image file name */

{
    struct WorldCoor *wcs, *GetWCSFITS();
    void TabWCS();

    wcs = GetWCSFITS (filename);
    if (nowcs (wcs))
	return;

    TabWCS (filename, wcs);

    free (wcs);
    return;
}


void
TabWCS (filename, wcs)

char	*filename;	/* FITS or IRAF image file name */

struct WorldCoor *wcs;

{
    int i;
    char rastr[32], decstr[32], fform[8];

    if (wcs->ctype[0][0] == (char) 0)
	return;
    if (tabout && nf == 1) {
	printf ("filename");
	for (i = 1; i < nchar - 8; i++) printf (" ");
	printf ("	naxis1	naxis2");
	printf ("	ctype1  	ctype2  ");
	printf ("	crval1	crval2	radecsys");
	printf ("	crpix1	crpix2");
	printf ("	cdelt1	cdelt2");
	printf ("	crota2\n");
	printf ("--------");
	for (i = 1; i < nchar - 8; i++) printf ("-");
	printf ("	------	------");
	printf ("	------	------	--------");
	printf ("	-------	-------");
	printf ("	------	------");
	printf ("	------	------");
	printf ("	--------\n");
	}

    sprintf (fform,"%%%d.%ds",nchar);
    if (tabout)
	printf (fform, filename);
    else
	printf (fform, filename);
    if (tabout)
	printf ("	%.0f	%.0f", wcs->nxpix, wcs->nypix);
    else
	printf (" %4.0f %4.0f", wcs->nxpix, wcs->nypix);

    if (tabout)
	printf ("	%s	%s", wcs->ctype[0], wcs->ctype[1]);
    else
	printf (" %s %s", wcs->ctype[0], wcs->ctype[1]);

    if (tabout) {
	if (hms) {
	    if (wcs->coorflip) {
		ra2str (rastr, wcs->yref, ndec);
		dec2str (decstr, wcs->xref, ndec-1);
		}
	    else {
		ra2str (rastr, wcs->xref, ndec);
		dec2str (decstr, wcs->yref, ndec-1);
		}
	    printf (" %s %s %s", rastr, decstr, wcs->radecsys);
	    }
	else
	    printf ("	%7.2f	%7.2f	%s", wcs->xref, wcs->yref, wcs->radecsys);
	}
    else {
	if (hms) {
	    if (wcs->coorflip) {
		ra2str (rastr, wcs->yref, ndec);
		dec2str (decstr, wcs->xref, ndec-1);
		}
	    else {
		ra2str (rastr, wcs->xref, ndec);
		dec2str (decstr, wcs->yref, ndec-1);
		}
	    printf (" %s %s %s", rastr, decstr, wcs->radecsys);
	    }
	else
	    printf (" %7.2f %7.2f %s", wcs->xref, wcs->yref, wcs->radecsys);
	}

    if (tabout)
	printf ("	%7.2f	%7.2f", wcs->xrefpix, wcs->yrefpix);
    else
	printf (" %7.2f %7.2f", wcs->xrefpix, wcs->yrefpix);

    if (tabout)
	printf ("	%7.4f	%7.4f", 3600.0*wcs->xinc, 3600.0*wcs->yinc);
    else
	printf (" %7.4f %7.4f", 3600.0*wcs->xinc, 3600.0*wcs->yinc);

    if (tabout)
	printf ("	%7.4f\n", wcs->rot);
    else
	printf (" %7.4f\n", wcs->rot);

    return;
}
/* Feb 18 1998	New program
 */
