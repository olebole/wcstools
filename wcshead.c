/* File wcshead.c
 * November 30, 1999
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
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"

static void usage();
static void ListWCS();

static int verbose = 0;		/* verbose/debugging flag */
static int tabout = 0;		/* tab table output flag */
static int ndec = 3;		/* maximum number of decimal places for output*/
static int nchar = 16;		/* maximum number of characters for filename */
static int hms = 0;		/* 1 for output in hh:mm:ss dd:mm:ss */
static int nf = 0;
static int version = 0;		/* If 1, print only program name and version */


main (ac, av)
int ac;
char **av;
{
    char *str;
    int readlist = 0;
    char *lastchar;
    char filename[128];
    FILE *flist;
    char *listfile;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

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

    	case 'z':	/* Use AIPS classic WCS */
    	    setdefwcs(1);
    	    break;

	case '@':	/* List of files to be read */
	    readlist++;
	    listfile = ++str;
	    str = str + strlen (str) - 1;
	    av++;
	    ac--;

	default:
	    usage();
	    break;
	}
    }

    /* Find number of images to search and leave listfile open for reading */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMHEAD: List file %s cannot be read\n",
		     listfile);
	    usage();
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
	    ListWCS (filename);
	    if (verbose)
		printf ("\n");
	    }
	fclose (flist);
	}

    /* If no arguments left, print usage */
    if (ac == 0)
	usage();

    if (verbose)

    nf = 0;
    while (ac-- > 0) {
	char *fn = *av++;
	nf++;
	ListWCS (fn);
	}

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Print WCS part of FITS or IRAF image header\n");
    fprintf (stderr,"usage: wcshead [-htv] file.fit ...\n");
    fprintf (stderr,"  -h: print CRVALs as hh:mm:ss dd:mm:ss\n");
    fprintf (stderr,"  -t: print tab table output\n");
    fprintf (stderr,"  -v: verbose\n");
    fprintf (stderr,"  -z: Use AIPS classic WCS subroutines\n");
    exit (1);
}

static void
ListWCS (filename)

char	*filename;	/* FITS or IRAF image file name */
{
    int i;
    char str[256], temp[80];
    char rastr[32], decstr[32], fform[8];
    struct WorldCoor *wcs, *GetWCSFITS();

    wcs = GetWCSFITS (filename);
    if (nowcs (wcs))
	return;

    if (wcs->ctype[0][0] == (char) 0)
	return;
    if (tabout && nf == 1) {
	strcpy (str, "filename");
	for (i = 1; i < nchar - 8; i++) strcat (str, " ");
	strcat (str, "	naxis1	naxis2");
	strcat (str, "	ctype1  	ctype2  ");
	strcat (str, "	crval1	crval2	radecsys");
	strcat (str, "	crpix1	crpix2");
	strcat (str, "	cdelt1	cdelt2");
	strcat (str, "	crota2\n");
	strcat (str, "--------");
	for (i = 1; i < nchar - 8; i++) strcat (str, "-");
	strcat (str, "	------	------");
	strcat (str, "	------	------	--------");
	strcat (str, "	-------	-------");
	strcat (str, "	------	------");
	strcat (str, "	------	------");
	strcat (str, "	--------\n");
	printf ("%s", str);
	}

    sprintf (fform,"%%%d.%ds",nchar, nchar);
    if (tabout)
	sprintf (str, fform, filename);
    else
	sprintf (str, fform, filename);

    if (tabout)
	sprintf (temp, "	%.0f	%.0f", wcs->nxpix, wcs->nypix);
    else
	sprintf (temp, " %4.0f %4.0f", wcs->nxpix, wcs->nypix);
    strcat (str, temp);

    if (tabout)
	sprintf (temp, "	%s	%s", wcs->ctype[0], wcs->ctype[1]);
    else
	sprintf (temp, " %s %s", wcs->ctype[0], wcs->ctype[1]);
    strcat (str, temp);

    if (tabout) {
	if (hms) {
	    if (wcs->coorflip) {
		ra2str (rastr, 32, wcs->yref, ndec);
		dec2str (decstr, 32, wcs->xref, ndec-1);
		}
	    else {
		ra2str (rastr, 32, wcs->xref, ndec);
		dec2str (decstr, 32, wcs->yref, ndec-1);
		}
	    sprintf (temp, " %s %s %s", rastr, decstr, wcs->radecsys);
	    }
	else
	    sprintf (temp, "	%7.2f	%7.2f	%s", wcs->xref, wcs->yref, wcs->radecsys);
	}
    else {
	if (hms) {
	    if (wcs->coorflip) {
		ra2str (rastr, 32, wcs->yref, ndec);
		dec2str (decstr, 32, wcs->xref, ndec-1);
		}
	    else {
		ra2str (rastr, 32, wcs->xref, ndec);
		dec2str (decstr, 32, wcs->yref, ndec-1);
		}
	    sprintf (temp, " %s %s %s", rastr, decstr, wcs->radecsys);
	    }
	else
	    sprintf (temp, " %7.2f %7.2f %s", wcs->xref, wcs->yref, wcs->radecsys);
	}
    strcat (str, temp);

    if (tabout)
	sprintf (temp, "	%7.2f	%7.2f", wcs->xrefpix, wcs->yrefpix);
    else
	sprintf (temp, " %7.2f %7.2f", wcs->xrefpix, wcs->yrefpix);
    strcat (str, temp);

    if (tabout)
	sprintf (temp, "	%7.4f	%7.4f", 3600.0*wcs->xinc, 3600.0*wcs->yinc);
    else
	sprintf (temp, " %7.4f %7.4f", 3600.0*wcs->xinc, 3600.0*wcs->yinc);
    strcat (str, temp);

    if (tabout)
	sprintf (temp, "	%7.4f\n", wcs->rot);
    else
	sprintf (temp, " %7.4f\n", wcs->rot);
    strcat (str, temp);

    printf ("%s", str);

    wcsfree (wcs);

    return;
}
/* Feb 18 1998	New program
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Jul 10 1998	Add option to use AIPS classic WCS subroutines
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Apr  7 1999	Print lines all at once instead of one variable at a time
 * Jun  3 1999	Change PrintWCS to ListWCS to avoid name conflict
 * Oct 15 1999	Free wcs using wcsfree()
 * Oct 22 1999	Drop unused variables after lint
 * Nov 30 1999	Fix declaration of ListWCS()
 */
