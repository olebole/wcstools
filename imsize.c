/* File imsize.c
 * October 22, 1999
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
static void PrintWCS ();
extern void setsys();
extern void setcenter();
extern void setsecpix();
extern struct WorldCoor *GetFITSWCS();
extern char *GetFITShead();

static char coorsys[4];
static double size = 0.0;
static double frac = 0.0;

static int verbose = 0;		/* verbose/debugging flag */
static int dss = 0;		/* Flag to drop extra stuff for DSS */
static int dssc = 0;		/* Flag to drop extra stuff for DSS */
static double eqout = 0.0;
static double eqim = 0.0;
static int sysout = 0;
static int sysim = 0;
static int printepoch = 0;
static int printrange = 0;	/* Flag to print range rather than center */
static int version = 0;		/* If 1, print only program name and version */

main (ac, av)
int ac;
char **av;
{
    char *str, *str1;
    char rastr[16];
    char decstr[16];
    int readlist = 0;
    char *lastchar;
    char filename[128];
    FILE *flist;
    char *listfile;

    coorsys[0] = 0;

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

	case 'v':	/* more verbosity */
	    verbose++;
	    break;

	case 'b':	/* ouput B1950 (B1950) coordinates */
	    eqout = 1950.0;
	    sysout = WCS_B1950;
	    str1 = *(av+1);
	    if (*(str+1) || !strchr (str1,':'))
		strcpy (coorsys, "B1950");
	    else if (ac < 3)
		usage ();
	    else {
		setsys (WCS_B1950);
		strcpy (coorsys, "B1950");
		strcpy (rastr, *++av);
		ac--;
		strcpy (decstr, *++av);
		ac--;
		setcenter (rastr, decstr);
		}
	    break;

	case 'c':	/* Change output for DSS */
	    strcpy (coorsys, "J2000");
	    eqout = 2000.0;
	    sysout = WCS_J2000;
	    dssc++;
	    str1 = *(av+1);
	    if (!*(str+1) && (strchr (str1,'-') || strchr (str1,'+')) ) {
		size = atof (*++av);
		ac--;
		}
	    else if (!*(str+1) && strchr (str1,'x') ) {
		frac = atof (*++av+1);
		ac--;
		}
	    break;

	case 'd':	/* Change output for DSS */
	    strcpy (coorsys, "J2000");
	    eqout = 2000.0;
	    sysout = WCS_J2000;
	    dss++;
	    str1 = *(av+1);
	    if (!*(str+1) && (strchr (str1,'-') || strchr (str1,'+')) ) {
		size = atof (*++av);
		ac--;
		}
	    else if (!*(str+1) && strchr (str1,'x') ) {
		frac = atof (*++av+1);
		ac--;
		}
	    break;

	case 'e':	/* ouput epoch of plate */
	    printepoch++;
	    break;

	case 'j':	/* ouput J2000 (J2000) coordinates */
	    str1 = *(av+1);
	    eqout = 2000.0;
	    sysout = WCS_J2000;
	    if (*(str+1) || !strchr (str1,':'))
		strcpy (coorsys, "J2000");
	    else if (ac < 3)
		usage ();
	    else {
		setsys (WCS_J2000);
		strcpy (coorsys, "J2000");
		strcpy (rastr, *++av);
		ac--;
		strcpy (decstr, *++av);
		ac--;
		setcenter (rastr, decstr);
		}
	    break;

    	case 'p':	/* Initial plate scale in arcseconds per pixel */
    	    if (ac < 2)
    		usage();
    	    setsecpix (atof (*++av));
    	    ac--;
    	    break;

	case 'r':
	    printrange++;
	    break;

	case 'z':       /* Use AIPS classic WCS */
	    setdefwcs (1);
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
	    fprintf (stderr,"IMSIZE: List file %s cannot be read\n",
		     listfile);
	    usage ();
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
	    PrintWCS (filename);
	    if (verbose)
		printf ("\n");
	    }
	fclose (flist);
	}

    /* If no arguments left, print usage */
    if (ac == 0)
	usage ();

    while (ac-- > 0) {
	char *fn = *av++;
	PrintWCS (fn);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Print size of image in WCS and pixels\n");
    fprintf (stderr,"Usage: [-vcd] [-p scale] [-b ra dec] [-j ra dec] FITS or IRAF file(s)\n");
    fprintf (stderr,"  -b: Output B1950 (B1950) coordinates (optional center)\n");
    fprintf (stderr,"  -c: Format output without pixel dimensions (optional size change)\n");
    fprintf (stderr,"  -d: Format output as input to DSS getimage (optional size change)\n");
    fprintf (stderr,"  -e: Add epoch of image to output line\n");
    fprintf (stderr,"  -j: Output J2000 (J2000) coordinates (optional center)\n");
    fprintf (stderr,"  -p: Initial plate scale in arcsec per pixel (default 0)\n");
    fprintf (stderr,"  -r: Print range in RA and Dec\n");
    fprintf (stderr,"  -v: Verbose\n");
    fprintf (stderr,"  -z: use AIPS classic projections instead of WCSLIB\n");
    exit (1);
}


static void
PrintWCS (name)

char *name;

{
    char *header;		/* FITS image header */
    int nc;
    char fileroot[64];
    char *filename, *ext;
    int nax;
    int hp, wp, i, lfroot;
    double cra, cdec, dra, ddec, secpix;
    double xmin, xmax, ymin, ymax, dx, dy;
    struct WorldCoor *wcs;
    char *colon;
    char rstr[32], dstr[32], blanks[64];
    char ramin[32], ramax[32], decmin[32], decmax[32];

    /* Find root file name */
    if (strrchr (name,'/'))
	filename = strrchr (name,'/') + 1;
    else
	filename = name;
    ext = strsrch (filename, ".imh");
    if (ext == NULL)
	ext = strsrch (filename,".fit");
    if (ext == NULL)
	ext = strsrch (filename,".fts");
    if (ext != NULL)
	nc = ext - filename;
    else
	nc = strlen (filename);
    strncpy (fileroot, filename, nc);
    fileroot[nc] = 0;

    if (verbose) {
	fprintf (stderr,"Print World Coordinate System from ");
	if (isiraf (name))
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}

    if ((header = GetFITShead (name)) == NULL)
	return;

    /* Set image dimensions */
    nax = 0;
    if (hgeti4 (header,"NAXIS",&nax) < 1)
	return;
    else {
	if (hgeti4 (header,"NAXIS1",&wp) < 1)
	    return;
	else {
	    if (hgeti4 (header,"NAXIS2",&hp) < 1)
		return;
	    }
	}

    /* Read world coordinate system information from the image header */
    wcs = GetFITSWCS (name, header, verbose, &cra, &cdec, &dra, &ddec, &secpix,
		      &wp, &hp, &sysim, &eqim);
    if (nowcs (wcs)) {
	printf ("%s: No WCS for file, cannot compute image size\n", name);
	return;
	}

    /* Convert to desired output coordinates */
    wcscon (sysim, sysout, eqim, eqout, &cra, &cdec, wcs->epoch);

    /* Image center */
    ra2str (rstr, 16, cra, 3);
    dec2str (dstr, 16, cdec, 2);

    /* Image size in arcminutes */
    dra = 2.0 * dra * 60.0 * cos (degrad(cdec));
    ddec = 2.0 * ddec * 60.0;

    if (coorsys[0] == 0)
	wcscstr (coorsys, wcs->syswcs, wcs->equinox, 0.0);
    else
	wcsoutinit (wcs, coorsys);

    /* Print information */
    if (frac > 0.0) {
	dra = dra * frac;
	ddec = ddec * frac;
	}
    else if (size != 0.0) {
	dra = dra + size;
	ddec = ddec + size;
	}

    /* Print coverage of image in right ascension and declination */
    if (printrange) {
	wcsrange (wcs, &xmin, &xmax, &ymin, &ymax);
	ra2str (ramin, 32, xmin, 3);
	ra2str (ramax, 32, xmax, 3);
	dec2str (decmin, 32, ymin, 3);
	dec2str (decmax, 32, ymax, 3);
	dx = wcs->xinc * 3600.0;
	dy = wcs->yinc * 3600.0;
	strcpy (blanks, "                                       ");
	lfroot = strlen (fileroot);
	blanks[lfroot-1] = 0;
	printf ("%s RA:  %s -  %s %.4f arcsec/pix \n",
		fileroot, ramin, ramax, dy);
	printf ("%s Dec: %s - %s %.4f arcsec/pix %s\n",
		blanks, decmin, decmax, dy, coorsys);
	}

    /* Input for DSS GETIMAGE program */
    else if (dss) {
	for (i = 1; i < 3; i++) {
	    colon = strchr (rstr,':');
	    if (colon)
		*colon = ' ';
	    colon = strchr (dstr,':');
	    if (colon)
		*colon = ' ';
	    }
	printf ("%s %s %s ", fileroot, rstr, dstr);
	if (secpix > 0.0)
	    printf (" %.3f %.3f\n", dra, ddec);
	else
	    printf (" 10.0 10.0\n");
	}
    else {
	printf ("%s %s %s %s", fileroot, rstr, dstr, coorsys);
	if (secpix > 0.0)
	    printf (" %.3f\'x%.3f\'", dra, ddec);
	else if (dssc)
	    printf (" 10.000\'x10.000\'");
	if (!dssc) {
	    dx = wcs->xinc * 3600.0;
	    dy = wcs->yinc * 3600.0;
	    if (secpix > 0.0) {
		if (dx == dy)
		    printf (" %.4f\"/pix", secpix);
		else
		    printf (" %.4f/%.4f\"/pix", dx, dy);
		}
	    printf ("  %dx%d pix", wp, hp);
	    }
	if (printepoch)
	    printf (" (%.5f)",wcs->epoch);
	printf ("\n");
	}

    wcsfree (wcs);
    free (header);
    return;
}
/* Jul  9 1996	New program
 * Jul 18 1996	Update header reading
 * Jul 19 1996	Add option to change coordinate system
 * Aug  8 1996	Fix coordinate change option
 * Aug  9 1996	Add file name to output string
 * Aug 15 1996	Clean up image header reading code
 * Aug 27 1996	Add output format for DSS getimage
 * Aug 28 1996	Allow size increase or decrease if DSS-format output
 * Sep  3 1996	Add option to print DSS format with colons
 * Oct 16 1996  Rewrite to allow optional new center and use GetWCSFITS
 * Oct 17 1996	Do not print angular size and scale if not set
 * Oct 30 1996	Make equinox default to J2000 if not in image header
 * Dec 10 1996	Change equinox in getfitswcs call to double
 *
 * Feb 21 1997  Add optional epoch of image to output
 * Feb 24 1997  Read header using GetFITShead()
 * May 28 1997  Add option to read a list of filenames from a file
 * Sep  8 1997	Add option to change size by fraction of size or constant
 * Oct 14 1997	Add option to print RA and Dec range
 * Nov 13 1997	Print both plate scales if they are different
 * Dec 15 1997	Read IRAF 2.11 image format; fix bug if no WCS
 *
 * Jan 27 1998  Implement Mark Calabretta's WCSLIB
 * Jan 29 1998  Add -z for AIPS classic WCS projections
 * Feb 18 1998	Version 2.0: Full Calabretta WCS
 * Apr 24 1998	change coordinate setting to setsys() from setfk4()
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Jun 25 1998	Set WCS subroutine choice with SETDEFWCS()
 * Jul 24 1998	Drop unused variable irafheader
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Sep 17 1998	Add coordinate system to GetFITSWCS()
 * Sep 17 1998	Set coordinate system string using wcscstr()
 * Sep 29 1998	Change call to GetFITSWCS()
 * Oct 14 1998	Use isiraf() to determine file type
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Apr  7 1999	Add file name argument to GetFITSWCS
 * Jun 17 1999	Fix coordinate conversion
 * Oct 15 1999	Free wcs using wcsfree()
 * Oct 22 1999	Drop unused variables after lint
 */
