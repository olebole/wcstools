/* File imsize.c
 * February 21, 1997
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
static void PrintWCS ();
extern void fk524e();
extern void fk425e();
extern void setfk4();
extern void setcenter();
extern void setsecpix();
extern struct WorldCoor *GetFITSWCS();
extern char *GetFITShead();

static char coorsys[4];
static double size = 0.0;

static int verbose = 0;		/* verbose/debugging flag */
static int dss = 0;		/* Flag to drop extra stuff for DSS */
static int dssc = 0;		/* Flag to drop extra stuff for DSS */
static int ieq = 0;
static double equinox = 0.0;
static int printepoch = 0;

main (ac, av)
int ac;
char **av;
{
    char *str, *str1;
    char rastr[16];
    char decstr[16];

    coorsys[0] = 0;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {

	case 'v':	/* more verbosity */
	    verbose++;
	    break;

	case 'b':	/* ouput B1950 (B1950) coordinates */
	    strcpy (coorsys, "B1950");
	    ieq = 1950;
	    equinox = 1950.0;
	    str1 = *(av+1);
	    if (*(str+1) || !strchr (str1,':'))
		setfk4 ();
	    else if (ac < 3)
		usage ();
	    else {
		setfk4 ();
		strcpy (rastr, *++av);
		ac--;
		strcpy (decstr, *++av);
		ac--;
		setcenter (rastr, decstr);
		}
	    break;

	case 'c':	/* Change output for DSS */
	    strcpy (coorsys, "J2000");
	    dssc++;
	    str1 = *(av+1);
	    if (!*(str+1) && (strchr (str1,'-') || strchr (str1,'+')) ) {
		size = atof (*++av);
		ac--;
		}
	    break;

	case 'd':	/* Change output for DSS */
	    strcpy (coorsys, "J2000");
	    dss++;
	    str1 = *(av+1);
	    if (!*(str+1) && (strchr (str1,'-') || strchr (str1,'+')) ) {
		size = atof (*++av);
		ac--;
		}
	    break;

	case 'e':	/* ouput epoch of plate */
	    printepoch++;
	    break;

	case 'j':	/* ouput J2000 (J2000) coordinates */
	    str1 = *(av+1);
	    ieq = 2000;
	    equinox = 2000.0;
	    if (*(str+1) || !strchr (str1,':'))
		strcpy (coorsys, "J2000");
	    else if (ac < 3)
		usage ();
	    else {
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

	default:
	    usage();
	    break;
	}
    }

    /* now there are ac remaining file names starting at av[0] */
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
    fprintf (stderr,"Print size of image in WCS and pixels\n");
    fprintf(stderr,"Usage: [-vcd] [-p scale] [-b ra dec] [-j ra dec] FITS or IRAF file(s)\n");
    fprintf(stderr,"  -b: Output B1950 (B1950) coordinates (optional center)\n");
    fprintf(stderr,"  -c: Format output without pixel dimensions (optional size change)\n");
    fprintf(stderr,"  -d: Format output as input to DSS getimage (optional size change)\n");
    fprintf(stderr,"  -e: Add epoch of image to output line\n");
    fprintf(stderr,"  -j: Output J2000 (J2000) coordinates (optional center)\n");
    fprintf(stderr,"  -p: Initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}


static void
PrintWCS (name)

char *name;

{
    char *header;		/* FITS image header */
    int lhead, nbhead, nc;
    int iraffile;		/* 1 if IRAF image */
    int *irafheader;		/* IRAF image header */
    char fileroot[64];
    char *filename, *ext;
    int nax;
    int hp, wp, i;
    double cra, cdec, dra, ddec, secpix;
    struct WorldCoor *wcs;
    char *colon;
    char rstr[64], dstr[64];

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
	if (strsrch (name,".imh") != NULL)
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
    wcs = GetFITSWCS (header, verbose, &cra, &cdec, &dra, &ddec, &secpix,
                &wp, &hp, equinox);

    /* Image center */
    ra2str (rstr, cra, 3);
    dec2str (dstr, cdec, 2);

    /* Image size in arcminutes */
    dra = 2.0 * dra * 60.0 * cos (degrad(cdec));
    ddec = 2.0 * ddec * 60.0;

    if (coorsys[0] == 0) {
	ieq = -1;
	if (iswcs (wcs)) {
	    if (wcs->equinox == 1950.0)
		strcpy (coorsys, "B1950");
	    else if (wcs->equinox == 2000.0)
		strcpy (coorsys, "J2000");
	    else
		sprintf (coorsys,"%6.1f",wcs->equinox);
	    }
	else if (hgetr8 (header,"EQUINOX",&equinox)) {
	    if (equinox == 0)
		equinox = 1950.0;
	    if (equinox == 1950.0)
		strcpy (coorsys, "B1950");
	    else if (ieq == 2000)
		strcpy (coorsys, "J2000");
	    else
		sprintf (coorsys,"%6.1f",equinox);
	    }
	else if (hgetr8 (header,"EPOCH",&equinox)) {
	    if (equinox == 0)
		equinox = 1950.0;
	    if (equinox == 1950.0)
		strcpy (coorsys, "B1950");
	    else if (equinox == 2000.0)
		strcpy (coorsys, "J2000");
	    else
		sprintf (coorsys,"%6.1f",equinox);
	    }
	else
	    strcpy (coorsys, "J2000");
	}

    /* Print information */
    if (dss) {
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
	    printf (" %.3f %.3f\n", dra+size, ddec+size);
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
	    if (secpix > 0.0)
		printf (" %.3f \"/pix", secpix);
	    printf ("  %dx%d pix", wp, hp);
	    }
	if (printepoch)
	    printf (" (%.5f)",wcs->epoch);
	printf ("\n");
	}

    free (wcs);
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
 * Feb 21 1997  Read header using GetFITShead()
 */
