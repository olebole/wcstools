/* File imsize.c
 * October 31, 1996
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

static char coorsys[4];
static double size = 0.0;

static int verbose = 0;		/* verbose/debugging flag */
static int dss = 0;		/* Flag to drop extra stuff for DSS */
static int dssc = 0;		/* Flag to drop extra stuff for DSS */
static int ieq = 0;

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

	case 'j':	/* ouput J2000 (J2000) coordinates */
	    str1 = *(av+1);
	    ieq = 2000;
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
    fprintf(stderr,"  -b: output B1950 (B1950) coordinates (optional center)\n");
    fprintf(stderr,"  -c: format output without pixel dimensions (optional size change)\n");
    fprintf(stderr,"  -d: format output as input to DSS getimage (optional size change)\n");
    fprintf(stderr,"  -j: output J2000 (J2000) coordinates (optional center)\n");
    fprintf(stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -v: verbose\n");
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

    /* Open IRAF header if .imh extension is present */
    ext = strsrch (name, ".imh");
    if (ext) {
	iraffile = 1;
	if ((irafheader = irafrhead (name, &lhead))) {
	    header = iraf2fits (name, irafheader, lhead, &nbhead);
	    free (irafheader);
	    if (!header) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",name);
		return;
		}
	    if (strrchr (name,'/'))
		filename = strrchr (name,'/') + 1;
	    else
		filename = name;
	    nc = ext - filename;
	    strncpy (fileroot, filename, nc);
	    fileroot[nc] = 0;
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", name);
	    return;
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	if ((header = fitsrhead (name, &lhead, &nbhead))) {
	    if (strrchr (name,'/'))
		filename = strrchr (name,'/') + 1;
	    else
		filename = name;
	    ext = strsrch (name,".fit");
	    if (!ext)
		ext = strsrch (name,".fts");
	    if (ext)
		nc = ext - filename;
	    strncpy (fileroot, filename, nc);
	    fileroot[nc] = 0;
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", name);
	    return;
	    }
	}
    if (verbose) {
	fprintf (stderr,"Print World Coordinate System from ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}

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
                &wp, &hp, ieq);

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
	else if (hgeti4 (header,"EQUINOX",&ieq)) {
	    if (ieq == 0)
		ieq = 1950;
	    if (ieq == 1950)
		strcpy (coorsys, "B1950");
	    else if (ieq == 2000)
		strcpy (coorsys, "J2000");
	    else
		sprintf (coorsys,"%4d",ieq);
	    }
	else if (hgeti4 (header,"EPOCH",&ieq)) {
	    if (ieq == 0)
		ieq = 1950;
	    if (ieq == 1950)
		strcpy (coorsys, "B1950");
	    else if (ieq == 2000)
		strcpy (coorsys, "J2000");
	    else
		sprintf (coorsys,"%4d",ieq);
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
	    printf ("  %dx%d pix\n", wp, hp);
	    }
	else
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
 */
