/* File imsize.c
 * September 3, 1996
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
static int PrintFITSWCS ();
static void PrintWCS ();
extern void fk524e();
extern void fk425e();

static char coorsys0[4];
static double cra0 = 0.0;
static double cdec0 = 0.0;
static double size = 0.0;

static int verbose = 0;		/* verbose/debugging flag */
static int dss = 0;		/* Flag to drop extra stuff for DSS */
static int dssc = 0;		/* Flag to drop extra stuff for DSS */

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str, *str1;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {

	case 'v':	/* more verbosity */
	    verbose++;
	    break;

	case 'b':	/* ouput FK4 (B1950) coordinates */
	    strcpy (coorsys0, "FK4");
	    str1 = *(av+1);
	    if (*(str+1) || !strchr (str1,':'))
		strcpy (coorsys0, "FK4");
	    else if (ac < 3)
		usage (progname);
	    else {
		cra0 = str2ra (*++av);
		ac--;
		cdec0 = str2dec (*++av);
		ac--;
		}
	    break;

	case 'c':	/* Change output for DSS */
	    strcpy (coorsys0, "FK5");
	    dssc++;
	    str1 = *(av+1);
	    if (!*(str+1) && (strchr (str1,'-') || strchr (str1,'+')) ) {
		size = atof (*++av);
		ac--;
		}
	    break;

	case 'd':	/* Change output for DSS */
	    strcpy (coorsys0, "FK5");
	    dss++;
	    str1 = *(av+1);
	    if (!*(str+1) && (strchr (str1,'-') || strchr (str1,'+')) ) {
		size = atof (*++av);
		ac--;
		}
	    break;

	case 'j':	/* ouput FK5 (J2000) coordinates */
	    str1 = *(av+1);
	    if (*(str+1) || !strchr (str1,':'))
		strcpy (coorsys0, "FK5");
	    else if (ac < 3)
		usage (progname);
	    else {
		strcpy (coorsys0, "FK5");
		cra0 = str2ra (*++av);
		ac--;
		cdec0 = str2dec (*++av);
		ac--;
		}
	    break;

	default:
	    usage(progname);
	    break;
	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage (progname);

    while (ac-- > 0) {
	char *fn = *av++;
	PrintWCS (fn);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"Print size of image in WCS and pixels\n");
    fprintf(stderr,"%s: usage: [-v] file.fit ...\n", progname);
    fprintf(stderr,"  -b: output B1950 (FK4) coordinates (optional center)\n");
    fprintf(stderr,"  -d: format output as input to DSS getimage (optional size change)\n");
    fprintf(stderr,"  -j: output J2000 (FK5) coordinates (optional center)\n");
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
	    fprintf (stderr,"IRAF image file %s", name);
	else
	    fprintf (stderr,"FITS image file %s", name);
	}

    if (!PrintFITSWCS (header, fileroot))
	fprintf (stderr,"%s: no WCS fields found\n", name);

    free (header);
    return;
}


static int
PrintFITSWCS (header, filename)

char	*header;	/* Image FITS header */
char	*filename;	/* Image file name */
{
    int nax;
    int hp, wp, i;
    double cra, cdec, dra, ddec, secpix;
    struct WorldCoor *wcs;
    char *colon;
    char rstr[64], dstr[64], coorsys[4];


    /* Set image dimensions */
    nax = 0;
    if (hgeti4 (header,"NAXIS",&nax) < 1)
	return (0);
    else {
	if (hgeti4 (header,"NAXIS1",&wp) < 1)
	    return (0);
	else {
	    if (hgeti4 (header,"NAXIS2",&hp) < 1)
		return (0);
	    }
	}

    /* Find center for pre-existing WCS, if there is one */
    wcs = wcsinit (header);
    if (iswcs (wcs))
	wcssize (wcs, &cra, &cdec, &dra, &ddec);
    else
	return (0);

    /* Reset center if asked for */
    if (cra0 > 0.0)
	cra = cra0;
    if (cdec0 != 0.0)
	cdec = cdec0;
    if (cra0 > 0.0 || cdec0 != 0.0) {
	if (coorsys0[1])
	    wcsshift (wcs,cra,cdec,coorsys0);
	else
	    wcsshift (wcs,cra,cdec,wcs->radecsys);
	}

    if (coorsys0[2] == '5' && wcs->radecsys[2] == '4') {
	(void) fk425e (&cra, &cdec, wcs->epoch);
	strcpy (coorsys, coorsys0);
	}
    else if (coorsys0[2] == '4' && wcs->radecsys[2] == '5') {
	(void) fk524e (&cra, &cdec, wcs->epoch);
	strcpy (coorsys, coorsys0);
	}
    else
	strcpy (coorsys, wcs->radecsys);

    /* Plate scale from WCS if it is present */
    secpix = 3600.0 * 2.0 * ddec / (double) hp;

    /* Image size from header */
    ra2str (rstr, cra, 3);
    dec2str (dstr, cdec, 2);
    dra = 2.0 * dra * 60.0 * cos (degrad(cdec));
    ddec = 2.0 * ddec * 60.0;

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
	printf ("%s %s %s ", filename, rstr, dstr);
	printf (" %.3f %.3f\n", dra+size, ddec+size);
	}
    else {
	printf ("%s %s %s %s", filename, rstr, dstr, coorsys);
	printf (" %.3f\'x%.3f\'", dra, ddec);
	if (!dssc) {
	    printf (" %.3f \"/pix", secpix);
	    printf ("  %dx%d pix\n", wp, hp);
	    }
	else
	    printf ("\n");
	}

    free (wcs);
    return (1);
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
 */
