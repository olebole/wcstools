/* File immatch.c
 * December 8, 1997
 * By Doug Mink, after Elwood Downey
 * (Harvard-Smithsonian Center for Astrophysics)
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

static void usage();
static void FitWCS ();

static int verbose = 0;		/* verbose/debugging flag */
static int overwrite = 0;	/* allow overwriting of input image file */
static int rot = 0;
static int mirror = 0;
static int bitpix = 0;
static int imsearch = 1;	/* set to 0 if image catalog provided */

extern int RotFITS ();
extern int SetWCSFITS ();
extern int DelWCSFITS();
extern int PrintWCS();
extern void settolerance ();
extern void setreflim ();
extern void setrot ();
extern void setnfit ();
extern void setsecpix ();
extern void setcenter ();
extern void setfk4 ();
extern void setminb ();
extern void setmaxcat ();
extern void setstarsig ();
extern void setclass();
extern void setplate();
extern void setrefcat();
extern void setimcat();
extern void setbmin();
extern void setfrac();
extern void setrefpix();
extern void setwcstype();
extern void setfitwcs();

main (ac, av)
int ac;
char **av;
{
    char *str;
    double bmin, maglim1, maglim2;
    char rastr[16];
    char decstr[16];
    int readlist = 0;
    char *lastchar;
    char filename[128];
    FILE *flist;
    char *listfile;
    double x, y;

    setfitwcs (0);

    /* Decode arguments */
    for (av++; --ac > 0 && (*(str = *av)=='-' || *str == '@'); av++) {
	char c;
	if (*str == '@')
	    str = str - 1;
	while ((c = *++str) != 0)
    	switch (c) {
    	case 'a':	/* Initial rotation angle in degrees */
    	    if (ac < 2)
    		usage();
    	    setrot (atof (*++av));
    	    ac--;
    	    break;

    	case 'b':	/* initial coordinates on command line in B1950 */
	    setfk4 ();
    	    if (ac < 3)
    		usage();
	    strcpy (rastr, *++av);
	    ac--;
	    strcpy (decstr, *++av);
	    ac--;
	    setcenter (rastr, decstr);
    	    break;

	case 'c':       /* Set reference catalog */
	    if (ac < 2)
		usage();
	    setrefcat (*++av);
	    ac--;
	    break;

	case 'd':	/* Read image star positions from DAOFIND file */
	    if (ac < 2)
		usage();
	    setimcat (*++av);
	    imsearch = 0;
	    ac--;
	    break;

	case 'e':	/* Set WCS projection
	    if (ac < 2)
		usage();
	    setwcsproj (*++av);
	    ac--;
	    break; */
	    
    	case 'j':	/* center coordinates on command line in J2000 */
    	    if (ac < 3)
    		usage();
	    strcpy (rastr, *++av);
	    ac--;
	    strcpy (decstr, *++av);
	    ac--;
	    setcenter (rastr, decstr);
    	    break;

	case 'g':	/* Guide Star object class */
    	    if (ac < 2)
    		usage();
    	    setclass ((int) atof (*++av));
    	    ac--;
    	    break;

	case 'h':	/* Maximum number of reference stars */
    	    if (ac < 2)
    		usage();
    	    setmaxcat ((int) atof (*++av));
    	    ac--;
    	    break;

	case 'i':       /* Image star minimum peak value */
	    if (ac < 2)
		usage();
	    bmin = atof (*++av);
	    if (bmin < 0)
		setstarsig (-bmin);
	    else
		setbmin (bmin);
	    ac--;
	    break;

    	case 'l':	/* Left-right reflection before rotating */
	    mirror = 1;
    	    break;

    	case 'm':	/* Limiting reference star magnitude */
    	    if (ac < 2)
    		usage();
	    maglim1 = -99.0;
    	    maglim2 = atof (*++av);
    	    ac--;
	    if (ac > 1 && isnum (*(av+1))) {
		maglim1 = maglim2;
		maglim2 = atof (*++av);
		ac--;
		}
    	    setreflim (maglim1, maglim2);
    	    break;

    	case 'p':	/* Initial plate scale in arcseconds per pixel */
    	    if (ac < 2)
    		usage();
    	    setsecpix (atof (*++av));
    	    ac--;
    	    break;

    	case 'r':	/* Angle in degrees to rotate before fitting */
    	    if (ac < 2)
    		usage();
    	    rot = (int) atof (*++av);
    	    ac--;
    	    break;

    	case 's':   /* this fraction more image stars than GSC or vice versa */
    	    if (ac < 2)
    		usage();
    	    setfrac (atof (*++av));
    	    ac--;
    	    break;

    	case 't':	/* +/- this many pixels is a hit */
    	    if (ac < 2)
    		usage();
    	    settolerance (atof (*++av));
    	    ac--;
    	    break;

	case 'u':	/* UJ Catalog plate number */
    	    if (ac < 2)
    		usage();
    	    setplate ((int) atof (*++av));
    	    ac--;
    	    break;

    	case 'v':	/* more verbosity */
    	    verbose++;
    	    break;

	case 'x':	/* X and Y coordinates of reference pixel */
	    if (ac < 3)
		usage();
	    x = atof (*++av);
	    ac--;
	    y = atof (*++av);
	    ac--;
    	    setrefpix (x, y);
    	    break;

	case 'y':	/* Multiply dimensions of image by fraction */
	    if (ac < 2)
		usage();
	    setimfrac (atof (*++av));
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
	    fprintf (stderr,"IMMATCH: List file %s cannot be read\n",
		     listfile);
	    usage ();
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
	    FitWCS (filename);
	    if (verbose)
		printf ("\n");
	    }
	fclose (flist);
	}

    /* If no arguments left, print usage */
    if (ac == 0)
	usage();

    while (ac-- > 0) {
    	char *fn = *av++;
    	if (verbose)
    	    printf ("%s:\n", fn);
    	FitWCS (fn);
    	if (verbose)
    	    printf ("\n");
	}

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"\n");
    fprintf (stderr,"Match catalog and image stars from WCS in image file\n");
    fprintf(stderr,"Usage: [-vl] [-m mag] [-n frac] [-s mode] [-g class] [-h maxref] [-i peak]\n");
    fprintf(stderr,"       [-c catalog] [-p scale] [-b ra dec] [-j ra dec] [-r deg] [-t tol] [-x x y] [-y frac]\n");
    fprintf(stderr,"       FITS or IRAF file(s)\n");
    fprintf(stderr,"  -a: initial rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b: initial center in B1950 (FK4) RA and Dec\n");
    fprintf(stderr,"  -c: reference catalog (gsc, uac, ujc, tab table file\n");
    fprintf(stderr,"  -d: Use following DAOFIND output catalog instead of search\n");
    fprintf(stderr,"  -e: WCS type (TAN default)\n");
    fprintf(stderr,"  -g: Guide Star Catalog class (-1=all,0,3 (default -1)\n");
    fprintf(stderr,"  -h: maximum number of reference stars to use (10-200, default 25\n");
    fprintf(stderr,"  -i: minimum peak value for star in image (<0=-sigma)\n");
    fprintf(stderr,"  -j: initial center in J2000 (FK5) RA and Dec\n");
    fprintf(stderr,"  -l: reflect left<->right before rotating and fitting\n");
    fprintf(stderr,"  -m: initial reference catalog magnitude limits\n");
    fprintf(stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -r: rotation angle in degrees before fitting (default 0)\n");
    fprintf(stderr,"  -s: use this fraction extra stars (default 1.0)\n");
    fprintf(stderr,"  -t: offset tolerance in pixels (default 20)\n");
    fprintf(stderr,"  -u: USNO catalog single plate number to accept\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -x: X and Y coordinates of reference pixel (default is center)\n");
    fprintf(stderr,"  -y: multiply image dimensions by this for search (default is 1)\n");
    exit (1);
}


static void
FitWCS (name)
char *name;
{
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    int iraffile;		/* 1 if IRAF image */
    char *image;		/* Image */
    char *header;		/* FITS header */
    int *irafheader;		/* IRAF image header */
    char newname[64];		/* Name for revised image */
    char pixname[64];		/* Pixel file name for revised image */
    char temp[16];
    char *ext;
    char *fname;
    int lext, lname;
    int newimage;

    image = NULL;

    /* Open IRAF image if .imh extension is present */
    if (strsrch (name,".imh") != NULL) {
	iraffile = 1;
	if ((irafheader = irafrhead (name, &lhead)) != NULL) {
	    header = iraf2fits (name, irafheader, lhead, &nbhead);
	    if (header == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",name);
		free (irafheader);
		return;
		}
	    if (imsearch || rot || mirror) {
		if ((image = irafrimage (header)) == NULL) {
		    hgets (header,"PIXFILE", 64, pixname);
		    fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
		    free (irafheader);
		    free (header);
		    return;
		    }
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", name);
	    return;
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	if ((header = fitsrhead (name, &lhead, &nbhead)) != NULL) {
	    if (imsearch || rot || mirror) {
		if ((image = fitsrimage (name, nbhead, header)) == NULL) {
		    fprintf (stderr, "Cannot read FITS image %s\n", name);
		    free (header);
		    return;
		    }
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", name);
	    return;
	    }
	}

    if (verbose) {
	fprintf (stderr,"Matching catalog to ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}

    /* Print existing WCS headers and check for permission to overwrite */
    (void) PrintWCS (header, verbose);

    /* Rotate and/or reflect image */
    if (imsearch  && (rot != 0 || mirror)) {
	if (RotFITS (name, header, &image, rot, mirror, bitpix, verbose)) {
	    fprintf (stderr,"Image %s could not be rotated\n", name);
	    return;
	    }
	}

    if (SetWCSFITS (name, header, image, verbose)) {
	if (verbose)
	    (void) PrintWCS (header, verbose);	/* print new WCS */
	}

    free (header);
    if (iraffile)
	free (irafheader);
    if (image != NULL)
	free (image);
    return;
}


char *
{
}
/* Nov  6 1997	New program based on IMWCS
 * Nov 17 1997	Add optional second magnitude limit
 * Dec  8 1997	Fix bug in setting nominal WCS
 */
