/* File newhead.c
 * April 7, 1999
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 */

/* Write a FITS header without any data.  Add information using edhead or sethead */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>

#include "fitsfile.h"
#include "wcs.h"

static void usage();
static void MakeHead ();

extern struct WorldCoor *GetFITSWCS();

static int verbose = 0;	/* verbose/debugging flag */
static int bitpix = 0;	/* number of bits per pixel (FITS code, 0=no image) */
static int fitsout = 0;	/* Output FITS file from IRAF input if 1 */
static int version = 0;		/* If 1, print only program name and version */
static int wcshead = 0;	/* If 1, add WCS information from command line */
static int nx = 100;	/* width of image in pixels */
static int ny = 100;	/* height of image in pixels */

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
    char rastr[32], decstr[32];
    double x, y;

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
    	case 'a':	/* Initial rotation angle in degrees */
    	    if (ac < 2)
    		usage();
    	    setrot (atof (*++av));
    	    ac--;
	    wcshead++;
    	    break;

    	case 'b':	/* initial coordinates on command line in B1950 */
    	    if (ac < 3)
    		usage();
	    setsys (WCS_B1950);
	    strcpy (rastr, *++av);
	    ac--;
	    strcpy (decstr, *++av);
	    ac--;
	    setcenter (rastr, decstr);
	    wcshead++;
    	    break;

    	case 'j':	/* center coordinates on command line in J2000 */
    	    if (ac < 3)
    		usage();
	    setsys (WCS_J2000);
	    strcpy (rastr, *++av);
	    ac--;
	    strcpy (decstr, *++av);
	    ac--;
	    setcenter (rastr, decstr);
	    wcshead++;
    	    break;

    	case 'o':	/* Number of bits per pixel */
    	    if (ac < 2)
    		usage ();
    	    bitpix = (int) atof (*++av);
    	    ac--;
    	    break;

    	case 'p':	/* Initial plate scale in arcseconds per pixel */
    	    if (ac < 2)
    		usage();
    	    setsecpix (atof (*++av));
    	    ac--;
	    if (ac > 1 && isnum (*(av+1))) {
		setsecpix2 (atof (*++av));
		ac--;
		}
	    wcshead++;
    	    break;

	case 's':	/* size of image in X and Y pixels */
	    if (ac < 3)
		usage();
	    nx = atoi (*++av);
	    ac--;
	    ny = atoi (*++av);
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
	    wcshead++;
    	    break;

	case '@':	/* List of files to be created */
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

    /* Process files in file of filenames */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMROT: List file %s cannot be read\n",
		     listfile);
	    usage ();
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
	    MakeHead (filename);
	    if (verbose)
		printf ("\n");
	    }
	fclose (flist);
	}

    /* If no arguments left, print usage */
    if (ac == 0)
	usage ();

    /* Process files on command line */
    else {
	while (ac-- > 0) {
    	    char *fn = *av++;
    	    MakeHead (fn);
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
    fprintf (stderr,"Make dataless FITS image header files\n");
    fprintf(stderr,"Usage: [-v] [-a degrees] [-p scale] [-b ra dec] [-j ra dec] [-s nx ny]\n");
    fprintf(stderr,"       [-x x y] file.fits ...\n");
    fprintf(stderr,"  -a: initial rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b: initial center in B1950 (FK4) RA and Dec\n");
    fprintf(stderr,"  -j: initial center in J2000 (FK5) RA and Dec\n");
    fprintf(stderr,"  -o: output pixel size in bits (FITS code, default=0)\n");
    fprintf(stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -s: size of image in x and y pixels (default 100x100)\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -x: X and Y coordinates of reference pixel (default is center)\n");
    exit (1);
}

static void
MakeHead (name)
char *name;
{
    char *image;	/* FITS image */
    char *header;	/* FITS header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    double cra;		/* Center right ascension in degrees (returned) */
    double cdec;	/* Center declination in degrees (returned) */
    double dra;		/* Right ascension half-width in degrees (returned) */
    double ddec;	/* Declination half-width in degrees (returned) */
    double secpix;	/* Arcseconds per pixel (returned) */
    int wp;		/* Image width in pixels (returned) */
    int hp;		/* Image height in pixels (returned) */
    int sysout;		/* Coordinate system to return (0=image, returned) */
    double  *eqout;	/* Equinox to return (0=image, returned) */
    int i;
    char temp[8];
    struct WorldCoor *wcs;
    FILE *diskfile;

    if (verbose) {
	fprintf (stderr,"Create header as ");
	fprintf (stderr,"FITS file %s\n", name);
	}

    /* Make sure that no existing file is overwritten */
    if ((diskfile = fopen (name, "r")) != NULL) {
	if (verbose)
	    fprintf (stderr,"NEWHEAD: FITS file %s exists, no new file written\n",
		     name);
	fclose (diskfile);
	return;
	}

    lhead = 14400;
    header = malloc (lhead);
    strcpy (header, "END ");
    for (i = 4; i < lhead; i++)
	header[i] = ' ';
    hputl (header, "SIMPLE", 1);
    hputi4 (header, "BITPIX", bitpix);
    hputi4 (header, "NAXIS", 2);
    hputi4 (header, "NAXIS1", nx);
    hputi4 (header, "NAXIS2", ny);

    if (bitpix != 0) {
	if (verbose)
	    fprintf (stderr, "NEWHEAD: %d bits/pixel\n", bitpix);
	}

    if (wcshead)
	wcs = GetFITSWCS (name,header,verbose,&cra,&cdec,&dra,&ddec,&secpix,
			  &wp,&hp,&sysout,&eqout);

    if (fitswimage (name, header, image) > 0 && verbose)
	printf ("%s: written successfully.\n", name);

    free (header);
    return;
}
/* Jan  4 1999	New program
 * Apr  7 1999	Add filename argument to GetFITSWCS
 */
