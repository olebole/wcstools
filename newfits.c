/* File newfits.c
 * November 1, 1999
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 */

/* Write a FITS header without any data or a blank FITS image.
 * Add information using edhead or sethead */

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
static void MakeFITS ();
extern void setcenter();
extern void setsys();
extern void setrot();
extern void setsecpix();
extern void setsecpix2();
extern void setrefpix();
extern void setcdelt();
extern struct WorldCoor *GetFITSWCS();

static int verbose = 0;	/* verbose/debugging flag */
static int bitpix = 0;	/* number of bits per pixel (FITS code, 0=no image) */
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
    for (av++; --ac > 0; av++) {

	/* Set RA, Dec, and equinox if WCS-generated argument */
	if (strsrch (*av,":") != NULL) {
	    if (ac < 3)
		usage();
	    else {
		strcpy (rastr, *av);
		ac--;
		strcpy (decstr, *++av);
		setcenter (rastr, decstr);
		ac--;
		setsys (wcscsys (*++av));
		}
	    }

	/* Set decimal degree center */
	else if (isnum (*av)) {
	    if (ac < 3)
		usage();
	    else {
		strcpy (rastr, *++av);
		ac--;
		strcpy (decstr, *++av);
		setcenter (rastr, decstr);
		ac--;
		setsys (wcscsys (*++av));
		}
	    }

	/* Read list of files to create from file */
	else if (*(str = *av) == '@') {
	    readlist = 1;
	    listfile = *av + 1;
	    }

	/* Otherwise, read command */
	else if (*(str = *av) == '-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

    		case 'a':	/* Initial rotation angle in degrees */
    		    if (ac < 2)
    			usage();
    		    setrot (atof (*++av));
    		    ac--;
		    wcshead++;
    		    break;
	
    		case 'b':	/* Reference pixel coordinates in B1950 */
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

    		case 'd':	/* Set CDELTn, CROTAn instead of CD matrix */
		    setcdelt();
    		    break;
	
    		case 'j':	/* Reference pixel coordinates in J2000 */
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
	
    		case 'p':	/* Initial plate scale in arcseconds / pixel */
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
	
		case 's':	/* Size of image in X and Y pixels */
		    if (ac < 3)
			usage();
		    nx = atoi (*++av);
		    ac--;
		    ny = atoi (*++av);
		    ac--;
    		    break;
	
    		case 'v':	/* More verbosity */
    		    verbose++;
    		    break;
	
		case 'x':	/* Reference pixel X and Y coordinates */
		    if (ac < 3)
			usage();
		    x = atof (*++av);
		    ac--;
		    y = atof (*++av);
		    ac--;
    		    setrefpix (x, y);
		    wcshead++;
    		    break;
	
    		default:
    		    usage();
    		    break;
    		}
	    }
	else {
    	    MakeFITS (*av);
    	    if (verbose)
    		printf ("\n");
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
	    MakeFITS (filename);
	    if (verbose)
		printf ("\n");
	    }
	fclose (flist);
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
    fprintf(stderr,"  -d: set CDELTn, CROTAn instead of CD matrix\n");
    fprintf(stderr,"  -j: initial center in J2000 (FK5) RA and Dec\n");
    fprintf(stderr,"  -o: output pixel size in bits (FITS code, default=0)\n");
    fprintf(stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -s: size of image in x and y pixels (default 100x100)\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -x: X and Y coordinates of reference pixel (default is center)\n");
    exit (1);
}

static void
MakeFITS (name)
char *name;
{
    char *image;	/* FITS image */
    char *header;	/* FITS header */
    int lhead;		/* Maximum number of bytes in FITS header */
    double cra;		/* Center right ascension in degrees (returned) */
    double cdec;	/* Center declination in degrees (returned) */
    double dra;		/* Right ascension half-width in degrees (returned) */
    double ddec;	/* Declination half-width in degrees (returned) */
    double secpix;	/* Arcseconds per pixel (returned) */
    int wp;		/* Image width in pixels (returned) */
    int hp;		/* Image height in pixels (returned) */
    int sysout=0;	/* Coordinate system to return (0=image, returned) */
    double eqout=0.0;	/* Equinox to return (0=image, returned) */
    int i, nbimage;
    struct WorldCoor *wcs;
    FILE *diskfile;

    if (verbose) {
	if (bitpix != 0)
	    fprintf (stderr,"Create ");
	else
	    fprintf (stderr,"Create header as ");
	fprintf (stderr,"FITS file %s\n", name);
	}

    /* Make sure that no existing file is overwritten */
    if ((diskfile = fopen (name, "r")) != NULL) {
	fprintf (stderr,"NEWFITS: FITS file %s exists, no new file written\n",
		     name);
	fclose (diskfile);
	return;
	}

    lhead = 14400;
    header = (char *) calloc (1, lhead);
    strcpy (header, "END ");
    for (i = 4; i < lhead; i++)
	header[i] = ' ';
    hlength (header, 14400);
    hputl (header, "SIMPLE", 1);
    hputi4 (header, "BITPIX", bitpix);
    hputi4 (header, "NAXIS", 2);
    hputi4 (header, "NAXIS1", nx);
    hputi4 (header, "NAXIS2", ny);


    /* Set up blank image, if bitpix is non-zero */
    if (bitpix != 0) {
	if (bitpix < 0)
	    nbimage = -bitpix * nx * ny;
	else
	    nbimage =  bitpix * nx * ny;
	image = calloc (nbimage, 1);
	if (verbose)
	    fprintf (stderr, "NEWFITS: %d bits/pixel\n", bitpix);
	}
    else
	image = NULL;

    /* Initialize header */
    if (wcshead) {
	wcs = GetFITSWCS (name,header,verbose,&cra,&cdec,&dra,&ddec,&secpix,
			  &wp,&hp,&sysout,&eqout);
	wcsfree (wcs);
	}

    if (fitswimage (name, header, image) > 0 && verbose)
	printf ("%s: written successfully.\n", name);

    free (header);
    if (image != NULL)
	free (image);
    return;
}
/* Jan  4 1999	New program
 * Apr  7 1999	Add filename argument to GetFITSWCS
 * May 13 1999	Change name to NEWFITS and add option to write blank image
 * Jun  3 1999	Allow center to be set as standard WCSTools coordinate string
 * Jun  3 1999	Move command line filenames into processing loop
 * Oct 15 1999	Free wcs using wcsfree()
 * Oct 22 1999	Drop unused variables after lint
 * Nov  1 1999	Set header length after creating it; add option for CDELTn
 */
