/* File imwcs.c
 * October 11, 1996
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
static int PrintWCS ();

static int verbose = 0;		/* verbose/debugging flag */
static int writeheader = 0;	/* write header fields; else read-only */
static int overwrite = 0;	/* allow overwriting of input image file */
static int rot = 0;
static int mirror = 0;
static int bitpix = 0;
static int fitsout = 0; /* Output FITS file from IRAF input if 1 */

extern int RotFITS ();
extern int SetWCSFITS ();
extern int DelWCSFITS();
extern void settolerance ();
extern void setreflim ();
extern void setrot ();
extern void setnfit ();
extern void setsecpix ();
extern void setcenter ();
extern void setfk4 ();
extern void setminb ();
extern void setmaxr ();
extern void setstarsig ();
extern void setclass();
extern void setplate();
extern void setrefcat();
extern void setbmin();
extern void setfrac();
extern void setwcstype();

main (ac, av)
int ac;
char **av;
{
    char *str;
    double bmin;
    char rastr[16];
    char decstr[16];

    /* Decode arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
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

    	case 'f':	/* Write FITS file */
	    fitsout = 1;
    	    break;

	case 'g':	/* Guide Star object class */
    	    if (ac < 2)
    		usage();
    	    setclass ((int) atof (*++av));
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
    	    setreflim (atof (*++av));
    	    ac--;
    	    break;
    	case 'n':	/* Number of parameters to fit */
    	    if (ac < 2)
    		usage();
    	    setnfit ((int) atof (*++av));
    	    ac--;
    	    break;
    	case 'o':	/* allow overwriting of existing image file */
    	    overwrite++;
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
    	case 'w':	/* update the fields */
    	    writeheader++;
    	    break;
    	default:
    	    usage();
    	    break;
    	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage();

    if (!writeheader && !verbose) {
	fprintf (stderr, "IMWCS: Must have either w or v argument.\n");
	usage();
	}

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
    fprintf (stderr,"Set WCS in FITS and IRAF image files (after UIowa SETWCS)\n");
    fprintf(stderr,"Usage: [-vowdfl] [-m mag] [-n frac] [-s mode] [-g class] [-i peak] [-c catalog]\n");
    fprintf(stderr,"       [-p scale] [-b ra dec] [-j ra dec] [-r deg] [-t tol] FITS or IRAF file(s)\n");
    fprintf(stderr,"  -a: initial rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b: initial center in B1950 (FK4) RA and Dec\n");
    fprintf(stderr,"  -c: reference catalog (gsc, ujc, or tab table file\n");
    fprintf(stderr,"  -e: WCS type (TAN default)\n");
    fprintf(stderr,"  -f: write FITS output no matter what input\n");
    fprintf(stderr,"  -g: Guide Star Catalog class (-1=all,0,3 (default -1)\n");
    fprintf(stderr,"  -i: minimum peak value for star in image (<0=-sigma)\n");
    fprintf(stderr,"  -j: initial center in J2000 (FK5) RA and Dec\n");
    fprintf(stderr,"  -l: reflect left<->right before rotating and fitting\n");
    fprintf(stderr,"  -m: initial GSC or UJC limiting magnitude (default 17)\n");
    fprintf(stderr,"  -n: number of parameters to fit (2-5)\n");
    fprintf(stderr,"  -o: allow overwriting of input image, else write new one\n");
    fprintf(stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -r: rotation angle in degrees before fitting (default 0)\n");
    fprintf(stderr,"  -s: use this fraction extra stars (default 1.0)\n");
    fprintf(stderr,"  -t: offset tolerance in pixels (default 20)\n");
    fprintf(stderr,"  -u: UJ catalog single plate number to accept\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -w: write header (default is read-only)\n");
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

    /* Open IRAF image if .imh extension is present */
    if (strsrch (name,".imh") != NULL) {
	iraffile = 1;
	irafheader = irafrhead (name, &lhead);
	if (irafheader) {
	    header = iraf2fits (name, irafheader, lhead, &nbhead);
	    if (header == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",name);
		free (irafheader);
		return;
		}
	    image = irafrimage (header);
	    if (image == NULL) {
		hgets (header,"PIXFILE", 64, &pixname);
		fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
		free (irafheader);
		free (header);
		return;
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
	fitsout = 1;
	header = fitsrhead (name, &lhead, &nbhead);
	if (header) {
	    image = fitsrimage (name, nbhead, header);
	    if (image == NULL) {
		fprintf (stderr, "Cannot read FITS image %s\n", name);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", name);
	    return;
	    }
	}

    if (verbose) {
	fprintf (stderr,"Set World Coordinate System in ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}

    /* Print existing WCS headers and check for permission to overwrite */
    if (PrintWCS (header) == 0)
	(void) DelWCSFITS (header, verbose);

    /* Rotate and/or reflect image */
    if (rot != 0 || mirror) {
	if (RotFITS (name, header, &image, rot, mirror, bitpix, verbose)) {
	    fprintf (stderr,"Image %s could not be rotated\n");
	    return;
	    }
	if (!overwrite)
	    newimage = 1;
	}
    else if (overwrite)
	newimage = 0;
    else
	newimage = 1;

    /* Make up name for new FITS or IRAF output file */
    if (newimage) {

    /* Remove directory path and extension from file name */
	ext = strrchr (name, '.');
	fname = strrchr (name, '/');
	if (fname)
	    fname = fname + 1;
	else
	    fname = name;
	lname = strlen (fname);
	if (ext) {
	    lext = strlen (ext);
	    strncpy (newname, fname, lname - lext);
	    *(newname + lname - lext) = 0;
	    }
	else
	    strcpy (newname, fname);

    /* Add rotation and reflection to image name */
	if (rot < 10 && rot > -1)
	    sprintf (temp,"%1d",rot);
	else if (rot < 100 && rot > -10)
	    sprintf (temp,"%2d",rot);
	else if (rot < 1000 && rot > -100)
	    sprintf (temp,"%3d",rot);
	else
	    sprintf (temp,"%4d",rot);
	if (rot != 0)
	    strcat (newname, temp);
	if (mirror)
	    strcat (newname, "m");

    /* Add file extension preceded by a w */
	if (fitsout)
	    strcat (newname, "w.fit");
	else {
	    strcpy (pixname, "HDR$");
	    strcat (pixname, newname);
	    strcat (pixname, "w.pix");
	    hputs (header, "PIXFILE", pixname);
	    strcat (newname, "w.imh");
	    }
	}
    else
	strcpy (newname, name);

    if (SetWCSFITS (header, image, verbose)) {
	if (writeheader) {
	    if (verbose)
		(void) PrintWCS (header);	/* print new WCS */

	/* Log WCS program version in the image header */

	    if (fitsout) {
		if (fitswimage (newname, header, image) > 0 && verbose)
		    printf ("%s: rewritten successfully.\n", newname);
		else if (verbose)
		    printf ("%s could not be written.\n", newname);
		}
	    else if (newimage) {
		if (irafwimage (newname,lhead,irafheader,header,image) > 0 && verbose)
		    printf ("%s rewritten successfully.\n", newname);
		else if (verbose)
		    printf ("%s could not be written.\n", newname);
		}
	    else {
		if (irafwhead (newname,lhead,irafheader,header) > 0 && verbose)
		    printf ("%s rewritten successfully.\n", newname);
		else if (verbose)
		    printf ("%s could not be written.\n", newname);
		}
	    }
	else if (verbose)
	    printf ("%s: file unchanged.\n", name);
	}
    else if (verbose)
	printf ("%s: file unchanged.\n", name);

    free (header);
    if (iraffile)
	free (irafheader);
    free (image);
    return;
}


/* check the WCS fields and print any that are found if verbose.
 * return 0 if all are found, else -1.
 */
static int
PrintWCS (header)
char *header;

{
    char str[80];
    double v;
    int n;

    n = 0;

    if (hgets (header,"CTYPE1",16,str)) {
	if (verbose) printf ("CTYPE1 = %s\n", str);
	n++;
	}
    if (hgetr8 (header, "CRVAL1", &v)) {
	if (verbose) printf ("CRVAL1 = %.8f\n", v);
	n++;
	}
    if (hgetr8 (header, "CDELT1", &v)) {
	if (verbose) printf ("CDELT1 = %.8f\n", v);
	n++;
	}
    if (hgetr8 (header, "CRPIX1", &v)) {
	if (verbose) printf ("CRPIX1 = %.8f\n", v);
	n++;
	}
    if (hgetr8 (header, "CROTA1", &v)) {
	if (verbose) printf ("CROTA1 = %.3f\n", v);
	n++;
	}

    if (hgets (header,"CTYPE2",16,str)) {
	if (verbose) printf ("CTYPE2 = %s\n", str);
	n++;
	}
    if (hgetr8 (header, "CRVAL2", &v)) {
	if (verbose) printf ("CRVAL2 = %.8f\n", v);
	n++;
	}
    if (hgetr8 (header, "CDELT2", &v)) {
	if (verbose) printf ("CDELT2 = %.8f\n", v);
	n++;
	}
    if (hgetr8 (header, "CRPIX2", &v)) {
	if (verbose) printf ("CRPIX2 = %.8f\n", v);
	n++;
	}
    if (hgetr8 (header, "CROTA2", &v)) {
	if (verbose) printf ("CROTA2 = %.3f\n", v);
	n++;
	}
    if (hgets (header,"IMWCS",80,str)) {
	if (verbose) printf ("IMWCS = %s\n", str);
	n++;
	}

    return (n == 11 ? 0 : -1);
}


char *
{
}
/* Feb 16 1996	New program
 * Apr 15 1996	Move delWCSFITS to libwcs
 * Apr 24 1996	Add optional initial plate center on command line
 * May  2 1996	Add option to rotate and/or reflect before fitting
 * May 14 1996	Change GSCLIM to REFLIM
 * May 22 1996	Rearrange commands; set initial plate center explicitly
 * May 31 1996	Rename subroutines; drop WCS deletion as an option
 * Jun  4 1996	Allow writing of FITS file even if input is IRAF
 * Jun 14 1996	Write IMWCS record in header
 * Jun 28 1996	Add option to set number of parameters to fit
 * Jul  3 1996	Always write to new file unless -o
 * Jul 22 1996	Add option to change WCS projection
 * Aug  6 1996	Force number of decimal places in PrintWCS
 * Aug 26 1996	Change HGETC call to HGETS; fix irafrhead call
 * Aug 28 1996	Declare undeclared variables after lint
 * Aug 29 1996	Allow writing of new IRAF files
 * Sep  1 1996	Move parameter defaults to lwcs.h
 * Sep  3 1996	Fix star finding
 * Sep 17 1996	Fix bug in GSC reading
 * Oct 11 1996	Fix DelWCS declaration and do not free strings in PrintWCS
 */
