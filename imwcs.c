/* File imwcs.c
 * June 2, 1998
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
#include "libwcs/wcs.h"

static void usage();
static void FitWCS ();

static int verbose = 0;		/* verbose/debugging flag */
static int writeheader = 0;	/* write header fields; else read-only */
static int overwrite = 0;	/* allow overwriting of input image file */
static int rot = 0;
static int mirror = 0;
static int bitpix = 0;
static int fitsout = 0; /* Output FITS file from IRAF input if 1 */
static int imsearch = 1;	/* set to 0 if image catalog provided */
static int erasewcs = 0;	/* Set to 1 to erase initial image WCS */
char outname[128];		/* Name for output image */

extern int RotFITS ();
extern int SetWCSFITS ();
extern int DelWCSFITS();
extern int PrintWCS();
extern void settolerance();
extern void setreflim();
extern void setrot();
extern void setnfit();
extern void setsecpix();
extern void setcenter();
extern void setsys();
extern void setminb();
extern void setmaxcat();
extern void setstarsig();
extern void setclass();
extern void setplate();
extern void setrefcat();
extern void setimcat();
extern void setbmin();
extern void setfrac();
extern void setrefpix();
extern void setwcstype();
extern void setoldwcs();	/* AIPS classic WCS flag */
extern void setfitplate();
extern void setproj();
extern void setiterate();
extern void setrecenter();

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

    outname[0] = 0;

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
    	    if (ac < 3)
    		usage();
	    setsys (WCS_B1950);
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

	case 'e':	/* Erase WCS projection in image header */
	    erasewcs++;
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

    	case 'j':	/* center coordinates on command line in J2000 */
    	    if (ac < 3)
    		usage();
	    setsys (WCS_J2000);
	    strcpy (rastr, *++av);
	    ac--;
	    strcpy (decstr, *++av);
	    ac--;
	    setcenter (rastr, decstr);
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

    	case 'n':	/* Number of parameters to fit */
    	    if (ac < 2)
    		usage();
    	    setnfit ((int) atof (*++av));
    	    ac--;
    	    break;

    	case 'o':	/* Specifiy output image filename */
    	    if (ac < 2)
    		usage();
	    strcpy (outname, *++av);
	    ac--;
	    if (outname[0] == '-')
		overwrite++;
	    else
		overwrite = 0;
    	    writeheader++;
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
    	    break;

    	case 'q':	/* fit residuals */
    	    if (ac < 2)
    		usage();
	    if (strchr (*++av,'i') != NULL)
    		setiterate (1);
	    if (strchr (*av,'r') != NULL)
    		setrecenter (1);
	    if (strchr (*av,'p') != NULL)
    		setfitplate (6);
	    if (strchr (*av,'8') != NULL)
    		setfitplate (8);
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

    	case 't':	/* tolerance in pixels for star match */
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
	    ac--;
    	    break;

	case 'z':       /* Use AIPS classic WCS */
	    setoldwcs (1);
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
	if (!writeheader && !verbose) {
	    fprintf (stderr, "IMWCS: Must have either w or v argument.\n");
	    usage();
	    }
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMWCS: List file %s cannot be read\n",
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
    fprintf(stderr,"Usage: [-vwdfl] [-o filename] [-m mag] [-n frac] [-s mode] [-g class]\n");
    fprintf(stderr,"       [-h maxref] [-i peak] [-c catalog] [-p scale] [-b ra dec] [-j ra dec]\n");
    fprintf(stderr,"       [-r deg] [-t tol] [-x x y] [-y frac] FITS or IRAF file(s)\n");
    fprintf(stderr,"  -a: initial rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b: initial center in B1950 (FK4) RA and Dec\n");
    fprintf(stderr,"  -c: reference catalog (gsc, uac, usac, ujc, tab table file\n");
    fprintf(stderr,"  -d: Use following DAOFIND output catalog instead of search\n");
    fprintf(stderr,"  -e: Erase image WCS keywords\n");
    fprintf(stderr,"  -f: write FITS output no matter what input\n");
    fprintf(stderr,"  -g: Guide Star Catalog class (-1=all,0,3 (default -1)\n");
    fprintf(stderr,"  -h: maximum number of reference stars to use (10-200, default 25\n");
    fprintf(stderr,"  -i: minimum peak value for star in image (<0=-sigma)\n");
    fprintf(stderr,"  -j: initial center in J2000 (FK5) RA and Dec\n");
    fprintf(stderr,"  -l: reflect left<->right before rotating and fitting\n");
    fprintf(stderr,"  -m: reference catalog magnitude limit(s) (default none)\n");
    fprintf(stderr,"  -n: list of parameters to fit (12345678; negate for refinement)\n");
    fprintf(stderr,"  -o: name for output image, - to overwrite\n");
    fprintf(stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -q: <i>terate, <r>ecenter, <p>olynomial fit\n");
    fprintf(stderr,"  -r: rotation angle in degrees before fitting (default 0)\n");
    fprintf(stderr,"  -s: use this fraction extra stars (default 1.0)\n");
    fprintf(stderr,"  -t: offset tolerance in pixels (default 20)\n");
    fprintf(stderr,"  -u: USNO catalog single plate number to accept\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -w: write header (default is read-only)\n");
    fprintf(stderr,"  -x: X and Y coordinates of reference pixel (default is center)\n");
    fprintf(stderr,"  -y: add this fraction to image dimensions for search (default is 0)\n");
    fprintf(stderr,"  -z: use AIPS classic projections instead of WCSLIB\n");
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
    char newname[128];		/* Name for revised image */
    char pixname[128];		/* Pixel file name for revised image */
    char temp[16];
    char *ext;
    char *fname;
    int lext, lname;
    int newimage;
    char *imext, *imext1;

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
	    if (imsearch || writeheader || rot || mirror) {
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
	fitsout = 1;
	if ((header = fitsrhead (name, &lhead, &nbhead)) != NULL) {
	    if (imsearch || writeheader || rot || mirror) {
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
	fprintf (stderr,"Set World Coordinate System in ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}

    /* Print existing WCS keywords and optionally erase them */
    if (PrintWCS (header, verbose) == 0) {
	if (erasewcs)
	    (void) DelWCSFITS (header, verbose);
	}

    /* Rotate and/or reflect image */
    if ((imsearch || writeheader) && (rot != 0 || mirror)) {
	if (RotFITS (name, header, &image, rot, mirror, bitpix, verbose)) {
	    fprintf (stderr,"Image %s could not be rotated\n", name);
	    return;
	    }
	if (!overwrite)
	    newimage = 1;
	}

    /* Check for permission to overwrite */
    else if (overwrite)
	newimage = 0;
    else
	newimage = 1;

    /* Use output filename if it is set on the command line */
    if (outname[0] > 0)
	strcpy (newname, outname);

    /* Make up name for new FITS or IRAF output file */
    else if (newimage) {

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

    /* Add image extension number or name to output file name */
	imext = strchr (fname, ',');
	imext1 = NULL;
	if (imext == NULL) {
	    imext = strchr (fname, '[');
	    if (imext != NULL) {
		imext1 = strchr (fname, ']');
		*imext1 = (char) 0;
		}
	    }
	if (imext != NULL) {
	    strcat (newname, "_");
	    strcat (newname, imext+1);
	    }

    /* Add rotation and reflection to image name */
	if (mirror)
	    strcat (newname, "m");
	else if (rot != 0)
	    strcat (newname, "r");
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

    /* Add file extension preceded by a w */
	if (fitsout)
	    strcat (newname, "w.fits");
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

    if (SetWCSFITS (name, header, image, verbose)) {
	if (writeheader) {
	    if (verbose)
		(void) PrintWCS (header, verbose);	/* print new WCS */

	/* Log WCS program version in the image header */

	    if (fitsout) {
		if (fitswimage (newname, header, image) > 0 && verbose) {
		    if (overwrite)
			printf ("%s: rewritten successfully.\n", newname);
		    else
			printf ("%s: written successfully.\n", newname);
		    }
		else if (verbose)
		    printf ("%s could not be written.\n", newname);
		}
	    else if (newimage) {
		if (irafwimage (newname,lhead,irafheader,header,image) > 0 && verbose) {
		    if (overwrite)
			printf ("%s: rewritten successfully.\n", newname);
		    else
			printf ("%s: written successfully.\n", newname);
		    }
		else if (verbose)
		    printf ("%s could not be written.\n", newname);
		}
	    else {
		if (irafwhead (newname,lhead,irafheader,header) > 0 && verbose) {
		    if (overwrite)
			printf ("%s: rewritten successfully.\n", newname);
		    else
			printf ("%s: written successfully.\n", newname);
		    }
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
    if (image != NULL)
	free (image);
    return;
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
 * Oct 17 1996	Fix bugs which Sun C ignored
 * Nov 19 1996	Revised search subroutines, USNO A catalog added
 * Dec 10 1996	Revised WCS initialization
 * Dec 10 1996	Add option to get image stars from DAOFIND output list
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 * Mar 20 1997	Fix bug in GetFITSWCS which affected odd equinoxes
 * Apr 25 1997	Fix bug in uacread
 * May 28 1997  Add option to read a list of filenames from a file
 * Jul 12 1997	Add option to center reference pixel ccords on the command line
 * Aug 20 1997	Add option to set maximum number of reference stars to try to match
 * Sep  3 1996	Add option to change dimensions of search by fraction
 * Sep  9 1997	Add option to turn on residual refinement by negating nfit
 * Sep  9 1997	Do not read image unless it is needed
 * Oct 31 1997	Specify parameters to fit with numeric string
 * Nov  6 1997	Move PrintWCS to library
 * Nov  7 1997	Specify output image filename
 * Nov 14 1997	Change image increase from multiple to fraction added
 * Nov 17 1997	Add optional second magnitude limit
 * Dec  8 1997	Fixed bug in setting nominal WCS
 * Dec 15 1997	Add capability of reading and writing IRAF 2.11 images
 *
 * Jan 27 1998  Implement Mark Calabretta's WCSLIB
 * Jan 29 1998  Add -z for AIPS classic WCS projections
 * Feb 18 1998	Version 2.0: Calabretta WCS
 * Feb 20 1998	Add -q to fit residuals
 * Mar  1 1998	Add optional second axis plate scale argument to -p
 * Mar  3 1998	Add option to use first WCS result to try again
 * Mar  6 1998	Add option to recenter on second pass
 * Mar  6 1998	Change default FITS extension from .fit to .fits
 * Mar  6 1998	Add option to set projection type
 * Mar 27 1998	Version 2.2: Drop residual fitting; add polynomial fit
 * Apr 14 1998	New coordinate conversion software; polynomial debugged
 * Apr 24 1998	change coordinate setting to setsys() from setfk4()
 * Apr 27 1998	Add image extension name/number to output file name
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * May  4 1998	Make erasure of original image WCS optional
 * Jun  2 1998	Fix bugs in tabread() and hput()
 */
