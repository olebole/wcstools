/* File immwcs.c
 * January 28, 2000
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
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"

static void usage();
static void FitWCS ();

static int verbose = 0;		/* verbose/debugging flag */
static int writeheader = 0;	/* write header fields; else read-only */
static int overwrite = 0;	/* allow overwriting of input image file */
static int rot = 0;		/* Set to counterclockwise angle in degrees to rotate */
static int mirror = 0;		/* Set to 1 to flip image right to left */
static int bitpix = 0;		/* Bits per pixel for output image, if different */
static int fitsout = 0;		/* Output FITS file from IRAF input if 1 */
static int imsearch = 1;	/* set to 0 if image catalog provided */
static int version = 0;		/* If 1, print only program name and version */
char outname[128];		/* Name for output image */

extern char *RotFITS ();
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
extern void setmatch();
extern void setstarsig();
extern void setclass();
extern void setplate();
extern void setimcat();
extern void setbmin();
extern void setrefpix();
extern void setwcsproj();
extern void setfitplate();
extern void setproj();
extern void setiterate();
extern void setrecenter();
extern void setsecpix2();

main (ac, av)
int ac;
char **av;
{
    char *str;
    double maglim1, maglim2;
    char rastr[16];
    char decstr[16];
    int readlist = 0;
    char *lastchar;
    char filename[128];
    FILE *flist;
    char *listfile;
    double x, y, arot, drot;

    outname[0] = 0;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

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
	    drot = atof (*++av);
	    arot = fabs (drot);
	    if (arot != 90.0 && arot != 180.0 && arot != 270.0) {
		setrot (drot);
		rot = 0;
		}
	    else
		rot = atoi (*av);
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

	case 'd':	/* Read image star positions from DAOFIND file */
	    if (ac < 2)
		usage();
	    setmatch (*++av);
	    imsearch = 0;
	    ac--;
	    break;

	case 'e':	/* Set WCS projection */
	    if (ac < 2)
		usage();
	    setproj (*++av);
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

    	case 'f':	/* Write FITS file */
	    fitsout = 1;
    	    break;

	case 'h':	/* Maximum number of reference stars */
    	    if (ac < 2)
    		usage();
    	    setmaxcat ((int) atof (*++av));
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

	case 'z':       /* Use AIPS classic WCS */
	    setdefwcs (WCS_ALT);
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
	    fprintf (stderr, "IMMWCS: Must have either w or v argument.\n");
	    usage();
	    }
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMMWCS: List file %s cannot be read\n",
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
	fprintf (stderr, "IMMWCS: Must have either w or v argument.\n");
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
    if (version)
	exit (-1);
    fprintf (stderr,"Set WCS in FITS and IRAF image files from star/image match file\n");
    fprintf(stderr,"Usage: [-vwdfl] [-o filename] [-m mag] [-n frac] [-s mode] [-g class]\n");
    fprintf(stderr,"       [-h maxref] [-i peak] [-c catalog] [-p scale] [-b ra dec] [-j ra dec]\n");
    fprintf(stderr,"       [-r deg] [-t tol] [-x x y] [-y frac] FITS or IRAF file(s)\n");
    fprintf(stderr,"  -a: initial rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b: initial center in B1950 (FK4) RA and Dec\n");
    fprintf(stderr,"  -d: Use following file of ra/dec/x/y matches\n");
    fprintf(stderr,"  -e: WCS projection type (TAN default)\n");
    fprintf(stderr,"  -f: write FITS output no matter what input\n");
    fprintf(stderr,"  -h: maximum number of reference stars to use (10-200, default all\n");
    fprintf(stderr,"  -j: initial center in J2000 (FK5) RA and Dec\n");
    fprintf(stderr,"  -l: reflect left<->right before rotating and fitting\n");
    fprintf(stderr,"  -n: list of parameters to fit (12345678)\n");
    fprintf(stderr,"  -o: name for output image, - to overwrite\n");
    fprintf(stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -q: <i>terate, <r>ecenter, <p>olynomial fit\n");
    fprintf(stderr,"  -r: rotation angle in degrees before fitting (default 0)\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -w: write header (default is read-only)\n");
    fprintf(stderr,"  -x: X and Y coordinates of reference pixel (default is center)\n");
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
    char *irafheader;		/* IRAF image header */
    char newname[128];		/* Name for revised image */
    char pixname[128];		/* Pixel file name for revised image */
    char temp[16];
    char *ext;
    char *fname;
    int lext, lname;
    char *newimage;
    int rename;
    char *imext, *imext1;

    image = NULL;

    /* Open IRAF image */
    if (isiraf (name)) {
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

    /* Open FITS file */
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

    /* Print existing WCS headers */
    if (PrintWCS (header, verbose) == 0)
	(void) DelWCSFITS (header, verbose);

    /* Rotate and/or reflect image */
    if ((imsearch || writeheader) && (rot != 0 || mirror)) {
	if ((newimage = RotFITS (name,header,image,rot,mirror,bitpix,verbose))
	    == NULL) {
	    fprintf (stderr,"Image %s could not be rotated\n", name);
	    if (iraffile)
		free (irafheader);
	    if (image != NULL)
		free (image);
	    free (header);
	    return;
	    }
	else {
	    if (image != NULL)
		free (image);
	    image = newimage;
	    }

	if (!overwrite)
	    rename = 1;
	}

    /* Check for permission to overwrite */
    else if (overwrite)
	rename = 0;
    else
	rename = 1;

    /* Use output filename if it is set on the command line */
    if (outname[0] > 0)
	strcpy (newname, outname);

    /* Make up name for new FITS or IRAF output file */
    else if (rename) {

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
	    else if (rename) {
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
/* Apr 30 1998	New program based on IMWCS
 * May 27 1998	Include fitsio.h instead of fitshead.h
   Jun  2 1998  Fix bugs in hput() and tabread()
 * Jun 11 1998	Change setwcstype() to setwcsproj() to avoid conflict
 * Jun 25 1998	Set WCS subroutine choice with SETDEFWCS()
 * Jul 24 1998	Make irafheader char instead of int
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Oct 14 1998	Use isiraf() to determine file type
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Jun 10 1999  If -a argument is multiple of 90, rotate image
 * Jul  7 1999	Fix bug setting rotation
 * Oct 22 1999	Drop unused variables after lint
 *
 * Jan 28 2000	Call setdefwcs() with WCS_ALT instead of 1
 */
