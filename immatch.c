/* File immatch.c
 * July 7, 1999
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

#include "fitsfile.h"
#include "wcs.h"

static void usage();
static void MatchCat();

static int verbose = 0;		/* verbose/debugging flag */
static int overwrite = 0;	/* allow overwriting of input image file */
static int rot = 0;
static int mirror = 0;
static int bitpix = 0;
static int imsearch = 1;	/* set to 0 if image catalog provided */
static char *refcatname;	/* Name of reference catalog to match */
static int version = 0;		/* If 1, print only program name and version */

extern char *RotFITS();
extern int SetWCSFITS();
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
extern void setimcat();
extern void setbmin();
extern void setfrac();
extern void setrefpix();
extern void setwcsproj();
extern void setfitwcs();

main (ac, av)
int ac;
char **av;
{
    char *str;
    double bmin, maglim1, maglim2, arot, drot;
    char rastr[16];
    char decstr[16];
    int readlist = 0;
    char *lastchar;
    char filename[128];
    FILE *flist;
    char *listfile;
    double x, y;
    int i;
    char *refcatn;
    char progpath[128];
    char *progname;		/* Name of program as executed */

    setfitwcs (0);

    /* Check name used to execute programe and set catalog name accordingly */
    strcpy (progpath, av[0]);
    progname = progpath;
    for (i = strlen (progpath); i > -1; i--) {
	if (progpath[i] > 63 && progpath[i] < 90)
	    progpath[i] = progpath[i] + 32;
	if (progpath[i] == '/') {
	    progname = progpath + i + 1;
	    break;
	    }
	}
    if (strsrch (progname,"gsc") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "gsc");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"uac") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "uac");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"ua1") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "ua1");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"ua2") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "ua2");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"usac") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "usac");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"usa1") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "usa1");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"usa2") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "usa2");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"ujc") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "ujc");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"sao") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "sao");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"ppm") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "ppm");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"ira") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "iras");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"tyc") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "tycho");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"hip") != NULL) {
	refcatn = (char *) calloc (1,16);
	strcpy (refcatn, "hipparcos");
	refcatname = refcatn;
	}
    else if (strsrch (progname,"act") != NULL) {
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "act");
	refcatname = refcatn;
	}

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage (progname);
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage (progname);
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
    		usage (progname);
	    drot = atof (*++av);
	    arot = fabs (drot);
	    if (arot != 90.0 && arot != 180.0 && arot != 270.0) {
		setrot (rot);
		rot = 0;
		}
	    else
		rot = atoi (*av);
    	    ac--;
    	    break;

    	case 'b':	/* initial coordinates on command line in B1950 */
    	    if (ac < 3)
    		usage (progname);
	    setsys (WCS_B1950);
	    strcpy (rastr, *++av);
	    ac--;
	    strcpy (decstr, *++av);
	    ac--;
	    setcenter (rastr, decstr);
    	    break;

	case 'c':       /* Set reference catalog */
	    if (ac < 2)
		usage (progname);
	    refcatname = *++av;
	    ac--;
	    break;

	case 'd':	/* Read image star positions from DAOFIND file */
	    if (ac < 2)
		usage (progname);
	    setimcat (*++av);
	    imsearch = 0;
	    ac--;
	    break;

	case 'e':	/* Set WCS projection
	    if (ac < 2)
		usage (progname);
	    setwcsproj (*++av);
	    ac--;
	    break; */

	case 'f':	/* Set IRAF output format */
	    setirafout();
	    break;

    	case 'j':	/* center coordinates on command line in J2000 */
    	    if (ac < 3)
    		usage (progname);
	    setsys (WCS_J2000);
	    strcpy (rastr, *++av);
	    ac--;
	    strcpy (decstr, *++av);
	    ac--;
	    setcenter (rastr, decstr);
    	    break;

	case 'g':	/* Guide Star object class */
    	    if (ac < 2)
    		usage (progname);
    	    setclass ((int) atof (*++av));
    	    ac--;
    	    break;

	case 'h':	/* Maximum number of reference stars */
    	    if (ac < 2)
    		usage (progname);
    	    setmaxcat ((int) atof (*++av));
    	    ac--;
    	    break;

	case 'i':       /* Image star minimum peak value */
	    if (ac < 2)
		usage (progname);
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
    		usage (progname);
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
    		usage (progname);
    	    setsecpix (atof (*++av));
    	    ac--;
    	    break;

    	case 'r':	/* Angle in degrees to rotate before fitting */
    	    if (ac < 2)
    		usage (progname);
    	    rot = (int) atof (*++av);
    	    ac--;
    	    break;

    	case 's':   /* this fraction more image stars than GSC or vice versa */
    	    if (ac < 2)
    		usage (progname);
    	    setfrac (atof (*++av));
    	    ac--;
    	    break;

    	case 't':	/* +/- this many pixels is a hit */
    	    if (ac < 2)
    		usage (progname);
    	    settolerance (atof (*++av));
    	    ac--;
    	    break;

	case 'u':	/* UJ Catalog plate number */
    	    if (ac < 2)
    		usage (progname);
    	    setplate ((int) atof (*++av));
    	    ac--;
    	    break;

    	case 'v':	/* more verbosity */
    	    verbose++;
    	    break;

	case 'x':	/* X and Y coordinates of reference pixel */
	    if (ac < 3)
		usage (progname);
	    x = atof (*++av);
	    ac--;
	    y = atof (*++av);
	    ac--;
    	    setrefpix (x, y);
    	    break;

	case 'y':	/* Multiply dimensions of image by fraction */
	    if (ac < 2)
		usage (progname);
	    setimfrac (atof (*++av));
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
    	    usage (progname);
    	    break;
    	}
    }

    /* Find number of images to search and leave listfile open for reading */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMMATCH: List file %s cannot be read\n",
		     listfile);
	    usage (progname);
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
	    MatchCat (progname, filename);
	    if (verbose)
		printf ("\n");
	    }
	fclose (flist);
	}

    /* If no arguments left, print usage */
    if (ac == 0)
	usage (progname);

    while (ac-- > 0) {
    	char *fn = *av++;
	if (fn == NULL)
	    break;
    	if (verbose)
    	    printf ("%s:\n", fn);
    	MatchCat (progname, fn);
    	if (verbose)
    	    printf ("\n");
	}

    return (0);
}

static void
usage (progname)

char	*progname;		/* Name of program being executed */

{
    if (version)
	exit (-1);
    if (strsrch (progname,"gsc") != NULL)
	fprintf (stderr,"Match HST Guide Star Catalog to image stars from WCS in image file\n");
    else if (strsrch (progname,"ujc") != NULL)
	fprintf (stderr,"Match USNO J Catalog to image stars from WCS in image file\n");
    else if (strsrch (progname,"uac") != NULL)
	fprintf (stderr,"Match USNO A Catalog to image stars from WCS in image file\n");
    else if (strsrch (progname,"ua1") != NULL)
	fprintf (stderr,"Match USNO A1.0 Catalog to image stars from WCS in image file\n");
    else if (strsrch (progname,"ua2") != NULL)
	fprintf (stderr,"Match USNO A2.0 Catalog to image stars from WCS in image file\n");
    else if (strsrch (progname,"usac") != NULL)
	fprintf (stderr,"Match USNO SA Catalog to image stars from WCS in image file\n");
    else if (strsrch (progname,"usa1") != NULL)
	fprintf (stderr,"Match USNO SA1.0 Catalog to image stars from WCS in image file\n");
    else if (strsrch (progname,"usa2") != NULL)
	fprintf (stderr,"Match USNO SA2.0 Catalog to image stars from WCS in image file\n");
    else if (strsrch (progname,"act") != NULL)
	fprintf (stderr,"Match USNO ACT Catalog to image stars from WCS in image file\n");
    else if (strsrch (progname,"iras") != NULL)
	fprintf (stderr,"Match IRAS Point Source Catalog to image stars from WCS in image file\n");
    else if (strsrch (progname,"sao") != NULL)
	fprintf (stderr,"Match SAO Catalog to image stars from WCS in image file\n");
    else if (strsrch (progname,"ppm") != NULL)
	fprintf (stderr,"Match PPM Catalog to image stars from WCS in image file\n");
    else if (strsrch (progname,"tycho") != NULL)
	fprintf (stderr,"Match Tycho Catalog to image stars from WCS in image file\n");
    else
	fprintf (stderr,"Match catalog to image stars from WCS in image file\n");
    fprintf(stderr,"Usage: [-vl] [-m mag] [-n frac] [-s mode] [-g class] [-h maxref] [-i peak]\n");
    fprintf(stderr,"       [-c catalog] [-p scale] [-b ra dec] [-j ra dec] [-r deg] [-t tol] [-x x y] [-y frac]\n");
    fprintf(stderr,"       FITS or IRAF file(s)\n");
    fprintf(stderr,"  -a: initial rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b: initial center in B1950 (FK4) RA and Dec\n");
    fprintf(stderr,"  -c: reference catalog (gsc, uac, ujc, tab table file\n");
    fprintf(stderr,"  -d: Use following DAOFIND output catalog instead of search\n");
    /* fprintf(stderr,"  -e: WCS type (TAN default)\n"); */
    fprintf(stderr,"  -f: Write output X Y RA Dec, instead of N RA Dec X Y\n");
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
    fprintf(stderr,"  -z: use AIPS classic projections instead of WCSLIB\n");
    exit (1);
}


static void
MatchCat (progname, name)

char	*progname;		/* Name of program being executed */
char	*name;			/* Name of FITS or IRAF image file */

{
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    int iraffile;		/* 1 if IRAF image */
    char *image;		/* Image */
    char *header;		/* FITS header */
    char *irafheader;		/* IRAF image header */
    char newname[64];		/* Name for revised image */
    char pixname[64];		/* Pixel file name for revised image */
    char temp[16];
    char *ext;
    char *fname;
    int lext, lname;
    char *newimage;

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

    /* Open FITS file */
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
	if ((newimage = RotFITS (name,header,image,rot,mirror,bitpix,verbose))
	    == NULL) {
	    fprintf (stderr,"Image %s could not be rotated\n", name);
	    free (header);
	    if (iraffile)
		free (irafheader);
	    if (image != NULL)
		free (image);
	    return;
	    }
	free (image);
	image = newimage;
	}

    (void) SetWCSFITS (name, header, image, refcatname, verbose);

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
 * Dec 15 1997	Add capability of reading and writing IRAF 2.11 images
 *
 * Jan 27 1998  Implement Mark Calabretta's WCSLIB
 * Jan 29 1998  Add -z for AIPS classic WCS projections
 * Feb 18 1998	Version 2.0: Full Calabretta implementation
 * Mar 27 1998	Version 2.1: Add IRAF TNX projection
 * Apr 13 1998	Version 2.2: Add polynomial plate fit
 * Apr 24 1998	change coordinate setting to setsys() from setfk4()
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jun  2 1998	Fix bugs in hput() and tabread()
 * Jun 11 1998	Change setwcstype() to setwcsproj() to avoid conflict
 * Jun 25 1998	Set WCS subroutine choice with SETDEFWCS()
 * Jul 24 1998	Make irafheader char instead of int
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Oct 14 1998	Use isiraf() to determine file type
 * Nov 13 1998	Pass reference catalog name to SetWCSFITS
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Jan 26 1999	Add option to format output for IRAF coord fitting task
 * Apr 13 1999	Fix progname to drop / when full pathname
 * Jun  8 1999	Return image pointer from RotFITS, not flag
 * Jun 10 1999	If -a argument is multiple of 90, rotate image
 * Jul  7 1999	Fix bug setting rotation
 */
