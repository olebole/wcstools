/* File imstar.c
 * September 20, 1999
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
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
#include "wcscat.h"

static int verbose = 0;		/* verbose flag */
static int debug = 0;		/* debugging flag */
static int version = 0;		/* If 1, print only program name and version */
static int rot = 0;		/* Angle to rotate image (multiple of 90 deg) */
static int mirror = 0;		/* If 1, flip image right-left before rotating*/

static void usage();
static void ListStars();
extern char *RotFITS();
extern void RASortStars();
extern void FluxSortStars();
extern void setstarsig();
extern void setbmin();
extern void setmaxrad();
extern void setborder();
extern void setimcat();
extern void setparm();
extern struct WorldCoor *GetFITSWCS();

static double magoff = 0.0;
static int rasort = 0;
static int printhead = 0;
static int tabout = 0;
static int tabfile = 1;
static int nstar = 0;
static double cra0 = 0.0;
static double cdec0 = 0.0;
static double eqout = 0.0;
static int sysout = -1;
static int daofile = 0;
static int ascfile = 0;
static int setuns = 0;		/* Change to unsigned integer flag */
static int imsearch = 1;	/* If 1, search for stars in image */
static int region_char;
static int region_radius;

main (ac, av)
int ac;
char **av;
{
    char *str;
    double bmin, arot, drot;
    char rastr[32], decstr[32];
    int readlist = 0;
    char *lastchar;
    char *cstr;
    char filename[128];
    FILE *flist;
    char *listfile;
    int maxrad;

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
   	str = *av; 
	if (*str == '@') {
	    readlist++;
	    listfile = ++str;
	    str = str + strlen (str) - 1;
	    av++;
	    ac--;
	    }
	else if (strchr (str, '='))
	    setparm (str);
	else if (*str == '-') {
	    char c;
	    while (c = *++str) {
		switch (c) {
		case 'v':	/* more verbosity */
		    verbose++;
		    break;
	
		case 'a':       /* Initial rotation angle in degrees */
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
	
		case 'b':	/* ouput FK4 (B1950) coordinates */
		    eqout = 1950.0;
		    sysout = WCS_B1950;
		    break;
	
		case 'c':	/* Set center RA and Dec */
		    if (ac < 3)
			usage();
		    strcpy (rastr, *++av);
		    ac--;
		    strcpy (decstr, *++av);
		    ac--;
		    setcenter (rastr, decstr);
		    break;
	
		case 'd':	/* Read image star positions from DAOFIND file */
		    if (ac < 2)
			usage();
		    setimcat (*++av);
		    imsearch = 0;
		    ac--;
		    break;
	
		case 'e':	/* Number of pixels to ignore around image edge */
		    if (ac < 2)
			usage();
		    setborder (atof (*++av));
		    ac--;
		    break;
	
		case 'f':	/* Write ASCII catalog format for SKYMAP */
		    ascfile = 1;
		    tabfile = 0;
		    break;
	
		case 'h':	/* ouput descriptive header */
		    printhead++;
		    break;
	
		case 'i':	/* Image star minimum peak value (or minimum sigma */
		    if (ac < 2)
			usage();
		    bmin = atof (*++av);
		    if (bmin < 0)
			setstarsig (-bmin);
		    else
			setbmin (bmin);
		    ac--;
		    break;
	
		case 'j':	/* ouput FK5 (J2000) coordinates */
		    eqout = 2000.0;
		    sysout = WCS_J2000;
		    break;
	
		case 'k':	/* Print each star as it is found */
		    debug++;
		    break;
	
    		case 'l':	/* Left-right reflection before rotating */
		    mirror = 1;
    		    break;
	
		case 'm':	/* Magnitude offset */
		    if (ac < 2)
			usage();
		    magoff = atof (*++av);
		    ac--;
		    break;
	
		case 'n':	/* Number of brightest stars to read */
		    if (ac < 2)
			usage();
		    nstar = atoi (*++av);
		    ac--;
		    break;
	
    		case 'p':	/* Plate scale in arcseconds per pixel */
    		    if (ac < 2)
    			usage();
    		    setsecpix (atof (*++av));
    		    ac--;
    		    break;
	
		case 'q':	/* Output region file shape for SAOimage */
    		    if (ac < 2)
    			usage();
		    cstr = *++av;
		    switch (cstr[0]){
			case 'c':
			    if (cstr[1] == 'i')
				region_char = WCS_CIRCLE;
			    else
				region_char = WCS_CROSS;
			    break;
			case 'd':
			    region_char = WCS_DIAMOND;
			    break;
			case 's':
			    region_char = WCS_SQUARE;
			    break;
			case 'x':
			    region_char = WCS_EX;
			    break;
			case 'v':
			    region_char = WCS_VAR;
			    break;
			case '+':
			    region_char = WCS_CROSS;
			    break;
			case 'o':
			default:
			    region_char = WCS_CIRCLE;
			}
		    if (region_radius == 0)
			region_radius = 10;
    		    ac--;
		    break;
	
		case 'r':	/* Maximum acceptable radius for a star */
		    if (ac < 2)
			usage();
		    maxrad = (int) atof (*++av);
		    region_radius = maxrad;
		    setmaxrad (maxrad);
		    ac--;
		    break;
	
		case 's':	/* sort by RA */
		    rasort = 1;
		    break;
	
		case 't':	/* tab table to stdout */
		    tabout = 1;
		    break;
	
		case 'u':	/* Set 16-bit int image file to unsigned */
		    setuns = 1;
		    break;
	
		case 'w':	/* Write DAOFIND-format output file */
		    daofile = 1;
		    tabfile = 0;
		    break;

		case 'x':	/* X and Y coordinates of reference pixel */
		    if (ac < 3)
			usage();
    		    setrefpix (atof (*++av), atof (*++av));
		    ac--;
		    ac--;
    		    break;

		case 'z':       /* Use AIPS classic WCS */
		    setdefwcs (1);
		    break;

		default:
		    usage();
		}
		}
	    }

	/* Otherwise assume that this argument is a FITS or IRAF file */
	else {
	    char *fn = str;
	    ListStars (fn);
	    if (verbose)
		printf ("\n");
	    }
	}

    /* Find number of images to search and leave listfile open for reading */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMSTAR: List file %s cannot be read\n",
		     listfile);
	    usage();
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
	    ListStars (filename);
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
    fprintf (stderr,"Find stars in FITS and IRAF image files\n");
    fprintf(stderr,"usage: imstar [-vbsjt] [-m mag_off] [-n num] [-c ra dec]file.fits ...\n");
    fprintf(stderr,"  -a: initial rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b: Output B1950 (FK4) coordinates \n");
    fprintf(stderr,"  -c: Use following RA and Dec as center \n");
    fprintf(stderr,"  -d: Use following DAOFIND output catalog instead of search \n");
    fprintf(stderr,"  -e: Number of pixels to ignore around image edge \n");
    fprintf(stderr,"  -f: Write simple ASCII catalog file, not tab table \n");
    fprintf(stderr,"  -h: Print heading, else do not \n");
    fprintf(stderr,"  -i: Minimum peak value for star in image (<0=-sigma)\n");
    fprintf(stderr,"  -j: Output J2000 (FK5) coordinates \n");
    fprintf(stderr,"  -k: Print each star as it is found for debugging \n");
    fprintf(stderr,"  -l: reflect left<->right before rotating and searching\n");
    fprintf(stderr,"  -m: Magnitude offset (set brightest to abs(offset) if < 0)\n");
    fprintf(stderr,"  -n: Number of brightest stars to print \n");
    fprintf(stderr,"  -p: Plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -q: Output region file shape for SAOimage (default o)\n");
    fprintf(stderr,"  -r: Maximum radius for star in pixels \n");
    fprintf(stderr,"  -s: Sort by RA instead of flux \n");
    fprintf(stderr,"  -t: Tab table format star list\n");
    fprintf(stderr,"  -u: Set BITPIX to -16 for unsigned integer\n");
    fprintf(stderr,"  -v: Verbose; print star list to stdout\n");
    fprintf(stderr,"  -w: write DAOFIND format star list to output file\n");
    fprintf(stderr,"  -x: X and Y coordinates of reference pixel (if not in header or center)\n");
    fprintf (stderr,"  -z: use AIPS classic projections instead of WCSLIB\n");
    fprintf(stderr,"  @listfile: file containing a list of filenames to search\n");
    exit (1);
}


extern int FindStars ();
struct WorldCoor *wcsinit();	
extern int pix2wcst();

static void
ListStars (filename)

char	*filename;	/* FITS or IRAF file filename */

{
    char *image;		/* FITS image */
    char *newimage;		/* Rotated FITS image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    char *irafheader;		/* IRAF image header */
    double *sx=0, *sy=0;	/* image stars, pixels */
    double *sb=0;		/* image star brightesses */
    double *sra=0, *sdec=0;	/* image star RA and Dec */
    int ns;			/* n image stars */
    double *smag;		/* image star magnitudes */
    int *sp;			/* peak flux in counts */
    double ra, dec;
    double cra,cdec,dra,ddec,secpix;
    int wp, hp;
    char rastr[32], decstr[32];
    int i, bitpix;
    char headline[160];
    char pixname[128];
    char outfile[64];
    char *ext;
    char temp[32];
    FILE *fd;
    struct WorldCoor *wcs;	/* World coordinate system structure */
    int iraffile = 0;

    /* Open IRAF header */
    if (isiraf (filename)) {
	if ((irafheader = irafrhead (filename, &lhead)) != NULL) {
	    if ((header = iraf2fits (filename,irafheader,lhead,&nbhead)) == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",filename);
		free (irafheader);
		return;
		}
	    if (imsearch) {
		if ((image = irafrimage (header)) == NULL) {
		    hgets (header,"PIXFILE", 64, pixname);
		    fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
		    free (irafheader);
		    free (header);
		    return;
		    }
		}
	    iraffile = 1;
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", filename);
	    free (header);
	    return;
	    }
	}

    /* Read FITS image header */
    else {
	if ((header = fitsrhead (filename, &lhead, &nbhead)) != NULL) {
	    if (imsearch) {
		if ((image = fitsrimage (filename, nbhead, header)) == NULL) {
		    fprintf (stderr, "Cannot read FITS image %s\n", filename);
		    free (header);
		    return;
		    }
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    return;
	    }
	}
    if (verbose && printhead)

    /* Set image to unsigned integer if 16-bit and flag set */
    if (setuns) {
	hgeti4 (header, "BITPIX", &bitpix);
	if (bitpix == 16)
	    hputi4 (header, "BITPIX", -16);
	}

    /* Rotate and/or reflect image */
    if (imsearch && (rot != 0 || mirror)) {
	if ((newimage=RotFITS(filename,header,image,rot,mirror,bitpix,verbose))
	    == NULL) {
	    fprintf (stderr,"Image %s could not be rotated\n", filename);
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
	}

/* Find the stars in an image and use the world coordinate system
 * information in the header to produce a plate catalog with right
 * ascension, declination, and a plate magnitude
 */

    wcs = GetFITSWCS (filename, header,verbose, &cra, &cdec, &dra, &ddec,
		      &secpix, &wp, &hp, &sysout, &eqout);
    if (wcs == NULL && !tabfile)
	daofile = 1;

    /* Discover star-like things in the image, in pixels */
    ns = FindStars (header, image, &sx, &sy, &sb, &sp, debug);
    if (ns < 1) {
	fprintf (stderr,"ListStars: no stars found in image %s\n", filename);
	return;
	}

    /* Save star positions */
    if (nstar > 0 && ns > nstar)
	ns = nstar;

    /* If no magnitude offset, set brightest star to 0 magnitude */
    /* If magnitude offset is given < 0, set brightest star to abs(magoff) */
    FluxSortStars (sx, sy, sb, sp, ns);
    if (ns > 0 && magoff <= 0.0) {
        magoff = 2.5 * log10 (sb[0]) - magoff;
        }

    /* Compute right ascension and declination for all stars to be listed */
    smag = (double *) malloc (ns * sizeof (double));
    sra = (double *) malloc (ns * sizeof (double));
    sdec = (double *) malloc (ns * sizeof (double));
    for (i = 0; i < ns; i++) {
	if (iswcs (wcs))
	    pix2wcs (wcs, sx[i], sy[i], &sra[i], &sdec[i]);
	else {
	    sra[i] = 0.0;
	    sdec[i] = 0.0;
	    }
	smag[i] = -2.5 * log10 (sb[i]) + magoff;
	}

    /* Sort star-like objects in image by right ascension */
    if (rasort && iswcs (wcs))
	RASortStars (0, sra, sdec, sx, sy, sb, 0, sp, ns);
    sprintf (headline, "IMAGE	%s", filename);

    /* Open plate catalog file */
    if (strcmp (filename,"stdin")) {
	if (strrchr (filename, '/'))
	    strcpy (outfile, strrchr (filename, '/')+1);
	else
	    strcpy (outfile,filename);
	if ((ext = strsrch (outfile, ".fit")) != NULL ||
	    (ext = strsrch (outfile, ".imh")) != NULL)
	    *ext = (char) 0;
	}
    else {
	strcpy (outfile,filename);
	(void) hgets (header,"OBJECT",64,outfile);
	}

    /* Add rotation and reflection to output file name */
    if (mirror)
	strcat (outfile, "m");
    else if (rot != 0)
	strcat (outfile, "r");
    if (rot != 0) {
	if (rot < 10 && rot > -1) {
	    sprintf (temp,"%1d",rot);
	    strcat (outfile, temp);
	    }
	else if (rot < 100 && rot > -10) {
	    sprintf (temp,"%2d",rot);
	    strcat (outfile, temp);
	    }
	else if (rot < 1000 && rot > -100) {
	    sprintf (temp,"%3d",rot);
	    strcat (outfile, temp);
	    }
	else {
	    sprintf (temp,"%4d",rot);
	    strcat (outfile, temp);
	    }
	}

    if (region_char)
	strcat (outfile, ".reg");
    else if (daofile || nowcs (wcs))
	strcat (outfile,".dao");
    else
	strcat (outfile,".imstar");
    if (verbose)
	printf ("%s\n", outfile);
		
    fd = fopen (outfile, "w");
    if (fd == NULL) {
	fprintf (stderr, "IMSTAR:  cannot write file %s\n", outfile);
        return;
        }

    /* Write file of positions for SAOimage regions */
    if (region_char) {
	int radius, ix, iy;
	char snum[32], rstr[16];
	fprintf (fd, "# stars in %s\n", filename);
	switch (region_char) {
	    case WCS_SQUARE:
		strcpy (rstr, "SQUARE");
		break;
	    case WCS_DIAMOND:
		strcpy (rstr, "DIAMOND");
		break;
	    case WCS_CROSS:
		strcpy (rstr, "CROSS");
		break;
	    case WCS_EX:
		strcpy (rstr, "EX");
		break;
	    case WCS_CIRCLE:
	    default:
		strcpy (rstr, "CIRCLE");
	    }
	radius = region_radius;
	for (i = 0; i < ns; i++) {
	    ix = (int)(sx[i] + 0.5);
	    iy = (int)(sy[i] + 0.5);
	    fprintf (fd, "%s(%d,%d,%d) # %s %d\n",
		     rstr, ix, iy, radius, filename, i);
	    }
	}
    else {

    /* Write header */
    if (tabfile)
	fprintf (fd,"%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (daofile)
	fprintf (fd,"#%s\n", headline);
    if (iswcs (wcs)) {
	if (rasort) {
	    if (tabfile)
		fprintf (fd, "RASORT	T\n");
	    if (tabout)
		printf ("RASORT	T\n");
	    }

	if (wcs->sysout == WCS_B1950)
	    sprintf (headline, "EQUINOX	1950.0");
	else
	    sprintf (headline, "EQUINOX	2000.0");
	if (ascfile) {
	    if (wcs->sysout == WCS_B1950)
		fprintf (fd, "%s.cat\n", filename);
	    else
		fprintf (fd, "%s.cat/j\n", filename);
	    }
	else if (tabfile)
	    fprintf (fd, "%s\n", headline);
	else if (daofile)
	    fprintf (fd, "#%s\n", headline);
	if (tabfile)
	    fprintf (fd, "EPOCH	%9.4f\n", wcs->epoch);
	if (tabout)
	    printf ("%s\n", headline);
	if (ascfile)
	else if (daofile)
	else if (tabfile)
	if (tabout)

	sprintf (headline,"ID 	RA      	DEC     	MAG   	X    	Y    	COUNTS   	PEAK");
	if (tabfile)
	    fprintf (fd, "%s\n", headline);
	if (tabout)
	    printf ("%s\n", headline);
	sprintf (headline,"---	------------	------------	------	-----	-----	--------	------");
	if (tabfile)
	    fprintf (fd, "%s\n", headline);
	if (tabout)
	    printf ("%s\n", headline);
	}

    for (i = 0; i < ns; i++) {
	ra2str (rastr, 32, sra[i], 3);
	dec2str (decstr, 32, sdec[i], 2);
	sprintf (headline, "%d	%s	%s	%.2f	%.2f	%.2f	%.2f	%d",
		     i+1, rastr,decstr, smag[i], sx[i], sy[i], sb[i], sp[i]);
	if (tabout)
	    printf ("%s\n", headline);
	if (tabfile)
	    fprintf (fd, "%s\n", headline);
	if (daofile) {
	    sprintf (headline, "%7.2f %7.2f %6.2f  %d",
		    sx[i],sy[i],smag[i],sp[i]);
	    if (iswcs (wcs))
		sprintf (headline, "%s %s %s", headline, rastr, decstr);
	    fprintf (fd, "%s\n", headline);
	    }
	sprintf (headline, "%3d %s %s %.2f", i+1,rastr,decstr,smag[i]);
	sprintf (headline, "%s  %.2f %.2f %.2f %d",
	    headline, sx[i],sy[i],sb[i], sp[i]);
	if (ascfile)
	    fprintf (fd, "%s\n", headline);
	if (verbose)
	    printf ("%s\n", headline);
	}
	}

    fclose (fd);
    if (sx) free ((char *)sx);
    if (sy) free ((char *)sy);
    if (sb) free ((char *)sb);
    if (sra) free ((char *)sra);
    if (sdec) free ((char *)sdec);
    if (smag) free ((char *)smag);
    free ((char *)wcs);
    free (header);
    if (imsearch)
	free (image);
    return;
}

/* Feb 29 1996	New program
 * Apr 30 1996	Add FOCAS-style catalog matching
 * May  1 1996	Add initial image center from command line
 * May  2 1996	Set up four star matching modes
 * May 14 1996	Pass verbose flag; allow other reference catalogs
 * May 21 1996	Sort by right ascension; allow output in FK4 or FK5
 * May 29 1996	Add optional new image center coordinates
 * Jun 10 1996	Drop 3 arguments flux sorting subroutine
 * Jul 16 1996	Update input code
 * Aug 26 1996	Change HGETC call to HGETS
 * Aug 27 1996	Remove unused variables after lint
 * Aug 30 1996	Allow border to be set
 * Sep  1 1996	Move parameter defaults to lwcs.h
 * Oct 17 1996	Drop unused variables
 * Dec 10 1996	Improve hot pixel rejection
 * Dec 11 1996	Allow reading from DAOFIND file instead of searching image
 * Dec 11 1996	Add WCS default rotation and use getfitswcs
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 * Mar 18 1997	Skip WCS calls if no WCS
 * May 28 1997	Add option to read a list of filenames from a file
 * Jul 12 1997  Add option to center reference pixel coords on the command line
 * Nov  7 1997	Print file in tab, DAO, or ASCII format, just like STDOUT
 * Dec 16 1997	Support IRAF 2.11 image headers
 * Dec 16 1997	Fix spacing in non-tab-table output
 *
 * Jan 27 1998  Implement Mark Calabretta's WCSLIB
 * Jan 29 1998  Add -z for AIPS classic WCS projections
 * Feb 18 1998	Version 2.0: Full Calabretta WCS
 * Mar  2 1998	Fix RA sorting bug
 * Mar 27 1998	Version 2.1: Add IRAF TNX projection
 * Mar 27 1998	Version 2.2: Add polynomial plate fit
 * Apr 27 1998	Drop directory from output file name
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * May 28 1998	Include fitsio.h instead of fitshead.h
 * Jun  2 1998  Fix bug in hput()
 * Jun 15 1998	Default to tab table file; ASCII table verbosee
 * Jun 15 1998	Write DAO-format file if -w flag set; ASCII table verbose
 * Jun 17 1998	Add option to set 16-bit files to unsigned int BITPIX=-16
 * Jul 24 1998	Make irafheader char instead of int
 * Jul 27 1998	Fix bug in ra2str() and dec2str() arguments
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Sep 17 1998	Add coordinate system to GetFITSWCS() argument list
 * Sep 29 1998	Changesystem and equinox arguments to GetFITSWCS()
 * Oct 14 1998	Use isiraf() to determine file type
 * Oct 27 1998	Add option to write region file to plot results over image
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Apr  7 1999	Add filename argument to GetFITSWCS
 * May 25 1999	Add epoch to output tab table header
 * Jun  9 1999	Set brightest magnitude for any magnitude offset (J-B Marquette)
 * Jun 10 1999	Add option to rotation and reflect image before searching
 * Jun 10 1999	Drop .fits or .imh file extension from output file name
 * Jun 11 1999	Add parameter setting on command line
 * Jul  1 1999	Only free image if it was allocated
 * Jul  7 1999	Fix bug setting rotation
 * Jul  7 1999	Do not add 0 to file name if no rotation
 * Jul  7 1999	If -n argument more than found stars, list only number found
 * Sep 20 1999	Drop second call to pix2wcs
 */
