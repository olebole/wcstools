/* File imstar.c
 * December 11, 1996
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
#include "libwcs/fitshead.h"
#include "libwcs/wcs.h"

static void usage();
static int verbose = 0;		/* verbose flag */
static int debug = 0;		/* debugging flag */

static void ListStars ();
extern void RASortStars ();
extern void FluxSortStars ();
extern void setstarsig ();
extern void setbmin ();
extern void setmaxrad ();
extern void setborder ();
extern void setimcat();
extern struct WorldCoor *GetFITSWCS();

static double magoff = 0.0;
static int rasort = 0;
static int printhead = 0;
static int tabout = 0;
static int nstar = 0;
static double cra0 = 0.0;
static double cdec0 = 0.0;
static double eqout = 0.0;

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    double bmin;
    char rastr[32], decstr[32];
    int maxrad;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {
	case 'v':	/* more verbosity */
	    verbose++;
	    break;

	case 'a':       /* Initial rotation angle in degrees */
	    if (ac < 2)
		usage();
	    setrot (atof (*++av));
	    ac--;
	    break;

	case 'b':	/* ouput FK4 (B1950) coordinates */
	    eqout = 1950.0;
	    break;

	case 'c':	/* Set center RA and Dec */
	    if (ac < 3)
		usage (progname);
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
	    ac--;
	    break;

	case 'e':	/* Number of pixels to ignore around image edge */
	    if (ac < 2)
		usage (progname);
	    setborder (atof (*++av));
	    ac--;
	    break;

	case 'h':	/* ouput descriptive header */
	    printhead++;
	    break;

	case 'i':	/* Image star minimum peak value (or minimum sigma */
	    if (ac < 2)
		usage (progname);
	    bmin = atof (*++av);
	    if (bmin < 0)
		setstarsig (-bmin);
	    else
		setbmin (bmin);
	    ac--;
	    break;

	case 'j':	/* ouput FK5 (J2000) coordinates */
	    eqout = 2000.0;
	    break;

	case 'k':	/* Print each star as it is found */
	    debug++;
	    break;

	case 'm':	/* Magnitude offset */
	    if (ac < 2)
		usage (progname);
	    magoff = atof (*++av);
	    ac--;
	    break;

	case 'n':	/* Number of brightest stars to read */
	    if (ac < 2)
		usage (progname);
	    nstar = atoi (*++av);
	    ac--;
	    break;

	case 'r':	/* Maximum acceptable radius for a star */
	    if (ac < 2)
		usage (progname);
	    maxrad = (int) atof (*++av);
	    setmaxrad (maxrad);
	    ac--;
	    break;

	case 's':	/* sort by RA */
	    rasort = 1;
	    break;

	case 't':	/* tab table to stdout */
	    tabout = 1;
	    break;

	default:
	    usage (progname);
	    break;
	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage (progname);

    while (ac-- > 0) {
	char *fn = *av++;
	ListStars (fn);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"Find stars in FITS and IRAF image files\n");
    fprintf(stderr,"%s: usage: [-vbsjt] [-m mag_off] [-n num] [-c ra dec]file.fts ...\n",
	    progname);
    fprintf(stderr,"  -a: initial rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b: Output B1950 (FK4) coordinates \n");
    fprintf(stderr,"  -c: Use following RA and Dec as center \n");
    fprintf(stderr,"  -d: Use following DAOFIND output catalog instead of search \n");
    fprintf(stderr,"  -e: Number of pixels to ignore around image edge \n");
    fprintf(stderr,"  -h: Print heading, else do not \n");
    fprintf(stderr,"  -i: Minimum peak value for star in image (<0=-sigma)\n");
    fprintf(stderr,"  -j: Output J2000 (FK5) coordinates \n");
    fprintf(stderr,"  -k: Print each star as it is found for debugging \n");
    fprintf(stderr,"  -m: Magnitude offset\n");
    fprintf(stderr,"  -n: Number of brightest stars to print \n");
    fprintf(stderr,"  -r: Maximum radius for star in pixels \n");
    fprintf(stderr,"  -s: Sort by RA instead of flux \n");
    fprintf(stderr,"  -t: Tab table to standard output as well as file\n");
    fprintf(stderr,"  -v: Verbose\n");
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
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    int *irafheader;		/* IRAF image header */
    double *sx=0, *sy=0;	/* image stars, pixels */
    double *sb=0;		/* image star brightesses */
    double *sra=0, *sdec=0;	/* image star RA and Dec */
    int ns;			/* n image stars */
    double *smag;		/* image star magnitudes */
    int *sp;			/* peak flux in counts */
    double ra, dec;
    double cra,cdec,dra,ddec,secpix;
    int wp, hp;
    char rastr[16], decstr[16];
    int i;
    char headline[160];
    char pixname[128];
    char outfile[64];
    FILE *fd;
    struct WorldCoor *wcs;	/* World coordinate system structure */

    /* Open IRAF header if .imh extension is present */
    if (strsrch (filename,".imh") != NULL) {
	irafheader = irafrhead (filename, &lhead);
	if (irafheader) {
	    header = iraf2fits (filename, irafheader, lhead, &nbhead);
	    image = irafrimage (header);
	    if (image == NULL) {
		hgets (header,"PIXFILE", 64, pixname);
		fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
		free (irafheader);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", filename);
	    free (header);
	    return;
	    }
	}

    /* Read FITS image header if .imh extension is not present */
    else {
	header = fitsrhead (filename, &lhead, &nbhead);
	if (header) {
	    image = fitsrimage (filename, nbhead, header);
	    if (image == NULL) {
		fprintf (stderr, "Cannot read FITS image %s\n", filename);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    return;
	    }
	}
    if (verbose && printhead)

/* Find the stars in an image and use the world coordinate system
 * information in the header to produce a plate catalog with right
 * ascension, declination, and a plate magnitude
 */

    wcs = GetFITSWCS (header,verbose, &cra, &cdec, &dra, &ddec, &secpix,
	  &wp, &hp, eqout);

    /* Discover star-like things in the image, in pixels */
    ns = FindStars (header, image, &sx, &sy, &sb, &sp, debug);
    if (ns < 1) {
	fprintf (stderr,"ListStars: no stars found in image %s\n", filename);
	return;
	}

    /* Save star positions */
    if (nstar > 0)
	ns = nstar;

    /* If no magnitude offset, set brightest star to 0 magnitude */
    if (ns > 0 && magoff == 0.0) {
	FluxSortStars (sx, sy, sb, sp, ns);
	magoff = 2.5 * log10 (sb[0]);
	}

    /* Compute right ascension and declination for all stars to be listed */
    sra = (double *) malloc (ns * sizeof (double));
    sdec = (double *) malloc (ns * sizeof (double));
    smag = (double *) malloc (ns * sizeof (double));
    for (i = 0; i < ns; i++) {
	pix2wcs (wcs, sx[i], sy[i], &sra[i], &sdec[i]);
	smag[i] = -2.5 * log10 (sb[i]) + magoff;
	}

    /* Sort star-like objects in image by right ascension */
    if (rasort)
	RASortStars (0, sra, sdec, sx, sy, sb, sp, ns);
    sprintf (headline, "IMAGE	%s", filename);

    /* Open plate catalog file */
    if (strcmp (filename,"stdin")) {
	strcpy (outfile,filename);
	strcat (outfile,".stars");
	}
    else {
	strcpy (outfile,filename);
	(void) hgets (header,"OBJECT",64,outfile);
	strcat (outfile,".stars");
	}
    if (verbose)
	printf ("%s\n", outfile);
		
    fd = fopen (outfile, "w");
    if (fd == NULL) {
	fprintf (stderr, "IMSTAR:  cannot write file %s\n", outfile);
        return;
        }

    /* Write header */
    fprintf (fd,"%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (rasort) {
	fprintf (fd, "RASORT	T\n");
	if (tabout)
	    printf ("RASORT	T\n");
	}

    if (strcmp (wcs->sysout,"FK4") == 0 || strcmp (wcs->sysout,"fk4") == 0)
	sprintf (headline, "EQUINOX	1950.0");
    else
	sprintf (headline, "EQUINOX	2000.0");
    fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (tabout)

    sprintf (headline,"ID 	RA      	DEC     	MAG   	X    	Y    	COUNTS   	PEAK");
    fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    sprintf (headline,"---	------------	------------	------	-----	-----	--------	------");
    fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    for (i = 0; i < ns; i++) {
	pix2wcs (wcs, sx[i], sy[i], &ra, &dec);
	ra2str (rastr, ra, 3);
	dec2str (decstr, dec, 2);
	sprintf (headline, "%d	%s	%s	%.2f	%.2f	%.2f	%.2f	%d",
		 i+1, rastr,decstr, smag[i], sx[i], sy[i], sb[i], sp[i]);
	fprintf (fd, "%s\n", headline);
	if (tabout)
	    printf ("%s\n", headline);
	else if (verbose)
	    printf ("%3d %s %s %6.2f %7.2f %7.2f %8.1f %d\n",
		    i+1, rastr, decstr, smag[i],sx[i],sy[i],sb[i], sp[i]);
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
 */
