/* File imstar.c
 * February 26, 1996
 * By Doug Mink
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include "fitshead.h"
#include "wcs.h"

#define MAXHEADLEN 14400

static void usage();

static int verbose = 0;		/* verbose/debugging flag */
static char *RevMsg = "IMSTAR version 1.0, 21 February 1996";

static void listStars ();
static double magoff = 0.0;

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    int nstar = 0;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {
	case 'v':	/* more verbosity */
	    verbose++;
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
	if (verbose)
	    printf ("%s:\n", fn);
		
	listStars (fn, nstar);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"%s\n",RevMsg);
    fprintf (stderr,"Find stars in FITS and IRAF image files\n");
    fprintf (stderr, "By D. Mink, SAO\n");
    fprintf(stderr,"%s: usage: [-v] [-m mag_off] [-n num] file.fts ...\n",
	    progname);
    fprintf(stderr,"  -m: magnitude offset\n");
    fprintf(stderr,"  -n: number of brightest stars to print \n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


extern int findStars ();
extern void sortStars ();
struct WorldCoor *wcsinit();	
extern int pix2wcst();

static void
listStars (filename, nstars)

char	*filename;	/* FITS or IRAF file filename */
int	nstars;		/* Number of brightest stars to list */

{
    char *image;		/* FITS image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int iraffile;		/* 1 if IRAF image */
    int *irafheader;		/* IRAF image header */
    double *sx=0, *sy=0;	/* image stars, pixels */
    double *sb=0;		/* image star brightesses */
    int ns;			/* n image stars */
    double mag;
    char wcstring[64];
    int i;
    FILE *fd;
    struct WorldCoor *wcs;	/* World coordinate system structure */

    /* Allocate FITS header */
    header = malloc (MAXHEADLEN);
    lhead = MAXHEADLEN;

    /* Open IRAF image if .imh extension is present */
    if (strsrch (filename,".imh") != NULL) {
	iraffile = 1;
	irafheader = irafrhead (filename, lhead, header);
	if (irafheader == NULL) {
	    free (header);
	    return;
	    }
	image = irafrimage (filename, irafheader, header);
	if (image == NULL) {
	    free (header);
	    free (irafheader);
	    return;
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	image = fitsrimage (filename, lhead, header);
	if (image == NULL) {
	    free (header);
	    return;
	    }
	}
    if (verbose) {
	fprintf (stderr,"%s\n",RevMsg);
	fprintf (stderr,"Find stars in FITS or IRAF image file\n");
	fprintf (stderr, "By D. Mink, SAO\n");
	}

/* Find the stars in an image and use the world coordinate system
 * information in the header to produce a plate catalog with right
 * ascension, declination, and a plate magnitude
 */

    wcs = wcsinit (header);
    wcs->printsys = 0;

    /* Discover star-like things in the image, in pixels */
    ns = findStars (header, image, &sx, &sy, &sb);
    if (ns < 1) {
	printf ("imWCSstars: no stars found in image %s\n", filename);
	return;
	}

    /* Sort star-like objects in image by brightness */
    sortStars (sx, sy, sb, ns);

    /* Open plate catalog file */
    strcat (filename,".stars");
    fd = fopen (filename, "w");
    if (fd == NULL) {
	fprintf (stderr, "IMSTAR:  cannot write file %s %d\n", filename, fd);
        return;
        }

    /* Write header */
    fprintf (fd, "IMAGE	%s\n", filename);
    fprintf (fd, "EQUINOX	2000.0\n");
    fprintf (fd,"ID 	RA      	DEC     	MAG   	X    	Y    	V\n");
    fprintf (fd,"---	------------	------------	------	-----	-----	--------\n");


    /* Save star positions */
    if (nstars > 0)
	ns = nstars;
    for (i = 0; i < ns; i++) {
	mag = -2.5 * log10 (sb[i]) + magoff;
	if (verbose) {
	    wcs->tabsys = 0;
	    pix2wcst (wcs, sx[i], sy[i], wcstring, 64);
	    printf ("%3d %s %6.2f %6.1f %6.1f %8.1f\n",
		    i+1, wcstring, mag,sx[i],sy[i],sb[i]);
	    }
	wcs->tabsys = 1;
	pix2wcst (wcs, sx[i], sy[i], wcstring, 64);
	fprintf (fd, "%d	%s	%.2f	%.1f	%.1f	%.1f\n",
		 i+1, wcstring, mag,sx[i],sy[i],sb[i]);
	}

    fclose (fd);
    if (sx)
	free ((char *)sx);
    if (sy)
	free ((char *)sy);
    if (sb)
	free ((char *)sb);
    free ((char *)wcs);

    free (header);
    if (iraffile)
	free (irafheader);
    free (image);
    return;
}
