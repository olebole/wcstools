/* File setwcs.c
 * February 16, 1996
 * By Elwood Downey, revised by Doug Mink
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>

#include "fitshead.h"

#define MAXHEADLEN 14400

static void usage();
static void delPos ();
static void calPos ();
static void checkEquinox ();
static int checkWCSFITS ();

static int verbose = 0;		/* verbose/debugging flag */
static int writeheader = 0;	/* write header fields; else read-only */
static int overwrite = 0;	/* allow overwriting if fields already exist */
static int delete = 0;		/* delete the WCD header fields */
static char *RevMsg = "SETWCS version 1.0, 8 February 1996";

extern int setWCSFITS ();
extern void settolerance ();
extern void setgsclim ();
extern void setrot ();
extern void setmatch ();
extern void setsecpix ();

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    double rot, scale, gsclim, frac;
    int tolerance;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
    	switch (c) {
    	case 'o':	/* allow overwriting existing WCS headers */
    	    overwrite++;
    	    break;
    	case 'v':	/* more verbosity */
    	    verbose++;
    	    break;
    	case 'f':	/* image flipped around N-S axis */
	    setflip (-1);
    	    break;
    	case 'j':	/* Use Joe Mohr's triangle-matching code */
	    setmatch ();
    	    break;
    	case 'm':	/* Limiting GSC magnitude */
    	    if (ac < 2)
    		usage (progname);
    	    gsclim = atof (*++av);
    	    setgsclim (gsclim);
    	    ac--;
    	    break;
    	case 'n':	/* this many more image stars than GSC or vice versa */
    	    if (ac < 2)
    		usage (progname);
    	    frac = atof (*++av);
    	    setfrac (frac);
    	    ac--;
    	    break;
    	case 'r':	/* Initial rotation angle in degrees */
    	    if (ac < 2)
    		usage (progname);
    	    rot = atof (*++av);
    	    setrot (rot);
    	    ac--;
    	    break;
    	case 's':	/* Initial scale factor in arcseconds per pixel */
    	    if (ac < 2)
    		usage (progname);
    	    scale = atof (*++av);
    	    setsecpix (scale);
    	    ac--;
    	    break;
    	case 't':	/* +/- this many pixels is a hit */
    	    if (ac < 2)
    		usage (progname);
    	    tolerance = atoi (*++av);
    	    settolerance (tolerance);
    	    ac--;
    	    break;
    	case 'w':	/* update the fields */
    	    writeheader++;
    	    break;
    	default:
    	    usage(progname);
    	    break;
    	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage (progname);

    if (delete && !writeheader) {
	fprintf (stderr, "%s: -d also requires -w.\n", progname);
	usage(progname);
	}

    if (!delete && !writeheader && !verbose) {
	fprintf (stderr, "%s: need at least one of -dwv.\n", progname);
	usage(progname);
	}

    if (delete)
	while (ac-- > 0) {
    	char *fn = *av++;
    	delPos (fn);
	}
    else
	while (ac-- > 0) {
    	char *fn = *av++;
    	if (verbose)
    	    printf ("%s:\n", fn);
    	calPos (fn);
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
    fprintf (stderr,"Set WCS in FITS and IRAF image files\n");
    fprintf (stderr, "By E. Downey, UIowa and D. Mink, SAO\n");
    fprintf(stderr,"%s: usage: [-vowd] [-m mag] file.fts ...\n", progname);
    fprintf(stderr,"  -d: delete any WCS header fields (requires -w)\n");
    fprintf(stderr,"  -f: image flipped around N-S axis\n");
    fprintf(stderr,"  -j: use Joe Mohr's star matching subroutine\n");
    fprintf(stderr,"  -m: initial GSC limiting magnitude (default 17)\n");
    fprintf(stderr,"  -o: allow overwriting any existing WCS fields\n");
    fprintf(stderr,"  -r: initial rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -s: initial scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -t: offset tolerance in pixels (default 20)\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -w: write header (default is read-only)\n");
    exit (1);
}

static void
calPos (name)
char *name;
{
    char *image;		/* FITS image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int iraffile;		/* 1 if IRAF image */
    int *irafheader;		/* IRAF image header */

    /* Allocate FITS header */
    header = malloc (MAXHEADLEN);
    lhead = MAXHEADLEN;

    /* Open IRAF image if .imh extension is present */
    if (strsrch (name,".imh") != NULL) {
	iraffile = 1;
	irafheader = irafrhead (name, lhead, header);
	if (irafheader == NULL) {
	    free (header);
	    return;
	    }
	image = irafrimage (name, irafheader, header);
	if (image == NULL) {
	    free (header);
	    free (irafheader);
	    return;
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	image = fitsrimage (name, lhead, header);
	if (image == NULL) {
	    free (header);
	    return;
	    }
	}
    if (verbose) {
	fprintf (stderr,"%s\n",RevMsg);
	fprintf (stderr,"Set World Coordinate System in FITS and IRAF image files\n");
	fprintf (stderr, "By E. Downey, UIowa and D. Mink, SAO\n");
	}

    /* print existing WCS headers and check for permission to overwrite */
    if (checkWCSFITS (header) == 0 && writeheader && !overwrite) {
	fprintf (stderr,
    	    "%s: use -o to overwrite existing WCS fields.\n", name);
	free (header);
	if (iraffile)
	    free (irafheader);
	free (image);
	return;
	}

    if (setWCSFITS (header, image, verbose)) {
	if (writeheader) {
	    if (verbose)
		(void) checkWCSFITS (header);	/* print new WCS */

	    /* insure there is at least one of EPOCH or EQUINOX fields */
	    checkEquinox (header);

	    if (iraffile) {
		if (irafwhead (name, irafheader, header) > 0 && verbose)
		    printf ("%s: rewritten successfully.\n", name);
		}
	    else {
		if (fitswimage (name, header, image) > 0 && verbose)
		    printf ("%s: rewritten successfully.\n", name);
		}
	    }
	}
    else if (verbose)
	printf ("%s: no -w -- file unchanged.\n", name);

    free (header);
    if (iraffile)
	free (irafheader);
    free (image);
    return;
}


static void
delPos (name)
char *name;
{
    char *header;
    char *image;
    int lhead;

    /* open the image */
    lhead = MAXHEADLEN;
    header = malloc (lhead);
    image = fitsrimage (name, lhead, header);
    if (image == NULL) {
	fprintf (stderr, "Cannot read from %s\n", name);
	free (header);
	return;
	}

    if (delWCSFITS (header, verbose) < 1) {
	if (verbose)
	    printf ("%s: no WCS fields found -- file unchanged\n", name);
	}
    else  {
	if (fitswimage (name, header, image) < 1)
	    fprintf (stderr, "%s: Could not write FITS file\n", name);
	else if (verbose)
	    printf ("%s: rewritten successfully without WCS\n", name);
	}

    free (header);
    free (image);
    return;
}


/* insure header has at least one of EPOCH or EQUINOX.
 * if neither, add both set to 2000.0.
 */
static void
checkEquinox (header)
char *header;
{
    static char ep[] = "EPOCH";
    static char eq[] = "EQUINOX";
    double p, q;

    if (hgetr8 (header, ep, &p) < 0 && hgetr8 (header, eq, &q) < 0) {
	hputr8 (header, ep, 2000.0);
	hputr8 (header, eq, 2000.0);
    }
}

/* check the WCS fields and print any that are found if verbose.
 * return 0 if all are found, else -1.
 */
static int
checkWCSFITS (header)
char *header;

{
    char str[16];
    double v;
    int n;

    n = 0;

    if (hgets (header,"CTYPE1",16,str)) {
	if (verbose) printf ("CTYPE1 = %s\n", str);
	free (str);
	n++;
	}
    if (hgetr8 (header, "CRVAL1", &v)) {
	if (verbose) printf ("CRVAL1 = %g\n", v);
	n++;
	}
    if (hgetr8 (header, "CDELT1", &v)) {
	if (verbose) printf ("CDELT1 = %g\n", v);
	n++;
	}
    if (hgetr8 (header, "CRPIX1", &v)) {
	if (verbose) printf ("CRPIX1 = %g\n", v);
	n++;
	}
    if (hgetr8 (header, "CROTA1", &v)) {
	if (verbose) printf ("CROTA1 = %g\n", v);
	n++;
	}

    if (hgets (header,"CTYPE2",16,str)) {
	if (verbose) printf ("CTYPE2 = %s\n", str);
	free (str);
	n++;
	}
    if (hgetr8 (header, "CRVAL2", &v)) {
	if (verbose) printf ("CRVAL2 = %g\n", v);
	n++;
	}
    if (hgetr8 (header, "CDELT2", &v)) {
	if (verbose) printf ("CDELT2 = %g\n", v);
	n++;
	}
    if (hgetr8 (header, "CRPIX2", &v)) {
	if (verbose) printf ("CRPIX2 = %g\n", v);
	n++;
	}
    if (hgetr8 (header, "CROTA2", &v)) {
	if (verbose) printf ("CROTA2 = %g\n", v);
	n++;
	}

    return (n == 10 ? 0 : -1);
}


/* delete all the C* fields.
 * return 0 if at least one such field is found, else -1.  */

int
delWCSFITS (header, verbose)

char *header;
int verbose;
{
    static char *flds[] = {
	"CTYPE1", "CRVAL1", "CDELT1", "CRPIX1", "CROTA1",
	"CTYPE2", "CRVAL2", "CDELT2", "CRPIX2", "CROTA2" };
    int i;
    int n;

    n = 0;

    for (i = 0; i < sizeof(flds)/sizeof(flds[0]); i++) {
	if (hdel (header, flds[i])) {
	    n++;
	    if (verbose)
    	    printf ("%s: deleted\n", flds[i]);
	}
	else if (verbose)
	    printf ("%s: not found\n", flds[i]);
	}

    return (n);
}

char *
getRevMsg ()
{
    return (RevMsg);
}
