/* File sethead.c
 * February 21, 1997
 * By Doug Mink Harvard-Smithsonian Center for Astrophysics)
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

#define MAXKWD 50

static void usage();
static void SetValues ();

static int verbose = 0;		/* verbose/debugging flag */
static int newimage = 1;
static void SetValues ();

main (ac, av)
int ac;
char **av;
{
    char *str;
    char *kwd[MAXKWD];
    int nkwd;
    char *fn;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {
	case 'v':	/* more verbosity */
	    verbose++;
	    break;
	case 'w':	/* overwrite old file */
	    newimage = 0;
	    break;
	default:
	    usage();
	    break;
	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage ();

    if (ac-- > 0) {
	fn = *av++;
	}
    nkwd = 0;
    while (ac-- > 0 && nkwd < MAXKWD) {
	kwd[nkwd] = *av++;
	nkwd++;
	}
    if (nkwd > 0)
	SetValues (fn, nkwd, kwd);

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"Set FITS or IRAF header keyword values\n");
    fprintf(stderr,"Usage: [-v] file.fit kw1=val1 kw2=val2 ... kwn=valuen\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -w: overwrite header\n");
    exit (1);
}


static void
SetValues (filename, nkwd, kwd)

char	*filename;	/* Name of FITS or IRAF image file */
int	nkwd;		/* Number of keywords for which to set values */
char	*kwd[];		/* Names and values of those keywords */

{
    char *header;	/* FITS image header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    int *irafheader;	/* IRAF image header */
    int iraffile;	/* 1 if IRAF image, 0 if FITS image */
    int i, lname, lext;
    char *image;
    char newname[128];
    char *ext, *fname;
    char *kw, *kwv, *kwl;
    char *v, *vq0, *vq1;
    int ikwd, lkwd, lkwv;
    int squote = 39;
    int dquote = 34;
    char cval[24];

    /* Open IRAF image if .imh extension is present */
    if (strsrch (filename,".imh") != NULL) {
	iraffile = 1;
	if ((irafheader = irafrhead (filename, &lhead)) != NULL) {
	    if ((header = iraf2fits (filename, irafheader, lhead, &nbhead)) == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",filename);
		free (irafheader);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF file %s\n", filename);
	    return;
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	if ((header = fitsrhead (filename, &lhead, &nbhead)) != NULL) {
	    if ((image = fitsrimage (filename, nbhead, header)) == NULL) {
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
    if (verbose) {
	fprintf (stderr,"Print Header Parameter Values from ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", filename);
	else
	    fprintf (stderr,"FITS image file %s\n", filename);
	}

    if (nkwd < 1)
	return;

    /* Set keywords one at a time */
    for (ikwd = 0; ikwd < nkwd; ikwd++) {
        strcpy (cval,"                    ");
	kwv = strchr (kwd[ikwd],'=');
	*kwv = 0;
	lkwd = kwv - kwd[ikwd];
	kwv = kwv + 1;
	lkwv = strlen (kwv);

	/* Make keyword all upper case */
	kwl = kwd[ikwd] + lkwd;
	for (kw = kwd[ikwd]; kw < kwl; kw++) {
	    if (*kw > 96 && *kw < 123)
		*kw = *kw - 32;
	    }

	/* Write value to keyword */
	if ((vq0 = strchr (kwv,dquote))) {
	    vq0 = vq0 + 1;
	    vq1 = strchr (vq0,dquote);
	    if (vq0 && vq1) {
		kwv = vq0;
		*vq1 = 0;
		hputs (header, kwd[ikwd], kwv);
		}
	    else
		hputs (header, kwd[ikwd], kwv);
	    }
	else if ((vq0 = strchr (kwv,squote))) {
	    vq0 = vq0 + 1;
	    vq1 = strchr (vq0,squote);
	    if (vq0 && vq1) {
		kwv = vq0;
		*vq1 = 0;
		hputs (header, kwd[ikwd], kwv);
		}
	    else
		hputs (header, kwd[ikwd], kwv);
	    }
	else if (isnum (kwv)) {
	    i = 21 - lkwv;
	    for (v = kwv; v < kwv+lkwv; v++)
		cval[i++] = *v;
	    cval[21] = 0;
	    hputc (header, kwd[ikwd], cval);
	    }
	else
	    hputs (header, kwd[ikwd], kwv);
	if (verbose)
	    printf ("%s = %s\n", kwd[ikwd], kwv);
	}

    /* Make up name for new FITS or IRAF output file */
    if (newimage) {

    /* Remove directory path and extension from file name */
	ext = strrchr (filename, '.');
	fname = strrchr (filename, '/');
	if (fname)
	    fname = fname + 1;
	else
	    fname = filename;
	lname = strlen (fname);
	if (ext) {
	    lext = strlen (ext);
	    strncpy (newname, fname, lname - lext);
	    *(newname + lname - lext) = 0;
	    }
	else
	    strcpy (newname, fname);

    /* Add file extension preceded by a e */
	if (iraffile)
	    strcat (newname, "e.imh");
	else
	    strcat (newname, "e.fit");
	}
    else
	strcpy (newname, filename);

    /* Write fixed header to output file */
    if (iraffile) {
	if (irafwhead (newname, lhead, irafheader, header) > 0 && verbose)
	    printf ("%s rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
	free (irafheader);
	}
    else {
	if (fitswimage (newname, header, image) > 0 && verbose)
	    printf ("%s: rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
	free (image);
	}

    free (header);
    return;
}

/* Oct 11 1996	New program
 * Dec 12 1996	Move ISNUM subroutine to hget.c
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 */
