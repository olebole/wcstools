/* File keyhead.c
 * May 27, 1998
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
#include "fitsio.h"

#define MAXKWD 50
#define MAXFILES 1000

static void usage();
static void ChangeKeyNames ();

static int verbose = 0;		/* verbose/debugging flag */
static int newimage = 0;
static int nfile = 0;
static int maxlfn = 0;


main (ac, av)
int ac;
char **av;
{
    char *str;
    char *kwd[MAXKWD];
    int ifile;
    int nkwd = 0;
    char *fn[MAXFILES];
    int readlist = 0;
    int lfn;
    char *lastchar;
    char filename[128];
    FILE *flist;
    char *listfile;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {

	case 'v':	/* more verbosity */
	    verbose++;
	    break;

	case 'n':	/* write to new file */
	    newimage = 1;
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

    /* now there are ac remaining file names and keywords starting at av[0] */
    if (ac == 0)
	usage ();

    /* Make keyword and filename lists from rest of arguments */
    nkwd = 0;
    nfile = 0;
    while (ac-- > 0  && nfile < MAXFILES && nkwd < MAXKWD) {
	if (strsrch (*av,".fit") != NULL ||
	    strsrch (*av,".fts") != NULL ||
	    strsrch (*av,".FIT") != NULL ||
	    strsrch (*av,".FTS") != NULL ||
	    strsrch (*av,".imh") != NULL) {
	    fn[nfile] = *av;
	    lfn = strlen (fn[nfile]);
	    if (lfn > maxlfn)
		maxlfn = lfn;
	    nfile++;
	    }
	else {
	    kwd[nkwd] = *av;
	    nkwd++;
	    }
	av++;
	}

    /* Change keyword names one file at a time */
    if (nkwd > 0) {

	/* Read through headers of images in listfile */
	if (readlist) {
	    if ((flist = fopen (listfile, "r")) == NULL) {
		fprintf (stderr,"GETHEAD: List file %s cannot be read\n",
		     listfile);
		usage ();
		}
	    while (fgets (filename, 128, flist) != NULL) {
		lastchar = filename + strlen (filename) - 1;
		if (*lastchar < 32) *lastchar = 0;
		ChangeKeyNames (filename, nkwd, kwd);
		if (verbose)
		    printf ("\n");
		}
	    fclose (flist);
	    }

	/* Read image headers from command line list */
	else {
	    for (ifile = 0; ifile < nfile; ifile++)
	    	ChangeKeyNames (fn[ifile], nkwd, kwd);
	    }
	}
    return (0);
}

static void
usage ()
{
    fprintf (stderr,"Change FITS or IRAF header keyword names\n");
    fprintf(stderr,"Usage: [-vn] file.fit [file.fits...] kw1=kw1a ... kwn=kwna\n");
    fprintf(stderr,"Usage: [-vn] @filelist kw1=kw1a kw2=kw2a ... kwn=kwna\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -n: write new file\n");
    exit (1);
}


static void
ChangeKeyNames (filename, nkwd, kwd)

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
    char *kw, *kwv, *kwl, *kwn;
    int ikwd, lkwd, lkwn;
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
	fprintf (stderr,"Change Header Keyword Names from ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", filename);
	else
	    fprintf (stderr,"FITS image file %s\n", filename);
	}

    if (nkwd < 1)
	return;

    /* Change keywords one at a time */
    for (ikwd = 0; ikwd < nkwd; ikwd++) {
        strcpy (cval,"                    ");
	kwv = strchr (kwd[ikwd],'=');
	*kwv = 0;
	lkwd = kwv - kwd[ikwd];
	kwn = kwv + 1;
	lkwn = strlen (kwn);

	/* Make keywords all upper case */
	kwl = kwd[ikwd] + lkwd;
	for (kw = kwd[ikwd]; kw < kwl; kw++) {
	    if (*kw > 96 && *kw < 123)
		*kw = *kw - 32;
	    }
	kwl = kwn + lkwn;
	for (kw = kwn; kw < kwl; kw++) {
	    if (*kw > 96 && *kw < 123)
		*kw = *kw - 32;
	    }

	/* Change keyword name */
	hchange (header, kwd[ikwd], kwn);
	if (verbose)
	    printf ("%s => %s\n", kwd[ikwd], kwn);
	*kwv = '=';
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

/* Dec 17 1997	New program
 *
 * May 28 1998	Include fitsio.h instead of fitshead.h
 */
