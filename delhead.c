/* File delhead.c
 * April 2, 1999
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
#include "fitsfile.h"

#define MAXKWD 100
#define MAXFILES 1000
static int maxnkwd = MAXKWD;
static int maxnfile = MAXFILES;

static void usage();
static void DelKeywords();

static int verbose = 0;		/* verbose/debugging flag */
static int newimage = 0;
static int version = 0;		/* If 1, print only program name and version */

main (ac, av)
int ac;
char **av;
{
    char *str;
    char **kwd;
    int nkwd = 0;
    int nxkwd = 0;
    char **fn;
    int nfile = 0;
    int nxfile = 0;
    int readlist = 0;
    int ifile;
    char *lastchar;
    char filename[128];
    char keyword[16];
    char *name;
    FILE *flist;
    char *listfile;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage ();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage ();
	}

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {

	case 'f':	/* maximum number of files, overrides MAXFILES */
	    if (ac < 2)
		usage();
	    maxnfile = (int) (atof (*++av));
	    ac--;
	    break;

	case 'm':	/* maximum number of keywords, overrides MAXKWD */
	    if (ac < 2)
		usage();
	    maxnkwd = (int) (atof (*++av));
	    ac--;
	    break;

	case 'n':	/* write new file */
	    newimage++;
	    break;

	case 'v':	/* more verbosity */
	    verbose++;
	    break;

	default:
	    usage();
	    break;
	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage ();

    nkwd = 0;
    nfile = 0;
    nxfile = 0;
    nxkwd = 0;
    fn = (char **)calloc (maxnfile, sizeof(char *));
    kwd = (char **)calloc (maxnkwd, sizeof(char *));
    while (ac-- > 0) {
	if (*av[0] == '@') {
	    readlist++;
	    listfile = *av + 1;
	    if (nfile == 0)
		nfile = 2;
	    }
	else if (isfits (*av) || isiraf (*av)) {
	    if (nfile < maxnfile) {
		fn[nfile] = *av;
		nfile++;
		}
	    else
		nxfile++;
	    }
	else if (nkwd < maxnkwd) {
	    kwd[nkwd] = *av;
	    nkwd++;
	    }
	else
	    nxkwd++;
	av++;
	}

    /* Delete keyword values one file at a time */
    if (nkwd > 0) {

	/* Read through headers of images in listfile */
	if (readlist) {
	    if ((flist = fopen (listfile, "r")) == NULL) {
		fprintf (stderr,"DELHEAD: List file %s cannot be read\n",
		     listfile);
		usage ();
		}
	    while (fgets (filename, 128, flist) != NULL) {
		lastchar = filename + strlen (filename) - 1;
		if (*lastchar < 32) *lastchar = 0;
		DelKeywords (filename, nkwd, kwd);
		if (verbose)
		    printf ("\n");
		}
	    fclose (flist);
	    }

	/* Read image headers from command line list */
	else {
	    for (ifile = 0; ifile < nfile; ifile++)
	    	DelKeywords (fn[ifile], nkwd, kwd);
	    }
	}

    /* Read keywords from list file */
    else if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"DELHEAD: List file %s cannot be read\n",
		     listfile);
	    usage ();
	    }
	while (fgets (keyword, 16, flist) != NULL) {
	    lastchar = keyword + strlen (keyword) - 1;
	    if (*lastchar < 32) *lastchar = 0;
	    kwd[nkwd] = (char *) calloc (1, 16);
	    strcpy (kwd[nkwd], keyword);
	    nkwd++;
	    }
	fclose (flist);

	for (ifile = 0; ifile < nfile; ifile++)
	    DelKeywords (fn[ifile], nkwd, kwd);
	}

    /* Print error message(s) if the command line got too long */
    if (nxfile > 0)
	fprintf (stderr, "DELHEAD limit of %d files exceeded by %d\n",
		 maxnfile, nxfile);
    if (nxkwd > 0)
	fprintf (stderr, "DELHEAD limit of %d keywords exceeded by %d\n",
		 maxnkwd, nxkwd);

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Delete FITS or IRAF header keyword entries\n");
    fprintf(stderr,"Usage: [-nv][-f num][-m num] file1.fits [ ... filen.fits] kw1 [... kwn]\n");
    fprintf(stderr,"Usage: [-nv][-f num][-m num] @listfile kw1 [... kwn]\n");
    fprintf(stderr,"Usage: [-nv][-f num][-m num] file1.fits [ ... filen.fits] @keylistfile\n");
    fprintf(stderr,"  -f: Max number of files if not %d\n", maxnfile);
    fprintf(stderr,"  -m: Max number of keywords if not %d\n", maxnkwd);
    fprintf(stderr,"  -n: write new file\n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
DelKeywords (filename, nkwd, kwd)

char	*filename;	/* Name of FITS or IRAF image file */
int	nkwd;		/* Number of keywords to delete */
char	*kwd[];		/* Names of those keywords */

{
    char *header;	/* FITS image header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    char *irafheader;	/* IRAF image header */
    int iraffile;	/* 1 if IRAF image, 0 if FITS image */
    int i, lname, lext, lroot, naxis;
    char *image;
    char newname[128];
    char *ext, *fname, *imext, *imext1;
    char *kw, *kwl;
    char echar;
    int ikwd, lkwd, lkwv;
    int fdr, fdw, ipos, nbr, nbw;

    /* Open IRAF image if .imh extension is present */
    if (isiraf (filename)) {
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
	    hgeti4 (header,"NAXIS",&naxis);
	    if (naxis > 0) {
		if ((image = fitsrimage (filename, nbhead, header)) == NULL) {
		    if (verbose)
			fprintf (stderr,"Cannot read FITS image %s\n",filename);
		    }
		}
	    else {
		if (verbose)
		    fprintf (stderr, "Rewriting primary header, copying rest\n");
		newimage = 1;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    return;
	    }
	}
    if (verbose) {
	fprintf (stderr,"Delete Header Parameter Entries from ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", filename);
	else
	    fprintf (stderr,"FITS image file %s\n", filename);
	}

    if (nkwd < 1)
	return;

    /* Delete keywords one at a time */
    for (ikwd = 0; ikwd < nkwd; ikwd++) {

	/* Make keyword all upper case */
	kwl = kwd[ikwd] + strlen (kwd[ikwd]);
	for (kw = kwd[ikwd]; kw < kwl; kw++) {
	    if (*kw > 96 && *kw < 123)
		*kw = *kw - 32;
	    }

	/* Delete keyword */
	if (hdel (header, kwd[ikwd]) && verbose)
	    printf ("%s: %s deleted\n", filename, kwd[ikwd]);
	}

    /* Make up name for new FITS or IRAF output file */
    if (newimage) {

    /* Remove directory path and extension from file name */
	fname = strrchr (filename, '/');
	if (fname)
	    fname = fname + 1;
	else
	    fname = filename;
	ext = strrchr (fname, '.');
	if (ext != NULL) {
	    lext = (fname + strlen (fname)) - ext;
	    lroot = ext - fname;
	    strncpy (newname, fname, lroot);
	    *(newname + lroot) = 0;
	    }
	else {
	    lext = 0;
	    lroot = strlen (fname);
	    strcpy (newname, fname);
	    }
	imext = strchr (fname, ',');
	imext1 = NULL;
	if (imext == NULL) {
	    imext = strchr (fname, '[');
	    if (imext != NULL) {
		imext1 = strchr (fname, ']');
		*imext1 = (char) 0;
		}
	    }
	if (imext != NULL && *(imext+1) != '0') {
	    strcat (newname, "_");
	    strcat (newname, imext+1);
	    }
	if (fname)
	    fname = fname + 1;
	else
	    fname = filename;
	strcat (newname, "e");
	if (lext > 0) {
	    if (imext != NULL) {
		echar = *imext;
		*imext = (char) 0;
		strcat (newname, ext);
		*imext = echar;
		if (imext1 != NULL)
		    *imext1 = ']';
		}
	    else
		strcat (newname, ext);
	    }
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
    else if (naxis > 0) {
	if (fitswimage (newname, header, image) > 0 && verbose)
	    printf ("%s: rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
	free (image);
	}
    else {
	if ((fdw = fitswhead (newname, header)) > 0) {
	    fdr = fitsropen (filename);
	    ipos = lseek (fdr, nbhead, SEEK_SET);
	    image = (char *) calloc (2880, 1);
	    while ((nbr = read (fdr, image, 2880)) > 0) {
		nbw = write (fdw, image, nbr);
		if (nbw < nbr)
		    fprintf (stderr,"SETHEAD: %d / %d bytes written\n",nbw,nbr);
		}
	    close (fdr);
	    close (fdw);
	    if (verbose)
		printf ("%s: rewritten successfully.\n", newname);
	    free (image);
	    }
	}
    free (header);
    return;
}

/* Jul 27 1998	New program
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Aug 14 1998	If changing primary header, write out entire input file
 * Oct  5 1998	Allow header changes even if no data is present
 * Oct  5 1998	Use isiraf() and isfits() to check for data file
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Mar  2 1999	Add option to delete list of keyword names from file
 * Apr  2 1999	Add warning if too many files or keywords on command line
 * Apr  2 1999	Add -f and -m to change maximum number of files or keywords
 */
