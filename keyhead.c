/* File keyhead.c
 * November 30, 1998
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

#define MAXKWD 50
#define MAXFILES 1000

static void usage();
static void ChangeKeyNames ();

static int verbose = 0;		/* verbose/debugging flag */
static int newimage = 0;	/* write new image with modified header */
static int replace = 0;		/* replace value of first keyword with second */
static int nfile = 0;
static int maxlfn = 0;
static int keyset = 0;
static int histset = 0;
static int version = 0;		/* If 1, print only program name and version */


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

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {

	case 'h':	/* set HISTORY */
	    histset++;
	    break;

	case 'k':	/* set SETHEAD keyword */
	    keyset++;
	    break;

	case 'n':	/* write to new file */
	    newimage = 1;
	    break;

	case 'r':	/* write to new file */
	    replace = 1;
	    break;

	case 'v':	/* more verbosity */
	    verbose++;
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
	if (strsrch (*av,"=") != NULL) {
	    kwd[nkwd] = *av;
	    nkwd++;
	    }
	else {
	    fn[nfile] = *av;
	    lfn = strlen (fn[nfile]);
	    if (lfn > maxlfn)
		maxlfn = lfn;
	    nfile++;
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
    if (version)
	exit (-1);
    fprintf (stderr,"Change FITS or IRAF header keyword names\n");
    fprintf(stderr,"Usage: [-vn] file.fit [file.fits...] kw1=kw1a ... kwn=kwna\n");
    fprintf(stderr,"Usage: [-nrv] @filelist kw1=kw1a kw2=kw2a ... kwn=kwna\n");
    fprintf(stderr,"  -n: write new file\n");
    fprintf(stderr,"  -r: replace value of 1st keyword with value of 2nd keyword\n");
    fprintf(stderr,"  -v: verbose\n");
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
    char *irafheader;	/* IRAF image header */
    int iraffile;	/* 1 if IRAF image, 0 if FITS image */
    int i, lname, lext, lroot, lhist;
    char *image, *imext, *imext1;
    char newname[128];
    char newkey[16];
    char *ext, *fname;
    char *kw, *kwv, *kwl, *kwn;
    char *value, *q, *line;
    char echar;
    int ikwd, lkwd, lkwn;
    int squote = 39;
    int dquote = 34;
    char cval[24];
    int fdr, fdw, ipos, nbr, nbw, naxis, nchange;
    char history[72];
    char *endchar;
    char *ltime;

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
	fprintf (stderr,"Change Header Keyword Names from ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", filename);
	else
	    fprintf (stderr,"FITS image file %s\n", filename);
	}

    if (nkwd < 1)
	return;

    /* Change keywords one at a time */
    nchange = 0;
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

	/* Change keyword value */
	if (replace) {
	    if ((line = ksearch (header, kwn)) == NULL)
		continue;
	    q = strchr (line, squote);
	    value = hgetc (header, kwn);
	    if (q != NULL && q < line+80) {
		hputs (header, kwd[ikwd], value);
		if (verbose)
		    printf ("%s = %s = '%s'\n", kwd[ikwd], kwn, value);
		}
	    else {
		hputc (header, kwd[ikwd], value);
		if (verbose)
		    printf ("%s = %s = %s\n", kwd[ikwd], kwn, value);
		}
	    }

	/* Change keyword name */
	else {
	    if ((line = ksearch (header, kwd[ikwd])) == NULL)
		continue;
	    hchange (header, kwd[ikwd], kwn);
	    if (verbose)
		printf ("%s => %s\n", kwd[ikwd], kwn);
	    }
	*kwv = '=';
	nchange++;
	}

    if (!nchange)
	return;

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

    /* Add history to header */
    if (keyset || histset) {
	if (hgets (header, "KEYHEAD", 72, history))
	    hputc (header, "HISTORY", history);
	endchar = strchr (history, ',');
	*endchar = (char) 0;
	strcat (history, " ");
	ltime = getltime ();
	strcat (history, ltime);
	endchar = strrchr (history,':');
	*endchar = (char) 0;
	strcat (history, " ");
	for (ikwd = 0; ikwd < nkwd; ikwd++) {
	    lhist = strlen (history);
	    lkwd = strlen (kwd[ikwd]);

	    /* If too may keywords, start a second history line */
	    if (lhist + lkwd > 71) {
		if (histset) {
		    strcat (history, " updated");
		    hputc (header, "HISTORY", history);
		    endchar = strchr (history, ',');
		    *endchar = (char) 0;
		    strcat (history, " ");
		    ltime = getltime ();
		    strcat (history, ltime);
		    endchar = strrchr (history,':');
		    *endchar = (char) 0;
		    strcat (history, " ");
		    }
		else
		    break;
		}
	    strcat (history, kwd[ikwd]);
	    if (nkwd == 2 && ikwd < nkwd-1)
		strcat (history, " and ");
	    else if (ikwd < nkwd-1)
		strcat (history, ", ");
	    }
	if (keyset)
	    hputs (header, "KEYHEAD", history);
	if (histset)
	    hputc (header, "HISTORY", history);
	}

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

/* Dec 17 1997	New program
 *
 * May 28 1998	Include fitsio.h instead of fitshead.h
 * Jul 24 1998	Make irafheader char instead of int
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Aug 14 1998	Preserve extension when creating new file name
 * Aug 14 1998	If changing primary header, write out entire input file
 * Oct  5 1998	Allow header changes even if no data is present
 * Oct  5 1998	Determine assignment arguments by presence of equal sign
 * Oct  5 1998	Use isiraf() to check for file type
 * Oct 28 1998	Add option to replace keyword value with that of other keyword
 * Nov 25 1998	Fix bug in keyword name change code
 * Nov 30 1998	Add version and help commands for consistency
 */
