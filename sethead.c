/* File sethead.c
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

#define MAXKWD 50
#define MAXFILES 1000
static int maxnkwd = MAXKWD;
static int maxnfile = MAXFILES;

static void usage();
static void SetValues ();

static int verbose = 0;		/* verbose/debugging flag */
static int newimage0 = 0;
static int listpath = 0;
static int keyset = 0;
static int histset = 0;
static int krename = 0;
static char prefix[2];
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
    char *name;
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

	case 'f':	/* maximum number of files, overrides MAXFILES */
	    if (ac < 2)
		usage();
	    maxnfile = (int) (atof (*++av));
	    ac--;
	    break;

	case 'h':	/* set HISTORY */
	    histset++;
	    break;

	case 'k':	/* set SETHEAD keyword */
	    keyset++;
	    break;

	case 'm':	/* maximum number of keywords, overrides MAXKWD */
	    if (ac < 2)
		usage();
	    maxnkwd = (int) (atof (*++av));
	    ac--;
	    break;

	case 'n':	/* write new file */
	    newimage0++;
	    break;

	case 'r':	/* set flag to krename keywords with replaced values */
	    krename++;
	    if (ac > 1) {
		strncpy (prefix, *++av, 1);
		ac--;
		}
	    else
		prefix[0] = 'X';
	    prefix[1] = (char) 0;
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
	if (strsrch (*av,"@") != NULL) {
	    readlist++;
	    listfile = *av + 1;
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
    /* Get keyword values one at a time */
    if (nkwd > 0) {

	/* Read through headers of images in listfile */
	if (readlist) {
	    if ((flist = fopen (listfile, "r")) == NULL) {
		fprintf (stderr,"SETHEAD: List file %s cannot be read\n",
		     listfile);
		usage ();
		}
	    while (fgets (filename, 128, flist) != NULL) {
		lastchar = filename + strlen (filename) - 1;
		if (*lastchar < 32) *lastchar = 0;
		SetValues (filename, nkwd, kwd);
		if (verbose)
		    printf ("\n");
		}
	    fclose (flist);
	    }

	/* Read image headers from command line list */
	else {
	    for (ifile = 0; ifile < nfile; ifile++)
	    	SetValues (fn[ifile], nkwd, kwd);
	    }
	}

    /* Print error message(s) if the command line got too long */
    if (nxfile > 0)
	fprintf (stderr, "SETHEAD limit of %d files exceeded by %d\n",
		 maxnfile, nxfile);
    if (nxkwd > 0)
	fprintf (stderr, "SETHEAD limit of %d keywords exceeded by %d\n",
		 maxnkwd, nxkwd);

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Set FITS or IRAF header keyword values\n");
    fprintf(stderr,"Usage: [-nhkv][-f num][-m num][-r char] file1.fits [... filen.fits] kw1=val1 [ ... kwn=valuen]\n");
    fprintf(stderr,"Usage: [-nhkv][-f num][-m num][-r char] @listfile kw1=val1 [ ... kwn=valuen]\n");
    fprintf(stderr,"  -f: Max number of files if not %d\n", maxnfile);
    fprintf(stderr,"  -h: Write HISTORY line\n");
    fprintf(stderr,"  -m: Max number of keywords if not %d\n", maxnkwd);
    fprintf(stderr,"  -k: Write SETHEAD keyword\n");
    fprintf(stderr,"  -n: Write a new file (add e before the extension)\n");
    fprintf(stderr,"  -r: Rename reset keywords with prefix following\n");
    fprintf(stderr,"  -v: Verbose\n");
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
    char *irafheader;	/* IRAF image header */
    int iraffile;	/* 1 if IRAF image, 0 if FITS image */
    int newimage;	/* 1 to awrite new image file, else 0 */
    int i, lname, lext, lroot;
    char *image;
    char newname[128];
    char *ext, *fname, *imext, *imext1;
    char *kw, *kwv, *kwl, *kwv0;
    char *v, *vq0, *vq1;
    char echar;
    int ikwd, lkwd, lkwv, lhist;
    int fdr, fdw, ipos, nbr, nbw;
    int squote = 39;
    int dquote = 34;
    int naxis = 1;
    int imageread = 0;
    char cval[24];
    char history[72];
    char *endchar;
    char *ltime;
    char newkey[10];

    newimage = newimage0;

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
			fprintf (stderr, "No FITS image in %s\n", filename);
		    imageread = 0;
		    }
		else
		    imageread = 1;
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
	fprintf (stderr,"Set Header Parameter Values in ");
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
	kwv0 = strchr (kwd[ikwd],'=');
	*kwv0 = 0;
	lkwd = kwv0 - kwd[ikwd];
	kwv = kwv0 + 1;
	lkwv = strlen (kwv);

	/* Make keyword all upper case */
	kwl = kwd[ikwd] + lkwd;
	for (kw = kwd[ikwd]; kw < kwl; kw++) {
	    if (*kw > 96 && *kw < 123)
		*kw = *kw - 32;
	    }

	/* If keyword is already in header, krename it if requested */
	if (krename && ksearch (header, kwd[ikwd])) {
	    strcpy (newkey, prefix);
	    strcat (newkey, kwd[ikwd]);
	    if (strlen (newkey) > 8)
		newkey[8] = (char) 0;
	    hchange (header, kwd[ikwd], newkey);
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
	else if (!strcmp (kwv,"T") || !strcmp (kwv,"t"))
	    hputl (header, kwd[ikwd], 1);
	else if (!strcmp (kwv,"YES") || !strcmp (kwv,"yes"))
	    hputl (header, kwd[ikwd], 1);
	else if (!strcmp (kwv,"F") || !strcmp (kwv,"f"))
	    hputl (header, kwd[ikwd], 0);
	else if (!strcmp (kwv,"NO") || !strcmp (kwv,"no"))
	    hputl (header, kwd[ikwd], 0);
	else
	    hputs (header, kwd[ikwd], kwv);
	if (verbose)
	    printf ("%s = %s\n", kwd[ikwd], kwv);
	*kwv0 = '=';
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

    /* Add history to header */
    if (keyset || histset) {
	if (hgets (header, "SETHEAD", 72, history))
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
	    kwv0 = strchr (kwd[ikwd],'=');
	    if (kwv0)
		*kwv0 = (char) 0;
	    lhist = strlen (history);
	    lkwd = strlen (kwd[ikwd]);

	    /* If too may keywords, start a second history line */
	    if (lhist + lkwd + 10 > 71) {
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
	    if (kwv0)
		*kwv0 = '=';
	    if (nkwd == 2 && ikwd < nkwd-1)
		strcat (history, " and ");
	    else if (ikwd < nkwd-1)
		strcat (history, ", ");
	    }
	strcat (history, " updated");
	if (keyset)
	    hputs (header, "SETHEAD", history);
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
    else if (naxis > 0 && imageread) {
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

/* Oct 11 1996	New program
 * Dec 12 1996	Move ISNUM subroutine to hget.c
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 *
 * May 28 1998	Include fitsio.h instead of fitshead.h
 * Jun  2 1998	Fix bug in hput()
 * Jul 24 1998	Deal coorectly with logical T or F
 * Jul 24 1998	Make irafheader char instead of int
 * Jul 30 1998	Allow use of list of files and multiple files
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Aug 14 1998	Preserve extension when creating new file name
 * Aug 14 1998	If changing primary header, write out entire input file
 * Aug 31 1998	Add options to add HISTORY and/or SETHEAD keyword
 * Sep  1 1998	Add option to keep changed keywords with new names
 * Oct  5 1998	Allow header changes even if no data is present
 * Oct  5 1998	Determine assignment arguments by presence of equal sign
 * Oct  5 1998	Use isiraf() to determine if file is IRAF or FITS
 * Oct 29 1998	Fix history setting
 * Nov 30 1998	Add version and help commands for consistency
 * Dec 30 1998	Write header without image if no image is present
 *
 * Mar  4 1999	Reset = for each keyword after setting value in header
 * Apr  2 1999	Add warning if too many files or keywords on command line
 * Apr  2 1999	Add -f and -m to change maximum number of files or keywords
 */
