/* File gethead.c
 * March 12, 1998
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

#define MAXKWD 100
#define MAXFILES 1000

static void usage();
static void PrintValues();
static void strclean();

static int verbose = 0;		/* verbose/debugging flag */
static int nfile = 0;
static int maxlfn = 0;
static int listall = 0;
static int tabout = 0;
static int printhead = 0;

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
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
    int ikwd, lkwd, i;
    char *kw, *kwe;
    char string[80];

    /* crack arguments */
    for (av++; --ac > 0 && (*(str = *av)=='-'); av++) {
	char c;
	while (c = *++str)
	switch (c) {

	case 'a':	/* list file even if the keywords are not found */
	    listall++;
	    break;

	case 'v':	/* more verbosity */
	    verbose++;
	    break;

	case 'h':	/* output column headings */
	    printhead++;
	    break;

	case 't':	/* output tab table */
	    tabout++;
	    break;

	default:
	    usage(progname);
	    break;
	}
    }

    /* If no arguments left, print usage */
    if (ac == 0)
	usage (progname);

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
	else if (strsrch (*av,"@") != NULL) {
	    readlist++;
	    listfile = *av + 1;
	    nfile = 2;
	    }
	else {
	    kwd[nkwd] = *av;
	    nkwd++;
	    }
	av++;
	}

    if (nkwd > 0) {

	/* Print column headings if tab table or headings requested */
	if (tabout || printhead) {
	    printf ("FILENAME");
	    if (maxlfn > 8) {
		for (i = 8; i < maxlfn; i++)
		    printf (" ");
		}
	    if (tabout)
	    	printf ("	");
	    else
		printf (" ");
	    for (ikwd = 0; ikwd < nkwd; ikwd++) {
		lkwd = strlen (kwd[ikwd]);
		kwe = kwd[ikwd] + lkwd;
		for (kw = kwd[ikwd]; kw < kwe; kw++) {
		    if (*kw > 96 && *kw < 123)
			*kw = *kw - 32;
		    }
		printf ("%s",kwd[ikwd]);
		if (verbose || ikwd == nkwd - 1)
	    	    printf ("\n");
		else if (tabout)
	    	    printf ("	");
		else
		    printf (" ");
		}
	    }
	/* Print field-defining hyphens if tab table output requested */
	if (tabout) {
	    printf ("--------");
	    if (maxlfn > 8) {
		for (i = 8; i < maxlfn; i++)
		    printf (" ");
		}
	    if (tabout)
	    	printf ("	");
	    else
		printf (" ");
	    for (ikwd = 0; ikwd < nkwd; ikwd++) {
		strcpy (string, "----------");
		lkwd = strlen (kwd[ikwd]);
		string[lkwd] = 0;
		printf ("%s",string);
		if (verbose || ikwd == nkwd - 1)
	    	    printf ("\n");
		else if (tabout)
	    	    printf ("	");
		else
		    printf (" ");
		}
	    }

    /* Get keyword values one at a time */

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
		PrintValues (filename, nkwd, kwd);
		if (verbose)
		    printf ("\n");
		}
	    fclose (flist);
	    }

	/* Read image headers from command line list */
	else {
	    for (ifile = 0; ifile < nfile; ifile++)
	    	PrintValues (fn[ifile], nkwd, kwd);
	    }
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"Print FITS or IRAF header keyword values\n");
    fprintf(stderr,"usage: [-ahtv] file1.fit ... filen.fits kw1 kw2 ... kwn\n");
    fprintf(stderr,"usage: [-ahtv] @filelist kw1 kw2 ... kwn\n");
    fprintf(stderr,"  -a: list file even if keywords are not found\n");
    fprintf(stderr,"  -h: print column headings\n");
    fprintf(stderr,"  -t: output in tab-separated table format\n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
PrintValues (name, nkwd, kwd)

char	*name;	/* Name of FITS or IRAF image file */
int	nkwd;	/* Number of keywords for which to print values */
char	*kwd[];	/* Names of keywords for which to print values */

{
    char *header;	/* FITS image header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    int *irafheader;	/* IRAF image header */
    int iraffile;
    char fnform[8];
    char string[80];
    char temp[80];
    char outline[1000];
    char mstring[800];
    char *kw, *kwe;
    int ikwd, lkwd, nfound;

    /* Open IRAF image if .imh extension is present */
    if (strsrch (name,".imh") != NULL) {
	iraffile = 1;
	if ((irafheader = irafrhead (name, &lhead)) != NULL) {
	    header = iraf2fits (name, irafheader, lhead, &nbhead);
	    free (irafheader);
	    if (header == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",name);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF file %s\n", name);
	    return;
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	if ((header = fitsrhead (name, &lhead, &nbhead)) == NULL) {
	    fprintf (stderr, "Cannot read FITS file %s\n", name);
	    return;
	    }
	}
    if (verbose) {
	fprintf (stderr,"Print Header Parameter Values from ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}
    if (nfile > 1) {
	if (tabout)
	    sprintf (fnform, "%%-%ds	", maxlfn);
	else
	    sprintf (fnform, "%%-%ds ", maxlfn);
	sprintf (outline, fnform, name);
	}
    nfound = 0;

    for (ikwd = 0; ikwd < nkwd; ikwd++) {
	lkwd = strlen (kwd[ikwd]);
	kwe = kwd[ikwd] + lkwd;
	for (kw = kwd[ikwd]; kw < kwe; kw++) {
	    if (*kw > 96 && *kw < 123)
		*kw = *kw - 32;
	    }
	if (hgets (header, kwd[ikwd], 80, string)) {
	    strclean (string);
	    if (verbose)
		printf ("%s = %s", kwd[ikwd], string);
	    else
		strcat (outline, string);
	    nfound++;
	    }
	else if (hgetm (header, kwd[ikwd], 600, mstring)) {
	    if (verbose)
		printf ("%s = %s", kwd[ikwd], mstring);
	    else
		strcat (outline, mstring);
	    nfound++;
	    }
	else if (verbose)
	    printf ("%s not found", kwd[ikwd]);
	else
	    strcat (outline, "___");

	if (verbose)
	    printf ("\n");
	else if (ikwd < nkwd-1) {
	    if (tabout)
		strcat (outline, "	");
	    else
		strcat (outline, " ");
	    }
	}
    if (!verbose && (nfile < 2 || nfound > 0 || listall))
	printf ("%s\n", outline);

    free (header);
    return;
}


/* Remove exponent and trailing zeroes, if reasonable */
static void
strclean (string)

char *string;

{
    char *sdot, *s;
    int ndek, lstr, i;

    /* Remove positive exponent if there are enough digits given */
    if (strsrch (string, "E+") != NULL) {
	lstr = strlen (string);
	ndek = (int) (string[lstr-1] - 48);
	ndek = ndek + (10 * ((int) (string[lstr-2] - 48)));
	if (ndek < lstr - 7) {
	    lstr = lstr - 4;
	    string[lstr] = (char) 0;
	    string[lstr+1] = (char) 0;
	    string[lstr+2] = (char) 0;
	    string[lstr+3] = (char) 0;
	    sdot = strchr (string, '.');
	    if (ndek > 0 && sdot != NULL) {
		for (i = 1; i <= ndek; i++) {
		    *sdot = *(sdot+1);
		    sdot++;
		    *sdot = '.';
		    }
		}
	    }
	}

    /* Remove trailing zeroes */
    if (strchr (string, '.') != NULL) {
	lstr = strlen (string);
	s = string + lstr - 1;
	while (*s == '0' && lstr > 1) {
	    if (*(s - 1) != '.') {
		*s = (char) 0;
		lstr --;
		}
	    s--;
	    }
	}

    /* Remove trailing decimal point */
    lstr = strlen (string);
    s = string + lstr - 1;
    if (*s == '.')
	*s = (char) 0;

    return;
	
}

/* Sep  4 1996	New program
 * Oct  8 1996	Add newline after file name in verbose mode
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 * Nov  2 1997	Allow search of multiple files from command line or list file
 * Nov 12 1997	Add option to print tab tables or column headings
 * Dec 12 1997	Read IRAF version 2 .imh files
 *
 * Jan  5 1998	Do not print line unless keyword is found
 * Jan 27 1998	Clean up scientific notation and trailing zeroes
 * Mar 12 1998	Read IRAF multi-line keyword values
 */
