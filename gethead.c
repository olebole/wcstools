/* File gethead.c
 * June 21, 1999
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
#define MAXFILES 2000
static int maxnkwd = MAXKWD;
static int maxnfile = MAXFILES;

static void usage();
static void PrintValues();
static void strclean();
extern char *GetFITShead();

static int verbose = 0;		/* verbose/debugging flag */
static int nfile = 0;
static int ndec = -9;
static int maxlfn = 0;
static int listall = 0;
static int listpath = 0;
static int tabout = 0;
static int printhead = 0;
static int version = 0;		/* If 1, print only program name and version */
static int printfill=0;		/* If 1, print ___ for unfound keyword values */
static char *rootdir=NULL;	/* Root directory for input files */

main (ac, av)
int ac;
char **av;
{
    char *str;
    char **kwd;
    int nkwd = 0;
    int nxkwd = 0;
    char **fn;
    int nxfile = 0;
    int ifile;
    int readlist = 0;
    int lfn;
    char *lastchar;
    char filename[128];
    char *name;
    FILE *flist;
    char *listfile;
    int ikwd, lkwd, i;
    char *kw, *kwe;
    char string[80];

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    /* crack arguments */
    for (av++; --ac > 0 && (*(str = *av)=='-'); av++) {
	char c;
	while (c = *++str)
	switch (c) {

	case 'a':	/* list file even if the keywords are not found */
	    listall++;
	    break;

	case 'd':	/* root directory for input */
	    if (ac < 2)
		usage();
	    rootdir = *++av;
	    ac--;
	    break;
	    break;

	case 'f':	/* maximum number of files, overrides MAXFILES */
	    if (ac < 2)
		usage();
	    maxnfile = (int) (atof (*++av));
	    ac--;
	    break;

	case 'h':	/* output column headings */
	    printhead++;
	    break;

	case 'm':	/* maximum number of keywords, overrides MAXKWD */
	    if (ac < 2)
		usage();
	    maxnkwd = (int) (atof (*++av));
	    ac--;
	    break;

	case 'n':	/* number of decimal places in output */
	    if (ac < 2)
		usage();
	    ndec = (int) (atof (*++av));
	    ac--;
	    break;

	case 'p':	/* output column headings */
	    listpath++;
	    break;

	case 't':	/* output tab table */
	    tabout++;
	    break;

	case 'v':	/* more verbosity */
	    verbose++;
	    break;

	default:
	    usage();
	    break;
	}
    }

    /* If no arguments left, print usage */
    if (ac == 0)
	usage();

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
	    nfile = 2;
	    }
	else if (isfits (*av) || isiraf (*av)) {
	    if (nfile < maxnfile) {
		fn[nfile] = *av;

		if (listpath || (name = strrchr (fn[nfile],'/')) == NULL)
		    name = fn[nfile];
		else
		    name = name + 1;
		lfn = strlen (name);
		if (lfn > maxlfn)
		    maxlfn = lfn;
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

    if (nkwd > 0) {
	if (nkwd > 1)
	    printfill = 1;

	/* Print column headings if tab table or headings requested */
	if (printhead) {
	    if (nfile > 1 || listall) {
		printf ("FILENAME");
		if (maxlfn > 8) {
		    for (i = 8; i < maxlfn; i++)
			printf (" ");
		    }
		if (tabout)
	    	    printf ("	");
		else
		    printf (" ");
		}
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

	    /* Print field-defining hyphens if tab table output requested */
	    if (nfile > 1 || listall) {
		if (tabout) {
		    printf ("--------");
		    if (maxlfn > 8) {
			for (i = 8; i < maxlfn; i++)
			    printf (" ");
			}
		    printf ("	");
		    }
		else
		    printf (" ");
		}

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

    /* Print error message(s) if the command line got too long */
    if (nxfile > 0)
	fprintf (stderr, "GETHEAD limit of %d files exceeded by %d\n",
		 maxnfile, nxfile);
    if (nxkwd > 0)
	fprintf (stderr, "GETHEAD limit of %d keywords exceeded by %d\n",
		 maxnkwd, nxkwd);

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Print FITS or IRAF header keyword values\n");
    fprintf(stderr,"usage: gethead [-ahtv][-d dir][-f num][-m num][-n num] file1.fit ... filen.fits kw1 kw2 ... kwn\n");
    fprintf(stderr,"       gethead [-ahptv][-d dir][-f num][-m num][-n num] @filelist kw1 kw2 ... kwn\n");
    fprintf(stderr,"  -a: List file even if keywords are not found\n");
    fprintf(stderr,"  -d: Root directory for input files (default is cwd)\n");
    fprintf(stderr,"  -f: Max number of files if not %d\n", maxnfile);
    fprintf(stderr,"  -h: Print column headings\n");
    fprintf(stderr,"  -m: Max number of keywords if not %d\n", maxnkwd);
    fprintf(stderr,"  -n: Number of decimal places in numeric output\n");
    fprintf(stderr,"  -p: Print full pathnames of files\n");
    fprintf(stderr,"  -t: Output in tab-separated table format\n");
    fprintf(stderr,"  -v: Verbose\n");
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
    int iraffile;
    char fnform[8];
    char string[80];
    char temp[80];
    char *filename;
    char outline[1000];
    char mstring[800];
    char *kw, *kwe, *filepath;
    int ikwd, lkwd, nfound, notfound, nch;

    if (rootdir) {
	nch = strlen (rootdir) + strlen (name) + 1;
	filepath = (char *) calloc (1, nch);
	strcat (filepath, rootdir);
	strcat (filepath, "/");
	strcat (filepath, name);
	}
    else
	filepath = name;

    /* Retrieve FITS header from FITS or IRAF .imh file */
    if ((header = GetFITShead (filepath)) == NULL)
	return;

    if (verbose) {
	fprintf (stderr,"Print Header Parameter Values from ");
	hgeti4 (header, "IMHVER", &iraffile );
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", name);
	else
	    fprintf (stderr,"FITS image file %s\n", name);
	}

    /* Find file name */
    if (listpath || (filename = strrchr (name,'/')) == NULL)
	filename = name;
    else if (rootdir)
	filename = name;
    else
	filename = filename + 1;

    if (nfile > 1 || listall) {
	if (tabout)
	    sprintf (fnform, "%%-%ds	", maxlfn);
	else
	    sprintf (fnform, "%%-%ds ", maxlfn);
	sprintf (outline, fnform, filename);
	}
    else
	outline[0] = (char) 0;
    nfound = 0;

    notfound = 0;
    for (ikwd = 0; ikwd < nkwd; ikwd++) {
	lkwd = strlen (kwd[ikwd]);
	kwe = kwd[ikwd] + lkwd;
	for (kw = kwd[ikwd]; kw < kwe; kw++) {
	    if (*kw > 96 && *kw < 123)
		*kw = *kw - 32;
	    }
	if (hgets (header, kwd[ikwd], 80, string)) {
	    strclean (string);
	    if (ndec > -9 && isnum (string) && strchr (string, '.'))
		num2str (string, atof(string), 0, ndec);
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
	else if (printfill)
	    strcat (outline, "___");
	else
	    notfound = 1;

	if (verbose)
	    printf ("\n");
	else if (ikwd < nkwd-1) {
	    if (tabout)
		strcat (outline, "	");
	    else
		strcat (outline, " ");
	    }
	}
    if (!verbose && (nfound > 0 || printfill) && (nfile < 2 || nfound > 0 || listall))
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
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jun  2 1998  Fix bug in hput()
 * Jun  3 1998	Add -p option
 * Jun 18 1998	Print tab table heading only if -h option used
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Sep  1 1998	Set number of decimal places in floating output with -n
 * Oct  5 1998	Check more carefully for FITS and IRAF files on command line
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Feb 12 1999	Print null string if single keyword is not found
 * Feb 17 1999	Add -d option to set root input directory
 * Apr  1 1999	Add warning if too many files or keywords on command line
 * Apr  2 1999	Add -f and -m to change maximum number of files or keywords
 * Jun  9 1999	Initialize outline so Linux doesn't print garbage
 * Jun 21 1999	Fix bug so that -a option works, always printing filename
 */
