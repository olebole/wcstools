/* File gettab.c
 * October 22, 1999
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
#include "libwcs/wcs.h"
#include "libwcs/wcscat.h"

#define MAXCOL 100
#define MAXCOND 10
#define MAXFILES 1000
#define MAXLINES 1000

static void usage();
static void PrintValues();
static char *strclean();

static int verbose = 0;		/* verbose/debugging flag */
static int nfile = 0;
static int ndec = -9;
static int maxlfn = 0;
static int listall = 0;
static int listpath = 0;
static int tabout = 0;
static int printhead = 0;
static int assign = 0;
static int version = 0;		/* If 1, print only program name and version */
static int nlines = 0;
static int *keeplines;		/* List of lines to keep if nlines > 0 */

main (ac, av)
int ac;
char **av;
{
    char *str;
    char *kwd[MAXCOL];
    char *alias[MAXCOL];
    int ifile;
    int nkwd = 0;
    char *fn[MAXFILES];
    int readlist = 0;
    int lfn;
    char *lastchar;
    char filename[128];
    char *name;
    FILE *flist;
    char *listfile;
    int ikwd, lkwd, i;
    int lfield;
    char *kw1;
    char string[80];
    int icond,ncond;
    char *vali, *calias;
    char *cond[MAXCOND];
    char opcond[MAXCOND]; /* Logical conditions <, =, > */
    char *vcond[MAXCOND]; /* ASCII values to be tested against */
    double xcond[MAXCOND]; /* Numeric values to be tested against */
    struct TabTable *tabtable;
    char *ranges = NULL;
    char *temp;
    int nldef = 1;
    struct Range *range; /* Range of sequence numbers to list */

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    nkwd = 0;
    ncond = 0;
    nfile = 0;

    /* crack arguments */
    for (av++; --ac > 0; av++) {
	if (*(str = *av)=='-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

		case 'a':	/* list file even if no keywords are found */
		    listall++;
		    break;

		case 'e':	/* output assignments */
		    assign++;
		    break;

		case 'v':	/* more verbosity */
		    verbose++;
		    break;

		case 'h':	/* output column headings */
		    printhead++;
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

		default:
		    usage();
		    break;
		}
	    continue;
	    }

	/* Set range and make a list of line numbers from it */
	else if (strchr (*av + 1, '-') || strchr (*av + 1, ',')) {
	    if (ranges) {
		temp = ranges;
		ranges = (char *) calloc (strlen(ranges) + strlen(*av) + 2, 1);
		strcpy (ranges, temp);
		strcat (ranges, ",");
		strcat (ranges, *av);
		free (temp);
		}
	    else {
		ranges = (char *) calloc (strlen(*av) + 1, 1);
		strcpy (ranges, *av);
		}
	    continue;
	    }

	/* If numeric argument, set line to be read */
	else if (isnum (str)) {
	    if (ranges) {
		temp = ranges;
		ranges = (char *)calloc (strlen(ranges)+strlen(*av)+2, 1);
		strcpy (ranges, temp);
		strcat (ranges, ",");
		strcat (ranges, *av);
		free (temp);
		}
	    else {
		ranges = (char *) calloc (strlen(*av) + 1, 1);
		strcpy (ranges, *av);
		}
	    continue;
	    }

	else if (*av[0] == '@') {
	    readlist++;
	    listfile = *av + 1;
	    nfile = 2;
	    continue;
	    }
	else if (!strcmp (*av, "stdin")) {
	    fn[nfile] = *av;
	    name = fn[nfile];
	    lfn = strlen (name);
	    if (lfn > maxlfn)
		maxlfn = lfn;
	    nfile++;
	    continue;
	    }

	/* Record condition */
	else if (strchr (*av, '=') || strchr (*av, '<') || strchr (*av, '>')) {
	    cond[ncond] = *av;
	    if (ncond < MAXCOND)
		ncond++;
	    continue;
	    }

	/* Record column name with output alias */
	else if (strchr (*av, '@')) {
	    kwd[nkwd] = *av;
	    calias = strchr (*av, '@');
	    alias[nkwd] = calias + 1;
	    *calias = (char) 0;
	    if (nkwd < MAXCOL)
		nkwd++;
	    continue;
	    }

	/* Record tab table file name */
	else if (istab (*av)) {
	    fn[nfile] = *av;

	    if (listpath || (name = strrchr (fn[nfile],'/')) == NULL)
		name = fn[nfile];
	    else
		name = name + 1;
	    lfn = strlen (name);
	    if (lfn > maxlfn)
		maxlfn = lfn;
	    if (nfile < MAXFILES)
		nfile++;
	    continue;
	    }

	/* Record column name */
	else {
	    kwd[nkwd] = *av;
	    alias[nkwd] = NULL;
	    if (nkwd < MAXCOL)
		nkwd++;
	    continue;
	    }
	}

    /* Read from standard input if no file is specified */
    if (nfile == 0 && (nkwd > 0 || ncond > 0)) {
	name = malloc (8);
	strcpy (name, "stdin");
	fn[nfile] = name;
	lfn = strlen (name);
	if (lfn > maxlfn)
	    maxlfn = lfn;
	nfile++;
	}

    /* Set conditions */
    for (icond = 0; icond < ncond; icond++) {
	vali = strchr (cond[icond], '=');
	if (vali != NULL) {
	    *vali = (char) 0;
	    vcond[icond] = vali + 1;
	    if (!isnum (vcond[icond]) && !strchr (vcond[icond],':'))
		opcond[icond] = 's';
	    else {
		opcond[icond] = '=';
		if (strchr (vcond[icond], ':'))
		    xcond[icond] = str2dec (vcond[icond]);
		else
		    xcond[icond] = atof (vcond[icond]);
		}
	    }
	vali = strchr (cond[icond], '<');
	if (vali != NULL) {
	    opcond[icond] = '<';
	    *vali = (char) 0;
	    vcond[icond] = vali + 1;
	    if (strchr (vcond[icond], ':')) {
		if (strsrch (cond[icond], "RA") || strsrch (cond[icond], "ra"))
		    xcond[icond] = str2ra (vcond[icond]);
		else
		    xcond[icond] = str2dec (vcond[icond]);
		}
	    else
		xcond[icond] = atof (vcond[icond]);
	    }
	vali = strchr (cond[icond], '>');
	if (vali != NULL) {
	    opcond[icond] = '>';
	    *vali = (char) 0;
	    vcond[icond] = vali + 1;
	    if (strchr (vcond[icond], ':'))
		xcond[icond] = str2dec (vcond[icond]);
	    else
		xcond[icond] = atof (vcond[icond]);
	    }
	}

    /* Decode ranges */
    if (ranges != NULL) {
	range = RangeInit (ranges, nldef);
	nlines = rgetn (range);
	keeplines = (int *) calloc (1, nlines);
	for (i = 0; i < nlines; i++)
	    keeplines[i] = rgeti4 (range);
	}

    if (nkwd > 0) {

	/* Print column headings if tab table or headings requested */
	if (printhead || tabout) {

	    /* For tab table output, keep input header information */
	    if (tabout) {
    		if ((tabtable = tabopen (name)) != NULL &&
		    tabtable->tabbuff != tabtable->tabhead) {
		    *(tabtable->tabhead-1) = (char) 0;
		    printf ("%s\n", tabtable->tabbuff);
		    }
		}

	    /* If multiple files, add filname at start of output line */
	    if (nfile > 1) {
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

	    /* Print column names */
	    for (ikwd = 0; ikwd < nkwd; ikwd++) {
		if (alias[ikwd])
		    kw1 = alias[ikwd];
		else
		    kw1 = kwd[ikwd];
		printf ("%s",kw1);
		if ((i = tabcol (tabtable, kwd[ikwd])) > 0) {
		    lfield = tabtable->lcfld[i-1];
		    if (lfield > 32)
			lfield = 32;
		    }
		if (alias[ikwd])
		    lkwd = strlen (alias[ikwd]);
		else
		    lkwd = strlen (kwd[ikwd]);
		if (tabout && lfield > lkwd) {
		    for (i = lkwd; i < lfield; i++)
			printf (" ");
		    }
		if (verbose || ikwd == nkwd - 1)
	    	    printf ("\n");
		else if (tabout)
	    	    printf ("	");
		else
		    printf (" ");
		}

	    /* Print field-defining hyphens if tab table output requested */
	    if (nfile > 1) {
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
		strcpy (string, "---------------------------------");
		if (tabout) {
		    if ((i = tabcol (tabtable, kwd[ikwd])) > 0) {
			lkwd = tabtable->lcfld[i-1];
			if (lkwd > 32)
			    lkwd = 32;
			}
		    else if (alias[ikwd])
			lkwd = strlen (alias[ikwd]);
		    else
			lkwd = strlen (kwd[ikwd]);
		    }
		else if (alias[ikwd])
		    lkwd = strlen (alias[ikwd]);
		else
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
	    if (tabout)
		 tabclose (tabtable);
	    }

    /* Get table values one at a time */

	/* Read through tables in listfile */
	if (readlist) {
	    if ((flist = fopen (listfile, "r")) == NULL) {
		fprintf (stderr,"GETHEAD: List file %s cannot be read\n",
		     listfile);
		usage ();
		}
	    while (fgets (filename, 128, flist) != NULL) {
		lastchar = filename + strlen (filename) - 1;
		if (*lastchar < 32) *lastchar = 0;
		PrintValues (filename,nkwd,kwd,alias,ncond,cond,opcond,vcond,xcond);
		if (verbose)
		    printf ("\n");
		}
	    fclose (flist);
	    }

	/* Read tables from command line list */
	else {
	    for (ifile = 0; ifile < nfile; ifile++)
	    	PrintValues (fn[ifile],nkwd,kwd,alias,ncond,cond,opcond,vcond,xcond);
	    }
	}

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Print FITS or IRAF header keyword values\n");
    fprintf(stderr,"usage: gethead [-ahtv][-n num] file1.tab ... filen.tab kw1 kw2 ... kwn\n");
    fprintf(stderr,"       gethead [-ahptv][-n num] @filelist kw1 kw2 ... kwn\n");
    fprintf(stderr,"       gethead [-ahptv][-n num] <file1.tab kw1 kw2 ... kwn\n");
    fprintf(stderr,"  -a: List file even if keywords are not found\n");
    fprintf(stderr,"  -e: Print keyword=value list\n");
    fprintf(stderr,"  -h: Print column headings\n");
    fprintf(stderr,"  -n: Number of decimal places in numeric output\n");
    fprintf(stderr,"  -p: Print full pathnames of files\n");
    fprintf(stderr,"  -t: Output in tab-separated table format\n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}


static void
PrintValues (name, nkwd, kwd, alias, ncond, kcond, opcond, vcond, xcond)

char	*name;	  /* Name of FITS or IRAF image file */
int	nkwd;	  /* Number of keywords for which to print values */
char	*kwd[];	  /* Names of keywords for which to print values */
char	*alias[]; /* Output names of keywords if different from input */
int	ncond;	  /* Number of conditions which must be met */
char	*kcond[]; /* Columns to be tested */
char	*opcond;  /* Logical conditions <, =, > */
char	*vcond[]; /* Values to be tested against */
double	*xcond;	  /* Values to be tested against */

{
    char *str;
    int drop;
    int icond;
    char fnform[8];
    char string[80];
    char *filename;
    char outline[1000];
    char *line;
    char newline = 10;
    int ikwd, nfound;
    int i, iline, keep;
    double xnum;
    struct TabTable *tabtable;
    int col[MAXCOL];
    int ccond[MAXCOND];

    /* Figure out conditions first, separating out keywords to check */

    /* Retrieve FITS header from FITS or IRAF .imh file */
    if ((tabtable = tabopen (name)) == NULL)
	return;

    if (verbose) {
	fprintf (stderr,"Print table Values from tab table file %s\n", name);
	}

    /* Find file name */
    if (listpath || (filename = strrchr (name,'/')) == NULL)
	filename = name;
    else
	filename = filename + 1;

    if (nfile > 1) {
	if (tabout)
	    sprintf (fnform, "%%-%ds	", maxlfn);
	else
	    sprintf (fnform, "%%-%ds ", maxlfn);
	sprintf (outline, fnform, filename);
	}
    nfound = 0;
    line = tabtable->tabdata;

    /* Find column numbers for column names to speed up extraction */
    ccond[0] = 0;
    for (icond = 0; icond < ncond; icond++)
	ccond[icond] = tabcol (tabtable, kcond[icond]);
    col[0] = 0;
    for (ikwd = 0; ikwd < nkwd; ikwd++)
	col[ikwd] = tabcol (tabtable, kwd[ikwd]);

    iline = 0;
    while (line != NULL && strlen (line) > 0) {
	outline[0] = (char) 0;

	/* Check line number if extracting specific lines */
	iline++;
	drop = 0;
	if (nlines > 0) {
	    keep = 0;
	    for (i = 0; i < nlines; i++) {
		if (iline == keeplines[i])
		    keep = 1;
		}
	    if (!keep)
		drop = 1;
	    }

	/* Check conditions */
	for (icond = 0; icond < ncond; icond++) {
	    if (ccond[icond] > 0 &&
		!tabgetc (tabtable, line, ccond[icond], string, 80)) {
		str = strclean (string);
		if (opcond[icond] == 's') {
		    if (!strcmp (str, vcond[icond])) {
			drop++;
			break;
			}
		    }
		else {
		    if (strchr (str, ':')) {
			if (strsrch (kcond[icond], "RA") ||
			    strsrch (kcond[icond], "ra"))
			    xnum = str2ra (string);
			else
			    xnum = str2dec (string);
			}
		    else if (!isnum (str)) {
			drop++;
			break;
			}
		    else
			xnum = atof (str);
		    if (opcond[icond] == '=' && xnum != xcond[icond]) {
			drop++;
			break;
			}
		    else if (opcond[icond] == '>' && xnum <= xcond[icond]) {
			drop++;
			break;
			}
		    else if (opcond[icond] == '<' && xnum >= xcond[icond]) {
			drop++;
			break;
			}
		    
		    }
		}
	    }

	/* Extract desired columns */
	if (!drop) {
	    for (ikwd = 0; ikwd < nkwd; ikwd++) {
		if (!tabgetc (tabtable, line, col[ikwd], string, 80)) {
		    str = strclean (string);
		    if (ndec > -9 && isnum (str) && strchr (str, '.'))
			num2str (str, atof(str), 0, ndec);
		    if (verbose) {
			if (alias[ikwd])
			    printf ("%s/%s = %s", kwd[ikwd], alias[ikwd], str);
			else
			    printf ("%s = %s", kwd[ikwd], str);
			}
		    else if (assign) {
			if (alias[ikwd])
			    strcat (outline, alias[ikwd]);
			else
			    strcat (outline, kwd[ikwd]);
			strcat (outline, "=");
			strcat (outline, str);
			}
		    else
			strcat (outline, str);
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
	    }

	line = strchr (line+1, newline);
	if (line != NULL)
	    line++;
	}

    tabclose (tabtable);
    return;
}


/* Remove exponent and trailing zeroes, if reasonable */

static char *
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

    /* Remove leading spaces */
    while (*string == ' ')
	string++;

    /* Remove trailing spaces */
    lstr = strlen (string);
    s = string + lstr - 1;
    while (*s == ' ') {
	*s = (char) 0;
	s--;
	}

    return (string);
}

/* Jan 22 1999	New program
 * Jan 25 1999	Keep header information
 * Mar  9 1999	Add range of lines; rework command line decoding logic
 * Oct 22 1999	Drop unused variables after lint
 */
