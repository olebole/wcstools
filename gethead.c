/* File gethead.c
 * February 21, 2001
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
#include "libwcs/fitsfile.h"
#include "libwcs/wcscat.h"

#define MAXKWD 100
#define MAXFILES 2000
static int maxnkwd = MAXKWD;
static int maxncond = MAXKWD;
static int maxnfile = MAXFILES;

#define FILE_FITS 1
#define FILE_IRAF 2
#define FILE_ASCII 3

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
static int printfile=1;		/* If 1, print filename first if >1 files */
static int fillblank=0;		/* If 1, replace blanks in strings with _ */
static int keyeqval=0;		/* If 1, print keyword=value, not just value */
static char *rootdir=NULL;	/* Root directory for input files */
static int ncond=0;		/* Number of keyword conditions to check */
static int condand=1;		/* If 1, AND comparisons, else OR */
static char **cond;		/* Conditions to check */
static char **ccond;		/* Condition characters */

main (ac, av)
int ac;
char **av;
{
    char *str;
    char **kwd;		/* Keywords to read */
    int nkwd = 0;
    char **fn;
    int *ft;
    int ifile;
    int lfn;
    char filename[256];
    char *name;
    FILE *flist, *fdk;
    char *listfile;
    char *ilistfile;
    char *klistfile;
    int ikwd, lkwd, i;
    char *kw, *kwe;
    char string[80];
    int nbytes;
    int filetype;
    int icond;

    ilistfile = NULL;
    klistfile = NULL;
    nkwd = 0;
    ncond = 0;
    nfile = 0;
    fn = (char **)calloc (maxnfile, sizeof(char *));
    ft = (int *)calloc (maxnfile, sizeof(int));
    kwd = (char **)calloc (maxnkwd, sizeof(char *));
    cond = (char **)calloc (maxncond, sizeof(char *));
    ccond = (char **)calloc (maxncond, sizeof(char *));

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    /* crack arguments */
    for (av++; --ac > 0; av++) {
	if (*(str = *av)=='-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

		case 'a': /* list file even if the keywords are not found */
		    listall++;
		    break;
	
		case 'b': /* Replace blanks with underscores */
		    fillblank++;
		    break;
	
		case 'd': /* Root directory for input */
		    if (ac < 2)
			usage();
		    rootdir = *++av;
		    ac--;
		    break;

		case 'e': /* list keyword=value for input to sethead */
		    keyeqval++;
		    break;
	
		case 'f': /* Do not print file names */
		    printfile = 0;
		    break;
	
		case 'h': /* Output column headings */
		    printhead++;
		    break;
	
		case 'n': /* Number of decimal places in output */
		    if (ac < 2)
			usage();
		    ndec = (int) (atof (*++av));
		    ac--;
		    break;
	
		case 'o': /* OR conditions insted of ANDing them */
		    condand = 0;
		    break;
	
		case 'p': /* List file pathnames, not just file names */
		    listall++;
		    listpath++;
		    break;
	
		case 't': /* Output tab table */
		    tabout++;
		    break;

		case 'u': /* Always print ___ if keyword not found */
		    printfill++;
		    break;
	
		case 'v': /* More verbosity */
		    verbose++;
		    break;

		default:
		    usage();
		    break;
		}
	    }

	/* File containing a list of keywords or files */
	else if (*av[0] == '@') {
	    listfile = *av + 1;
	    if (isimlist (listfile)) {
		ilistfile = listfile;
		nfile = getfilelines (ilistfile);
		}
	    else {
		klistfile = listfile;
		nkwd = getfilelines (klistfile);
		if (nkwd > 0) {
		    if (nkwd > maxnkwd) {
			kwd = (char **) realloc ((void *)kwd, nkwd);
			maxnkwd = nkwd;
			}
		    if ((fdk = fopen (klistfile, "r")) == NULL) {
			fprintf (stderr,"GETHEAD: File %s cannot be read\n",
				 klistfile);
			nkwd = 0;
			}
		    else {
			for (ikwd = 0; ikwd < nkwd; ikwd++) {
			    kwd[ikwd] = (char *) calloc (32, 1);
			    first_token (fdk, 31, kwd[ikwd]);
			    }
			fclose (fdk);
			}
		    }
		}
	    }

	/* Image file */
	else if (isfits (*av) || isiraf (*av)) {
	    if (nfile >= maxnfile) {
		maxnfile = maxnfile * 2;
		nbytes = maxnfile * sizeof (char *);
		fn = (char **) realloc ((void *)fn, nbytes);
		nbytes = maxnfile * sizeof (int);
		ft = (int *) realloc ((void *)ft, nbytes);
		}
	    fn[nfile] = *av;
	    if (isfits (*av))
		ft[nfile] = FILE_FITS;
	    else
		ft[nfile] = FILE_IRAF;

	    if (listpath || (name = strrchr (fn[nfile],'/')) == NULL)
		name = fn[nfile];
	    else
		name = name + 1;
	    lfn = strlen (name);
	    if (lfn > maxlfn)
		maxlfn = lfn;
	    nfile++;
	    }

	/* Image file */
	else if (isfile (*av)) {
	    if (nfile >= maxnfile) {
		maxnfile = maxnfile * 2;
		nbytes = maxnfile * sizeof (char *);
		fn = (char **) realloc ((void *)fn, nbytes);
		nbytes = maxnfile * sizeof (int);
		ft = (int *) realloc ((void *)ft, nbytes);
		}
	    fn[nfile] = *av;
	    ft[nfile] = FILE_ASCII;

	    if (listpath || (name = strrchr (fn[nfile],'/')) == NULL)
		name = fn[nfile];
	    else
		name = name + 1;
	    lfn = strlen (name);
	    if (lfn > maxlfn)
		maxlfn = lfn;
	    nfile++;
	    }

	/* Condition */
	else if (strchr (*av, '=') != NULL || strchr (*av, '#') != NULL ||
		 strchr (*av, '>') != NULL || strchr (*av, '<') != NULL ) {
	    if (ncond >= maxncond) {
		maxncond = maxncond * 2;
		cond = (char **)realloc((void *)cond, maxncond*sizeof(void *));
		ccond = (char **)realloc((void *)cond, maxncond*sizeof(void *));
		}
	    cond[ncond] = *av;
	    ccond[ncond] = strchr (*av, '=');
	    if (ccond[ncond] == NULL)
		ccond[ncond] = strchr (*av, '#');
	    if (ccond[ncond] == NULL)
		ccond[ncond] = strchr (*av, '>');
	    if (ccond[ncond] == NULL)
		ccond[ncond] = strchr (*av, '<');
	    kwe = ccond[ncond];
	    if (kwe != NULL) {
		for (kw = cond[ncond]; kw < kwe; kw++) {
		    if (*kw > 96 && *kw < 123)
			*kw = *kw - 32;
		    }
		}
	    ncond++;
	    }

	/* Keyword */
	else {
	    if (nkwd >= maxnkwd) {
		maxnkwd = maxnkwd * 2;
		kwd = (char **) realloc ((void *)kwd, maxnkwd);
		}
	    kwd[nkwd] = *av;
	    nkwd++;
	    }
	}

    if (nkwd <= 0 && nfile <= 0 )
	usage ();
    else if (nkwd <= 0) {
	fprintf (stderr, "GETHEAD: no keywords specified\n");
	exit (1);
	}
    else if (nfile <= 0 ) {
	fprintf (stderr, "GETHEAD: no files specified\n");
	exit (1);
	}

    if (nkwd > 1)
	printfill = 1;
    if (nfile < 2 && !listall)
	printfile = 0;

    /* Print column headings if tab table or headings requested */
    if (printhead) {

	/* Print conditions in header */
	for (icond = 0; icond < ncond; icond++) {
	    if (verbose) {
		if (condand || icond == 0)
		    printf ("%s\n",cond[icond]);
		else
		    printf (" or %s\n",cond[icond]);
		}
	    else if (tabout) {
		if (condand || icond == 0)
		    printf ("condition	%s\n", cond[icond]);
		else
		    printf ("condition	or %s\n", cond[icond]);
		}
	    }

	if (printfile) {
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

	/* Make keyword names upper case and print keyword names in header */
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
	if (printfile) {
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

    /* Open file containing a list of images, if there is one */
    if (ilistfile != NULL) {
	if ((flist = fopen (ilistfile, "r")) == NULL) {
	    fprintf (stderr,"GETHEAD: Image list file %s cannot be read\n",
		     ilistfile);
	    usage ();
	    }
	}

    /* Read through headers of images */
    for (ifile = 0; ifile < nfile; ifile++) {
	if (ilistfile != NULL) {
	    first_token (flist, 254, filename);
	    if (isiraf (filename))
		filetype = FILE_IRAF;
	    else
		filetype = FILE_FITS;
	    PrintValues (filename, filetype, nkwd, kwd);
	    }
	else
	    PrintValues (fn[ifile], ft[ifile], nkwd, kwd);

	if (verbose)
	    printf ("\n");
	}
    if (ilistfile != NULL)
	fclose (flist);

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Print FITS or IRAF header keyword values\n");
    fprintf(stderr,"usage: gethead [-abhoptv][-d dir][-f num][-m num][-n num] file1.fit ... filen.fits kw1 kw2 ... kwn\n");
    fprintf(stderr,"       gethead [-abhoptv][-d dir][-f num][-m num][-n num] @filelist kw1 kw2 ... kwn\n");
    fprintf(stderr,"usage: gethead [-abhoptv][-d dir][-f num][-m num][-n num] file1.fit ... filen.fits @keywordlist\n");
    fprintf(stderr,"       gethead [-abhoptv][-d dir][-f num][-m num][-n num] @filelist @keywordlist\n");
    fprintf(stderr,"  -a: List file even if keywords are not found\n");
    fprintf(stderr,"  -b: Replace blanks in strings with underscores\n");
    fprintf(stderr,"  -d: Root directory for input files (default is cwd)\n");
    fprintf(stderr,"  -e: Output keyword=value's on one line per file\n");
    fprintf(stderr,"  -f: Never print filenames (default is print if >1)\n");
    fprintf(stderr,"  -h: Print column headings\n");
    fprintf(stderr,"  -n: Number of decimal places in numeric output\n");
    fprintf(stderr,"  -o: OR conditions instead of ANDing them\n");
    fprintf(stderr,"  -p: Print full pathnames of files\n");
    fprintf(stderr,"  -t: Output in tab-separated table format\n");
    fprintf(stderr,"  -u: Always print ___ if keyword not found\n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}


static void
PrintValues (name, filetype, nkwd, kwd)

char	*name;		/* Name of FITS or IRAF image file */
int	filetype;	/* Type of file (FILE_FITS, FILE_IRAF, FILE_ASCII) */
int	nkwd;		/* Number of keywords for which to print values */
char	*kwd[];		/* Names of keywords for which to print values */

{
    char *header;	/* FITS image header or contents of ASCII file */
    char *cstr, *str;
    int iraffile;
    char fnform[8];
    char string[80];
    char temp[1028];
    char keyword[16];
    char *filename;
    char outline[1000];
    char mstring[800];
    char *kw, *kwe, *filepath;
    int ikwd, lkwd, nfound, notfound, nch;
    int jval, jcond, icond, i, lstr;
    double dval, dcond, dnum;
    char tcond;
    char cvalue[64], *cval;
    char numstr[32], numstr1[32];
    int pass, iwcs;

    if (rootdir) {
	nch = strlen (rootdir) + strlen (name) + 1;
	filepath = (char *) calloc (1, nch);
	strcat (filepath, rootdir);
	strcat (filepath, "/");
	strcat (filepath, name);
	}
    else
	filepath = name;

    /* Read ASCII file into buffer */
    if (filetype == FILE_ASCII) {
	if ((header = getfilebuff (filepath)) == NULL)
	    return;
	}

    /* Retrieve FITS header from FITS or IRAF .imh file */
    else if ((header = GetFITShead (filepath)) == NULL)
	return;

    if (verbose) {
	fprintf (stderr,"Print Header Parameter Values from ");
	hgeti4 (header, "IMHVER", &iraffile );
	if (filetype == FILE_ASCII)
	    fprintf (stderr,"ASCII file %s\n", name);
	else if (filetype == FILE_IRAF)
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

    if (printfile) {
	if (tabout)
	    sprintf (fnform, "%%-%ds	", maxlfn);
	else
	    sprintf (fnform, "%%-%ds ", maxlfn);
	sprintf (outline, fnform, filename);
	}
    else
	outline[0] = (char) 0;

    /* Check conditions */
    pass = 0;
    if (ncond > 0) {
	for (icond = 0; icond < ncond; icond++) {
	    if (condand)
		pass = 0;
	    tcond = *ccond[icond];

	    /* Extract test value from comparison string */
	    *ccond[icond] = (char) 0;
	    cstr = ccond[icond]+1;
	    if (strchr (cstr, ':')) {
		dnum = str2dec (cstr);
		num2str (numstr, dnum, 0, 7);
		cstr = numstr;
		}
	    strclean (cstr);

	    /* Read comparison value from header */
	    if (!hgets (header, cond[icond], 64, cvalue))
		continue;
	    cval = cvalue;
	    if (strchr (cval, ':')) {
		dnum = str2dec (cval);
		num2str (numstr1, dnum, 0, 7);
		cval = numstr1;
		}
	    strclean (cval);

	    /* Compare floating point numbers */
	    if (isnum (cstr) == 2 && isnum (cval)) {
		*ccond[icond] = tcond;
		dcond = atof (cstr);
		dval = atof (cval);
		if (tcond == '=' && dval == dcond)
		    pass = 1;
		if (tcond == '#' && dval != dcond)
		    pass = 1;
		if (tcond == '>' && dval > dcond)
		    pass = 1;
		if (tcond == '<' && dval < dcond)
		    pass = 1;
		}

	    /* Compare integers */
	    else if (isnum (cstr) == 1 && isnum (cval)) {
		*ccond[icond] = tcond;
		jcond = atoi (cstr);
		jval = atoi (cval);
		if (tcond == '=' && jval == jcond)
		    pass = 1;
		if (tcond == '#' && jval != jcond)
		    pass = 1;
		if (tcond == '>' && jval > jcond)
		    pass = 1;
		if (tcond == '<' && jval < jcond)
		    pass = 1;
		}

	    /* Compare strings (only equal or not equal */
	    else {
		*ccond[icond] = tcond;
		if (tcond == '=' && !strcmp (cstr, cval))
		    pass = 1;
		if (tcond == '#' && strcmp (cstr, cval))
		    pass = 1;
		}
	    if (condand && !pass)
		return;
	    }
	if (!pass)
	    return;
	}

    /* Read keywords from header */
    nfound = 0;
    notfound = 0;
    for (ikwd = 0; ikwd < nkwd; ikwd++) {

	/* IF FITS header, look for upper case keywords only */
	if (filetype != FILE_ASCII) {
	    lkwd = strlen (kwd[ikwd]);
	    kwe = kwd[ikwd] + lkwd;
	    for (kw = kwd[ikwd]; kw < kwe; kw++) {
		if (*kw > 96 && *kw < 123)
		    *kw = *kw - 32;
		}
	    }
	strcpy (keyword, kwd[ikwd]);

	/* Read keyword value from ASCII file */
	if (filetype == FILE_ASCII &&
	    agets (header, keyword, 80, string)) {
	    str = string;
	    strclean (str);
	    if (ndec > -9 && isnum (str) && strchr (str, '.'))
		num2str (str, atof(str), 0, ndec);
	    if (verbose)
		printf ("%s = %s\n", keyword, str);
	    else if (keyeqval) {
		sprintf (temp, " %s=%s", keyword, str);
		strcat (outline, temp);
		}
	    else
		strcat (outline, str);
	    nfound++;
	    }

	/* Read all keywords from multiple WCS's if wildcarded with @ */
	else if (keyword[lkwd-1] == '@') {
	    keyword[lkwd-1] = (char) 0;
	    for (iwcs = -1; iwcs < 26; iwcs++) {
		if (iwcs < 0)
		    keyword[lkwd-1] = (char)0;
		else
		    keyword[lkwd-1] = (char)(iwcs+65);
		if (hgets (header, keyword, 80, string)) {
		    if (string[0] == '#' && isnum (string+1)) {
			lstr = strlen (string);
			for (i = 0; i < lstr; i++)
			    string[i] = string[i+1];
			}
		    str = string;
		    strclean (str);
		    if (ndec > -9 && isnum (str) && strchr (str, '.'))
			num2str (string, atof(str), 0, ndec);
		    if (verbose)
			printf ("%s = %s\n", keyword, str);
		    else if (keyeqval) {
			sprintf (temp, " %s=%s", keyword, str);
			strcat (outline, temp);
			}
		    else
			strcat (outline, str);
		    nfound++;
		    }
		}
	    }

	/* Read one FITS keyword value */
	else if (hgets (header, keyword, 80, string)) {
	    if (string[0] == '#' && isnum (string+1)) {
		lstr = strlen (string);
		for (i = 0; i < lstr; i++)
		    string[i] = string[i+1];
		}
	    str = string;
	    strclean (str);
	    if (ndec > -9 && isnum (str) && strchr (str, '.'))
		num2str (string, atof(str), 0, ndec);
	    if (verbose)
		printf ("%s = %s", keyword, str);
	    else if (keyeqval) {
		sprintf (temp, " %s=%s", keyword, str);
		strcat (outline, temp);
		}
	    else
		strcat (outline, str);
	    nfound++;
	    }

	/* Read IRAF-style multiple-line keyword value */
	else if (hgetm (header, keyword, 600, mstring)) {
	    if (verbose)
		printf ("%s = %s\n", keyword, mstring);
	    else if (keyeqval) {
		sprintf (temp, " %s=%s", keyword, mstring);
		strcat (outline, temp);
		}
	    else
		strcat (outline, mstring);
	    nfound++;
	    }

	else if (verbose)
	    printf ("%s not found\n", keyword);
	else if (printfill)
	    strcat (outline, "___");
	else
	    notfound = 1;

	if (!verbose && ikwd < nkwd-1) {
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


/* Remove exponent, leading #, and/or trailing zeroes, if reasonable */
static void
strclean (string)

char *string;

{
    char *sdot, *s, *strend, *str, ctemp, *slast;
    int ndek, lstr, i;

    /* If number, ignore leading # and remove trailing non-numeric character */
    if (string[0] == '#') {
	strend = string + strlen (string);
	str = string + 1;
	strend = str + strlen (str) - 1;
	ctemp = *strend;
	if (!isnum (strend))
	    *strend = (char) 0;
	if (isnum (str)) {
	    strend = string + strlen (string);
	    for (str = string; str < strend; str++)
		*str = *(str + 1);
	    }
	else
	    *strend = ctemp;
	}

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

    /* Replace embedded blanks with underscores, if requested to */
    if (fillblank) {
	lstr = strlen (string);
	slast = string + lstr;
	for (s = string; s < slast; s++) {
	    if (*s == ' ') *s = '_';
	    }
	}

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
 * Jul 13 1999	Use only first token from line of list as filename
 * Jul 14 1999	Read lists of BOTH keywords and files simultaneously
 * Jul 15 1999	Reallocate keyword and file lists if default limits exceeded
 * Oct 22 1999	Drop unused variables after lint
 * Nov 16 1999	Add -u to always print underscores if keyword not found
 * Nov 19 1999	Add -f to never print filenames
 * Nov 30 1999	Fix so no file name is printed if only one file unless -a
 * Nov 30 1999	Cast realloc's
 *
 * Feb 24 2000	Add option to output keyword=value
 * Mar  1 2000	Add option to read ASCII files with keyword=value in them
 * Mar 17 2000	Add conditions
 * Mar 20 2000	Drop leading # from numbers
 * Mar 21 2000	Add -b option to replace blanks with underscores
 * Jun  8 2000	If no keywords or files specified, say so
 * Jun 12 2000	If -p is set, print all file names (-a)
 *
 * Feb 21 2001	Add @ wildcard option for multiple WCS keywords
 */
