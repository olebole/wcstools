/* File getcol.c
 * March 20, 2000
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include "libwcs/wcscat.h"

static void usage();

static int maxncond = 100;
static int verbose = 0;		/* Verbose/debugging flag */
static int wfile = 0;		/* True to print output file */
static int debug = 0;		/* True for extra information */
static int sumcol = 0;		/* True to sum columns */
static int meancol = 0;		/* True to compute mean of columns */
static int countcol = 0;	/* True to count entries in columns */
static char *objname = NULL;	/* Object name for output */
static char *keyword = NULL;	/* Column to add to tab table output */
static char progpath[128];
static int ncat = 0;
static int version = 0;		/* If 1, print only program name and version */
static int nread = 0;		/* Number of lines to read (0=all) */
static int nskip = 0;		/* Number of lines to skip */
static int tabout = 0;		/* If 1, separate output fields with tabs */
static int counttok = 0;	/* If 1, print number of columns on line */
static int printhead = 0;	/* If 1, print Starbase tab table header */
static int intcompare();
static int ncond=0;		/* Number of keyword conditions to check */
static int condand=1;		/* If 1, AND comparisons, else OR */
static char **cond;		/* Conditions to check */
static char **ccond;		/* Condition characters */
static void strclean();

main (ac, av)
int ac;
char **av;
{
    char *str;
    char *temp;
    char *filename;
    int systemp = 0;		/* Input search coordinate system */
    char *ranges = NULL;
    char *lfile = NULL;
    char *lranges = NULL;
    int match, newbytes, nrbytes;
    int icond;

    if (ac == 1)
        usage ();

    cond = (char **)calloc (maxncond, sizeof(char *));
    ccond = (char **)calloc (maxncond, sizeof(char *));

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage ();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage ();
	}

    /* crack arguments */
    for (av++; --ac > 0; av++) {

	/* Set column number */
	if (isnum (*av)) {
	    if (strchr (*av,'.'))
		match = 1;
	    if (ranges) {
		newbytes = strlen(ranges) + strlen(*av) + 2;
		newbytes = ((newbytes / 16) + 1) * 16;
		if (newbytes > nrbytes) {
		    temp = ranges;
		    ranges = (char *)calloc (newbytes, 1);
		    strcpy (ranges, temp);
		    free (temp);
		    }
		strcat (ranges, ",");
		strcat (ranges, *av);
		}
	    else {
		nrbytes = strlen(*av) + 2;
		nrbytes = ((nrbytes / 16) + 1) * 16;
		ranges = (char *) calloc (nrbytes, 1);
		strcpy (ranges, *av);
		}
	    }

	/* Set range and make a list of column numbers from it */
	else if (isrange (*av)) {
	    if (ranges) {
		newbytes = strlen(ranges) + strlen(*av) + 2;
		newbytes = ((newbytes / 16) + 1) * 16;
		if (newbytes > nrbytes) {
		    temp = ranges;
		    ranges = (char *) calloc (newbytes, 1);
		    strcpy (ranges, temp);
		    free (temp);
		    }
		strcat (ranges, ",");
		strcat (ranges, *av);
		}
	    else {
		ranges = (char *) calloc (strlen(*av) + 2, 1);
		if (strchr (*av,'.'))
		    match = 1;
		strcpy (ranges, *av);
		}
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
	    ncond++;
	    }

	/* Otherwise, read command */
	else if (*(str = *av) == '-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

	    case 'a':	/* Sum values in each numeric column */
		sumcol++;
		break;

	    case 'c':	/* Count entries in each column */
		countcol++;
		break;

	    case 'h':	/* Print Starbase tab table header */
		printhead++;
		break;

	    case 'k':	/* Count columns on first line */
		counttok++;
		break;

	    case 'm':	/* Compute mean of each numeric column */
		meancol++;
		break;

	    case 'n':	/* Number of lines to read */
		if (ac < 2)
		    usage ();
		nread = atoi (*++av);
		ac--;
		break;

	    case 'o': /* OR conditions insted of ANDing them */
		condand = 0;
		break;

	    case 'r':	/* Range of lines to read */
		if (ac < 2)
		    usage ();
		if (*(av+1)[0] == '@') {
		    lfile = *++av + 1;
		    ac--;
		    }
		else if (isrange (*(av+1))) {
		    lranges = (char *) calloc (strlen(*av) + 1, 1);
		    strcpy (lranges, *++av);
		    ac--;
		    }
		break;

	    case 's':	/* Number of lines to skip */
		if (ac < 2)
		    usage ();
		nskip = atoi (*++av);
		ac--;
		break;

	    case 't':	/* Tab-separated output */
		tabout++;
		break;

	    case 'v':	/* More verbosity */
		verbose++;
		break;

	    default:
		usage ();
		break;
	    }
	    }
	else
	    filename = *av;
	}
    ListFile (filename, ranges, lranges, lfile);

    free (lranges);
    free (ranges);
    return (0);
}

static void
usage ()

{
    if (version)
	exit (-1);
    fprintf (stderr,"Extract specified columns from an ASCII table file\n");
    fprintf (stderr,"Usage: [-amv][-n num][-r lines][-s num] filename [column number range]\n");
    fprintf(stderr,"  -a: Sum numeric colmuns\n");
    fprintf(stderr,"  -c: Add count of number of lines in each column at end\n");
    fprintf(stderr,"  -h: Print Starbase tab table header\n");
    fprintf(stderr,"  -k: Print number of columns on first line\n");
    fprintf(stderr,"  -m: Compute mean of numeric columns\n");
    fprintf(stderr,"  -n: Number of lines to read, if not all\n");
    fprintf(stderr,"  -r: Range or @file of lines to read, if not all\n");
    fprintf(stderr,"  -s: Number of lines to skip\n");
    fprintf(stderr,"  -t: Starbase tab table output\n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}

static int
ListFile (filename, ranges, lranges, lfile)

char	*filename;	/* File name */
char	*ranges;	/* String with range of column numbers to list */
char	*lranges;	/* String with range of lines to list */
char	*lfile;		/* Name of file with lines to list */

{
    int i, j, il, nbytes;
    char line[1024];
    char fline[1024];
    char *lastchar;
    FILE *fd;
    FILE *lfd;
    int nlog;
    struct Tokens tokens;  /* Token structure */
    struct Range *range;
    struct Range *lrange;
    int *iline;
    int nline;
    int iln;
    int nfdef = 9;
    double *sum;
    int *nsum;
    int *nent;
    int nlmax;
    int nfind, ntok, nt, ltok;
    int *inum;
    int skipline, icond, itok;
    char tcond, *cstr, *cval;
    char numstr[32], numstr1[32];
    double dcond, dval, dnum;
    int pass;
    int jcond, jval;
    char *cwhite;
    char token[80];

    cwhite = NULL;
    lrange = NULL;
    iline = NULL;

    if (verbose)

    if (debug)
	nlog = 1;
    else if (verbose)
	nlog = 100;
    else
	nlog = 0;
    if (nread < 1)
	nread = 100000;

    /* Make list of line numbers to read from list or range on command line */
    if (lranges != NULL) {
	lrange = RangeInit (lranges, nfdef);
	nline = rgetn (lrange);
	nbytes = nline * sizeof (int);
	if (!(iline = (int *) calloc (nline, sizeof(int))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for iline\n", nbytes);
	    return;
	    }
	for (i = 0; i < nline; i++)
	    iline[i] = rgeti4 (lrange);
	qsort (iline, nline, sizeof(int), intcompare);
	}

    /* Make list of line numbers to read from file specified on command line */
    if (lfile != NULL) {
	if (!(lfd = fopen (lfile, "r")))
            return (0);
	nlmax = 100;
	nline = 0;
	nbytes = nlmax * sizeof(int);
	if (!(iline = (int *) calloc (nlmax, sizeof(int))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for iline\n", nbytes);
	    fclose (lfd);
	    return;
	    }
	for (il = 0; il < nread; il++) {
	    if (fgets (line, 1024, lfd) == NULL)
		break;

	    /* Drop linefeeds */
	    lastchar = line + strlen(line) - 1;
	    if (*lastchar < 32)
		*lastchar = (char) 0;

	    ntok = setoken (&tokens, line, cwhite);
	    nt = 0;
	    if (il > nlmax) {
		nlmax = nlmax + 100;
		nbytes = nlmax * sizeof(int);
		if (!(iline = (int *) realloc (iline, nbytes))) {
		    fprintf (stderr, "Could not realloc %d bytes for iline\n",
			     nbytes);
		    fclose (lfd);
		    return;
		    }
		}
	    for (i = 0; i < ntok; i++) {
		if (getoken (tokens, i+1, token)) {
		    iline[il] = atoi (token);
		    if (iline[il] > 0) {
			nline++;
			break;
			}
		    }
		}
	    }
	fclose (lfd);
	qsort (iline, nline, sizeof(int), intcompare);
	}

    /* Open input file */
    if (!strcmp (filename, "stdin"))
	fd = stdin;
    else if (!(fd = fopen (filename, "r"))) {
        return (0);
	}

    /* Skip lines into input file */
    if (nskip > 0) {
	for (i = 0; i < nskip; i++) {
	    if (fgets (line, 80, fd) == NULL)
		break;
	    }
	}

    /* Print entire selected lines */
    if (ranges == NULL) {
	iln = 0;
	for (il = 0; il < nread; il++) {
	    if (fgets (line, 1024, fd) == NULL)
		break;

	    /* Skip if line is not on list, if there is one */
	    if (iline != NULL) {
		if (il+1 < iline[iln])
		    continue;
		else if (il+1 > iline[nline-1])
		    break;
		else
		    iln++;
		}

	    /* Turn linefeed into end of string */
	    lastchar = line + strlen(line) - 1;
	    if (*lastchar < 32)
		*lastchar = (char) 0;

	    /* Check conditions */
	    ntok = setoken (&tokens, line, cwhite);
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
		    itok = atoi (cond[icond]);
		    getoken (tokens, itok, token);
		    cval = token;
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
		    }
		if (condand && !pass)
		    break;
		}
	    if (pass)
		printf ("%s\n", line);
	    }
	}

    /* Find columns specified by number */
    else {
	range = RangeInit (ranges, nfdef);
	nfind = rgetn (range);
	nbytes = nfind * sizeof (double);
	if (!(sum = (double *) calloc (nfind, sizeof (double))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for sum\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return;
	    }
	nbytes = nfind * sizeof (int);
	if (!(nsum = (int *) calloc (nfind, sizeof(int))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for nsum\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return;
	    }
	if (!(nent = (int *) calloc (nfind, sizeof(int))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for nent\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return;
	    }
	if (!(inum = (int *) calloc (nfind, sizeof(int))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for inum\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return;
	    }
	else {
	    for (i = 0; i < nfind; i++)
		inum[i] = rgeti4 (range);
	    }
	wfile = 0;
	iln = 0;
	for (il = 0; il < nread; il++) {
	    if (fgets (line, 1024, fd) == NULL)
		break;

	    /* Skip if line is not on list, if there is one */
	    if (iline != NULL) {
		if (il+1 < iline[iln])
		    continue;
		else if (il > iline[nline-1])
		    break;
		else
		    iln++;
		}

	    /* Turn linefeed into end of string */
	    lastchar = line + strlen(line) - 1;
	    if (*lastchar < 32)
		*lastchar = (char) 0;

	    ntok = setoken (&tokens, line, cwhite);
	    if (counttok) {
		printf ("%d", ntok);
		if (verbose)
		    printf (" columns in %s", filename);
		else
		    printf ("\n");
		return (0);
		}

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
		    itok = atoi (cond[icond]);
		    getoken (tokens, itok, token);
		    cval = token;
		    if (strchr (cval, ':')) {
			dnum = str2dec (cval);
			num2str (numstr, dnum, 0, 7);
			cval = numstr;
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
			break;
		    }
		if (!pass)
		    continue;
		}

	    nt = 0;
	    if (il == 0 && printhead && tabout) {

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
			    printf ("condition  %s\n", cond[icond]);
			else
			    printf ("condition  or %s\n", cond[icond]);
			}
		    }

		for (i = 0; i < nfind; i++) {
		    if (getoken (tokens, inum[i], token)) {
			ltok = strlen (token);
			printf ("%03d", inum[i]);
			for (j = 3; j < ltok; j++)
			    printf (" ");
			}
		    else
			printf ("%03d", inum[i]);
		    printf ("	");
		    }
		printf ("\n");
		for (i = 0; i < nfind; i++) {
		    if (getoken (tokens, inum[i], token)) {
			ltok = strlen (token);
			for (j = 0; j < ltok; j++)
			    printf ("-");
			}
		    else
			printf ("---");
		    printf ("	");
		    }
		printf ("\n");
		}
	    for (i = 0; i < nfind; i++) {
		if (getoken (tokens, inum[i], token)) {
		    if (inum[i] > tokens.ntok || inum[i] < 1)
			printf ("___");
		    else
			printf ("%s", token);
		    if (i < nfind-1) {
			if (tabout)
			    printf ("	");
			else
			    printf (" ");
			}
		    nt++;
		    nent[i]++;
		    if (isnum (token)) {
			sum[i] = sum[i] + atof (token);
			nsum[i]++;
			}
		    }
		}
	    if (nt > 0)
		printf ("\n");
	    }
        }

    /* Print sums of numeric columns */
    if (sumcol) {
	for (i = 0; i < nfind; i++) {
	    if (nsum[i] > 0) {
		if (i < nfind-1)
		    printf ("%f ", sum[i]);
		else
		    printf ("%f", sum[i]);
		}
	    else if (i < nfind-1)
		printf ("___ ");
	    else
		printf ("___");
	    }
	if (nfind > 0)
	    printf ("\n");
	}

    /* Print means of numeric columns */
    if (meancol) {
	for (i = 0; i < nfind; i++) {
	    if (nsum[i] > 0) {
		if (i < nfind-1)
		    printf ("%f ", sum[i]/(double)nsum[i]);
		else
		    printf ("%f", sum[i]/(double)nsum[i]);
		}
	    else if (i < nfind-1)
		printf ("___ ");
	    else
		printf ("___");
	    }
	if (nfind > 0)
	    printf ("\n");
	}

    /* Print count for each column */
    if (countcol) {
	for (i = 0; i < nfind; i++) {
	    if (nent[i] > 0) {
		if (i < nfind-1)
		    printf ("%d ", nent[i]);
		else
		    printf ("%d", nent[i]);
		}
	    }
	if (nfind > 0)
	    printf ("\n");
	}

    /* Free memory used for search results */
    if (inum) free ((char *)inum);
    if (nsum) free ((char *)nsum);
    if (nent) free ((char *)nent);
    if (sum) free ((char *)sum);

    if (fd != stdin) fclose (fd);
    return (nfind);
}

static int
intcompare (int *i, int *j)
{
    if (*i > *j)
	return (1);
    if (*i < *j)
	return (-1);
    return (0);
}


/* Remove exponent and/or trailing zeroes, if reasonable */
static void
strclean (string)

char *string;

{
    char *sdot, *s, *strend, *str;
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


/* Nov  2 1999	New program
 * Nov  3 1999	Add option to read from stdin as input filename
 * Dec  1 1999	Add options to print counts, means, and sums of columns
 * Dec 14 1999	Add option for tab-separated output
 *
 * Jan  7 2000	Add option to list range of lines or filed list of lines
 * Jan 26 2000	Add documentation of entry count and tab output
 * Jan 26 2000	Add option to print tab table header
 * Feb 11 2000	Fix reallocation of range variables
 * Mar 20 2000	Add conditional line printing
 */
