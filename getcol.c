/* File getcol.c
 * April 12, 2002
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
#include "libwcs/fitshead.h"

#define	MAX_LTOK	256
#define	MAX_NTOK	256

static void usage();
static int ListFile();
static int iscolop();
static int iscol();
static double median();

static double badval = 0.0;	/* Value to ignore */
static int isbadval = 0;	/* 1 if badval is set */
static int maxncond = 100;
static int maxnop = 100;
static int verbose = 0;		/* Verbose/debugging flag */
static int debug = 0;		/* True for extra information */
static int sumcol = 0;		/* True to sum column values */
static int ameancol = 0;	/* True for absolute mean of column values */
static int meancol = 0;		/* True to compute mean of column values */
static int amedcol = 0;		/* True for absolute median of column values */
static int medcol = 0;		/* True to compute median of column values */
static int countcol = 0;	/* True to count entries in columns */
static int rangecol = 0;	/* True to print range of column values */
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
static int nop=0;		/* Number of keyword values to operate on */
static char **op;		/* Operations to check */
static char **cop;		/* Operation characters */
static void strclean();
static int napp=0;		/* Number of lines to append */
static int ndec=-1;		/* Number of decimal places in f.p. output */
static int printcol = 1;	/* Flag to print extracted data */
static char *cwhite;
static int *frstchar;		/* First character of column to use (1-n) */
static int *lastchar;		/* Last character of column to use (1-n) */
static int qmeancol = 0;	/* If 1, print mean of columns added in quadrature */
				/* 0 = no mean of quadruature */
main (ac, av)
int ac;
char **av;
{
    char *str;
    char *temp;
    char *filename;
    char *ranges = NULL;
    char *lfile = NULL;
    char *lranges = NULL;
    char *cdot, *ccol;
    int icol;
    int match, newbytes;
    int nrbytes = 0;

    cwhite = NULL;

    if (ac == 1)
        usage ();

    cond = (char **)calloc (maxncond, sizeof(char *));
    ccond = (char **)calloc (maxncond, sizeof(char *));
    op = (char **)calloc (maxnop, sizeof(char *));
    cop = (char **)calloc (maxnop, sizeof(char *));
    frstchar = (int *) calloc (MAX_NTOK, sizeof (int));
    lastchar = (int *) calloc (MAX_NTOK, sizeof (int));

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
	if (iscol (*av)) {
	    cdot = strchr (*av,'.');
	    if (cdot) {
		ccol = strchr (*av, ':');
		if (ccol) {
		    *cdot = (char) 0;
		    icol = atoi (*av);
		    *ccol = (char) 0;
		    if (icol > 0 && icol < MAX_NTOK) {
			if (strlen (cdot+1) > 0)
			    frstchar[icol] = atoi (cdot+1);
			else
			    frstchar[icol] = 1;
			if (strlen (ccol+1) > 0)
			    lastchar[icol] = atoi (ccol+1);
			else
			    lastchar[icol] = 0;
			}
		    }
		else
		    match = 1;
		}
		
	    if (ranges) {
		newbytes = strlen(ranges) + strlen(*av) + 2;
		newbytes = ((newbytes / 16) + 1) * 16;
		if (newbytes > nrbytes) {
		    temp = ranges;
		    ranges = (char *)calloc (newbytes, 1);
		    strcpy (ranges, temp);
		    nrbytes = newbytes;
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
		    nrbytes = newbytes;
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
		ccond = (char **)realloc((void *)ccond, maxncond*sizeof(void *));
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

	    case 'b':	/* bar-separated input */
		cwhite = (char *) calloc (4,1);
		strcpy (cwhite, "bar");
		break;

	    case 'c':	/* Count entries in each column */
		countcol++;
		break;

	    case 'd':	/* Number of decimal places in f.p. output */
		if (ac < 2)
		    usage ();
		ndec = atoi (*++av);
		ac--;
		break;

	    case 'e':	/* Compute median of each numeric column */
		medcol++;
		break;

	    case 'f':	/* Range of values in column */
		rangecol++;
		break;

	    case 'g':	/* Compute absolute median of each numeric column */
		amedcol++;
		break;

	    case 'h':	/* Print Starbase tab table header */
		printhead++;
		break;

	    case 'i':	/* tab-separated input */
		cwhite = (char *) calloc (4,1);
		strcpy (cwhite, "tab");
		break;

	    case 'j':	/* Compute absolute mean of each numeric column */
		ameancol++;
		break;

	    case 'k':	/* Count columns on first line */
		counttok++;
		break;

	    case 'l':	/* Number of lines to append */
		if (ac < 2)
		    usage ();
		napp = atoi (*++av);
		ac--;
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

	    case 'p': /* Print sum or mean only */
		printcol = 0;
		break;

	    case 'q':	/* Compute mean in quadrature of selected numeric columns */
		qmeancol = 1;
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

	    case 'x':	/* Value to ignore */
		if (ac < 2)
		    usage ();
		badval = atof (*++av);
		isbadval++;
		ac--;
		break;

	    default:
		usage ();
		break;
	    }
	    }

	/* Operation */
	else if (iscolop (*av)) {
	    if (nop >= maxnop) {
		maxnop = maxnop * 2;
		op = (char **)realloc((void *)op, maxncond*sizeof(void *));
		cop = (char **)realloc((void *)cop, maxncond*sizeof(void *));
		}
	    op[nop] = *av;
	    cop[nop] = strchr (*av, '*');
	    if (cop[nop] == NULL)
		cop[nop] = strchr (*av, '+');
	    if (cop[nop] == NULL)
		cop[nop] = strchr (*av, '/');
	    if (cop[nop] == NULL)
		cop[nop] = strchr (*av, 's');
	    if (cop[nop] == NULL)
		cop[nop] = strchr (*av, 'm');
	    if (cop[nop] == NULL)
		cop[nop] = strchr (*av, 'p');
	    nop++;
	    }

	/* File to read */
	else
	    filename = *av;
	}
    (void)ListFile (filename, ranges, lranges, lfile);

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
    fprintf (stderr,"Usage: [-abcefghijkmopqtv][-d num][-l num][-n num][-r lines][-s num] filename [col] [col] ...\n");
    fprintf(stderr," col: Number range (n1-n2,n3-n4...) or col.c1:c2\n");
    fprintf(stderr,"  -a: Sum selected numeric column(s)\n");
    fprintf(stderr,"  -b: Input is bar-separate table file\n");
    fprintf(stderr,"  -c: Add count of number of lines in each column at end\n");
    fprintf(stderr,"  -d num: Number of decimal places in f.p. output\n");
    fprintf(stderr,"  -e: Median values of selected numeric column(s)\n");
    fprintf(stderr,"  -f: Print range of values in selected column(s)\n");
    fprintf(stderr,"  -g: Median absolute values of selected column(s)\n");
    fprintf(stderr,"  -h: Print Starbase tab table header\n");
    fprintf(stderr,"  -i: Input is tab-separate table file\n");
    fprintf(stderr,"  -j: Print means of absolute values of selected column(s)\n");
    fprintf(stderr,"  -k: Print number of columns on first line\n");
    fprintf(stderr,"  -l num: Number of lines to add to each line\n");
    fprintf(stderr,"  -m: Print means of selected numeric column(s)\n");
    fprintf(stderr,"  -n num: Number of lines to read, if not all\n");
    fprintf(stderr,"  -o: OR conditions insted of ANDing them\n");
    fprintf(stderr,"  -p: Print only sum or mean, not individual values\n");
    fprintf(stderr,"  -q: Compute mean of selected columns added in quadrature\n");
    fprintf(stderr,"  -r: Range or @file of lines to read, if not all\n");
    fprintf(stderr,"  -s: Number of lines to skip\n");
    fprintf(stderr,"  -t: Starbase tab table output\n");
    fprintf(stderr,"  -v: Verbose\n");
    fprintf(stderr,"  -x num: Set value to ignore\n");
    exit (1);
}

static int
ListFile (filename, ranges, lranges, lfile)

char	*filename;	/* File name */
char	*ranges;	/* String with range of column numbers to list */
char	*lranges;	/* String with range of lines to list */
char	*lfile;		/* Name of file with lines to list */

{
    int i, j, il, ir, nbytes;
    char line[1024];
    char *nextline;
    char *lastchar;
    FILE *fd;
    FILE *lfd;
    int nlog;
    struct Tokens tokens;  /* Token structure */
    struct Range *range;
    struct Range *lrange;
    int *iline;
    int nline;
    int idnum;
    int iln;
    int nfdef = 9;
    double *sum, *asum, *colmin, *colmax, **med, **amed;
    double *qmed, qsum = 0.0;
    double qsum1;
    int *nsum;
    int *nent;
    int *hms;		/* Flag for hh:mm:ss or dd:mm:ss format */
    int *limset;	/* Flag for range initialization */
    int nlmax;
    double dtok, dnum;
    int nfind, ntok, nt, ltok,iop;
    int *inum;
    int icond, itok;
    char tcond, *cstr, *cval, top;
    char numstr[32], numstr1[32];
    double dcond, dval;
    int pass;
    int nchar, k;
    int iapp;
    int jcond, jval;
    int unset;
    char token[MAX_LTOK];
    int nq, nqsum = 0;

    nent = NULL;
    sum = NULL;
    inum = NULL;
    nsum = NULL;
    colmin = NULL;
    colmax = NULL;
    lrange = NULL;
    iline = NULL;
    hms = NULL;
    med = NULL;
    limset = NULL;

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
	    return (-1);
	    }
	for (i = 0; i < nline; i++)
	    iline[i] = rgeti4 (lrange);
	qsort (iline, nline, sizeof(int), intcompare);
	}

    /* Make list of line numbers to read from file specified on command line */
    if (lfile != NULL) {
	if (!(lfd = fopen (lfile, "r")))
            return (-1);
	nlmax = 99;
	nline = 0;
	nbytes = nlmax * sizeof(int);
	if (!(iline = (int *) calloc (nlmax, sizeof(int))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for iline\n", nbytes);
	    fclose (lfd);
	    return (-1);
	    }
	il = 0;
	nextline = line;
	for (ir = 0; ir < nread; ir++) {
	    if (fgets (nextline, 1023, lfd) == NULL)
		break;

	    /* Skip lines with comments */
	    if (line[0] == '#')
		continue;

	    /* Drop linefeeds */
	    lastchar = nextline + strlen(nextline) - 1;
	    if (*lastchar < 32)
		*lastchar = (char) 0;

	    /* Add lines with escaped linefeeds */
	    lastchar = nextline + strlen(nextline) - 1;
	    if (*lastchar == (char) 92) {
		nextline = lastchar;
		continue;
		}
	    else
		nextline = line;

	    ntok = setoken (&tokens, line, cwhite);
	    nt = 0;
	    il++;
	    if (il > nlmax) {
		nlmax = nlmax + 100;
		nbytes = nlmax * sizeof(int);
		if (!(iline = (int *) realloc ((void *) iline, nbytes))) {
		    fprintf (stderr, "Could not realloc %d bytes for iline\n",
			     nbytes);
		    fclose (lfd);
		    return (-1);
		    }
		}
	    for (i = 0; i < ntok; i++) {
		if (getoken (tokens, i+1, token, MAX_LTOK)) {
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
        return (-1);
	}

    /* Skip lines into input file */
    if (nskip > 0) {
	for (i = 0; i < nskip; i++) {
	    if (fgets (line, 1023, fd) == NULL)
		break;
	    }
	}

    /* Print entire selected lines */
    if (ranges == NULL) {
	iln = 0;
	il = 0;
	iapp = 0;
	nextline = line;
	for (ir = 0; ir < nread; ir++) {
	    if (fgets (nextline, 1023, fd) == NULL)
		break;

	    /* Skip lines with comments */
	    if (nextline[0] == '#')
		continue;

	    /* Add lines with escaped linefeeds */
	    lastchar = nextline + strlen(nextline) - 1;
	    if (*lastchar == (char) 92) {
		nextline = lastchar;
		continue;
		}
	    else if (iapp++ < napp) {
		*lastchar = ' ';
		nextline = lastchar + 1;
		continue;
		}
	    else {
		iapp = 0;
		nextline = line;
		}

	    il++;

	    /* Skip if line is not on list, if there is one */
	    if (iline != NULL) {
		if (il+1 < iline[iln])
		    continue;
		else if (il+1 > iline[nline-1])
		    break;
		else
		    iln++;
		}

	    /* Drop control character at end of string */
	    lastchar = line + strlen(line) - 1;
	    if (*lastchar < 32)
		*lastchar = (char) 0;

	    /* Drop second control character at end of string */
	    lastchar = line + strlen(line) - 1;
	    if (*lastchar < 32)
		*lastchar = (char) 0;

	    /* Echo line if it is a comment */
	    if (line[0] == '#') {
		printf ("%s\n", line);
		continue;
		}

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
		    getoken (tokens, itok, token, MAX_LTOK);
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
		    if (condand && !pass)
			break;
		    }
		if (pass)
		    printf ("%s\n", line);
		}
	    else
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
	    return (-1);
	    }
	if (!(asum = (double *) calloc (nfind, sizeof (double))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for asum\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return (-1);
	    }
	if (!(colmin = (double *) calloc (nfind, sizeof (double))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for colmin\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return (-1);
	    }
	if (!(colmax = (double *) calloc (nfind, sizeof (double))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for colmax\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return (-1);
	    }
	if (!(med = (double **) calloc (nfind, sizeof (double *))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for med\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return (-1);
	    }
	if (!(amed = (double **) calloc (nfind, sizeof (double *))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for amed\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return (-1);
	    }
	else {
	    nbytes = nread * sizeof (double);
	    if (!(qmed = calloc (nread, sizeof(double)))) {
		fprintf (stderr, "Could not calloc %d bytes for qmed\n", nbytes);
		if (fd != stdin) fclose (fd);
		return (-1);
		}
	    for (i = 0; i < nfind; i++) {
		if (!(med[i] = calloc (nread, sizeof(double)))) {
		    fprintf (stderr, "Could not calloc %d bytes for med%d\n",
			     nbytes, i);
		    if (fd != stdin) fclose (fd);
		    return (-1);
		    }
		if (!(amed[i] = calloc (nread, sizeof(double)))) {
		    fprintf (stderr, "Could not calloc %d bytes for amed%d\n",
			     nbytes, i);
		    if (fd != stdin) fclose (fd);
		    return (-1);
		    }
		}
	    }
	nbytes = nfind * sizeof (int);
	if (!(nsum = (int *) calloc (nfind, sizeof(int))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for nsum\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return (-1);
	    }
	if (!(hms = (int *) calloc (nfind, sizeof(int))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for hms\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return (-1);
	    }
	if (!(limset = (int *) calloc (nfind, sizeof(int))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for limset\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return (-1);
	    }
	if (!(nent = (int *) calloc (nfind, sizeof(int))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for nent\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return (-1);
	    }
	if (!(inum = (int *) calloc (nfind, sizeof(int))) ) {
	    fprintf (stderr, "Could not calloc %d bytes for inum\n", nbytes);
	    if (fd != stdin) fclose (fd);
	    return (-1);
	    }
	else {
	    for (i = 0; i < nfind; i++)
		inum[i] = rgeti4 (range);
	    }
	iln = 0;
	iapp = 0;
	nextline = line;
	for (il = 0; il < nread; il++) {
	    if (fgets (nextline, 1023, fd) == NULL)
		break;

	    /* Clear control character at end of string */
	    lastchar = line + strlen(line) - 1;
	    if (*lastchar < 32)
		*lastchar = (char) 0;

	    /* Clear second control character at end of string */
	    lastchar = line + strlen(line) - 1;
	    if (*lastchar < 32)
		*lastchar = (char) 0;

	    /* Add lines with escaped linefeeds */
	    lastchar = nextline + strlen(nextline) - 1;
	    if (*lastchar == (char) 92) {
		nextline = lastchar;
		continue;
		}
	    else if (iapp++ < napp) {
		*++lastchar = ' ';
		nextline = lastchar + 1;
		continue;
		}
	    else
		nextline = line;
		iapp = 0;

	    /* Skip if line is not on list, if there is one */
	    if (iline != NULL) {
		if (il+1 < iline[iln])
		    continue;
		else if (il > iline[nline-1])
		    break;
		else
		    iln++;
		}

	    /* Echo line if it is a comment */
	    if (line[0] == '#') {
		printf ("%s\n", line);
		continue;
		}

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
		    getoken (tokens, itok, token, MAX_LTOK);
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
		    if (getoken (tokens, inum[i], token, MAX_LTOK)) {
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
		    if (getoken (tokens, inum[i], token, MAX_LTOK)) {
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

	    /* Print requested columns */
	    qsum1 = 0;
	    nq = 0;
	    for (i = 0; i < nfind; i++) {
		if (getoken (tokens, inum[i], token, MAX_LTOK)) {

		    /* Get substring of column, if requested */
		    if (frstchar[i] > 0) {
			ltok = strlen (token);
			if (lastchar[i] > 0)
			    nchar = lastchar[i] - frstchar[i] + 1;
			else
			    nchar = ltok - frstchar[i] + 1;
			k = frstchar[i];
			for (j = 0; j < nchar; j++)
			    token[j] = token[k++];
			for (j = nchar; j < ltok; j++);
			    token[j] = (char) 0;
			}

		    if (isnum (token)) {
			dval = atof (token);
			hms[i] = 0;
			idnum = 1;
			}
		    else if (strchr (token, ':')) {
			dval = str2dec (token);
			hms[i] = 1;
			idnum = 1;
			}
		    else {
			dval = 0.0;
			idnum = 0;
			}
		    if (idnum) {
			if (!isbadval || (isbadval && dval != badval)) {
			    qsum1 = qsum1 + dval * dval;
			    nq++;
			    sum[i] = sum[i] + dval;
			    asum[i] = asum[i] + fabs (dval);
			    if (!limset[i]) {
				colmin[i] = dval;
				colmax[i] = dval;
				limset[i]++;
				}
			    else if (dval < colmin[i])
				colmin[i] = dval;
			    else if (dval > colmax[i])
				colmax[i] = dval;
			    med[i][nsum[i]] = dval;
			    amed[i][il] = fabs (dval);
			    nsum[i]++;
			    }
			}
		    if (printcol) {
			if (i > 0) {
			    if (tabout)
				printf ("	");
			    else
				printf (" ");
			    }
			if (inum[i] > tokens.ntok || inum[i] < 1)
			    printf ("___");
			else if (ndec > -1 && isnum (token) == 2) {
			    num2str (numstr, atof (token), 0, ndec);
			    printf ("%s", numstr);
			    }
			else
			    printf ("%s", token);
			}
		    nt++;
		    nent[i]++;
		    }
		if (nq > 1) {
		    qmed[nqsum] = sqrt (qsum1);
		    nqsum++;
		    qsum = qsum + sqrt (qsum1);
		    }
		}

	    /* Print columns being operated on */
	    for (iop = 0; iop < nop; iop++) {
		if (i > 0 || nfind > 0) {
		    if (tabout)
			printf ("	");
		    else
			printf (" ");
		    }

		/* Extract test value from comparison string */
		top = *cop[iop];
		*cop[iop] = (char) 0;
		cstr = cop[iop]+1;
		if (isnum (cstr) > 1)
		    dnum = atof (cstr);
		else {
		    itok = atoi (cstr);
		    if (getoken (tokens, itok, token, MAX_LTOK))
			dnum = atof (token);
		    else {
			printf ("___");
			continue;
			}
		    }
		*cop[iop] = top;

		/* Extract token from input line */
		itok = atoi (op[iop]);
		if (getoken (tokens, itok, token, MAX_LTOK)) {
		    dtok = atof (token);
		    if (top == '+' || top == 'a')
			printf ("%f", dtok + dnum);
		    else if (top == '-' || top == 's')
			printf ("%f", dtok - dnum);
		    else if (top == '*' || top == 'm')
			printf ("%f", dtok * dnum);
		    else if (top == '/' || top == 'd')
			printf ("%f", dtok / dnum);
		    else
			printf ("___");
		    }
		}
	    if (nt > 0 && printcol)
		printf ("\n");
	    }
        }

    /* Print sums of values in numeric columns */
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

    /* Print means of absolute values in numeric columns */
    if (ameancol) {
	for (i = 0; i < nfind; i++) {
	    if (nsum[i] > 0) {
		dval = asum[i] / (double)nsum[i];
		if (hms[i])
		    dec2str (numstr, 32, dval, 3);
		else
		    sprintf (numstr, "%f", dval);
		strclean (numstr);
		if (i < nfind-1)
		    printf ("%s ", numstr);
		else
		    printf ("%s", numstr);
		}
	    else if (i < nfind-1)
		printf ("___ ");
	    else
		printf ("___");
	    }
	if (nfind > 0)
	    printf ("\n");
	}

    /* Print means of values in numeric columns */
    if (meancol) {
	for (i = 0; i < nfind; i++) {
	    if (nsum[i] > 0) {
		dval = sum[i] / (double)nsum[i];
		if (hms[i])
		    dec2str (numstr, 32, dval, 3);
		else
		    sprintf (numstr, "%f", dval);
		strclean (numstr);
		if (i < nfind-1)
		    printf ("%s ", numstr);
		else
		    printf ("%s", numstr);
		}
	    else if (i < nfind-1)
		printf ("___ ");
	    else
		printf ("___");
	    }
	if (nfind > 0)
	    printf ("\n");
	}

    /* Print medians of absolute values in numeric columns */
    if (amedcol) {
	for (i = 0; i < nfind; i++) {
	    if (nsum[i] > 0) {
		dval = median (amed[i], nsum[i]);
		if (hms[i])
		    dec2str (numstr, 32, dval, 3);
		else
		    sprintf (numstr, "%f", dval);
		strclean (numstr);
		if (i < nfind-1)
		    printf ("%s ", numstr);
		else
		    printf ("%s", numstr);
		}
	    else if (i < nfind-1)
		printf ("___ ");
	    else
		printf ("___");
	    }
	if (nfind > 0)
	    printf ("\n");
	}

    /* Print medians of values in numeric columns */
    if (medcol) {
	for (i = 0; i < nfind; i++) {
	    if (nsum[i] > 0) {
		dval = median (med[i], nsum[i]);
		if (hms[i])
		    dec2str (numstr, 32, dval, 3);
		else
		    sprintf (numstr, "%f", dval);
		strclean (numstr);
		if (i < nfind-1)
		    printf ("%s ", numstr);
		else
		    printf ("%s", numstr);
		}
	    else if (i < nfind-1)
		printf ("___ ");
	    else
		printf ("___");
	    }
	if (nfind > 0)
	    printf ("\n");
	}

    /* Print ranges of values in numeric columns */
    if (rangecol) {
	for (i = 0; i < nfind; i++) {
	    if (hms[i])
		dec2str (numstr, 32, colmin[i], 3);
	    else
		sprintf (numstr, "%f", colmin[i]);
	    strclean (numstr);
	    printf ("%s-", numstr);
	    if (hms[i])
		dec2str (numstr, 32, colmax[i], 6);
	    else
		sprintf (numstr, "%f", colmax[i]);
	    strclean (numstr);
	    if (i < nfind-1)
		printf ("%s ", numstr);
	    else
		printf ("%s", numstr);
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

    /* Print mean of all numeric columns added in quadrature */
    if (qmeancol) {
	if (nqsum > 0) {
	    dval = qsum / (double)nqsum;
	    sprintf (numstr, "%f", dval);
	    strclean (numstr);
	    printf ("%s", numstr);
	    dval = median (qmed, nqsum);
	    sprintf (numstr, "%f", dval);
	    strclean (numstr);
	    printf (" %s\n", numstr);
	    }
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
intcompare (i, j)

int *i, *j;
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

static int
iscolop (string)

char *string;
{
    /* Check for presence of operation */
    if (strchr (string, '+') != NULL || strchr (string, '-') != NULL ||
	strchr (string, '*') != NULL || strchr (string, '/') != NULL ||
	strchr (string, 'a') != NULL || strchr (string, 's') != NULL ||
	strchr (string, 'm') != NULL || strchr (string, 'd') != NULL) {

	/* Check to see if string is a file */
	if (access(string,0) && strcmp (string, "stdin"))
	    return (1);
	else
	    return (0);
	}
    else
	return (0);
}

static double
median (x, n)

int	n;
double	*x;
{
    int rhs, lhs;
    int NComp();

    if (n <= 0)
	return (0.0);
   else if (n == 1)
	return (x[0]);

    qsort (x, n, sizeof (double), NComp);

    lhs = (n - 1) / 2;
    rhs = n / 2;

    if (lhs == rhs)
	return (x[lhs]);
    else
	return ((x[lhs] + x[rhs]) / 2.0);
}

iscol (string)

char *string;   /* Character string */
{
    int lstr, i, nd;
    char cstr, cstr1;
    int fpcode;

    /* Return 0 if string is NULL */
    if (string == NULL)
        return (0);

    lstr = strlen (string);
    nd = 0;
    fpcode = 1;

    /* Remove trailing spaces */
    while (string[lstr-1] == ' ')
        lstr--;

    /* Column strings contain 0123456789 and . or : for subranges */
    for (i = 0; i < lstr; i++) {
        cstr = string[i];
        if (cstr == '\n')
            break;

        /* Ignore leading spaces */
        if (cstr == ' ' && nd == 0)
            continue;

        if ((cstr < 48 || cstr > 57) &&
            cstr != ':' && cstr != '.' )
            return (0);
	else if (cstr >= 47 && cstr <= 57)
	    nd++;
        }
    if (nd > 0)
        return (1);
    else
        return (0);
}


int
NComp (pd1, pd2)

void *pd1, *pd2;
{
    double d1 = *((double *)pd1);
    double d2 = *((double *)pd2);

    if (d1 > d2)
	return (1);
    else if (d1 < d2)
	return (-1);
    else
	return (0);
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
 * Apr  4 2000	Add option to operate on keyword values
 * Apr  5 2000	Simply echo lines starting with #
 * Jul 21 2000	Link lines with escaped linefeeds at end to next line
 * Aug  8 2000	Add -l option to append lines without backslashes
 * Oct 23 2000	Declare ListFile(); fix column arithmetic command line options
 * Nov 20 2000	Clean up code using lint
 * Dec 11 2000	Include fitshead.h for string search
 *
 * Jan 17 2001	Add -d option to set number of output decimal places
 * Jan 17 2001	Add a, s, m, d for add, subtract, multiply, divide const or col
 * Mar 19 2001	Drop type declarations from intcompare argument list
 * Jun 18 2001	Add maximum length of returned string to getoken()
 * Jul 17 2001	Check operations for stdin as well as file
 * Oct 10 2001	Add sum, mean, sigma for hh:mm:ss and dd:mm:ss entries
 * Oct 10 2001	Add -p option to print only sum, mean, sigma, not entries
 * Oct 11 2001	Add -f option to print range of values in selected columns
 * Oct 16 2001	Add -e option to compute medians
 * Oct 16 2001	Ignore non-numeric values for sums, means, and medians
 * Nov 13 2001	Add option to specifiy characters of column to use
 * Dec 28 2001	Clear second control character at the end of a line (CR/LF)
 *
 * Apr  9 2002	NComp() cannot be static as it is passed to qsort()
 * Apr 10 2002	Add option to print means and medians of absolute values
 * Apr 11 2002	Add option to print mean of selected columns added in quadrature
 * Apr 12 2002	Add option to print median of selected columns added in quadrature
 * Apr 12 2002	Fix bug in computing median of filtered file
 * Apr 12 2002	Add -x option to set ignorable value
 */
