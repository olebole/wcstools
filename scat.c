/* File scat.c
 * July 5, 2001
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
#include "libwcs/wcs.h"
#include "libwcs/lwcs.h"
#include "libwcs/wcscat.h"
#include "libwcs/fitsfile.h"

static void usage();
static int scatparm();
static void scatcgi();

static int ListCat ();
extern void setcenter();
static void SearchHead();
static int GetArea();

static int verbose = 0;		/* Verbose/debugging flag */
static int afile = 0;		/* True to append output file */
static int wfile = 0;		/* True to print output file */
static int classd = -1;		/* Guide Star Catalog object classes */
static double maglim1 = MAGLIM1; /* Catalog bright magnitude limit */
static double maglim2 = MAGLIM2; /* Catalog faint magnitude limit */
static int sysout0 = 0;		/* Output coordinate system */
static double eqcoor = 2000.0;	/* Equinox of search center */
static int degout0 = 0;		/* 1 if degrees output instead of hms */
static double ra0 = -99.0;	/* Initial center RA in degrees */
static double dec0 = -99.0;	/* Initial center Dec in degrees */
static double rad0 = 0.0;	/* Search box radius */
static double dra0 = 0.0;	/* Search box width */
static double ddec0 = 0.0;	/* Search box height */
static double epoch0 = 0.0;	/* Epoch for coordinates */
static int syscoor = 0;		/* Input search coordinate system */
static int nstars = 0;		/* Number of brightest stars to list */
static int printhead = 0;	/* 1 to print table heading */
static int nohead = 1;		/* 1 to print table heading */
static searchcenter = 0;	/* 1 to print simpler format */
static int tabout = 0;		/* 1 for tab table to standard output */
static int catsort = SORT_MAG;	/* Default to sort stars by magnitude */
static int debug = 0;		/* True for extra information */
static char *objname;		/* Object name for output */
static char *keyword;		/* Column to add to tab table output */
static char *progname;		/* Name of program as executed */
static int version = 0;		/* If 1, print only program name and version */
static int match = 0;		/* If 1, match num exactly in BIN or ASC cats*/
static int printobj = 0;	/* If 1, print object name instead of number */
static char cpname[16];		/* Name of program for error messages */
static int oneline = 0;		/* If 1, print center and closest on 1 line */
static struct Star *srch;	/* Search center structure for catalog search */
static struct StarCat *srchcat; /* Search catalog structure */
static int readlist = 0;	/* If 1, search centers are from a list */
static int notprinted = 1;	/* If 1, print header */
static char *listfile;		/* Name of catalog file with search centers */
static int printxy = 0;		/* If 1, print X Y instead of object number */
static char *xstr, *ystr;	/* X and Y strings if printxy */
static int closest;		/* 1 if printing only closest star */
static double *gnum;		/* Catalog star numbers */
static double *gra;		/* Catalog star right ascensions */
static double *gdec;		/* Catalog star declinations */
static double *gpra;		/* Catalog star RA proper motions */
static double *gpdec;		/* Catalog star declination proper motions */
static double *gm;		/* Catalog magnitudes */
static double *gmb;		/* Catalog B magnitudes */
static double *gx;		/* Catalog star X positions on image */
static double *gy;		/* Catalog star Y positions on image */
static int *gc;			/* Catalog star object classes */
static char **gobj;		/* Catalog star object names */
static char **gobj1;		/* Catalog star object names */
static int nalloc = 0;
static struct StarCat *starcat[5]; /* Star catalog data structure */
static double eqout = 0.0;	/* Equinox for output coordinates */
static int ncat = 0;		/* Number of reference catalogs to search */
static char *refcatname[5];	/* reference catalog names */
static char *ranges;		/* Catalog numbers to print */
static int http=0;		/* Set to one if http header needed on output */
static int padspt = 0;		/* Set to one to pad out long spectral type */

main (ac, av)
int ac;
char **av;
{
    FILE *fd;
    char *str;
    char rastr[16];
    char decstr[16];
    char line[200];
    int i, lcat;
    char *refcatn;
    char cs;
    int srchtype;
    char *temp;
    int systemp = 0;		/* Input search coordinate system */
    int istar;
    char *blank;
    double epoch;
    int nmag, mprop;
    char title[80];
    char *query;
    int ndcat;
    int sysref;		/* Coordinate system of reference catalog */
    double eqref;	/* Equinox of catalog to be searched */
    double epref;	/* Epoch of catalog to be searched */
    int refcat;
    int lrange;
    char *rstr, *dstr, *cstr;

    ranges = NULL;
    keyword = NULL;
    objname = NULL;
    for (i = 0; i < 5; i++)
	starcat[i] = NULL;

    /* Null out buffers before starting */
    gnum = NULL;
    gra = NULL;
    gdec = NULL;
    gpra = NULL;
    gpdec = NULL;
    gm = NULL;
    gmb = NULL;
    gx = NULL;
    gy = NULL;
    gc = NULL;
    gobj = NULL;
    gobj1 = NULL;

    /* Check name used to execute programe and set catalog name accordingly */
    progname = ProgName (av[0]);
    for (i = 0; i < strlen (progname); i++) {
	if (progname[i] > 95 && progname[i] < 123)
	    cpname[i] = progname[i] - 32;
	else
	    cpname[i] = progname[i];
	}
    refcatn = ProgCat (progname);
    if (refcatn != NULL) {
	refcatname[ncat] = refcatn;
	ncat++;
	refcat = RefCat (refcatn,title,&sysref,&eqref,&epref,&mprop,&nmag);
	ndcat = CatNdec (refcat);
	}
    else
	ndcat = -1;

    /* Set parameters from keyword=value arguments */
    if ((query = getenv ("QUERY_STRING")) != NULL) {
        scatcgi (query);
	http++;
	}

    /* If not http and no arguments, print command list */
    else if (ac == 1)
        usage (progname, NULL);

    /* Check for help or version command first */
    if (!http) {
	str = *(av+1);
	if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	    usage (progname, NULL);
	if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	    version = 1;
	    usage (progname, NULL);
	    }
	if (!strncmp (str, "band", 4) || !strncmp (str, "filt", 4)) {
	    fprintf (stderr, "HST Guide Star Catalog Bandpass Codes\n");
	    fprintf (stderr, "Code Bandpass Emulsion/Filter Notes\n");
	    fprintf (stderr, " 0      J     IIIaJ+GG395     SERC-J/EJ\n");
	    fprintf (stderr, " 1      V     IIaD+W12        Pal Quick-V\n");
	    fprintf (stderr, " 3      B     -               Johnson\n");
	    fprintf (stderr, " 4      V     -               Johnson\n");
	    fprintf (stderr, " 5      R     IIIaF+RG630     -\n");
	    fprintf (stderr, " 6      V495  IIaD+GG495      Pal QV/AAO XV\n");
	    fprintf (stderr, " 7      O     103aO+no filt   POSS-I Blue\n");
	    fprintf (stderr, " 8      E     103aE+redplex   POSS-I Red\n");
	    fprintf (stderr, " 9      R     IIIaF+RG630     -\n");
	    fprintf (stderr, "10      -     IIaD+GG495+yel  GPO Astrograph\n");
	    fprintf (stderr, "11      -     103aO+blue      Black Birch Astrograph\n");
	    fprintf (stderr, "12      -     103aO+blue      Black Birch Astrograph (GSC cal)\n");
	    fprintf (stderr, "13      -     103aG+GG495+yel Black Birch Astrograph\n");
	    fprintf (stderr, "14      -     103aG+GG495+yel Black Birch Astrograph (GSC cal)\n");
	    fprintf (stderr, "16      J     IIIaJ+GG495     -\n");
	    fprintf (stderr, "18      V     IIIaJ+GG385     POSS-II Blue\n");
	    fprintf (stderr, "19      U     -               Johnson\n");
	    fprintf (stderr, "20      R     -               Johnson\n");
	    fprintf (stderr, "21      I     -               Johnson\n");
	    fprintf (stderr, "22      U     -               Cape\n");
	    fprintf (stderr, "23      R     -               Kron\n");
	    fprintf (stderr, "24      I     -               Kron\n");
	    exit (1);
	    }
	if (!strncmp (str, "clas", 4) || !strncmp (str, "obje", 4)) {
	    fprintf (stderr, "HST Guide Star Catalog Object Classes\n");
	    fprintf (stderr, "0: Stellar\n");
	    fprintf (stderr, "3: Non-Stellar\n");
	    fprintf (stderr, "5: Not really an object\n");
	    exit (1);
	    }
	}

    /* crack arguments */
    for (av++; --ac > 0; av++) {

	/* Set parameters from keyword=value& arguments */
	if (strchr (*av, '&')) {
            scatcgi (*av);
	    }

	/* Set parameters from keyword=value arguments */
	else if (strchr (*av, '=')) {
            if (scatparm (*av))
		fprintf (stderr, "SCAT: %s is not a parameter.\n", *av);
	    }

	/* Set search RA, Dec, and equinox if colon in argument */
	else if (strsrch (*av,":") != NULL) {
	    if (ac < 2)
		usage(progname, *av);
	    else {
		strcpy (rastr, *av);
		ac--;
		strcpy (decstr, *++av);
		ra0 = str2ra (rastr);
		dec0 = str2dec (decstr);
		ac--;
		if (ac < 1) {
		    syscoor = WCS_J2000;
		    eqcoor = 2000.0;
		    }
		else if ((syscoor = wcscsys (*(av+1))) >= 0)
		    eqcoor = wcsceq (*++av);
		else {
		    syscoor = WCS_J2000;
		    eqcoor = 2000.0;
		    }
		}
	    }

	/* Set range and make a list of star numbers from it */
	else if (isrange (*av)) {
	    if (ranges) {
		temp = ranges;
		lrange = strlen(ranges) + strlen(*av) + 16;
		ranges = (char *) calloc (lrange, 1);
		strcpy (ranges, temp);
		strcat (ranges, ",");
		strcat (ranges, *av);
		free (temp);
		}
	    else {
		lrange = strlen(*av) + 16;
		ranges = (char *) calloc (lrange, 1);
		if (strchr (*av,'.'))
		    match = 1;
		strcpy (ranges, *av);
		}
	    if (strchr (*av, 'x') == NULL && ndcat > 0) {
		int n = ndcat;
		strcat (ranges, "x0.");
		while (n-- > 0)
		    strcat (ranges, "0");
		strcat (ranges, "1");
		}
	    }

	/* Set decimal degree center or star number */
	else if (isnum (*av)) {

	    if (ac > 1 && isnum (*(av+1))) {
		int ndec1, ndec2;

		/* Check for second number and coordinate system */
		if (ac > 2 && (systemp = wcscsys (*(av + 2))) > 0) {
		    rstr = *av++;
		    ac--;
		    dstr = *av++;
		    ac--;
		    cstr = *av;
		    }

		/* Check for two numbers which aren't catalog numbers */
		else {
		    ndec1 = StrNdec (*av);
		    ndec2 = StrNdec (*(av+1));
		    if (ndcat > -1 && ndec1 == ndcat && ndec2 == ndcat)
			cstr = NULL;
		    else {
			rstr = *av++;
			ac--;
			dstr = *av;
			cstr = (char *) malloc (8);
			strcpy (cstr, "J2000");
			systemp = WCS_J2000;
			}
		    }
		}
	    else
		cstr = NULL;

	    /* Set decimal degree center */
	    if (cstr != NULL) {
		ra0 = atof (rstr);
		dec0 = atof (dstr);
		syscoor = systemp;
		eqcoor = wcsceq (cstr);
		}

	    /* Assume number to be star number if no coordinate system */
	    else {
		if (strchr (*av,'.'))
		    match = 1;
		if (ranges) {
		    temp = ranges;
		    lrange = strlen(ranges) + strlen(*av) + 2;
		    ranges = (char *)calloc (lrange, 1);
		    strcpy (ranges, temp);
		    strcat (ranges, ",");
		    strcat (ranges, *av);
		    free (temp);
		    }
		else {
		    lrange = strlen(*av) + 2;
		    ranges = (char *) calloc (lrange, 1);
		    strcpy (ranges, *av);
		    }
		}
	    }

	else if (*(str = *av) == '@') {
	    readlist = 1;
	    listfile = *av + 1;
	    }

	/* Otherwise, read command */
	else if (*(str = *av) == '-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

	    case 'v':	/* more verbosity */
		if (verbose)
		    debug++;
		else
		    verbose++;
		break;

	    case 'a':	/* Get closest source */
		catsort = SORT_DIST;
		closest++;
		break;

    	    case 'b':	/* output coordinates in B1950 */
		sysout0 = WCS_B1950;
		eqout = 1950.0;
    		break;

	    case 'c':       /* Set reference catalog */
		if (ac < 2)
		    usage(progname, *av);
		lcat = strlen (*++av);
		refcatn = (char *) calloc (1, lcat + 2);
		strcpy (refcatn, *av);
		refcatname[ncat] = refcatn;
		ncat = ncat + 1;
		ac--;
		break;

	    case 'd':	/* output in degrees instead of sexagesimal */
		degout0++;
		break;

	    case 'e':	/* Set ecliptic coordinate output */
		sysout0 = WCS_ECLIPTIC;
		break;

	    case 'f':	/* output in centering coordinates */
		searchcenter++;
		break;

	    case 'g':	/* Set galactic coordinate output and optional center */
		sysout0 = WCS_GALACTIC;
		break;

	    case 'h':	/* ouput descriptive header */
		printhead++;
		break;

	    case 'i':	/* ouput catalog object name instead of number */
		printobj++;
		break;

    	    case 'j':	/* center coordinates on command line in J2000 */
		sysout0 = WCS_J2000;
		eqout = 2000.0;
    		break;

	    case 'k':	/* Keyword (column) to add to output from tab table */
		if (ac < 2)
		    usage (progname, *av);
		keyword = *++av;
		settabkey (keyword);
		ac--;
		break;

	    case 'l':	/* Print center and closest star on one line */
		oneline++;
		catsort = SORT_DIST;
		closest++;
		nstars = 1;
		break;

	    case 'm':	/* Magnitude limit */
		if (ac < 2)
		    usage (progname, *av);
		maglim2 = atof (*++av);
		ac--;
		if (ac > 1 && isnum (*(av+1))) {
		    maglim1 = maglim2;
		    maglim2 = atof (*++av);
		    ac--;
		    }
		else if (MAGLIM1 == MAGLIM2)
		    maglim1 = -2.0;
		break;

	    case 'n':	/* Number of brightest stars to read */
		if (ac < 2)
		    usage (progname, *av);
		nstars = atoi (*++av);
		ac--;
		break;

	    case 'o':	/* Object name */
		if (ac < 2)
		    usage (progname, *av);
		objname = *++av;
		ac--;
		break;

	    case 'p':	/* Sort by distance from center */
		catsort = SORT_DIST;
		break;

    	    case 'q':	/* Output equinox in years */
    		if (ac < 2)
    		    usage(progname, *av);
    		eqout = fd2ep (*++av);
    		ac--;
    		break;

    	    case 'r':	/* Box radius in arcseconds */
    		if (ac < 2)
    		    usage(progname, *av);
		av++;
		if ((dstr = strchr (*av, ',')) != NULL) {
		    *dstr = (char) 0;
		    dstr++;
		    }
		if (strchr (*av,':'))
		    rad0 = 3600.0 * str2dec (*av);
		else
		    rad0 = atof (*av);
		if (dstr != NULL) {
		    dra0 = rad0;
		    rad0 = 0.0;
		    if (strchr (dstr, ':'))
			ddec0 = 3600.0 * str2dec (dstr);
		    else
			ddec0 = atof (dstr);
		    if (ddec0 <= 0.0)
			ddec0 = dra0;
		    /* rad0 = sqrt (dra0*dra0 + ddec0*ddec0); */
		    }
    		ac--;
    		break;

	    case 's':	/* sort by RA  or anything else */
		catsort = SORT_RA;
		if (ac > 1) {
		    cs = *(av+1)[0];
		    if (strchr ("dmnpr",(int)cs)) {
			av++;
			ac--;
			}
		    else
			cs = 'r';
		    }
		else
		    cs = 'r';
		if (cs) {

		    /* Declination */
		    if (cs == 'd')
			catsort = SORT_DEC;

		    /* Magnitude (brightest first) */
		    else if (cs == 'm')
			catsort = SORT_MAG;

		    /* No sorting */
		    else if (cs == 'n')
			catsort = NOSORT;

		    /* Distance from search center (closest first) */
		    else if (cs == 'p')
			catsort = SORT_DIST;

		    /* Right ascension */
		    else if (cs == 'r')
			catsort = SORT_RA;
		    else
			catsort = SORT_RA;
		    }
		else
		    catsort = SORT_RA;
		break;

	    case 't':	/* tab table to stdout */
		tabout = 1;
		break;

	    case 'u':       /* Print following 2 numbers at start of line */
		if (ac > 2) {
		    printxy = 1;
		    xstr = *++av;
		    ac--;
		    ystr = *++av;
		    ac--;
		    }
		break;

    	    case 'w':	/* write output file */
    		wfile++;
    		break;

	    case 'x':       /* Guide Star object class */
		if (ac < 2)
		    usage(progname, *av);
		classd = (int) atof (*++av);
		setgsclass (classd);
		ac--;
		break;

	    case 'y':	/* Set output coordinate epoch */
		if (ac < 2)
		    usage(progname, *av);
		epoch0 = fd2ep (*++av);
		ac--;
		break;

	    case 'z':	/* Set append flag */
		afile++;
		wfile++;
		break;

	    default:
		usage (progname, *av);
		break;
	    }
	    }
	else {
	    lcat = strlen (*av);
	    refcatn = (char *) calloc (1, lcat + 2);
	    strcpy (refcatn, *av);
	    refcatname[ncat] = refcatn;
	    ncat = ncat + 1;
	    }
	}

    /* Set output epoch appropriately if output system is specified */
    if (sysout0 == 0 && syscoor != 0)
	sysout0 = syscoor;
    if (epoch0 == 0.0 && sysout0 != 0) {
	if (sysout0 == WCS_J2000)
	    epoch0 = 2000.0;
	if (sysout0 == WCS_B1950)
	    epoch0 = 1950.0;
	}

    /* Set output equinox appropriately if output system is specified */
    if (eqout == 0.0 && sysout0 != 0) {
	if (sysout0 == WCS_J2000)
	    eqout = 2000.0;
	if (sysout0 == WCS_B1950)
	    eqout = 1950.0;
	}

    /* If http output, send header */
    if (http)
	printf ("Content-type: text/plain\n\n");

    if (readlist) {

	/* Read search center list from starbase tab table catalog */
	if (istab (listfile)) {
	    ranges = NULL;
	    srchcat = tabcatopen (listfile, NULL, 10000);
	    if (srchcat != NULL) {
		srch = (struct Star *) calloc (1, sizeof (struct Star));
		for (istar = 1; istar <= srchcat->nstars; istar ++) {
		    if (tabstar (istar, srchcat, srch, verbose)) {
			if (verbose)
			    fprintf (stderr,"%s: Cannot read star %d\n",
				     cpname, istar);
                	break;
                	}
		    ra0 = srch->ra;
		    dec0 = srch->dec;
		    if (eqout > 0.0)
			eqcoor = eqout;
		    else
			eqcoor = srch->equinox;
		    if (epoch0 != 0.0)
			epoch = epoch0;
		    else
			epoch = srch->epoch;
		    if (sysout0)
			syscoor = sysout0;
		    else
			syscoor = srch->coorsys;
		    wcsconp (srch->coorsys, syscoor, srch->equinox, eqcoor,
			     srch->epoch,epoch,
			     &ra0,&dec0,&srch->rapm,&srch->decpm);
		    ListCat (ranges, eqout);
		    }
		tabcatclose (srchcat);
		}
	    }

	/* Read search center list from SAOTDC ASCII table catalog */
	else if (isacat (listfile)) {
	    ranges = NULL;
	    if (!(srchtype = RefCat (listfile,title,&syscoor,&eqcoor,
				     &epoch,&mprop,&nmag))) {
		fprintf (stderr,"List catalog '%s' is missing\n", listfile);
		return (0);
		}
	    srchcat = ctgopen (listfile, srchtype);
	    if (srchcat != NULL) {
		srch = (struct Star *) calloc (1, sizeof (struct Star));
		for (istar = 1; istar <= srchcat->nstars; istar ++) {
		    if (ctgstar (istar, srchcat, srch)) {
			if (verbose)
			    fprintf (stderr,"%s: Cannot read star %d\n",
				 cpname, istar);
                	}
		    else {
			ra0 = srch->ra;
			dec0 = srch->dec;
			if (eqout > 0.0)
			    eqcoor = eqout;
			else
			    eqcoor = srch->equinox;
			if (epoch0 != 0.0)
			    epoch = epoch0;
			else
			    epoch = srch->epoch;
			if (sysout0)
			    syscoor = sysout0;
			else
			    syscoor = srch->coorsys;
			wcsconp (srch->coorsys, syscoor, srch->equinox, eqcoor,
				 srch->epoch,epoch,
				 &ra0,&dec0,&srch->rapm,&srch->decpm);
			ListCat (ranges, eqout);
			}
		    }
		ctgclose (srchcat);
		}
	    }
	else {
	    if (strcmp (listfile,"STDIN")==0 || strcmp (listfile,"stdin")==0)
                fd = stdin;
            else
                fd = fopen (listfile, "r");
            if (fd != NULL) {
                while (fgets (line, 200, fd)) {
		    blank = strchr (line, ' ');
		    if (blank)
			*blank = (char) 0;
		    blank = strchr (line, '\n');
		    if (blank)
			*blank = (char) 0;
		    ranges = (char *) calloc (strlen(line) + 1, 1);
		    strcpy (ranges, line);
		    ListCat (ranges, eqout);
		    }
		fclose (fd);
		}
	    }
	}
    else {
	if (sysout0 && !syscoor)
	    syscoor = sysout0;
	if (syscoor) {
	    if (epoch0 == 0.0) {
		if (syscoor == WCS_B1950)
		    epoch0 = 1950.0;
		else
		    epoch0 = 2000.0;
		}
	    if (eqout == 0.0) {
		if (syscoor == WCS_B1950)
		    eqout = 1950.0;
		else
		    eqout = 2000.0;
		}
	    }
	ListCat (ranges, eqout);
	}

    for (i = 0; i < ncat; i++) {
	free (refcatname[i]);
	ctgclose (starcat[i]);
	}

    /* Free memory used for search results and return */
    if (gx) free ((char *)gx);
    if (gy) free ((char *)gy);
    if (gm) free ((char *)gm);
    if (gmb) free ((char *)gmb);
    if (gra) free ((char *)gra);
    if (gdec) free ((char *)gdec);
    if (gnum) free ((char *)gnum);
    if (gc) free ((char *)gc);
    if (gobj) free ((char *)gobj);

    return (0);
}

static void
usage (progname, command)
char *progname;
char *command;
{
    FILE *dev;

    if (http)
	dev = stdout;
    else
	dev = stderr;
    if (version)
	exit (-1);
    if (strsrch (progname,"gsc2") != NULL)
	fprintf (dev,"Find GSC II Stars in a region on the sky\n");
    else if (strsrch (progname,"gsca") != NULL)
	fprintf (dev,"Find GSC-ACT Stars in a region on the sky\n");
    else if (strsrch (progname,"gsc") != NULL)
	fprintf (dev,"Find HST Guide Stars in a region on the sky\n");
    else if (strsrch (progname,"tmc") != NULL ||
    	strsrch (progname,"2mp") != NULL)
	fprintf (dev,"Find 2MASS Point Sources in a region on the sky\n");
    else if (strsrch (progname,"ujc") != NULL)
	fprintf (dev,"Find USNO J Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"uac") != NULL)
	fprintf (dev,"Find USNO A Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"ua1") != NULL)
	fprintf (dev,"Find USNO A-1.0 Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"ua2") != NULL)
	fprintf (dev,"Find USNO A-2.0 Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"usac") != NULL)
	fprintf (dev,"Find USNO SA Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"usa1") != NULL)
	fprintf (dev,"Find USNO SA-1.0 Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"usa2") != NULL)
	fprintf (dev,"Find USNO SA-2.0 Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"tmc") != NULL)
	fprintf (dev,"Find 2MASS Point Source Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"sao") != NULL)
	fprintf (dev,"Find SAO Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"ppm") != NULL)
	fprintf (dev,"Find PPM Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"ira") != NULL)
	fprintf (dev,"Find IRAS Point Sources in a region on the sky\n");
    else if (strsrch (progname,"tyc") != NULL)
	fprintf (dev,"Find Tycho Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"hip") != NULL)
	fprintf (dev,"Find Hipparcos Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"act") != NULL)
	fprintf (dev,"Find ACT Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"bsc") != NULL)
	fprintf (dev,"Find Bright Star Catalog stars in a region on the sky\n");
    else
	fprintf (dev,"Find catalog stars in a region on the sky\n");

    if (command)
	fprintf (dev, "Bad command: %s\n", command);
    fprintf (dev,"Usage: %s [arguments] ra dec system (J2000, B1950, etc.)\n",
	progname);
    fprintf (dev,"  or : %s [arguments] list of catalog number ranges\n",
	progname);
    fprintf (dev,"  or : %s [arguments] @file of either positions or numbers)\n",
	progname);
    fprintf(dev,"  -a: List single closest catalog source\n");
    fprintf(dev,"  -b: Output B1950 (FK4) coordinates\n");
    if (!strcmp (progname, "scat"))
    fprintf(dev,"  -c name: Reference catalog (act, gsc, ua2, usa2, or local file\n");
    fprintf(dev,"  -d: Output RA and Dec in degrees instead of hms dms\n");
    fprintf(dev,"  -e: Output ecliptic coordinates\n");
    fprintf(dev,"  -f: Output search center for other programs\n");
    fprintf(dev,"  -g: Output galactic coordinates\n");
    fprintf(dev,"  -h: Print heading, else do not \n");
    fprintf(dev,"  -i: Print catalog object name, not catalog number\n");
    fprintf(dev,"  -j: Output J2000 (FK5) coordinates\n");
    fprintf(dev,"  -k kwd: Add this keyword to output from tab table search\n");
    fprintf(dev,"  -l: Print center and closest star on one line\n");
    fprintf(dev,"  -m mag1 [mag2]: Magnitude limit(s)\n");
    fprintf(dev,"  -n num: Number of brightest stars to print \n");
    fprintf(dev,"  -o name: Object name \n");
    fprintf(dev,"  -q year: Equinox of output positions in years\n");
    fprintf(dev,"  -r rad [dy]: Search radius (<0=-half-width) in arcsec\n");
    fprintf(dev,"  -s d|m|n|p|r: Sort by r=RA d=Dec m=Mag n=none p=distance\n");
    fprintf(dev,"  -t: Tab table to standard output as well as file\n");
    fprintf(dev,"  -u x y: Print x y instead of number in front of non-tab entry\n");
    fprintf(dev,"  -v: Verbose\n");
    fprintf(dev,"  -w: Write output file search[objname].[catalog]\n");
    if (!strcmp (progname, "scat") || !strcmp (progname, "sgsc"))
	fprintf(dev,"  -x type: GSC object type (0=stars 3=galaxies -1=all -2=bands)\n");
    fprintf(dev,"  -y year: Epoch of output positions in years\n");
    fprintf(dev,"  -z: Append to output file search[objname].[catalog]\n");
    exit (1);
}

#define TABMAX 64

static int
ListCat (ranges, eqout)

char	*ranges;	/* String with range of catalog numbers to list */
double	eqout;		/* Equinox for output coordinates */

{
    double cra, cdec;
    double epout = 0.0;
    int sysref;		/* Coordinate system of reference catalog */
    double eqref;	/* Equinox of catalog to be searched */
    double epref;	/* Epoch of catalog to be searched */
    int ng;		/* Number of catalog stars */
    int ns;		/* Number of brightest catalog stars actually used */
    struct Range *range; /* Range of catalog numbers to list */
    int nfind;		/* Number of stars to find */
    int i, ngmax, mprop, nc;
    double das, dds, drs;
    double gnmax;
    int degout;
    double xmag, xmag1;
    double maxnum;
    FILE *fd;
    char rastr[32], decstr[32];	/* coordinate strings */
    char numstr[32];	/* Catalog number */
    char cstr[32];	/* Coordinate system */
    char *catalog;
    double drad, dra, ddec, mag1, mag2;
    double gdist, da, dd, dec, gdmax;
    int nlog;
    int nmag;
    int typecol;
    int band;
    int sysout;
    char headline[160];
    char filename[80];
    char title[80];
    char string[TABMAX];
    char temp[80];
    char isp[4];
    int ngsc;
    int refcat;		/* reference catalog switch */
    int icat, nndec, nnfld, nsfld;
    double date, time;
    int gcset;
    int ndist;
    int distsort;
    double flux1, flux2, flux3, flux4;
    char lim1, lim2, lim3, lim4;
    double pra, pdec;
    void ep2dt();
    void PrintNum();
    int LenNum();

    /* Drop out if no catalog is specified */
    if (ncat < 1) {
	fprintf (stderr, "No catalog specified\n");
	exit (-1);
	}

    /* Allocate space for returned catalog information */
    if (ranges != NULL) {
	int nfdef = 9;

	/* Allocate and fill list of numbers to read */
	range = RangeInit (ranges, nfdef);
	ngmax = rgetn (range) * 4;
	}
    else if (nstars != 0)
	ngmax = nstars;
    else
	ngmax = MAXCAT;

    if (ngmax > nalloc) {

	/* Free currently allocated buffers if more entries are needed */
	if (nalloc > 0) {
	    if (gm) free ((char *)gm);
	    if (gmb) free ((char *)gmb);
	    if (gra) free ((char *)gra);
	    if (gdec) free ((char *)gdec);
	    if (gpra) free ((char *)gpra);
	    if (gpdec) free ((char *)gpdec);
	    if (gnum) free ((char *)gnum);
	    if (gc) free ((char *)gc);
	    if (gx) free ((char *)gx);
	    if (gy) free ((char *)gy);
	    if (gobj) free ((char *)gobj);
	    }
	gm = NULL;
	gmb = NULL;
	gra = NULL;
	gdec = NULL;
	gpra = NULL;
	gpdec = NULL;
	gnum = NULL;
	gc = NULL;
	gx = NULL;
	gy = NULL;
	gobj = NULL;

	if (!(gnum = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gnum\n",
		    ngmax*sizeof(double));
	if (!(gra = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gra\n",
		     ngmax*sizeof(double));
	if (!(gdec = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gdec\n",
		     ngmax*sizeof(double));
	if (!(gm = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gm\n",
		     ngmax*sizeof(double));
	if (!(gmb = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gmb\n",
		     ngmax*sizeof(double));
	if (!(gc = (int *) calloc (ngmax, sizeof(int))))
	    fprintf (stderr, "Could not calloc %d bytes for gc\n",
		     ngmax*sizeof(int));
	if (!(gx = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gx\n",
		     ngmax*sizeof(double));
	if (!(gy = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gy\n",
		     ngmax*sizeof(double));
	if (!(gobj = (char **) calloc (ngmax, sizeof(char *))))
	    fprintf (stderr, "Could not calloc %d bytes for obj\n",
		     ngmax*sizeof(char *));
	else {
	    for (i = 0; i < ngmax; i++)
		gobj[i] = NULL;
	    }
	if (!(gpra = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gpra\n",
		     ngmax*sizeof(double));
	if (!(gpdec = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gpdec\n",
		     ngmax*sizeof(double));
	if (gnum==NULL || gra==NULL || gdec==NULL || gm==NULL || gmb==NULL ||
	    gc==NULL || gx==NULL || gy==NULL || gobj==NULL || gpra == NULL ||
	    gpdec == NULL) {
	    if (gm) free ((char *)gm);
	    gm = NULL;
	    if (gmb) free ((char *)gmb);
	    gmb = NULL;
	    if (gra) free ((char *)gra);
	    gra = NULL;
	    if (gdec) free ((char *)gdec);
	    gdec = NULL;
	    if (gpra) free ((char *)gpra);
	    gpra = NULL;
	    if (gpdec) free ((char *)gpdec);
	    gpdec = NULL;
	    if (gnum) free ((char *)gnum);
	    gnum = NULL;
	    if (gc) free ((char *)gc);
	    gc = NULL;
	    if (gx) free ((char *)gx);
	    gx = NULL;
	    if (gy) free ((char *)gy);
	    gy = NULL;
	    if (gobj) free ((char *)gobj);
	    gobj = NULL;
	    nalloc = 0;
	    return (0);
	    }

	/* Initialize catalog entry values */
	for (i = 0; i < ngmax; i++) {
	    gm[i] = 99.0;
	    gmb[i] = 99.0;
	    gra[i] = 0.0;
	    gdec[i] = 0.0;
	    gpra[i] = 0.0;
	    gpdec[i] = 0.0;
	    gnum[i] = 0.0;
	    gc[i] = 0;
	    gx[i] = 0.0;
	    gy[i] = 0.0;
	    }
	nalloc = ngmax;
	}

    /* Start of per catalog loop */
    for (icat = 0; icat < ncat; icat++) {
	nndec = 0;
	isp[2] = (char) 0;
	isp[3] = (char) 0;

	/* Skip this catalog if no name is given */
	if (refcatname[icat] == NULL || strlen (refcatname[icat]) == 0) {
	    fprintf (stderr, "Catalog %d not specified\n", icat);
	    continue;
	    }

        if (verbose || (printhead && notprinted)) {
	    if (closest)
	    else
	    }

	if (debug)
	    nlog = 1;
	else if (verbose) {
	    if (refcat == UAC  || refcat == UA1  || refcat == UA2 ||
		refcat == USAC || refcat == USA1 || refcat == USA2 ||
		refcat == GSC  || refcat == GSCACT || refcat == TMPSC)
		nlog = 1000;
	    else
		nlog = 100;
	    }
	else
	    nlog = 0;

	/* Figure out which catalog we are searching */
	if (!(refcat = RefCat (refcatname[icat],title,&sysref,&eqref,
			       &epref,&mprop,&nmag))) {
	    fprintf (stderr,"ListCat: Catalog '%s' is missing\n", refcatname[icat]);
	    return (0);
	    }

	/* Set output coordinate system from command line or catalog */
	if (sysout0)
	    sysout = sysout0;
	else if (srch!= NULL && srch->epoch != 0.0)
	    sysout = srch->coorsys;
	if (!sysout)
	    sysout = sysref;
	if (!sysout)
	    sysout = WCS_J2000;

	/* Set epoch from command line, search catalog, or searched catalog */
	epout = epoch0;
	if (epout == 0.0) {
	    if (srch!= NULL && srch->epoch != 0.0)
		epout = srch->epoch;
	    else if (sysout0 == WCS_J2000)
		epout = 2000.0;
	    else if (sysout0 == WCS_B1950)
		epout = 1950.0;
	    else if (!mprop)
		epout = epref;
	    else {
		if (sysout == WCS_B1950)
		    epout = 1950.0;
		else
		    epout = 2000.0;
		}
	    }

	/* Set equinox from command line, search catalog, or searched catalog */
	if (eqout == 0.0) {
	    if (srch!= NULL && srch->equinox != 0.0)
		eqout = srch->equinox;
	    else if (sysout0 == WCS_J2000)
		eqout = 2000.0;
	    else if (sysout0 == WCS_B1950)
		eqout = 1950.0;
	    else
		eqout = epref;
	    if (eqout == 0.0)
		eqout = 2000.0;
	    }
	if (epout == 0.0)
	    epout = eqout;

	/* Set degree flag for output */
	if (sysout == WCS_ECLIPTIC || sysout == WCS_GALACTIC)
	    degout = 1;
	else
	    degout = degout0;

	/* Find stars specified by number */
	if (ranges != NULL) {
	    nfind = rgetn (range);
	    for (i = 0; i < nfind; i++)
		gnum[i] = rgetr8 (range);
	    wfile = 0;

	    /* Find the specified catalog stars */
	    ng = ctgrnum (refcatname[icat], refcat,
		      nfind, sysout, eqout, epout, match, &starcat[icat],
		      gnum, gra, gdec, gpra, gpdec, gm, gmb, gc, gobj, debug);

	    if (gobj[0] == NULL)
		gobj1 = NULL;
	    else
		gobj1 = gobj;
	    /* if (ng > nfind)
		ns = nfind;
	    else */
		ns = ng;

	    for (i = 0; i < ns; i++ ) {
		gx[i] = 0.0;
		gy[i] = 1.0;
		}

	    /* Find largest catalog number printed */
	    maxnum = 0.0;
	    for (i = 0; i < ns; i++ ) {
		if (gnum[i] > maxnum)
		    maxnum = gnum[i];
		}
	    nnfld = CatNumLen (refcat, maxnum, nndec);

	    /* Check to see whether gc is set at all */
	    gcset = 0;
	    for (i = 0; i < ns; i++ ) {
		if (gc[i] != 0) {
		    gcset = 1;
		    break;
		    }
		}

	    /* Set flag for plate, class, or type column */
	    if (refcat == BINCAT || refcat == SAO  || refcat == PPM ||
		refcat == BSC)
		typecol = 1;
	    else if ((refcat == GSC || refcat == GSCACT) && classd < -1)
		typecol = 3;
	    else if (refcat == GSC || refcat == GSCACT ||
		refcat == UJC ||
		refcat == USAC || refcat == USA1   || refcat == USA2 ||
		refcat == UAC  || refcat == UA1    || refcat == UA2 ||
		refcat == BSC  || (refcat == TABCAT&&gcset))
		typecol = 2;
	    else
		typecol = 0;

	    /* Write out entries for use as image centers */
	    if (searchcenter) {
		gnmax = gnum[0];
		for (i = 1; i < ns; i++ ) {
		    if (gnum[i] > gnmax) gnmax = gnum[i];
		    }
		for (i = 0; i < ns; i++ ) {
		    if (sysref == WCS_XY) {
			num2str (rastr, gra[i], 10, 5);
			num2str (decstr, gdec[i], 10, 5);
			}
		    else if (degout) {
			deg2str (rastr, 32, gra[i], 6);
			deg2str (decstr, 32, gdec[i], 6);
			}
		    else {
			ra2str (rastr, 32, gra[i], 3);
			dec2str (decstr, 32, gdec[i], 2);
			}
		    wcscstr (cstr, sysout, eqout, epout);
		    if (printobj && gobj1 != NULL)
			printf ("%s %s %s %s\n",
				gobj[i], rastr, decstr, cstr);
		    else {
			CatNum (refcat, -nnfld, nndec, gnum[i], numstr);
			printf ("%s_%s %s %s %s\n",
				refcatname[icat],numstr,rastr,decstr,cstr);
			}
		    }
		return (ns);
		}
	    }

	/* Find stars specified by location */
	else {

	    /* Set search radius if finding closest star */
	    if (rad0 == 0.0 && dra0 == 0.0) {
		if (closest) {
		    if (refcat == GSC || refcat == GSCACT || refcat == UJC ||
			refcat == USAC || refcat == USA1 || refcat == USA2)
			rad0 = 900.0;
		    else if ( refcat == UAC  || refcat == UA1  || refcat == UA2)
			rad0 = 120.0;
		    else if ( refcat == GSC2 || refcat == TMPSC)
			rad0 = 120.0;
		    else if (refcat == SAO || refcat == PPM || refcat == IRAS)
			rad0 = 5000.0;
		    else
			rad0 = 1800.0;
		    }
		else
		    rad0 = 10.0;
		}

	    /* Set limits from defaults and command line information */
	    if (GetArea (verbose,syscoor,sysout,eqout,epout,
		     &cra,&cdec,&dra,&ddec,&drad))
		return (0);

	    if (srch != NULL) {
		if ((srchcat->stnum == 0 || srchcat->stnum > 4) &&
		    strlen (srch->objname) > 0) {
		    if (objname == NULL)
			objname = (char *) malloc (32);
		    strcpy (objname, srch->objname);
		    }
		else {
		    if (objname == NULL)
			objname = (char *) calloc (1,32);
		    nsfld = CatNumLen (TXTCAT, srch->num, srchcat->nndec);
		    CatNum (TXTCAT, nsfld, srchcat->nndec, srch->num, objname);
		    }
		}
	    nnfld = CatNumLen (refcat, 0.0, nndec);

	    /* Print search center and size in input and output coordinates */
	    if (verbose || (printhead && !oneline)) {
		SearchHead (icat,sysout,sysout,eqout,eqout,epout,epout,
			cra,cdec,dra,ddec,drad,nnfld);
		if (!closest) {
		    if (sysout != syscoor)
			SearchHead (icat,sysout,syscoor,eqout,eqcoor,
				    epout,epout,cra,cdec,dra,ddec,drad,nnfld);
		    if (sysref != syscoor && sysref != sysout)
			SearchHead (icat,sysout,sysref,eqout,eqref,
				    epout,epref,cra,cdec,dra,ddec,drad,nnfld);
		    }
		}

	    /* Set the magnitude limits for the catalog search */
	    if (maglim2 == 0.0) {
		mag1 = 0.0;
		mag2 = 0.0;
		}
	    else {
		mag1 = maglim1;
		mag2 = maglim2;
		}

	    /* Find the nearby reference stars, in ra/dec */
	    if (catsort == SORT_DIST)
		distsort = 1;
	    else
		distsort = 0;
	    ng = ctgread (refcatname[icat], refcat, distsort,
		      cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,mag2,
		      ngmax,&starcat[icat],
		      gnum,gra,gdec,gpra,gpdec,gm,gmb,gc,gobj,nlog);
	    if (ngmax < 1)
		return (ng);

	    /* Set flag if any proper motions are non-zero
	    mprop = 0;
	    for (i = 0; i < ng; i++) {
		if (gpra[i] != 0.0 || gpdec[i] != 0.0) {
		    mprop = 1;
		    break;
		    }
		} */

	    if ((verbose || printhead) && ncat == 1 && ng < 1) {
		fprintf (stderr, "No stars found in %s\n",refcatname[icat]);
		return (0);
		}

	    if (gobj[0] == NULL)
		gobj1 = NULL;
	    else
		gobj1 = gobj;
	    if (ng > ngmax)
		ns = ngmax;
	    else
		ns = ng;

	    /* Find largest catalog number to be printed */
	    maxnum = 0.0;
	    for (i = 0; i < ns; i++ ) {
		if (gnum[i] > maxnum)
		    maxnum = gnum[i];
		}
	    nnfld = CatNumLen (refcat, maxnum, nndec);

	    /* Check to see whether gc is set at all */
	    gcset = 0;
	    for (i = 0; i < ng; i++ ) {
		if (gc[i] != 0) {
		    gcset = 1;
		    break;
		    }
		}

	    /* Set flag for plate, class, or type column */
	    if (refcat == BINCAT || refcat == SAO  || refcat == PPM ||
		refcat == BSC)
		typecol = 1;
	    else if ((refcat == GSC || refcat == GSCACT) && classd < -1)
		typecol = 3;
	    else if (refcat == GSC || refcat == GSCACT ||
		refcat == UJC || 
		refcat == USAC || refcat == USA1   || refcat == USA2 ||
		refcat == UAC  || refcat == UA1    || refcat == UA2 ||
		refcat == BSC  || (refcat == TABCAT&&gcset))
		typecol = 2;
	    else
		typecol = 0;

	    /* Compute distance from star to search center */
	    for (i = 0; i < ns; i++ ) {
		gx[i] = wcsdist (cra, cdec, gra[i], gdec[i]);
		gy[i] = 1.0;
		}

	    if (ns > 1) {

		/* Sort reference stars from closest to furthest */
		if (catsort == SORT_DIST)
		    XSortStars (gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gmb,gc,gobj1,ns);

		/* Sort star-like objects in image by right ascension */
		else if (catsort == SORT_RA)
		    RASortStars (gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gmb,gc,gobj1,ns);

		/* Sort star-like objects in image by declination */
		else if (catsort == SORT_DEC)
		    DecSortStars (gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gmb,gc,gobj1,ns);

		/* Sort reference stars from brightest to faintest */
		else if (catsort == SORT_MAG)
		    MagSortStars (gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gmb,gc,gobj1,ns);
		}

	    /* Print one line with search center and found star */
	    if (oneline) {
		if (ns > 0) {
		    das = dra * cos (degrad(cdec)) * 3600.0;
		    dds = ddec * 3600.0;
		    drs = drad * 3600.0;
		    if (drs <= 0.0)
			drs = sqrt (das*das + dds*dds);
		    if (das <= 0.0) {
			das = drs;
			dds = drs;
			}
		    if (nohead & tabout) {

			/* Write tab table heading */
			catalog = CatName (refcat, refcatname[icat]);
			printf ("catalog	%s\n", catalog);
			if (listfile != NULL)
			    printf ("search	%s\n", listfile);
			printf ("equinox	%.1f\n", eqout);
			printf ("epoch	%.1f\n", epout);
			if (dra0 > 0.0) {
			    printf ("drasec	%.2f\n", dra0);
			    printf ("ddecsec	%.2f\n", ddec0);
			    }
			else if (rad0 > 0)
			    printf ("radsec	%.2f\n", drs);
			else if (rad0 < 0)
			    printf ("boxsec	%.2f\n", dds);
			if (mprop) {
			    printf ("pmunit	mas/yr\n");
			    }

		        /* Write column headings */
		        if (srch != NULL) {
			    if (srchcat->keyid[0] > 0) {
				printf ("%s", srchcat->keyid);
				nc = strlen (srchcat->keyid);
				}
			    else {
				printf ("srch_id");
				nc = 7;
				}
			    if (srchcat->nnfld > nc) {
				for (i = nc; i < srchcat->nnfld; i++)
				    printf (" ");
				}
			    printf ("	");
			    }
		        printf ("srch_ra     	srch_dec    	");
		        if (srchcat != NULL) {
			    if (srchcat->nmag > 1) {
				if (strlen (srchcat->keymag1) >0)
				    printf ("%s	", srchcat->keymag1);
				else
				    printf ("smag1	");
				if (strlen (srchcat->keymag2) >0)
				    printf ("%s	", srchcat->keymag2);
				else
				    printf ("smag2	");
				if (srchcat->nmag > 2) {
				    if (strlen (srchcat->keymag3) >0)
					printf ("%s	", srchcat->keymag3);
				    else
					printf ("smag3	");
				    }
				}
			    else if (srchcat->nmag > 0) {
				if (strlen (srchcat->keymag1) >0)
				    printf ("%s	", srchcat->keymag1);
				else
				    printf ("s_mag	");
				}
			    if (srchcat->nepoch)
				printf ("epoch    	");
			    if (srchcat->sptype) {
				if (strlen (srch->isp) > 2) {
				    padspt = 1;
				    printf ("srch_spt       	");
				    }
				else {
				    padspt = 0;
				    printf ("sspt	");
				    }
				}
			    }
			if (refcat == TABCAT && starcat[icat]->keyid[0] >0) {
			    strcpy (headline, starcat[icat]->keyid);
			    strcat (headline, "                ");
			    }
			else
			    CatID (headline, refcat);
			headline[nnfld] = (char) 0;
			printf ("%s", headline);
			printf ("	ra          	dec        	");
			if (refcat == TABCAT) {
			    if (starcat[icat]->nmag > 1) {
				if (strlen (starcat[icat]->keymag1) >0)
				    printf ("%s	", starcat[icat]->keymag1);
				else
				    printf ("mag1	");
				if (strlen (starcat[icat]->keymag2) > 0)
				    printf ("%s	", starcat[icat]->keymag2);
				else
				    printf ("mag2	");
				if (starcat[icat]->nmag > 2) {
				    if (strlen (starcat[icat]->keymag3) >0)
				        printf ("%s	", starcat[icat]->keymag3);
				    else
					printf ("mag3	");
				    }
				}
			    else if (starcat[icat]->nmag > 0) {
				if (starcat[icat]->keymag1[0] >0)
				    printf ("%s	", starcat[icat]->keymag1);
				else
				    printf ("mag	");
				}
			    }
			else if (refcat == GSC2)
			    printf ("magf 	magj 	magv 	magn 	");
			else if (refcat == HIP)
			    printf ("magb 	magv 	parlx	parer	");
			else if (refcat == IRAS)
			    printf ("f10m 	f25m 	f60m 	f100m	");
			else if (refcat == TMPSC)
			    printf ("magj  	magh  	magk  	");
			else if (nmag > 1)
			    printf ("magb  	magr  	");
			else
			    printf ("mag   	");
			if (typecol == 1)
			    printf ("spt   	");
			printf ("n    	");
			printf ("dra");
			for (i = 3; i < LenNum(das,2); i++)
			    printf (" ");
			printf ("	ddec");
			for (i = 4; i < LenNum(dds,2); i++)
			    printf (" ");
			printf ("	drad");
			for (i = 4; i < LenNum(drs,2); i++)
			    printf (" ");
			printf ("\n");
			if (srch != NULL) {
			    printf ("-------");
			    if (srchcat->nnfld > 7) {
				for (i = 8; i < srchcat->nnfld; i++)
				    printf ("-");
				}
			    printf ("	");
			    }
			printf ("------------	------------	");
			if (srchcat != NULL) {
			    if (srchcat->nmag > 0)
				printf ("-----	");
			    if (srchcat->nmag > 1)
				printf ("-----	");
			    if (srchcat->nmag > 2)
				printf ("-----	");
			    if (srchcat->nepoch)
				printf ("---------	");
			    if (srchcat->sptype != 0) {
				if (padspt)
				    printf ("---------------	");
				else
				    printf ("----	");
				}
			    }
			strcpy (headline,"--------------------");
			headline[nnfld] = (char) 0;
			printf ("%s", headline);
			printf ("	------------	------------	");
			if (refcat == GSC2)
			    printf ("-----	-----	-----	-----	");
			else if (refcat == HIP || refcat == IRAS)
			    printf ("-----	-----	-----	-----	");
			else if (refcat == TMPSC)
			    printf ("------	------	------	");
			else if (nmag > 1)
			    printf ("-----	-----	");
			else
			    printf ("-----	");
			if (typecol == 1)
			    printf ("---   	");
			printf ("--	");
			for (i = 0; i < LenNum(das,2); i++)
			    printf ("-");
			printf ("	");
			for (i = 0; i < LenNum(dds,2); i++)
			    printf ("-");
			printf ("	");
			for (i = 0; i < LenNum(drs,2); i++)
			    printf ("-");
			printf ("\n");
			nohead = 0;
			}
		    if (srch != NULL) {
			if ((srchcat->stnum == 0 || srchcat->stnum > 4) &&
			    strlen (srch->objname) > 0)
			    strcat (numstr, srch->objname);
			else
			    CatNum (TXTCAT,-srchcat->nnfld,srchcat->nndec,
				    srch->num,numstr);
		        if (tabout)
			    printf ("%s	", numstr);
		        else
			    printf ("%s ", numstr);
			}
		    if (degout) {
			num2str (rastr, cra, 10, 5);
			num2str (decstr, cdec, 10, 5);
			}
		    else {
			ra2str (rastr, 32, cra, 3);
			dec2str (decstr, 32, cdec, 2);
			}
		    if (tabout)
			printf ("%s	%s", rastr, decstr);
		    else
			printf ("%s %s", rastr, decstr);
		    if (srch != NULL) {
			if (srchcat->nmag > 0) {
			    if (tabout)
				printf ("	%5.2f", srch->xmag[0]);
			    else
				printf (" %5.2f", srch->xmag[0]);
			    }
			if (srchcat->nmag > 1) {
			    if (tabout)
				printf ("	%5.2f", srch->xmag[1]);
			    else
				printf (" %5.2f", srch->xmag[1]);
			    }
			if (srchcat->nmag > 2) {
			    if (tabout)
				printf ("	%5.2f", srch->xmag[2]);
			    else
				printf (" %5.2f", srch->xmag[2]);
			    }
			if (srchcat->nepoch) {
			    ep2dt (srch->epoch, &date, &time);
			    if (tabout)
				printf ("	%9.4f", date);
			    else
				printf (" %9.4f", date);
			    }
			if (srchcat->sptype != 0) {
			    if (padspt) {
				for (i = 0; i < 14; i++) {
				    if (srch->isp[i] == (char) 0)
					srch->isp[i] = ' ';
				    }
				srch->isp[14] = (char) 0;
				if (tabout)
				    printf ("	%14s", srch->isp);
				else
				    printf ("  %14s ", srch->isp);
				}
			    else if (tabout)
				printf ("	%s", srch->isp);
			    else
				printf ("  %s ", srch->isp);
			    }
			}
		    if (gobj1 != NULL) {
			if (strlen (gobj1[0]) > 0)
			    strcpy (numstr, gobj1[0]);
			}
		    else if (starcat[icat] != NULL)
			CatNum (refcat,-nnfld,starcat[icat]->nndec,gnum[0],numstr);
		    else
			CatNum (refcat,-nnfld,nndec,gnum[0],numstr);
		    if (degout) {
			num2str (rastr, gra[0], 10, 5);
			num2str (decstr, gdec[0], 10, 5);
			}
		    else {
			ra2str (rastr, 32, gra[0], 3);
			dec2str (decstr, 32, gdec[0], 2);
			}
		    if (tabout)
			printf ("	%s	%s	%s",
			        numstr, rastr, decstr);
		    else
			printf (" %s %s %s",
			        numstr, rastr, decstr);
		    if (refcat == GSC2 || refcat == HIP) {
			xmag = 0.01 * (double) (gc[0] / 10000);
			xmag1 = 0.01 * (double) (gc[0] % 10000);
			if (tabout)
			    printf ("	%5.2f	%5.2f	%5.2f	%5.2f",
				    gmb[0], gm[0], xmag, xmag1);
			else
			    printf (" %5.2f %5.2f %5.2f %5.2f",
				    gmb[0], gm[0], xmag, xmag1);
			}
		    else if (refcat == IRAS) {
			xmag = 0.01 * (double) (gc[0] / 100000);
			xmag1 = 0.01 * (double) (gc[0] % 100000);
			if (gmb[0] > 100.0) {
			    gmb[0] = gmb[0] - 100.0;
			    lim1 = 'L';
			    }
			else
			    lim1 = ' ';
			if (gm[0] > 100.0) {
			    gm[0] = gm[0] - 100.0;
			    lim2 = 'L';
			    }
			else
			    lim2 = ' ';
			if (xmag > 100.0) {
			    xmag = xmag - 100.0;
			    lim3 = 'L';
			    }
			else
			    lim3 = ' ';
			if (xmag1 > 100.0) {
			    xmag1 = xmag1 - 100.0;
			    lim4 = 'L';
			    }
			else
			    lim4 = ' ';
			flux1 = 1000.0 * pow (10.0, -gmb[0] / 2.5);
			flux2 = 1000.0 * pow (10.0, -gm[0] / 2.5);
			flux3 = 1000.0 * pow (10.0, -xmag / 2.5);
			flux4 = 1000.0 * pow (10.0, -xmag1 / 2.5);
			if (tabout)
			    printf ("	%.2f%c	%.2f%c	%.2f%c	%.2f%c",
				   flux1,lim1,flux2,lim2,flux3,lim3,flux4,lim4);
			else
			    printf (" %5.2f%c %5.2f%c %5.2f%c %5.2f%c",
				   flux1,lim1,flux2,lim2,flux3,lim3,flux4,lim4);
			}
		    else if (refcat == TMPSC) {
			xmag = 0.001 * (double) gc[0];
			if (tabout)
			    printf ("	%6.3f	%6.3f	%6.3f",
				    gm[0], gmb[0], xmag);
			else
			    printf (" %6.3f %6.3f %6.3f",
				    gm[0], gmb[0], xmag);
			}
		    else if (nmag > 1) {
			if (tabout)
			    printf ("	%5.2f	%5.2f", gmb[0], gm[0]);
			else
			    printf (" %5.2f %5.2f", gmb[0], gm[0]);
			}
		    else {
			if (tabout)
			    printf ("	%5.2f", gm[0]);
			else
			    printf (" %5.2f", gm[0]);
			}
		    if (typecol == 1) {
			isp[0] = gc[0] / 1000;
			isp[1] = gc[0] % 1000;
			if (isp[0] == ' ' && isp[1] == ' ') {
			    isp[0] = '_';
			    isp[1] = '_';
			    }
			if (tabout)
			    printf ("	%2s ", isp);
			else
			    printf (" %2s ", isp);
			}
		    if (tabout)
			printf ("	%d", ng);
		    else
			printf (" %d", ng);
		    dec = (gdec[0] + cdec) * 0.5;
		    if (degout) {
			if ((gra[0] - cra) > 180.0)
			    da = gra[0] - cra - 360.0;
			else if ((gra[0] - cra) < -180.0)
			    da = gra[0] - cra + 360.0;
			else
			    da = gra[0] - cra;
			dd = gdec[0] - cdec;
			gdist = sqrt (da*da + dd*dd);
			ndist = 5;
			}
		    else {
			if ((gra[0] - cra) > 180.0)
			    da = 3600.0*(gra[0]-cra-360.0)*cos(degrad(dec));
			else if ((gra[0] - cra) < -180.0)
			    da = 3600.0*(gra[0]+360.0-cra)*cos(degrad(dec));
			else
			    da = 3600.0 * (gra[0] - cra) * cos (degrad (dec));
			dd = 3600.0 * (gdec[0] - cdec);
			gdist = 3600.0 * gx[0];
			ndist = 2;
			}
		    if (tabout)
			printf ("	");
		    else
			printf (" ");
		    PrintNum (das, da, ndist);
		    if (tabout)
			printf ("	");
		    else
			printf (" ");
		    PrintNum (dds, dd, ndist);
		    if (tabout)
			printf ("	");
		    else
			printf (" ");
		    PrintNum (drs, gdist, ndist);
		    printf ("\n");
		    }
		notprinted = 0;
		continue;
		}

	    /* List the brightest or closest MAXSTARS reference stars */
	    if (ng > ngmax) {
		if ((verbose || printhead) && !closest) {
		    if (distsort) {
			if (ng > 1)
			    printf ("Closest %d / %d %s (closer than %.2f arcsec)",
				ns, ng, title, 3600.0*gx[ns-1]);
		        else
			    printf ("Closest of %d %s",ng, title);
			}
		    else if (maglim1 > 0.0)
			printf ("%d / %d %s (between %.2f and %.2f)",
			    ns, ng, title, gm[0], gm[ns-1]);
		    else
			printf ("%d / %d %s (brighter than %.2f)",
		 	    ns, ng, title, gm[ns-1]);
		    printf ("\n");
		    }
		}
	    else {
	        if (verbose || printhead) {
		    if (maglim1 > 0.0)
			printf ("%d %s between %.2f and %.2f\n",
			    ng, title, maglim1, maglim2);
		    else if (maglim2 > 0.0)
			printf ("%d %s brighter than %.2f\n",
			    ng, title, maglim2);
		    else if (verbose)
			printf ("%d %s\n", ng, title);
		    }
		}

	    /* Open result catalog file */
	    if (wfile && icat == 0) {
		if (objname)
		    strcpy (filename,objname);
		else
		    strcpy (filename,"search");
		for (i = 0; i < ncat; i++) {
		    strcat (filename,".");
		    strcat (filename,refcatname[icat]);
		    }
		if (printxy)
		    strcat (filename,".match");

		if (afile)
		    fd = fopen (filename, "a");
		else
		    fd = fopen (filename, "w");

		/* Free result arrays and return if cannot write file */
		if (fd == NULL) {
		    if (afile)
			fprintf (stderr, "%s:  cannot append to file %s\n",
				 cpname, filename);
		    else
			fprintf (stderr, "%s:  cannot write file %s\n",
				 cpname, filename);
        	    return (0);
		    }
		}
            }

    /* Write heading */
    if (tabout && nohead) {
	catalog = CatName (refcat, refcatname[icat]);
	sprintf (headline, "catalog	%s", catalog);
	if (wfile)
	    fprintf (fd, "%s\n", headline);
	else
	    printf ("%s\n", headline);

	if (!ranges) {
	    ra2str (rastr, 32, cra, 3);
	    if (wfile)
		fprintf (fd, "ra	%s\n", rastr);
	    else
		printf ("ra	%s\n", rastr);
	    dec2str (decstr, 32, cdec, 2);
	    if (wfile)
		fprintf (fd, "dec	%s\n", decstr);
	    else
		printf ("dec	%s\n", decstr);
	    }

	if (wfile)
	    fprintf (fd, "equinox	%.1f\n", eqout);
	else
	    printf ("equinox	%.1f\n", eqout);

	if (wfile)
	    fprintf (fd, "epoch	%.1f\n", epout);
	else
	    printf ("epoch	%.1f\n", epout);

	if (mprop) {
	    if (degout) {
		if (wfile)
		    fprintf (fd, "rpmunit	arcsec/century\n");
		else
		    printf ("rpmunit	arcsec/century\n");
		}
	    else {
		if (wfile)
		    fprintf (fd, "rpmunit	tsec/century\n");
		else
		    printf ("rpmunit	tsec/century\n");
		}
	    if (wfile)
		fprintf (fd, "dpmunit	arcsec/century\n");
	    else
		printf ("dpmunit	arcsec/century\n");
	    }

	das = dra * cos (degrad(cdec)) * 3600.0;
	dds = ddec * 3600.0;
	drs = drad * 3600.0;
	if (drs <= 0.0)
	    drs = sqrt (das*das + dds*dds);
	if (das <= 0.0) {
	    das = drs;
	    dds = drs;
	    }
 	if (dra0 > 0.0) {
	    if (wfile) {
		printf ("drasec	%.2f\n", dra0);
		printf ("ddecsec	%.2f\n", ddec0);
		}
	    else {
		fprintf (fd, "drasec	%.2f\n", dra0);
		fprintf (fd, "ddecsec	%.2f\n", ddec0);
		}
	    }
	else if (rad0 > 0) {
	    if (wfile)
		fprintf (fd, "radsec	%.2f\n", drs);
	    else
		printf ("radsec	%.2f\n", drs);
	    }
	else if (rad0 < 0) {
	    if (wfile)
		fprintf (fd, "boxsec	%.2f\n",dds);
	    else
		printf ("boxsec	%.2f\n",dds);
	    }

	if (catsort > 0) {
	    switch (catsort) {
		case SORT_DEC:
		    if (wfile)
			fprintf (fd, "catsort	dec\n");
		    else
			printf ("catsort	dec\n");
		    break;
		case SORT_DIST:
		    if (wfile)
			fprintf (fd, "catsort	dist\n");
		    else
			printf ("catsort	dist\n");
		    break;
		case SORT_MAG:
		    if (wfile)
			fprintf (fd, "catsort	mag\n");
		    else
			printf ("catsort	mag\n");
		    break;
		case SORT_RA:
		    if (wfile)
			fprintf (fd, "catsort	ra\n");
		    else
			printf ("catsort	ra\n");
		    break;
		case SORT_X:
		    if (wfile)
			fprintf (fd, "catsort	x\n");
		    else
			printf ("catsort	x\n");
		    break;
		default:
		    break;
		}
	    }

	if (wfile)
	else

	/* Print column headings */
	if (refcat == TABCAT && strlen(starcat[icat]->keyid) > 0)
	    sprintf (headline,"%s          ", starcat[icat]->keyid);
	else
	    CatID (headline, refcat);
	headline[nnfld] = (char) 0;

	if (sysout == WCS_GALACTIC)
	    strcat (headline,"	long_gal   	lat_gal  ");
	else if (sysout == WCS_ECLIPTIC)
	    strcat (headline,"	long_ecl   	lat_ecl  ");
	else if (sysout == WCS_B1950)
	    strcat (headline,"	ra1950      	dec1950  ");
	else
	    strcat (headline,"	ra      	dec      ");
	if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
	    refcat == UAC  || refcat == UA1  || refcat == UA2)
	    strcat (headline,"	magb	magr	plate");
	else if (refcat == TMPSC)
	    strcat (headline,"	magj	magh	magk");
	else if (refcat == GSC2)
	    strcat (headline,"	magf	magj	magv  	magn ");
	else if (refcat == IRAS)
	    strcat (headline,"	f10m  	f25m  	f60m   	f100m ");
	else if (refcat == HIP)
	    strcat (headline,"	magb	magv	parlx 	parer");
	else if (refcat==TYCHO || refcat==TYCHO2 || refcat==ACT)
	    strcat (headline,"	magb	magv");
	else if (refcat==HIP)
	    strcat (headline,"	magb 	magv 	prllx	parer");
	else if (refcat == GSC || refcat == GSCACT)
	    strcat (headline,"	mag	class	band	N");
	else if (refcat == UJC)
	    strcat (headline,"	mag	plate");
	else
	    strcat (headline,"	mag");
	if (typecol == 1)
	    strcat (headline,"	type");
	if (mprop)
	    strcat (headline,"	Ura    	Udec  ");
	if (ranges == NULL)
	    strcat (headline,"	arcsec");
	if (refcat == TABCAT && keyword != NULL) {
	    strcat (headline,"	");
	    strcat (headline, keyword);
	    }
	if (gobj1 != NULL)
	    strcat (headline,"	object");
	if (printxy)
	    strcat (headline, "	X      	Y      ");

	if (wfile)
	    fprintf (fd, "%s\n", headline);
	else
	    printf ("%s\n", headline);

	strcpy (headline, "---------------------");
	headline[nnfld] = (char) 0;
	strcat (headline,"	------------	------------");
	if (refcat == TMPSC)
	    strcat (headline,"	------	------	------");
	else if (refcat == IRAS)
	    strcat (headline,"	-----	-----	-----	-----");
	else if (refcat==HIP || refcat == GSC2)
	    strcat (headline,"	-----	-----	-----	-----");
	else if (nmag == 2)
	    strcat (headline,"	-----	-----");
	else
	    strcat (headline,"	-----");
	if (refcat == GSC || refcat == GSCACT)
	    strcat (headline,"	-----	----	-");
	else if (typecol == 1)
	    strcat (headline,"	----");
	else if (typecol == 2)
	    strcat (headline,"	-----");
	if (mprop)
	    strcat (headline,"	-------	------");
	if (ranges == NULL)
	    strcat (headline, "	------");
	if (refcat == TABCAT && keyword != NULL)
	    strcat (headline,"	------");
	if (printxy)
	    strcat (headline, "	-------	-------");

	if (wfile)
	    fprintf (fd, "%s\n", headline);
	else
	    printf ("%s\n", headline);
	nohead = 0;
	}

    else if (printhead && nohead) {
	if (ng == 0)
	    printf ("No %s Stars Found\n", title);
	else {
	    if (printxy)
		strcpy (headline, "  X     Y   ");
	    else {
		if (refcat == GSC || refcat == GSCACT)
		    strcpy (headline, "GSC number ");
		else if (refcat == GSC2)
		    strcpy (headline, "GSC II number ");
		else if (refcat == USAC)
		    strcpy (headline, "USNO SA number ");
		else if (refcat == USA1)
		    strcpy (headline, "USNO SA1 number");
		else if (refcat == USA2)
		    strcpy (headline, "USNO SA2 number");
		else if (refcat == UAC)
		    strcpy (headline, "USNO A number  ");
		else if (refcat == UA1)
		    strcpy (headline, "USNO A1 number ");
		else if (refcat == UA2)
		    strcpy (headline, "USNO A2 number ");
		else if (refcat == TMPSC)
		    strcpy (headline, "2MASS num. ");
		else if (refcat == UJC)
		    strcpy (headline, " UJ number    ");
		else if (refcat == SAO)
		    strcpy (headline, "SAO number ");
		else if (refcat == PPM)
		    strcpy (headline, "PPM number ");
		else if (refcat == BSC)
		    strcpy (headline, "BSC number ");
		else if (refcat == IRAS)
		    strcpy (headline, "IRAS num ");
		else if (refcat == TYCHO)
		    strcpy (headline, "Tycho number ");
		else if (refcat == TYCHO2)
		    strcpy (headline, "Tycho2 num  ");
		else if (refcat == HIP)
		    strcpy (headline, "Hip num ");
		else if (refcat == ACT)
		    strcpy (headline, "ACT number  ");
		else
		    strcpy (headline, "Number   ");
		}
	    if (sysout == WCS_B1950) {
		if (degout) {
		    if (eqout == 1950.0)
			strcat (headline, "  RA1950   Dec1950  ");
		    else {
			sprintf (temp, "RAB%7.2f DecB%7.2f  ", eqout, eqout);
			strcat (headline, temp);
			}
		    }
		else {
		    if (eqout == 1950.0)
			strcat (headline, "RAB1950      DecB1950    ");
		    else {
			sprintf (temp, "RAB%7.2f   DecB%7.2f  ", eqout, eqout);
			strcat (headline, temp);
			}
		    }
		}
	    else if (sysout == WCS_ECLIPTIC)
		strcat (headline, "Ecl Lon    Ecl Lat  ");
	    else if (sysout == WCS_GALACTIC)
		strcat (headline, "Gal Lon    Gal Lat  ");
	    else {
		if (degout) {
		    if (eqout == 2000.0)
			strcat (headline, "  RA2000   Dec2000  ");
		    else {
			sprintf (temp,"RAJ%7.2f  DecJ%7.2f ", eqout, eqout);
			strcat (headline, temp);
			}
		    }
		else {
		    if (eqout == 2000.0)
			strcat (headline, " RA2000       Dec2000   ");
		    else {
			sprintf (temp,"RAJ%7.2f   DecJ%7.2f  ", eqout, eqout);
			strcat (headline, temp);
			}
		    }
		}
	    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		refcat == UAC  || refcat == UA1  || refcat == UA2)
		strcat (headline, "  MagB  MagR Plate");
	    else if (refcat == UJC)
		strcat (headline, "  Mag  Plate");
	    else if (refcat == GSC || refcat == GSCACT)
		strcat (headline, "   Mag  Class Band N");
	    else if (refcat == SAO || refcat == PPM || refcat == BSC)
		strcat (headline, "   Mag  Type");
	    else if (refcat==TMPSC)
		strcat (headline, "   MagJ  MagH  MagK");
	    else if (refcat==IRAS)
		strcat (headline, "f10m   f25m   f60m   f100m");
	    else if (refcat==GSC2)
		strcat (headline, "   MagF  MagJ  MagV  MagN");
	    else if (refcat==HIP)
		strcat (headline, "  MagB  MagV  Parlx Parer");
	    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==ACT)
		strcat (headline, "  MagB   MagV ");
	    else if (refcat == TABCAT && gcset) {
		if (strlen(starcat[icat]->keymag1) > 0) {
		    strcat (headline," ");
		    strcat (headline, starcat[icat]->keymag1);
		    }
		else if (starcat[icat]->nmag > 1)
		    strcat (headline, " Mag1");
		else
		    strcat (headline, " Mag");
		if (starcat[icat]->nmag > 1) {
		    if (strlen(starcat[icat]->keymag2) > 0) {
			strcat (headline,"   ");
			strcat (headline, starcat[icat]->keymag2);
			}
		    else
			strcat (headline, "   Mag2");
		    }
		if (gcset)
		    strcat (headline, "     Peak");
		}
	    else
		strcat (headline, "  Mag");
	    if (ranges == NULL)
		strcat (headline, "  Arcsec");
	    if (wfile)
		fprintf (fd, "%s\n", headline);
	    else
		printf ("%s\n", headline);
	    }
	nohead = 0;
	}

    /* Find maximum separation for formatting */
    gdmax = 0.0;
    for (i = 0; i < ns; i++) {
	if (gx[i] > 0.0) {
	    gdist = 3600.0 * gx[i];
	    if (gdist > gdmax)
		gdmax = gdist;
	    }
	}

    string[0] = (char) 0;
    if (closest) ns = 1;
    for (i = 0; i < ns; i++) {
	if (typecol == 1) {
	    isp[0] = gc[i] / 1000;
	    isp[1] = gc[i] % 1000;
	    if (isp[0] == ' ' && isp[1] == ' ') {
		isp[0] = '_';
		isp[1] = '_';
		}
	    }
	if (refcat == GSC || refcat == GSCACT) {
	    ngsc = gc[i] / 10000;
	    gc[i] = gc[i] - (ngsc * 10000);
	    band = gc[i] / 100;
	    gc[i] = gc[i] - (band * 100);
	    }
	if (refcat == GSC2 || refcat == HIP) {
	    xmag = 0.01 * (double) (gc[i] / 10000);
	    xmag1 = 0.01 * (double) (gc[i] % 10000);
	    }
	else if (refcat == IRAS) {
	    xmag = 0.01 * (double) (gc[i] / 100000);
	    xmag1 = 0.01 * (double) (gc[i] % 100000);
	    }
	else if (refcat == TMPSC)
	    xmag = 0.001 * (double) gc[i];
	if (refcat == IRAS) {
	    if (gmb[i] > 100.0) {
		gmb[i] = gmb[i] - 100.0;
		lim1 = 'L';
		}
	    else
		lim1 = ' ';
	    if (gm[i] > 100.0) {
		gm[i] = gm[i] - 100.0;
		lim2 = 'L';
		}
	    else
		lim2 = ' ';
	    if (xmag > 100.0) {
		xmag = xmag - 100.0;
		lim3 = 'L';
		}
	    else
		lim3 = ' ';
	    if (xmag1 > 100.0) {
		xmag1 = xmag1 - 100.0;
		lim4 = 'L';
		}
	    else
		lim4 = ' ';
	    flux1 = 1000.0 * pow (10.0, -gmb[i] / 2.5);
	    flux2 = 1000.0 * pow (10.0, -gm[i] / 2.5);
	    flux3 = 1000.0 * pow (10.0, -xmag / 2.5);
	    flux4 = 1000.0 * pow (10.0, -xmag1 / 2.5);
	    }
	if (gy[i] > 0.0) {
	    if (degout) {
		deg2str (rastr, 32, gra[i], 5);
		deg2str (decstr, 32, gdec[i], 5);
		}
	    else {
		ra2str (rastr, 32, gra[i], 3);
		dec2str (decstr, 32, gdec[i], 2);
		}
	    if (gx[i] > 0.0)
		gdist = 3600.0 * gx[i];
	    else
		gdist = 0.0;

	    /* Convert proper motion to milliarcsec/year from deg/year */
	    if (mprop) {
		if (degout)
		    pra = (gpra[i] * 3600.0 * 100.0);
		else
		    pra = (gpra[i] * 240.0 * 100.0);
		pdec = gpdec[i] * 3600.0 * 100.0;
		}
	    if (starcat[icat] != NULL)
		CatNum (refcat,-nnfld,starcat[icat]->nndec,gnum[i],numstr);
	    else
		CatNum (refcat, -nnfld, nndec, gnum[i], numstr);

	    /* Print or write tab-delimited output line for one star */
	    if (tabout) {
	    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		refcat == UAC  || refcat == UA1  || refcat == UA2)
		sprintf (headline, "%s	%s	%s	%.1f	%.1f	%d",
		 numstr, rastr, decstr, gmb[i], gm[i], gc[i]);
	    else if (refcat == GSC || refcat == GSCACT)
	        sprintf (headline,
			 "%s	%s	%s	%.2f	%d	%d	%d",
			 numstr, rastr, decstr, gm[i], gc[i], band, ngsc);
	    else if (refcat == TMPSC)
		sprintf (headline, "%s	%s	%s	%.3f	%.3f	%.3f",
		 numstr, rastr, decstr, gm[i], gmb[i], xmag);
	    else if (refcat == GSC2)
		sprintf (headline, "%s	%s	%s	%.2f	%.2f	%.2f	%.2f",
		 numstr, rastr, decstr, gm[i], gmb[i], xmag, xmag1);
	    else if (refcat == IRAS)
		sprintf (headline, "%s	%s	%s	%.2f%c	%.2f%c	%.2f%c	%.2f%c",
		 numstr,rastr,decstr,flux1,lim1,flux2,lim2,flux3,lim3,flux4,lim4);
	    else if (refcat == HIP)
		sprintf (headline, "%s	%s	%s	%.2f	%.2f	%.2f	%.2f",
		 numstr, rastr, decstr, gmb[i], gm[i], xmag, xmag1);
	    else if (refcat == UJC)
	        sprintf (headline, "%s	%s	%s	%.2f	%d",
		 numstr, rastr, decstr, gm[i], gc[i]);
	    else if (refcat==SAO || refcat==PPM || refcat == BSC) {
	        sprintf (headline, "%s	%s	%s	%.2f	%2s",
		 numstr, rastr, decstr, gm[i], isp);
		}
	    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==ACT) {
	        sprintf (headline, "%s	%s	%s	%.2f	%.2f",
		 numstr, rastr, decstr, gmb[i], gm[i]);
		}
	    else if (refcat == TABCAT)
		sprintf (headline, "%s	%s	%s	%.2f",
			 numstr, rastr, decstr, gm[i]);
	    else if (refcat == BINCAT) {
		sprintf (headline, "%s	%s	%s	%.2f	%2s",
			 numstr, rastr, decstr, gm[i], isp);
		}
	    else
	        sprintf (headline, "%s	%s	%s	%.2f",
		 numstr, rastr, decstr, gm[i]);
	    if (mprop) {
	        sprintf (temp, "	%7.3f	%6.2f", pra, pdec);
	        strcat (headline, temp);
		}
	    if (refcat == TABCAT && gcset) {
	        sprintf (temp, "	%d", gc[i]);
	        strcat (headline, temp);
		}

	    if (ranges == NULL) {
	        sprintf (temp, "	%.2f", gdist);
	        strcat (headline, temp);
		}
	    if (refcat == TABCAT && keyword != NULL) {
		strcat (headline, "	");
		strcat (headline, gobj[i]);
		}
	    if ((refcat == BINCAT || refcat == TXTCAT) &&
		 gobj1 != NULL && gobj[i] != NULL) {
		strcat (headline, "	");
		strcat (headline, gobj[i]);
		}
	    if (printxy) {
		strcat (headline, "	");
		strcat (numstr, xstr);
		strcat (headline, "	");
		strcat (numstr, ystr);
		}
		if (wfile)
		    fprintf (fd, "%s\n", headline);
		else
		    printf ("%s\n", headline);
		}

	    /* Print or write space-delimited output line for one star */
	    else {
		if (printxy) {
		    strcpy (numstr, xstr);
		    strcat (numstr, " ");
		    strcat (numstr, ystr);
		    }
		if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		    refcat == UAC  || refcat == UA1  || refcat == UA2)
		    sprintf (headline,"%s %s %s %5.1f %5.1f  %d ",
			     numstr,rastr,decstr,gmb[i],gm[i],gc[i]);
		else if (refcat == TMPSC)
		    sprintf (headline, "%s %s %s %6.3f %6.3f %6.3f",
			     numstr, rastr, decstr, gm[i], gmb[i], xmag);
		else if (refcat == IRAS)
		    sprintf (headline, "%s %s %s %5.2f%c %5.2f%c %5.2f%c %5.2f%c",
		 numstr,rastr,decstr,flux1,lim1,flux2,lim2,flux3,lim3,flux4,lim4);
		else if (refcat == GSC2 || refcat == HIP)
		    sprintf (headline, "%s %s %s %5.2f %5.2f %5.2f %5.2f",
			     numstr, rastr, decstr, gm[i], gmb[i], xmag,xmag1);
		else if (refcat == GSC || refcat == GSCACT)
	            sprintf (headline, "%s %s %s %6.2f %4d %4d %2d",
			     numstr, rastr, decstr, gm[i], gc[i], band, ngsc);
		else if (refcat == UJC)
		    sprintf (headline,"%s %s %s %6.2f %4d",
			     numstr, rastr, decstr, gm[i],gc[i]);
		else if (refcat==SAO || refcat==PPM ||
			 refcat==IRAS || refcat == BSC) {
		    sprintf (headline,"  %s  %s %s %6.2f  %2s",
			     numstr,rastr,decstr,gm[i],isp);
		    }
		else if (refcat==TYCHO || refcat==TYCHO2 || refcat==ACT) {
		    sprintf (headline,"%s %s %s %6.2f %6.2f ",
			     numstr,rastr,decstr,gmb[i],gm[i]);
		    }
		else if (refcat == TABCAT) {
		    if (gcset)
			sprintf (headline,"%s %s %s %6.2f %7d",
				 numstr, rastr, decstr, gm[i],gc[i]);
		    else
			sprintf (headline,"%s %s %s %6.2f",
				 numstr, rastr, decstr, gm[i]);
		    }
		else if (refcat == BINCAT) {
		    sprintf (headline,"%s %s %s %6.2f %2s",
			     numstr, rastr, decstr, gm[i], isp);
		    }
		else {
		    sprintf (headline, "%s %s %s %6.2f",
			     numstr, rastr, decstr, gm[i]);
		    }
		if (ranges == NULL) {
		    if (gdmax < 100.0)
			sprintf (temp, "  %5.2f", gdist);
		    else if (gdmax < 1000.0)
			sprintf (temp, "  %6.2f", gdist);
		    else if (gdmax < 10000.0)
			sprintf (temp, "  %7.2f", gdist);
		    else if (gdmax < 100000.0)
			sprintf (temp, "  %8.2f", gdist);
		    else
			sprintf (temp, "  %.2f", gdist);
		    strcat (headline, temp);
		    }
		if (refcat == TABCAT && keyword != NULL) {
		    sprintf (temp, " %s", gobj[i]);
		    strcat (headline, temp);
		    }
		else if ((refcat == BINCAT || refcat == TXTCAT) &&
			 gobj1 != NULL && gobj[i] != NULL) {
		    sprintf (temp, " %s", gobj[i]);
		    strcat (headline, temp);
		    }
		if (wfile)
		    fprintf (fd, "%s\n", headline);
		else
		    printf ("%s\n", headline);
		}
	    }
	}

	/* If searching more than one catalog, separate them with blank line */
	if (ncat > 0 && icat < ncat-1)
	    printf ("\n");

	/* Free memory used for object names in current catalog */
	if (gobj1 != NULL) {
	    for (i = 0; i < ns; i++) {
		if (gobj[i] != NULL) free (gobj[i]);
		gobj[i] = NULL;
		}
	    }
	}

    /* Close output file */
    if (wfile)
	fclose (fd);

    return (ns);
}


/* Get a center and radius for a search area.  If the image center is not
 * given in the system of the reference catalog, convert it.
 * Return 0 if OK, else -1
 */

static int
GetArea (verbose, syscoor, sysout, eqout, epout, cra, cdec, dra, ddec, drad)

int	verbose;	/* Extra printing if =1 */
int	syscoor;	/* Coordinate system of input search coordinates */
int	sysout;		/* Coordinate system of output coordinates */
double	eqout;		/* Equinox in years of output coordinates */
double	epout;		/* Epoch in years of output coordinates (0=eqcoor */
double	*cra;		/* Center longitude/right ascension (degrees returned)*/
double	*cdec;		/* Center latitude/declination (degrees returned) */
double	*dra;		/* Longitude/RA half-width (degrees returned) */
double	*ddec;		/* Latitude/Declination half-width (degrees returned) */
double	*drad;		/* Radius to search in degrees (0=box) (returned) */
{
    char rstr[32], dstr[32], cstr[32];

    *cra = ra0;
    *cdec = dec0;
    if (verbose) {
	if (syscoor == WCS_XY) {
	    num2str (rstr, *cra, 10, 5);
            num2str (dstr, *cdec, 10, 5);
	    }
	else if (syscoor==WCS_ECLIPTIC || syscoor==WCS_GALACTIC || degout0) {
	    deg2str (rstr, 32, *cra, 5);
            deg2str (dstr, 32, *cdec, 5);
	    }
	else {
	    ra2str (rstr, 32, *cra, 3);
            dec2str (dstr, 32, *cdec, 2);
	    }
	wcscstr (cstr, syscoor, 0.0, 0.0);
	fprintf (stderr,"Center:  %s   %s %s\n", rstr, dstr, cstr);
	}

    if (syscoor != sysout) {
	wcscon (syscoor, sysout, 0.0, 0.0, cra, cdec, epout);
	if (verbose) {
	    if (syscoor == WCS_ECLIPTIC || syscoor == WCS_GALACTIC || degout0) {
		deg2str (rstr, 32, *cra, 5);
        	deg2str (dstr, 32, *cdec, 5);
		}
	    else {
		ra2str (rstr, 32, *cra, 3);
        	dec2str (dstr, 32, *cdec, 2);
		}
	    wcscstr (cstr, syscoor, 0.0, 0.0);
	    fprintf (stderr,"Center:  %s   %s %s\n", rstr, dstr, cstr);
	    }
	}

    /* Set search box radius from command line, if it is there */
    if (dra0 > 0.0) {
	*drad = 0.0;
	if (syscoor == WCS_XY) {
	    *dra = dra0;
	    *ddec = ddec0;
	    }
	else {
	    *dra = dra0 / 3600.0;
	    *ddec = ddec0 / 3600.0;
	    }
	}

    /* Search box */
    else if (rad0 < 0.0) {
	*drad = 0.0;
	if (syscoor == WCS_XY) {
	    *dra = -rad0;
	    *ddec = -rad0;
	    }
	else {
	    *ddec = -rad0 / 3600.0;
	    if (*cdec < 90.0 && *cdec > -90.0)
		*dra = *ddec / cos (degrad (*cdec));
	    else
		*dra = 180.0;
	    }
	}

    /* Search circle */
    else if (rad0 > 0.0) {
	if (syscoor == WCS_XY) {
	    *drad = rad0;
	    *dra = *drad;
	    *ddec = *drad;
	    }
	else {
	    *drad = rad0 / 3600.0;
	    *ddec = *drad;
	    if (*cdec < 90.0 && *cdec > -90.0)
		*dra = *ddec / cos (degrad (*cdec));
	    else
		*dra = 180.0;
	    }
	}
    else {
	if (verbose)
	    fprintf (stderr, "GetArea: Illegal radius, rad= %.5f\n",rad0);
	return (-1);
	}

    if (verbose) {
	if (syscoor == WCS_XY) {
	    num2str (rstr, *dra * 2.0, 10, 5); 
	    num2str (dstr, *ddec * 2.0, 10, 5); 
	    }
	else {
	    ra2str (rstr, 32, *dra * 2.0, 2); 
	    dec2str (dstr, 32, *ddec * 2.0, 2); 
	    }
	fprintf (stderr,"Area:    %s x %s\n", rstr, dstr);
	}

    return (0);
}


static void
SearchHead (icat,sys1,sys2,eq1,eq2,ep1,ep2,cra,cdec,dra,ddec,drad,nnfld)

int	icat;		/* Number of catalog in list */
int	sys1;		/* Input coordinate system */
int	sys2;		/* Output coordinate system */
double	eq1, eq2;	/* Input and output equinoxes */
double	ep1,ep2;	/* Input and output epochs */
double	cra, cdec;	/* Center coordinates of search box in degrees */
double	dra, ddec;	/* Width and height of search box in degrees */
double	drad;		/* Radius of search region in degrees */
int	nnfld;		/* Number of characters in ID field */
{
    double ra, dec;
    char rastr[32];
    char decstr[32];
    char cstr[16];
    char oform[16];
    int refcat;
    char title[80];
    int sysref;
    int nmag, mprop;
    double eqref, epref;

    ra = cra;
    dec = cdec;
    wcscon (sys1, sys2, eq1, eq2, &ra, &dec, ep2);
    if (sys1 == WCS_XY) {
	num2str (rastr, ra, 10, 5);
	num2str (decstr, dec, 10, 5);
	}
    else if (sys2 == WCS_ECLIPTIC || sys2 == WCS_GALACTIC || degout0) {
	deg2str (rastr, 32, ra, 5);
	deg2str (decstr, 32, dec, 5);
	}
    else {
	ra2str (rastr, 32, ra, 3);
	dec2str (decstr, 32, dec, 2);
	}

    /* Set type of catalog being searched */
    if (!(refcat = RefCat (refcatname[icat], title, &sysref, &eqref, &epref,
			   &mprop, &nmag))) {
	fprintf (stderr,"ListCat: Catalog '%s' is missing\n",refcatname[icat]);
	return;
	}

    /* Label search center */
    if (objname) {
	sprintf (oform, "%%%ds", nnfld);
	printf (oform, objname);
	}
    else {
	if (refcat == ACT)
	    printf ("USNO ACT ");
	else if (refcat == GSC)
	    printf ("HST GSC  ");
	else if (refcat == USAC)
	    printf ("USNO SA      ");
	else if (refcat == USA1)
	    printf ("USNO SA-1.0  ");
	else if (refcat == USA2)
	    printf ("USNO SA-2.0  ");
	else if (refcat == UAC)
	    printf ("USNO A       ");
	else if (refcat == UA1)
	    printf ("USNO A-1.0   ");
	else if (refcat == UA2)
	    printf ("USNO A-2.0   ");
	else if (refcat == UJC)
	    printf ("USNO UJ1.0  ");
	else if (refcat == SAO)
	    printf ("SAO      ");
	else if (refcat == PPM)
	    printf ("PPM      ");
	else if (refcat == IRAS)
	    printf ("IRAS     ");
	else if (refcat == TYCHO)
	    printf ("Tycho    ");
	else if (refcat == TYCHO2)
	    printf ("Tycho2   ");
	else if (refcat == HIP)
	    printf ("Hip      ");
	else if (refcat == ACT)
	    printf ("ACT      ");
	else
	    printf ("%9s ", refcatname[icat]);
	}
    wcscstr (cstr, sys2, eq2, ep2);
    printf (" %s %s %s", rastr, decstr, cstr);
    if (rad0 != 0.0)
	printf (" r= %.2f", rad0);
    else
	printf (" +- %.2f %.2f", dra0, ddec0);
    if (classd == 0)
	printf (" stars");
    else if (classd == 3)
	printf (" nonstars");
    if (epoch0 != 0.0)
	printf (" at epoch %7.2f\n", epoch0);
    else
	printf ("\n");
    return;
}


void
PrintNum (maxnum, num, ndec)

double	maxnum;	/* Maximum value of number to print */
double	num;	/* Number to print */
int	ndec;	/* Number of decimal places in output */
{
    char nform[8];

    if (ndec > 0)
	sprintf (nform, "%%%d.%df", LenNum (maxnum,ndec), ndec);
    else
	sprintf (nform, "%%%dd", LenNum (maxnum,ndec));
    printf (nform, num);
}


int
LenNum (maxnum, ndec)

double	maxnum;	/* Maximum value of number to print */
int	ndec;	/* Number of decimal places in output */
{
    if (ndec <= 0)
	ndec = -1;
    if (maxnum < 9.999)
	return (ndec + 3);
    else if (maxnum < 99.999)
	return (ndec + 4);
    else if (maxnum < 999.999)
	return (ndec + 5);
    else if (maxnum < 9999.999)
	return (ndec + 6);
    else
	return (ndec + 7);
}


/* Set parameter values from cgi query line as kw=val&kw=val... */

static void
scatcgi (qstring)

char *qstring;
{
    char *pstring, *parend;

    pstring = qstring;
    while ((parend = strchr (pstring, '&')) != NULL) {
	*parend = (char) 0;
	if (scatparm (pstring))
	    fprintf (stderr, "SCATCGI: %s is not a parameter.\n", pstring);
	pstring = parend + 1;
	}
    if (scatparm (pstring))
	fprintf (stderr, "SCATCGI: %s is not a parameter.\n", pstring);
    return;
}


/* Set parameter values from the command line as keyword=value
 * Return 1 if successful, else 0 */

static int
scatparm (parstring)

char *parstring;
{
    char *parname;
    char *parvalue;
    char *parequal;
    char *temp;
    char *refcatn;
    int lcat, lrange;

    tabout = 1;

    /* First, check for request for help */
    if (!strncasecmp (parstring,"help", 4)) {
	printf ("Content-type: text/plain\n\n");
	printf ("Results will be returned as a tab-separated Starbase table\n");
	printf ("Enter a sequence of keyword=value, separated by & in one line\n");
	printf ("\n");
	printf ("keyword  value description\n");
	printf ("catalog  gsc, ua2(-USNO-A2.0), gsc2(GSC II), etc\n");
	printf ("ra       right ascension in degrees or hh:mm:ss.sss\n");
	printf ("dec      declination in degrees or [+/-]dd:mm:ss.sss\n");
	printf ("sys      coordinate system (B1950, J2000, Ecliptic, Galactic\n");
	printf ("outsys   output coordinate system if not same as above\n");
	printf ("nstar    maximum number of stars to be returned\n");
	printf ("rad      search radius in arcseconds or dd:mm:ss.sss\n");
	printf ("         (negate for square box)\n");
	printf ("dra      search halfwidth in RA arcseconds or dd:mm:ss.sss\n");
	printf ("ddec     search halfheight in Dec arcseconds or dd:mm:ss.sss\n");
	printf ("sort     dist   distance from search center\n");
	printf ("         ra     right ascension\n");
	printf ("         dec    declination\n");
	printf ("         mag    magnitude (brightest first)\n");
	printf ("         none   no sort\n");
	printf ("epoch    output epoch in fractional years\n");
	printf ("equinox  output equinox in fractional years\n");
	printf ("format   col    output coordinates hh:mm:ss dd:mm:ss\n");
	printf ("         deg    output coordinates in degrees\n");
	printf ("         rad    output coordinates in radians\n");
	printf ("min      minimum magnitude\n");
	printf ("max      maximum magnitude\n");
	printf ("num      range of catalog numbers (n-nxn or n,n,n)\n");
	exit (0);
	}

    /* Separate parameter name and value */
    parname = parstring;
    if ((parequal = strchr (parname,'=')) == NULL)
	return (0);
    *parequal = (char) 0;
    parvalue = parequal + 1;

    /* Get closest source */
    if (!strcasecmp (parname, "closest")) {
	if (!strncasecmp (parvalue, "y", 1)) {
	    catsort = SORT_DIST;
	    nstars = 1;
	    closest++;
	    }
	}

    /* Set range of source numbers to print */
    else if (!strncasecmp (parname,"num",3)) {
	if (ranges) {
	    temp = ranges;
	    lrange = strlen(ranges) + strlen(parvalue) + 2;
	    ranges = (char *) malloc (lrange);
	    strcpy (ranges, temp);
	    strcat (ranges, ",");
	    strcat (ranges, parvalue);
	    free (temp);
	    }
	else {
	    lrange = strlen (parvalue) + 2;
	    ranges = (char *) malloc (lrange);
	    if (strchr (parvalue,'.'))
		match = 1;
	    strcpy (ranges, parvalue);
	    }
	}

    /* Radius in arcseconds */
    else if (!strncasecmp (parname,"rad",3)) {
	if (strchr (parvalue,':'))
	    rad0 = 3600.0 * str2dec (parvalue);
	else
	    rad0 = atof (parvalue);
	}

    /* Search center right ascension */
    else if (!strcasecmp (parname,"ra"))
	ra0 = str2ra (parvalue);

    /* Search center declination */
    else if (!strcasecmp (parname,"dec"))
	dec0 = str2dec (parvalue);

    /* Search center coordinate system */
    else if (!strncasecmp (parname,"sys",3)) {
	syscoor = wcscsys (parvalue);
	eqcoor = wcsceq (parvalue);
	}

    /* Output coordinate system */
    else if (!strcasecmp (parname, "outsys")) {

	/* B1950 (FK4) coordinates */
	if (!strcasecmp (parvalue, "B1950") ||
	    !strcasecmp (parvalue, "FK4")) {
	    sysout0 = WCS_B1950;
	    eqout = 1950.0;
    	    }

	/* J2000 (FK5) coordinates */
	else if (!strcasecmp (parvalue, "J2000") ||
	    !strcasecmp (parvalue, "FK5")) {
	    sysout0 = WCS_J2000;
	    eqout = 2000.0;
    	    }

	/* Galactic coordinates */
	else if (!strncasecmp (parvalue, "GAL", 3))
	    sysout0 = WCS_GALACTIC;

	/* Ecliptic coordinates */
	else if (!strncasecmp (parvalue, "ECL", 3))
	    sysout0 = WCS_ECLIPTIC;
	}

    /* Set reference catalog */
    else if (!strcasecmp (parname, "catalog")) {
	lcat = strlen (parvalue) + 2;
	refcatn = (char *) malloc (lcat);
	strcpy (refcatn, parvalue);
	refcatname[ncat] = refcatn;
	ncat = ncat + 1;
	}

    /* Set output coordinate epoch */
    else if (!strcasecmp (parname, "epoch"))
	epoch0 = fd2ep (parvalue);

    /* Output equinox in years */
    else if (!strcasecmp (parname, "equinox"))
    	eqout = fd2ep (parvalue);

    /* Output in degrees instead of sexagesimal */
    else if (!strcasecmp (parname, "format")) {
	if (!strncasecmp (parvalue, "deg", 3))
	    degout0 = 1;
	else if (!strncasecmp (parvalue, "rad", 3))
	    degout0 = 2;
	else
	    degout0 = 0;
	}

    /* Print center and closest star on one line */
    else if (!strcasecmp (parname, "oneline")) {
	if (parvalue[0] == 'y' || parvalue[0] == 'Y') {
	    oneline++;
	    catsort = SORT_DIST;
	    closest++;
	    nstars = 1;
	    }
	}

    /* Magnitude limit */
    else if (!strncasecmp (parname,"mag",3)) {
	maglim2 = atof (parvalue);
	if (MAGLIM1 == MAGLIM2)
	    maglim1 = -2.0;
	}
    else if (!strncasecmp (parname,"max",3))
	maglim2 = atof (parvalue);
    else if (!strncasecmp (parname,"min",3))
	maglim1 = atof (parvalue);

    /* Number of brightest stars to read */
    else if (!strncasecmp (parname,"nstar",5))
	nstars = atoi (parvalue);

    /* Object name */
    else if (!strcmp (parname, "object") || !strcmp (parname, "OBJECT")) {
	lcat = strlen (parvalue) + 2;
	objname = (char *) malloc (lcat);
	strcpy (objname, parvalue);
	}

    /* Output sorting */
    else if (!strcasecmp (parname, "sort")) {

	/* Sort by distance from center */
	if (!strncasecmp (parvalue,"di",2))
	    catsort = SORT_DIST;

	/* Sort by RA */
	if (!strncasecmp (parvalue,"r",1))
	    catsort = SORT_RA;

	/* Sort by Dec */
	if (!strncasecmp (parvalue,"de",2))
	    catsort = SORT_DEC;

	/* Sort by magnitude */
	if (!strncasecmp (parvalue,"m",1))
	    catsort = SORT_MAG;

	/* Sort by magnitude */
	if (!strncasecmp (parvalue,"n",1))
	    catsort = NOSORT;
	}

    /* Search box half-width in RA */
    else if (!strcasecmp (parname,"dra")) {
	if (strchr (parvalue,':'))
	    dra0 = 3600.0 * str2ra (parvalue);
	else
	    dra0 = atof (parvalue);
	if (ddec0 <= 0.0)
	    ddec0 = dra0;
	/* rad0 = sqrt (dra0*dra0 + ddec0*ddec0); */
	}

    /* Search box half-height in Dec */
    else if (!strcasecmp (parname,"ddec")) {
	if (strchr (parvalue,':'))
	    ddec0 = 3600.0 * str2dec (parvalue);
	else
	    ddec0 = atof (parvalue);
	if (dra0 <= 0.0)
	    dra0 = ddec0;
	rad0 = sqrt (dra0*dra0 + ddec0*ddec0);
	}

    /* Guide Star object class */
    else if (!strncasecmp (parname,"cla",3)) {
	classd = (int) atof (parvalue);
	setgsclass (classd);
	}

    /* Catalog to be searched */
    else if (!strncasecmp (parname,"cat",3)) {
	lcat = strlen (parvalue) + 2;
	refcatn = (char *) malloc (lcat);
	strcpy (refcatn, parvalue);
	refcatname[ncat] = refcatn;
	ncat = ncat + 1;
	}
    else {
	*parequal = '=';
	return (1);
	}
    *parequal = '=';
    return (0);
}

/* Oct 18 1996	New program based on imtab
 * Nov 13 1996	Set maximum nstar from command line if greater than default
 * Nov 14 1996	Set limits from subroutine
 * Nov 19 1996	Fix usage
 * Dec 12 1996	Allow bright as well as faint magnitude limit
 * Dec 12 1996	Fix header for UAC
 * Dec 12 1996	Fix header for UAC magnitudes
 * Dec 13 1996	Write plate into header if selected
 * Dec 18 1996	Allow WCS sky coordinate format as input argument for center
 * Dec 18 1996	Add option to print entries for specified catalog numbers
 * Dec 30 1996	Clean up closest star message
 * Dec 30 1996	Print message instead of heading if no stars are found
 *
 * Jan 10 1997	Fix bug in RASort Stars which did not sort magnitudes
 * Mar 12 1997	Add USNO SA-1.0 catalog as USAC
 * Apr 25 1997	Fix bug in uacread
 * May 29 1997	Add option to add keyword to tab table output
 * Nov 12 1997	Fix DEC in header to print Dec string instead of RA string
 * Nov 17 1997	Initialize both magnitude limits
 * Dec  8 1997	Set up program to be called by various names
 * Dec 12 1997	Fix center coordinate printing in heading
 *
 * Apr 10 1998	Fix bug search USNO A-1.0 catalog created by last revision
 * Apr 10 1998	Set search radius only if argument negative
 * Apr 14 1998	Version 2.2: to match other software
 * Apr 23 1998	Add output in galactic or ecliptic coordinates; use wcscon()
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * Jun  2 1998	Fix bug in tabread()
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Jun 30 1998	Fix declaration of GetArea()
 * Jul  1 1998	Allow specification of center different from output system
 * Jul  9 1998	Adjust all report headings
 * Jul 30 1998	Realign heading and fix help
 * Aug  6 1998	Do not include fitshead.h; it is in wcs.h
 * Sep 10 1998	Add SAOTDC binary format catalogs
 * Sep 15 1998	Adjust output format for binary format catalogs
 * Sep 16 1998	Fix bug creating output filename
 * Sep 21 1998	Do not print distance to search center if not searching
 * Sep 21 1998	Print epoch if not that of equinox
 * Sep 24 1998	Increase search radius for closest star
 * Sep 24 1998	Add second magnitude for Tycho Catalogue
 * Oct 15 1998	Add ability to read TDC ASCII catalog files
 * Oct 16 1998	Add ability to read any TDC binary catalog file
 * Oct 21 1998	Add object name to TDC binary catalogs
 * Oct 21 1998	Use wcscat.h common
 * Oct 23 1998	Allow up to 10 catalogs to be searched at once
 * Oct 26 1998	Return object name in same operation as object position
 * Oct 27 1998	Fix RefCat() calls
 * Oct 29 1998	Add GSC class selection argument x; it should have been there
 * Oct 30 1998	Read object name if accessing catalogs by number, too
 * Nov 20 1998	Implement USNO A-2.0 and SA-2.0 catalogs; differentiate from A1
 * Nov 30 1998	Add version and help commands for consistency
 * Dec  8 1998	Add Hipparcos and ACT catalogs
 * Dec 14 1998	Fix format for UJC
 * Dec 21 1998	Fix format for BINCAT and TXTCAT format catalogs

 * Jan 20 1999	Add option to print search center as output
 * Jan 21 1999	Drop option of adding coordinates to -b -e -g -j
 * Jan 21 1999	Improve command parser to accept fractional degree centers
 * Jan 21 1999	Set output coordinate system for TDC ASCII and Binary catalogs
 * Jan 26 1999	Add option of range of star numbers
 * Feb  1 1999	Add switch to get sequence number unless . in number
 * Feb  2 1999	Vary number of decimal places according to input ASCII catalog
 * Feb  2 1999	Set output equinox, epoch, coorsys from catalog if not set
 * Feb 10 1999	Increase search radius for closest star in smaller catalogs
 * Feb 12 1999	Add ACT catalog from CDROM
 * Feb 18 1999	Make printing name in search centers optional
 * Apr 13 1999	Fix progname to drop / when full pathname
 * May 12 1999	Add option to search from a TDC ASCII catalog
 * May 19 1999	Add option to print search center and closest star on 1 line
 * May 19 1999	Format catalog number using CatNum()
 * May 21 1999	Allow option of setting epoch from search catalog
 * May 28 1999	Allow search from starbase table catalog as well as ASCII
 * May 28 1999	Write epoch of coordinates into tab header if not = equinox
 * Jun  4 1999	Improve labelling for search from file
 * Jun  4 1999	Use calloc() instead of malloc() for proper initialization
 * Jun  4 1999	Allow rectangular in addition to square and circular searches
 * Jun  7 1999	Allow radius input in sexagesimal
 * Jun 30 1999	Use isrange() to check for a range of source numbers
 * Jul  1 1999	Allow any legal FITS date format for epoch
 * Aug 16 1999	Be more careful checking for 2nd -r argument
 * Aug 20 1999	Change u from USNO plate to print X Y 
 * Aug 23 1999	Fix closest star search by setting search radius not box
 * Aug 25 1999	Add the Bright Star Catalog
 * Sep  8 1999	Clean up code
 * Sep 10 1999	Do all searches through catread() and catrnum()
 * Sep 16 1999	Fix galactic coordinate header
 * Sep 16 1999	Add distsort argument to catread() call
 * Sep 21 1999	Change description of search area from square to region
 * Oct 15 1999	Fix calls to catopen() and catclose()
 * Oct 22 1999	Drop unused variables after lint
 * Oct 22 1999	Change catread() to ctgread() to avoid system conflict
 * Oct 22 1999	Increase default r for closest search to 1800 arcsec
 * Nov 19 1999	Use CatNum when concocting names from numbers
 * Nov 29 1999	Include fitsfile.h for date conversion
 *
 * Jan 11 2000	Get nndec for Starbase catalogs
 * Feb 10 2000	Use reference catalog coordinate system by default
 * Feb 15 2000	Use MAXCAT from lwcs.h instead of MAXREF
 * Feb 28 2000	Drop Peak column if not set in TAB catalogs
 * Mar 10 2000	Move catalog selection from executable name to subroutine
 * Mar 14 2000	Lowercase all header keywords
 * Mar 14 2000	Add proper motion, if present, to Starbase output
 * Mar 27 2000	Drop unused subroutine setradius() declaration
 * Mar 28 2000	Clean up output for catalog IDs and GSC classes
 * Apr  3 2000	Allocate search return buffers only once per execution
 * May 26 2000	Add Tycho 2 catalog
 * May 26 2000	Always use CatNumLen() to get ID number field size
 * May 31 2000	Do not sort if only one star
 * Jun  2 2000	Set to NULL when freeing
 * Jun 23 2000	Add degree output for one-line matches (-l)
 * Jun 26 2000	Add XY output
 * Jul 13 2000	Add star catalog data structure to ctgread() argument list
 * Jul 13 2000	Precess search catalog sources to specified system and epoch
 * Jul 25 2000	Pass address of star catalog data structure address
 * Sep  1 2000	Call CatNum with -nnfld to print leading zeroes on search ctrs
 * Sep 21 2000	Print spectral type instead of plate number of USNO A catalogs
 * Sep 25 2000	Print spectral type and second magnitude if present in one-line
 * Oct 26 2000	Print proper motion in msec/year and masec/year unless degout
 * Nov  3 2000	Print proper motion in sec/century and arcsec/century
 * Nov  9 2000	Set output equinox and epoch from -j and -b flags
 * Nov 17 2000	Add keyword=value command line parameter decoding
 * Nov 22 2000	Pass starcat to ctgrnum()
 * Nov 28 2000	Add CGI query decoding
 * Dec 15 2000	Add code to deal with coordinates without systems
 * Dec 18 2000	Change two-argument box size separator from space to comma
 * Dec 18 2000	Always allocate proper motion arrays
 * Dec 29 2000	Set debug flag if two v's encountered as command line arguments
 *
 * Jan  2 2001	Fix proper motion test; fix box size in heading
 * Mar  1 2001	If printing x and y, add .match extension to output file
 * Mar  1 2001	Print output file as tab/Starbase only if -t
 * Mar  1 2001	Add -z option to append to output file
 * Mar 23 2001	If catalog system, equinox, and epoch not set, set to J2000
 * Mar 23 2001	Set epoch and equinox to match search coords if not set
 * Apr 24 2001	Add HST GSC band output
 * May 22 2001	Rewrite sort options to include declination and no sort
 * May 23 2001	Add GSC-ACT catalog (updated HST GSC)
 * May 30 2001	Add 2MASS Point Source Catalog
 * Jun  5 2001	Read tab table search catalog one line at a time
 * Jun  6 2001	Set output equinox like output epoch; fix one-line tab bugs
 * Jun  7 2001	Add arguments to RefCat()
 * Jun 13 2001	Make -r argument radius if positive box half-width if negative
 * Jun 13 2001	Print id and magnitude column headings for tab table searches
 * Jun 13 2001	Print headings only of first of list of star numbers
 * Jun 14 2001	Initialize gobj elements to NULL
 * Jun 18 2001	Fix printing of long spectral type from search catalog
 * Jun 22 2001	Add GSC2 and modify Hipparcos output
 * Jun 25 2001	Add fifth digit to IRAS 3rd and 4th magnitudes
 * Jun 26 2001	If nstars < 0, print from catalog subroutine
 * Jul  3 2001	Add help to CGI response
 * Jul  5 2001	Add hundredths digit to magnitude limit
 */
