/* File scat.c
 * July 1, 1999
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
#include "wcs.h"
#include "libwcs/lwcs.h"
#include "libwcs/wcscat.h"

#define MAXREF 100

static void usage();

static int ListCat ();
extern void setcenter();
extern void setradius();
static void SearchHead();
static int GetArea();

static int verbose = 0;		/* Verbose/debugging flag */
static int wfile = 0;		/* True to print output file */
static int classd = -1;		/* Guide Star Catalog object classes */
static int uplate = 0;		/* UJ Catalog plate number to use */
static double maglim1 = MAGLIM1; /* Catalog bright magnitude limit */
static double maglim2 = MAGLIM2; /* Catalog faint magnitude limit */
static int sysout0 = 0;		/* Output coordinate system */
static double eqcoor = 2000.0;	/* Equinox of search center */
static double eqout = 0.0;	/* Equinox for output coordinates */
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
static int rasort = 0;		/* 1 to sort stars by brighness */
static int distsort = 0;	/* 1 to sort stars by distance from center */
static int debug = 0;		/* True for extra information */
static char *objname;		/* Object name for output */
static char *keyword;		/* Column to add to tab table output */
static char *progname;		/* Name of program as executed */
static char progpath[128];
static int ncat = 0;
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

main (ac, av)
int ac;
char **av;
{
    char *str;
    char rastr[16];
    char decstr[16];
    int i, lcat;
    char *refcatname[5];	/* reference catalog names */
    char *refcatn;
    char *temp;
    int systemp = 0;		/* Input search coordinate system */
    char *ranges;
    int istar;

    ranges = NULL;
    keyword = NULL;
    objname = NULL;

    /* Check name used to execute programe and set catalog name accordingly */
    strcpy (progpath, av[0]);
    progname = progpath;
    for (i = strlen (progpath); i > -1; i--) {
	if (progpath[i] > 63 && progpath[i] < 90)
	    progpath[i] = progpath[i] + 32;
	if (progpath[i] == '/') {
	    progname = progpath + i + 1;
	    break;
	    }
	}
    for (i = 0; i < strlen (progname); i++) {
	if (progname[i] > 95 && progname[i] < 123)
	    cpname[i] = progname[i] - 32;
	else
	    cpname[i] = progname[i];
	}
    if (strsrch (progname,"gsc") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "gsc");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"uac1") != NULL ||
	strsrch (progname,"ua1") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "ua1");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"uac2") != NULL ||
	strsrch (progname,"ua2") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "ua2");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"uac") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "uac");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"usac1") != NULL ||
	strsrch (progname,"usa1") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "usa1");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"usac2") != NULL ||
	strsrch (progname,"usa2") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "usa2");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"usac") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "usac");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"ujc") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "ujc");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"sao") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "sao");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"ppm") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "ppm");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"ira") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "iras");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"tyc") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "tycho");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"hip") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,16);
	strcpy (refcatn, "hipparcos");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"act") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "act");
	refcatname[0] = refcatn;
	}

    if (ac == 1)
        usage (progname);

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage (progname);
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage (progname);
	}

    /* crack arguments */
    for (av++; --ac > 0; av++) {

	/* Set RA, Dec, and equinox if WCS-generated argument */
	if (strsrch (*av,":") != NULL) {
	    if (ac < 3)
		usage(progname);
	    else {
		strcpy (rastr, *av);
		ac--;
		strcpy (decstr, *++av);
		ra0 = str2ra (rastr);
		dec0 = str2dec (decstr);
		ac--;
		syscoor = wcscsys (*++av);
		eqcoor = wcsceq (*av);
		}
	    }

	/* Set range and make a list of star numbers from it */
	else if (isrange (*av)) {
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
		if (strchr (*av,'.'))
		    match = 1;
		strcpy (ranges, *av);
		}
	    }

	/* Set decimal degree center or star number */
	else if (isnum (*av)) {

	    /* Set decimal degree center */
	    if (ac > 2 && (systemp = wcscsys (*(av + 2))) > 0) {
		ra0 = atof (*av);
		ac--;
		dec0 = atof (*++av);
		ac--;
		av++;
		syscoor = systemp;
		eqcoor = wcsceq (*av);
		}

	    /* Assume number to be star number if no coordinate system */
	    else {
		if (strchr (*av,'.'))
		    match = 1;
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
		verbose++;
		break;

	    case 'a':	/* Get closest source */
		distsort++;
		nstars = 1;
		break;

    	    case 'b':	/* output coordinates in B1950 */
		sysout0 = WCS_B1950;
		eqout = 1950.0;
    		break;

	    case 'c':       /* Set reference catalog */
		if (ac < 2)
		    usage(progname);
		lcat = strlen (*++av);
		refcatn = (char *) calloc (1, lcat + 1);
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
		    usage (progname);
		keyword = *++av;
		ac--;
		break;

	    case 'l':	/* Print center and closest star on one line */
		oneline++;
		distsort++;
		nstars = 1;
		break;

	    case 'm':	/* Magnitude limit */
		if (ac < 2)
		    usage (progname);
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
		    usage (progname);
		nstars = atoi (*++av);
		ac--;
		break;

	    case 'o':	/* Object name */
		if (ac < 2)
		    usage (progname);
		objname = *++av;
		ac--;
		break;

	    case 'p':	/* Sort by distance from center */
		distsort++;
		break;

    	    case 'q':	/* Output equinox in years */
    		if (ac < 2)
    		    usage(progname);
    		eqout = fd2ep (*++av);
    		ac--;
    		break;

    	    case 'r':	/* Box radius in arcseconds */
    		if (ac < 2)
    		    usage(progname);
		if (strchr (*++av,':'))
		    rad0 = 3600.0 * str2dec (*av);
		else
		    rad0 = atof (*av);
    		ac--;
		if (ac > 1 && isnum (*(av+1))) {
		    dra0 = rad0;
		    rad0 = 0.0;
		    if (strchr (*++av,':'))
			ddec0 = 3600.0 * str2dec (*av);
		    else
			ddec0 = atof (*av);
		    if (ddec0 <= 0.0)
			ddec0 = dra0;
		    rad0 = sqrt (dra0*dra0 + ddec0*ddec0);
		    ac--;
		    }
    		break;

	    case 's':	/* sort by RA */
		rasort = 1;
		break;

	    case 't':	/* tab table to stdout */
		tabout = 1;
		break;

	    case 'u':       /* USNO Catalog plate number */
		if (ac < 2)
		    usage(progname);
		uplate = (int) atof (*++av);
		ac--;
		break;

    	    case 'w':	/* write output file */
    		wfile++;
    		break;

	    case 'x':       /* Guide Star object class */
		if (ac < 2)
		    usage(progname);
		classd = (int) atof (*++av);
		ac--;
		break;

	    case 'y':	/* Set output coordinate epoch */
		if (ac < 2)
		    usage(progname);
		epoch0 = fd2ep (*++av);
		ac--;
		break;

	    default:
		usage (progname);
		break;
	    }
	    }
	else {
	    lcat = strlen (*av);
	    refcatn = (char *) calloc (1, lcat + 1);
	    strcpy (refcatn, *av);
	    refcatname[ncat] = refcatn;
	    ncat = ncat + 1;
	    }
	}

    if (readlist) {
	ranges = NULL;

	/* Read search center list from starbase tab table catalog */
	if (istab (listfile)) {
	    srchcat = tabcatopen (listfile);
	    if (srchcat != NULL) {
		srch = (struct Star *) calloc (1, sizeof (struct Star));
		for (istar = 1; istar <= srchcat->nstars; istar ++) {
		    if (tabstar (istar, srchcat, srch)) {
			fprintf (stderr,"%s: Cannot read star %d\n",
				 cpname, istar);
                	break;
                	}
		    ra0 = srch->ra;
		    dec0 = srch->dec;
		    syscoor = srch->coorsys;
		    eqcoor = srch->equinox;
		    if (epoch0 != 0.0) {
			ra0 = ra0 + (srch->rapm * (epoch0 - srch->epoch));
			dec0 = dec0 + (srch->decpm * (epoch0 - srch->epoch));
			}
		    ListCat (ncat, refcatname, ranges);
		    }
		tabcatclose (srchcat);
		}
	    }
	else {
	    srchcat = catopen (listfile);
	    if (srchcat != NULL) {
		srch = (struct Star *) calloc (1, sizeof (struct Star));
		for (istar = 1; istar <= srchcat->nstars; istar ++) {
		    if (catstar (istar, srchcat, srch)) {
			fprintf (stderr,"%s: Cannot read star %d\n",
				 cpname, istar);
                	break;
                	}
		    ra0 = srch->ra;
		    dec0 = srch->dec;
		    syscoor = srch->coorsys;
		    eqcoor = srch->equinox;
		    if (epoch0 != 0.0) {
			ra0 = ra0 + (srch->rapm * (epoch0 - srch->epoch));
			dec0 = dec0 + (srch->decpm * (epoch0 - srch->epoch));
			}
		    ListCat (ncat, refcatname, ranges);
		    }
		catclose (srchcat);
		}
	    }
	}
    else
	ListCat (ncat, refcatname, ranges);

    for (i = 0; i < ncat; i++)
	free (refcatname[i]);
    return (0);
}

static void
usage (progname)

char *progname;
{
    if (version)
	exit (-1);
    if (strsrch (progname,"gsc") != NULL)
	fprintf (stderr,"Find HST Guide Stars in a square on the sky\n");
    else if (strsrch (progname,"ujc") != NULL)
	fprintf (stderr,"Find USNO J Catalog stars in a square on the sky\n");
    else if (strsrch (progname,"uac") != NULL)
	fprintf (stderr,"Find USNO A Catalog stars in a square on the sky\n");
    else if (strsrch (progname,"ua1") != NULL)
	fprintf (stderr,"Find USNO A-1.0 Catalog stars in a square on the sky\n");
    else if (strsrch (progname,"ua2") != NULL)
	fprintf (stderr,"Find USNO A-2.0 Catalog stars in a square on the sky\n");
    else if (strsrch (progname,"usac") != NULL)
	fprintf (stderr,"Find USNO SA Catalog stars in a square on the sky\n");
    else if (strsrch (progname,"usa1") != NULL)
	fprintf (stderr,"Find USNO SA-1.0 Catalog stars in a square on the sky\n");
    else if (strsrch (progname,"usa2") != NULL)
	fprintf (stderr,"Find USNO SA-2.0 Catalog stars in a square on the sky\n");
    else if (strsrch (progname,"sao") != NULL)
	fprintf (stderr,"Find SAO Catalog stars in a square on the sky\n");
    else if (strsrch (progname,"ppm") != NULL)
	fprintf (stderr,"Find PPM Catalog stars in a square on the sky\n");
    else if (strsrch (progname,"ira") != NULL)
	fprintf (stderr,"Find IRAS Point Sources in a square on the sky\n");
    else if (strsrch (progname,"tyc") != NULL)
	fprintf (stderr,"Find Tycho Catalog stars in a square on the sky\n");
    else if (strsrch (progname,"hip") != NULL)
	fprintf (stderr,"Find Hipparcos Catalog stars in a square on the sky\n");
    else if (strsrch (progname,"act") != NULL)
	fprintf (stderr,"Find ACT Catalog stars in a square on the sky\n");
    else
	fprintf (stderr,"Find catalog stars in a square on the sky\n");
    fprintf (stderr,"Usage: [-abdehjlstvw] [-m [mag1] mag2] [-n num] [-r arcsec] ra dec sys or @list\n");
    fprintf(stderr,"  -a: List single closest catalog source\n");
    fprintf(stderr,"  -b: Output B1950 (FK4) coordinates\n");
    fprintf(stderr,"  -c: Reference catalog (act, gsc, ua2, usa2, or local file\n");
    fprintf(stderr,"  -d: Output RA and Dec in degrees instead of hms dms\n");
    fprintf(stderr,"  -e: Output ecliptic coordinates\n");
    fprintf(stderr,"  -f: Output search center for other programs\n");
    fprintf(stderr,"  -g: Output galactic coordinates\n");
    fprintf(stderr,"  -h: Print heading, else do not \n");
    fprintf(stderr,"  -i: Print catalog object name, not catalog number\n");
    fprintf(stderr,"  -j: Output J2000 (FK5) coordinates\n");
    fprintf(stderr,"  -k: Add this keyword to output from tab table search\n");
    fprintf(stderr,"  -l: Print center and closest star on one line\n");
    fprintf(stderr,"  -m: Magnitude limit(s)\n");
    fprintf(stderr,"  -n: Number of brightest stars to print \n");
    fprintf(stderr,"  -o: Object name \n");
    fprintf(stderr,"  -p: Sort by distance from center instead of flux\n");
    fprintf(stderr,"  -q: Equinox of output positions in years\n");
    fprintf(stderr,"  -r: Search half-width (<0=-radius) in arcsec (def 10)\n");
    fprintf(stderr,"  -s: Sort by RA instead of flux \n");
    fprintf(stderr,"  -t: Tab table to standard output as well as file\n");
    fprintf(stderr,"  -u: USNO catalog single plate number to accept\n");
    fprintf(stderr,"  -v: Verbose\n");
    fprintf(stderr,"  -w: Write tab table output file search[objname].[catalog]\n");
    fprintf(stderr,"  -x: GSC object type (0=stars 3=galaxies -1=all)\n");
    fprintf(stderr,"  -y: Epoch of output positions in years\n");
    exit (1);
}

#define TABMAX 64

static int
ListCat (ncat, refcatname, ranges)

int	ncat;		/* Number of catalogs to search */
char	**refcatname;	/* Reference catalog name */
char	*ranges;	/* String with range of catalog numbers to list */

{
    double *gnum;	/* Catalog star numbers */
    double *gra;	/* Catalog star right ascensions, rads */
    double *gdec;	/* Catalog star declinations rads */
    double *gm;		/* Catalog magnitudes */
    double *gmb;	/* Catalog B magnitudes */
    double *gx;		/* Catalog star X positions on image */
    double *gy;		/* Catalog star Y positions on image */
    int *gc;		/* Catalog star object classes */
    char **gobj;	/* Catalog star object names */
    char **gobj1;	/* Catalog star object names */
    double cra, cdec;
    double epout;
    int sysref;		/* Coordinate system of reference catalog */
    double eqref;	/* Equinox of catalog to be searched */
    double epref;	/* Epoch of catalog to be searched */
    int ng;		/* Number of catalog stars */
    int nbg;		/* Number of brightest catalog stars actually used */
    struct Range *range; /* Range of catalog numbers to list */
    int nfind;		/* Number of stars to find */
    int i, ngmax, nbytes, nbobj;
    double das, dds, drs;
    int degout;
    char nform[64];
    FILE *fd;
    char rastr[32], decstr[32];	/* coordinate strings */
    char numstr[32];	/* Catalog number */
    double drad, dra, ddec, mag1, mag2;
    double gdist, dr, da, dd, dec;
    int nlog, closest;
    int sysout;
    char headline[160];
    char filename[80];
    char title[80];
    char string[TABMAX];
    char temp[80];
    char isp[4];
    int refcat;		/* reference catalog switch */
    int icat, nndec;
    struct StarCat *starcat;
    double date, time;
    void ep2dt();
    void PrintNum();
    int LenNum();

    gnum = NULL;
    gra = NULL;
    gdec = NULL;
    gm = NULL;
    gmb = NULL;
    gx = NULL;
    gy = NULL;
    gc = NULL;
    gobj = NULL;
    gobj1 = NULL;

    /* Start of per catalog loop */
    for (icat = 0; icat < ncat; icat++) {
    nndec = 0;
    isp[2] = (char) 0;
    isp[3] = (char) 0;
    if (distsort && nstars == 1)
	closest = 1;
    else
	closest = 0;

    if (verbose || (printhead && notprinted)) {
	if (closest)
	else
	}

    if (sysout0) {
	sysout = sysout0;
	}
    else {
	sysout = syscoor;
	if (syscoor == WCS_B1950)
	    eqout = 1950.0;
	else if (ranges == NULL)
	    eqout = 2000.0;
	}
    if (debug)
	nlog = 1;
    else if (verbose)
	nlog = 100;
    else
	nlog = 0;

    if (!(refcat = RefCat (refcatname[icat], title, &sysref, &eqref, &epref))) {
	fprintf (stderr,"ListCat: Catalog '%s' is missing\n", refcatname[icat]);
	return (0);
	}

    /* Set epoch from command line, search catalog, or searched catalog */
    epout = epoch0;
    if (epout == 0.0) {
	if (srch!= NULL && srch->epoch != 0.0)
	    epout = srch->epoch;
	else
	    epout = epref;
	}


    /* Set degree flag for output */
    if (sysout == WCS_ECLIPTIC || sysout == WCS_GALACTIC)
	degout = 1;
    else
	degout = degout0;

    /* Find stars specified by number */
    if (ranges != NULL) {
	int nfdef = 9;
	range = RangeInit (ranges, nfdef);
	nfind = rgetn (range);
	nbytes = nfind * sizeof (double);
	nbobj = nfind * sizeof (char *);
	gnum = (double *) calloc (nfind, sizeof(double));
	if (!(gnum = (double *) calloc (nfind, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gnum\n", nbytes);
	else {
	    for (i = 0; i < nfind; i++)
		gnum[i] = rgetr8 (range);
	    }
	if (!(gra = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gra\n", nbytes);
	if (!(gdec = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gdec\n", nbytes);
	if (!(gm = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gm\n", nbytes);
	if (!(gmb = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gmb\n", nbytes);
	if (!(gc = (int *) calloc (ngmax, sizeof(int))))
	    fprintf (stderr, "Could not calloc %d bytes for gc\n", nbytes);
	if (!(gx = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gx\n", nbytes);
	if (!(gy = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gy\n", nbytes);
	if (!(gobj = (char **) calloc (nfind, sizeof(char *))))
	    fprintf (stderr, "Could not calloc %d bytes for obj\n", nbobj);
	if (!gnum || !gra || !gdec || !gm || !gmb || !gc || !gx || !gy || !gobj) {
	    if (gm) free ((char *)gm);
	    if (gmb) free ((char *)gmb);
	    if (gra) free ((char *)gra);
	    if (gdec) free ((char *)gdec);
	    if (gnum) free ((char *)gnum);
	    if (gc) free ((char *)gc);
	    if (gx) free ((char *)gx);
	    if (gy) free ((char *)gy);
	    if (gobj) free ((char *)gobj);
	    return (0);
	    }
	wfile = 0;

	/* Set output coordinate system if not set on the command line*/
	if (refcat == BINCAT) {
	    starcat = binopen (refcatname[icat]);
	    if (!sysout)
		sysout = starcat->coorsys;
	    if (!eqout)
		eqout = starcat->equinox;
	    if (!epout)
		epout = starcat->epoch;
	    nndec = starcat->nndec;
	    binclose (starcat);
	    }
	else if (refcat == TXTCAT) {
	    starcat = catopen (refcatname[icat]);
	    if (!sysout)
		sysout = starcat->coorsys;
	    if (!eqout)
		eqout = starcat->equinox;
	    if (!epout)
		epout = starcat->epoch;
	    nndec = starcat->nndec;
	    catclose (starcat);
	    }
	else {
	    if (!sysout)
		sysout = WCS_J2000;
	    if (!eqout)
		eqout = 2000.0;
	    if (!epout)
		epout = eqout;
	    }

	/* Find the specified catalog stars */
	if (refcat == GSC)
	    nbg = gscrnum (nfind,sysout,eqout,epout,
			   gnum,gra,gdec,gm,gc,nlog);
	else if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		 refcat == UAC  || refcat == UA1  || refcat == UA2)
	    nbg = uacrnum (refcatname[icat],nfind,sysout,eqout,epout,
			   gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == UJC)
	    nbg = ujcrnum (nfind,sysout,eqout,epout,
			   gnum,gra,gdec,gm,gc,nlog);
	else if (refcat == SAO)
	    nbg = binrnum ("SAO",nfind,sysout,eqout,epout,match,
			   gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == PPM)
	    nbg = binrnum ("PPM",nfind,sysout,eqout,epout,match,
			   gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == IRAS)
	    nbg = binrnum ("IRAS",nfind,sysout,eqout,epout,match,
			   gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == TYCHO)
	    nbg = binrnum ("tycho",nfind,sysout,eqout,epout,match,
			   gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == HIP)
	    nbg = binrnum ("hipparcos",nfind,sysout,eqout,epout,match,
			   gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == ACT)
	    nbg = actrnum (nfind,sysout,eqout,epout,
			   gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == TABCAT) {
	    nbg = tabrnum (refcatname[icat],nfind,sysout,eqout,epout,
			   gnum,gra,gdec,gm,gc,debug);
	    if (keyword != NULL)
		tabrkey (refcatname[icat], nfind, gnum, keyword, gobj);
	    }
	else if (refcat == BINCAT) {
	    nbg = binrnum (refcatname[icat],nfind,sysout,eqout,epout,match,
			   gnum,gra,gdec,gm,gmb,gc,gobj,nlog);
	    starcat = binopen (refcatname[icat]);
	    if (!sysout)
		sysout = starcat->coorsys;
	    if (!eqout)
		eqout = starcat->equinox;
	    if (!epout)
		epout = starcat->epoch;
	    nndec = starcat->nndec;
	    binclose (starcat);
	    }
	else {
	    nbg = catrnum (refcatname[icat],nfind,sysout,eqout,epout,match,
			   gnum,gra,gdec,gm,gobj,debug);
	    starcat = catopen (refcatname[icat]);
	    if (!sysout)
		sysout = starcat->coorsys;
	    if (!eqout)
		eqout = starcat->equinox;
	    if (!epout)
		epout = starcat->epoch;
	    nndec = starcat->nndec;
	    catclose (starcat);
	    }

	for (i = 0; i < nbg; i++ ) {
	    gx[i] = 0.0;
	    gy[i] = 1.0;
	    }

	/* Write out entries for use as image centers */
	if (searchcenter) {
	    for (i = 0; i < nbg; i++ ) {
		if (printobj) {
		    if (refcat == TABCAT && keyword != NULL)
			printf ("%s ", gobj[i]);
		    else if ((refcat == BINCAT || refcat == TXTCAT) &&
			     gobj != NULL && gobj[i] != NULL)
			printf ("%s ", gobj[i]);
		    else if (nndec > 0) {
			sprintf (nform, "%%s_%%%d.%d ", nndec+5, nndec);
			printf (nform, refcatname[icat], gnum[i]);
			}
		    else if (gnum[i] < 100.0)
			printf ("%s_%03d ",refcatname[icat], (int)gnum[i]);
		    else
			printf ("%s_%d ",refcatname[icat], (int)gnum[i]);
		    }
		else if (nndec > 0) {
		    sprintf (nform, "%%s_%%%d.%d ", nndec+5, nndec);
		    printf (nform, refcatname[icat], gnum[i]);
		    }
		else if (gnum[i] < 100.0)
		    printf ("%s_%03d ",refcatname[icat], (int)gnum[i]);
		else
		    printf ("%s_%d ",refcatname[icat], (int)gnum[i]);
		if (degout) {
		    deg2str (rastr, 32, gra[i], 5);
		    deg2str (decstr, 32, gdec[i], 5);
		    }
		else {
		    ra2str (rastr, 32, gra[i], 3);
		    dec2str (decstr, 32, gdec[i], 2);
		    }
		if (sysout == WCS_J2000)
		    printf ("%s %s J2000\n", rastr, decstr);
		else if (sysout == WCS_B1950)
		    printf ("%s %s B1950\n", rastr, decstr);
		else if (sysout == WCS_GALACTIC)
		    printf ("%s %s galactic\n", rastr, decstr);
		else if (sysout == WCS_ECLIPTIC)
		    printf ("%s %s ecliptic\n", rastr, decstr);
		else
		    printf ("%s %s none\n", rastr, decstr);
		}
	    return (nbg);
	    }
	}

    /* Find stars specified by location */
    else {

	/* Set search radius if finding closest star */
	if (rad0 == 0.0 && dra0 == 0.0) {
	    if (closest) {
		if (refcat == GSC || refcat == UJC ||
		    refcat == USAC || refcat == USA1 || refcat == USA2 ||
		    refcat == UAC  || refcat == UA1  || refcat == UA2)
		    rad0 = 1800.0;
		else
		    rad0 = 360.0;
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
		strlen (srch->objname) > 0)
		strcat (objname, srch->objname);
	    else {
		if (objname == NULL)
		    objname = (char *) calloc (1,32);
		CatNum (TXTCAT, srchcat->nndec, srch->num, objname);
		}
	    }

	/* Print search center and size in input and output coordinates */
	if (verbose || (printhead && !oneline)) {
	    SearchHead (refcatname[icat],sysout,sysout,eqout,eqout,epout,epout,
			cra,cdec,dra,ddec,drad,nndec);
	    if (sysout != syscoor && !closest)
		SearchHead (refcatname[icat],sysout,syscoor,eqout,eqcoor,
			    epout,epout,cra,cdec,dra,ddec,drad,nndec);
	    if (sysref != syscoor && sysref != sysout && !closest)
		SearchHead (refcatname[icat],sysout,sysref,eqout,eqref,
			    epout,epref,cra,cdec,dra,ddec,drad,nndec);
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

	/* Allocate memory for results of catalog search over image region */
	if (nstars > MAXREF)
	    ngmax = nstars;
	else
	    ngmax = MAXREF;
	nbytes = ngmax * sizeof (double);
	nbobj = ngmax * sizeof (char *);
	if (icat == 0) {
	if (!(gnum = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gnum\n", nbytes);
	if (!(gra = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gra\n", nbytes);
	if (!(gdec = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gdec\n", nbytes);
	if (!(gm = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gm\n", nbytes);
	if (!(gmb = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gmb\n", nbytes);
	if (!(gc = (int *) calloc (ngmax, sizeof(int))))
	    fprintf (stderr, "Could not calloc %d bytes for gc\n", nbytes);
	if (!(gx = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gx\n", nbytes);
	if (!(gy = (double *) calloc (ngmax, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gy\n", nbytes);
	if (!(gobj = (char **) calloc (ngmax, sizeof(char *))))
	    fprintf (stderr, "Could not calloc %d bytes for obj\n", nbobj);
	if (!gnum || !gra || !gdec || !gm || !gmb || !gc || !gx || !gy || !gobj) {
	    if (gm) free ((char *)gm);
	    if (gmb) free ((char *)gmb);
	    if (gra) free ((char *)gra);
	    if (gdec) free ((char *)gdec);
	    if (gnum) free ((char *)gnum);
	    if (gc) free ((char *)gc);
	    if (gx) free ((char *)gx);
	    if (gy) free ((char *)gy);
	    if (gobj) free ((char *)gobj);
	    return (0);
	    }
	    }

	/* Find the nearby reference stars, in ra/dec */
	if (refcat == GSC)
	    ng = gscread (cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,mag2,
			  classd,ngmax,gnum,gra,gdec,gm,gc,nlog);
	else if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		 refcat == UAC  || refcat == UA1  || refcat == UA2)
	    ng = uacread (refcatname[icat],cra,cdec,dra,ddec,drad,sysout,eqout,epout,
			  mag1,mag2,uplate,ngmax,gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == UJC)
	    ng = ujcread (cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,mag2,
			  uplate,ngmax,gnum,gra,gdec,gm,gc,nlog);
	else if (refcat == SAO)
	    ng = binread ("SAOra",cra,cdec,dra,ddec,drad,sysout,eqout,epout,
			  mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,NULL,nlog);
	else if (refcat == PPM)
	    ng = binread ("PPMra",cra,cdec,dra,ddec,drad,sysout,eqout,epout,
			  mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,NULL,nlog);
	else if (refcat == IRAS)
	    ng = binread ("IRAS",cra,cdec,dra,ddec,drad,sysout,eqout,epout,
			  mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,NULL,nlog);
	else if (refcat == TYCHO)
	    ng = binread ("tychora",cra,cdec,dra,ddec,drad,sysout,eqout,epout,
			  mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,NULL,nlog);
	else if (refcat == HIP)
	    ng = binread ("hipparcosra",cra,cdec,dra,ddec,drad,sysout,eqout,
		      epout,mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,NULL,nlog);
	else if (refcat == ACT)
	    ng = actread (cra,cdec,dra,ddec,drad,sysout,eqout,epout,
			  mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == BINCAT) {
	    ng = binread (refcatname[icat],cra,cdec,dra,ddec,drad,sysout,eqout,epout,
			  mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,gobj,nlog);
	    starcat = binopen (refcatname[icat]);
	    if (!sysout)
		sysout = starcat->coorsys;
	    if (!eqout)
		eqout = starcat->equinox;
	    if (!epout)
		epout = starcat->epoch;
	    nndec = starcat->nndec;
	    binclose (starcat);
	    }
	else if (refcat == TABCAT) {
	    ng = tabread (refcatname[icat],cra,cdec,dra,ddec,drad,sysout,eqout,epout,
			  mag1,mag2,ngmax,gnum,gra,gdec,gm,gc,nlog);
	    if (keyword != NULL)
		tabrkey (refcatname[icat], ng, gnum, keyword, gobj);
	    }
	else {
	    ng = catread (refcatname[icat],cra,cdec,dra,ddec,drad,sysout,eqout,epout,
			  mag1,mag2,ngmax,gnum,gra,gdec,gm,gobj,nlog);
	    starcat = catopen (refcatname[icat]);
	    if (!sysout)
		sysout = starcat->coorsys;
	    if (!eqout)
		eqout = starcat->equinox;
	    if (!epout)
		epout = starcat->epoch;
	    nndec = starcat->nndec;
	    catclose (starcat);
	    }
	if (gobj[0] == NULL)
	    gobj1 = NULL;
	else
	    gobj1 = gobj;

	/* Compute distance from star to search center */
	for (i = 0; i < ng; i++ ) {
	    gx[i] = wcsdist (cra, cdec, gra[i], gdec[i]);
	    gy[i] = 1.0;
	    }

	/* Sort reference stars from closest to furthest */
	if (distsort)
	    XSortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, gobj1, ng);

	/* Sort star-like objects in image by right ascension */
	else if (rasort)
	    RASortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, gobj1, ng);

	/* Sort reference stars from brightest to faintest */
	else
	    MagSortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, gobj1, ng);

	/* Print one line with search center and found star */
	if (oneline) {
	    if (ng > 0) {
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
		    int lid;

		    /* Write tab table heading */
		    if (refcat == GSC)
			printf ("CATALOG	HST GSC 1.1\n");
		    else if (refcat == USAC)
			printf ("CATALOG	USNO SA-1.0\n");
		    else if (refcat == USA1)
			printf ("CATALOG	USNO SA-1.0\n");
		    else if (refcat == USA2)
			printf ("CATALOG	USNO SA-2.0\n");
		    else if (refcat == UAC)
			printf ("CATALOG	USNO A-1.0\n");
		    else if (refcat == UA1)
			printf ("CATALOG	USNO A-1.0\n");
		    else if (refcat == UA2)
			printf ("CATALOG	USNO A-2.0\n");
		    else if (refcat == UJC)
			printf ("CATALOG	USNO UJ1.0\n");
		    else if (refcat == SAO)
			printf ("CATALOG	SAO\n");
		    else if (refcat == PPM)
			printf ("CATALOG	PPM\n");
		    else if (refcat == IRAS)
			printf ("CATALOG	IRAS Point Source\n");
		    else if (refcat == TYCHO)
			printf ("CATALOG	Tycho\n");
		    else if (refcat == HIP)
			printf ("CATALOG	Hipparcos\n");
		    else if (refcat == ACT)
			printf ("CATALOG	ACT\n");
		    else
			printf ("CATALOG	%s\n", refcatname[icat]);
		    printf ("SEARCH	%s\n", listfile);
		    printf ("EQUINOX	%.1f\n", eqout);
		    if (epout != eqout)
			printf ("EPOCH	%.1f\n", epout);
		    if (dra0 > 0.0) {
			printf ("DRASEC	%.2f\n", das);
			printf ("DDECSEC	%.2f\n", dds);
			}
		    else if (rad0 < 0)
			printf ("RADSEC	%.2f\n", drs);
		    else if (rad0 > 0)
			printf ("BOXSEC	%.2f\n", dds);

		    /* Write column headings */
		    if (srch != NULL)
			printf ("srch_id	");
		    printf ("srch_ra     	srch_dec    	");
		    if (srchcat->nmag != 0)
			printf ("srch_mag	");
		    if (srchcat->nepoch)
			printf ("epoch    	");
		    printf ("id");
		    lid = CatNumLen(refcat, nndec);
		    for (i = 0; i < lid-2; i++)
			printf (" ");
		    printf ("	ra          	dec        	");
		    printf ("mag    	n 	");
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
		    if (srch != NULL)
			printf ("-------	");
		    printf ("------------	------------	");
		    if (srchcat->nmag != 0)
			printf ("-----	");
		    if (srchcat->nepoch)
			printf ("---------	");
		    for (i = 0; i < lid; i++)
			printf ("-");
		    printf ("	------------	------------	");
		    printf ("-------	--	");
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
			CatNum (TXTCAT, srchcat->nndec, srch->num, numstr);
		    if (tabout)
			printf ("%s	", numstr);
		    else
			printf ("%s ", numstr);
		    }
		ra2str (rastr, 32, cra, 3);
		dec2str (decstr, 32, cdec, 2);
		if (tabout)
		    printf ("%s	%s", rastr, decstr);
		else
		    printf ("%s %s", rastr, decstr);
		if (srch != NULL) {
		    if (srchcat->nmag != 0) {
			if (tabout)
			    printf ("	%5.2f", srch->xmag[0]);
			else
			    printf (" %5.2f", srch->xmag[0]);
			}
		    if (srchcat->nepoch) {
			ep2dt (srch->epoch, &date, &time);
			if (tabout)
			    printf ("	%9.4f", date);
			else
			    printf (" %9.4f", date);
			}
		    }
		if (gobj1 != NULL) {
		    if (strlen (gobj1[0]) > 0)
			strcpy (numstr, gobj1[0]);
		    }
		else
		    CatNum (refcat, nndec, gnum[0], numstr);
		ra2str (rastr, 32, gra[0], 3);
		dec2str (decstr, 32, gdec[0], 2);
		if (tabout)
		    printf ("	%s	%s	%s	%5.2f	%d",
			    numstr, rastr, decstr, gm[0], ng);
		else
		    printf (" %s %s %s %5.2f %d",
			    numstr, rastr, decstr, gm[0], ng);
		dec = (gdec[0] + cdec) * 0.5;
		da = 3600.0 * (gra[0] - cra) * cos (degrad (dec));
		dd = 3600.0 * (gdec[0] - cdec);
		gdist = 3600.0 * gx[0];
		if (tabout)
		    printf ("	");
		else
		    printf (" ");
		PrintNum (das, da, 2);
		if (tabout)
		    printf ("	");
		else
		    printf (" ");
		PrintNum (dds, dd, 2);
		if (tabout)
		    printf ("	");
		else
		    printf (" ");
		PrintNum (drs, gdist, 2);
		printf ("\n");
		}
	    notprinted = 0;

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
	    return;
	    }

	/* List the brightest or closest MAXSTARS reference stars */
	if (nstars > 0 && ng > nstars) {
	    nbg = nstars;
	    if ((verbose || printhead) && !closest) {
		if (distsort) {
		    if (ng > 1)
			printf ("Closest %d / %d %s (closer than %.2f arcsec)",
				nbg, ng, title, 3600.0*gx[nbg-1]);
		    else
			printf ("Closest of %d %s",ng, title);
		    }
		else if (maglim1 > 0.0)
		    printf ("%d / %d %s (between %.1f and %.1f)",
			    nbg, ng, title, gm[0], gm[nbg-1]);
		else
		    printf ("%d / %d %s (brighter than %.1f)",
		 	    nbg, ng, title, gm[nbg-1]);
		printf ("\n");
		}
	    }
	else {
	    nbg = ng;
	    if (verbose || printhead) {
		if (maglim1 > 0.0)
		    printf ("%d %s between %.1f and %.1f\n",
			    ng, title, maglim1, maglim2);
		else if (maglim2 > 0.0)
		    printf ("%d %s brighter than %.1f\n",
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
		strcat (filename,refcatname[i]);
		}

	    fd = fopen (filename, "w");

	    /* Free result arrays and return if cannot write file */
	    if (fd == NULL) {
		fprintf (stderr, "%s:  cannot write file %s\n", cpname, filename);
		if (gx) free ((char *)gx);
		if (gy) free ((char *)gy);
		if (gm) free ((char *)gm);
		if (gmb) free ((char *)gmb);
		if (gra) free ((char *)gra);
		if (gdec) free ((char *)gdec);
		if (gnum) free ((char *)gnum);
		if (gc) free ((char *)gc);
        	return (0);
		}
	    }
        }

    /* Write heading */
    if (refcat == GSC)
	sprintf (headline, "CATALOG	HSTGSC1.1");
    else if (refcat == USAC)
	sprintf (headline, "CATALOG	USNO SA-1.0");
    else if (refcat == USA1)
	sprintf (headline, "CATALOG	USNO SA-1.0");
    else if (refcat == USA2)
	sprintf (headline, "CATALOG	USNO SA-2.0");
    else if (refcat == UAC)
	sprintf (headline, "CATALOG	USNO A-1.0");
    else if (refcat == UA1)
	sprintf (headline, "CATALOG	USNO A-1.0");
    else if (refcat == UA2)
	sprintf (headline, "CATALOG	USNO A-2.0");
    else if (refcat == UJC)
	sprintf (headline, "CATALOG	USNO UJ1.0");
    else if (refcat == SAO)
	sprintf (headline, "CATALOG	SAO");
    else if (refcat == PPM)
	sprintf (headline, "CATALOG	PPM");
    else if (refcat == IRAS)
	sprintf (headline, "CATALOG	IRAS Point Source");
    else if (refcat == TYCHO)
	sprintf (headline, "CATALOG	Tycho");
    else if (refcat == HIP)
	sprintf (headline, "CATALOG	Hipparcos");
    else if (refcat == ACT)
	sprintf (headline, "CATALOG	ACT");
    else
	sprintf (headline, "CATALOG	%s", refcatname[icat]);
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (!ranges) {
	ra2str (rastr, 32, cra, 3);
	if (wfile)
	    fprintf (fd, "RA	%s\n", rastr);
	if (tabout)
	    printf ("RA	%s\n", rastr);

	dec2str (decstr, 32, cdec, 2);
	if (wfile)
	    fprintf (fd, "DEC	%s\n", decstr);
	if (tabout)
	    printf ("DEC	%s\n", decstr);

	}
    if (wfile)
	fprintf (fd, "EQUINOX	%.1f\n", eqout);
    if (tabout)
	printf ("EQUINOX	%.1f\n", eqout);

    if (epout != eqout) {
	if (wfile)
	    fprintf (fd, "EPOCH	%.1f\n", epout);
	if (tabout)
	    printf ("EPOCH	%.1f\n", epout);
	}

    if (tabout) {
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
	    printf ("DRASEC	%.2f\n", das);
	    printf ("DDECSEC	%.2f\n", dds);
	    }
	else if (rad0 < 0)
	    printf ("RADSEC	%.2f\n", drs);
	else if (rad0 > 0)
	    printf ("BOXSEC	%.2f\n",dds);
	}

    if (uplate > 0) {
	if (wfile)
            fprintf (fd, "PLATE     %d\n", uplate);
        if (tabout)
            printf ("PLATE      %d\n", uplate);
        }

    if (distsort) {
	if (wfile)
	    fprintf (fd, "DISTSORT	T\n");
	if (tabout)
	    printf ("DISTSORT	T\n");
	}
    else if (rasort) {
	if (wfile)
	    fprintf (fd, "RASORT	T\n");
	if (tabout)
	    printf ("RASORT	T\n");
	}
    if (wfile)
    if (tabout)

    /* Print column headings */
    if (refcat == ACT)
	strcpy (headline, "act_id  	");
    else if (refcat == GSC)
	strcpy (headline, "gsc_id  	");
    else if (refcat == USAC)
	strcpy (headline,"usac_id  	");
    else if (refcat == USA1)
	strcpy (headline,"usa1_id  	");
    else if (refcat == USA2)
	strcpy (headline,"usa2_id  	");
    else if (refcat == UAC)
	strcpy (headline,"usnoa_id  	");
    else if (refcat == UA1)
	strcpy (headline,"usnoa1_id  	");
    else if (refcat == UA2)
	strcpy (headline,"usnoa2_id  	");
    else if (refcat == UJC)
	strcpy (headline,"usnoj_id  	");
    else if (refcat == SAO)
	strcpy (headline,"sao_id  	");
    else if (refcat == PPM)
	strcpy (headline,"ppm_id  	");
    else if (refcat == IRAS)
	strcpy (headline,"iras_id  	");
    else if (refcat == TYCHO)
	strcpy (headline,"tycho_id  	");
    else if (refcat == HIP)
	strcpy (headline,"hip_id  	");
    else if (refcat == ACT)
	strcpy (headline,"act_id  	");
    else
	strcpy (headline,"id       	");
    if (sysout == WCS_GALACTIC)
	strcat (headline,"long.ecl   	lat.ecl  	");
    else if (sysout == WCS_ECLIPTIC)
	strcat (headline,"long.ecl   	lat.ecl  	");
    else if (sysout == WCS_B1950)
	strcat (headline,"ra1950      	dec1950  	");
    else
	strcat (headline,"ra      	dec      	");
    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
	refcat == UAC  || refcat == UA1  || refcat == UA2)
	strcat (headline,"magb	magr	plate");
    else if (refcat == TYCHO || refcat == HIP || refcat == ACT)
	strcat (headline,"magb	magv	type");
    else if (refcat == TABCAT)
	strcat (headline,"mag	type");
    else
	strcat (headline,"mag");
    if (ranges == NULL)
	strcat (headline,"	arcsec");
    if (refcat == TABCAT && keyword != NULL) {
	strcat (headline,"	");
	strcat (headline, keyword);
	}
    if (gobj1 != NULL)
	strcat (headline,"	object");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (refcat == USAC || refcat == USA1 || refcat == USA2 || refcat == HIP ||
	refcat == UAC  || refcat == UA1  || refcat == UA2 || refcat == TYCHO ||
	refcat == ACT)
	sprintf(headline,"----------	--------	---------	----	-----	-----");
    else if (refcat == BINCAT)
        sprintf (headline,"----------	------------	------------	------	----");
    else
        sprintf (headline,"----------	------------	------------	------");
    if (ranges == NULL)
	strcat (headline, "	------");
    if (refcat == TABCAT && keyword != NULL)
	strcat (headline,"	------");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (printhead) {
	if (nbg == 0) {
	    if (refcat == GSC)
		printf ("No Guide Stars Found\n");
	    else if (refcat == USAC)
		printf ("No USNO SA Stars Found\n");
	    else if (refcat == USA1)
		printf ("No USNO SA-1.0 Stars Found\n");
	    else if (refcat == USA2)
		printf ("No USNO SA-2.0 Stars Found\n");
	    else if (refcat == UAC)
		printf ("No USNO A Stars Found\n");
	    else if (refcat == UA1)
		printf ("No USNO A-1.0 Stars Found\n");
	    else if (refcat == UA2)
		printf ("No USNO A-2.0 Stars Found\n");
	    else if (refcat == UJC)
		printf ("No UJ 1.0 Stars Found\n");
	    else if (refcat == SAO)
		printf ("No SAO Stars Found\n");
	    else if (refcat == PPM)
		printf ("No PPM Stars Found\n");
	    else if (refcat == IRAS)
		printf ("No IRAS Point Sources Found\n");
	    else if (refcat == TYCHO)
		printf ("No Tycho Stars Found\n");
	    else if (refcat == HIP)
		printf ("No Hipparcos Stars Found\n");
	    else if (refcat == ACT)
		printf ("No ACT Stars Found\n");
	    else
		printf ("No Stars Found\n");
	    }
	else {
	    if (refcat == GSC)
		printf ("GSC number ");
	    else if (refcat == USAC)
		printf ("USNO SA number ");
	    else if (refcat == USA1)
		printf ("USNO SA1 number");
	    else if (refcat == USA2)
		printf ("USNO SA2 number");
	    else if (refcat == UAC)
		printf ("USNO A number  ");
	    else if (refcat == UA1)
		printf ("USNO A1 number ");
	    else if (refcat == UA2)
		printf ("USNO A2 number ");
	    else if (refcat == UJC)
		printf (" UJ number    ");
	    else if (refcat == SAO)
		printf ("SAO number  ");
	    else if (refcat == PPM)
		printf ("PPM number  ");
	    else if (refcat == IRAS)
		printf ("IRAS number  ");
	    else if (refcat == TYCHO)
		printf ("Tycho number ");
	    else if (refcat == HIP)
		printf ("Hip number  ");
	    else if (refcat == ACT)
		printf ("ACT number  ");
	    else
		printf ("  Number  ");
	    if (sysout == WCS_B1950) {
		if (degout) {
		    if (eqout == 1950.0)
			printf ("  RA1950   Dec1950  ");
		    else
			printf ("RAB%7.2f DecB%7.2f  ", eqout, eqout);
		    }
		else {
		    if (eqout == 1950.0)
			printf ("RAB1950      DecB1950    ");
		    else
			printf ("RAB%7.2f   DecB%7.2f  ", eqout, eqout);
		    }
		}
	    else if (sysout == WCS_ECLIPTIC)
		printf ("Ecl Lon    Ecl Lat  ");
	    else if (sysout == WCS_GALACTIC)
		printf ("Gal Lon    Gal Lat  ");
	    else {
		if (degout) {
		    if (eqout == 2000.0)
			printf ("  RAJ2000   DecJ2000  ");
		    else
			printf ("RAJ%7.2f  DecJ%7.2f ", eqout, eqout);
		    }
		else {
		    if (eqout == 2000.0)
			printf (" RAJ2000       DecJ2000   ");
		    else
			printf ("RAJ%7.2f   DecJ%7.2f  ", eqout, eqout);
		    }
		}
	    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		refcat == UAC  || refcat == UA1  || refcat == UA2)
		printf ("MagB  MagR Plate");
	    else if (refcat == UJC)
		printf ("  Mag  Plate");
	    else if (refcat == GSC)
		printf (" Mag  Type");
	    else if (refcat == SAO || refcat == PPM || refcat == IRAS)
		printf (" Mag  Type");
	    else if (refcat == TYCHO || refcat == HIP || refcat == ACT)
		printf (" MagB  MagV  Type");
	    else if (refcat == TABCAT)
		printf (" Mag     Peak");
	    else
		printf ("  Mag");
	    if (ranges == NULL)
		printf ("   Arcsec\n");
	    else
		printf ("\n");
	    }
	}

    string[0] = (char) 0;
    for (i = 0; i < nbg; i++) {
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
	    CatNum (refcat, nndec, gnum[i], numstr);
	    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		refcat == UAC  || refcat == UA1  || refcat == UA2)
		sprintf (headline, "%s	%s	%s	%.1f	%.1f	%d",
		 numstr, rastr, decstr, gmb[i], gm[i], gc[i]);
	    else if (refcat == UJC)
	        sprintf (headline, "%s	%s	%s	%.2f	%d",
		 numstr, rastr, decstr, gm[i], gc[i]);
	    else if (refcat==SAO || refcat==PPM || refcat==IRAS ) {
		isp[0] = gc[i] / 1000;
		isp[1] = gc[i] % 1000;
	        sprintf (headline, "%s	%s	%s	%.2f	%2s",
		 numstr, rastr, decstr, gm[i], isp);
		}
	    else if (refcat == TYCHO || refcat == HIP || refcat == ACT) {
		isp[0] = gc[i] / 1000;
		isp[1] = gc[i] % 1000;
	        sprintf (headline, "%s	%s	%s	%.2f	%.2f	%2s",
		 numstr, rastr, decstr, gmb[i], gm[i], isp);
		}
	    else if (refcat == TABCAT)
	        sprintf (headline, "%s	%s	%s	%.2f	%d",
		 numstr, rastr, decstr, gm[i], gc[i]);
	    else if (refcat == BINCAT) {
		isp[0] = gc[i] / 1000;
		isp[1] = gc[i] % 1000;
	        sprintf (headline, "%s	%s	%s	%.2f	%2s",
			 numstr, rastr, decstr, gm[i], isp);
		}
	    else
	        sprintf (headline, "%s	%s	%s	%.2f",
		 numstr, rastr, decstr, gm[i]);

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
	    if (wfile)
		fprintf (fd, "%s\n", headline);
	    else if (tabout)
		printf ("%s\n", headline);
	    else {
		if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		    refcat == UAC  || refcat == UA1  || refcat == UA2)
		    sprintf (headline,"%s %s %s %5.1f %5.1f %4d ",
			numstr,rastr,decstr,gmb[i],gm[i],gc[i]);
		else if (refcat == UJC)
		    sprintf (headline,"%s %s %s %6.2f %4d",
			numstr, rastr, decstr, gm[i],gc[i]);
		else if (refcat == GSC)
		    sprintf (headline,"%s %s %s %6.2f %2d",
			numstr, rastr, decstr, gm[i],gc[i]);
		else if (refcat==SAO || refcat==PPM || refcat==IRAS ) {
		    isp[0] = gc[i] / 1000;
		    isp[1] = gc[i] % 1000;
		    sprintf (headline,"%s  %s %s %6.2f  %2s",
			numstr,rastr,decstr,gm[i],isp);
		    }
		else if (refcat == TYCHO || refcat == HIP || refcat == ACT) {
		    isp[0] = gc[i] / 1000;
		    isp[1] = gc[i] % 1000;
		    sprintf (headline,"%s %s %s %6.2f %6.2f  %2s",
			numstr,rastr,decstr,gmb[i],gm[i],isp);
		    }
		else if (refcat == TABCAT)
		    sprintf (headline,"%s %s %s %6.2f %7d",
			numstr, rastr, decstr, gm[i],gc[i]);
		else if (refcat == BINCAT) {
		    isp[0] = gc[i] / 1000;
		    isp[1] = gc[i] % 1000;
		    sprintf (headline,"%s %s %s %6.2f %2s",
			numstr, rastr, decstr, gm[i], isp);
		    }
		else {
		    sprintf (headline, "%s %s %s %6.2f",
			numstr, rastr, decstr, gm[i]);
		    }
		if (ranges == NULL) {
		    sprintf (temp, "  %7.2f", gdist);
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
		printf ("%s\n", headline);
		}
	    }
	}

	/* If searching more than one catalog, separate them with blank line */
	if (ncat > 0 && icat < ncat-1)
	    printf ("\n");

	/* Free memory used for object names in current catalog */
	if (gobj1 != NULL) {
	    for (i = 0; i < nbg; i++)
		if (gobj[i] != NULL) free (gobj[i]);
	    }
	}

    /* Close output file */
    if (wfile)
	fclose (fd);

    /* Free memory used for search results */
    if (gx) free ((char *)gx);
    if (gy) free ((char *)gy);
    if (gm) free ((char *)gm);
    if (gmb) free ((char *)gmb);
    if (gra) free ((char *)gra);
    if (gdec) free ((char *)gdec);
    if (gnum) free ((char *)gnum);
    if (gc) free ((char *)gc);
    if (gobj) free ((char *)gobj);

    return (nbg);
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
	*dra = (dra0 / cos (degrad (*cdec))) / 3600.0;
	*ddec = ddec0 / 3600.0;
	}
    else if (rad0 > 0.0) {
	*drad = 0.0;
	*ddec = rad0 / 3600.0;
	if (*cdec < 90.0 && *cdec > -90.0)
	    *dra = *ddec / cos (degrad (*cdec));
	else
	    *dra = 180.0;
	}
    else if (rad0 < 0.0) {
	*drad = -rad0 / 3600.0;
	*ddec = *drad;
	if (*cdec < 90.0 && *cdec > -90.0)
	    *dra = *ddec / cos (degrad (*cdec));
	else
	    *dra = 180.0;
	}
    else {
	if (verbose)
	    fprintf (stderr, "GetArea: Illegal radius, rad= %.5f\n",rad0);
	return (-1);
	}

    if (verbose) {
	ra2str (rstr, 32, *dra * 2.0, 2); 
	dec2str (dstr, 32, *ddec * 2.0, 2); 
	fprintf (stderr,"Area:    %s x %s\n", rstr, dstr);
	}

    return (0);
}


static void
SearchHead (refcatname,sys1,sys2,eq1,eq2,ep1,ep2,cra,cdec,dra,ddec,drad,nndec)

char	*refcatname;
int	sys1;
int	sys2;
double	eq1, eq2;
double	ep1, ep2;
double	cra, cdec;
double	dra, ddec;
double	drad;
int	nndec;
{
    double ra, dec;
    char rastr[32];
    char decstr[32];
    char cstr[16];
    char oform[16];
    int refcat;
    char *title[80];
    int sysref;
    double eqref, epref;

    ra = cra;
    dec = cdec;
    wcscon (sys1, sys2, eq1, eq2, &ra, &dec, ep2);
    if (sys2 == WCS_ECLIPTIC || sys2 == WCS_GALACTIC || degout0) {
	deg2str (rastr, 32, ra, 5);
	deg2str (decstr, 32, dec, 5);
	}
    else {
	ra2str (rastr, 32, ra, 3);
	dec2str (decstr, 32, dec, 2);
	}

    /* Set type of catalog being searched */
    if (!(refcat = RefCat (refcatname, title, &sysref, &eqref, &epref))) {
	fprintf (stderr,"ListCat: Catalog '%s' is missing\n",refcatname);
	return;
	}

    /* Label search center */
    if (objname) {
	sprintf (oform, "%%%ds", CatNumLen(refcat, nndec));
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
	else if (refcat == HIP)
	    printf ("Hip      ");
	else if (refcat == ACT)
	    printf ("ACT      ");
	else
	    printf ("%9s %s %s ", refcatname);
	}
    wcscstr (cstr, sys2, eq2, ep2);
    printf (" %s %s %s", rastr, decstr, cstr);
    if (drad != 0.0)
	printf ("r= %.2f", drad*3600.0);
    else
	printf ("+- %.2f", ddec*3600.0);
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
 */
