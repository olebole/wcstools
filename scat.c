/* File scat.c
 * September 25, 2000
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

static int ListCat ();
extern void setcenter();
static void SearchHead();
static int GetArea();

static int verbose = 0;		/* Verbose/debugging flag */
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
static int rasort = 0;		/* 1 to sort stars by brighness */
static int distsort = 0;	/* 1 to sort stars by distance from center */
static int debug = 0;		/* True for extra information */
static char *objname;		/* Object name for output */
static char *keyword;		/* Column to add to tab table output */
static char *progname;		/* Name of program as executed */
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
static struct StarCat *starcat[5]; /* Star catalog data structure */

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
    char *refcatname[5];	/* reference catalog names */
    char *refcatn;
    int srchtype;
    char *temp;
    int systemp = 0;		/* Input search coordinate system */
    char *ranges;
    int istar;
    char *blank;
    double epoch;
    char title[80];
    double eqout = 0.0;		/* Equinox for output coordinates */

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
		closest++;
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
		settabkey (keyword);
		ac--;
		break;

	    case 'l':	/* Print center and closest star on one line */
		oneline++;
		distsort++;
		closest++;
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
		if (ac > 1 && isnum (*(av+1)) &&
		    (ac < 4 || wcscsys (*(av+3)) < 1)) {
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
		    usage(progname);
		classd = (int) atof (*++av);
		setgsclass (classd);
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

	/* Read search center list from starbase tab table catalog */
	if (istab (listfile)) {
	    ranges = NULL;
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
		    ListCat (ncat, refcatname, ranges, eqout);
		    }
		tabcatclose (srchcat);
		}
	    }

	/* Read search center list from SAOTDC ASCII table catalog */
	else if (isacat (listfile)) {
	    ranges = NULL;
	    if (!(srchtype = RefCat (listfile, title, &syscoor, &eqcoor, &epoch))) {
		fprintf (stderr,"List catalog '%s' is missing\n", listfile);
		return (0);
		}
	    srchcat = ctgopen (listfile, srchtype);
	    if (srchcat != NULL) {
		srch = (struct Star *) calloc (1, sizeof (struct Star));
		for (istar = 1; istar <= srchcat->nstars; istar ++) {
		    if (ctgstar (istar, srchcat, srch)) {
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
		    ListCat (ncat, refcatname, ranges, eqout);
		    }
		ctgclose (srchcat, srchtype);
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
		    ListCat (ncat, refcatname, ranges, eqout);
		    }
		fclose (fd);
		}
	    }
	}
    else
	ListCat (ncat, refcatname, ranges, eqout);

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
usage (progname)

char *progname;
{
    if (version)
	exit (-1);
    if (strsrch (progname,"gsc") != NULL)
	fprintf (stderr,"Find HST Guide Stars in a region on the sky\n");
    else if (strsrch (progname,"ujc") != NULL)
	fprintf (stderr,"Find USNO J Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"uac") != NULL)
	fprintf (stderr,"Find USNO A Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"ua1") != NULL)
	fprintf (stderr,"Find USNO A-1.0 Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"ua2") != NULL)
	fprintf (stderr,"Find USNO A-2.0 Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"usac") != NULL)
	fprintf (stderr,"Find USNO SA Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"usa1") != NULL)
	fprintf (stderr,"Find USNO SA-1.0 Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"usa2") != NULL)
	fprintf (stderr,"Find USNO SA-2.0 Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"sao") != NULL)
	fprintf (stderr,"Find SAO Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"ppm") != NULL)
	fprintf (stderr,"Find PPM Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"ira") != NULL)
	fprintf (stderr,"Find IRAS Point Sources in a region on the sky\n");
    else if (strsrch (progname,"tyc") != NULL)
	fprintf (stderr,"Find Tycho Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"hip") != NULL)
	fprintf (stderr,"Find Hipparcos Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"act") != NULL)
	fprintf (stderr,"Find ACT Catalog stars in a region on the sky\n");
    else if (strsrch (progname,"bsc") != NULL)
	fprintf (stderr,"Find Bright Star Catalog stars in a region on the sky\n");
    else
	fprintf (stderr,"Find catalog stars in a region on the sky\n");
    fprintf (stderr,"Usage: [-abdehjlstvw] [-m [mag1] mag2] [-n num] [-r arcsec] [-u x y] ra dec sys or @list\n");
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
    fprintf(stderr,"  -u: Print X Y instead of number in front of non-tab entry\n");
    fprintf(stderr,"  -v: Verbose\n");
    fprintf(stderr,"  -w: Write tab table output file search[objname].[catalog]\n");
    fprintf(stderr,"  -x: GSC object type (0=stars 3=galaxies -1=all)\n");
    fprintf(stderr,"  -y: Epoch of output positions in years\n");
    exit (1);
}

#define TABMAX 64
static int nalloc = 0;

static int
ListCat (ncat, refcatname, ranges, eqout)

int	ncat;		/* Number of catalogs to search */
char	**refcatname;	/* Reference catalog name */
char	*ranges;	/* String with range of catalog numbers to list */
double	eqout;		/* Equinox for output coordinates */

{
    double cra, cdec;
    double epout;
    int sysref;		/* Coordinate system of reference catalog */
    double eqref;	/* Equinox of catalog to be searched */
    double epref;	/* Epoch of catalog to be searched */
    int ng;		/* Number of catalog stars */
    int ns;		/* Number of brightest catalog stars actually used */
    struct Range *range; /* Range of catalog numbers to list */
    int nfind;		/* Number of stars to find */
    int i, ngmax, mprop;
    double das, dds, drs;
    double gnmax;
    int degout;
    double maxnum;
    FILE *fd;
    char rastr[32], decstr[32];	/* coordinate strings */
    char numstr[32];	/* Catalog number */
    char cstr[32];	/* Coordinate system */
    double drad, dra, ddec, mag1, mag2;
    double gdist, da, dd, dec, gdmax;
    int nlog;
    int nmag;
    int typecol;
    int sysout;
    char headline[160];
    char filename[80];
    char title[80];
    char string[TABMAX];
    char temp[80];
    char isp[4];
    int refcat;		/* reference catalog switch */
    int icat, nndec, nnfld, nsfld;
    double date, time;
    int gcset;
    int ndist;
    double pra, pdec;
    void ep2dt();
    void PrintNum();
    int LenNum();

    /* Allocate space for returned catalog information */
    if (ranges != NULL) {
	int nfdef = 9;

	/* Allocate and fill list of numbers to read */
	range = RangeInit (ranges, nfdef);
	ngmax = rgetn (range);
	}
    else if (nstars > 0)
	ngmax = nstars;
    else
	ngmax = MAXCAT;

    /* Free currently allocated buffers if more entries are needed */
    if (nalloc > 0 && ngmax > nalloc) {
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
	}
    else if (nalloc == 0) {
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
	}
    else {
	for (i = 0; i < nalloc; i++) {
	    gm[i] = 99.0;
	    gmb[i] = 99.0;
	    gra[i] = 0.0;
	    gdec[i] = 0.0;
	    if (gpra != NULL) gpra[i] = 0.0;
	    if (gpdec != NULL) gpdec[i] = 0.0;
	    gnum[i] = 0.0;
	    gc[i] = 0;
	    gx[i] = 0.0;
	    gy[i] = 0.0;
	    }
	}

    nalloc = ngmax;
    if (gnum == NULL) {
	if (!(gnum = (double *) calloc (nalloc, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gnum\n",
		     nalloc*sizeof(double));
	}
    if (gra == NULL) {
	if (!(gra = (double *) calloc (nalloc, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gra\n",
		     nalloc*sizeof(double));
	}
    if (gdec == NULL) {
	if (!(gdec = (double *) calloc (nalloc, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gdec\n",
		     nalloc*sizeof(double));
	}
    if (gm == NULL) {
	if (!(gm = (double *) calloc (nalloc, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gm\n",
		     nalloc*sizeof(double));
	}
    if (gmb == NULL) {
	if (!(gmb = (double *) calloc (nalloc, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gmb\n",
		     nalloc*sizeof(double));
	}
    if (gc == NULL) {
	if (!(gc = (int *) calloc (nalloc, sizeof(int))))
	    fprintf (stderr, "Could not calloc %d bytes for gc\n",
		     nalloc*sizeof(int));
	}
    if (gx == NULL) {
	if (!(gx = (double *) calloc (nalloc, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gx\n",
		     nalloc*sizeof(double));
	}
    if (gy == NULL) {
	if (!(gy = (double *) calloc (nalloc, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gy\n",
		     nalloc*sizeof(double));
	}
    if (gobj == NULL) {
	if (!(gobj = (char **) calloc (nalloc, sizeof(char *))))
	    fprintf (stderr, "Could not calloc %d bytes for obj\n",
		     nalloc*sizeof(char *));
	}
    if (gnum==NULL || gra==NULL || gdec==NULL || gm==NULL || gmb==NULL ||
	gc==NULL || gx==NULL || gy==NULL || gobj==NULL) {
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

    /* Start of per catalog loop */
    for (icat = 0; icat < ncat; icat++) {
	nndec = 0;
	isp[2] = (char) 0;
	isp[3] = (char) 0;

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

	/* Figure out which catalog we are searching */
	if (!(refcat = RefCat (refcatname[icat], title, &sysref, &eqref, &epref))) {
	    fprintf (stderr,"ListCat: Catalog '%s' is missing\n", refcatname[icat]);
	    return (0);
	    }

	/* Find out whether catalog has proper motions */
	mprop = PropCat ();

	/* If it does, allocate memory for them, if not already done */
	if (mprop) {
	    if (gpra == NULL) {
		if (!(gpra = (double *) calloc (ngmax, sizeof(double))))
		    fprintf (stderr, "Could not calloc %d bytes for gpra\n",
			     ngmax*sizeof(double));
		}
	    if (gpdec == NULL) {
		if (!(gpdec = (double *) calloc (ngmax, sizeof(double))))
		    fprintf (stderr, "Could not calloc %d bytes for gpdec\n",
			     ngmax*sizeof(double));
		}
	    if (gpra == NULL || gpdec == NULL) {
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
	    }

	/* Set number of magnitudes */
	if (refcat==USAC || refcat==USA1 || refcat==USA2 ||
	    refcat==HIP || refcat==UAC  || refcat==UA1  ||
	    refcat==UA2 || refcat==TYCHO || refcat==TYCHO2 || refcat==ACT)
	    nmag = 2;
	else
	    nmag = 1;

	/* Set epoch from command line, search catalog, or searched catalog */
	epout = epoch0;
	if (epout == 0.0) {
	    if (srch!= NULL && srch->epoch != 0.0)
		epout = srch->epoch;
	    else if (sysout0 == WCS_J2000)
		epout = 2000.0;
	    else if (sysout0 == WCS_B1950)
		epout = 1950.0;
	    else
		epout = epref;
	    }

	/* Set output coordinate system from searched catalog if not set yet */
	if (!sysout)
	    sysout = sysref;
	if (!eqout)
	    eqout = eqref;
	if (!epout)
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
		      nfind, sysout, eqout, epout, match,
		      gnum, gra, gdec, gpra, gpdec, gm, gmb, gc, gobj, debug);

	    if (gobj[0] == NULL)
		gobj1 = NULL;
	    else
		gobj1 = gobj;
	    if (ng > nfind)
		ns = nfind;
	    else
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
		refcat == USAC || refcat == USA1   || refcat == USA2 ||
		refcat == UAC  || refcat == UA1    || refcat == UA2 ||
		refcat == ACT  || refcat == TYCHO2 || refcat == BSC)
		typecol = 1;
	    else if (refcat == GSC  || refcat == UJC || refcat == IRAS ||
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
		    if (refcat == GSC || refcat == UJC ||
			refcat == USAC || refcat == USA1 || refcat == USA2)
			rad0 = -900.0;
		    else if ( refcat == UAC  || refcat == UA1  || refcat == UA2)
			rad0 = -180.0;
		    else
			rad0 = -1800.0;
		    }
		else
		    rad0 = -10.0;
		}

	    /* Set limits from defaults and command line information */
	    if (GetArea (verbose,syscoor,sysout,eqout,epout,
		     &cra,&cdec,&dra,&ddec,&drad))
		return (0);

	    if (srch != NULL) {
		if ((srchcat->stnum == 0 || srchcat->stnum > 4) &&
		    strlen (srch->objname) > 0) {
		    if (objname == NULL)
			objname = (char *) calloc (1,32);
		    strcat (objname, srch->objname);
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
		SearchHead (refcatname[icat],sysout,sysout,eqout,eqout,epout,epout,
			cra,cdec,dra,ddec,drad,nnfld);
		if (sysout != syscoor && !closest)
		    SearchHead (refcatname[icat],sysout,syscoor,eqout,eqcoor,
			    epout,epout,cra,cdec,dra,ddec,drad,nnfld);
		if (sysref != syscoor && sysref != sysout && !closest)
		    SearchHead (refcatname[icat],sysout,sysref,eqout,eqref,
			    epout,epref,cra,cdec,dra,ddec,drad,nnfld);
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
	    ng = ctgread (refcatname[icat], refcat, distsort,
		      cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,mag2,
		      ngmax,&starcat[icat],
		      gnum,gra,gdec,gpra,gpdec,gm,gmb,gc,gobj,nlog);

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
		refcat == USAC || refcat == USA1   || refcat == USA2 ||
		refcat == UAC  || refcat == UA1    || refcat == UA2 ||
		refcat == ACT  || refcat == TYCHO2 || refcat == BSC)
		typecol = 1;
	    else if (refcat == GSC  || refcat == UJC || refcat == IRAS ||
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
		if (distsort)
		    XSortStars (gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gmb,gc,gobj1,ns);

		/* Sort star-like objects in image by right ascension */
		else if (rasort)
		    RASortStars (gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gmb,gc,gobj1,ns);

		/* Sort reference stars from brightest to faintest */
		else
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
			if (refcat == GSC)
			    printf ("catalog	HST GSC 1.1\n");
			else if (refcat == USAC)
			    printf ("catalog	USNO SA-1.0\n");
			else if (refcat == USA1)
			    printf ("catalog	USNO SA-1.0\n");
			else if (refcat == USA2)
			    printf ("catalog	USNO SA-2.0\n");
			else if (refcat == UAC)
			    printf ("catalog	USNO A-1.0\n");
			else if (refcat == UA1)
			    printf ("catalog	USNO A-1.0\n");
			else if (refcat == UA2)
			    printf ("catalog	USNO A-2.0\n");
			else if (refcat == UJC)
			    printf ("catalog	USNO UJ1.0\n");
			else if (refcat == SAO)
			    printf ("catalog	SAO\n");
			else if (refcat == PPM)
			    printf ("catalog	PPM\n");
			else if (refcat == IRAS)
			    printf ("catalog	IRAS Point Source\n");
			else if (refcat == TYCHO)
			    printf ("catalog	Tycho\n");
			else if (refcat == TYCHO2)
			    printf ("catalog	Tycho 2\n");
			else if (refcat == HIP)
			    printf ("catalog	Hipparcos\n");
			else if (refcat == ACT)
			    printf ("catalog	ACT\n");
    			else if (refcat == BSC)
			    printf ("catalog	Yale Bright Star Catalog");
			else
			    printf ("catalog	%s\n", refcatname[icat]);
			printf ("search	%s\n", listfile);
			printf ("equinox	%.1f\n", eqout);
			if (epout != eqout)
			    printf ("epoch	%.1f\n", epout);
			if (dra0 > 0.0) {
			    printf ("drasec	%.2f\n", das);
			    printf ("ddecsec	%.2f\n", dds);
			    }
			else if (rad0 < 0)
			    printf ("radsec	%.2f\n", drs);
			else if (rad0 > 0)
			    printf ("boxsec	%.2f\n", dds);
			if (mprop) {
			    printf ("pmunit	mas/yr\n");
			    }

		        /* Write column headings */
		        if (srch != NULL)
			    printf ("srch_id	");
		        printf ("srch_ra     	srch_dec    	");
		        if (srchcat->nmag != 0)
			    printf ("srch_mag	");
		        if (srchcat->nepoch)
			    printf ("epoch    	");
			strcpy (headline,"id                  ");
			headline[nnfld] = (char) 0;
			printf ("%s", headline);
			printf ("	ra          	dec        	");
			if (srchcat->nmag > 0)
			    printf ("mag    	");
			if (srchcat->nmag > 1)
			    printf ("mag2 	");
			if (srchcat->sptype != 0)
			    printf ("SpT 	");
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
			if (srch != NULL)
			    printf ("-------	");
			printf ("------------	------------	");
			if (srchcat->nmag > 0)
			    printf ("-----	");
			if (srchcat->nmag > 1)
			    printf ("-----	");
			if (srchcat->sptype != 0)
			    printf ("----	");
			if (srchcat->nepoch)
			    printf ("---------	");
			strcpy (headline,"--------------------");
			headline[nnfld] = (char) 0;
			printf ("%s", headline);
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
			else {
			    nsfld = CatNumLen (TXTCAT,srch->num,srchcat->nndec);
			    CatNum (TXTCAT,nsfld,srchcat->nndec,srch->num,numstr);
			    }
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
			if (srchcat->sptype != 0) {
			    if (tabout)
				printf ("	%2s", srch->isp);
			    else
				printf ("  %2s ", srch->isp);
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
		    else if (starcat[icat] != NULL)
			CatNum (refcat,nnfld,starcat[icat]->nndec,gnum[0],numstr);
		    else
			CatNum (refcat,nnfld,nndec,gnum[0],numstr);
		    if (degout) {
			num2str (rastr, gra[0], 10, 5);
			num2str (decstr, gdec[0], 10, 5);
			}
		    else {
			ra2str (rastr, 32, gra[0], 3);
			dec2str (decstr, 32, gdec[0], 2);
			}
		    if (tabout)
			printf ("	%s	%s	%s	%5.2f",
			        numstr, rastr, decstr, gm[0]);
		    else
			printf (" %s %s %s %5.2f",
			        numstr, rastr, decstr, gm[0]);
		    if (nmag > 1) {
			if (tabout)
			    printf ("	%5.2f", gmb[0]);
			else
			    printf (" %5.2f", gmb[0]);
			}
		    if (typecol == 1) {
			isp[0] = gc[i] / 1000;
			isp[1] = gc[i] % 1000;
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
			da = gra[0] - cra;
			dd = gdec[0] - cdec;
			gdist = sqrt (da*da + dd*dd);
			ndist = 5;
			}
		    else {
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
			printf ("%d / %d %s (between %.1f and %.1f)",
			    ns, ng, title, gm[0], gm[ns-1]);
		    else
			printf ("%d / %d %s (brighter than %.1f)",
		 	    ns, ng, title, gm[ns-1]);
		    printf ("\n");
		    }
		}
	    else {
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
		    fprintf (stderr, "%s:  cannot write file %s\n",
			     cpname, filename);
        	    return (0);
		    }
		}
            }

    /* Write heading */
    if (refcat == GSC)
	sprintf (headline, "catalog	HSTGSC1.1");
    else if (refcat == USAC)
	sprintf (headline, "catalog	USNO SA-1.0");
    else if (refcat == USA1)
	sprintf (headline, "catalog	USNO SA-1.0");
    else if (refcat == USA2)
	sprintf (headline, "catalog	USNO SA-2.0");
    else if (refcat == UAC)
	sprintf (headline, "catalog	USNO A-1.0");
    else if (refcat == UA1)
	sprintf (headline, "catalog	USNO A-1.0");
    else if (refcat == UA2)
	sprintf (headline, "catalog	USNO A-2.0");
    else if (refcat == UJC)
	sprintf (headline, "catalog	USNO UJ1.0");
    else if (refcat == SAO)
	sprintf (headline, "catalog	SAO");
    else if (refcat == PPM)
	sprintf (headline, "catalog	PPM");
    else if (refcat == IRAS)
	sprintf (headline, "catalog	IRAS Point Source");
    else if (refcat == TYCHO)
	sprintf (headline, "catalog	Tycho");
    else if (refcat == TYCHO2)
	sprintf (headline, "catalog	Tycho 2");
    else if (refcat == HIP)
	sprintf (headline, "catalog	Hipparcos");
    else if (refcat == ACT)
	sprintf (headline, "catalog	ACT");
    else if (refcat == BSC)
	sprintf (headline, "catalog	Yale Bright Star Catalog");
    else
	sprintf (headline, "catalog	%s", refcatname[icat]);
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (!ranges) {
	ra2str (rastr, 32, cra, 3);
	if (wfile)
	    fprintf (fd, "ra	%s\n", rastr);
	if (tabout)
	    printf ("ra	%s\n", rastr);

	dec2str (decstr, 32, cdec, 2);
	if (wfile)
	    fprintf (fd, "dec	%s\n", decstr);
	if (tabout)
	    printf ("dec	%s\n", decstr);

	}
    if (wfile)
	fprintf (fd, "equinox	%.1f\n", eqout);
    if (tabout)
	printf ("equinox	%.1f\n", eqout);

    if (epout != eqout) {
	if (wfile)
	    fprintf (fd, "epoch	%.1f\n", epout);
	if (tabout)
	    printf ("epoch	%.1f\n", epout);
	}
    if (mprop) {
	if (wfile) {
	    fprintf (fd, "rpmunit	mas/yr\n");
	    fprintf (fd, "dpmunit	mas/yr\n");
	    }
	if (tabout) {
	    printf ("rpmunit	mas/yr\n");
	    printf ("dpmunit	mas/yr\n");
	    }
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
	    printf ("drasec	%.2f\n", das);
	    printf ("ddecsec	%.2f\n", dds);
	    }
	else if (rad0 < 0)
	    printf ("radsec	%.2f\n", drs);
	else if (rad0 > 0)
	    printf ("boxsec	%.2f\n",dds);
	}

    if (distsort) {
	if (wfile)
	    fprintf (fd, "distsort	T\n");
	if (tabout)
	    printf ("distsort	T\n");
	}
    else if (rasort) {
	if (wfile)
	    fprintf (fd, "rasort	T\n");
	if (tabout)
	    printf ("rasort	T\n");
	}
    if (wfile)
    if (tabout)

    /* Print column headings */
    if (refcat == ACT)
	strcpy (headline, "act_id       ");
    else if (refcat == BSC)
	strcpy (headline, "bsc_id       ");
    else if (refcat == GSC)
	strcpy (headline, "gsc_id       ");
    else if (refcat == USAC)
	strcpy (headline,"usac_id       ");
    else if (refcat == USA1)
	strcpy (headline,"usa1_id       ");
    else if (refcat == USA2)
	strcpy (headline,"usa2_id       ");
    else if (refcat == UAC)
	strcpy (headline,"usnoa_id      ");
    else if (refcat == UA1)
	strcpy (headline,"usnoa1_id     ");
    else if (refcat == UA2)
	strcpy (headline,"usnoa2_id     ");
    else if (refcat == UJC)
	strcpy (headline,"usnoj_id      ");
    else if (refcat == SAO)
	strcpy (headline,"sao_id        ");
    else if (refcat == PPM)
	strcpy (headline,"ppm_id        ");
    else if (refcat == IRAS)
	strcpy (headline,"iras_id       ");
    else if (refcat == TYCHO)
	strcpy (headline,"tycho_id      ");
    else if (refcat == TYCHO2)
	strcpy (headline,"tycho2_id     ");
    else if (refcat == HIP)
	strcpy (headline,"hip_id        ");
    else
	strcpy (headline,"id            ");
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
	strcat (headline,"	magb	magr	type");
    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==HIP || refcat==ACT)
	strcat (headline,"	magb	magv	type");
    else if (refcat == GSC)
	strcat (headline,"	mag	class");
    else if (refcat == UJC)
	strcat (headline,"	mag	plate");
    else
	strcat (headline,"	mag");
    if (refcat == SAO || refcat == PPM || (refcat == TABCAT && gcset))
	strcat (headline,"	type");
    if (mprop)
	strcat (headline,"	Ux    	Uy    ");
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
    if (tabout)
	printf ("%s\n", headline);

    strcpy (headline, "---------------------");
    headline[nnfld] = (char) 0;
    strcat (headline,"	------------	------------");
    if (nmag == 2)
	strcat (headline,"	-----	-----");
    else
	strcat (headline,"	-----");
    if (typecol)
	strcat (headline,"	-----");
    if (refcat == SAO || refcat == PPM || (refcat == TABCAT && gcset))
	strcat (headline,"	----");
    if (mprop)
	strcat (headline,"	------	------");
    if (ranges == NULL)
	strcat (headline, "	------");
    if (refcat == TABCAT && keyword != NULL)
	strcat (headline,"	------");
    if (printxy)
	strcat (headline, "	-------	-------");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (ranges == NULL)
	strcat (headline,"	------");
    if (printhead) {
	if (ng == 0)
	    printf ("No %s Stars Found\n", title);
	else {
	    if (!printxy) {
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
		    printf ("SAO number ");
		else if (refcat == PPM)
		    printf ("PPM number ");
		else if (refcat == BSC)
		    printf ("BSC number ");
		else if (refcat == IRAS)
		    printf ("IRAS number ");
		else if (refcat == TYCHO)
		    printf ("Tycho number ");
		else if (refcat == TYCHO2)
		    printf ("Tycho2 num  ");
		else if (refcat == HIP)
		    printf ("Hip number  ");
		else if (refcat == ACT)
		    printf ("ACT number  ");
		else
		    printf ("Number   ");
		}
	    else
		printf ("  X     Y   ");
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
			printf ("  RA2000   Dec2000  ");
		    else
			printf ("RAJ%7.2f  DecJ%7.2f ", eqout, eqout);
		    }
		else {
		    if (eqout == 2000.0)
			printf (" RA2000       Dec2000   ");
		    else
			printf ("RAJ%7.2f   DecJ%7.2f  ", eqout, eqout);
		    }
		}
	    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		refcat == UAC  || refcat == UA1  || refcat == UA2)
		printf ("  MagB  MagR Type");
	    else if (refcat == UJC)
		printf ("  Mag  Plate");
	    else if (refcat == GSC)
		printf ("   Mag Class");
	    else if (refcat == SAO || refcat == PPM ||
		     refcat == IRAS || refcat == BSC)
		printf ("   Mag  Type");
	    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==HIP || refcat==ACT)
		printf ("  MagB   MagV Type");
	    else if (refcat == TABCAT && gcset)
		printf (" Mag     Peak");
	    else
		printf ("  Mag");
	    if (ranges == NULL)
		printf ("  Arcsec");
	    printf ("\n");
	    }
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
    for (i = 0; i < ns; i++) {
	if (typecol == 1) {
	    isp[0] = gc[i] / 1000;
	    isp[1] = gc[i] % 1000;
	    if (isp[0] == ' ' && isp[1] == ' ') {
		isp[0] = '_';
		isp[1] = '_';
		}
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
		pra = (gpra[i] * 3600.0 * 1000.0) * cosdeg (gdec[i]);
		pdec = gpdec[i] * 3600.0 * 1000.0;
		}
	    CatNum (refcat, nnfld, nndec, gnum[i], numstr);
	    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		refcat == UAC  || refcat == UA1  || refcat == UA2)
		sprintf (headline, "%s	%s	%s	%.1f	%.1f	%2s",
		 numstr, rastr, decstr, gmb[i], gm[i], isp);
	    else if (refcat == UJC || refcat == GSC)
	        sprintf (headline, "%s	%s	%s	%.2f	%d",
		 numstr, rastr, decstr, gm[i], gc[i]);
	    else if (refcat==SAO || refcat==PPM ||
		     refcat==IRAS || refcat == BSC) {
	        sprintf (headline, "%s	%s	%s	%.2f	%2s",
		 numstr, rastr, decstr, gm[i], isp);
		}
	    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==HIP || refcat==ACT) {
	        sprintf (headline, "%s	%s	%s	%.2f	%.2f	%2s",
		 numstr, rastr, decstr, gmb[i], gm[i], isp);
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
	        sprintf (temp, "	%6.2f	%6.2f", pra, pdec);
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
	    else if (tabout)
		printf ("%s\n", headline);
	    else {
		if (printxy) {
		    strcpy (numstr, xstr);
		    strcat (numstr, " ");
		    strcat (numstr, ystr);
		    }
		if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		    refcat == UAC  || refcat == UA1  || refcat == UA2)
		    sprintf (headline,"%s %s %s %5.1f %5.1f  %2s ",
			numstr,rastr,decstr,gmb[i],gm[i],isp);
		else if (refcat == UJC || refcat == GSC)
		    sprintf (headline,"%s %s %s %6.2f %4d",
			numstr, rastr, decstr, gm[i],gc[i]);
		else if (refcat == GSC)
		    sprintf (headline,"%s %s %s %6.2f %2d",
			numstr, rastr, decstr, gm[i],gc[i]);
		else if (refcat==SAO || refcat==PPM ||
			 refcat==IRAS || refcat == BSC) {
		    sprintf (headline,"  %s  %s %s %6.2f  %2s",
			numstr,rastr,decstr,gm[i],isp);
		    }
		else if (refcat==TYCHO || refcat==TYCHO2 || refcat==HIP || refcat==ACT) {
		    sprintf (headline,"%s %s %s %6.2f %6.2f  %2s",
			numstr,rastr,decstr,gmb[i],gm[i],isp);
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
	    *dra = (dra0 / cos (degrad (*cdec))) / 3600.0;
	    *ddec = ddec0 / 3600.0;
	    }
	}
    else if (rad0 > 0.0) {
	*drad = 0.0;
	if (syscoor == WCS_XY) {
	    *dra = rad0;
	    *ddec = rad0;
	    }
	else {
	    *ddec = rad0 / 3600.0;
	    if (*cdec < 90.0 && *cdec > -90.0)
		*dra = *ddec / cos (degrad (*cdec));
	    else
		*dra = 180.0;
	    }
	}
    else if (rad0 < 0.0) {
	if (syscoor == WCS_XY) {
	    *drad = -rad0;
	    *dra = *drad;
	    *ddec = *drad;
	    }
	else {
	    *drad = -rad0 / 3600.0;
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
SearchHead (refcatname,sys1,sys2,eq1,eq2,ep1,ep2,cra,cdec,dra,ddec,drad,nnfld)

char	*refcatname;	/* Name of catalog */
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
    char *title[80];
    int sysref;
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
    if (!(refcat = RefCat (refcatname, title, &sysref, &eqref, &epref))) {
	fprintf (stderr,"ListCat: Catalog '%s' is missing\n",refcatname);
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
	    printf ("%9s ", refcatname);
	}
    wcscstr (cstr, sys2, eq2, ep2);
    printf (" %s %s %s", rastr, decstr, cstr);
    if (drad != 0.0)
	printf (" r= %.2f", drad*3600.0);
    else
	printf (" +- %.2f", ddec*3600.0);
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
 */
