/* File scat.c
 * July 9, 1998
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
#include "fitshead.h"
#include "wcs.h"
#include "libwcs/lwcs.h"

#define GSC	1	/* refcat value for HST Guide Star Catalog */
#define UJC	2	/* refcat value for USNO UJ Star Catalog */
#define UAC	3	/* refcat value for USNO A-1.0 Star Catalog */
#define USAC	4	/* refcat value for USNO SA-1.0 Star Catalog */

#define MAXREF 100

static void usage();

extern int gscread();
extern int uacread();
extern int ujcread();
extern int gscrnum();
extern int uacrnum();
extern int ujcrnum();
static int ListCat ();
extern void XSortStars ();
extern void RASortStars ();
extern void MagSortStars ();
extern void setcenter();
extern void setradius();
static void SearchHead();
static int GetArea();
extern double wcsdist();
extern int tabopen();
extern int tabgetk();
extern void tabclose();
extern char *tabstar();

static int verbose = 0;		/* Verbose/debugging flag */
static int wfile = 0;		/* True to print output file */
static int classd = -1;		/* Guide Star Catalog object classes */
static int uplate = 0;		/* UJ Catalog plate number to use */
static double maglim1 = MAGLIM1; /* Catalog bright magnitude limit */
static double maglim2 = MAGLIM2; /* Catalog faint magnitude limit */
static int sysout = -1;		/* Output coordinate system */
static double eqref = 2000.0;	/* Equinox of catalog to be searched */
static double eqcoor = 2000.0;	/* Equinox of search center */
static double eqout = 2000.0;	/* Equinox for output coordinates */
static int degout0 = 0;		/* 1 if degrees output instead of hms */
static double ra0 = -99.0;	/* Initial center RA in degrees */
static double dec0 = -99.0;	/* Initial center Dec in degrees */
static double rad0 = 10.0;	/* Search box radius */
static double epoch = 2000.0;	/* Epoch for coordinates */
static int syscoor = -1;	/* Input search coordinate system */
static int nstars = 0;		/* Number of brightest stars to list */
static int printhead = 0;	/* 1 to print table heading */
static int tabout = 0;		/* 1 for tab table to standard output */
static int rasort = 0;		/* 1 to sort stars by brighness */
static int distsort = 0;	/* 1 to sort stars by distance from center */
static int debug = 0;		/* True for extra information */
static char *objname = NULL;	/* Object name for output */
static int refcat = GSC;	/* reference catalog switch */
static char refcatname[32]="GSC";	/* reference catalog name */
static int refcat2 = 0;		/* Second reference catalog switch */
static char *keyword = NULL;	/* Column to add to tab table output */
static char *progname;		/* Name of program as executed */
static char progpath[128];

main (ac, av)
int ac;
char **av;
{
    char *str, *str1;
    char rastr[16];
    char decstr[16];
    int nfind = 0;
    int i, nc;
    double *snum;
    snum = NULL;

    /* Check name used to execute programe and set catalog name accordingly */
    strcpy (progpath, av[0]);
    progname = progpath;
    for (i = strlen (progpath); i > -1; i--) {
	if (progpath[i] > 95 && progpath[i] < 122)
	    progpath[i] = progpath[i] - 32;
	if (progpath[i] == '/') {
	    progname = progpath + i;
	    break;
	    }
	}
    if (strsrch (progname,"GSC") != NULL)
	refcat = GSC;
    else if (strsrch (progname,"UAC") != NULL)
	refcat = UAC;
    else if (strsrch (progname,"USAC") != NULL)
	refcat = USAC;
    else if (strsrch (progname,"UJC") != NULL)
	refcat = UJC;
    else
	refcat = -1;

    if (ac == 1)
        usage ();

    /* crack arguments */
    for (av++; --ac > 0; av++) {

	/* Set RA, Dec, and equinox if WCS-generated argument */
	if (strsrch (*av,":") != NULL) {
	    if (ac < 3)
		usage();
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

	/* Set star number if getting just one star */
	else if (strsrch (*av,".") != NULL) {
	    if (snum == NULL)
		snum = (double *) malloc (sizeof(double));
	    else
		snum = (double *) realloc (snum, (nfind+1)*sizeof(double));
	    snum[nfind] = atof (*av);
	    nfind++;
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
    		setradius (60.0);
		break;

    	    case 'b':	/* initial coordinates on command line in B1950 */
		str1 = *(av+1);
		nc = (int)str1[0];
		if (syscoor >= 0 || *(str+1) || nc < 47 || nc > 58) {
		    sysout = WCS_B1950;
		    eqout = 1950.0;
		    }
		else if (ac < 3)
		    usage ();
		else {
		    syscoor = WCS_B1950;
		    eqcoor = 1950.0;
		    sysout = WCS_B1950;
		    eqout = 1950.0;
		    strcpy (rastr, *++av);
		    ac--;
		    strcpy (decstr, *++av);
		    ac--;
		    ra0 = str2ra (rastr);
		    dec0 = str2dec (decstr);
		    if (ac > 1) {
			if ((syscoor = wcscsys (*(av+1))) > -1) {
			    ac--;
			    av++;
			    eqcoor = wcsceq (*av);
			    }
			}
		    }
    		break;

	    case 'c':       /* Set reference catalog */
		if (ac < 2)
		    usage();
		strcpy (refcatname, *++av);
		if (strncmp(refcatname,"gs",2)==0 ||
		    strncmp (refcatname,"GS",2)== 0)
		    refcat = GSC;
		else if (strncmp(refcatname,"us",2)==0 ||
		    strncmp(refcatname,"US",2)==0)
		    refcat = USAC;
		else if (strncmp(refcatname,"ua",2)==0 ||
		    strncmp(refcatname,"UA",2)==0)
		    refcat = UAC;
		else if (strncmp(refcatname,"uj",2)==0 ||
		    strncmp(refcatname,"UJ",2)==0)
		    refcat = UJC;
		else
		    refcat = 0;
		ac--;
		break;

	    case 'd':	/* output in degrees instead of sexagesimal */
		degout0++;
		break;

	    case 'e':	/* Set ecliptic coordinate output and optional center */
		str1 = *(av+1);
		nc = (int)str1[0];
		if (syscoor >= 0 || *(str+1) || nc < 47 || nc > 58)
		    sysout = WCS_ECLIPTIC;
		else if (ac < 3)
		    usage ();
		else {
		    syscoor = WCS_ECLIPTIC;
		    strcpy (rastr, *++av);
		    ac--;
		    strcpy (decstr, *++av);
		    ac--;
		    ra0 = atof (rastr);
		    dec0 = atof (decstr);
		    if (ac > 1) {
			if ((syscoor = wcscsys (*(av+1))) > -1) {
			    ac--;
			    av++;
			    }
			else
			    syscoor = WCS_ECLIPTIC;
			}
		    }
		break;

	    case 'g':	/* Set galactic coordinate output and optional center */
		if (ac < 2 || *(str+1) != 0 || nc < 47 || nc > 58)
		    sysout = WCS_GALACTIC;
		else if (ac < 3)
		    usage ();
		else {
		    syscoor = WCS_GALACTIC;
		    strcpy (rastr, *++av);
		    ac--;
		    strcpy (decstr, *++av);
		    ac--;
		    ra0 = atof (rastr);
		    dec0 = atof (decstr);
		    if (ac > 1) {
			if ((syscoor = wcscsys (*(av+1))) > -1) {
			    ac--;
			    av++;
			    }
			else
			    syscoor = WCS_GALACTIC;
			}
		    }
		break;

	    case 'h':	/* ouput descriptive header */
		printhead++;
		break;

    	    case 'j':	/* center coordinates on command line in J2000 */
		str1 = *(av+1);
		nc = (int)str1[0];
		if (syscoor >= 0 || *(str+1) || nc < 47 || nc > 58) {
		    sysout = WCS_J2000;
		    eqout = 2000.0;
		    }
		else if (ac < 3)
		    usage ();
		else {
		    sysout = WCS_J2000;
		    eqout = 2000.0;
		    syscoor = WCS_J2000;
		    eqcoor = 2000.0;
		    strcpy (rastr, *++av);
		    ac--;
		    strcpy (decstr, *++av);
		    ac--;
		    ra0 = str2ra (rastr);
		    dec0 = str2dec (decstr);
		    if (ac > 1) {
			if ((syscoor = wcscsys (*(av+1))) > -1) {
			    ac--;
			    av++;
			    eqcoor = wcsceq (*av);
			    }
			}
		    }
    		break;

	    case 'k':	/* Keyword (column) to add to output from tab table */
		if (ac < 2)
		    usage ();
		keyword = *++av;
		ac--;
		break;

	    case 'm':	/* Magnitude limit */
		if (ac < 2)
		    usage ();
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
		    usage ();
		nstars = atoi (*++av);
		ac--;
		break;

	    case 'o':	/* Object name */
		if (ac < 2)
		    usage ();
		objname = *++av;
		ac--;
		break;

	    case 'p':	/* Sort by distance from center */
		distsort++;
		break;

    	    case 'q':	/* Output equinox in years */
    		if (ac < 2)
    		    usage();
    		eqout = atof (*++av);
    		ac--;
    		break;

    	    case 'r':	/* Box radius in arcseconds */
    		if (ac < 2)
    		    usage();
    		rad0 = atof (*++av);
    		ac--;
    		break;

	    case 's':	/* sort by RA */
		rasort = 1;
		break;

	    case 't':	/* tab table to stdout */
		tabout = 1;
		break;

	    case 'u':       /* USNO Catalog plate number */
		if (ac < 2)
		    usage();
		uplate = (int) atof (*++av);
		ac--;
		break;

    	    case 'w':	/* write output file */
    		wfile++;
    		break;

	    case 'y':	/* Set output coordinate epoch */
		if (ac < 2)
		    usage();
		epoch = atof (*++av);
		ac--;
		break;

    	    case '2':	/* check second catalog */
		if (ac < 2)
		    usage();
		strcpy (refcatname, *++av);
		if (strncmp(refcatname,"gs",2)==0 ||
		    strncmp (refcatname,"GS",2)== 0)
		    refcat2 = GSC;
		else if (strncmp(refcatname,"us",2)==0 ||
		    strncmp(refcatname,"US",2)==0)
		    refcat2 = USAC;
		else if (strncmp(refcatname,"ua",2)==0 ||
		    strncmp(refcatname,"UA",2)==0)
		    refcat2 = UAC;
		else if (strncmp(refcatname,"uj",2)==0 ||
		    strncmp(refcatname,"UJ",2)==0)
		    refcat2 = UJC;
		else
		    refcat2 = 0;
		ac--;
    		break;

	    default:
		usage ();
		break;
	    }
	}
    }

    if (ListCat (nfind, snum) < 1 && refcat2 > 0) {
	refcat = refcat2;
	(void)ListCat (nfind, snum);
	}

    return (0);
}

static void
usage ()
{
    if (refcat == GSC)
	fprintf (stderr,"Find HST Guide Stars in a square on the sky\n");
    else if (refcat == UJC)
	fprintf (stderr,"Find USNO J Catalog stars in a square on the sky\n");
    else if (refcat == UAC)
	fprintf (stderr,"Find USNO A 1.0 Catalog stars in a square on the sky\n");
    else if (refcat == USAC)
	fprintf (stderr,"Find USNO SA 1.0 Catalog stars in a square on the sky\n");
    else
	fprintf (stderr,"Find catalog stars in a square on the sky\n");
    fprintf (stderr,"Usage: [-adhstvw] [-m [mag1] mag2] [-e sys] [-n num] [-r arcsec] [-b][-j] ra dec\n");
    fprintf(stderr,"  -a: List single closest catalog source\n");
    fprintf(stderr,"  -b: Output B1950 (FK4) coordinates with optional center\n");
    fprintf(stderr,"  -c: Reference catalog (gsc, uac, usac, or tab table file\n");
    fprintf(stderr,"  -d: Sort by distance from center instead of flux\n");
    fprintf(stderr,"  -e: Output ecliptic coordinates with optional center)\n");
    fprintf(stderr,"  -g: Output galactic coordinates with optional center)\n");
    fprintf(stderr,"  -h: Print heading, else do not \n");
    fprintf(stderr,"  -j: Output J2000 (FK5) coordinates with optional center\n");
    fprintf(stderr,"  -k: Add this keyword to output from tab table search\n");
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
    fprintf(stderr,"  -w: Write tab table output file imagename.gsc\n");
    fprintf(stderr,"  -x: GSC object type (0=stars 3=galaxies -1=all)\n");
    fprintf(stderr,"  -y: Epoch of output positions in years\n");
    exit (1);
}

#define TABMAX 64

static int
ListCat (nfind, snum)

int nfind;		/* Number of stars to find */
double *snum;		/* Catalog numbers */

{
    double *gnum=0;	/* Catalog star numbers */
    double *gra=0;	/* Catalog star right ascensions, rads */
    double *gdec=0;	/* Catalog star declinations rads */
    double *gm=0;	/* Catalog magnitudes */
    double *gmb=0;	/* Catalog B magnitudes */
    double *gx=0;	/* Catalog star X positions on image */
    double *gy=0;	/* Catalog star Y positions on image */
    int *gc=0;		/* Catalog star object classes */
    double cra, cdec;
    int sysref;		/* Coordinate system of reference catalog */
    int ng;		/* Number of catalog stars */
    int nbg;		/* Number of brightest catalog stars actually used */
    int i, ngmax, nbytes;
    int degout;
    FILE *fd;
    char rastr[32], decstr[32];	/* coordinate strings */
    double drad, dra, ddec, ra1, dec1, mag1, mag2;
    double epref;
    double mag;
    int offscale, nlog;
    char headline[160];
    char filename[80];
    char title[80];
    char *tabline;
    char string[TABMAX];
    int ntab;

    if (verbose || printhead) {
	if (nstars == 1)
	else
	}
    tabline = 0;
    if (sysout < 0) {
	sysout = syscoor;
	eqout = eqcoor;
	}
    sysref = WCS_J2000;
    epref = 2000.0;
    if (debug)
	nlog = 1;
    else if (verbose)
	nlog = 100;
    else
	nlog = 0;

    if (refcat == GSC)
	strcpy (title, "HST Guide Stars");
    else if (refcat == UAC)
	strcpy (title, "USNO A Catalog Stars");
    else if (refcat == USAC)
	strcpy (title, "USNO SA Catalog Stars");
    else if (refcat == UJC)
	strcpy (title, "USNO J Catalog Stars");
    else if (refcatname[0] > 0)
	sprintf (title, "%s Catalog Stars", refcatname);
    else {
	fprintf (stderr,"ListCat: Catalog name is missing\n");
	return (0);
	}

    /* Find stars specified by number */
    if (snum != NULL) {
	nbytes = nfind * sizeof (double);
	if (!(gnum = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gnum\n", nbytes);
	else {
	    for (i = 0; i < nfind; i++)
		gnum[i] = snum[i];
	    }
	if (!(gra = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gra\n", nbytes);
	if (!(gdec = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gdec\n", nbytes);
	if (!(gm = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gm\n", nbytes);
	if (!(gmb = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gmb\n", nbytes);
	if (!(gc = (int *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gc\n", nbytes);
	if (!(gx = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gx\n", nbytes);
	if (!(gy = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gy\n", nbytes);
	if (!gnum || !gra || !gdec || !gm || !gmb || !gc || !gx || !gy) {
	    if (gm) free ((char *)gm);
	    if (gmb) free ((char *)gmb);
	    if (gra) free ((char *)gra);
	    if (gdec) free ((char *)gdec);
	    if (gnum) free ((char *)gnum);
	    if (gc) free ((char *)gc);
	    if (gx) free ((char *)gx);
	    if (gy) free ((char *)gy);
	    return (0);
	    }
	wfile = 0;

	/* Find the specified catalog stars */
	if (refcat == GSC)
	    nbg = gscrnum (nfind,gnum,gra,gdec,gm,gc,nlog);
	else if (refcat == USAC)
	    nbg = usarnum (nfind,gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == UAC)
	    nbg = uacrnum (nfind,gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == UJC)
	    nbg = ujcrnum (nfind,gnum,gra,gdec,gm,gc,nlog);
	else
	    nbg = tabrnum (refcatname,nfind,gnum,gra,gdec,gm,gc,debug);
	for (i = 0; i < nbg; i++ ) {
	    wcscon (sysref, sysout, eqref, eqout, &gra[i],&gdec[i], 2000.0);
	    gx[i] = 0.0;
	    gy[i] = 1.0;
	    }
	}

    /* Find stars specified by location */
    else {


	/* Set limits from defaults and command line information */
	if (GetArea (verbose,syscoor,sysref,epoch,&cra,&cdec,&dra,&ddec,&drad))
	    return (0);

	/* Print search center and size in input and output coordinates */
	if (verbose || printhead) {
	    SearchHead (sysref,syscoor,eqref,eqcoor,cra,cdec,dra,ddec,drad);
	    if (sysref != syscoor)
		SearchHead (sysref,sysref,eqref,eqref,cra,cdec,dra,ddec,drad);
	    if (sysout != syscoor && sysout != sysref)
		SearchHead (sysref,sysout,eqref,eqout,cra,cdec,dra,ddec,drad);
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

	if (nstars > MAXREF)
	    ngmax = nstars;
	else
	    ngmax = MAXREF;
	nbytes = ngmax * sizeof (double);
	if (!(gnum = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gnum\n", nbytes);
	if (!(gra = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gra\n", nbytes);
	if (!(gdec = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gdec\n", nbytes);
	if (!(gm = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gm\n", nbytes);
	if (!(gmb = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gmb\n", nbytes);
	if (!gm)
	    fprintf (stderr, "Could not malloc %d bytes for gm\n", nbytes);
	if (!(gc = (int *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gc\n", nbytes);
	if (!gnum || !gra || !gdec || !gm || !gmb || !gc) {
	    if (gm) free ((char *)gm);
	    if (gmb) free ((char *)gmb);
	    if (gra) free ((char *)gra);
	    if (gdec) free ((char *)gdec);
	    if (gnum) free ((char *)gnum);
	    if (gc) free ((char *)gc);
	    return (0);
	    }

	/* Find the nearby reference stars, in ra/dec */
	if (refcat == GSC)
	    ng = gscread (cra,cdec,dra,ddec,drad,mag1,mag2,classd,ngmax,
			  gnum,gra,gdec,gm,gc,nlog);
	else if (refcat == USAC)
	    ng = usaread (cra,cdec,dra,ddec,drad,mag1,mag2,uplate,ngmax,
			  gnum,gra,gdec,gm,gmb,gc,debug);
	else if (refcat == UAC)
	    ng = uacread (cra,cdec,dra,ddec,drad,mag1,mag2,uplate,ngmax,
			  gnum,gra,gdec,gm,gmb,gc,debug);
	else if (refcat == UJC)
	    ng = ujcread (cra,cdec,dra,ddec,drad,mag1,mag2,uplate,ngmax,
			  gnum,gra,gdec,gm,gc,debug);
	else
	    ng = tabread (refcatname,cra,cdec,dra,ddec,drad,mag1,mag2,ngmax,
			  gnum,gra,gdec,gm,gc,debug);

	if (!(gx = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gx\n", nbytes);
	if (!(gy = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gy\n", nbytes);
	if (!gx || !gy) {
	    if (gm) free ((char *)gm);
	    if (gmb) free ((char *)gmb);
	    if (gra) free ((char *)gra);
	    if (gdec) free ((char *)gdec);
	    if (gnum) free ((char *)gnum);
	    if (gc) free ((char *)gc);
	    if (gx) free ((char *)gx);
	    if (gy) free ((char *)gy);
	    return (0);
	    }

	/* Compute distance from star to search center */
	for (i = 0; i < ng; i++ ) {
	    offscale = 0;
	    gx[i] = wcsdist (cra, cdec, gra[i], gdec[i]);
	    gy[i] = 1.0;
	    }

	/* Sort reference stars from closest to furthest */
	if (distsort)
	    XSortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, ng);

	/* Sort star-like objects in image by right ascension */
	else if (rasort)
	    RASortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, nbg);

	/* Sort reference stars from brightest to faintest */
	else
	    MagSortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, ng);

	/* List the brightest or closest MAXSTARS reference stars */
	if (nstars > 0 && ng > nstars) {
	    nbg = nstars;
	    if (verbose || printhead) {
		if (distsort) {
		    if (ng > 1)
			printf ("Closest %d / %d %s (closer than %.2f arcsec)\n",
				nbg, ng, title, 3600.0*gx[nbg-1]);
		    else
			printf ("Closest of %d %s\n",ng, title);
		    }
		else if (maglim1 > 0.0)
		    printf ("%d / %d %s (between %.1f and %.1f)\n",
			    nbg, ng, title, gm[0], gm[nbg-1]);
		else
		    printf ("%d / %d %s (brighter than %.1f)\n",
		 	    nbg, ng, title, gm[nbg-1]);
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
	if (wfile) {
	    if (objname)
		strcpy (filename,objname);
	    else
		strcpy (filename,"search");
	    if (refcat == GSC)
		strcat (filename,".gsc");
	    else if (refcat == UJC)
		strcat (filename,".ujc");
	    else if (refcat == USAC)
		strcat (filename,".usac");
	    else if (refcat == UAC)
		strcat (filename,".uac");
	    else
		strcat (filename,".");
	    strcat (filename,refcatname);

	    fd = fopen (filename, "w");

	    /* Free result arrays and return if cannot write file */
	    if (fd == NULL) {
		fprintf (stderr, "SCAT:  cannot write file %s\n", filename);
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

    /* Set degree flag for output */
    if (sysout == WCS_ECLIPTIC || sysout == WCS_GALACTIC)
	degout = 1;
    else
	degout = degout0;

    /* Write heading */
    if (refcat == GSC)
	sprintf (headline, "CATALOG	HSTGSC1.1");
    else if (refcat == USAC)
	sprintf (headline, "CATALOG	USNO SA 1.0");
    else if (refcat == UAC)
	sprintf (headline, "CATALOG	USNO A 1.0");
    else if (refcat == UJC)
	sprintf (headline, "CATALOG	USNO UJ1.0");
    else
	sprintf (headline, "CATALOG	%s", refcatname);
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (cra >= 0.0) {
	ra2str (rastr, 32, cra, 3);
	if (wfile)
	    fprintf (fd, "RA	%s\n", rastr);
	if (tabout)
	    printf ("RA	%s\n", rastr);
	}

    if (cdec >= -90.0) {
	dec2str (decstr, 32, cdec, 2);
	if (wfile)
	    fprintf (fd, "DEC	%s\n", decstr);
	if (tabout)
	    printf ("DEC	%s\n", decstr);
	}

    if (ddec > 0.0) {
	dec2str (decstr, 32, ddec, 2);
	if (wfile)
	    fprintf (fd, "RADIUS	%s\n", decstr);
	if (tabout)
	    printf ("RADIUS	%s\n", decstr);
	}

    if (syscoor == WCS_B1950)
	sprintf (headline, "EQUINOX	1950.0");
    else
	sprintf (headline, "EQUINOX	2000.0");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

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
    if (refcat == GSC)
	strcpy (headline, "gsc_id  	");
    else if (refcat == USAC)
	strcpy (headline,"usac_id  	");
    else if (refcat == UAC)
	strcpy (headline,"usnoa_id  	");
    else if (refcat == UJC)
	strcpy (headline,"usnoj_id  	");
    else
	strcpy (headline,"id    	");
    if (sysout == WCS_GALACTIC)
	strcat (headline,"long.ecl   	lat.ecl  	");
    else if (sysout == WCS_ECLIPTIC)
	strcat (headline,"long.ecl   	lat.ecl  	");
    else if (sysout == WCS_B1950)
	strcat (headline,"ra1950      	dec1950  	");
    else
	strcat (headline,"ra      	dec      	");
    if (refcat == USAC || refcat == UAC)
	strcat (headline,"magb	magr	plate  	arcsec");
    else
	strcat (headline,"mag	plate	arcsec");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (refcat == UAC || refcat == USAC)
	sprintf(headline,"----------	--------	---------	----	-----	-----	-----	------");
    else
        sprintf (headline,"----------	------------	------------	------	----	------");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (printhead) {
	if (nbg == 0) {
	    if (refcat == GSC)
		printf ("No Guide Stars Found\n");
	    else if (refcat == USAC)
		printf ("No USNO A 1.0 Stars Found\n");
	    else if (refcat == USAC)
		printf ("No USNO SA 1.0 Stars Found\n");
	    else if (refcat == UJC)
		printf ("No UJ 1.0 Stars Found\n");
	    else
		printf ("No Stars Found\n");
	    }
	else {
	    if (refcat == GSC)
		printf ("GSC number ");
	    else if (refcat == USAC)
		printf ("USNO SA number ");
	    else if (refcat == UAC)
		printf ("USNO A number  ");
	    else if (refcat == UJC)
		printf (" UJ number    ");
	    else
		printf (" Number    ");
	    if (sysout == WCS_B1950) {
		if (degout) {
		    if (eqout == 1950.0)
			printf ("  RA1950   Dec1950  ");
		    else
			printf ("RAB%7.2f DecB%7.2f  ", eqout, eqout);
		    }
		else {
		    if (eqout == 1950.0)
			printf ("RAB1950      DecB1950     ");
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
			printf (" RA2000       Dec2000     ");
		    else
			printf ("RAJ%7.2f   DecJ%7.2f  ", eqout, eqout);
		    }
		}
	    if (refcat == UAC || refcat == USAC)
		printf ("MagB  MagR Plate  Arcsec\n");
	    else if (refcat == GSC)
		printf (" Mag  Type  Arcsec\n");
	    else
		printf (" Mag  Plate Arcsec\n");
	    }
	}

    if (keyword != NULL)
	ntab = tabopen (refcat);

    for (i = 0; i < nbg; i++) {
	if (gx[i] > 0.0 && gy[i] > 0.0) {
	    wcscon (sysref, sysout, eqref, eqout, &gra[i],&gdec[i], epref);
	    if (degout) {
		deg2str (rastr, 32, gra[i], 5);
		deg2str (decstr, 32, gdec[i], 5);
		}
	    else {
		ra2str (rastr, 32, gra[i], 3);
		dec2str (decstr, 32, gdec[i], 2);
		}
	    if (refcat == UAC || refcat == USAC)
		sprintf (headline, "%13.8f	%s	%s	%.1f	%.1f	%d	%.2f",
		 gnum[i], rastr, decstr, gmb[i], gm[i], gc[i], 3600.0*gx[i]);
	    else if (refcat == UJC)
	        sprintf (headline, "%12.7f	%s	%s	%.2f	%d	%.2f",
		 gnum[i], rastr, decstr, gm[i], gc[i], 3600.0*gx[i]);
	    else {
	        sprintf (headline, "%9.4f	%s	%s	%.2f	%d	%.2f",
		 gnum[i], rastr, decstr, gm[i], gc[i], 3600.0*gx[i]);
		if (keyword != NULL) {
		    tabline = tabstar ((int)gnum[i], tabline);
		    (void) tabgetk (tabline, keyword, string, TABMAX);
		    strcat (headline, "	");
		    strcat (headline, string);
		    }
		}
	    if (wfile)
		fprintf (fd, "%s\n", headline);
	    else if (tabout)
		printf ("%s\n", headline);
	    else if (refcat == UAC || refcat == USAC)
		printf ("%13.8f %s %s %5.1f %5.1f %4d %8.2f\n",
			gnum[i],rastr,decstr,gmb[i],gm[i],gc[i], 3600.0*gx[i]);
	    else if (refcat == UJC)
		printf ("%12.7f %s %s %6.2f %4d %7.2f\n",
			gnum[i], rastr, decstr, gm[i],gc[i], 3600.0*gx[i]);
	    else if (refcat == GSC)
		printf ("%9.4f %s %s %6.2f %2d %7.2f\n",
			gnum[i], rastr, decstr, gm[i],gc[i], 3600.0*gx[i]);
	    else {
		if (keyword != NULL) {
		    tabline = tabstar ((int)gnum[i], tabline);
		    (void) tabgetk (tabline, keyword, string, TABMAX);
		    strcat (headline, "	");
		    strcat (headline, string);
		    printf ("%9.4f %s %s %6.2f %7d %7.2f %s\n",
		     gnum[i], rastr, decstr, gm[i],gc[i], 3600.0*gx[i], string);
		    }
		else
		    printf ("%9.4f %s %s %6.2f %7d %7.2f\n",
			gnum[i], rastr, decstr, gm[i],gc[i], 3600.0*gx[i]);
		}
	    }
	}

    if (keyword != NULL)
	tabclose();
    if (wfile)
	fclose (fd);
    if (gx) free ((char *)gx);
    if (gy) free ((char *)gy);
    if (gm) free ((char *)gm);
    if (gmb) free ((char *)gmb);
    if (gra) free ((char *)gra);
    if (gdec) free ((char *)gdec);
    if (gnum) free ((char *)gnum);
    if (gc) free ((char *)gc);

    return (nbg);
}


/* Get a center and radius for a search area.  If the image center is not
 * given in the system of the reference catalog, convert it.
 * Return 0 if OK, else -1
 */

static int
GetArea (verbose, syscoor, sysref, epoch, cra, cdec, dra, ddec, drad)

int	verbose;	/* Extra printing if =1 */
int	syscoor;	/* Coordinate system of search area */
int	sysref;		/* Coordinate system of reference catalog */
double	epoch;		/* Epoch in years of search area (ignored if zero) */
double	*cra;		/* Center right ascension in degrees (returned) */
double	*cdec;		/* Center declination in degrees (returned) */
double	*dra;		/* Right ascension half-width in degrees (returned) */
double	*ddec;		/* Declination half-width in degrees (returned) */
double	*drad;		/* Radius to search in degrees (0=box) (returned) */
{
    char rstr[32], dstr[32];

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
	if (syscoor == WCS_B1950)
	    fprintf (stderr,"Center:  %s   %s (B1950)\n", rstr, dstr);
	else if (syscoor == WCS_GALACTIC)
	    fprintf (stderr,"Center:  %s   %s (galactic)\n", rstr, dstr);
	else if (syscoor == WCS_ECLIPTIC)
	    fprintf (stderr,"Center:  %s   %s (ecliptic)\n", rstr, dstr);
	else
	    fprintf (stderr,"Center:  %s   %s (J2000)\n", rstr, dstr);
	}

    /* Convert coordinate system to match reference catalog system */
    wcscon (syscoor, sysref, eqcoor, eqref, cra, cdec, epoch);

    /* Set search box radius from command line, if it is there */
    if (rad0 > 0.0) {
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


    if (verbose && sysref != syscoor) {
	if (degout0) {
	    deg2str (rstr, 32, *cra, 5);
            deg2str (dstr, 32, *cdec, 5);
	    }
	else {
	    ra2str (rstr, 32, *cra, 3);
            dec2str (dstr, 32, *cdec, 2);
	    }
	if (sysref == WCS_B1950)
	    fprintf (stderr,"Center:  %s   %s (B1950)\n", rstr, dstr);
	else
	    fprintf (stderr,"Center:  %s   %s (J2000)\n", rstr, dstr);
	}
    if (verbose) {
	ra2str (rstr, 32, *dra * 2.0, 2); 
	dec2str (dstr, 32, *ddec * 2.0, 2); 
	fprintf (stderr,"Area:    %s x %s\n", rstr, dstr);
	}

    return (0);
}


static void
SearchHead (sys1, sys2, eq1, eq2, cra, cdec, dra, ddec, drad)

int	sys1;
int	sys2;
double	eq1, eq2;
double	cra, cdec;
double	dra, ddec;
double	drad;
{
    double ra, dec;
    char rastr[32];
    char decstr[32];

    ra = cra;
    dec = cdec;
    wcscon (sys1, sys2, eq1, eq2, &ra, &dec, epoch);
    if (sys2 == WCS_ECLIPTIC || sys2 == WCS_GALACTIC || degout0) {
	deg2str (rastr, 32, ra, 5);
	deg2str (decstr, 32, dec, 5);
	}
    else {
	ra2str (rastr, 32, ra, 3);
	dec2str (decstr, 32, dec, 2);
	}
    if (objname)
	printf ("%12s %s %s ", objname, rastr, decstr);
    else {
	if (refcat == GSC)
	    printf ("HST GSC   %s %s ", rastr, decstr);
	else if (refcat == USAC)
	    printf ("USNO SA 1.0   %s %s ", rastr, decstr);
	else if (refcat == UAC)
	    printf ("USNO A 1.0    %s %s ", rastr, decstr);
	else if (refcat == UJC)
	    printf ("USNO UJ1.0   %s %s ", rastr, decstr);
	else
	    printf ("%9s %s %s ", refcatname, rastr, decstr);
	}
    if (sys2 == WCS_B1950) {
	if (eq2 == 1950.0)
	    printf ("(B1950) ");
	else
	    printf ("(B%7.2f) ", eq2);
	}
    else if (sys2 == WCS_GALACTIC)
	printf ("(galactic) ");
    else if (sys2 == WCS_ECLIPTIC)
	printf ("(ecliptic) ");
    else {
	if (eq2 == 2000.0)
	    printf ("(J2000) ");
	else
	    printf ("(J%7.2f) ", eq2);
	}
    if (drad != 0.0)
	printf ("r= %.2f\n", drad*3600.0);
    else
	printf ("+- %.2f\n", ddec*3600.0);
    return;
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
 * Mar 12 1997	Add USNO SA 1.0 catalog as USAC
 * Apr 25 1997	Fix bug in uacread
 * May 29 1997	Add option to add keyword to tab table output
 * Nov 12 1997	Fix DEC in header to print Dec string instead of RA string
 * Nov 17 1997	Initialize both magnitude limits
 * Dec  8 1997	Set up program to be called by various names
 * Dec 12 1997	Fix center coordinate printing in heading
 *
 * Apr 10 1998	Fix bug search USNO A 1.0 catalog created by last revision
 * Apr 10 1998	Set search radius only if argument negative
 * Apr 14 1998	Version 2.2: to match other software
 * Apr 23 1998	Add output in galactic or ecliptic coordinates; use wcscon()
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * Jun  2 1998	Fix bug in tabread()
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Jun 30 1998	Fix declaration of GetArea()
 * Jul  1 1998	Allow specification of center different from output system
 * Jul  9 1998	Adjust all report headings
 */
