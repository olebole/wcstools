/* File scat.c
 * November 30, 1998
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
static int sysout = 0;		/* Output coordinate system */
static double eqcoor = 2000.0;	/* Equinox of search center */
static double eqout = 2000.0;	/* Equinox for output coordinates */
static int degout0 = 0;		/* 1 if degrees output instead of hms */
static double ra0 = -99.0;	/* Initial center RA in degrees */
static double dec0 = -99.0;	/* Initial center Dec in degrees */
static double rad0 = 10.0;	/* Search box radius */
static double epoch0 = 0.0;	/* Epoch for coordinates */
static int syscoor = 0;		/* Input search coordinate system */
static int nstars = 0;		/* Number of brightest stars to list */
static int printhead = 0;	/* 1 to print table heading */
static int tabout = 0;		/* 1 for tab table to standard output */
static int rasort = 0;		/* 1 to sort stars by brighness */
static int distsort = 0;	/* 1 to sort stars by distance from center */
static int debug = 0;		/* True for extra information */
static char *objname = NULL;	/* Object name for output */
static char *keyword = NULL;	/* Column to add to tab table output */
static char *progname;		/* Name of program as executed */
static char progpath[128];
static int ncat = 0;
static int version = 0;		/* If 1, print only program name and version */

main (ac, av)
int ac;
char **av;
{
    char *str, *str1;
    char rastr[16];
    char decstr[16];
    int nfind = 0;
    int i, nc, lcat;
    double *snum;
    char *refcatname[5];	/* reference catalog names */
    char *refcatn;
    snum = NULL;

    /* Check name used to execute programe and set catalog name accordingly */
    strcpy (progpath, av[0]);
    progname = progpath;
    for (i = strlen (progpath); i > -1; i--) {
	if (progpath[i] > 63 && progpath[i] < 90)
	    progpath[i] = progpath[i] + 32;
	if (progpath[i] == '/') {
	    progname = progpath + i;
	    break;
	    }
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

	/* Set star number if getting just one star */
	else if (isnum (*av)) {
	    if (snum == NULL)
		snum = (double *) calloc (1, sizeof(double));
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
    		rad0 = 120.0;
		break;

    	    case 'b':	/* initial coordinates on command line in B1950 */
		str1 = *(av+1);
		nc = (int)str1[0];
		if (syscoor || *(str+1) || nc < 47 || nc > 58) {
		    sysout = WCS_B1950;
		    eqout = 1950.0;
		    }
		else if (ac < 3)
		    usage (progname);
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

	    case 'e':	/* Set ecliptic coordinate output and optional center */
		str1 = *(av+1);
		nc = (int)str1[0];
		if (syscoor || *(str+1) || nc < 47 || nc > 58)
		    sysout = WCS_ECLIPTIC;
		else if (ac < 3)
		    usage (progname);
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
		    usage (progname);
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
		if (syscoor || *(str+1) || nc < 47 || nc > 58) {
		    sysout = WCS_J2000;
		    eqout = 2000.0;
		    }
		else if (ac < 3)
		    usage (progname);
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
		    usage (progname);
		keyword = *++av;
		ac--;
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
    		eqout = atof (*++av);
    		ac--;
    		break;

    	    case 'r':	/* Box radius in arcseconds */
    		if (ac < 2)
    		    usage(progname);
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
		epoch0 = atof (*++av);
		ac--;
		break;

	    default:
		usage (progname);
		break;
	    }
	}
    }

    ListCat (nfind, snum, ncat, refcatname);

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
    else
	fprintf (stderr,"Find catalog stars in a square on the sky\n");
    fprintf (stderr,"Usage: [-adhstvw] [-m [mag1] mag2] [-e sys] [-n num] [-r arcsec] [-b][-j] ra dec\n");
    fprintf(stderr,"  -a: List single closest catalog source\n");
    fprintf(stderr,"  -b: Output B1950 (FK4) coordinates with optional center\n");
    fprintf(stderr,"  -c: Reference catalog (gsc, uac, usac, or tab table file\n");
    fprintf(stderr,"  -d: Output RA and Dec in degrees instead of hms dms\n");
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
ListCat (nfind, snum, ncat, refcatname)

int	nfind;		/* Number of stars to find */
double	*snum;		/* Catalog numbers */
char	**refcatname;	/* reference catalog name */

{
    double *gnum=0;	/* Catalog star numbers */
    double *gra=0;	/* Catalog star right ascensions, rads */
    double *gdec=0;	/* Catalog star declinations rads */
    double *gm=0;	/* Catalog magnitudes */
    double *gmb=0;	/* Catalog B magnitudes */
    double *gx=0;	/* Catalog star X positions on image */
    double *gy=0;	/* Catalog star Y positions on image */
    int *gc=0;		/* Catalog star object classes */
    char **gobj=0;	/* Catalog star object names */
    char **gobj1=NULL;	/* Catalog star object names */
    double cra, cdec;
    double epout;
    int sysref;		/* Coordinate system of reference catalog */
    double eqref;	/* Equinox of catalog to be searched */
    double epref;	/* Epoch of catalog to be searched */
    int ng;		/* Number of catalog stars */
    int nbg;		/* Number of brightest catalog stars actually used */
    int i, ngmax, nbytes;
    int degout;
    FILE *fd;
    char rastr[32], decstr[32];	/* coordinate strings */
    double drad, dra, ddec, ra1, dec1, mag1, mag2;
    double gdist;
    double mag;
    int offscale, nlog, closest;
    char headline[160];
    char filename[80];
    char title[80];
    char string[TABMAX];
    char temp[80];
    char isp[4];
    int ntab;
    int refcat;		/* reference catalog switch */
    int icat;

    if (verbose || printhead) {
	if (closest)
	else
	}

    /* Start of per ctalog loop */
    for (icat = 0; icat < ncat; icat++) {
    isp[2] = (char) 0;
    isp[3] = (char) 0;
    if (distsort && nstars == 1)
	closest = 1;
    else
	closest = 0;
    if (!sysout) {
	sysout = syscoor;
	if (syscoor == WCS_B1950)
	    eqout = 1950.0;
	else
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
    epout = epoch0;
    if (epout == 0.0)
	epout = epref;

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
	if (!(gobj = malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for obj\n", nbytes);
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
	    nbg = binrnum ("SAO",nfind,sysout,eqout,epout,
			   gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == PPM)
	    nbg = binrnum ("PPM",nfind,sysout,eqout,epout,
			   gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == IRAS)
	    nbg = binrnum ("IRAS",nfind,sysout,eqout,epout,
			   gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == TYCHO)
	    nbg = binrnum ("tycho",nfind,sysout,eqout,epout,
			   gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == TABCAT) {
	    nbg = tabrnum (refcatname[icat],nfind,sysout,eqout,epout,
			   gnum,gra,gdec,gm,gc,debug);
	    if (keyword != NULL)
		tabrkey (refcatname[icat], nfind, gnum, keyword, gobj);
	    }
	else if (refcat == BINCAT)
	    nbg = binrnum (refcatname[icat],nfind,sysout,eqout,epout,
			   gnum,gra,gdec,gm,gmb,gc,gobj,nlog);
	else
	    nbg = catrnum (refcatname[icat],nfind,sysout,eqout,epout,
			   gnum,gra,gdec,gm,gobj,debug);

	for (i = 0; i < nbg; i++ ) {
	    gx[i] = 0.0;
	    gy[i] = 1.0;
	    }
	}

    /* Find stars specified by location */
    else {

	/* Set limits from defaults and command line information */
	if (GetArea (verbose,syscoor,sysout,eqout,epout,
		     &cra,&cdec,&dra,&ddec,&drad))
	    return (0);

	/* Print search center and size in input and output coordinates */
	if ((verbose || printhead) && !closest) {
	    SearchHead (refcatname[icat],sysout,sysout,eqout,eqout,epout,epout,
			cra,cdec,dra,ddec,drad);
	    if (sysout != syscoor)
		SearchHead (refcatname[icat],sysout,syscoor,eqout,eqcoor,
			    epout,epout,cra,cdec,dra,ddec,drad);
	    if (sysref != syscoor && sysref != sysout)
		SearchHead (refcatname[icat],sysout,sysref,eqout,eqref,
			    epout,epref,cra,cdec,dra,ddec,drad);
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
	if (icat == 0) {
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
	if (!(gx = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gx\n", nbytes);
	if (!(gy = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gy\n", nbytes);
	if (!(gobj = malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for obj\n", nbytes);
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
	else if (refcat == BINCAT)
	    ng = binread (refcatname[icat],cra,cdec,dra,ddec,drad,sysout,eqout,epout,
			  mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,gobj,nlog);
	else if (refcat == TABCAT) {
	    ng = tabread (refcatname[icat],cra,cdec,dra,ddec,drad,sysout,eqout,epout,
			  mag1,mag2,ngmax,gnum,gra,gdec,gm,gc,nlog);
	    if (keyword != NULL)
		tabrkey (refcatname[icat], ng, gnum, keyword, gobj);
	    }
	else
	    ng = catread (refcatname[icat],cra,cdec,dra,ddec,drad,sysout,eqout,epout,
			  mag1,mag2,ngmax,gnum,gra,gdec,gm,gobj,nlog);
	if (gobj[0] == NULL)
	    gobj1 = NULL;
	else
	    gobj1 = gobj;

	/* Compute distance from star to search center */
	for (i = 0; i < ng; i++ ) {
	    offscale = 0;
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
    else
	sprintf (headline, "CATALOG	%s", refcatname[icat]);
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
    sprintf (headline, "EQUINOX	%.1f", eqout);
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
    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
	refcat == UAC  || refcat == UA1  || refcat == UA2)
	strcat (headline,"magb	magr	plate");
    else if (refcat == TYCHO)
	strcat (headline,"magb	magv	type");
    else
	strcat (headline,"mag	type");
    if (snum == NULL)
	strcat (headline,"	arcsec");
    if (refcat == TABCAT && keyword != NULL) {
	strcat (headline,"	");
	strcat (headline, keyword);
	}
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
	refcat == UAC  || refcat == UA1  || refcat == UA2 || refcat == TYCHO)
	sprintf(headline,"----------	--------	---------	----	-----	-----");
    else
        sprintf (headline,"----------	------------	------------	------	----");
    if (snum == NULL)
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
	    else
		printf ("Number    ");
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
	    else if (refcat == GSC)
		printf (" Mag  Type");
	    else if (refcat == SAO || refcat == PPM || refcat == IRAS)
		printf (" Mag  Type");
	    else if (refcat == TYCHO)
		printf (" MagB  MagV  Type");
	    else if (refcat == TABCAT)
		printf (" Mag     Peak");
	    else
		printf ("    Mag");
	    if (snum == NULL)
		printf ("  Arcsec\n");
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
	    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		refcat == UAC  || refcat == UA1  || refcat == UA2)
		sprintf (headline, "%13.8f	%s	%s	%.1f	%.1f	%d",
		 gnum[i], rastr, decstr, gmb[i], gm[i], gc[i]);
	    else if (refcat == UJC)
	        sprintf (headline, "%12.7f	%s	%s	%.2f	%d",
		 gnum[i], rastr, decstr, gm[i], gc[i]);
	    else if (refcat==SAO || refcat==PPM || refcat==IRAS ) {
		isp[0] = gc[i] / 1000;
		isp[1] = gc[i] % 1000;
	        sprintf (headline, "%9d	%s	%s	%.2f	%2s",
		 (int)(gnum[i]+0.5), rastr, decstr, gm[i], isp);
		}
	    else if (refcat == TYCHO) {
		isp[0] = gc[i] / 1000;
		isp[1] = gc[i] % 1000;
	        sprintf (headline, "%10.5f	%s	%s	%.2f	%.2f	%2s",
		 gnum[i], rastr, decstr, gmb[i], gm[i], isp);
		}
	    else if (refcat == TABCAT)
	        sprintf (headline, "%9.4f	%s	%s	%.2f	%d",
		 gnum[i], rastr, decstr, gm[i], gc[i]);
	    else
	        sprintf (headline, "%9d	%s	%s	%.2f	%d",
		 (int)gnum[i], rastr, decstr, gm[i], gc[i]);
	    if (snum == NULL) {
	        sprintf (temp, "	%.2f", gdist);
	        strcat (headline, temp);
		}
	    if (refcat == TABCAT && keyword != NULL) {
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
		    printf ("%13.8f %s %s %5.1f %5.1f %4d ",
			gnum[i],rastr,decstr,gmb[i],gm[i],gc[i]);
		else if (refcat == UJC)
		    printf ("%12.7f %s %s %6.2f %4d",
			gnum[i], rastr, decstr, gm[i],gc[i]);
		else if (refcat == GSC)
		    printf ("%9.4f %s %s %6.2f %2d",
			gnum[i], rastr, decstr, gm[i],gc[i]);
		else if (refcat==SAO || refcat==PPM || refcat==IRAS ) {
		    isp[0] = gc[i] / 1000;
		    isp[1] = gc[i] % 1000;
		    printf ("%8d  %s %s %6.2f  %2s",
			(int)(gnum[i]+0.5),rastr,decstr,gm[i],isp);
		    }
		else if (refcat == TYCHO) {
		    isp[0] = gc[i] / 1000;
		    isp[1] = gc[i] % 1000;
		    printf ("%10.5f %s %s %6.2f %6.2f  %2s",
			gnum[i],rastr,decstr,gmb[i],gm[i],isp);
		    }
		else if (refcat == TABCAT)
		    printf ("%9.4f %s %s %6.2f %7d",
			gnum[i], rastr, decstr, gm[i],gc[i]);
		else if (refcat == BINCAT) {
		    isp[0] = gc[i] / 1000;
		    isp[1] = gc[i] % 1000;
		    printf ("%9.4f %s %s %6.2f %2s",
			gnum[i], rastr, decstr, gm[i], isp);
		    }
		else
		    printf ("%9.4f %s %s %6.2f",
			gnum[i], rastr, decstr, gm[i]);
		if (snum == NULL)
		    printf ("  %7.2f", gdist);
		if (refcat == TABCAT && keyword != NULL)
		    printf (" %s\n", gobj[i]);
		else if (refcat == BINCAT && gobj != NULL && gobj[i] != NULL)
		    printf (" %s\n", gobj[i]);
		else if (refcat == TXTCAT && gobj != NULL && gobj[i] != NULL)
		    printf (" %s\n", gobj[i]);
		else
		    printf ("\n");
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

    if (verbose) {
	ra2str (rstr, 32, *dra * 2.0, 2); 
	dec2str (dstr, 32, *ddec * 2.0, 2); 
	fprintf (stderr,"Area:    %s x %s\n", rstr, dstr);
	}

    return (0);
}


static void
SearchHead (refcatname,sys1,sys2, eq1,eq2, ep1,ep2, cra, cdec, dra, ddec, drad)

char	*refcatname;
int	sys1;
int	sys2;
double	eq1, eq2;
double	ep1, ep2;
double	cra, cdec;
double	dra, ddec;
double	drad;
{
    double ra, dec;
    char rastr[32];
    char decstr[32];
    char cstr[16];
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
    if (objname)
	printf ("%12s %s %s ", objname, rastr, decstr);
    else {
        if (!(refcat = RefCat (refcatname, title, &sysref, &eqref, &epref))) {
	    fprintf (stderr,"ListCat: Catalog '%s' is missing\n",refcatname);
	    return;
	    }
	if (refcat == GSC)
	    printf ("HST GSC   %s %s ", rastr, decstr);
	else if (refcat == USAC)
	    printf ("USNO SA       %s %s ", rastr, decstr);
	else if (refcat == USA1)
	    printf ("USNO SA-1.0   %s %s ", rastr, decstr);
	else if (refcat == USA2)
	    printf ("USNO SA-2.0   %s %s ", rastr, decstr);
	else if (refcat == UAC)
	    printf ("USNO A        %s %s ", rastr, decstr);
	else if (refcat == UA1)
	    printf ("USNO A-1.0    %s %s ", rastr, decstr);
	else if (refcat == UA2)
	    printf ("USNO A-2.0    %s %s ", rastr, decstr);
	else if (refcat == UJC)
	    printf ("USNO UJ1.0   %s %s ", rastr, decstr);
	else if (refcat == SAO)
	    printf ("SAO       %s %s ", rastr, decstr);
	else if (refcat == PPM)
	    printf ("PPM       %s %s ", rastr, decstr);
	else if (refcat == IRAS)
	    printf ("IRAS      %s %s ", rastr, decstr);
	else if (refcat == TYCHO)
	    printf ("Tycho     %s %s ", rastr, decstr);
	else
	    printf ("%9s %s %s ", refcatname, rastr, decstr);
	}
    wcscstr (cstr, sys2, eq2, ep2);
    printf ("%s ", cstr);
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
 */
