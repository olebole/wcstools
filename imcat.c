/* File imcat.c
 * May 12, 1999
 * By Doug Mink
 * (Harvard-Smithsonian Center for Astrophysics)
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
#include "wcscat.h"
#include "libwcs/lwcs.h"

#define MAXREF 100

static void usage();
static void ListCat();
extern void fk524e();
extern struct WorldCoor *GetFITSWCS();
extern char *GetFITShead();
extern void setsys();
extern void setcenter();
extern void setsecpix();
extern void setrefpix();


static int verbose = 0;		/* verbose/debugging flag */
static int wfile = 0;		/* True to print output file */
static int refcat = GSC;	/* reference catalog switch */
static int classd = -1;		/* Guide Star Catalog object classes */
static int uplate = 0;		/* UJ Catalog plate number to use */
static double maglim1 = MAGLIM1; /* reference catalog bright magnitude limit */
static double maglim2 = MAGLIM2; /* reference catalog faint magnitude limit */
static int nstars = 0;		/* Number of brightest stars to list */
static int printhead = 0;	/* 1 to print table heading */
static int tabout = 0;		/* 1 for tab table to standard output */
static int rasort = 0;		/* 1 to sort stars by brighness */
static int debug = 0;		/* True for extra information */
static int degout0 = 0;		/* True for RA and Dec in fractional degrees */
static char *keyword = NULL;	/* Column to add to tab table output */
static int sysout = 0;		/* Output coordinate system */
static int sysref = WCS_J2000;	/* Output coordinate system */
static double eqout = 0.0;	/* Equinox for output coordinates */
static int version = 0;		/* If 1, print only program name and version */

main (ac, av)
int ac;
char **av;
{
    char *str, *str1;
    char rastr[16];
    char decstr[16];
    int readlist = 0;
    char *lastchar;
    char filename[128];
    FILE *flist;
    char *listfile;
    char *cstr;
    int i, ic;
    double x, y;
    char *refcatname[5];	/* reference catalog name */
    int ncat = 0;
    int region_radius[5];	/* Flag for SAOimage region file output */
    int rcat = 0;
    int region_char[5];		/* Character for SAOimage region file output */
    char *refcatn;
    int lcat;
    int scat = 0;
    char progpath[128];
    char *progname;		/* Name of program as executed */

    for (i = 0; i < 5; i++) {
	region_radius[i] = 0;
	region_char[i] = 0;
	}

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
    if (strsrch (progname,"gsc") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "gsc");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"uac") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "uac");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"ua1") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "ua1");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"ua2") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "ua2");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"usac") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "usac");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"usa1") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "usa1");
	refcatname[0] = refcatn;
	}
    else if (strsrch (progname,"usa2") != NULL) {
	ncat = 1;
	refcatn = (char *) calloc (1,8);
	strcpy (refcatn, "usa2");
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

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage (progname);
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage (progname);
	}

    /* Decode arguments */
    for (av++; --ac > 0 && (*(str = *av)=='-' || *str == '@'); av++) {
	char c;
	if (*str == '@')
	    str = str - 1;
	while (c = *++str)
    	switch (c) {

    	case 'b':	/* initial coordinates on command line in B1950 */
	    str1 = *(av+1);
	    ic = (int)str1[0];
	    if (*(str+1) || ic < 48 || ic > 58) {
		setsys(WCS_B1950);
		sysout = WCS_B1950;
		eqout = 1950.0;
		}
	    else if (ac < 3)
		usage (progname);
	    else {
		setsys(WCS_B1950);
		sysout = WCS_B1950;
		eqout = 1950.0;
		strcpy (rastr, *++av);
		ac--;
		strcpy (decstr, *++av);
		ac--;
		setcenter (rastr, decstr);
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

	case 'd':
	    degout0++;
	    break;

	case 'e':
	    sysout = WCS_ECLIPTIC;
	    break;

	case 'g':
	    sysout = WCS_GALACTIC;
	    break;

	case 'h':	/* ouput descriptive header */
	    printhead++;
	    break;

    	case 'j':	/* center coordinates on command line in J2000 */
	    str1 = *(av+1);
	    ic = (int)str1[0];
	    if (*(str+1) || ic < 48 || ic > 58) {
		setsys(WCS_J2000);
		sysout = WCS_J2000;
		eqout = 2000.0;
		}
	    else if (ac < 3)
		usage (progname);
	    else {
		setsys(WCS_J2000);
		sysout = WCS_J2000;
		eqout = 2000.0;
		strcpy (rastr, *++av);
		ac--;
		strcpy (decstr, *++av);
		ac--;
		setcenter (rastr, decstr);
		}
    	    break;

	case 'k':	/* Keyword (column) to add to output from tab table */
	    if (ac < 2)
		usage (progname);
	    keyword = *++av;
	    ac--;
	    break;

    	case 'm':	/* Limiting reference star magnitude */
    	    if (ac < 2)
    		usage(progname);
    	    maglim2 =  atof (*++av);
    	    ac--;
	    if (ac > 1 && isnum (*(av+1))) {
		maglim1 = maglim2;
		maglim2 = atof (*++av);
		ac--;
		}
	    else if (maglim1 == 0.0)
		maglim1 = -2.0;
    	    break;

	case 'n':	/* Number of brightest stars to read */
	    if (ac < 2)
		usage (progname);
	    nstars = atoi (*++av);
	    ac--;
	    break;

	case 'o':	/* Guide Star object class */
    	    if (ac < 2)
    		usage(progname);
    	    classd = (int) atof (*++av);
    	    ac--;
    	    break;

    	case 'p':	/* Initial plate scale in arcseconds per pixel */
    	    if (ac < 2)
    		usage(progname);
    	    setsecpix (atof (*++av));
    	    ac--;
    	    break;

	case 'q':	/* Output region file shape for SAOimage */
    	    if (ac < 2)
    		usage(progname);
	    cstr = *++av;
	    switch (cstr[0]){
		case 'c':
		    if (cstr[1] == 'i')
			region_char[scat] = WCS_CIRCLE;
		    else
			region_char[scat] = WCS_CROSS;
		    break;
		case 'd':
		    region_char[scat] = WCS_DIAMOND;
		    break;
		case 's':
		    region_char[scat] = WCS_SQUARE;
		    break;
		case 'x':
		    region_char[scat] = WCS_EX;
		    break;
		case 'v':
		    region_char[scat] = WCS_VAR;
		    break;
		case '+':
		    region_char[scat] = WCS_CROSS;
		    break;

		case 'o':
		default:
		    region_char[scat] = WCS_CIRCLE;
		}
	    if (region_radius[scat] == 0)
		region_radius[scat] = -1;
	    scat++;
    	    wfile++;
    	    ac--;
	    break;

	case 'r':	/* Output region file with shape radius for SAOimage */
    	    if (ac < 2)
    		usage(progname);
	    region_radius[rcat] = atoi (*++av);
	    if (region_radius[rcat] == 0)
		region_radius[rcat] = -1;
	    rcat++;
    	    wfile++;
    	    ac--;
	    break;

	case 's':	/* sort by RA */
	    rasort = 1;
	    break;

	case 't':	/* tab table to stdout */
	    tabout = 1;
	    break;

	case 'u':	/* UJ Catalog plate number */
    	    if (ac < 2)
    		usage(progname);
    	    uplate = (int) atof (*++av);
    	    ac--;
    	    break;

    	case 'v':	/* more verbosity */
    	    verbose++;
    	    break;

    	case 'w':	/* write output file */
    	    wfile++;
    	    break;

	case 'x':	/* X and Y coordinates of reference pixel */
	    if (ac < 3)
		usage(progname);
	    x = atof (*++av);
	    ac--;
	    y = atof (*++av);
	    ac--;
    	    setrefpix (x, y);
    	    break;

	case 'z':       /* Use AIPS classic WCS */
	    setdefwcs (1);
	    break;

	case '@':       /* List of files to be read */
	    readlist++;
	    listfile = ++str;
	    str = str + strlen (str) - 1;
	    av++;
	    ac--;


    	default:
    	    usage(progname);
    	    break;
    	}
    }
    /* if (!verbose && !wfile)
	verbose = 1; */

    /* Find number of images to search and leave listfile open for reading */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMCAT: List file %s cannot be read\n",
		     listfile);
	    usage (progname);
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
    	    if (debug)
    		printf ("%s:\n", filename);
	    ListCat (progname,filename,ncat,refcatname,region_radius,region_char);
	    }
	fclose (flist);
	}

    /* If no arguments left, print usage */
    /* if (ac == 0 || !refcat) */
    if (ac == 0)
	usage(progname);

    while (ac-- > 0) {
    	char *fn = *av++;
    	if (debug)
    	    printf ("%s:\n", fn);
    	ListCat (progname,fn, ncat, refcatname, region_radius, region_char);
    	if (debug)
    	    printf ("\n");
	}

    return (0);
}

static void
usage (progname)
char	*progname;
{
    if (version)
	exit (1);
    if (strsrch (progname,"gsc") != NULL) {
	fprintf (stderr,"List HST Guide Stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-g class]\n");
	}
    else if (strsrch (progname,"ujc") != NULL) {
	fprintf (stderr,"List USNO J Catalog stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-u plate]\n");
	}
    else if (strsrch (progname,"uac") != NULL) {
	fprintf (stderr,"List USNO A stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-u plate]\n");
	}
    else if (strsrch (progname,"ua1") != NULL) {
	fprintf (stderr,"List USNO A-1.0 stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-u plate]\n");
	}
    else if (strsrch (progname,"ua2") != NULL) {
	fprintf (stderr,"List USNO A-2.0 stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-u plate]\n");
	}
    else if (strsrch (progname,"usac") != NULL) {
	fprintf (stderr,"List USNO SA stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-u plate]\n");
	}
    else if (strsrch (progname,"usa1") != NULL) {
	fprintf (stderr,"List USNO SA-1.0 stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-u plate]\n");
	}
    else if (strsrch (progname,"usa2") != NULL) {
	fprintf (stderr,"List USNO SA-2.0 stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-u plate]\n");
	}
    else if (strsrch (progname,"act") != NULL) {
	fprintf (stderr,"List ACT Catalog Stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2]\n");
	}
    else if (strsrch (progname,"iras") != NULL) {
	fprintf (stderr,"List IRAS Point Sources in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2]\n");
	}
    else if (strsrch (progname,"sao") != NULL) {
	fprintf (stderr,"List SAO Catalog stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2]\n");
	}
    else if (strsrch (progname,"ppm") != NULL) {
	fprintf (stderr,"List PPM Catalog stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2]\n");
	}
    else if (strsrch (progname,"tycho") != NULL) {
	fprintf (stderr,"List Tycho Catalog stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2]\n");
	}
    else {
	fprintf (stderr,"List catalog stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vwhst][-m [mag1] mag2][-c or f catalog][-x x y]\n");
	}
    fprintf (stderr,"       [-p scale][-q osd+x][-b ra dec][-j ra dec][-r arcsec] FITS or IRAF file(s)\n");
    fprintf (stderr,"  -b: Output, (center) in B1950 (FK4) RA and Dec\n");
    fprintf (stderr,"  -c: Reference catalog (gsc, uac, or tab table or ASCII file\n");
    fprintf (stderr,"  -d: Output RA,Dec positions in fractional degrees\n");
    fprintf (stderr,"  -e: Output in ecliptic longitude and latitude\n");
    fprintf (stderr,"  -f: Reference catalog (gsc, uac, or binary table file\n");
    fprintf (stderr,"  -g: Output in galactic longitude and latitude\n");
    fprintf (stderr,"  -h: Print heading, else do not \n");
    fprintf (stderr,"  -j: Output (center) in J2000 (FK5) RA and Dec\n");
    fprintf (stderr,"  -k: Add this keyword to output from tab table search\n");
    fprintf (stderr,"  -m: Limiting catalog magnitude(s) (default none)\n");
    fprintf (stderr,"  -n: Number of brightest stars to print \n");
    fprintf (stderr,"  -o: Set HST Guide Star object class to print \n");
    fprintf (stderr,"  -p: Initial plate scale in arcsec per pixel (default 0)\n");
    fprintf (stderr,"  -q: Write SAOimage region file of this shape (filename.cat)\n");
    fprintf (stderr,"  -r: Write SAOimage region file of this radius (filename.cat)\n");
    fprintf (stderr,"  -s: Sort by RA instead of flux \n");
    fprintf (stderr,"  -t: Tab table to standard output as well as file\n");
    fprintf (stderr,"  -u: USNO catalog single plate number to accept\n");
    fprintf (stderr,"  -w: Write tab table output file [imagename].[catalog]\n");
    fprintf (stderr,"  -v: Verbose\n");
    fprintf (stderr,"  -x: X and Y coordinates of reference pixel (default is center)\n");
    fprintf (stderr,"  -z: Use AIPS classic projections instead of WCSLIB\n");
    exit (1);
}


struct WorldCoor *wcsinit();	

static void
ListCat (progname, filename, ncat, refcatname, region_radius, region_char)

char	*progname;	/* Name of program being executed */
char	*filename;	/* FITS or IRAF file filename */
int	ncat;		/* Number oc catalogs to search */
char	**refcatname;	/* reference catalog name */
int	*region_radius;	/* Flag for SAOimage region file output */
int	*region_char;	/* Character for SAOimage region file output */

{
    char *header;	/* FITS image header */
    double *gnum=0;	/* Catalog numbers */
    double *gra=0;	/* Catalog right ascensions, rads */
    double *gdec=0;	/* Catalog declinations rads */
    double *gm=0;	/* Catalog star magnitudes */
    double *gmb=0;	/* Catalog star blue magnitudes */
    double *gx, *gy;	/* Catalog star positions on image */
    char **gobj=NULL;	/* Catalog object names */
    char **gobj1=NULL;	/* Catalog object names */
    int *gc;		/* Catalog object classes, plates, etc. */
    int ng;		/* Number of catalog stars */
    int nbg;		/* Number of brightest catalog stars actually used */
    int imh, imw;	/* Image height and width in pixels */
    int i, ngmax, nbytes;
    int degout;
    FILE *fd;
    struct WorldCoor *wcs;	/* World coordinate system structure */
    double eqref;	/* Equinox of catalog to be searched */
    double epref;	/* Epoch of catalog to be searched */
    double epout;	/* Epoch of catalog to be searched */
    int sysref;		/* Coordinate system of catalog to be searched */
    char rastr[16], decstr[16];	/* coordinate strings */
    double cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, mag1, mag2,secpix;
    double mag, gdist, drad;
    int offscale, nlog;
    char headline[160];
    char temp[80];
    char title[80];
    char outfile[80];
    char *fname;
    char isp[4];
    int icat;
    int printobj = 0;
    int nndec;
    char nform[64];
    struct StarCat *starcat;

    isp[2] = 0;
    isp[3] = 0;
    if (verbose || printhead)

    /* Loop through catalogs */
    for (icat = 0; icat < ncat; icat++) {

    /* Find title and coordinate system for catalog */
    if (!(refcat = RefCat (refcatname[icat],title,&sysref,&eqref,&epref))) {
	fprintf (stderr,"ListCat: No catalog named %s\n", refcatname);
	return;
	}
    if (classd == 0)
	strcat (title, " stars");
    else if (classd == 3)
	strcat (title, " nonstars");

    /* Read world coordinate system information from the image header */
    if ((header = GetFITShead (filename)) == NULL)
	return;
    wcs = GetFITSWCS (filename, header, verbose, &cra, &cdec, &dra, &ddec,
		      &secpix, &imw, &imh, &sysout, &eqout);
    free (header);
    if (nowcs (wcs))
	return;

    /* Set up limits for search */
    ra1 = cra - dra;
    ra2 = cra + dra;
    dec1 = cdec - ddec;
    dec2 = cdec + ddec;
    epout = wcs->epoch;
    if (verbose || printhead) {
	char rastr1[16],rastr2[16],decstr1[16],decstr2[16], cstr[16];
	wcscstr (cstr, sysout, eqout, epout);
	ra2str (rastr1, 16, ra1, 3);
	ra2str (rastr2, 16, ra2, 3);
	printf ("%s: RA:  %s - %s %s\n",filename,rastr1, rastr2, cstr);
	dec2str (decstr1, 16, dec1, 2);
	dec2str (decstr2, 16, dec2, 2);
	printf ("%s: Dec: %s - %s %s\n",filename, decstr1, decstr2, cstr);
	}

/* Set the magnitude limits for the search */
    if (maglim2 == 0.0) {
	mag1 = 0.0;
	mag2 = 0.0;
	}
    else {
	mag1 = maglim1;
	mag2 = maglim2;
	}
    if (mag2 < mag1) {
	mag = mag1;
	mag1 = mag2;
	mag2 = mag;
	}

    if (nstars > 0)
	ngmax = nstars;
    else
	ngmax = MAXREF;
    nbytes = MAXREF * sizeof (double);

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
	free ((char *)wcs);
	return;
	}
    for (i = 0; i < ngmax; i++)
	gmb[i] = 0.0;
    if (debug)
	nlog = 100;
    else
	nlog = 0;

    /* Find the nearby reference stars, in ra/dec */
    drad = 0.0;
    if (refcat == GSC)
	ng = gscread (cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,mag2,
		      classd,ngmax,gnum,gra,gdec,gm,gc,nlog);
    else if (refcat == UAC  || refcat == UA1  || refcat == UA2 ||
	     refcat == USAC || refcat == USA1 || refcat == USA2)
	ng = uacread (refcatname[icat],cra,cdec,dra,ddec,drad,sysout,eqout,
		     epout,mag1,mag2,uplate,ngmax,gnum,gra,gdec,gm,gmb,gc,nlog);
    else if (refcat == UJC)
	ng = ujcread (cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,mag2,
		      uplate,ngmax,gnum,gra,gdec,gm,gc,nlog);
    else if (refcat == ACT)
	ng = actread (cra,cdec,dra,ddec,drad,sysout,eqout,epout,
		      mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,nlog);
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
	ng = binread ("hipparcosra",cra,cdec,dra,ddec,drad,sysout,eqout,epout,
		      mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,NULL,nlog);
    else if (refcat == BINCAT) {
	ng = binread (refcatname[icat],cra,cdec,dra,ddec,drad,sysout,eqout,
		      epout,mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,gobj,nlog);
	starcat = binopen (refcatname[icat]);
	nndec = starcat->nndec;
	binclose (starcat);
	}
    else if (refcat == TABCAT) {
	ng = tabread (refcatname[icat],cra,cdec,dra,ddec,drad,sysout,eqout,
		      epout,mag1,mag2,ngmax,gnum,gra,gdec,gm,gc,nlog);
	if (keyword != NULL)
	    tabrkey (refcatname[icat], ng, gnum, keyword, gobj);
	}
    else {
	ng = catread (refcatname[icat],cra,cdec,dra,ddec,drad,sysout,eqout,
		      epout,mag1,mag2,ngmax,gnum,gra,gdec,gm,gobj,nlog);
	starcat = catopen (refcatname[icat]);
	nndec = starcat->nndec;
	catclose (starcat);
	}
    if (gobj[0] == NULL)
	gobj1 = NULL;
    else
	gobj1 = gobj;

    /* Get image pixel coordinates for each star found in reference catalog */
    for (i = 0; i < ng; i++ ) {
	offscale = 0;
	wcs2pix (wcs, gra[i], gdec[i], &gx[i], &gy[i], &offscale);
	}

    /* Sort reference stars by brightness (magnitude) */
    MagSortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, gobj1, ng);

    /* List the brightest MAXSTARS reference stars */
    if (nstars > 0 && ng > nstars) {
	nbg = nstars;
	if (verbose || printhead)
	    if (mag2 > 0.0)
		printf ("%d / %d %s between %.1f and %.1f",
			nbg, ng, title, gm[0], gm[nbg-1]);
	    else
		printf ("%d / %d %s brighter than %.1f",
			nbg, ng, title, gm[nbg-1]);
	}
    else {
	nbg = ng;
	if (verbose || printhead) {
	    if (maglim1 > 0.0)
		printf ("%d %s between %.1f and %.1f",ng,title,maglim1,maglim2);
	    else if (maglim2 > 0.0)
		printf ("%d %s brighter than %.1f",ng, title, maglim2);
	    else
		printf ("%d %s", ng, title);
	    }
	}
    if (printhead || verbose) {
	if (wcs->epoch != wcs->equinox)
	    printf (" at %7.2f",wcs->epoch);
	if (isiraf (filename))
	    printf (" in IRAF image %s\n",filename);
	else
	    printf (" in FITS image %s\n", filename);
	}

    /* Sort star-like objects in image by right ascension */
    if (rasort)
	RASortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, gobj1, nbg);

    sprintf (headline, "IMAGE	%s", filename);

    /* Open plate catalog file */
    if (wfile && icat == 0) {
	fname = strrchr (filename, '/');
	if (fname != NULL)
	    strcpy (outfile, fname+1);
	else
	    strcpy (outfile, filename);
	for (i = 0; i < ncat; i++) {
	    strcat (outfile,".");
	    strcat (outfile,refcatname[i]);
	    }
	fd = fopen (outfile, "w");
	if (fd == NULL) {
	    fprintf (stderr, "IMCAT:  cannot write file %s\n", outfile);
	    if (gx) free ((char *)gx);
	    if (gy) free ((char *)gy);
	    if (gm) free ((char *)gm);
	    if (gmb) free ((char *)gmb);
	    if (gra) free ((char *)gra);
	    if (gdec) free ((char *)gdec);
	    if (gnum) free ((char *)gnum);
	    if (gc) free ((char *)gc);
	    if (gx) free ((char *)gx);
	    if (gy) free ((char *)gy);
	    if (gobj) free ((char *)gobj);
	    free ((char *)wcs);
            return;
	    }
        }

    /* Write region file for SAOimage overplotting */
    if (region_radius[0]) {
	double x, y, ddec, rmax, min_mag, max_mag, magscale;
	int radius, ix, iy;
	char snum[32], rstr[16];
	if (region_radius[icat] == 0) {
	    if (icat > 0)
		region_radius[icat] = region_radius[icat - 1];
	    else
		region_radius[icat] = 15;
	    }
	if (region_char[icat] == WCS_VAR)
	    strcat (title, " (+ stars, x nonstars)");
	else if (region_radius[icat] > 0) {
	    sprintf (temp, " (%d\" radius)", region_radius[icat]);
	    strcat (title, temp);
	    }
	fprintf (fd, "# %s\n", title);
	ddec = (double)region_radius[icat] / 3600.0;
	if (region_radius[icat] < 0) {
	    max_mag = gm[0];
	    min_mag = gm[0];
	    for (i = 0; i < nbg; i++) {
		if (gm[i] > max_mag) max_mag = gm[i];
		if (gm[i] < min_mag) min_mag = gm[i];
		}
	    if (max_mag == min_mag)
		rmax = 0;
	    else {
		rmax = (wcs->nxpix + wcs->nypix) / 50;
		if (rmax < 10)
		    rmax = 5;
		else
		    rmax = rmax - 5;
		magscale = max_mag - min_mag;
		}
	    }
	if (region_char[icat] == 0) {
	    if (icat > 0)
		region_char[icat] = region_char[icat - 1] + 1;
	    else
		region_char[icat] = WCS_CIRCLE;
	    }
	switch (region_char[icat]) {
	    case WCS_SQUARE:
		strcpy (rstr, "SQUARE");
		break;
	    case WCS_DIAMOND:
		strcpy (rstr, "DIAMOND");
		break;
	    case WCS_CROSS:
		strcpy (rstr, "CROSS");
		break;
	    case WCS_EX:
		strcpy (rstr, "EX");
		break;
	    case WCS_CIRCLE:
	    default:
		strcpy (rstr, "CIRCLE");
	    }

	for (i = 0; i < nbg; i++) {
	    if (gx[i] > 0.0 && gy[i] > 0.0) {
		if (region_radius[icat] > 0) {
		    wcs2pix (wcs, gra[i], gdec[i]+ddec, &x, &y, &offscale);
		    radius = (int) (sqrt ((x-gx[i])*(x-gx[i]) +
				          (y-gy[i])*(y-gy[i])) + 0.5);
		    }
		else if (rmax == 0)
		    radius = 20;
		else
		    radius = 5 + (int) (rmax * (max_mag - gm[i]) / magscale);
		ix = (int)(gx[i] + 0.5);
		iy = (int)(gy[i] + 0.5);
		printobj = 0;
		if (refcat == UAC  || refcat == UA1  || refcat == UA2 ||
		    refcat == USAC || refcat == USA1 || refcat == USA2)
		    sprintf (snum,"%13.8f", gnum[i]);
		else if (refcat == UJC)
		    sprintf (snum,"%12.7f", gnum[i]);
		else if (refcat == GSC)
		    sprintf (snum,"%9.4f", gnum[i]);
		else if (refcat==SAO || refcat==PPM || refcat==IRAS )
		    sprintf (snum,"%d", (int)gnum[i]);
		else if (refcat == TYCHO || refcat == HIP || refcat == ACT)
		    sprintf (snum,"%10.5f", gnum[i]);
		else {
		    if (gobj != NULL && gobj[i] != NULL) {
			if (strlen (gobj[i]) < 32)
			    strcpy (snum, gobj[i]);
			else
			    strncpy (snum, gobj[i], 31);
			printobj = 1;
			}
		    else
			sprintf (snum,"%d", (int)gnum[i]);
		    }
		if (region_char[icat] == WCS_VAR) {
		    if (gc[i] == 0)
			strcpy (rstr, "CROSS");
		    else
			strcpy (rstr, "EX");
		    }
		if (printobj)
		    fprintf (fd, "%s(%d,%d,%d) # %s\n",
			     rstr, ix, iy, radius, snum);
		else
		    fprintf (fd, "%s(%d,%d,%d) # %s %s\n",
			     rstr, ix, iy, radius, refcatname[icat], snum);
		}
	    }
	if (icat == ncat-1)
	    printf ("%s\n", outfile);
	continue;
	}

    /* Write heading */
    if (wfile)
	fprintf (fd,"%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    /* Set degree flag for output */
    if (sysout == WCS_ECLIPTIC || sysout == WCS_GALACTIC)
	degout = 1;
    else
	degout = degout0;

    if (refcat == GSC) 
	sprintf (headline, "CATALOG     HSTGSC1.1");
    else if (refcat == UAC)
	sprintf (headline, "CATALOG     USNO A");
    else if (refcat == UA1)
	sprintf (headline, "CATALOG     USNO A-1.0");
    else if (refcat == UA2)
	sprintf (headline, "CATALOG     USNO A-2.0");
    else if (refcat == USAC)
	sprintf (headline, "CATALOG     USNO SA");
    else if (refcat == USA1)
	sprintf (headline, "CATALOG     USNO SA-1.0");
    else if (refcat == USA2)
	sprintf (headline, "CATALOG     USNO SA-2.0");
    else if (refcat == UJC)
	sprintf (headline, "CATALOG     USNO UJ 1.0");
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
	sprintf (headline, "CATALOG     %s", refcatname[icat]);
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

    if (rasort) {
	if (wfile)
	    fprintf (fd, "RASORT	T\n");
	if (tabout)
	    printf ("RASORT	T\n");
	}

    if (wcs->sysout == WCS_B1950)
	sprintf (headline, "EQUINOX	1950.0");
    else
	sprintf (headline, "EQUINOX	2000.0");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (wfile)
    if (tabout)

    if (refcat == GSC)
	strcpy (headline,"gsc_id    	");
    else if (refcat == UAC)
	strcpy (headline,"uac_id    	");
    else if (refcat == UA1)
	strcpy (headline,"ua1_id    	");
    else if (refcat == UA2)
	strcpy (headline,"ua2_id    	");
    else if (refcat == USAC)
	strcpy (headline,"usac_id    	");
    else if (refcat == USA1)
	strcpy (headline,"usa1_id    	");
    else if (refcat == USA2)
	strcpy (headline,"usa2_id    	");
    else if (refcat == UJC)
	strcpy (headline,"ujc_id    	");
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
	strcpy (headline,"id    	");
    if (sysout == WCS_B1950)
	strcat (headline,"ra1950   	dec1950      	");
    else if (sysout == WCS_ECLIPTIC)
	strcat (headline,"long_ecl 	lat_ecl       	");
    else if (sysout == WCS_GALACTIC)
	strcat (headline,"long_gal  	lat_gal       	");
    else
	strcat (headline,"ra      	dec           	");
    if (refcat == UAC  || refcat == UA1  || refcat == UA2 ||
	refcat == USAC || refcat == USA1 || refcat == USA2)
	strcat (headline,"magb  	magr  	x    	y    	");
    else if (refcat == TYCHO || refcat == HIP || refcat == ACT)
	strcat (headline,"magb  	magv  	x    	y    	");
    else
	strcat (headline,"mag   	x    	y    	");
    if (refcat == GSC)
	strcat (headline,"type	");
    else if (refcat == UAC  || refcat == UA1  || refcat == UA2 ||
	    refcat == USAC || refcat == USA1 || refcat == USA2 || refcat == UJC)
	strcat (headline,"plate	");
    else if (refcat==SAO || refcat==PPM || refcat==IRAS || refcat==TYCHO ||
	     refcat==HIP || refcat==ACT)
	strcat (headline,"type 	");
    else
	strcat (headline,"peak	");
    strcat (headline,"dist");
    if (refcat == TABCAT && keyword != NULL) {
	strcat (headline,"	");
	strcat (headline, keyword);
	}
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (refcat == UAC  || refcat == UA1  || refcat == UA2 || refcat == HIP ||
	refcat == USAC || refcat == USA1 || refcat == USA2 || refcat == TYCHO ||
	refcat == ACT)
	sprintf(headline,"----------	--------	---------	----	-----	-----	-----	-----	----");
    else
        sprintf (headline,"----------	------------	------------	------	----	-------	----");
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
	    else if (refcat == UAC)
		printf ("No USNO A Stars Found\n");
	    else if (refcat == UA1)
		printf ("No USNO A-1.0 Stars Found\n");
	    else if (refcat == UA2)
		printf ("No USNO A-2.0 Stars Found\n");
	    else if (refcat == USAC)
		printf ("No USNO SA Stars Found\n");
	    else if (refcat == USA1)
		printf ("No USNO SA-1.0 Stars Found\n");
	    else if (refcat == USA2)
		printf ("No USNO SA-2.0 Stars Found\n");
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
		printf ("SAO number ");
	    else if (refcat == PPM)
		printf ("PPM number ");
	    else if (refcat == IRAS)
		printf ("IRAS number");
	    else if (refcat == TYCHO)
		printf ("Tycho number ");
	    else if (refcat == HIP)
		printf ("Hip number   ");
	    else if (refcat == ACT)
		printf ("ACT number   ");
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
			printf (" RA2000       Dec2000    ");
		    else
			printf ("RAJ%7.2f   DecJ%7.2f  ", eqout, eqout);
		    }
		}
	    if (refcat == UAC  || refcat == UA1  || refcat == UA2 ||
		refcat == USAC || refcat == USA1 || refcat == USA2)
		printf ("MagB  MagR    X      Y   Plate Arcsec\n");
	    else if (refcat == UJC)
		printf ("  Mag     X      Y   Plate Arcsec\n");
	    else if (refcat == GSC)
		printf (" Mag     X      Y   Class Arcsec\n");
	    else if (refcat == SAO || refcat == PPM || refcat == IRAS)
		printf (" Mag     X      Y   Type  Arcsec\n");
	    else if (refcat == TYCHO || refcat == HIP || refcat == ACT)
		printf (" MagR MagV     X      Y   Type  Arcsec\n");
	    else if (refcat == TABCAT && keyword != NULL)
		printf (" Mag     X      Y   Peak  Arcsec  %s\n", keyword);
	    else if (refcat == BINCAT)
		printf (" Mag     X      Y   Type  Arcsec   Object\n");
	    else if (refcat == TXTCAT)
		printf (" Mag     X      Y   Arcsec   Object\n");
	    else
		printf (" Mag     X      Y   Peak  Arcsec\n");
	    }
	}

    /* Print positions from reference catalog */
    for (i = 0; i < nbg; i++) {
	if (gx[i] > 0.0 && gy[i] > 0.0) {
	    gdist = 3600.0 * wcsdist (cra, cdec, gra[i], gdec[i]);
	    if (degout) {
		deg2str (rastr, 32, gra[i], 5);
		deg2str (decstr, 32, gdec[i], 5);
		}
	    else {
		ra2str (rastr, 32, gra[i], 3);
		dec2str (decstr, 32, gdec[i], 2);
		}
	    if (tabout || wfile) {
		if (refcat == GSC)
		    sprintf (headline, "%9.4f	%s	%s	%.2f	%.1f	%.1f	%d	%.2f",
		     gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i], gdist);
		else if (refcat == UAC  || refcat == UA1  || refcat == UA2 ||
			 refcat == USAC || refcat == USA1 || refcat == USA2)
		    sprintf (headline, "%13.8f	%s	%s	%.1f	%.1f	%.1f	%.1f	%d	%.2f",
		     gnum[i],rastr,decstr,gmb[i],gm[i],gx[i],gy[i],gc[i],gdist);
		else if (refcat == UJC)
		    sprintf (headline, "%12.7f	%s	%s	%.2f	%.1f	%.1f	%d	%.2f",
		     gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i], gdist);
		else if (refcat==SAO || refcat==PPM || refcat==IRAS ) {
		    isp[0] = gc[i] / 1000;
		    isp[1] = gc[i] % 1000;
		    sprintf (headline, "%9d	%s	%s	%.2f	%.1f	%.1f	%2s	%.2f",
		     (int)(gnum[i]+0.5), rastr, decstr, gm[i], gx[i], gy[i], isp, gdist);
		    }
		else if (refcat==TYCHO || refcat == HIP || refcat == ACT) {
		    isp[0] = gc[i] / 1000;
		    isp[1] = gc[i] % 1000;
		    sprintf (headline, "%10.5f	%s	%s	%.2f	%.2f	%.1f	%.1f	%2s	%.2f",
		     gnum[i], rastr, decstr, gmb[i], gm[i], gx[i], gy[i], isp, gdist);
		    }
		else if (refcat == BINCAT && nndec > 0) {
		    sprintf (nform, "%%%d.%df	%%s	%%s	%%.2f	%%.1f	%%.1f	%%d	%%.2f",
			     nndec+5, nndec);
		    sprintf (headline, nform, gnum[i], rastr, decstr, gm[i],
			     gx[i], gy[i], gc[i], gdist);
		    }
		else if (refcat == BINCAT)
		    sprintf (headline, "%10.0f	%s	%s	%.2f	%.1f	%.1f	%d	%.2f",
		     gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i], gdist);
		else if (nndec > 0) {
		    sprintf (nform, "%%%d.%df	%%s	%%s	%%.2f	%%.1f	%%.1f	%%d	%%.2f",
			     nndec+5, nndec);
		    sprintf (headline, nform, gnum[i], rastr, decstr, gm[i],
			     gx[i], gy[i], gc[i], gdist);
		    }
		else
		    sprintf (headline, "%9.4f	%s	%s	%.2f	%.1f	%.1f	%d	%.2f",
		     gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i], gdist);
		if (wfile)
		    fprintf (fd, "%s\n", headline);
		else if (tabout)
		    printf ("%s\n", headline);
		}
	    else if (!tabout) {
		if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		    refcat == UAC  || refcat == UA1  || refcat == UA2)
		    sprintf (headline,"%13.8f %s %s %5.1f %5.1f %.1f %.1f %4d ",
			gnum[i],rastr,decstr,gmb[i],gm[i],gx[i],gy[i],gc[i]);
		else if (refcat == UJC)
		    sprintf (headline,"%12.7f %s %s %6.2f %.1f %.1f %4d",
			gnum[i], rastr, decstr, gm[i],gx[i],gy[i],gc[i]);
		else if (refcat == GSC)
		    sprintf (headline,"%9.4f %s %s %6.2f %.1f %.1f %2d",
			gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i]);
		else if (refcat==SAO || refcat==PPM || refcat==IRAS ) {
		    isp[0] = gc[i] / 1000;
		    isp[1] = gc[i] % 1000;
		    sprintf (headline,"%8d  %s %s %6.2f  %2s",
			(int)(gnum[i]+0.5),rastr,decstr,gm[i],isp);
		    sprintf (temp, "  %.1f  %.1f", gx[i], gy[i]);
		    strcat (headline, temp);
		    }
		else if (refcat == TYCHO || refcat == HIP || refcat == ACT) {
		    isp[0] = gc[i] / 1000;
		    isp[1] = gc[i] % 1000;
		    sprintf (headline,"%10.5f %s %s %6.2f %6.2f  %2s",
			gnum[i],rastr,decstr,gmb[i],gm[i],isp);
		    sprintf (temp, "  %.1f  %.1f", gx[i], gy[i]);
		    strcat (headline, temp);
		    }
		else if (refcat == TABCAT)
		    sprintf (headline,"%9.4f %s %s %6.2f %.1f %.1f%7d",
			gnum[i], rastr, decstr, gm[i],gx[i],gy[i],gc[i]);
		else if (refcat == BINCAT) {
		    isp[0] = gc[i] / 1000;
		    isp[1] = gc[i] % 1000;
		    if (nndec > 0) {
			sprintf (nform,"%%%d.%df %%s %%s %%6.2f %%2s",
				 nndec+5,nndec);
			sprintf (headline,nform,gnum[i],rastr,decstr,gm[i],isp);
			}
		    else
			sprintf (headline,"%8d %s %s %6.2f %2s",
				(int)(gnum[i]+0.5), rastr, decstr, gm[i], isp);
		    sprintf (temp, "  %.1f  %.1f", gx[i], gy[i]);
		    strcat (headline, temp);
		    }
		else {
		    if (nndec > 0) {
			sprintf (nform,"%%%d.%df %%s %%s %%6.2f",
				 nndec+5,nndec);
			sprintf (headline,nform, gnum[i], rastr, decstr, gm[i]);
			}
		    else
			sprintf (headline, "%8d %s %s %6.2f",
				(int)(gnum[i]+0.5), rastr, decstr, gm[i]);
		    sprintf (temp, "  %.1f  %.1f", gx[i], gy[i]);
		    strcat (headline, temp);
		    }
		sprintf (temp, "  %7.2f", gdist);
		strcat (headline, temp);
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
    if (gobj) free ((char *)gobj);
    free ((char *)wcs);

    return;
}

/* May 21 1996	New program
 * Jul 11 1996	Update file reading
 * Jul 16 1996	Remove unused image pointer; do not free unallocated header
 * Aug 15 1996	Clean up file reading code
 * Sep  4 1996	Free header immediately after use
 * Oct  1 1996	Set file extension according to the catalog which is used
 * Oct  1 1996	Write output file only if flag is set
 * Oct 15 1996	Use GetFITSWCS instead of local code
 * Oct 16 1996	Write list of stars to standard output by default
 * Nov 13 1996	Add UA catalog reading capabilities
 * Nov 15 1996	Change catalog reading subroutine arguments
 * Dec 10 1996	Change equinox in getfitswcs call to double
 * Dec 12 1996	Version 1.2
 * Dec 12 1996	Add option for bright as well as faint magnitude limits
 * Dec 12 1996	Fix header for UAC magnitudes
 * Dec 13 1996	Write plate into header if selected
 *
 * Jan 10 1997	Fix bug in RASort Stars which did not sort magnitudes
 * Feb 21 1997  Get image header from GetFITSWCS()
 * Mar 14 1997	Add support for USNO SA-1.0 catalog
 * Apr 25 1997	Fix bug in uacread
 * May 28 1997  Add option to read a list of filenames from a file
 * Nov 17 1997	Initialize both magnitude limits
 * Dec  8 1997	Set up program to be called by various names
 * Dec 15 1997	Add capability of reading and writing IRAF 2.11 images
 *
 * Jan 27 1998  Implement Mark Calabretta's WCSLIB
 * Jan 27 1998  Add -z for AIPS classic WCS projections
 * Feb 18 1998	Version 2.0: Full Calabretta WCS implementation
 * Mar 27 1998	Version 2.1: Add IRAF TNX projection
 * Apr 14 1998	Version 2.2: Add polynomial plate fit
 * Apr 24 1998	change coordinate setting to setsys() from setfk4()
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * May 13 1998	If nstars is set use it as a limit no matter how small
 * May 27 1998	Do not include fitshead.h
 * Jun  2 1998	Fix bug in tabread()
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Jun 25 1998	Set WCS subroutine choice with SETDEFWCS()
 * Jul  8 1998	Add other coordinate types
 * Jul  9 1998	Adjust report headings
 * Sep 10 1998	Add SAOTDC binary format catalogs
 * Sep 16 1998	Add coordinate system and equinox to binread()
 * Sep 17 1998	Add coordinate system to GetFITSWCS()
 * Sep 21 1998	Add epoch to heading
 * Sep 25 1998	Add system, equinox, and epoch to catalog search calls
 * Oct  9 1998	Add option to write SAOimage region file
 * Oct  9 1998	Add option to read arbitrary SAO binary catalog file
 * Oct  9 1998	Add source ID in comment after SAOimage region output
 * Oct 13 1998	Make sure refcatname is always set
 * Oct 14 1998	Use isiraf() to determine file type
 * Oct 15 1998	Add TDC ASCII catalog access
 * Oct 19 1998	Add variable SAOimage region shapes
 * Oct 19 1998	Add magnitude-scaled SAOimage region size
 * Oct 21 1998	Add object name to binary catalogs
 * Oct 22 1998	Use RefCat() to set type of catalog name
 * Oct 23 1998	Allow searches of multiple catalogs into one output file
 * Oct 26 1998	Return object name in same operation as object position
 * Oct 27 1998	Move region shape codes to wcscat.h
 * Oct 29 1998	Add GSC class to output header
 * Oct 29 1998	Add tab table keyword to output
 * Nov 20 1998	Add support for USNO A-2.0 and SA-2.0 catalogs
 * Nov 30 1998	Add x command for new reference pixel
 * Nov 30 1998	Add version and help commands for consistency
 * Dec  8 1998	Add support for Hipparcos and ACT catalogs
 * Dec 21 1998	Fix formats for text catalogs
 * Dec 21 1998	Write output file to current working directory
 *
 * Jan 25 1999	Add -i for IRAF formatted output (X Y RA Dec)
 * Jan 26 1999	Drop -i; add similar feature to immatch
 * Feb 12 1999	Finish adding support for ACT catalog
 * Feb 18 1998	Add variable number of decimal places to TDC catalog output
 * Mar  2 1999	Add x and y to non-tab output (bug fix)
 * Apr  7 1999	Add filename argument to GetFITSWCS
 * Apr 13 1999	Fix progname to drop / when full pathname
 * Apr 20 1999	Fix minor bug in character assignment code
 * May 12 1999	Adjust command listing
 */
