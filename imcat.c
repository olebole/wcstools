/* File imcat.c
 * April 2, 2003
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

#include "libwcs/wcs.h"
#include "libwcs/wcscat.h"
#include "libwcs/lwcs.h"
#include "libwcs/fitsfile.h"

#define MAXFILES 1000
static int maxnfile = MAXFILES;

static void PrintUsage();
static void ListCat();
extern void fk524e();
extern struct WorldCoor *GetFITSWCS();
extern char *GetFITShead();
extern void setsys();
extern void setcenter();
extern void setsecpix();
extern void setrefpix();
extern void setdateobs();
extern void setparm();


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
static int catsort = SORT_MAG;	/* Default to sort stars by magnitude */
static int debug = 0;		/* True for extra information */
static int degout0 = 0;		/* True for RA and Dec in fractional degrees */
static char *keyword = NULL;	/* Column to add to tab table output */
static int sysout = 0;		/* Output coordinate system */
static double eqout = 0.0;	/* Equinox for output coordinates */
static int version = 0;		/* If 1, print only program name and version */
static int obname[5];		/* If 1, print object name, else number */
static struct StarCat *starcat[5]; /* Star catalog data structure */
static int nmagmax = 4;
static int sortmag = 0;		/* Magnitude by which to sort stars */
static webdump = 0;
static char *progname;		/* Name of program as executed */
static int minid = 0;		/* Minimum number of plate IDs for USNO-B1.0 */
static int minpmqual = 0;	/* Minimum USNO-B1.0 proper motion quality */
extern int getminpmqual();
extern int getminid();
static void ImageLim();

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
    char errmsg[256];
    FILE *flist;
    char *listfile;
    char *cstr;
    char cs, cs1;
    char **fn;
    int i, ic;
    double x, y;
    char *refcatname[5];	/* reference catalog name */
    int ncat = 0;
    int region_radius[5];	/* Flag for SAOimage region file output */
    int rcat = 0;
    int region_char[5];		/* Character for SAOimage region file output */
    int region_pixel;
    char *refcatn;
    int ifile, nfile;
    int lcat;
    int scat = 0;
    char c1, c, *ccom;
    double drot;

    nfile = 0;
    fn = (char **)calloc (maxnfile, sizeof(char *));

    for (i = 0; i< 5; i++)
	starcat[i] = NULL;

    for (i = 0; i < 5; i++) {
	region_radius[i] = 0;
	region_char[i] = 0;
	obname[i] = 0;
	}

    /* Check name used to execute programe and set catalog name accordingly */
    progname = ProgName (av[0]);
    refcatn = ProgCat (progname);
    if (refcatn != NULL) {
	refcatname[ncat] = refcatn;
	ncat++;
	}

    if (ac == 1)
	PrintUsage (NULL);

    /* Loop through the arguments */
    for (av++; --ac > 0; av++) {
	str = *av;

	/* Check for help command first */
	if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	    PrintUsage (NULL);

	/* Check for version command */
	else if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	    version = 1;
	    PrintUsage ("version");
	    }

	else if (strchr (str, '=')) {
	    setparm (str);
	    minid = getminid ();
	    minpmqual = getminpmqual ();
	    }

	/* Image list file */
	else if (str[0] == '@') {
	    readlist++;
	    listfile = ++str;
	    str = str + strlen (str) - 1;
	    av++;
	    ac--;
	    }

	/* Decode arguments */
	else if (str[0] == '-') {

	    while ((c = *++str) != 0) {
    		switch (c) {

		case 'a':       /* Initial rotation angle in degrees */
		    if (ac < 2)
			PrintUsage (str);
		    drot = atof (*++av);
                    setrot (drot);
                    ac--;
                    break;

		case 'b':	/* initial coordinates on command line in B1950 */
		    str1 = *(av+1);
		    ic = (int)str1[0];
		    if (*(str+1) || ic < 48 || ic > 58) {
			setsys(WCS_B1950);
			sysout = WCS_B1950;
			eqout = 1950.0;
			}
		    else if (ac < 3)
			PrintUsage (str);
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
			PrintUsage (str);
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

		case 'i':	/* Label region with name, not number */
		    if (ncat > 0)
			obname[ncat-1]++;
		    else
			obname[0]++;
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
			PrintUsage (str);
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
			PrintUsage (str);
		    keyword = *++av;
		    settabkey (keyword);
		    if (ncat > 0)
			obname[ncat-1]++;
		    else
			obname[0]++;
		    ac--;
		    break;

		case 'm':	/* Limiting reference star magnitude */
		    if (ac < 2)
			PrintUsage (str);
		    cs1 = *(str+1);
		    if (cs1 != (char) 0) {
			++str;
			if (cs1 > '9')
			    sortmag = (int) cs1;
			else
			    sortmag = (int) cs1 - 48;
			}
		    ac--;
		    av++;
		    if ((ccom = strchr (*av, ',')) != NULL) {
			*ccom = (char) 0;
			maglim1 = atof (*av);
			maglim2 = atof (ccom+1);
			}
		    else {
			maglim2 = atof (*av);
			if (ac > 1 && isnum (*(av+1))) {
			    av++;
			    ac--;
			    maglim1 = maglim2;
			    maglim2 = atof (*av);
			    }
		        else if (MAGLIM1 == MAGLIM2)
			    maglim1 = -2.0;
			}
		    break;

		case 'n':	/* Number of brightest stars to read */
		    if (ac < 2)
			PrintUsage (str);
		    nstars = atoi (*++av);
		    ac--;
		    break;

		case 'o':	/* Guide Star object class */
		    if (ac < 2)
			PrintUsage (str);
		    classd = (int) atof (*++av);
		    setgsclass (classd);
		    ac--;
		    break;

		case 'p':	/* Initial plate scale in arcseconds per pixel */
		    if (ac < 2)
			PrintUsage (str);
		    setsecpix (atof (*++av));
		    ac--;
		    break;

		case 'q':	/* Output region file shape for SAOimage */
		    if (ac < 2)
			PrintUsage (str);
		    cstr = *++av;
		    c1 = cstr[0];
		    region_pixel = 0;
		    if (strlen(cstr) > 1) {
			if (cstr[0] == 'p') {
			    c1 = cstr[1];
			    region_pixel = 10;
			    }
			else if (cstr[1] == 'p')
			    region_pixel = 10;
			}
		    switch (c1){
			case 'c':
			    if (cstr[1] == 'i')
				region_char[scat] = WCS_CIRCLE;
			    else
				region_char[scat] = WCS_CROSS;
			    break;
			case 'd':
			    region_char[scat] = WCS_DIAMOND;
			    break;
			case 'p':
			    region_char[scat] = WCS_PCIRCLE;
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
		    region_char[scat] = region_char[scat] + region_pixel;
		    if (region_radius[scat] == 0)
			region_radius[scat] = -1;
		    scat++;
		    wfile++;
		    ac--;
		    break;

		case 'r':	/* Output region file with shape radius for SAOimage */
		    if (ac < 2)
			PrintUsage (str);
		    region_radius[rcat] = atoi (*++av);
		    if (region_radius[rcat] == 0)
			region_radius[rcat] = -1;
		    rcat++;
		    wfile++;
		    ac--;
		    break;

		case 's':	/* sort by RA, Dec, magnitude or nothing */
		    catsort = SORT_RA;
		    if (ac > 1) {
			str1 = *(av + 1);
			cs = str1[0];
			if (strchr ("dmnrxy",(int)cs)) {
			    cs1 = str1[1];
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
			else if (cs == 'm') {
			    catsort = SORT_MAG;
			    if (cs1 != (char) 0) {
				if (cs1 > '9')
				    sortmag = (int) cs1;
				else
				    sortmag = (int) cs1 - 48;
				}
			    }

			/* No sorting */
			else if (cs == 'n')
			    catsort = NOSORT;

			/* X coordinate */
			else if (cs == 'x')
			    catsort = SORT_X;

			/* Y coordinate */
			else if (cs == 'y')
			    catsort = SORT_Y;

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

		case 'u':	/* UJ Catalog plate number */
		    if (ac < 2)
			PrintUsage (str);
		    uplate = (int) atof (*++av);
		    setuplate (uplate);
		    ac--;
		    break;

		case 'v':	/* more verbosity */
		    if (debug) {
			webdump++;
			debug = 0;
			verbose = 0;
			}
		    else if (verbose)
			debug++;
		    else
			verbose++;
		    break;

		case 'w':	/* write output file */
		    wfile++;
		    break;

		case 'x':	/* X and Y coordinates of reference pixel */
		    if (ac < 3)
			PrintUsage (str);
		    x = atof (*++av);
		    ac--;
		    y = atof (*++av);
		    ac--;
		    setrefpix (x, y);
		    break;

		case 'y':	/* Epoch of image in FITS date format */
		    if (ac < 2)
			PrintUsage (str);
		    setdateobs (*++av);
		    ac--;
		    break;

		case 'z':       /* Use AIPS classic WCS */
		    setdefwcs (WCS_ALT);
		    break;

		default:
		    sprintf (errmsg, "* Illegal command -%c-", c);
		    PrintUsage (errmsg);
		    break;
		}
		}
	    }

	/* Image file */
	else if (isfits (str) || isiraf (str)) {
	    if (nfile >= maxnfile) {
		maxnfile = maxnfile * 2;
		fn = (char **) realloc ((void *)fn, maxnfile);
		}
	    fn[nfile] = str;
	    nfile++;
	    }

	else {
	    sprintf (errmsg, "* %s is not a FITS or IRAF file.", str);
	    PrintUsage (errmsg);
	    }
	}

    /* if (!verbose && !wfile)
	verbose = 1; */

    /* Process image files from list file */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    sprintf (errmsg,"* List file %s cannot be read\n", listfile);
	    PrintUsage (errmsg);
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

    /* Process image files */
    else if (nfile > 0) {
	for (ifile = 0; ifile < nfile; ifile++) {
	    if ( verbose)
		printf ("%s:\n", fn[ifile]);
	    ListCat (progname,fn[ifile], ncat, refcatname, region_radius, region_char);
	    if (verbose)
		printf ("\n");
	    }
	}

    /* Print error message if no image files to process */
    else
	PrintUsage ("* No files to process.");

    /* Close source catalogs */
    for (i = 0; i < 5; i++)
	if (starcat[i] != NULL) ctgclose (starcat[i]);

    return (0);
}

static void
PrintUsage (command)
char	*command;
{
    if (version)
	exit (0);

    if (command != NULL) {
	if (command[0] == '*')
	    fprintf (stderr, "%s\n", command);
	else
	    fprintf (stderr, "* Missing argument for command: %c\n", command[0]);
	exit (1);
	}

    if (strsrch (progname,"gsca") != NULL) {
	fprintf (stderr,"List GSC-ACT Stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-g class]\n");
	}
    else if (strsrch (progname,"gsc2") != NULL) {
	fprintf (stderr,"List GSC II Stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-g class]\n");
	}
    else if (strsrch (progname,"gsc") != NULL) {
	fprintf (stderr,"List HST Guide Stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-g class]\n");
	}
    else if (strsrch (progname,"tmc") != NULL ||
	strsrch (progname,"2mp") != NULL) {
	fprintf (stderr,"List 2MASS Point Sources in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2]\n");
	}
    else if (strsrch (progname,"ujc") != NULL) {
	fprintf (stderr,"List USNO J Catalog stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-u plate]\n");
	}
    else if (strsrch (progname,"uac") != NULL) {
	fprintf (stderr,"List USNO A stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-u plate]\n");
	}
    else if (strsrch (progname,"ub1") != NULL) {
	fprintf (stderr,"List USNO B-1.0 stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-mn [mag1] mag2]\n");
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
    else if (strsrch (progname,"ub1") != NULL) {
	fprintf (stderr,"List USNO B-1.0 stars in FITS and IRAF image files\n");
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
	fprintf (stderr,"Usage: [-vwhst][-a deg][-m [mag1] mag2][-c catalog][-x x y]\n");
	}
    fprintf (stderr,"       [-p scale][-q osd+x][-b ra dec][-j ra dec][-r arcsec] FITS or IRAF file(s)\n");
    fprintf (stderr,"  -a: initial rotation angle in degrees (default 0)\n");
    fprintf (stderr,"  -b [RA Dec]: Output, (center) in B1950 (FK4) RA and Dec\n");
    fprintf (stderr,"  -c name: Reference catalog (gsc, ua2, local file, etc.\n");
    fprintf (stderr,"  -d: Output RA,Dec positions in fractional degrees\n");
    fprintf (stderr,"  -e: Output in ecliptic longitude and latitude\n");
    fprintf (stderr,"  -g: Output in galactic longitude and latitude\n");
    fprintf (stderr,"  -h: Print heading, else do not \n");
    fprintf (stderr,"  -i: Print name instead of number in region file \n");
    fprintf (stderr,"  -j [RA Dec]: Output (center) in J2000 (FK5) RA and Dec\n");
    fprintf (stderr,"  -k keyword: Add this keyword to output from tab table search\n");
    fprintf (stderr,"  -mx m1 [m2]: Catalog magnitude #x limit(s) (only one set allowd, default none)\n");
    fprintf (stderr,"  -n num: Number of brightest stars to print \n");
    fprintf (stderr,"  -o name: Set HST Guide Star object class to print \n");
    fprintf (stderr,"  -p num: Initial plate scale in arcsec per pixel (default 0)\n");
    fprintf (stderr,"  -q code: Write SAOimage region file of this shape (filename.cat)\n");
    fprintf (stderr,"  -r num: Write SAOimage region file of this radius (filename.cat)\n");
    fprintf (stderr,"  -s d|mx|n|r|x|y: Sort by r=RA d=Dec mx=Mag#x n=none x=X y=Y\n");
    fprintf (stderr,"  -t: Tab table to standard output as well as file\n");
    fprintf (stderr,"  -u num: USNO catalog single plate number to accept\n");
    fprintf (stderr,"  -v: Verbose\n");
    fprintf (stderr,"  -w: Write tab table output file [imagename].[catalog]\n");
    fprintf (stderr,"  -x x y: X and Y coordinates of reference pixel (default is center)\n");
    fprintf (stderr,"  -y date: Epoch of image in FITS date format or year\n");
    fprintf (stderr,"  -z: Use AIPS classic projections instead of WCSLIB\n");
    exit (1);
    fprintf (stderr,"   x: Number of magnitude must be same for sort and limits\n");
    fprintf (stderr,"      and x may be omitted from either or both -m and -s m\n");
}


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
    double *gnum;	/* Catalog numbers */
    double *gra;	/* Catalog right ascensions, degrees */
    double *gdec;	/* Catalog declinations, degrees */
    double *gpra;	/* Catalog right ascension proper motions, degrees/year */
    double *gpdec;	/* Catalog declination proper motions, degrees/year */
    double **gm;		/* Catalog star magnitudes */
    double *gx, *gy;	/* Catalog star positions on image */
    char **gobj;	/* Catalog object names */
    char **gobj1;	/* Catalog object names */
    int *gc;		/* Catalog object classes, plates, etc. */
    int ng;		/* Number of catalog stars */
    int nbg;		/* Number of brightest catalog stars actually used */
    int imh, imw;	/* Image height and width in pixels */
    int i, ngmax;
    int degout;
    FILE *fd;
    struct WorldCoor *wcs;	/* World coordinate system structure */
    double eqref;	/* Equinox of catalog to be searched */
    double epref;	/* Epoch of catalog to be searched */
    double epout;	/* Epoch of catalog to be searched */
    int sysref;		/* Coordinate system of catalog to be searched */
    char rastr[16], decstr[16];	/* coordinate strings */
    char numstr[32];	/* Catalog number */
    double cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, mag1, mag2,secpix;
    double mag, drad, flux;
    int offscale, nlog, imag;
    char headline[160];
    char temp[80];
    char title[80];
    char outfile[80];
    char *fname;
    char magname[8];
    char isp[4];
    int icat;
    int printobj = 0;
    int nndec, nnfld;
    char nform[64];
    char blanks[256];
    int lfn;
    int band;
    int ngsc;
    int gcset;
    int mprop;
    int nmag;
    int sptype;
    int magsort;
    double gxmax, gymax;
    double pra, pdec;
    double maxnum;
    double xmag, xmag1;
    char *catalog;

    /* Drop out if no catalog is specified */
    if (ncat < 1) {
	fprintf (stderr, "No catalog specified\n");
	exit (-1);
	}

    gnum = NULL;
    gra = NULL;
    gdec = NULL;
    gpra = NULL;
    gpdec = NULL;
    gm = NULL;
    gx = NULL;
    gy = NULL;
    gc = NULL;
    gobj = NULL;
    gobj1 = NULL;

    isp[2] = 0;
    isp[3] = 0;
    if (verbose || printhead)
    for (i = 0; i < 255; i++)
	blanks[i] = ' ';
    blanks[255] = (char) 0;

    /* Loop through catalogs */
    for (icat = 0; icat < ncat; icat++) {

	/* Skip this catalog if no name is given */
	if (refcatname[icat] == NULL || strlen (refcatname[icat]) == 0) {
	    fprintf (stderr, "Catalog %d not specified\n", icat);
	    continue;
	    }

    /* Find title and coordinate system for catalog */
    if (!(refcat = RefCat (refcatname[icat],title,&sysref,&eqref,&epref,&mprop,&nmag))) {
	fprintf (stderr,"ListCat: No catalog named %s\n", refcatname[icat]);
	return;
	}

    /* If more magnitudes are needed, allocate space for them */
    if (nmag > nmagmax)
	nmagmax = nmag;

    if (classd == 0)
	strcat (title, " stars");
    else if (classd == 3)
	strcat (title, " nonstars");

    /* Read world coordinate system information from the image header */
    if ((header = GetFITShead (filename, verbose)) == NULL)
	return;
    wcs = GetFITSWCS (filename, header, verbose, &cra, &cdec, &dra, &ddec,
		      &secpix, &imw, &imh, &sysout, &eqout);
    gxmax = (double) imw + 0.5;
    gymax = (double) imh + 0.5;
    free (header);
    if (nowcs (wcs)) {
	wcsfree (wcs);
	return;
	}

    /* Set up limits for search, taking into account image rotation */
    ImageLim (wcs,&cra, &cdec, &dra, &ddec, &ra1, &ra2, &dec1, &dec2);
    epout = wcs->epoch;
    if (verbose || printhead) {
	char rastr1[16],rastr2[16],decstr1[16],decstr2[16], cstr[16];
	wcscstr (cstr, sysout, eqout, epout);
	ra2str (rastr1, 16, ra1, 3);
	ra2str (rastr2, 16, ra2, 3);
	printf ("%s: RA:  %s - %s %s\n",filename,rastr1, rastr2, cstr);
	dec2str (decstr1, 16, dec1, 2);
	dec2str (decstr2, 16, dec2, 2);
	lfn = strlen (filename);
	blanks[lfn] = (char) 0;
	printf ("%s  Dec: %s - %s %s\n", blanks, decstr1, decstr2, cstr);
	blanks[lfn] = ' ';
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
    if (sortmag > 9)
	sortmag = CatMagNum (sortmag, refcat);

    if (nstars > 0)
	ngmax = nstars;
    else
	ngmax = MAXCAT;

    if (!(gnum = (double *) calloc (ngmax, sizeof (double))))
	fprintf (stderr, "Could not calloc %d bytes for gnum\n",
		 ngmax * sizeof (double));
    if (!(gra = (double *) calloc (ngmax, sizeof (double))))
	fprintf (stderr, "Could not calloc %d bytes for gra\n",
		 ngmax * sizeof (double));
    if (!(gdec = (double *) calloc (ngmax, sizeof (double))))
	fprintf (stderr, "Could not calloc %d bytes for gdec\n",
		 ngmax * sizeof (double));
    if (!(gm = (double **) calloc (nmagmax, sizeof(double *))))
	fprintf (stderr, "Could not calloc %d bytes for gm\n",
		 nmagmax*sizeof(double *));
    else {
	for (imag = 0; imag < nmagmax; imag++) {
	    if (!(gm[imag] = (double *) calloc (ngmax, sizeof(double))))
		fprintf (stderr, "Could not calloc %d bytes for gm\n",
			 ngmax*sizeof(double));
	    }
	}
    if (!(gc = (int *) calloc (ngmax, sizeof (int))))
	fprintf (stderr, "Could not calloc %d bytes for gc\n",
		 ngmax * sizeof (double));
    if (!(gx = (double *) calloc (ngmax, sizeof (double))))
	fprintf (stderr, "Could not calloc %d bytes for gx\n",
		 ngmax * sizeof (double));
    if (!(gy = (double *) calloc (ngmax, sizeof (double))))
	fprintf (stderr, "Could not calloc %d bytes for gy\n",
		 ngmax * sizeof (double));
    if (!(gobj = (char **) calloc (ngmax, sizeof (void *))))
	fprintf (stderr, "Could not calloc %d bytes for obj\n",
		 ngmax*sizeof(void *));
    if (!(gpra = (double *) calloc (ngmax, sizeof(double))))
	fprintf (stderr, "Could not calloc %d bytes for gpra\n",
		 ngmax*sizeof(double));
    if (!(gpdec = (double *) calloc (ngmax, sizeof(double))))
	fprintf (stderr, "Could not calloc %d bytes for gpdec\n",
		 ngmax*sizeof(double));

    if (!gnum || !gra || !gdec || !gm || !gc || !gx || !gy || !gobj ||
	!gpra || !gpdec) {
	if (gm) {
	    for (imag = 0; imag < nmagmax; imag++)
		free ((char *) gm[imag]);
	    free ((char *)gm);
	    gm = NULL;
	    }
	if (gra) free ((char *)gra);
	if (gdec) free ((char *)gdec);
	if (gpra) free ((char *)gpra);
	if (gpdec) free ((char *)gpdec);
	if (gnum) free ((char *)gnum);
	if (gc) free ((char *)gc);
	if (gx) free ((char *)gx);
	if (gy) free ((char *)gy);
	if (gobj) free ((char *)gobj);
	wcsfree (wcs);
	return;
	}
    if (webdump)
	nlog = -1;
    else if (verbose) {
	if (refcat == UAC  || refcat == UA1  || refcat == UA2 || refcat==UB1 ||
	    refcat == USAC || refcat == USA1 || refcat == USA2 ||
	    refcat == GSC  || refcat == GSCACT || refcat == TMPSC)
	    nlog = 1000;
	else
	    nlog = 100;
	}
    else
	nlog = 0;

    /* Find the nearby reference stars, in ra/dec */
    drad = 0.0;
    ng = ctgread (refcatname[icat], refcat, 0,
		  cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,mag2,
		  sortmag,ngmax,&starcat[icat],
		  gnum,gra,gdec,gpra,gpdec,gm,gc,gobj,nlog);

    /* Set flag if any proper motions are non-zero
    mprop = 0;
    for (i = 0; i < ng; i++) {
	if (gpra[i] != 0.0 || gpdec[i] != 0.0) {
	    mprop = 1;
	    break;
	    }
	} */

    /* Set flag if spectral type is present */
    if (refcat==SAO || refcat==PPM || refcat==IRAS || refcat==TYCHO ||
	refcat==HIP || refcat==BSC)
	sptype = 1;
    else if (refcat == TMPSC)
	sptype = 2;
    else if (starcat[icat] != NULL && starcat[icat]->sptype > 0)
	sptype = 1;
    else
	sptype = 0;

    if (refcat == BINCAT || refcat == TABCAT || refcat == TXTCAT)
	nndec = starcat[icat]->nndec;

    /* Find out whether object names are set */
    if (gobj[0] == NULL)
	gobj1 = NULL;
    else
	gobj1 = gobj;

    if (ng > ngmax)
	nbg = ngmax;
    else
	nbg = ng;

    /* Find largest catalog number printed */
    maxnum = 0.0;
    for (i = 0; i < nbg; i++ ) {
	if (gnum[i] > maxnum)
	    maxnum = gnum[i];
	}
    nnfld = CatNumLen (refcat, maxnum, nndec);

    /* Get image pixel coordinates for each star found in reference catalog */
    for (i = 0; i < nbg; i++ ) {
	offscale = 0;
	wcs2pix (wcs, gra[i], gdec[i], &gx[i], &gy[i], &offscale);
	if (offscale) {
	    gx[i] = 0.0;
	    gy[i] = 0.0;
	    }
	}

    /* Check to see whether gc is set at all */
    gcset = 0;
    for (i = 0; i < nbg; i++ ) {
	if (gc[i] != 0) {
	    gcset = 1;
	    break;
	    }
	}

    /* Sort reference stars by brightness (magnitude) */
    MagSortStars (gnum, gra, gdec, gpra, gpdec, gx, gy, gm, gc, gobj1, nbg,
		  nmag, sortmag);

    /* List the brightest reference stars */
    CatMagName (sortmag, refcat, magname);
    if (sortmag > 0 && sortmag <= nmag)
	magsort = sortmag - 1;
    else
	magsort = 0;
    if (ng > ngmax) {
	if (verbose || printhead) {
	    if (mag2 > 0.0)
		printf ("%d / %d %s %.1f < %s < %.1f",
			nbg,ng,title,gm[magsort][0],magname,gm[magsort][nbg-1]);
	    else
		printf ("%d / %d %s %s <= %.1f",
			nbg, ng, title, magname, gm[magsort][nbg-1]);
	    }
	}
    else {
	if (verbose || printhead) {
	    if (maglim1 > 0.0)
		printf ("%d %s %.1f < %s < %.1f",
			ng,title,maglim1,magname,maglim2);
	    else if (maglim2 > 0.0)
		printf ("%d %s %s < %.1f",ng, title, magname, maglim2);
	    else
		printf ("%d %s", ng, title);
	    }
	}
    if (printhead || verbose) {
	if (mprop) {
	    if (wcs->epoch != wcs->equinox)
		printf (" at %7.2f",wcs->epoch);
	    if (isiraf (filename))
		printf (" in IRAF image %s\n",filename);
	    else
		printf (" in FITS image %s\n", filename);
	    }
	else
	    printf ("\n");
	}

    /* Sort catalogued objects, if requested */
    if (nbg > 1) {

	/* Sort reference stars by image X coordinate */
	if (catsort == SORT_X)
	    XSortStars (gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gc,gobj1,nbg,
			 nmagmax);

	/* Sort reference stars by image Y coordinate */
	if (catsort == SORT_Y)
	    YSortStars (gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gc,gobj1,nbg,
			 nmagmax);

	/* Sort star-like objects in image by right ascenbgion */
	else if (catsort == SORT_RA)
	    RASortStars (gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gc,gobj1,nbg,
			 nmagmax);

	/* Sort star-like objects in image by declination */
	else if (catsort == SORT_DEC)
	    DecSortStars(gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gc,gobj1,nbg,
			 nmagmax);

	/* Sort reference stars from brightest to faintest */
	else if (catsort == SORT_MAG) {
	    MagSortStars(gnum,gra,gdec,gpra,gpdec,gx,gy,gm,gc,gobj1,nbg,
			 nmagmax,sortmag);
	    }
	}

    sprintf (headline, "image	%s", filename);

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
	    if (gm) {
		for (imag = 0; i < nmagmax; i++)
		    free ((char *)gm[imag]);
		free ((char *)gm);
		}
	    if (gra) free ((char *)gra);
	    if (gdec) free ((char *)gdec);
	    if (gpra) free ((char *)gpra);
	    if (gpdec) free ((char *)gpdec);
	    if (gnum) free ((char *)gnum);
	    if (gc) free ((char *)gc);
	    if (gx) free ((char *)gx);
	    if (gy) free ((char *)gy);
	    if (gobj) free ((char *)gobj);
	    wcsfree (wcs);
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
		region_radius[icat] = 20.0 * wcs->xinc / 3600.0;
	    }
	else if (region_radius[icat] < -1)
	    region_radius[icat] = -region_radius[icat] * wcs->xinc / 3600.0;
	if (region_char[icat] == WCS_VAR)
	    strcat (title, " (+ stars, x nonstars)");
	else if (region_radius[icat] > 0 && region_char[icat] > 10) {
	    sprintf (temp, " (%d pixel radius)", region_radius[icat]);
	    strcat (title, temp);
	    }
	else if (region_radius[icat] > 0) {
	    sprintf (temp, " (%d\" radius)", region_radius[icat]);
	    strcat (title, temp);
	    }
	fprintf (fd, "# %s\n", title);
	ddec = (double)region_radius[icat] / 3600.0;
	if (region_radius[icat] < 0) {
	    max_mag = gm[magsort][0];
	    min_mag = gm[magsort][0];
	    for (i = 0; i < nbg; i++) {
		if (gm[magsort][i] > max_mag) max_mag = gm[magsort][i];
		if (gm[magsort][i] < min_mag) min_mag = gm[magsort][i];
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
	    case WCS_PSQUARE:
		strcpy (rstr, "SQUARE");
		break;
	    case WCS_DIAMOND:
	    case WCS_PDIAMOND:
		strcpy (rstr, "DIAMOND");
		break;
	    case WCS_CROSS:
	    case WCS_PCROSS:
		strcpy (rstr, "CROSS");
		break;
	    case WCS_EX:
	    case WCS_PEX:
		strcpy (rstr, "EX");
		break;
	    case WCS_CIRCLE:
	    case WCS_PCIRCLE:
	    default:
		strcpy (rstr, "CIRCLE");
	    }

	for (i = 0; i < nbg; i++) {
	    if (gx[i] > 0.0 && gy[i] > 0.0) {
		if (region_radius[icat] > 0) {
		    if (region_char[icat] > 10)
			radius = region_radius[icat];
		    else {
			wcs2pix (wcs, gra[i], gdec[i]+ddec, &x, &y, &offscale);
			radius = (int) (sqrt ((x-gx[i])*(x-gx[i]) +
					      (y-gy[i])*(y-gy[i])) + 0.5);
			}
		    }
		else if (rmax == 0)
		    radius = 20;
		else
		    radius = 5 + (int) (rmax * (max_mag - gm[magsort][i]) / magscale);
		ix = (int)(gx[i] + 0.5);
		iy = (int)(gy[i] + 0.5);
		printobj = 0;
		if (obname[icat] && gobj1 != NULL) {
		    if (gobj1[i] != NULL) {
			if (strlen (gobj1[i]) < 32)
			    strcpy (snum, gobj1[i]);
			else
			    strncpy (snum, gobj1[i], 31);
			printobj = 1;
			}
		    else
			CatNum (refcat, 0, -1, gnum[i], snum);
		    }
		else
		    CatNum (refcat, 0, -1, gnum[i], snum);
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

    catalog = CatName (refcat, refcatname[icat]);
    if (wfile)
	fprintf (fd, "catalog	%s\n", catalog);
    if (tabout)
	printf ("catalog	%s\n", catalog);

    if (uplate > 0) {
	sprintf (headline, "plate	%d", uplate);
	if (wfile)
            fprintf (fd, "%s\n", headline);
        if (tabout)
            printf ("%s\n", headline);
        }

    /* Minimum number of plate IDs for USNO-B1.0 catalog */
    if (minid != 0) {
	sprintf (headline, "minid	%d", minid);
	if (wfile)
	    fprintf (fd, "%s\n", headline);
	if (tabout)
	    printf ("%s\n", headline);
	}

    /* Minimum proper motion quality for USNO-B1.0 catalog */
    if (minpmqual > 0) {
	sprintf (headline, "minpmq	%d", minpmqual);
	if (wfile)
	    fprintf (fd, "%s\n", headline);
	if (tabout)
	    printf ("%s\n", headline);
	}

    if (catsort == SORT_RA) {
	sprintf (headline, "rasort	T");
	if (wfile)
	    fprintf (fd, "%s\n", headline);
	if (tabout)
	    printf ("%s\n", headline);
	}

    /* Equinox of output coordinates */
    if (sysout == WCS_J2000)
	sprintf (headline, "radecsys	FK5");
    else if (sysout == WCS_B1950)
	sprintf (headline, "radecsys	FK4");
    else if (sysout == WCS_ECLIPTIC)
	sprintf (headline, "radecsys	ecliptic");
    else if (sysout == WCS_GALACTIC)
	sprintf (headline, "radecsys	galactic");
    else
	sprintf (headline, "radecsys	unknown");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    /* Equinox of output coordinates */
    sprintf (headline, "equinox	%.2f", eqout);
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    /* Epoch of output coordinates from image header time of observation */
    sprintf (headline, "epoch	%.2f", epout);
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    /* Proper Motion units, if proper motion is in this catalog */
    if (mprop == 1) {
	if (degout) {
	    if (wfile)
		fprintf (fd, "rpmunit	arcsec/century\n");
	    else if (tabout)
		printf ("rpmunit	arcsec/century\n");
	    }
	else {
	    if (wfile)
		fprintf (fd, "rpmunit	tsec/century\n");
	    else if (tabout)
		printf ("rpmunit	tsec/century\n");
	    }
	if (wfile)
	    fprintf (fd, "dpmunit	arcsec/century\n");
	else if (tabout)
	    printf ("dpmunit	arcsec/century\n");
	}

    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (refcat == TABCAT && starcat[icat]->keyid[0] >0) {
	strcpy (headline, starcat[icat]->keyid);
	strcat (headline, "                ");
	}
    else
	CatID (headline, refcat);
    headline[nnfld] = (char) 0;
    strcat (headline, "	");
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
	strcat (headline,"magb  	magr  	");
    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==HIP || refcat==ACT)
	strcat (headline,"magb  	magv  	");
    else if (refcat ==IRAS)
	strcat (headline,"f10m  	f25m  	f60m  	f100m 	");
    else if (refcat == GSC2)
	strcat (headline,"magf  	magj 	magv	magn	");
    else if (refcat == UB1)
	strcat (headline,"magb1 	magr1	magb2	magr2	magn 	");
    else if (refcat == TMPSC)
	strcat (headline,"magj   	magh   	magk   	");
    else {
	for (imag = 0; imag < nmag; imag++) {
	    if (starcat[icat] != NULL &&
		strlen (starcat[icat]->keymag[imag]) > 0)
		sprintf (temp, "    %s ", starcat[icat]->keymag[imag]);
	    else if (nmag > 1)
		sprintf (temp, "    mag%d ", imag);
	    else
		sprintf (temp, "    mag  ");
	    strcat (headline, temp);
	    }
	}

    if (refcat == HIP)
	strcat (headline,"parlx	parer	");
    else if (refcat == IRAS)
	strcat (headline,"f60m 	f100m	");
    else if (refcat == UAC  || refcat == UA1  || refcat == UA2 ||
	refcat == USAC || refcat == USA1 || refcat == USA2 || refcat == UJC)
	strcat (headline,"plate	");
    else if (refcat == GSC || refcat == GSCACT)
	strcat (headline,"class 	band	N	");
    else if (refcat == UB1)
	strcat (headline,"pm 	ni	");
    else if (sptype == 1)
	strcat (headline,"type	");
    else if (gcset)
	strcat (headline,"peak	");
    if (mprop)
	strcat (headline, "pmra 	pmdec	");
    strcat (headline,"x    	y    ");
    if (refcat == TABCAT && keyword != NULL) {
	strcat (headline,"	");
	strcat (headline, keyword);
	}
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    strcpy (headline,"----------------------");		/* ID number */
    headline[nnfld] = (char) 0;
    strcat (headline, "	-----------	------------	-----");/* RA Dec Mag */
    if (refcat == UAC  || refcat == UA1  || refcat == UA2 || 
	refcat == USAC || refcat == USA1 || refcat == USA2 || refcat == TYCHO ||
	refcat == TYCHO2 || refcat == ACT)
	strcat (headline,"	-----");		/* Second magnitude */
    else if (refcat == TMPSC)
	strcat (headline,"--	-------	-------"); /* JHK Magnitudes */
    else if (refcat == IRAS)
	strcat (headline,"-	------	------	------"); /* 4 fluxes */
    else if (refcat == GSC2 || refcat == HIP)
	strcat (headline,"	-----	-----	-----"); /* 4 magnitudes */
    else if (refcat == UB1)
	strcat (headline,"	-----	-----	-----	-----");
    else {
	for (imag = 1; imag < nmag; imag++) {
	    sprintf (temp, "	-----");
	    strcat (headline, temp);
	    }
	}

    if (refcat == GSC || refcat == GSCACT)
	strcat (headline,"	-----	----	-");	/* class, band, n */
    else if (refcat == UB1)
	strcat (headline,"	--	--");
    else if (gcset)
	strcat (headline, "	-----");		/* plate or peak */
    if (mprop)
	strcat (headline, "	------	------");	/* Proper motion */
    strcat (headline, "	------	------");		/* X and Y */
    if (refcat == TABCAT && keyword != NULL)
	strcat (headline,"	------");		/* Additional keyword */
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (printhead) {
	if (nbg == 0)
	    printf ("No %s Found\n", title);
	else {
	    if (refcat == TABCAT && strlen(starcat[icat]->keyid) > 0)
		printf ("%s          ", starcat[icat]->keyid);
	    else
		CatID (headline, refcat);
	    headline[nnfld] = (char) 0;
	    printf ("%s", headline);

	    if (sysout == WCS_B1950) {
		if (degout) {
		    if (eqout == 1950.0)
			printf ("  RA1950   Dec1950  ");
		    else
			printf (" RAB%7.2f DecB%7.2f  ", eqout, eqout);
		    }
		else {
		    if (eqout == 1950.0)
			printf ("  RAB1950      DecB1950    ");
		    else
			printf (" RAB%7.2f   DecB%7.2f  ", eqout, eqout);
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
			printf (" RAJ%7.2f  DecJ%7.2f ", eqout, eqout);
		    }
		else {
		    if (eqout == 2000.0)
			printf ("  RA2000        Dec2000    ");
		    else
			printf (" RAJ%7.2f   DecJ%7.2f  ", eqout, eqout);
		    }
		}
	    if (refcat == UAC  || refcat == UA1  || refcat == UA2 ||
		refcat == USAC || refcat == USA1 || refcat == USA2)
		printf ("MagB  MagR Plate   X      Y   \n");
	    else if (refcat == UJC)
		printf ("  Mag Plate   X      Y   \n");
	    else if (refcat == GSC || refcat == GSCACT)
		printf ("  Mag Class Band N    X       Y   \n");
	    else if (refcat == GSC2)
		printf ("MagF  MagJ  MagV  MagN    X       Y   \n");
	    else if (refcat == UB1)
		printf ("MagB1 MagR1 MagB2 MagR2 MagN  PM NI    X       Y   \n");
	    else if (refcat == GSC2)
		printf ("MagB1 MagR1 MagB2 MagR2 MagN    X       Y   \n");
	    else if (refcat == IRAS)
		printf ("f10m  f25m  f60m  f100m   X       Y   \n");
	    else if (refcat == HIP)
		printf ("MagB  MagV  parlx parer   X       Y   \n");
	    else if (refcat == TMPSC)
		printf ("MagJ    MagH    MagK      X       Y   \n");
	    else if (refcat == SAO || refcat == PPM || refcat == BSC)
		printf ("  Mag  Type   X       Y     \n");
	    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==ACT)
		printf ("MagR   MagV    X       Y     \n");
	    else if (refcat == TABCAT) {
		for (imag = 0; imag < nmag; imag++) {
		    if (starcat[icat] != NULL &&
			strlen (starcat[icat]->keymag[imag]) > 0)
			sprintf (temp, "    %s ", starcat[icat]->keymag[imag]);
		    else if (nmag > 1)
			sprintf (temp, "    mag%d ", imag);
		    else
			sprintf (temp, "    mag  ");
		    strcat (headline, temp);
		    }
		if (gcset)
		    printf (" Peak ");
		printf (" X      Y  ");
		if (keyword != NULL)
		    printf ("   %s\n", keyword);
		}
	    else if (refcat == BINCAT)
		printf (" Mag   Type   X      Y     Object\n");
	    else if (refcat == TXTCAT)
		printf (" Mag     X      Y     Object\n");
	    else if (gcset)
		printf (" Mag  Peak     X       Y   \n");
	    else
		printf (" Mag     X       Y   \n");
	    }
	}

    /* Print positions from reference catalog */
    for (i = 0; i < nbg; i++) {
	if (gx[i] > 0.0 && gy[i] > 0.0 && gx[i] < gxmax && gy[i] < gymax) {
	    if (sptype == 1) {
	    	isp[0] = gc[i] / 1000;
		isp[1] = gc[i] % 1000;
		}
	    if (refcat == GSC || refcat == GSCACT) {
		ngsc = gc[i] / 10000;
		gc[i] = gc[i] - (ngsc * 10000);
		band = gc[i] / 100;
		gc[i] = gc[i] - (band * 100);
		}
	    CatNum (refcat, -nnfld, nndec, gnum[i], numstr);
	    if (degout) {
		deg2str (rastr, 32, gra[i], 5);
		deg2str (decstr, 32, gdec[i], 5);
		}
	    else {
		ra2str (rastr, 32, gra[i], 3);
		dec2str (decstr, 32, gdec[i], 2);
		}
	    if (tabout || wfile) {
		if (refcat == GSC || refcat == GSCACT)
		    sprintf (headline, "%s	%s	%s	%5.2f	%d	%d	%d",
		     numstr, rastr, decstr, gm[0][i], gc[i], band, ngsc);
		else if (refcat == GSC2 || refcat == HIP)
		    sprintf (headline, "%s	%s	%s	%5.2f	%5.2f	%5.2f	%5.2f",
		     numstr,rastr,decstr,gm[0][i],gm[1][i],gm[2][i],gm[3][i]);
		else if (refcat == UB1)
		    sprintf (headline, "%s	%s	%s	%5.2f	%5.2f	%5.2f	%5.2f	%5.2f	%2d	%2d",
		     numstr,rastr,decstr,gm[0][i],gm[1][i],gm[2][i],gm[3][i],
		     gm[4][i],gc[i]/100,gc[i]%100);
		else if (refcat == TMPSC) {
		    sprintf (headline, "%s	%s	%s", numstr, rastr, decstr);
		    for (imag = 0; imag < 3; imag++) {
			if (gm[imag][i] > 100.0)
			    sprintf (temp, "	%6.3fL", gm[imag][i]-100.0);
			else
			    sprintf (temp, "	%6.3f ", gm[imag][i]);
			strcat (headline, temp);
			}
		    }
		else if (refcat == IRAS) {
		    sprintf (headline, "%s	%s	%s", numstr, rastr, decstr);
		    for (imag = 0; imag < 4; imag++) {
			if (gm[imag][i] > 100.0) {
			    flux = 1000.0 * pow (10.0,-(gm[imag][i]-100.0)/2.5);
			    sprintf (temp, "	%.2fL", flux);
			    }
			else {
			    flux = 1000.0 * pow (10.0, -gm[imag][i] / 2.5);
			    sprintf (temp, "	%.2f ", flux);
			    }
			strcat (headline, temp);
			}
		    }
		else if (refcat == UAC  || refcat == UA1  || refcat == UA2 ||
			 refcat == USAC || refcat == USA1 || refcat == USA2)
		    sprintf (headline, "%s	%s	%s	%5.1f	%5.1f	%d",
		     numstr,rastr,decstr,gm[0][i],gm[1][i],gc[i]);
		else if (refcat == UJC)
		    sprintf (headline, "%s	%s	%s	%5.2f	%d",
		     numstr, rastr, decstr, gm[0][i], gc[i]);
		else if (refcat==SAO || refcat==PPM || refcat== BSC ) {
		    sprintf (headline, "%s	%s	%s	%5.2f	%2s",
		     numstr,rastr,decstr,gm[0][i],isp);
		    }
		else if (refcat==TYCHO || refcat==TYCHO2 || refcat==ACT) {
		    sprintf (headline, "%s	%s	%s	%5.2f	%5.2f",
		     numstr,rastr,decstr,gm[0][i],gm[1][i],isp);
		    }
		else {
		    sprintf (headline, "%s	%s	%s",
			     numstr, rastr, decstr);
		    for (imag = 0; imag < nmag; imag++) {
			sprintf (temp, "	%5.2f",gm[imag][i]);
			strcat (headline, temp);
			}
		    if (gcset) {
			sprintf (temp, "	%d", gc[i]);
			strcat (headline, temp);
			}
		    if (sptype == 1) {
			sprintf (temp, "	%s", isp);
			strcat (headline, temp);
			}
		    }
		if (mprop) {
		    if (degout)
			pra = gpra[i] * 360000.0;
		    else
			pra = gpra[i] * 24000.0;
		    pdec = gpdec[i] * 360000.0;
		    sprintf (temp, "	%.2f	%.2f", pra,pdec);
		    strcat (headline, temp);
		    }
		sprintf (temp, "	%.2f	%.2f",
			 gx[i],gy[i]);
		strcat (headline, temp);
		if (wfile)
		    fprintf (fd, "%s\n", headline);
		if (tabout)
		    printf ("%s\n", headline);
		}
	    else if (!tabout) {
		if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
		    refcat == UAC  || refcat == UA1  || refcat == UA2)
		    sprintf (headline,"%s %s %s %5.1f %5.1f %4d",
			numstr,rastr,decstr,gm[0][i],gm[1][i],gc[i]);
		else if (refcat == UJC)
		    sprintf (headline,"%s %s %s %6.2f %4d",
			numstr, rastr, decstr, gm[0][i], gc[i]);
		else if (refcat == GSC || refcat == GSCACT)
		    sprintf (headline,"%s %s %s %6.2f %4d %4d %2d",
			numstr, rastr, decstr, gm[0][i], gc[i], band, ngsc);
		else if (refcat == TMPSC) {
		    sprintf (headline, "%s %s %s", numstr, rastr, decstr);
		    for (imag = 0; imag < 3; imag++) {
			if (gm[imag][i] > 100.0)
			    sprintf (temp, " %6.3fL", gm[imag][i]-100.0);
			else
			    sprintf (temp, " %6.3f ", gm[imag][i]);
			strcat (headline, temp);
			}
		    }
		else if (refcat == GSC2 || refcat == HIP)
		    sprintf (headline,"%s %s %s %5.2f %5.2f %5.2f %5.2f",
			     numstr,rastr,decstr,gm[0][i],gm[1][i],gm[2][i],gm[3][i]);
		else if (refcat == UB1)
		    sprintf (headline,"%s %s %s %5.2f %5.2f %5.2f %5.2f %5.2f %2d %2d",
			     numstr,rastr,decstr,gm[0][i],gm[1][i],gm[2][i],
			     gm[3][i],gm[4][i],gc[i]/100,gc[i]%100);
		else if (refcat==IRAS) {
		    sprintf (headline, "%s %s %s", numstr, rastr, decstr);
		    for (imag = 0; imag < 3; imag++) {
			if (gm[imag][i] > 100.0) {
			    flux = 1000.0 * pow (10.0, -(gm[imag][i]-100.0) / 2.5);
			    sprintf (temp, " %5.2fL", flux);
			    }
			else {
			    flux = 1000.0 * pow (10.0, -gm[imag][i] / 2.5);
			    sprintf (temp, " %5.2f ", flux);
			    }
			strcat (headline, temp);
			}
		    }
		else if (refcat==SAO || refcat==PPM || refcat== BSC)
		    sprintf (headline,"%s  %s %s %6.2f  %2s",
			numstr,rastr,decstr,gm[0][i],isp);
		else if (refcat==TYCHO || refcat==TYCHO2 || refcat==ACT)
		    sprintf (headline,"%s %s %s %6.2f %6.2f",
			     numstr,rastr,decstr,gm[0][i],gm[1][i],isp);
		else {
		    sprintf (headline,"%s %s %s",
			     numstr,rastr,decstr);
		    for (imag = 0; imag < nmag; imag++) {
			sprintf (temp, " %5.2f",gm[imag][i]);
			strcat (headline, temp);
			}
		    if (sptype == 1) {
			sprintf (temp, " %2s", isp);
			strcat (headline, temp);
			}
		    if (refcat == TABCAT && gcset) {
			sprintf(temp,"	%d", gc);
			strcat (headline, temp);
			}
		    }

		/* Add image pixel coordinates to output line */
		if (wcs->nxpix < 1000.0 && wcs->nypix < 1000.0)
		    sprintf (temp, " %6.2f %6.2f", gx[i], gy[i]);
		else if (wcs->nxpix < 10000.0 && wcs->nypix < 10000.0)
		    sprintf (temp, " %6.1f %6.1f", gx[i], gy[i]);
		else
		    sprintf (temp, " %.1f %.1f", gx[i], gy[i]);
		strcat (headline, temp);

		/* Add object name to output line */
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
    if (gm) {
	for (imag = 0; i < nmagmax; i++)
	    free ((char *)gm[imag]);
	free ((char *)gm);
	}
    if (gra) free ((char *)gra);
    if (gdec) free ((char *)gdec);
    if (gnum) free ((char *)gnum);
    if (gc) free ((char *)gc);
    if (gobj) free ((char *)gobj);
    wcsfree (wcs);

    return;
}

/* Set up limits for search */
static void
ImageLim (wcs, cra, cdec, dra, ddec, ramin, ramax, decmin, decmax)

struct WorldCoor *wcs;		/* WCS parameter structure */
double	*cra, *cdec;		/* Center of search area  in degrees (returned) */
double	*dra, *ddec;		/* Horizontal and vertical half-widths in degrees (returned) */
double	*ramin, *ramax;		/* Right ascension limits in degrees (returned) */
double	*decmin, *decmax;	/* Declination limits in degrees (returned) */

{
    double xmin = 0.5;
    double xmax = wcs->nxpix + 0.5;
    double xcen = 0.5 + (wcs->nxpix * 0.5);
    double ymin = 0.5;
    double ymax = wcs->nypix + 0.5;
    double ycen = 0.5 + (wcs->nypix * 0.5);
    double ra[8], dec[8];
    int i;

    /* Find sky coordinates of corners and middles of sides */
    pix2wcs (wcs, xmin, ymin, &ra[0], &dec[0]);
    pix2wcs (wcs, xmin, ycen, &ra[1], &dec[1]);
    pix2wcs (wcs, xmin, ymax, &ra[2], &dec[2]);
    pix2wcs (wcs, xcen, ymin, &ra[3], &dec[3]);
    pix2wcs (wcs, xcen, ymax, &ra[4], &dec[4]);
    pix2wcs (wcs, xmax, ymin, &ra[5], &dec[5]);
    pix2wcs (wcs, xmax, ycen, &ra[6], &dec[6]);
    pix2wcs (wcs, xmax, ymax, &ra[7], &dec[7]);

    /* Find minimum and maximum right ascensions and declinations */
    *ramin = ra[0];
    *ramax = ra[0];
    *decmin = dec[0];
    *decmax = dec[0];
    for (i = 0; i < 8; i++) {
	if (ra[i] < *ramin)
	   *ramin = ra[i];
	if (ra[i] > *ramax)
	   *ramax = ra[i];
	if (dec[i] < *decmin)
	   *decmin = dec[i];
	if (dec[i] > *decmax)
	   *decmax = dec[i];
	}

    /* Set center and extent */
    *cra = 0.5 * (*ramin + *ramax);
    *cdec = 0.5 * (*decmin + *decmax);
    if (*ramax - *ramin > 180.0)
	*dra = 360.0 + *ramin  - *ramax;
    else
	*dra = 0.5 * (*ramax - *ramin);
    *ddec = 0.5 * (*decmax - *decmin);

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
 * Feb 21 1997	Get image header from GetFITSWCS()
 * Mar 14 1997	Add support for USNO SA-1.0 catalog
 * Apr 25 1997	Fix bug in uacread
 * May 28 1997	Add option to read a list of filenames from a file
 * Nov 17 1997	Initialize both magnitude limits
 * Dec  8 1997	Set up program to be called by various names
 * Dec 15 1997	Add capability of reading and writing IRAF 2.11 images
 *
 * Jan 27 1998	Implement Mark Calabretta's WCSLIB
 * Jan 27 1998	Add -z for AIPS classic WCS projections
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
 * Jun 17 1999	Use SearchLim() to compute search limits
 * Jul  8 1999	Fix bug when noll object name list
 * Aug 24 1999	If radius not set for region mode, use 20 pixels
 * Aug 25 1999	Add Bright Star Catalog, BSC
 * Aug 25 1999	Allocate using calloc(ngmax, ) instead of malloc(nbytes)
 * Aug 25 1999	Add option to set circle radius in pixels if < -1
 * Sep 10 1999	Do all searches through catread() and catrnum()
 * Sep 16 1999	Add zero distsort argument to catread() call
 * Oct 15 1999	Free wcs using wcsfree()
 * Oct 22 1999	Drop unused variables after lint
 * Oct 22 1999	Change catread() to ctgread() to avoid system conflict
 *
 * Jan 11 2000	Get nndec for Starbase catalogs
 * Jan 28 2000	Call setdefwcs() with WCS_ALT instead of 1
 * Feb 15 2000	Use MAXCAT from lwcs.h instead of MAXREF
 * Feb 28 2000	Drop Peak column if not set in TAB catalogs
 * Mar 10 2000	Move catalog selection from executable name to subroutine
 * Mar 15 2000	Add proper motions to catalog calls and Starbase output
 * Mar 28 2000	Clean up output for catalog IDs and GSC classes
 * May 26 2000	Add Tycho 2 catalog
 * May 26 2000	Always use CatNumLen() to get ID number field size
 * Jul 12 2000	Add star catalog data structure to ctgread() argument list
 * Jul 25 2000	Add coordinate system to SearchLim() call
 * Jul 25 2000	Fix star catalog structure initialization bug
 * Jul 25 2000	Pass address of star catalog data structure address
 * Sep 21 2000	Print spectral type instead of plate number of USNO-A catalogs
 * Dec  1 2000	Print plate, not type for USNO catalogs
 * Dec 15 2000	Deal with missing catalog names
 * Dec 18 2000	Always allocate proper motion arrays
 *
 * Jan 22 2001	Drop declaration of wcsinit()
 * Feb 23 2001	Drop distance from search center output everywhere
 * Mar 27 2001	Add option to set size of overplotted stars in pixels
 * May 23 2001	Add support for GSC-ACT catalog
 * May 24 2001	Add support for 2MASS Point Source Catalog
 * Jun  8 2001	Add proper motion flag and number of magnitudes to RefCat()
 * Jun 29 2001	Add support for GSC II catalog
 * Sep 13 2001	Allow sort by RA, Dec, X, Y, any magnitude
 * Sep 13 2001	Use 2-D array of magnitudes, rather than multiple vectors
 * Sep 14 2001	Add option to print catalog as returned over the web
 * Sep 18 2001	Add flags to IRAS and 2MASS point sources
 * Sep 20 2001	Get magnitude name for limits from CatMagName()
 * Oct 19 2001	Add -y to set epoch of observation
 * Oct 25 2001	Allow arbitrary argument order on command line
 * Oct 31 2001	Print complete help message if no arguments
 * Nov  6 2001	Declare undeclared subroutine setparm()
 *
 * Feb  1 2002	Print spectral type for TDC format catalogs, if present
 * Apr  3 2002	Add magnitude number to sort options
 * Apr  8 2002	Add magnitude number to magnitude limit setting
 * Apr  8 2002	Fix bug so that characters other than circles can be plotted
 * Apr 10 2002	Fix magnitude number bug and add magnitude letter
 * May  1 2002	Add -a command to set initial rotation angle
 * Jun 19 2002	Add verbose argument to GetFITShead()
 * Aug  6 2002	Print all magnitudes for BINARY, TABTABLE, or ASCII catalogs
 *
 * Jan 26 2003	Add support for USNO-B1.0 catalog
 * Jan 28 2003	Fix bug printing proper motion
 * Jan 29 2003	Add header lines if USNO-B1.0 ID or PM quality limits
 * Mar  4 2003	If star is offscale, set x and y to 0.0
 * Mar 25 2003	Deal correctly with rotated images
 * Apr  2 2003	Try rotated images again
 */
