/* File imcat.c
 * July 9, 1998
 * By Doug Mink, after Elwood Downey
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
#include "libwcs/lwcs.h"

#define GSC	1	/* refcat value for HST Guide Star Catalog */
#define UJC	2	/* refcat value for USNO UJ Star Catalog */
#define UAC	3	/* refcat value for USNO A-1.0 Star Catalog */
#define USAC	4	/* refcat value for USNO SA-1.0 Star Catalog */

extern int gscread();
extern int ujcread();
extern int uacread();
extern int usaread();
extern int tabread();
static void usage();
static void ListCat();
extern void RASortStars();
extern void MagSortStars();
extern void fk524e();
extern struct WorldCoor *GetFITSWCS();
extern char *GetFITShead();
extern void setsys();
extern void setcenter();
extern void setsecpix();


static int verbose = 0;		/* verbose/debugging flag */
static int wfile = 0;		/* True to print output file */
static int refcat = GSC;	/* reference catalog switch */
static char refcatname[32]="GSC"; /* reference catalog name */
static int classd = -1;		/* Guide Star Catalog object classes */
static int uplate = 0;		/* UJ Catalog plate number to use */
static double maglim1 = MAGLIM1; /* reference catalog bright magnitude limit */
static double maglim2 = MAGLIM2; /* reference catalog faint magnitude limit */
static char coorsys[4];		/* Output coordinate system */
static int nstars = 0;		/* Number of brightest stars to list */
static int printhead = 0;	/* 1 to print table heading */
static int tabout = 0;		/* 1 for tab table to standard output */
static int rasort = 0;		/* 1 to sort stars by brighness */
static int debug = 0;		/* True for extra information */
static int degout0 = 0;		/* True for RA and Dec in fractional degrees */
static char *progname;		/* Name of program as executed */
static int sysout = -1;		/* Output coordinate system */
static int sysref = WCS_J2000;	/* Output coordinate system */
static double eqout = 2000.0;	/* Equinox for output coordinates */
static double eqref = 2000.0;	/* Equinox of catalog to be searched */
static char progpath[128];

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
    int i, ic;

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
		strcpy (coorsys, "FK4");
		sysout = WCS_B1950;
		eqout = 1950.0;
		}
	    else if (ac < 3)
		usage ();
	    else {
		setsys(WCS_B1950);
		strcpy (coorsys, "FK4");
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
		usage();
	    strcpy (refcatname, *++av);
	    if (strncmp(refcatname,"gs",2)==0 ||
		strncmp (refcatname,"GS",2)== 0)
		refcat = GSC;
	    else if (strncmp(refcatname,"ua",2)==0 ||
		     strncmp(refcatname,"UA",2)==0)
		refcat = UAC;
	    else if (strncmp(refcatname,"us",2)==0 ||
		     strncmp(refcatname,"US",2)==0)
		refcat = USAC;
	    else if (strncmp(refcatname,"uj",2)==0 ||
		     strncmp(refcatname,"UJ",2)==0)
		refcat = UJC;
	    else
		refcat = 0;
	    ac--;
	    break;

	case 'd':
	    degout0++;
	    break;

	case 'e':
	    sysout = WCS_ECLIPTIC;
	    strcpy (coorsys, "ecliptic");
	    break;

	case 'g':
	    sysout = WCS_GALACTIC;
	    strcpy (coorsys, "galactic");
	    break;

	case 'h':	/* ouput descriptive header */
	    printhead++;
	    break;

    	case 'j':	/* center coordinates on command line in J2000 */
	    str1 = *(av+1);
	    ic = (int)str1[0];
	    if (*(str+1) || ic < 48 || ic > 58) {
		strcpy (coorsys, "FK5");
		setsys(WCS_J2000);
		sysout = WCS_J2000;
		eqout = 2000.0;
		}
	    else if (ac < 3)
		usage ();
	    else {
		setsys(WCS_J2000);
		sysout = WCS_J2000;
		eqout = 2000.0;
		strcpy (coorsys, "FK5");
		strcpy (rastr, *++av);
		ac--;
		strcpy (decstr, *++av);
		ac--;
		setcenter (rastr, decstr);
		}
    	    break;

    	case 'm':	/* Limiting reference star magnitude */
    	    if (ac < 2)
    		usage();
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
		usage ();
	    nstars = atoi (*++av);
	    ac--;
	    break;

    	case 'p':	/* Initial plate scale in arcseconds per pixel */
    	    if (ac < 2)
    		usage();
    	    setsecpix (atof (*++av));
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
    		usage();
    	    uplate = (int) atof (*++av);
    	    ac--;
    	    break;

    	case 'v':	/* more verbosity */
    	    verbose++;
    	    break;

    	case 'w':	/* write output file */
    	    wfile++;
    	    break;

	case 'x':	/* Guide Star object class */
    	    if (ac < 2)
    		usage();
    	    classd = (int) atof (*++av);
    	    ac--;
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
    	    usage();
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
	    usage ();
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
    	    if (debug)
    		printf ("%s:\n", filename);
	    ListCat (filename);
	    }
	fclose (flist);
	}

    /* If no arguments left, print usage */
    if (ac == 0 || refcat < 0)
	usage();

    while (ac-- > 0) {
    	char *fn = *av++;
    	if (debug)
    	    printf ("%s:\n", fn);
    	ListCat (fn);
    	if (debug)
    	    printf ("\n");
	}

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"\n");
    if (refcat == GSC) {
	fprintf (stderr,"List HST Guide Stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-g class]\n");
	}
    else if (refcat == UJC) {
	fprintf (stderr,"List USNO J Catalog stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-u plate]\n");
	}
    else if (refcat == UAC) {
	fprintf (stderr,"List USNO A 1.0 stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-u plate]\n");
	}
    else if (refcat == USAC) {
	fprintf (stderr,"List USNO SA 1.0 stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vhst] [-m [mag1] mag2] [-u plate]\n");
	}
    else {
	fprintf (stderr,"List catalog stars in FITS and IRAF image files\n");
	fprintf (stderr,"Usage: [-vwhst] [-m [mag1] mag2] [-c catalog]\n");
	}
    fprintf (stderr,"       [-p scale] [-b ra dec] [-j ra dec] FITS or IRAF file(s)\n");
    fprintf (stderr,"  -b: Output, (center) in B1950 (FK4) RA and Dec\n");
    fprintf (stderr,"  -c: Reference catalog (gsc, ujc, or tab table file\n");
    fprintf (stderr,"  -d: Output RA,Dec positions in fractional degrees\n");
    fprintf (stderr,"  -e: Output in ecliptic longitude and latitude\n");
    fprintf (stderr,"  -g: Output in galactic longitude and latitude\n");
    fprintf (stderr,"  -h: Print heading, else do not \n");
    fprintf (stderr,"  -j: Output (center) in J2000 (FK5) RA and Dec\n");
    fprintf (stderr,"  -m: Limiting catalog magnitude(s) (default none)\n");
    fprintf (stderr,"  -n: Number of brightest stars to print \n");
    fprintf (stderr,"  -p: Initial plate scale in arcsec per pixel (default 0)\n");
    fprintf (stderr,"  -s: Sort by RA instead of flux \n");
    fprintf (stderr,"  -t: Tab table to standard output as well as file\n");
    fprintf (stderr,"  -u: USNO catalog single plate number to accept\n");
    fprintf (stderr,"  -w: Write tab table output file imagename.cat\n");
    fprintf (stderr,"  -v: Verbose\n");
    fprintf (stderr,"  -x: Guide Star Catalog class (-1=all,0,3 (default -1)\n");
    fprintf (stderr,"  -z: Use AIPS classic projections instead of WCSLIB\n");
    exit (1);
}


struct WorldCoor *wcsinit();	
extern int pix2wcst();

static void
ListCat (filename)

char	*filename;	/* FITS or IRAF file filename */

{
    char *header;	/* FITS image header */
    double *gnum=0;	/* Catalog numbers */
    double *gra=0;	/* Catalog right ascensions, rads */
    double *gdec=0;	/* Catalog declinations rads */
    double *gm=0;	/* Catalog star magnitudes */
    double *gmb=0;	/* Catalog star blue magnitudes */
    double *gx, *gy;	/* Catalog star positions on image */
    int *gc;		/* Catalog object classes, plates, etc. */
    int ng;		/* Number of catalog stars */
    int nbg;		/* Number of brightest catalog stars actually used */
    int imh, imw;	/* Image height and width in pixels */
    int i, ngmax, nbytes;
    int degout;
    FILE *fd;
    struct WorldCoor *wcs;	/* World coordinate system structure */
    char rastr[16], decstr[16];	/* coordinate strings */
    double cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, mag1, mag2,secpix;
    double mag, gdist;
    int offscale, nlog;
    char headline[160];
    char title[80];
    char outfile[80];

    if (verbose || printhead)

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
	fprintf (stderr,"ListCat: No tab table catalog name\n");
	return;
	}

    /* Read world coordinate system information from the image header */
    if ((header = GetFITShead (filename)) == NULL)
	return;
    wcs = GetFITSWCS (header, verbose, &cra, &cdec, &dra, &ddec, &secpix,
		      &imw, &imh, eqref);
    free (header);
    if (nowcs (wcs))
	return;

    if (sysout < 0) {
        sysout = wcs->syswcs;
        eqout = wcs->equinox;
        }

    /* Set up limits for search */
    ra1 = cra - dra;
    ra2 = cra + dra;
    dec1 = cdec - ddec;
    dec2 = cdec + ddec;
    if (verbose || printhead) {
	char rastr1[16],rastr2[16],decstr1[16],decstr2[16];
	ra2str (rastr1, 16, ra1, 3);
	ra2str (rastr2, 16, ra2, 3);
	printf ("%s: RA:  %s - %s (J2000)\n",filename,rastr1,rastr2);
	dec2str (decstr1, 16, dec1, 2);
	dec2str (decstr2, 16, dec2, 2);
	printf ("%s: Dec: %s - %s (J2000)\n",filename, decstr1,decstr2);
	}

/* Set the output coordinate system */
    wcs->printsys = 0;
    if (coorsys[1])
	wcsoutinit (wcs, coorsys);

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
    if (!gnum || !gra || !gdec || !gm || !gc) {
	if (gm) free ((char *)gm);
	if (gmb) free ((char *)gmb);
	if (gra) free ((char *)gra);
	if (gdec) free ((char *)gdec);
	if (gnum) free ((char *)gnum);
	if (gc) free ((char *)gc);
	return;
	}
    if (debug)
	nlog = 100;
    else
	nlog = 0;

    /* Find the nearby reference stars, in ra/dec */
    if (refcat == GSC)
	ng = gscread (cra,cdec,dra,ddec,0.0,mag1,mag2,classd,ngmax,
		      gnum,gra,gdec,gm,gc,nlog);
    else if (refcat == UAC)
	ng = uacread (cra,cdec,dra,ddec,0.0,mag1,mag2,uplate,ngmax,
		      gnum,gra,gdec,gm,gmb,gc,debug);
    else if (refcat == USAC)
	ng = usaread (cra,cdec,dra,ddec,0.0,mag1,mag2,uplate,ngmax,
		      gnum,gra,gdec,gm,gmb,gc,debug);
    else if (refcat == UJC)
	ng = ujcread (cra,cdec,dra,ddec,0.0,mag1,mag2,uplate,ngmax,
		      gnum,gra,gdec,gm,gc,debug);
    else
	ng = tabread (refcatname,cra,cdec,dra,ddec,0.0,mag1,mag2,ngmax,
		      gnum,gra,gdec,gm,gc,debug);

    /* Project the reference stars into pixels on a plane at ra0/dec0 */
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
	free ((char *)wcs);
	return;
	}

    /* Get image pixel coordinates for each star found in reference catalog */
    for (i = 0; i < ng; i++ ) {
	offscale = 0;
	wcsc2pix (wcs, gra[i], gdec[i], "J2000", &gx[i], &gy[i], &offscale);
	}

    /* Sort reference stars by brightness (magnitude) */
    MagSortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, ng);

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
	if (strsrch (filename,".imh") != NULL)
	    printf (" in IRAF image %s\n",filename);
	else
	    printf (" in FITS image %s\n", filename);
	}

    /* Sort star-like objects in image by right ascension */
    if (rasort)
	RASortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, nbg);
    sprintf (headline, "IMAGE	%s", filename);

    /* Open plate catalog file */
    if (wfile) {
	strcpy (outfile,filename);
	if (refcat == GSC)
	    strcat (outfile,".gsc");
	else if (refcat == UAC)
	    strcat (outfile,".uac");
	else if (refcat == USAC)
	    strcat (outfile,".usac");
	else if (refcat == UJC)
	    strcat (outfile,".ujc");
	else {
	    strcat (outfile,".");
	    strcat (outfile,"refcatname");
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
	    free ((char *)wcs);
            return;
	    }
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
	sprintf (headline, "CATALOG     USNO A 1.0");
    else if (refcat == USAC)
	sprintf (headline, "CATALOG     USNO SA 1.0");
    else if (refcat == UJC)
	sprintf (headline, "CATALOG     USNO UJ1.0");
    else
	sprintf (headline, "CATALOG     %s", refcatname);
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
    else if (refcat == USAC)
	strcpy (headline,"usac_id    	");
    else if (refcat == UJC)
	strcpy (headline,"ujc_id    	");
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
    if (refcat == UAC || refcat == USAC)
	strcat (headline,"magb  	magr  	x    	y    	");
    else
	strcat (headline,"mag   	x    	y    	");
    if (refcat == GSC)
	strcat (headline,"type	");
    else if (refcat == UAC || refcat == USAC || refcat == UJC)
	strcat (headline,"plate	");
    else
	strcat (headline,"peak	");
    strcat (headline,"dist");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (refcat == UAC)
	sprintf(headline,"----------	--------	---------	----	-----	-----	-----	-----	----");
    else
        sprintf (headline,"----------	------------	------------	------	----	-------	----");
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
		printf ("MagB  MagR    X      Y   Plate Arcsec\n");
	    else if (refcat == GSC)
		printf (" Mag     X      Y   Class Arcsec\n");
	    else
		printf (" Mag     X      Y   Plate Arcsec\n");
	    }
	}

    /* Print positions from reference catalog */
    for (i = 0; i < nbg; i++) {
	if (gx[i] > 0.0 && gy[i] > 0.0) {
	    gdist = 3600.0 * wcsdist (cra, cdec, gra[i], gdec[i]);
	    wcscon (sysref, sysout, eqref, eqout, &gra[i],&gdec[i], wcs->epoch);
	    if (degout) {
		deg2str (rastr, 32, gra[i], 5);
		deg2str (decstr, 32, gdec[i], 5);
		}
	    else {
		ra2str (rastr, 32, gra[i], 3);
		dec2str (decstr, 32, gdec[i], 2);
		}
	    if (refcat == GSC)
		sprintf (headline, "%9.4f	%s	%s	%.2f	%.1f	%.1f	%d	%.2f",
		 gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i], gdist);
	    else if (refcat == UAC || refcat == USAC)
		sprintf (headline, "%13.8f	%s	%s	%.1f	%.1f	%.1f	%.1f	%d	%.2f",
		 gnum[i],rastr,decstr,gmb[i],gm[i],gx[i],gy[i],gc[i],gdist);
	    else if (refcat == UJC)
		sprintf (headline, "%12.7f	%s	%s	%.2f	%.1f	%.1f	%d	%.2f",
		 gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i], gdist);
	    else
		sprintf (headline, "%9.4f	%s	%s	%.2f	%.1f	%.1f	%d	%.2f",
		 gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i], gdist);
	    if (wfile)
		fprintf (fd, "%s\n", headline);
	    else if (tabout)
		printf ("%s\n", headline);
	    else if (refcat == UAC || refcat == USAC)
		printf ("%13.8f %s %s %5.1f %5.1f %6.1f %6.1f %4d %7.2f\n",
			gnum[i],rastr,decstr,gmb[i],gm[i],gx[i],gy[i],gc[i],gdist);
	    else if (refcat == UJC)
		printf ("%12.7f %s %s %6.2f %6.1f %6.1f %4d %7.2f\n",
			gnum[i], rastr, decstr, gm[i],gx[i],gy[i],gc[i],gdist);
	    else if (refcat == GSC)
		printf ("%9.4f %s %s %6.2f %6.1f %6.1f  %3d %7.2f\n",
			gnum[i], rastr, decstr, gm[i],gx[i],gy[i],gc[i],gdist);
	    else
		printf ("%9.4f %s %s %6.2f %6.1f %6.1f  %d %7.2f\n",
			gnum[i], rastr, decstr, gm[i],gx[i],gy[i],gc[i],gdist);
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
 */
