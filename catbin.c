/* File catbin.c
 * September 10, 1999
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

static void usage();

extern int catopen();
extern int catstar();
extern int binwopen();
extern int binwstar();
static int BinCat ();

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
static int syscoor = -1;	/* Input search coordinate system */
static int nstars = 0;		/* Number of brightest stars to list */
static int debug = 0;		/* True for extra information */
static int refcat2 = 0;		/* Second reference catalog switch */
static char *keyword = NULL;	/* Column to add to tab table output */
static int bincat = 0;		/* Flag for specific binary catalog */
static int version = 0;		/* If 1, print only program name and version */

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
    char *refcatname;	/* Input catalog name */

    refcatname = NULL;
    if (ac == 1)
        usage ();

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    /* crack arguments */
    for (av++; --ac > 0; av++) {
	if (*(str = *av) == '-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

	    case 'v':	/* more verbosity */
		verbose++;
		break;

    	    case 'b':	/* output catalog in B1950 coordinates */
		sysout = WCS_B1950;
		eqout = 1950.0;
		epout = 1950.0;
    		break;

	    case 'c':       /* Set input catalog */
		if (ac < 2)
		    usage();
		refcatname = *++av;
		refcat = -1;
		ac--;
		break;

    	    case 'j':	/* output catalog in J2000 coordinates */
		sysout = WCS_J2000;
		eqout = 2000.0;
		epout = 2000.0;
    		break;

	    case 'o':	/* Set number of characters in binary object name */
		if (ac < 2)
		    usage();
		ncobj = atoi (*++av);
		break;

	    default:
		if (refcatname) {
		    outname = *++av;
		    BinCat (refcatname);
		    }
		else
		    refcatname = *++av;
		break;
	    }
	    }
	}

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Usage: [-bjtvw] [-m [mag1] mag2] [-e sys] [-n num] [-r arcsec] [-b][-j] ra dec\n");
    fprintf(stderr,"  -b: Output B1950 (FK4) coordinates\n");
    fprintf(stderr,"  -c: Input catalog\n");
    fprintf(stderr,"  -j: Output J2000 (FK5) coordinates\n");
    fprintf(stderr,"  -o: Number of characters in output object name \n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}

#define TABMAX 64

static int
BinCat (refcatname)

char *refcatname;	/* Name of catalog to translate */

{
    double ra, dec, rapm, decpm;
    double epout;
    int sysref;		/* Coordinate system of reference catalog */
    double eqref;	/* Equinox of catalog to be searched */
    double epref;	/* Epoch of catalog to be searched */
    int i, ngmax, nbytes;
    char rastr[32], decstr[32];	/* coordinate strings */
    double drad, dra, ddec, ra1, dec1, mag1, mag2;
    double gdist;
    char **gobj;	/* Catalog star object names */
    char **gobj1;	/* Catalog star object names */
    double mag;
    int nlog;
    char headline[160];
    char filename[80];
    char title[80];
    char *tabline;
    char string[TABMAX];
    char temp[80];
    int ntab;
    int tabcat = 0;
    int refcat = 0;	/* reference catalog switch */

    if (verbose || printhead) {
	}
    tabline = 0;
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

    eqref = 2000.0;
    epref = 2000.0;
    sysref = WCS_J2000;
    if (refcat == GSC)
	strcpy (title, "HST Guide Stars");
    else if (refcat == UAC)
	strcpy (title, "USNO A Catalog Stars");
    else if (refcat == USAC)
	strcpy (title, "USNO SA Catalog Stars");
    else if (refcat == UJC)
	strcpy (title, "USNO J Catalog Stars");
    else if (refcat == SAO) {
	strcpy (title, "SAO Catalog Stars");
	sysref = WCS_B1950;
	eqref = 1950.0;
	epref = 1950.0;
	if (closest) rad0 = 3600.0;
	}
    else if (refcat == PPM) {
	strcpy (title, "PPM Catalog Stars");
	sysref = WCS_B1950;
	eqref = 1950.0;
	epref = 1950.0;
	if (closest) rad0 = 3600.0;
	}
    else if (refcat == IRAS) {
	strcpy (title, "IRAS Point Sources");
	sysref = WCS_B1950;
	eqref = 1950.0;
	epref = 1950.0;
	if (closest) rad0 = 3600.0;
	}
    else if (refcat == TYCHO) {
	strcpy (title, "Tycho Catalog Stars");
	if (closest) rad0 = 1200.0;
	}
    else if (refcatname[0] > 0) {
	if (!bincat && istab (refcatname))
	    tabcat = 1;
	sprintf (title, "%s Catalog Stars", refcatname);
	}
    else {
	fprintf (stderr,"ListCat: Catalog name is missing\n");
	return (0);
	}
    epout = epoch0;
    if (epout == 0.0)
	epout = epref;

    /* Find stars specified by number */
    if (snum != NULL) {
	nbytes = nfind * sizeof (double);
	if (!(gnum = (double *) calloc (nfind, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gnum\n",
		     nfind*sizeof(double));
	else {
	    for (i = 0; i < nfind; i++)
		gnum[i] = snum[i];
	    }
	if (!(gra = (double *) calloc (nfind, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gra\n",
		     nfind*sizeof(double));
	if (!(gdec = (double *) calloc (nfind, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gdec\n",
		     nfind*sizeof(double));
	if (!(gm = (double *) calloc (nfind, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gm\n",
		     nfind*sizeof(double));
	if (!(gmb = (double *) calloc (nfind, sizeof (double))))
	    fprintf (stderr, "Could not calloc %d bytes for gmb\n",
		     nfind*sizeof(double));
	if (!(gc = (int *) calloc (nfind, sizeof(int))))
	    fprintf (stderr, "Could not calloc %d bytes for gc\n",
		     nfind*sizeof(int));
	if (!(gx = (double *) calloc (nfind, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gx\n",
		     nfind*sizeof(double));
	if (!(gy = (double *) calloc (nfind, sizeof(double))))
	    fprintf (stderr, "Could not calloc %d bytes for gy\n",
		     nfind*sizeof(double));
	if (!(gobj = (char **) calloc (nfind, sizeof (void *))))
	    fprintf (stderr, "Could not calloc %d bytes for obj\n", nbytes);
		     nfind*sizeof(char *));
	if (!gnum || !gra || !gdec || !gm || !gmb || !gc || !gobj) {
	    if (gm) free ((char *)gm);
	    if (gmb) free ((char *)gmb);
	    if (gra) free ((char *)gra);
	    if (gdec) free ((char *)gdec);
	    if (gnum) free ((char *)gnum);
	    if (gc) free ((char *)gc);
	    if (gobj) free ((char *)gobj);
	    return (0);
	    }
	wfile = 0;

	/* Find the specified catalog stars */
	nbg = catrnum (refcatname, refcat, nfind,sysout,eqout,epout,
		       gnum,gra,gdec,gm,gmb,gc,gobj,nlog);

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
	    SearchHead (sysout,sysout,eqout,eqout,epout,epout,cra,cdec,
			dra,ddec,drad);
	    if (sysout != syscoor)
		SearchHead (sysout,syscoor,eqout,eqcoor,epout,epout,cra,cdec,
			    dra,ddec,drad);
	    if (sysref != syscoor && sysref != sysout)
		SearchHead (sysout,sysref,eqout,eqref,epout,epref,cra,cdec,
			    dra,ddec,drad);
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
	if (!(gobj = (char **) calloc (nfind, sizeof (void *))))
		     ngmax*sizeof(char *));
	    fprintf (stderr, "Could not calloc %d bytes for gobj\n",
	if (!gnum || !gra || !gdec || !gm || !gmb || !gc || !gobj) {
	    if (gm) free ((char *)gm);
	    if (gmb) free ((char *)gmb);
	    if (gra) free ((char *)gra);
	    if (gdec) free ((char *)gdec);
	    if (gnum) free ((char *)gnum);
	    if (gc) free ((char *)gc);
	    if (gobj) free ((char *)gobj);
	    return (0);
	    }

	/* Find the nearby reference stars, in ra/dec */
	ng = catread (refcatname, refcat
		      cra,cdec,dra,ddec,drad,sysout,eqout,epout,
		      mag1,mag2,ngmax,gnum,gra,gdec,gm,gmb,gc,gobj,nlog);

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
	    RASortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, ng);

	/* Sort reference stars from brightest to faintest */
	else
	    MagSortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, ng);

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
		    printf ("%d %s between %.1f and %.1f",
			    ng, title, maglim1, maglim2);
		else if (maglim2 > 0.0)
		    printf ("%d %s brighter than %.1f",
			    ng, title, maglim2);
		else if (verbose)
		    printf ("%d %s", ng, title);
		printf ("\n");
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
	    else if (refcat == SAO)
		strcat (filename,".sao");
	    else if (refcat == PPM)
		strcat (filename,".ppm");
	    else if (refcat == IRAS)
		strcat (filename,".iras");
	    else if (refcat == TYCHO)
		strcat (filename,".tycho");
	    else {
		strcat (filename,".");
		strcat (filename,refcatname);
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
	sprintf (headline, "CATALOG	USNO SA 1.0");
    else if (refcat == UAC)
	sprintf (headline, "CATALOG	USNO A 1.0");
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
    else if (refcat == UAC)
	strcpy (headline,"usnoa_id  	");
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
    if (refcat == USAC || refcat == UAC)
	strcat (headline,"magb	magr	plate");
    else if (refcat == TYCHO)
	strcat (headline,"magb	magv	type");
    else
	strcat (headline,"mag	type");
    if (snum == NULL)
	strcat (headline,"	arcsec");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (refcat == UAC || refcat == USAC || refcat == TYCHO)
	sprintf(headline,"----------	--------	---------	----	-----	-----	-----");
    else
        sprintf (headline,"----------	------------	------------	------	----");
    if (snum == NULL)
	strcat (headline, " ------");
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
	    else if (refcat == UAC)
		printf ("USNO A number  ");
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
	    if (refcat == UAC || refcat == USAC)
		printf ("MagB  MagR Plate");
	    else if (refcat == GSC)
		printf (" Mag  Type");
	    else if (refcat == SAO || refcat == PPM || refcat == IRAS)
		printf (" Mag  Type");
	    else if (refcat == TYCHO)
		printf (" MagB  MagV  Type");
	    else if (tabcat)
		printf (" Mag     Peak");
	    else
		printf ("    Mag");
	    if (snum == NULL)
		printf ("  Arcsec\n");
	    else
		printf ("\n");
	    }
	}

    if (keyword != NULL)
	ntab = tabcatopen (refcat);

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
	    if (refcat == UAC || refcat == USAC)
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
	    else {
	        sprintf (headline, "%9.4f	%s	%s	%.2f	%d",
		 gnum[i], rastr, decstr, gm[i], gc[i]);
		if (keyword != NULL) {
		    tabline = tabstar ((int)gnum[i], tabline);
		    (void) tabgetk (tabline, keyword, string, TABMAX);
		    }
		}
	    if (snum == NULL) {
	        sprintf (temp, "	%.2f", gdist);
	        strcat (headline, temp);
		}
	    if (keyword != NULL) {
		strcat (headline, "	");
		strcat (headline, string);
		}
	    if (wfile)
		fprintf (fd, "%s\n", headline);
	    else if (tabout)
		printf ("%s\n", headline);
	    else if (refcat == UAC || refcat == USAC)
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
	    else if (tabcat) {
		printf ("%9.4f %s %s %6.2f %7d",
			gnum[i], rastr, decstr, gm[i],gc[i]);
		if (keyword != NULL) {
		    tabline = tabstar ((int)gnum[i], tabline);
		    (void) tabgetk (tabline, keyword, string, TABMAX);
		    }
		}
	    else if (tabcat)
		printf ("%9.4f %s %s %6.2f %7d",
			gnum[i], rastr, decstr, gm[i],gc[i]);
	    else
		printf ("%9.4f %s %s %6.2f",
			gnum[i], rastr, decstr, gm[i]);
	    if (snum == NULL)
		printf ("  %7.2f", gdist);
	    if (keyword != NULL)
		printf (" %s\n", string);
	    else if (refcat < 0 && !tabcat && gobj != NULL && gobj[i] != NULL)
		printf (" %s\n", gobj[i]);
	    else
		printf ("\n");
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
SearchHead (sys1, sys2, eq1, eq2, ep1, ep2, cra, cdec, dra, ddec, drad)

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
	if (refcat == GSC)
	    printf ("HST GSC   %s %s ", rastr, decstr);
	else if (refcat == USAC)
	    printf ("USNO SA 1.0   %s %s ", rastr, decstr);
	else if (refcat == UAC)
	    printf ("USNO A 1.0    %s %s ", rastr, decstr);
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
    if (sys2 == WCS_B1950) {
	if (eq2 == 1950.0)
	    printf ("B1950 ");
	else
	    printf ("B%7.2f ", eq2);
	}
    else if (sys2 == WCS_GALACTIC)
	printf ("galactic ");
    else if (sys2 == WCS_ECLIPTIC)
	printf ("ecliptic ");
    else {
	if (eq2 == 2000.0)
	    printf ("J2000 ");
	else
	    printf ("J%7.2f ", eq2);
	}
    if (drad != 0.0)
	printf ("r= %.2f", drad*3600.0);
    else
	printf ("+- %.2f", ddec*3600.0);
    if (epoch0 != 0.0)
	printf (" at epoch %7.2f\n", epoch0);
    else
	printf ("\n");
    return;
}

/* Oct 21 1998	New program based on imcat
 * Nov 30 1998	Add version and help commands for consistency

 * Jan 19 1999	Update USNO A and SA catalog reading
 * Sep 10 1999	Do all searches through catread() and catrnum()
 */
