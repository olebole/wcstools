/* File imcat.c
 * October 1, 1996
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

#include "libwcs/fitshead.h"
#include "libwcs/wcs.h"
#include "libwcs/lwcs.h"

#define GSC	1	/* refcat value for HST Guide Star Catalog */
#define UJC	2	/* refcat value for USNO UJ Star Catalog */

extern int gscread();
extern int ujcread();
extern int tabread();
static void usage();
static void ListCat();
extern void RASortStars();
extern void MagSortStars();
extern void fk524e();
extern void fk425e();


static int verbose = 0;			/* verbose/debugging flag */
static int wfile = 0;			/* True to print output file */
static int refcat = GSC;		/* reference catalog switch */
static char refcatname[32]="GSC";	/* reference catalog name */
static int fk4 = 0;			/* Command line center is FK4 */
static int classd = -1;			/* Guide Star Catalog object classes */
static int uplate = 0;			/* UJ Catalog plate number to use */
static double maglim = MAGLIM;		/* reference catalog magnitude limit */
static double cra0 = 0.0;		/* Initial center RA in degrees */
static double cdec0 = 0.0;		/* Initial center Dec in degrees */
static char coorsys[4];			/* Output coordinate system */
static int nstars;			/* Number of brightest stars to list */
static int printhead = 0;		/* 1 to print table heading */
static int tabout = 0;			/* 1 for tab table to standard output */
static int rasort = 0;			/* 1 to sort stars by brighness */
static double secpix0 = 0;		/* Image plate scale in arcsec/pixel */
static int debug = 0;			/* True for extra information */

main (ac, av)
int ac;
char **av;
{
    char *str;

    /* Decode arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
    	switch (c) {

    	case 'b':	/* initial coordinates on command line in B1950 */
	    setfk4 ();
    	    if (ac < 3)
    		usage();
	    cra0 = str2ra (*++av);
	    ac--;
	    cdec0 = str2dec (*++av);
	    ac--;
	    strcpy (coorsys, "FK4");
    	    break;

	case 'c':       /* Set reference catalog */
	    if (ac < 2)
		usage();
	    strcpy (refcatname, *++av);
	    if (refcatname[0]=='G' || refcatname[0]=='g' ||
		strcmp(refcatname,"gsc")==0 || strcmp (refcatname,"GSC")== 0)
		refcat = GSC;
	    else if (refcatname[0]=='U' || refcatname[0]=='u' ||
		strcmp(refcatname,"ujc")==0 || strcmp(refcatname,"UJC")==0)
		refcat = UJC;
	    else
		refcat = 0;
	    ac--;
	    break;

	case 'd':
	    debug++;
	    break;

	case 'h':	/* ouput descriptive header */
	    printhead++;
	    break;

    	case 'j':	/* center coordinates on command line in J2000 */
    	    if (ac < 3)
    		usage();
	    cra0 = str2ra (*++av);
	    ac--;
	    cdec0 = str2dec (*++av);
	    ac--;
	    strcpy (coorsys, "FK5");
    	    break;

	case 'g':	/* Guide Star object class */
    	    if (ac < 2)
    		usage();
    	    classd = (int) atof (*++av);
    	    ac--;
    	    break;

    	case 'm':	/* Limiting reference star magnitude */
    	    if (ac < 2)
    		usage();
    	    maglim =  atof (*++av);
    	    ac--;
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
    	    secpix0 = atof (*++av);
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
    	default:
    	    usage();
    	    break;
    	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage();
    if (!verbose && !wfile)
	verbose = 1;

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
    fprintf (stderr,"List catalog stars in FITS and IRAF image files\n");
    fprintf (stderr,"Usage: [-vhst] [-m mag] [-g class] [-c catalog]\n");
    fprintf (stderr,"       [-p scale] [-b ra dec] [-j ra dec] FITS or IRAF file(s)\n");
    fprintf (stderr,"  -b: initial center in B1950 (FK4) RA and Dec\n");
    fprintf (stderr,"  -c: reference catalog (gsc, ujc, or tab table file\n");
    fprintf (stderr,"  -g: Guide Star Catalog class (-1=all,0,3 (default -1)\n");
    fprintf (stderr,"  -h: print heading, else do not \n");
    fprintf (stderr,"  -j: initial center in J2000 (FK5) RA and Dec\n");
    fprintf (stderr,"  -m: initial GSC or UJC limiting magnitude (default 17)\n");
    fprintf (stderr,"  -n: number of brightest stars to print \n");
    fprintf (stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf (stderr,"  -s: sort by RA instead of flux \n");
    fprintf (stderr,"  -t: tab table to standard output as well as file\n");
    fprintf (stderr,"  -u: UJ catalog single plate number to accept\n");
    fprintf (stderr,"  -v: verbose\n");
    exit (1);
}


struct WorldCoor *wcsinit();	
extern int pix2wcst();

static void
ListCat (filename)

char	*filename;	/* FITS or IRAF file filename */

{
    char *header;	/* FITS header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    int iraffile;	/* 1 if IRAF image */
    int *irafheader;	/* IRAF image header */
    double *gnum=0;	/* GSC star numbers */
    double *gra=0;	/* GSC star right ascensions, rads */
    double *gdec=0;	/* GSC star declinations rads */
    double *gm=0;	/* GCS magnitudes */
    double *gx, *gy;	/* GSC positions on image */
    int *gc;		/* GSC object classes */
    int ng;		/* Number of GSC stars */
    int nbg;		/* Number of brightest GCS stars actually used */
    int i, ngmax, nbytes;
    FILE *fd;
    struct WorldCoor *wcs;	/* World coordinate system structure */
    char rastr[16], decstr[16];	/* coordinate strings */
    double cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, mag1, mag2;
    int offscale, nlog;
    char headline[160];

    /* Open IRAF header if .imh extension is present */
    if (strsrch (filename,".imh")) {
	iraffile = 1;
	if ((irafheader = irafrhead (filename, &lhead))) {
	    header = iraf2fits (filename, irafheader, lhead, &nbhead);
	    free (irafheader);
	    if (!header) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",filename);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", filename);
	    return;
	    }
	}

    /* Read FITS image header if .imh extension is not present */
    else {
	iraffile = 0;
	if (!(header = fitsrhead (filename, &lhead, &nbhead))) {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    return;
	    }
	}
    if (verbose && printhead)

    /* Set plate scale from command line */
    if (secpix0 > 0)
	hputr8 (header, "SECPIX", secpix0);

    /* Read world coordinate system information from the image header */
    wcs = wcsinit (header);
    free (header);

    /* Set the RA and Dec limits in degrees for reference star search */
    wcssize (wcs, &cra, &cdec, &dra, &ddec);

    /* Reset the center if new center is on command line */
    if (cra0 > 0.0 || cdec0 != 0.0) {
	if (cra0 > 0.0)
	    cra = cra0;
	if (cdec0 != 0.0)
	    cdec = cdec0;
	if (coorsys[1])
	    wcsshift (wcs,cra,cdec,coorsys);
	else
	    wcsshift (wcs,cra,cdec,wcs->radecsys);
	}
    if (strcmp (wcs->radecsys,"FK4") == 0)
	fk425e (&cra, &cdec, wcs->epoch);
    ra1 = cra - dra;
    ra2 = cra + dra;
    dec1 = cdec - ddec;
    dec2 = cdec + ddec;
    if (debug) {
	char rastr1[16],rastr2[16],decstr1[16],decstr2[16];
	ra2str (rastr1,ra1,3);
	ra2str (rastr2,ra2,3);
	printf ("%s: RA:  %s - %s (J2000)\n",filename,rastr1,rastr2);
	dec2str (decstr1,dec1,2);
	dec2str (decstr2,dec2,2);
	printf ("%s: Dec: %s - %s (J2000)\n",filename, decstr1,decstr2);
	}

/* Set the output coordinate system */
    wcs->printsys = 0;
    if (coorsys[1])
	wcsoutinit (wcs, coorsys);

/* Set the magnitude limits for the GSC search */
    if (maglim == 0.0) {
	mag1 = 0.0;
	mag2 = 0.0;
	}
    else {
	mag1 = -2.0;
	mag2 = maglim;
	}

    ngmax = MAXREF;
    nbytes = MAXREF * sizeof (double);
    gnum = (double *) malloc (nbytes);
    gra = (double *) malloc (nbytes);
    gdec = (double *) malloc (nbytes);
    gm = (double *) malloc (nbytes);
    gc = (int *) malloc (nbytes);
    if (debug)
	nlog = 100;
    else
	nlog = 0;

    /* Find the nearby reference stars, in ra/dec */
    if (refcat == GSC)
	ng = gscread (ra1,ra2,dec1,dec2,mag1,mag2,classd,ngmax,
		      gnum,gra,gdec,gm,gc,nlog);
    else if (refcat == UJC)
	ng = ujcread (ra1,ra2,dec1,dec2,mag1,mag2,uplate,ngmax,
		      gnum,gra,gdec,gm,gc,debug);
    else if (refcatname[0] > 0)
	ng = tabread (refcatname,ra1,ra2,dec1,dec2,mag1,mag2,ngmax,
		      gnum,gra,gdec,gm,gc,debug);
    else {
        fprintf (stderr,"No reference star catalog specified\n");
	free ((char *)wcs);
        return;
        }

    /* Project the reference stars into pixels on a plane at ra0/dec0 */
    gx = (double *) malloc (ng * sizeof (double));
    gy = (double *) malloc (ng * sizeof (double));
    if (!gx || !gy) {
        fprintf (stderr, "Could not malloc temp space of %d bytes\n",
			 ng*sizeof(double)*2);
	return;
	}
    for (i = 0; i < ng; i++ ) {
	offscale = 0;
	if (strcmp (wcs->sysout,"FK4") == 0 || strcmp (wcs->sysout,"fk4") == 0)
	    fk524e (&gra[i],&gdec[i], wcs->epoch);
	wcs2pix (wcs, gra[i], gdec[i], &gx[i], &gy[i], &offscale);
	}

    /* Sort reference stars by brightness (magnitude) */
    MagSortStars (gnum, gra, gdec, gx, gy, gm, gc, ng);

    /* List the brightest MAXSTARS reference stars */
    if (nstars > 0 && ng > nstars) {
	nbg = nstars;
	if (debug)
	    printf ("using %d / %d HST Guide Stars brighter than %.1f\n",
		     nbg,ng,gm[nbg-1]);
	}
    else {
	nbg = ng;
	if (debug) {
	    if (maglim > 0.0) {
		if (refcat == GSC)
		    printf ("%d HST Guide Stars brighter than %.1f",ng,maglim);
		else if (refcat == UJC)
		    printf ("%d USNO J Catalog Stars brighter than %.1f",ng,maglim);
		else
		    printf ("%d %s stars brighter than %.1f",ng,refcatname,maglim);
		}
	    else {
		if (refcat == GSC)
		    printf ("%d HST Guide Stars", ng);
		else if (refcat == UJC)
		    printf ("%d USNO J Catalog Stars", ng);
		else
		    printf ("%d %s Catalog Stars", ng,refcatname);
		}
	    }
	}
    if (debug) {
	if (iraffile)
	    printf (" in IRAF image %s\n",filename);
	else
	    printf (" in FITS image %s\n", filename);
	}

    /* Sort star-like objects in image by right ascension */
    if (rasort)
	RASortStars (gnum, gra, gdec, gx, gy, gm, gc, nbg);
    sprintf (headline, "IMAGE	%s", filename);

    /* Open plate catalog file */
    if (refcat == GSC)
	strcat (filename,".gsc");
    else if (refcat == UJC)
	strcat (filename,".ujc");
    else
	strcat (filename,".cat");
    if (wfile) {
	fd = fopen (filename, "w");
	if (fd == NULL) {
	    fprintf (stderr, "IMGSC:  cannot write file %s %d\n", filename, fd);
            return;
	    }
        }

    /* Write heading */
    if (wfile)
	fprintf (fd,"%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    sprintf (headline, "CATALOG	HSTGSC1.1");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (rasort) {
	if (wfile)
	    fprintf (fd, "RASORT	T\n");
	if (tabout)
	    printf ("RASORT	T\n");
	}

    if (strcmp (wcs->sysout,"FK4") == 0 || strcmp (wcs->sysout,"fk4") == 0)
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
	sprintf (headline,"GSC_NUMBER	RA      	DEC      	MAG   	X    	Y	Type");
    else if (refcat == UJC)
	sprintf (headline,"UJC_NUMBER	RA      	DEC      	MAG   	X    	Y	Plate");
    else
	sprintf (headline,"NUMBER	RA      	DEC      	MAG   	X    	Y	Peak");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    sprintf (headline,"----------	------------	------------	------	-----	-----	----");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    for (i = 0; i < nbg; i++) {
	if (gx[i] > 0.0 && gy[i] > 0.0) {
	    ra2str (rastr, gra[i], 3);
	    dec2str (decstr, gdec[i], 2);
	    if (refcat == GSC)
		sprintf (headline, "%9.4f	%s	%s	%.2f	%.1f	%.1f	%d",
		 gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i]);
	    else if (refcat == UJC)
		sprintf (headline, "%12.7f	%s	%s	%.2f	%.1f	%.1f	%d",
		 gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i]);
	    else
		sprintf (headline, "%9.4f	%s	%s	%.2f	%.1f	%.1f	%d",
		 gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i]);
	    if (wfile)
		fprintf (fd, "%s\n", headline);
	    if (tabout)
		printf ("%s\n", headline);
	    else if (verbose) {
		if (refcat == GSC)
		    printf ("%9.4f %s %s %6.2f %6.1f %6.1f  %d\n",
			gnum[i], rastr, decstr, gm[i],gx[i],gy[i],gc[i]);
		else if (refcat == UJC)
		    printf ("%12.7f %s %s %6.2f %6.1f %6.1f  %d\n",
			gnum[i], rastr, decstr, gm[i],gx[i],gy[i],gc[i]);
		else
		    printf ("%9.4f %s %s %6.2f %6.1f %6.1f  %d\n",
			gnum[i], rastr, decstr, gm[i],gx[i],gy[i],gc[i]);
		}
	    }
	}

    if (wfile)
	fclose (fd);
    if (gm) free ((char *)gm);
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
 */
