/* File imgsc.c
 * July 16, 1996
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
#include "libwcs/fitshead.h"
#include "libwcs/wcs.h"

#define MAXREF 100

static void usage();
static int verbose = 0;		/* verbose/debugging flag */
static char *RevMsg = "IMGSC 1.0, 8 August 1996, Doug Mink SAO";

extern int gscread();
static void ListGSC ();
extern void RASortStars ();
extern void MagSortStars ();

static double magoff = 0.0;
static double maglim = 0.0;
static int rasort = 0;
static int classd = -1;
static char coorsys[4];
static int printhead = 0;
static int tabout = 0;
static double cra0 = 0.0;
static double cdec0 = 0.0;

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    int nstar = 0;
    double starsig, bmin;
    int maxrad;

    *coorsys = 0;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {
	case 'v':	/* more verbosity */
	    verbose++;
	    break;
	case 'b':	/* ouput FK4 (B1950) coordinates */
	    strcpy (coorsys, "FK4");
	    break;
	case 'c':	/* Set center RA and Dec */
	    if (ac < 3)
		usage (progname);
	    cra0 = str2ra (*++av);
	    ac--;
	    cdec0 = str2dec (*++av);
	    ac--;
	    break;
	case 'h':	/* ouput descriptive header */
	    printhead++;
	    break;
	case 'j':	/* ouput FK5 (J2000) coordinates */
	    strcpy (coorsys, "FK5");
	    break;
	case 'o':	/* Object type */
	    if (ac < 2)
		usage (progname);
	    classd = (int) (atof (*++av) + 0.1);
	    ac--;
	    break;
	case 'm':	/* Magnitude limit */
	    if (ac < 2)
		usage (progname);
	    maglim = atof (*++av);
	    ac--;
	    break;
	case 'n':	/* Number of brightest stars to read */
	    if (ac < 2)
		usage (progname);
	    nstar = atoi (*++av);
	    ac--;
	    break;
	case 's':	/* sort by RA */
	    rasort = 1;
	    break;
	case 't':	/* tab table to stdout */
	    tabout = 1;
	    break;
	default:
	    usage (progname);
	    break;
	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage (progname);

    while (ac-- > 0) {
	char *fn = *av++;

	ListGSC (fn, nstar);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"%s\n",RevMsg);
    fprintf (stderr,"Find HST GSC stars in FITS or IRAF image files\n");
    fprintf(stderr,"%s: usage: [-v] [-m mag_off] [-n num] [-c ra dec] file.fts ...\n",
	    progname);
    fprintf(stderr,"  -b: output B1950 (FK4) coordinates \n");
    fprintf(stderr,"  -c: Use following RA and Dec as center \n");
    fprintf(stderr,"  -h: print heading, else do not \n");
    fprintf(stderr,"  -j: output J2000 (FK5) coordinates \n");
    fprintf(stderr,"  -m: magnitude limit\n");
    fprintf(stderr,"  -n: number of brightest stars to print \n");
    fprintf(stderr,"  -o: object type (0=stars 3=galaxies -1=all)\n");
    fprintf(stderr,"  -s: sort by RA instead of flux \n");
    fprintf(stderr,"  -t: tab table to standard output as well as file\n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


struct WorldCoor *wcsinit();	
extern int pix2wcst();

static void
ListGSC (filename, nstars)

char	*filename;	/* FITS or IRAF file filename */
int	nstars;		/* Number of brightest stars to list */

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
    double cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, mag1, mag2, ra, dec;
    int offscale, nlog;
    char headline[160];

    /* Open IRAF header if .imh extension is present */
    if (strsrch (filename,".imh") != NULL) {
	iraffile = 1;
	irafheader = irafrhead (filename, &lhead);
	if (irafheader) {
	    header = iraf2fits (filename, irafheader, lhead, &nbhead);
	    free (irafheader);
	    if (header == NULL) {
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
	header = fitsrhead (filename, &lhead, &nbhead);
	if (header == NULL) {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    return;
	    }
	}
    if (verbose && printhead)
	printf ("%s\n",RevMsg);

    /* Read world coordinate system information from the image header */
    wcs = wcsinit (header);

    /* Set the RA and Dec limits in degrees for reference star search */
    wcssize (wcs, &cra, &cdec, &dra, &ddec);
    if (cra0 > 0.0)
	cra = cra0;
    if (cdec0 != 0.0)
	cdec = cdec0;
    if (cra0 > 0.0 || cdec0 != 0.0) {
	if (coorsys[1])
	    wcsshift (wcs,cra,cdec,coorsys);
	else
	    wcsshift (wcs,cra,cdec,wcs->radecsys);
	}
    if (strcmp (wcs->radecsys,"FK4") == 0)
	fk425 (&cra, &cdec);
    ra1 = cra - dra;
    ra2 = cra + dra;
    dec1 = cdec - ddec;
    dec2 = cdec + ddec;
    if (verbose && printhead) {
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
    nlog = 0;

    /* Find the nearby reference stars, in ra/dec */
    ng = gscread (ra1,ra2,dec1,dec2,mag1,mag2,classd,ngmax,gnum,gra,gdec,gm,gc,
		  nlog);

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
	if (verbose && printhead)
	    printf ("using %d / %d HST Guide Stars brighter than %.1f\n",
		     nbg,ng,gm[nbg-1]);
	}
    else {
	nbg = ng;
	if (verbose && printhead) {
	    if (maglim > 0.0)
		printf ("%d HST Guide Stars brighter than %.1f",
			ng, maglim);
	    else
		printf ("%d HST Guide Stars", ng);
	    }
	}
    if (verbose && printhead) {
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
    strcat (filename,".gsc");
    fd = fopen (filename, "w");
    if (fd == NULL) {
	fprintf (stderr, "IMGSC:  cannot write file %s %d\n", filename, fd);
        return;
        }

    /* Write header */
    fprintf (fd,"%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    sprintf (headline, "CATALOG	HSTGSC1.1");
    fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (rasort) {
	fprintf (fd, "RASORT	T\n");
	if (tabout)
	    printf ("RASORT	T\n");
	}

    if (strcmp (wcs->sysout,"FK4") == 0 || strcmp (wcs->sysout,"fk4") == 0)
	sprintf (headline, "EQUINOX	1950.0");
    else
	sprintf (headline, "EQUINOX	2000.0");
    fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    fprintf (fd, "PROGRAM	%s\n", RevMsg);
    if (tabout)
	printf ("PROGRAM	%s\n", RevMsg);

    sprintf (headline,"GSC_NUMBER	RA      	DEC      	MAG   	X    	Y	Type");
    fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    sprintf (headline,"----------	------------	------------	------	-----	-----	----");
    fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    for (i = 0; i < nbg; i++) {
	if (gx[i] > 0.0 && gy[i] > 0.0) {
	    ra2str (rastr, gra[i], 3);
	    dec2str (decstr, gdec[i], 2);
	    sprintf (headline, "%9.4f	%s	%s	%.2f	%.1f	%.1f	%d",
		 gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i]);
	    fprintf (fd, "%s\n", headline);
	    if (tabout)
		printf ("%s\n", headline);
	    else if (verbose)
		printf ("%9.4f %s %s %6.2f %6.1f %6.1f  %d\n",
			gnum[i], rastr, decstr, gm[i],gx[i],gy[i],gc[i]);
	    }
	}

    fclose (fd);
    if (gm) free ((char *)gm);
    if (gra) free ((char *)gra);
    if (gdec) free ((char *)gdec);
    if (gnum) free ((char *)gnum);
    if (gc) free ((char *)gc);
    free ((char *)wcs);

    free (header);
    return;
}

/* May 21 1996	New program
 * Jul 11 1996	Update file reading
 * Jul 16 1996	Remove unused image pointer; do not free unallocated header
 */
