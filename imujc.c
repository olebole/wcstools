/* File imujc.c
 * July 17, 1996
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
static char *RevMsg = "IMUJC 1.0, 8 August 1996, Doug Mink, SAO";

extern int gscread();
static void ListUJC ();
extern void RASortStars ();
extern void MagSortStars ();

static double magoff = 0.0;
static double maglim = 0.0;
static int rasort = 0;
static int printhead = 0;
static int tabout = 0;
static int plate = 0;
static char coorsys[4];
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
	case 'p':	/* Plate number to use for objects */
	    if (ac < 2)
		usage (progname);
	    plate = atoi (*++av);
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

	ListUJC (fn, nstar);
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
    fprintf (stderr,"Find UJC stars in FITS or IRAF image files\n");
    fprintf(stderr,"IMUJC: usage: [-v] [-m mag_off] [-n num] [-p plate] file.fts ...\n",
	    progname);
    fprintf(stderr,"  -b: output B1950 (FK4) coordinates \n");
    fprintf(stderr,"  -c: Use following RA and Dec as center \n");
    fprintf(stderr,"  -h: print heading, else do not \n");
    fprintf(stderr,"  -j: output J2000 (FK5) coordinates \n");
    fprintf(stderr,"  -m: magnitude limit\n");
    fprintf(stderr,"  -n: number of brightest stars to print \n");
    fprintf(stderr,"  -p: plate number for catalog sources (0=all)\n");
    fprintf(stderr,"  -s: sort by RA instead of flux \n");
    fprintf(stderr,"  -t: tab table to standard output as well as file\n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


extern int findStars ();
struct WorldCoor *wcsinit();	
extern int pix2wcst();

static void
ListUJC (filename, nstars)

char	*filename;	/* FITS or IRAF file filename */
int	nstars;		/* Number of brightest stars to list */

{
    char *header;	/* FITS header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    int iraffile;	/* 1 if IRAF image */
    int *irafheader;	/* IRAF image header */
    double *unum=0;	/* UJC star numbers */
    double *ura=0;	/* UJC star right ascensions, rads */
    double *udec=0;	/* UJC star declinations rads */
    double *um=0;	/* UJC magnitudes */
    double *ux, *uy;	/* UJC positions on image */
    int *up;		/* UJC plate numbers */
    int nu;		/* Number of UJC stars */
    int nbu;		/* Number of brightest GCS stars actually used */
    int i, numax, nbytes;
    FILE *fd;
    struct WorldCoor *wcs;	/* World coordinate system structure */
    char rastr[16], decstr[16];	/* coordnate strings */
    double cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, mag1, mag2;
    int offscale, nlog, lfn;
    char headline[160];

    /* Open IRAF image if .imh extension is present */
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

    /* Open FITS file if .imh extension is not present */
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

/* Set magnitude limits for the UJ Catalog search */
    if (maglim == 0.0) {
	mag1 = 0.0;
	mag2 = 0.0;
	}
    else {
	mag1 = -2.0;
	mag2 = maglim;
	}

    numax = MAXREF;
    nbytes = MAXREF * sizeof (double);
    unum = (double *) malloc (nbytes);
    ura = (double *) malloc (nbytes);
    udec = (double *) malloc (nbytes);
    um = (double *) malloc (nbytes);
    up = (int *) malloc (nbytes);
    nlog = 0;

    /* Find the nearby reference stars, in ra/dec */
    nu = ujcread (ra1,ra2,dec1,dec2,mag1,mag2,plate,numax,unum,ura,udec,
		  um,up,nlog);

    /* Project the reference stars into pixels on a plane at ra0/dec0 */
    ux = (double *) malloc (nu * sizeof (double));
    uy = (double *) malloc (nu * sizeof (double));
    if (!ux || !uy) {
	fprintf (stderr, "Could not malloc temp space of %d bytes\n",
					    nu*sizeof(double)*2);
	return;
	}
    for (i = 0; i < nu; i++ ) {
	offscale = 0;
	if (strcmp (wcs->sysout,"FK4") == 0 || strcmp (wcs->sysout,"fk4") == 0)
	    fk524e (&ura[i], &udec[i], wcs->epoch);
	wcs2pix (wcs, ura[i], udec[i], &ux[i], &uy[i], &offscale);
	}

    /* Sort reference stars by brightness (magnitude) */
    MagSortStars (unum, ura, udec, ux, uy, um, up, nu);

    if (nstars > 0 && nu > nstars) {
	nbu = nstars;
	if (verbose && printhead)
	    printf ("using %d / %d UJ Catalog Stars brighter than %.1f",
		    nbu, nu, um[nbu-1]);
	}
    else {
	nbu = nu;
	if (verbose && printhead) {
	    if (maglim > 0.0)
		printf ("%d UJ Catalog Stars brighter than %.1f",
			nu, maglim);
	    else
		printf ("%d UJ Catalog Stars", nu);
	    }
	if (iraffile)
	    printf (" in IRAF image %s\n",filename);
	else
	    printf (" in FITS image %s\n", filename);
	}
    sprintf (headline, "IMAGE	%s", filename);

    /* Sort star-like objects in image by right ascension */
    if (rasort)
	RASortStars (unum, ura, udec, ux, uy, um, up, nbu);

    /* Open plate catalog file */
    strcat (filename,".ujc");
    fd = fopen (filename, "w");
    if (fd == NULL) {
	fprintf (stderr, "IMUJC:  cannot write file %s %d\n", filename, fd);
        return;
        }

    /* Write header */
    fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    sprintf (headline, "CATALOG	UJ1.0");
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
    if (rasort)

    fprintf (fd, "PROGRAM	%s\n", RevMsg);
    if (tabout)
	printf ("PROGRAM	%s\n", RevMsg);

    sprintf (headline,"UJC_NUMBER	RA      	DEC      	MAG   	X    	Y	Plate");
    fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    sprintf (headline,"----------	------------	------------	------	-----	-----	----");
    fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    for (i = 0; i < nbu; i++) {
	if (ux[i] > 0.0 && uy[i] > 0.0) {
	    ra2str (rastr, ura[i], 3);
	    dec2str (decstr, udec[i], 2);
	    sprintf (headline, "%12.7f	%s	%s	%.2f	%.1f	%.1f	%d",
		 unum[i], rastr, decstr, um[i], ux[i], uy[i], up[i]);
	    fprintf (fd, "%s\n", headline);
	    if (tabout)
		printf ("%s\n", headline);
	    else if (verbose)
		printf ("%12.7f %s %s %6.2f %6.1f %6.1f %5.2f %d\n",
		    unum[i], rastr, decstr, um[i],ux[i],uy[i],um[i],up[i]);
	    }
	}

    fclose (fd);
    if (ux) free ((char *)ux);
    if (uy) free ((char *)uy);
    if (um) free ((char *)um);
    if (ura) free ((char *)ura);
    if (udec) free ((char *)udec);
    if (unum) free ((char *)unum);
    if (up) free ((char *)up);
    free ((char *)wcs);

    free (header);
    return;
}

/* May 21 1996	New program
 * May 23 1996	Add optional selection by plate number
 * May 28 1996	Clean up coordinate conversions
 * Jul 17 1996	Update image header reading
 */
