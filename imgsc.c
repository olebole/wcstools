/* File imgsc.c
 * December 12, 1996
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
#include "libwcs/lwcs.h"

#define MAXREF 100

static void usage();

extern int gscread();
static void ListGSC ();
extern void RASortStars ();
extern void MagSortStars ();
extern void fk524e();
extern void setfk4();
extern void setcenter();
extern void setsecpix();
extern struct WorldCoor *GetFITSWCS();

static int verbose = 0;		/* Verbose/debugging flag */
static int wfile = 0;		/* True to print output file */
static int classd = -1;		/* Guide Star Catalog object classes */
static double maglim1 = MAGLIM1; /* Guide Star Catalog bright magnitude limit */
static double maglim2 = MAGLIM;	/* Guide Star Catalog faint magnitude limit */
static char coorsys[4];		/* Output coordinate system */
static int nstars = 0;		/* Number of brightest stars to list */
static int printhead = 0;	/* 1 to print table heading */
static int tabout = 0;		/* 1 for tab table to standard output */
static int rasort = 0;		/* 1 to sort stars by brighness */
static int debug = 0;		/* True for extra information */

main (ac, av)
int ac;
char **av;
{
    char *str, *str1;
    char rastr[16];
    char decstr[16];

    *coorsys = 0;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {
	case 'v':	/* more verbosity */
	    verbose++;
	    break;
    	case 'b':	/* initial coordinates on command line in B1950 */
	    strcpy (coorsys, "FK4");
	    str1 = *(av+1);
	    if (*(str+1) || (str1[0] > 47 && str[0] < 58))
		setfk4();
	    else if (ac < 3)
		usage ();
	    else {
		setfk4 ();
		strcpy (rastr, *++av);
		ac--;
		strcpy (decstr, *++av);
		ac--;
		setcenter (rastr, decstr);
		}
    	    break;

	case 'd':
	    debug++;
	    break;

	case 'g':	/* Object type */
	    if (ac < 2)
		usage ();
	    classd = (int) (atof (*++av) + 0.1);
	    ac--;
	    break;

	case 'h':	/* ouput descriptive header */
	    printhead++;
	    break;

    	case 'j':	/* center coordinates on command line in J2000 */
	    str1 = *(av+1);
	    if (*(str+1) || (str1[0] > 47 && str[0] < 58))
		strcpy (coorsys, "FK5");
	    else if (ac < 3)
		usage ();
	    else {
		strcpy (coorsys, "FK5");
		strcpy (rastr, *++av);
		ac--;
		strcpy (decstr, *++av);
		ac--;
		setcenter (rastr, decstr);
		}
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

    	case 'w':	/* write output file */
    	    wfile++;
    	    break;

	default:
	    usage ();
	    break;
	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage ();

    while (ac-- > 0) {
	char *fn = *av++;

	ListGSC (fn);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"Find HST GSC stars in FITS or IRAF image files\n");
    fprintf (stderr,"Usage: [-vhstw] [-m [mag1] mag2] [-n num] [-p scale]\n");
    fprintf (stderr,"       [-b ra dec] [-j ra dec] FITS or IRAF file(s)\n");
    fprintf(stderr,"  -b: output B1950 (FK4) coordinates (optional center)\n");
    fprintf(stderr,"  -g: object type (0=stars 3=galaxies -1=all)\n");
    fprintf(stderr,"  -h: print heading, else do not \n");
    fprintf(stderr,"  -j: output J2000 (FK5) coordinates (optional center)\n");
    fprintf(stderr,"  -m: magnitude limit\n");
    fprintf(stderr,"  -n: number of brightest stars to print \n");
    fprintf(stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -s: sort by RA instead of flux \n");
    fprintf(stderr,"  -t: tab table to standard output as well as file\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -w: Write tab table output file imagename.gsc\n");
    exit (1);
}


struct WorldCoor *wcsinit();	
extern int pix2wcst();

static void
ListGSC (filename)

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
    double *gm=0;	/* GSC magnitudes */
    double *gmb=0;	/* Unused magnitudes */
    double *gx=0;	/* GSC X positions on image */
    double *gy=0;	/* GSC Y positions on image */
    int *gc=0;		/* GSC object classes */
    int ng;		/* Number of GSC stars */
    int nbg;		/* Number of brightest GCS stars actually used */
    int imh, imw;	/* Image height and width in pixels */
    int i, ngmax, nbytes;
    FILE *fd;
    struct WorldCoor *wcs;	/* World coordinate system structure */
    char rastr[16], decstr[16];	/* coordinate strings */
    double cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, mag1, mag2, secpix;
    double mag;
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
    if (verbose || printhead)

    /* Read world coordinate system information from the image header */
    wcs = GetFITSWCS (header, verbose, &cra, &cdec, &dra, &ddec, &secpix,
                &imw, &imh, 2000.0);
    free (header);
    if (nowcs (wcs))
        return;

    /* Set limits for search */
    ra1 = cra - dra;
    ra2 = cra + dra;
    dec1 = cdec - ddec;
    dec2 = cdec + ddec;
    if (verbose || printhead) {
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
    if (verbose)
	nlog = 100;
    else
	nlog = 0;

    /* Find the nearby reference stars, in ra/dec */
    ng = gscread (cra,cdec,dra,ddec,0.0,mag1,mag2,classd,ngmax,
		  gnum,gra,gdec,gm,gc,nlog);

    /* Project the reference stars into pixels on a plane at ra0/dec0 */
    nbytes = ng * sizeof (double);
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
    for (i = 0; i < ng; i++ ) {
	offscale = 0;
	if (strcmp (wcs->sysout,"FK4") == 0 || strcmp (wcs->sysout,"fk4") == 0)
	    fk524e (&gra[i],&gdec[i], wcs->epoch);
	wcs2pix (wcs, gra[i], gdec[i], &gx[i], &gy[i], &offscale);
	}

    /* Sort reference stars by brightness (magnitude) */
    MagSortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, ng);

    /* List the brightest MAXSTARS reference stars */
    if (nstars > 0 && ng > nstars) {
	nbg = nstars;
	if (verbose || printhead)
	    if (maglim1 > 0.0)
		printf ("%d / %d HST Guide Stars (between %.1f and %.1f)",
			nbg, ng, gm[0], gm[nbg-1]);
	    else
		printf ("%d / %d HST Guide Stars (brighter than %.1f)",
			nbg, ng, gm[nbg-1]);
	}
    else {
	nbg = ng;
	if (verbose || printhead) {
	    if (maglim1 > 0.0)
		printf ("%d HST Guide Stars between %.1f and %.1f",
			 ng, maglim1, maglim2);
	    else if (maglim2 > 0.0)
		printf ("%d HST Guide Stars brighter than %.1f", ng, maglim2);
	    else
		printf ("%d HST Guide Stars", ng);
	    }
	}
    if (verbose || printhead) {
	if (iraffile)
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
	strcat (filename,".gsc");
	fd = fopen (filename, "w");
	if (fd == NULL) {
	    if (gx) free ((char *)gx);
	    if (gy) free ((char *)gy);
	    if (gm) free ((char *)gm);
	    if (gmb) free ((char *)gmb);
	    if (gra) free ((char *)gra);
	    if (gdec) free ((char *)gdec);
	    if (gnum) free ((char *)gnum);
	    if (gc) free ((char *)gc);
	    free ((char *)wcs);
	    fprintf (stderr, "IMGSC:  cannot write file %s\n", filename);
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

    sprintf (headline,"GSC_NUMBER	RA      	DEC      	MAG   	X    	Y	Type");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    sprintf (headline,"----------	------------	------------	------	-----	-----	----");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (printhead)
	printf ("GSC number RA           Dec           Mag    X      Y   Type\n");

    for (i = 0; i < nbg; i++) {
	if (gx[i] > 0.0 && gy[i] > 0.0) {
	    ra2str (rastr, gra[i], 3);
	    dec2str (decstr, gdec[i], 2);
	    sprintf (headline, "%9.4f	%s	%s	%.2f	%.1f	%.1f	%d",
		 gnum[i], rastr, decstr, gm[i], gx[i], gy[i], gc[i]);
	    if (wfile)
		fprintf (fd, "%s\n", headline);
	    else if (tabout)
		printf ("%s\n", headline);
	    else
		printf ("%9.4f %s %s %6.2f %6.1f %6.1f %3d\n",
			gnum[i], rastr, decstr, gm[i],gx[i],gy[i],gc[i]);
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
 * Sep 17 1996	Set up command line center like IMWCS
 * Oct 16 1996	Rewrite to allow optional new center and use GetWCSFITS
 * Oct 16 1996	Write list of stars to stdout by default, -w to write file
 * Nov 14 1996	Add second magnitude for sort subroutines
 * Nov 15 1996	Change arguments for search routine
 * Dec 10 1996	Change equinox in getfitswcs call to double
 * Dec 12 1996	Allow bright as well as faint magnitude limit
 */
