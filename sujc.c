/* File sujc.c
 * November 19, 1996
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

extern int ujcread();
static void ListUJC ();
extern void XSortStars ();
extern void RASortStars ();
extern void MagSortStars ();
extern void fk524e();
extern void setfk4();
extern void setcenter();
extern void setradius();
extern int GetArea();
extern double wcsdist();

static int verbose = 0;		/* Verbose/debugging flag */
static int wfile = 0;		/* True to print output file */
static int plate = 0;           /* If not 0, use only stars from this plate */
static double maglim = MAGLIM;	/* Guide Star Catalog magnitude limit */
static char coorsys[4];		/* Output coordinate system */
static int nstars = 0;		/* Number of brightest stars to list */
static int printhead = 0;	/* 1 to print table heading */
static int tabout = 0;		/* 1 for tab table to standard output */
static int rasort = 0;		/* 1 to sort stars by brighness */
static int distsort = 0;	/* 1 to sort stars by distance from center */
static int debug = 0;		/* True for extra information */
static char *objname = NULL;	/* Object name for output */

main (ac, av)
int ac;
char **av;
{
    char *str, *str1;
    char rastr[16];
    char decstr[16];

    *coorsys = 0;

    if (ac == 1)
        usage ();

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
	    if (*(str+1) || (str1[0] < 47 && str[0] > 58))
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

	case '1':	/* Get closest source */
	    distsort++;
	    nstars = 1;
    	    setradius (600.0);
	    break;

	case 'd':	/* Sort by distance from center */
	    distsort++;
	    break;

	case 'h':	/* ouput descriptive header */
	    printhead++;
	    break;

    	case 'j':	/* center coordinates on command line in J2000 */
	    str1 = *(av+1);
	    if (*(str+1) || (str1[0] < 47 && str[0] > 58))
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
	    maglim = atof (*++av);
	    ac--;
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

    	case 'r':	/* Box radius in arcseconds */
    	    if (ac < 2)
    		usage();
    	    setradius (atof (*++av));
    	    ac--;
    	    break;

	case 's':	/* sort by RA */
	    rasort = 1;
	    break;

	case 't':	/* tab table to stdout */
	    tabout = 1;
	    break;

	case 'u':	/* Object type */
	    if (ac < 2)
		usage ();
	    plate = (int) (atof (*++av) + 0.1);
	    ac--;
	    break;

    	case 'w':	/* write output file */
    	    wfile++;
    	    break;

	default:
	    usage ();
	    break;
	}
    }

    ListUJC ();

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"Find USNO J Catalog stars in a square on the sky\n");
    fprintf (stderr,"Usage: [-1dhstvw] [-m faint mag] [-n num] [-r arcsec] [-b][-j] ra dec\n");
    fprintf(stderr,"  -1: list single closest catalog source\n");
    fprintf(stderr,"  -b: output B1950 (FK4) coordinates around this center\n");
    fprintf(stderr,"  -d: sort by distance from center instead of flux\n");
    fprintf(stderr,"  -h: print heading, else do not \n");
    fprintf(stderr,"  -j: output J2000 (FK5) coordinates around this center\n");
    fprintf(stderr,"  -m: magnitude limit\n");
    fprintf(stderr,"  -n: number of brightest stars to print \n");
    fprintf(stderr,"  -o: object name \n");
    fprintf(stderr,"  -r: search radius in arcsec (default 10)\n");
    fprintf(stderr,"  -s: sort by RA instead of flux \n");
    fprintf(stderr,"  -t: tab table to standard output as well as file\n");
    fprintf(stderr,"  -u: plate number for catalog sources (0=all)\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -w: Write tab table output file imagename.ujc\n");
    exit (1);
}


static void
ListUJC ()

{
    char *header;	/* FITS header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    int iraffile;	/* 1 if IRAF image */
    int *irafheader;	/* IRAF image header */
    double *unum=0;	/* UJC star numbers */
    double *ura=0;	/* UJC star right ascensions, rads */
    double *udec=0;	/* UJC star declinations rads */
    double *um=0;	/* GCS magnitudes */
    double *umb=0;	/* unused magnitudes */
    double *ux=0;	/* UJC X positions on image */
    double *uy=0;	/* UJC Y positions on image */
    int *up=0;		/* UJC object classes */
    int ng;		/* Number of UJC stars */
    int nbg;		/* Number of brightest GCS stars actually used */
    int imh, imw;	/* Image height and width in pixels */
    int i, ngmax, nbytes;
    FILE *fd;
    char rastr[16], decstr[16];	/* coordinate strings */
    double cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, mag1, mag2;
    int offscale, nlog;
    char headline[160];
    char filename[80];

    if (verbose || printhead)

    /* Read world coordinate system information from the image header */
    if (GetArea (verbose, 2000, &cra, &cdec, &dra, &ddec))
	return;

    /* Set limits for search */
    ra1 = cra - dra;
    ra2 = cra + dra;
    dec1 = cdec - ddec;
    dec2 = cdec + ddec;
    if (verbose || printhead) {
	ra2str (rastr, cra, 3);
	dec2str (decstr, cdec, 2);
	if (objname)
	    printf ("%12s %s %s ", objname, rastr, decstr);
	else
	    printf ("USNO J1.0    %s %s ", rastr, decstr);
	if (strcmp (coorsys,"FK4") == 0)
	    printf ("(B1950) ");
	else
	    printf ("(J2000) ");
	printf ("rad=%7.2f\n", ddec*3600.0);
	}

/* Set the magnitude limits for the UJC search */
    if (maglim == 0.0) {
	mag1 = 0.0;
	mag2 = 0.0;
	}
    else {
	mag1 = -2.0;
	mag2 = maglim;
	}

    if (nstars > MAXREF)
	ngmax = nstars;
    else
	ngmax = MAXREF;
    nbytes = ngmax * sizeof (double);
    if (!(unum = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for unum\n", nbytes);
    if (!(ura = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for ura\n", nbytes);
    if (!(udec = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for udec\n", nbytes);
    if (!(um = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for um\n", nbytes);
    if (!(umb = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for umb\n", nbytes);
    if (!(up = (int *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for up\n", nbytes);
    if (!unum || !ura || !udec || !um || !umb || !up) {
	if (um) free ((char *)um);
	if (umb) free ((char *)umb);
	if (ura) free ((char *)ura);
	if (udec) free ((char *)udec);
	if (unum) free ((char *)unum);
	if (up) free ((char *)up);
	return;
	}
    if (verbose)
	nlog = 10;
    else
	nlog = 0;

    /* Find the nearby reference stars, in ra/dec */
    ng = ujcread (cra,cdec,dra,ddec,ddec,mag1,mag2,plate,ngmax,
		  unum,ura,udec,um,up,nlog);

    /* Project the reference stars into pixels on a plane at ra0/dec0 */
    nbytes = ng * sizeof (double);
    if (!(ux = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for ux\n", nbytes);
    if (!(uy = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for uy\n", nbytes);
    if (!ux || !uy) {
	if (um) free ((char *)um);
	if (umb) free ((char *)umb);
	if (ura) free ((char *)ura);
	if (udec) free ((char *)udec);
	if (unum) free ((char *)unum);
	if (up) free ((char *)up);
	if (ux) free ((char *)ux);
	if (uy) free ((char *)uy);
	return;
	}
    for (i = 0; i < ng; i++ ) {
	offscale = 0;
	if (strcmp (coorsys,"FK4") == 0)
	    fk524 (&ura[i],&udec[i]);
	ux[i] = wcsdist (cra, cdec, ura[i], udec[i]);
	uy[i] = 1.0;
	}

    /* Sort reference stars by brightness (magnitude) */
    if (distsort)
	XSortStars (unum, ura, udec, ux, uy, um, umb, up, ng);
    else
	MagSortStars (unum, ura, udec, ux, uy, um, umb, up, ng);

    /* List the brightest or closest MAXSTARS reference stars */
    if (nstars > 0 && ng > nstars) {
	nbg = nstars;
	if (verbose || printhead) {
	    if (distsort)
		printf ("Closest %d / %d UJ Catalog stars (closer than %.2f arcsec)\n",
		     nbg, ng, 3600.0*ux[nbg-1]);
	    else
		printf ("Brightest %d / %d UJ Catalog stars (brighter than %.1f)\n",
		     nbg, ng, um[nbg-1]);
	    }
	}
    else {
	nbg = ng;
	if (verbose || printhead) {
	    if (maglim > 0.0)
		printf ("%d UJ Catalog stars brighter than %.1f\n",
			ng, maglim);
	    else if (verbose)
		printf ("%d UJ Catalog stars\n", ng);
	    }
	}

    /* Sort star-like objects in image by right ascension */
    if (rasort)
	RASortStars (unum, ura, udec, ux, uy, um, umb, up, nbg);

    /* Open plate catalog file */
    if (wfile) {
	if (objname)
	    strcpy (filename,objname);
	else
	    strcpy (filename,"search");
	strcat (filename,".ujc");
	fd = fopen (filename, "w");
	if (fd == NULL) {
	    if (ux) free ((char *)ux);
	    if (uy) free ((char *)uy);
	    if (um) free ((char *)um);
	    if (ura) free ((char *)ura);
	    if (udec) free ((char *)udec);
	    if (unum) free ((char *)unum);
	    if (up) free ((char *)up);
	    fprintf (stderr, "SUJC:  cannot write file %s\n", filename);
            return;
	    }
        }

    /* Write heading */
    sprintf (headline, "CATALOG	UJ1.0");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    ra2str (rastr, cra, 3);
    if (wfile)
	fprintf (fd, "RA	%s\n", rastr);
    if (tabout)
	printf ("RA	%s\n", rastr);
    dec2str (decstr, cdec, 2);
    if (wfile)
	fprintf (fd, "DEC	%s\n", rastr);
    if (tabout)
	printf ("DEC	%s\n", rastr);
    if (strcmp (coorsys,"FK4") == 0)
	sprintf (headline, "EQUINOX	1950.0");
    else
	sprintf (headline, "EQUINOX	2000.0");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (plate > 0) {
        fprintf (fd, "PLATE     %d\n", plate);
        if (tabout)
            printf ("PLATE      %d\n", plate);
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

    sprintf (headline,"UJC_NUMBER	RA      	DEC      	MAG   	Plate	Distance");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    sprintf (headline,"----------	------------	------------	------	----	_______");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (printhead)
	printf (" UJ number    RA           Dec           Mag  Plate Arcsec\n");

    for (i = 0; i < nbg; i++) {
	if (ux[i] > 0.0 && uy[i] > 0.0) {
	    ra2str (rastr, ura[i], 3);
	    dec2str (decstr, udec[i], 2);
	    sprintf (headline, "%12.7f	%s	%s	%.2f	%d	%.2f",
		 unum[i], rastr, decstr, um[i], up[i], 3600.0*ux[i]);
	    if (wfile)
		fprintf (fd, "%s\n", headline);
	    else if (tabout)
		printf ("%s\n", headline);
	    else
		printf ("%12.7f %s %s %6.2f %4d %7.2f\n",
			unum[i], rastr, decstr, um[i],up[i], 3600.0*ux[i]);
	    }
	}

    if (wfile)
	fclose (fd);
    if (ux) free ((char *)ux);
    if (uy) free ((char *)uy);
    if (um) free ((char *)um);
    if (ura) free ((char *)ura);
    if (udec) free ((char *)udec);
    if (unum) free ((char *)unum);
    if (up) free ((char *)up);

    return;
}

/* Oct 18 1996	New prouram based on imujc
 * Nov 13 1996	Use n argument to set maximum number of stars returned
 * Nov 15 1996	Change arguments to ujcread
 * Nov 19 1996	Fix usage()
 */
