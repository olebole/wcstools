/* File sujc.c
 * November 12, 1997
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
extern int ujcrnum();
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
static double maglim1 = MAGLIM1; /* Guide Star Catalog magnitude limit */
static double maglim2 = MAGLIM;	/* Guide Star Catalog magnitude limit */
static char coorsys[4];		/* Input coordinate system */
static char coorout[4];		/* Output coordinate system */
static int nstars = 0;		/* Number of brightest stars to list */
static int printhead = 0;	/* 1 to print table heading */
static int tabout = 0;		/* 1 for tab table to standard output */
static int rasort = 0;		/* 1 to sort stars by brighness */
static int distsort = 0;	/* 1 to sort stars by distance from center */
static int debug = 0;		/* True for extra information */
static char *objname = NULL;	/* Object name for output */
static char rastr[16];
static char decstr[16];

main (ac, av)
int ac;
char **av;
{
    char *str, *str1;
    double *ujcnum;
    int nfind = 0;

    ujcnum = NULL;
    *coorsys = 0;
    *coorout = 0;

    if (ac == 1)
        usage ();

    /* crack arguments */
    for (av++; --ac > 0; av++) {

	/* Set RA, Dec, and equinox if WCS-generated argument */
	if (strsrch (*av,":") != NULL) {
	    if (ac < 3)
		usage();
	    else {
		strcpy (rastr, *av);
		ac--;
		strcpy (decstr, *++av);
		setcenter (rastr, decstr);
		ac--;
		strcpy (coorsys, *++av);
		if (coorsys[0] == 'B' || coorsys[0] == 'b' ||
		    !strcmp (coorsys,"FK4") || !strcmp (coorsys,"fk4")) {
		    setfk4();
		    strcpy (coorsys, "FK4");
		    }
		else
		    strcpy (coorsys, "FK5");
		}
	    }

	/* Set star number if getting just one star */
	else if (strsrch (*av,".") != NULL) {
	    if (ujcnum == NULL)
		ujcnum = (double *) malloc (sizeof(double));
	    else
		ujcnum = (double *) realloc (ujcnum, (nfind+1)*sizeof(double));
	    ujcnum[nfind] = atof (*av);
	    nfind++;
	    }

	else if (*(str = *av) == '-') {
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
	}

    ListUJC (nfind, ujcnum);

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"Find USNO J Catalog stars in a square on the sky\n");
    fprintf (stderr,"Usage: [-1dhstvw] [-m [mag1] mag2] [-n num] [-r arcsec] [-b][-j] ra dec\n");
    fprintf(stderr,"  -1: list single closest catalog source\n");
    fprintf(stderr,"  -b: output B1950 (FK4) coordinates around this center\n");
    fprintf(stderr,"  -d: sort by distance from center instead of flux\n");
    fprintf(stderr,"  -h: print heading, else do not \n");
    fprintf(stderr,"  -j: output J2000 (FK5) coordinates around this center\n");
    fprintf(stderr,"  -m: magnitude limit(s)\n");
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
ListUJC (nfind, ujcnum)

int nfind;		/* Number of stars to find */
double *ujcnum;		/* USNO J Catalog zone.number */

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
    double cra = -1.0;
    double cdec = -99.0;
    double ddec = 0.0;
    double dra, ra1, ra2, dec1, dec2, mag1, mag2, mag;
    int offscale, nlog;
    char headline[160];
    char filename[80];

    if (verbose || printhead) {
	if (nstars == 1 && nfind == 0)
	else
	}
    if (verbose)
	nlog = 10;
    else
	nlog = 0;

    /* Find stars specified by number */
    if (ujcnum != NULL) {
	nbytes = nfind * sizeof (double);
	if (!(unum = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for unum\n", nbytes);
	else {
	    for (i = 0; i < nfind; i++)
		unum[i] = ujcnum[i];
	    }
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
	if (!(ux = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for ux\n", nbytes);
	if (!(uy = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for uy\n", nbytes);
	if (!unum || !ura || !udec || !um || !umb || !up || !ux || !uy) {
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
	wfile = 0;
	nbg = ujcrnum (nfind, unum, ura, udec, um, umb, up, nlog);
	for (i = 0; i < nbg; i++ ) {
	    if (strcmp (coorout,"FK4") == 0)
		fk524 (&ura[i],&udec[i]);
	    ux[i] = 0.0;
	    uy[i] = 1.0;
	    }
	}

    /* Find stars specified by location */
    else {

    /* Set limits from defaults and command line information */
    if (GetArea (verbose, 2000, &cra, &cdec, &dra, &ddec))
	return;

    if (coorout[0] == 0)
	strcpy (coorout, coorsys);
    if (strcmp (coorout, coorsys)) {
	if (strcmp (coorout, "FK5") == 0) {
	    ra2str (rastr, cra, 3);
	    dec2str (decstr, cdec, 2);
	    }
	else {
	    ra1 = cra;
	    dec1 = cdec;
	    fk524 (&ra1,&dec1);
	    ra2str (rastr, ra1, 3);
	    dec2str (decstr, dec1, 2);
	    }
	}
    if (verbose || printhead) {
	if (objname)
	    printf ("%12s %s %s ", objname, rastr, decstr);
	else
	    printf ("USNO J 1.0   %s %s ", rastr, decstr);
	if (strcmp (coorout,"FK4") == 0)
	    printf ("(B1950)  ");
	else
	    printf ("(J2000)  ");
	printf ("radius = %7.2f\n", ddec*3600.0);
	}

/* Set the magnitude limits for the UJC search */
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
	    if (distsort) {
		if (nbg > 1)
		    printf ("Closest %d / %d UJ Catalog stars (closer than %.2f arcsec)\n",
			     nbg, ng, 3600.0*ux[nbg-1]);
		else
		    printf ("Closest of %d UJ Catalog stars\n", ng);
		}
	    else if (mag1 > 0.0)
		printf ("Brightest %d / %d UJ Catalog stars (between %.1f and %.1f)\n",
		     nbg, ng, um[0], um[nbg-1]);
	    else
		printf ("Brightest %d / %d UJ Catalog stars (brighter than %.1f)\n",
		     nbg, ng, um[nbg-1]);
	    }
	}
    else {
	nbg = ng;
	if (verbose || printhead) {
	    if (maglim1 > 0.0)
		printf ("%d UJ Catalog stars between %.1f and %.1f\n",
			ng, maglim1, maglim2);
	    else if (maglim2 > 0.0)
		printf ("%d UJ Catalog stars brighter than %.1f\n",
			ng, maglim2);
	    else if (verbose)
		printf ("%d UJ Catalog stars\n", ng);
	    }
	}

    /* Sort star-like objects in image by right ascension */
    if (rasort)
	RASortStars (unum, ura, udec, ux, uy, um, umb, up, nbg);
    }

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

    if (cra >= 0.0) {
	ra2str (rastr, cra, 3);
	if (wfile)
	    fprintf (fd, "RA	%s\n", rastr);
	if (tabout)
	    printf ("RA	%s\n", rastr);
	}

    if (cdec >= -90.0) {
	dec2str (decstr, cdec, 2);
	if (wfile)
	fprintf (fd, "DEC	%s\n", decstr);
	if (tabout)
	    printf ("DEC	%s\n", decstr);
	}

    if (ddec > 0.0) {
	dec2str (decstr, ddec, 2);
	if (wfile)
	    fprintf (fd, "RADIUS	%s\n", decstr);
	if (tabout)
	    printf ("RADIUS	%s\n", decstr);
	}

    if (strcmp (coorsys,"FK4") == 0)
	sprintf (headline, "EQUINOX	1950.0");
    else
	sprintf (headline, "EQUINOX	2000.0");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (plate > 0) {
	if (wfile)
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

    if (strcmp (coorout,"FK4") == 0)
	sprintf (headline,"UJC_NUMBER	RA1950  	DEC1950  	MAG   	PLATE	ARCSEC");
    else
	sprintf (headline,"UJC_NUMBER	RA2000  	DEC2000  	MAG   	PLATE	ARCSEC");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    sprintf (headline,"----------	------------	------------	------	----	------");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (printhead) {
	if (nbg == 0)
	    printf ("No UJ Stars Found\n");
	else if (strcmp (coorout,"FK4") == 0)
	    printf (" UJ number    RA1950       Dec1950       Mag  Plate Arcsec\n");
	else
	    printf (" UJ number    RA2000       Dec2000       Mag  Plate Arcsec\n");
	}

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
 * Dec 12 1996	Allow bright as wella s faint magnitude limits
 * Dec 13 1996	Fix bug writing plate into header
 * Dec 18 1996	Allow WCS sky coordinate format as input argument for center
 * Dec 18 1996	Add option to print entries for specified catalog numbers
 * Dec 30 1996	Clean up closest star message
 * Dec 30 1996	Print message instead of heading if no stars are found
 * Jan 10 1997	Fix bug in RASort Stars which did not sort magnitudes
 * Mar 12 1997	Fix header alignment
 * Nov 12 1997	Fix DEC in header to print Dec string instead of RA string
 */
