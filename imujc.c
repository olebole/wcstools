/* File imujc.c
 * February 21, 1997
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
extern void setfk4();
extern void setcenter();
extern void setsecpix();
extern struct WorldCoor *GetFITSWCS();
extern char *GetFITShead();
extern void fk524e();

static int verbose = 0;		/* Verbose/debugging flag */
static int wfile = 0;		/* True to print output file */
static double maglim1 = MAGLIM1; /* USNO J Catalog bright magnitude limit */
static double maglim2 = MAGLIM;	/* USNO J Catalog faint magnitude limit */
static char coorsys[4];		/* Output coordinate system */
static int nstars = 0;		/* Number of brightest stars to list */
static int printhead = 0;	/* 1 to print table heading */
static int tabout = 0;		/* 1 for tab table to standard output */
static int rasort = 0;		/* 1 to sort stars by brighness */
static int plate = 0;		/* If not 0, use only stars from this plate*/

main (ac, av)
int ac;
char **av;
{
    char *str, *str1;
    char rastr[16], decstr[16];	/* coordinate strings */

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
		setfk4();
		strcpy (rastr, *++av);
		ac--;
		strcpy (decstr, *++av);
		ac--;
		setcenter (rastr, decstr);
		}
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

	case 'u':	/* Plate number to use for objects */
	    if (ac < 2)
		usage ();
	    plate = atoi (*++av);
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

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage ();

    while (ac-- > 0) {
	char *fn = *av++;

	ListUJC (fn);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"Find UJC stars in FITS or IRAF image files\n");
    fprintf (stderr,"Usage: [-vwhst] [-m [mag1] mag2] [-n num] [-u plate]\n");
    fprintf (stderr,"       [-p scale] [-b ra dec] [-j ra dec] FITS or IRAF file(s)\n");
    fprintf (stderr,"       [-b ra dec] [-j ra dec] file.fts ...\n");
    fprintf (stderr,"  -b: output B1950 (FK4) coordinates (optional center)\n");
    fprintf (stderr,"  -c: Use following RA and Dec as center \n");
    fprintf (stderr,"  -h: print heading, else do not \n");
    fprintf (stderr,"  -j: output J2000 (FK5) coordinates (optional center)\n");
    fprintf (stderr,"  -m: magnitude limit(s)\n");
    fprintf (stderr,"  -n: number of brightest stars to print \n");
    fprintf (stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf (stderr,"  -s: sort by RA instead of flux \n");
    fprintf (stderr,"  -t: tab table to standard output as well as file\n");
    fprintf (stderr,"  -u: plate number for catalog sources (0=all)\n");
    fprintf (stderr,"  -v: verbose\n");
    fprintf (stderr,"  -w: Write tab table output file imagename.ujc\n");
    exit (1);
}


static void
ListUJC (filename)

char	*filename;	/* FITS or IRAF file filename */

{
    char *header;	/* FITS header */
    double *unum=0;	/* UJC star numbers */
    double *ura=0;	/* UJC star right ascensions, rads */
    double *udec=0;	/* UJC star declinations rads */
    double *um=0;	/* UJC magnitudes */
    double *umb=0;	/* Unused magnitudes */
    double *ux=0;	/* UJC X positions on image */
    double *uy=0;	/* UJC Y positions on image */
    int *up=0;		/* UJC plate numbers */
    int nu;		/* Number of UJC stars */
    int nbu;		/* Number of brightest GCS stars actually used */
    int imw, imh;	/* Image width and height in pixels */
    int i, numax, nbytes;
    FILE *fd;
    struct WorldCoor *wcs;	/* World coordinate system structure */
    char rastr[16], decstr[16];	/* coordinate strings */
    double cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, mag1, mag2, secpix;
    double mag;
    int offscale, nlog;
    char headline[160];

    if (verbose || printhead)

    /* Read world coordinate system information from the image header */
    if ((header = GetFITShead (filename)) == NULL)
	return;
    wcs = GetFITSWCS (header, verbose, &cra, &cdec, &dra, &ddec, &secpix,
                &imw, &imh, 2000.0);
    free (header);
    if (nowcs (wcs))
        return;

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

/* Set magnitude limits for the UJ Catalog search */
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

    numax = MAXREF;
    nbytes = MAXREF * sizeof (double);
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
	free ((char *)wcs);
	return;
	}
    nlog = 0;

    /* Find the nearby reference stars, in ra/dec */
    nu = ujcread (cra,cdec,dra,ddec,0.0,mag1,mag2,plate,numax,
		  unum,ura,udec,um,up,nlog);

    /* Project the reference stars into pixels on a plane at ra0/dec0 */
    nbytes = nu * sizeof (double);
    ux = (double *) malloc (nbytes);
    uy = (double *) malloc (nbytes);
    if (!ux)
	fprintf (stderr, "Could not malloc %d bytes for ux\n", nbytes);
    if (!ux)
	fprintf (stderr, "Could not malloc %d bytes for uy\n", nbytes);
    if (!ux || !uy) {
	if (ux) free ((char *)ux);
	if (uy) free ((char *)uy);
	if (um) free ((char *)um);
	if (umb) free ((char *)umb);
	if (ura) free ((char *)ura);
	if (udec) free ((char *)udec);
	if (unum) free ((char *)unum);
	if (up) free ((char *)up);
	free ((char *)wcs);
	return;
	}
    for (i = 0; i < nu; i++ ) {
	offscale = 0;
	if (strcmp (wcs->sysout,"FK4") == 0 || strcmp (wcs->sysout,"fk4") == 0)
	    fk524e (&ura[i], &udec[i], wcs->epoch);
	wcs2pix (wcs, ura[i], udec[i], &ux[i], &uy[i], &offscale);
	}

    /* Sort reference stars by brightness (magnitude) */
    MagSortStars (unum, ura, udec, ux, uy, um, umb, up, nu);

    if (nstars > 0 && nu > nstars) {
	nbu = nstars;
	if (verbose || printhead)
	    if (maglim1 > 0.0)
		printf ("%d / %d UJ Catalog Stars between %.1f and %.1f",
			nbu, nu, um[0], um[nbu-1]);
	    else
		printf ("%d / %d UJ Catalog Stars brighter than %.1f",
			nbu, nu, um[nbu-1]);
	}
    else {
	nbu = nu;
	if (verbose || printhead) {
	    if (maglim1 > 0.0)
		printf ("%d UJ Catalog Stars between %.1f and %.1f",
		    nu, maglim1, maglim2);
	    else if (maglim2 > 0.0)
		printf ("%d UJ Catalog Stars brighter than %.1f", nu, maglim2);
	    else
		printf ("%d UJ Catalog Stars", nu);
	    }
	}
    if (verbose || printhead) {
	if (strsrch (filename,".imh") != NULL)
	    printf (" in IRAF image %s\n",filename);
	else
	    printf (" in FITS image %s\n", filename);
	}
    sprintf (headline, "IMAGE	%s", filename);

    /* Sort star-like objects in image by right ascension */
    if (rasort)
	RASortStars (unum, ura, udec, ux, uy, um, umb, up, nbu);

    /* Open plate catalog file */
    if (wfile) {
	strcat (filename,".ujc");
	fd = fopen (filename, "w");
	if (fd == NULL) {
	    fprintf (stderr, "IMUJC:  cannot write file %s\n", filename);
	    if (ux) free ((char *)ux);
	    if (uy) free ((char *)uy);
	    if (um) free ((char *)um);
	    if (umb) free ((char *)umb);
	    if (ura) free ((char *)ura);
	    if (udec) free ((char *)udec);
	    if (unum) free ((char *)unum);
	    if (up) free ((char *)up);
	    free ((char *)wcs);
            return;
            }
        }

    /* Write header */
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    sprintf (headline, "CATALOG	UJ1.0");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (strcmp (wcs->sysout,"FK4") == 0 || strcmp (wcs->sysout,"fk4") == 0)
	sprintf (headline, "EQUINOX	1950.0");
    else
	sprintf (headline, "EQUINOX	2000.0");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (plate > 0) {
	if (wfile)
	    fprintf (fd, "PLATE	%d\n", plate);
	if (tabout)
	    printf ("PLATE	%d\n", plate);
	}
    if (rasort) {
	if (wfile)
	    fprintf (fd, "RASORT	T\n");
	if (tabout)
	    printf ("RASORT	T\n");
	}

    if (wfile)
    if (tabout)

    sprintf (headline,"UJC_NUMBER	RA      	DEC      	MAG   	X    	Y	Plate");
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
	printf (" UJ number    RA           Dec           Mag    X      Y   Plate\n");

    for (i = 0; i < nbu; i++) {
	if (ux[i] > 0.0 && uy[i] > 0.0) {
	    ra2str (rastr, ura[i], 3);
	    dec2str (decstr, udec[i], 2);
	    sprintf (headline, "%12.7f	%s	%s	%.2f	%.1f	%.1f	%d",
		 unum[i], rastr, decstr, um[i], ux[i], uy[i], up[i]);
	    if (wfile)
		fprintf (fd, "%s\n", headline);
	    else if (tabout)
		printf ("%s\n", headline);
	    else
		printf ("%12.7f %s %s %6.2f %6.1f %6.1f %4d\n",
		    unum[i], rastr, decstr, um[i],ux[i],uy[i],up[i]);
	    }
	}

    if (wfile)
	fclose (fd);
    if (ux) free ((char *)ux);
    if (uy) free ((char *)uy);
    if (um) free ((char *)um);
    if (umb) free ((char *)umb);
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
 * Aug 15 1996	Clean up image header reading code
 * Aug 27 1996	Drop unused variables after lint
 * Sep 17 1996	Set up command line center like IMWCS
 * Sep 30 1996	Declare undeclared variables RASTR and DECSTR in main()
 * Oct 16 1996  Rewrite to allow optional new center and use GetWCSFITS
 * Oct 16 1996  Write list of stars to stdout by default, -w to write file
 * Oct 18 1996	Write selected plate to tab table header
 * Nov 15 1996	Change arguments to UJCREAD
 * Dec 10 1996	Change equinox in getfitswcs call to double
 * Dec 12 1996	Allow bright as well as faint magnitude limit
 * Dec 13 1996	Fix bug writing sorted header
 *
 * Jan 10 1997	Fix bug in RASort Stars which did not sort magnitudes
 * Feb 21 1997  Get FITS header from subroutine
 */
