/* File imtab.c
 * November 15, 1996
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

extern int tabread();
static void ListTab();
static int tabint();
extern void RASortStars();
extern void MagSortStars();
extern void fk524e();
extern void fk425e();
extern void setfk4();
extern void setcenter();
extern void setsecpix();
extern struct WorldCoor *GetFITSWCS();

static int verbose = 0;		/* verbose/debugging flag */
static int wfile = 0;		/* True to print output file */
static double maglim = MAGLIM;	/* Guide Star Catalog magnitude limit */
static char coorsys[4];		/* Output coordinate system */
static int nstars = 0;		/* Number of brightest stars to list */
static int printhead = 0;	/* 1 to print table heading */
static int tabout = 0;		/* 1 for tab table to standard output */
static int rasort = 0;		/* 1 to sort stars by brighness */
static char *tabcat;

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

	case 'b':	/* ouput FK4 (B1950) coordinates */
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

	case 'c':	/* Set tab table catalog name */
	    if (ac < 2)
		usage ();
	    tabcat = *++av;
	    ac--;
	    break;

	case 'h':	/* ouput descriptive header */
	    printhead++;
	    break;

	case 'j':	/* ouput FK5 (J2000) coordinates */
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
	    maglim = atof (*++av);
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

	ListTab (fn);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"Find tab-table catalog stars in FITS or IRAF image files\n");
    fprintf(stderr,"Usage: [-vhstw] [-m mag_off] [-n num] [-p scale]\n");
    fprintf(stderr,"       [-b ra dec] [-j ra dec] FITS or IRAF file(s)\n");
    fprintf(stderr,"  -b: output B1950 (FK4) coordinates (optional center)\n");
    fprintf(stderr,"  -c: Use following tab table catalog \n");
    fprintf(stderr,"  -h: print heading, else do not \n");
    fprintf(stderr,"  -j: output J2000 (FK5) coordinates (optional center)\n");
    fprintf(stderr,"  -m: magnitude limit\n");
    fprintf(stderr,"  -n: number of brightest stars to print \n");
    fprintf(stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -s: sort by RA instead of fltx \n");
    fprintf(stderr,"  -t: tab table to standard output as well as file\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -w: Write tab table output file imagename.cat\n");
    exit (1);
}


static void
ListTab (filename)

char	*filename;	/* FITS or IRAF file filename */

{
    char *header;	/* FITS header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    int iraffile;	/* 1 if IRAF image */
    int *irafheader;	/* IRAF image header */
    double *tnum=0;	/* Tab star numbers */
    double *tra=0;	/* Tab star right ascensions, rads */
    double *tdec=0;	/* Tab star declinations rads */
    double *tm=0;	/* Tab magnitudes */
    double *tmb=0;	/* Unused magnitudes */
    double *tx=0;	/* Tab X positions on image */
    double *ty=0;	/* Tab Y positions on image */
    int *tp=0;		/* Tab plate numbers */
    int nt;		/* Number of Tab stars */
    int nbt;		/* Number of brightest Tab stars actually used */
    int imh, imw;	/* Image height and width in pixels */
    int i, ntmax, nbytes;
    FILE *fd;
    struct WorldCoor *wcs;	/* World coordinate system structure */
    char rastr[16], decstr[16];	/* coordnate strings */
    double cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, mag1, mag2, secpix;
    int offscale, nlog;
    char headline[160];
    char starfile[128];
    char idnum[16];

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

    if (verbose || printhead)

    /* Read world coordinate system information from the image header */
    wcs = GetFITSWCS (header, verbose, &cra, &cdec, &dra, &ddec, &secpix,
		      &imw, &imh, 2000);
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

/* Set magnitude limits for the tab table catalog search */
    if (maglim == 0.0) {
	mag1 = 0.0;
	mag2 = 0.0;
	}
    else {
	mag1 = -2.0;
	mag2 = maglim;
	}

    if (nstars > MAXREF)
	ntmax = nstars;
    else
	ntmax = MAXREF;
    nbytes = ntmax * sizeof (double);
    if (!(tnum = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for tnum\n", nbytes);
    if (!(tra = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for tra\n", nbytes);
    if (!(tdec = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for tdec\n", nbytes);
    if (!(tm = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for tm\n", nbytes);
    if (!(tmb = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for tmb\n", nbytes);
    if (!(tp = (int *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for tp\n", nbytes);
    if (!tnum || !tra || !tdec || !tm || !tp) {
	if (tm) free ((char *)tm);
	if (tra) free ((char *)tra);
	if (tdec) free ((char *)tdec);
	if (tnum) free ((char *)tnum);
	if (tp) free ((char *)tp);
	free ((char *)wcs);
	}
    nlog = 1;

    /* Find the nearby reference stars, in ra/dec */
    nt = tabread (tabcat, cra,cdec,dra,ddec,0.0,mag1,mag2,ntmax,tnum,
		  tra,tdec,tm,tp,nlog);

    /* Project the reference stars into pixels on a plane at ra0/dec0 */
    if (!(tx = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d-bytes for dx\n", nbytes);
    if (!(ty = (double *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d-bytes for ty\n", nbytes);
    if (!tx || !ty) {
	if (tx) free ((char *)tx);
	if (ty) free ((char *)ty);
	if (tm) free ((char *)tm);
	if (tra) free ((char *)tra);
	if (tdec) free ((char *)tdec);
	if (tnum) free ((char *)tnum);
	if (tp) free ((char *)tp);
	free ((char *)wcs);
	return;
	}
    for (i = 0; i < nt; i++ ) {
	offscale = 0;
	if (strcmp (wcs->sysout,"FK4") == 0 || strcmp (wcs->sysout,"fk4") == 0)
	    fk524e (&tra[i], &tdec[i], wcs->epoch);
	wcs2pix (wcs, tra[i], tdec[i], &tx[i], &ty[i], &offscale);
	}

    /* Sort reference stars by brightness (magnitude) */
    MagSortStars (tnum, tra, tdec, tx, ty, tm, tmb, tp, nt);

    if (nstars > 0 && nt > nstars) {
	nbt = nstars;
	if (verbose || printhead)
	    printf ("%d / %d %s stars brighter than %.1f",
		    nbt, nt, tabcat, tm[nbt-1]);
	}
    else {
	nbt = nt;
	if (verbose || printhead) {
	    if (maglim > 0.0)
		printf ("%d %s stars brighter than %.1f",
			nt, tabcat, maglim);
	    else
		printf ("%d %s Catalog Stars", nt, tabcat);
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
	RASortStars (tnum, tra, tdec, tx, ty, tm, tmb, tp, nbt);

    /* Open plate catalog file */
    if (wfile) {
	strcpy (starfile, filename);
	strcat (starfile,".tabstars");
	fd = fopen (starfile, "w");
	if (fd == NULL) {
	    fprintf (stderr, "IMTAB:  cannot write file %s\n", starfile);
	    if (tx) free ((char *)tx);
	    if (ty) free ((char *)ty);
	    if (tm) free ((char *)tm);
	    if (tmb) free ((char *)tmb);
	    if (tra) free ((char *)tra);
	    if (tdec) free ((char *)tdec);
	    if (tnum) free ((char *)tnum);
	    if (tp) free ((char *)tp);
	    free ((char *)wcs);
            return;
	    }
        }

    /* Write header */
    sprintf (headline, "IMAGE	%s", filename);
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    sprintf (headline, "CATALOG	%s",tabcat);
    if (wfile)
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
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (wfile)
    if (tabout)

    sprintf (headline,"ID     	RA      	DEC      	MAG   	X    	Y	Peak");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (tabint (tnum[1]))
	sprintf (idnum,"------");
    else if (tnum[1] < 100.0 && tnum[1] > -10.0)
	sprintf (idnum,"----------");
    else
	sprintf (idnum,"-__---------");

    sprintf (headline,"%s	------------	------------	------	-----	-----	----", idnum);
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (printhead)
	printf (" Number    RA           Dec           Mag    X      Y   Peak\n");

    for (i = 0; i < nbt; i++) {
	if (tx[i] > 0.0 && ty[i] > 0.0) {
	    ra2str (rastr, tra[i], 3);
	    dec2str (decstr, tdec[i], 2);
	    if (tabint (tnum[i]))
		sprintf (idnum, "%d", (int)tnum[i]);
	    else if (tnum[i] < 10.0 && tnum[i] >= 0.0)
		sprintf (idnum, "%9.7f", tnum[i]);
	    else if (tnum[i] < 100.0 && tnum[i] > -10.0)
		sprintf (idnum, "%10.7f", tnum[i]);
	    else if (tnum[i] < 1000.0 && tnum[i] > -100.0)
		sprintf (idnum, "%11.7f", tnum[i]);
	    else
		sprintf (idnum, "%11.4f", tnum[i]);
	    sprintf (headline, "%s	%s	%s	%.2f	%.1f	%.1f	%d",
		 idnum, rastr, decstr, tm[i], tx[i], ty[i], tp[i]);
	    if (wfile)
		fprintf (fd, "%s\n", headline);
	    if (tabout)
		printf ("%s\n", headline);
	    else
		printf ("%s %s %s %6.2f %6.1f %6.1f %d\n",
		    idnum, rastr, decstr, tm[i],tx[i],ty[i],tp[i]);
	    }
	}

    if (wfile)
	fclose (fd);
    if (tx) free ((char *)tx);
    if (ty) free ((char *)ty);
    if (tm) free ((char *)tm);
    if (tmb) free ((char *)tmb);
    if (tra) free ((char *)tra);
    if (tdec) free ((char *)tdec);
    if (tnum) free ((char *)tnum);
    if (tp) free ((char *)tp);
    free ((char *)wcs);

    free (header);
    return;
}


static int
tabint (number)

double number;
{
    int num = (int) (number * 1000000.0);

    if (num % 1 == 0)
	return (1);
    else
	return (0);
}

/* Jul 19 1996	New program
 * Aug 16 1996	Clean up code
 * Aug 27 1996	Drop unused variables after lint
 * Oct 15 1996	Drop unused progname argument
 * Oct 16 1996  Rewrite to allow optional new center and use GetWCSFITS
 * Oct 16 1996  Write list of stars to stdout by default, -w to write file
 * Nov 15 1996	Change arguments in call to TABREAD
 */
