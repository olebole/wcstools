/* File imtab.c
 * August 27, 1996
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

extern int tabread();
static void ListTab();
static int tabint();
extern void RASortStars();
extern void MagSortStars();
extern void fk524e();
extern void fk425e();

static double maglim = 0.0;
static int rasort = 0;
static int printhead = 0;
static char *tabcat;
static int tabout = 0;
static char coorsys[4];
static double cra0 = 0.0;
static double cdec0 = 0.0;

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    char *str1;
    int nstar = 0;

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
		strcpy (coorsys, "FK4");
	    else if (ac < 3)
		usage (progname);
	    else {
		cra0 = str2ra (*++av);
		ac--;
		cdec0 = str2dec (*++av);
		ac--;
		}
	    break;
	case 'c':	/* Set tab table catalog name */
	    if (ac < 2)
		usage (progname);
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
		usage (progname);
	    else {
		strcpy (coorsys, "FK5");
		cra0 = str2ra (*++av);
		ac--;
		cdec0 = str2dec (*++av);
		ac--;
		}
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

	ListTab (fn, nstar);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"Find tab-table catalog stars in FITS or IRAF image files\n");
    fprintf(stderr,"IMTAB: usage: [-v] [-m mag_off] [-n num] file.fts ...\n",
	    progname);
    fprintf(stderr,"  -b: output B1950 (FK4) coordinates (optional center)\n");
    fprintf(stderr,"  -c: Use following tab table catalog \n");
    fprintf(stderr,"  -h: print heading, else do not \n");
    fprintf(stderr,"  -j: output J2000 (FK5) coordinates (optional center)\n");
    fprintf(stderr,"  -m: magnitude limit\n");
    fprintf(stderr,"  -n: number of brightest stars to print \n");
    fprintf(stderr,"  -s: sort by RA instead of fltx \n");
    fprintf(stderr,"  -t: tab table to standard output as well as file\n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


extern int findStars ();
struct WorldCoor *wcsinit();	
extern int pix2wcst();

static void
ListTab (filename, nstars)

char	*filename;	/* FITS or IRAF file filename */
int	nstars;		/* Number of brightest stars to list */

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
    double *tx, *ty;	/* Tab positions on image */
    int *tp;		/* Tab plate numbers */
    int nt;		/* Number of Tab stars */
    int nbt;		/* Number of brightest Tab stars actually used */
    int i, ntmax, nbytes;
    FILE *fd;
    struct WorldCoor *wcs;	/* World coordinate system structure */
    char rastr[16], decstr[16];	/* coordnate strings */
    double cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, mag1, mag2;
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

    if (verbose && printhead)

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
	fk425e (&cra, &cdec, wcs->epoch);
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

/* Set magnitude limits for the tab table catalog search */
    if (maglim == 0.0) {
	mag1 = 0.0;
	mag2 = 0.0;
	}
    else {
	mag1 = -2.0;
	mag2 = maglim;
	}

    ntmax = MAXREF;
    nbytes = MAXREF * sizeof (double);
    tnum = (double *) malloc (nbytes);
    tra = (double *) malloc (nbytes);
    tdec = (double *) malloc (nbytes);
    tm = (double *) malloc (nbytes);
    tp = (int *) malloc (nbytes);
    nlog = 1;

    /* Find the nearby reference stars, in ra/dec */
    nt = tabread (tabcat, ra1,ra2,dec1,dec2,mag1,mag2,ntmax,tnum,
		  tra,tdec,tm,tp,nlog);

    /* Project the reference stars into pixels on a plane at ra0/dec0 */
    tx = (double *) malloc (nt * sizeof (double));
    ty = (double *) malloc (nt * sizeof (double));
    if (!tx || !ty) {
	fprintf (stderr, "Could not malloc temp space of %d bytes\n",
					    nt*sizeof(double)*2);
	return;
	}
    for (i = 0; i < nt; i++ ) {
	offscale = 0;
	if (strcmp (wcs->sysout,"FK4") == 0 || strcmp (wcs->sysout,"fk4") == 0)
	    fk524e (&tra[i], &tdec[i], wcs->epoch);
	wcs2pix (wcs, tra[i], tdec[i], &tx[i], &ty[i], &offscale);
	}

    /* Sort reference stars by brightness (magnitude) */
    MagSortStars (tnum, tra, tdec, tx, ty, tm, tp, nt);

    if (nstars > 0 && nt > nstars) {
	nbt = nstars;
	if (verbose && printhead)
	    printf ("using %d / %d %s stars brighter than %.1f",
		    nbt, nt, tabcat, tm[nbt-1]);
	}
    else {
	nbt = nt;
	if (verbose && printhead) {
	    if (maglim > 0.0)
		printf ("%d %s stars brighter than %.1f",
			nt, tabcat, maglim);
	    else
		printf ("%d %s Catalog Stars", nt, tabcat);
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
	RASortStars (tnum, tra, tdec, tx, ty, tm, tp, nbt);

    /* Open plate catalog file */
    strcpy (starfile, filename);
    strcat (starfile,".tabstars");
    fd = fopen (starfile, "w");
    if (fd == NULL) {
	fprintf (stderr, "IMTAB:  cannot write file %s %d\n", starfile, fd);
        return;
        }

    /* Write header */
    sprintf (headline, "IMAGE	%s", filename);
    fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    sprintf (headline, "CATALOG	%s",tabcat);
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

    if (tabout)

    sprintf (headline,"ID     	RA      	DEC      	MAG   	X    	Y	Peak");
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
    fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

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
	    fprintf (fd, "%s\n", headline);
	    if (tabout)
		printf ("%s\n", headline);
	    else if (verbose)
		printf ("%s %s %s %6.2f %6.1f %6.1f %d\n",
		    idnum, rastr, decstr, tm[i],tx[i],ty[i],tp[i]);
	    }
	}

    fclose (fd);
    if (tx) free ((char *)tx);
    if (ty) free ((char *)ty);
    if (tm) free ((char *)tm);
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
 */
