/* File scat.c
 * January 10, 1997
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

#define GSC	1	/* refcat value for HST Guide Star Catalog */
#define UJC	2	/* refcat value for USNO UJ Star Catalog */
#define UAC	3	/* refcat value for USNO A Star Catalog */

#define MAXREF 100

static void usage();

extern int gscread();
static int ListCat ();
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
static int classd = -1;		/* Guide Star Catalog object classes */
static int uplate = 0;		/* UJ Catalog plate number to use */
static double maglim1 = MAGLIM1; /* Catalog bright magnitude limit */
static double maglim2 = MAGLIM;	/* Catalog faint magnitude limit */
static char coorsys[4];		/* Input coordinate system */
static char coorout[4];		/* Output coordinate system */
static int nstars = 0;		/* Number of brightest stars to list */
static int printhead = 0;	/* 1 to print table heading */
static int tabout = 0;		/* 1 for tab table to standard output */
static int rasort = 0;		/* 1 to sort stars by brighness */
static int distsort = 0;	/* 1 to sort stars by distance from center */
static int debug = 0;		/* True for extra information */
static char *objname = NULL;	/* Object name for output */
static int refcat = GSC;	/* reference catalog switch */
static char refcatname[32]="GSC";	/* reference catalog name */
static int refcat2 = 0;		/* Second reference catalog switch */

main (ac, av)
int ac;
char **av;
{
    char *str, *str1;
    char rastr[16];
    char decstr[16];
    int nfind = 0;
    double *snum;
    snum = NULL;

    *coorsys = 0;
    *coorout = 0;

    if (ac == 1)
        usage ();

    /* crack arguments */
    for (av++; --ac > 0; av++) {

	/* Set RA, Dec, and equinox if WCS-generated argument */
	if (strsrch (*av,":")) {
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
	    if (snum == NULL)
		snum = (double *) malloc (sizeof(double));
	    else
		snum = (double *) realloc (snum, (nfind+1)*sizeof(double));
	    snum[nfind] = atof (*av);
	    nfind++;
	    }

	/* Otherwise, read command */
	else if (*(str = *av) == '-') {
	    char c;
	    while (c = *++str)
	    switch (c) {
	    case 'v':	/* more verbosity */
		verbose++;
		break;
    	    case 'b':	/* initial coordinates on command line in B1950 */
		if (strlen (coorsys) != 0)
		    strcpy (coorout, "FK4");
		else {
		    strcpy (coorsys, "FK4");
		    strcpy (coorout, "FK4");
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
		    }
    		break;

	    case '1':	/* Get closest source */
		distsort++;
		nstars = 1;
    		setradius (60.0);
		break;

	    case 'c':       /* Set reference catalog */
		if (ac < 2)
		    usage();
		strcpy (refcatname, *++av);
		if (strncmp(refcatname,"gs",2)==0 ||
		    strncmp (refcatname,"GS",2)== 0)
		    refcat = GSC;
		else if (strncmp(refcatname,"ua",2)==0 ||
		     strncmp(refcatname,"UA",2)==0)
		    refcat = UAC;
		else if (strncmp(refcatname,"uj",2)==0 ||
		     strncmp(refcatname,"UJ",2)==0)
		    refcat = UJC;
		else
		    refcat = 0;
		ac--;
		break;

	    case 'd':	/* Sort by distance from center */
		distsort++;
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
		if (strlen (coorsys) != 0)
		    strcpy (coorout, "FK5");
		else {
		    strcpy (coorout, "FK5");
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

	    case 'u':       /* UJ Catalog plate number */
		if (ac < 2)
		    usage();
		uplate = (int) atof (*++av);
		ac--;
		break;

    	    case 'w':	/* write output file */
    		wfile++;
    		break;

    	    case '2':	/* check second catalog */
    		refcat2 = UAC;
    		break;

	    default:
		usage ();
		break;
	    }
	}
    }

    if (refcat2 > 0 && ListCat (nfind, snum) < 1) {
	refcat = refcat2;
	(void)ListCat (nfind, snum);
	}

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"Find catalog stars in a square on the sky\n");
    fprintf (stderr,"Usage: [-1dhstvw] [-m [mag1] mag2] [-n num] [-r arcsec] [-b][-j] ra dec\n");
    fprintf(stderr,"  -1: list single closest catalog source\n");
    fprintf(stderr,"  -b: output B1950 (FK4) coordinates around this center\n");
    fprintf(stderr,"  -c: reference catalog (gsc, ujc, or tab table file\n");
    fprintf(stderr,"  -d: sort by distance from center instead of flux\n");
    fprintf(stderr,"  -g: object type (0=stars 3=galaxies -1=all)\n");
    fprintf(stderr,"  -h: print heading, else do not \n");
    fprintf(stderr,"  -j: output J2000 (FK5) coordinates around this center\n");
    fprintf(stderr,"  -m: magnitude limit(s)\n");
    fprintf(stderr,"  -n: number of brightest stars to print \n");
    fprintf(stderr,"  -o: object name \n");
    fprintf(stderr,"  -r: search radius in arcsec (default 10)\n");
    fprintf(stderr,"  -s: sort by RA instead of flux \n");
    fprintf(stderr,"  -t: tab table to standard output as well as file\n");
    fprintf(stderr,"  -u: UJ catalog single plate number to accept\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -w: Write tab table output file imagename.gsc\n");
    exit (1);
}


static int
ListCat (nfind, snum)

int nfind;		/* Number of stars to find */
double *snum;		/* Catalog numbers */

{
    char *header;	/* FITS header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    int iraffile;	/* 1 if IRAF image */
    int *irafheader;	/* IRAF image header */
    double *gnum=0;	/* Catalog star numbers */
    double *gra=0;	/* Catalog star right ascensions, rads */
    double *gdec=0;	/* Catalog star declinations rads */
    double *gm=0;	/* Catalog magnitudes */
    double *gmb=0;	/* Catalog B magnitudes */
    double *gx=0;	/* Catalog star X positions on image */
    double *gy=0;	/* Catalog star Y positions on image */
    int *gc=0;		/* Catalog star object classes */
    int ng;		/* Number of catalog stars */
    int nbg;		/* Number of brightest catalog stars actually used */
    int imh, imw;	/* Image height and width in pixels */
    int i, ngmax, nbytes;
    FILE *fd;
    char rastr[16], decstr[16];	/* coordinate strings */
    double cra = -1.0;
    double cdec = -99.0;
    double dra, ddec, ra1, ra2, dec1, dec2, mag1, mag2, box;
    double mag;
    int offscale, nlog;
    char headline[160];
    char filename[80];
    char title[80];

    if (verbose || printhead) {
	if (nstars == 1)
	else
	}
    if (debug)
	nlog = 1;
    else if (verbose)
	nlog = 100;
    else
	nlog = 0;

    if (refcat == GSC)
	strcpy (title, "HST Guide Stars");
    else if (refcat == UAC)
	strcpy (title, "USNO A Catalog Stars");
    else if (refcat == UJC)
	strcpy (title, "USNO J Catalog Stars");
    else if (refcatname[0] > 0)
	sprintf (title, "%s Catalog Stars", refcatname);
    else {
	fprintf (stderr,"ListCat: Catalog name is missing\n");
	return (0);
	}

    /* Find stars specified by number */
    if (snum != NULL) {
	nbytes = nfind * sizeof (double);
	if (!(gnum = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gnum\n", nbytes);
	else {
	    for (i = 0; i < nfind; i++)
		gnum[i] = snum[i];
	    }
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
	if (!(gx = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gx\n", nbytes);
	if (!(gy = (double *) malloc (nbytes)))
	    fprintf (stderr, "Could not malloc %d bytes for gy\n", nbytes);
	if (!gnum || !gra || !gdec || !gm || !gmb || !gc || !gx || !gy) {
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
	wfile = 0;

	/* Find the specified catalog stars */
	if (refcat == GSC)
	    nbg = gscrnum (nfind,gnum,gra,gdec,gm,gc,nlog);
	else if (refcat == UAC)
	    nbg = uacrnum (nfind,gnum,gra,gdec,gm,gmb,gc,nlog);
	else if (refcat == UJC)
	    nbg = ujcrnum (nfind,gnum,gra,gdec,gm,gc,nlog);
	else
	    nbg = tabrnum (refcatname,nfind,gnum,gra,gdec,gm,gc,debug);
	for (i = 0; i < nbg; i++ ) {
	    if (strcmp (coorout,"FK4") == 0)
		fk524 (&gra[i],&gdec[i]);
	    gx[i] = 0.0;
	    gy[i] = 1.0;
	    }
	}

    /* Find stars specified by location */
    else {


    /* Set limits from defaults and command line information */
    if (GetArea (verbose,2000,&cra,&cdec,&dra,&ddec))
	return (0);

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
	else {
	    if (refcat == GSC)
		printf ("HST GSC   %s %s  ", rastr, decstr);
	    else if (refcat == UAC)
		printf ("USNO A 1.0     %s %s  ", rastr, decstr);
	    else if (refcat == UJC)
		printf ("USNO UJ1.0    %s %s  ", rastr, decstr);
	    else
		printf ("%9s %s %s  ", refcatname, rastr, decstr);
	    }
	if (strcmp (coorout,"FK4") == 0)
	    printf ("(B1950)  ");
	else
	    printf ("(J2000)  ");
	printf ("radius = %7.2f\n", ddec*3600.0);
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
    if (!gm)
	fprintf (stderr, "Could not malloc %d bytes for gm\n", nbytes);
    if (!(gc = (int *) malloc (nbytes)))
	fprintf (stderr, "Could not malloc %d bytes for gc\n", nbytes);
    if (!gnum || !gra || !gdec || !gm || !gmb || !gc) {
	if (gm) free ((char *)gm);
	if (gmb) free ((char *)gmb);
	if (gra) free ((char *)gra);
	if (gdec) free ((char *)gdec);
	if (gnum) free ((char *)gnum);
	if (gc) free ((char *)gc);
	return (0);
	}

    /* Find the nearby reference stars, in ra/dec */
    if (refcat == GSC)
	ng = gscread (cra,cdec,dra,ddec,ddec,mag1,mag2,classd,ngmax,
		      gnum,gra,gdec,gm,gc,nlog);
    else if (refcat == UAC)
	ng = uacread (cra,cdec,dra,ddec,ddec,mag1,mag2,uplate,ngmax,
		      gnum,gra,gdec,gm,gmb,gc,debug);
    else if (refcat == UJC)
	ng = ujcread (cra,cdec,dra,ddec,ddec,mag1,mag2,uplate,ngmax,
		      gnum,gra,gdec,gm,gc,debug);
    else
	ng = tabread (refcatname,cra,cdec,dra,ddec,ddec,mag1,mag2,ngmax,
		      gnum,gra,gdec,gm,gc,debug);

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
    for (i = 0; i < ng; i++ ) {
	offscale = 0;
	if (strcmp (coorsys,"FK4") == 0)
	    fk524 (&gra[i],&gdec[i]);
	gx[i] = wcsdist (cra, cdec, gra[i], gdec[i]);
	gy[i] = 1.0;
	}

    /* Sort reference stars by brightness (magnitude) */
    if (distsort)
	XSortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, ng);
    else
	MagSortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, ng);

    /* List the brightest or closest MAXSTARS reference stars */
    if (nstars > 0 && ng > nstars) {
	nbg = nstars;
	if (verbose || printhead) {
	    if (distsort) {
		if (ng > 1)
		    printf ("Closest %d / %d %s (closer than %.2f arcsec)\n",
			    nbg, ng, title, 3600.0*gx[nbg-1]);
		else
		    printf ("Closest of %d %s\n",ng, title);
		}
	    else if (maglim1 > 0.0)
		printf ("%d / %d %s (between %.1f and %.1f)\n",
		     nbg, ng, title, gm[0], gm[nbg-1]);
	    else
		printf ("%d / %d %s (brighter than %.1f)\n",
		     nbg, ng, title, gm[nbg-1]);
	    }
	}
    else {
	nbg = ng;
	if (verbose || printhead) {
	    if (maglim1 > 0.0)
		printf ("%d %s between %.1f and %.1f\n",
			ng, title, maglim1, maglim2);
	    else if (maglim2 > 0.0)
		printf ("%d %s brighter than %.1f\n",
			ng, title, maglim2);
	    else if (verbose)
		printf ("%d %s\n", ng, title);
	    }
	}

    /* Sort star-like objects in image by right ascension */
    if (rasort)
	RASortStars (gnum, gra, gdec, gx, gy, gm, gmb, gc, nbg);
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
	else if (refcat == UAC)
	    strcat (filename,".uac");
	else {
	    strcat (filename,".");
	    strcat (filename,refcatname);
	    }
	fd = fopen (filename, "w");
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

    /* Write heading */
    if (refcat == GSC)
	sprintf (headline, "CATALOG	HSTGSC1.1");
    else if (refcat == UAC)
	sprintf (headline, "CATALOG	USNO A 1.0");
    else if (refcat == UJC)
	sprintf (headline, "CATALOG	USNO UJ1.0");
    else
	sprintf (headline, "CATALOG	%s", refcatname);
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
	    fprintf (fd, "DEC	%s\n", rastr);
	if (tabout)
	    printf ("DEC	%s\n", rastr);
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
	if (strcmp (coorout,"FK4") == 0)
	    sprintf (headline,"GSC_NUMBER	RA1950    	DEC1950    	MAG   	TYPE	DISTANCE");
	else
	    sprintf (headline,"GSC_NUMBER	RA2000    	DEC2000    	MAG   	TYPE	DISTANCE");
    else if (refcat == UAC)
	if (strcmp (coorout,"FK4") == 0)
	    sprintf (headline,"UAC_NUMBER	RA1950  	DEC1950  	MAGB	MAGR	X    	Y    	Plate");
	else
	    sprintf (headline,"UAC_NUMBER	RA2000  	DEC2000  	MAGB	MAGR	X    	Y    	Plate");
    else if (refcat == UJC)
	if (strcmp (coorout,"FK4") == 0)
	    sprintf (headline,"UJC_NUMBER	RA1950  	DEC1950  	MAG   	Type	Distance");
	else
	    sprintf (headline,"UJC_NUMBER	RA2000  	DEC2000  	MAG   	Type	Distance");
    else
	if (strcmp (coorout,"FK4") == 0)
	    sprintf (headline,"NUMBER	RA1950  	DEC1950  	MAG   	Type	Distance");
	else
	    sprintf (headline,"NUMBER	RA2000  	DEC2000  	MAG   	Type	Distance");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);

    if (refcat == UAC)
	sprintf(headline,"----------	--------	---------	----	-----	-----	-----	-----");
    else
        sprintf (headline,"----------	------------	------------	------	----	_______");
    if (wfile)
	fprintf (fd, "%s\n", headline);
    if (tabout)
	printf ("%s\n", headline);
    if (printhead) {
	if (refcat == GSC)
	    if (nbg == 0)
		printf ("No Guide Stars Found\n");
	    else if (strcmp (coorout,"FK4") == 0)
		printf ("GSC number RA1950       Dec1950       Mag Type Arcsec\n");
	    else
		printf ("GSC number RA2000       Dec2000       Mag Type Arcsec\n");
	else if (refcat == UAC)
	    if (nbg == 0)
		printf ("No USNO A 1.0 Stars Found\n");
	    else if (strcmp (coorout,"FK4") == 0)
		printf ("USNO A number  RA1950       Dec1950      MagB  MagR Plate  Arcsec\n");
	    else
		printf ("USNO A number  RA2000       Dec2000      MagB  MagR Plate  Arcsec\n");
	else if (refcat == UJC)
	    if (nbg == 0)
		printf ("No UJ 1.0 Stars Found\n");
	    else if (strcmp (coorout,"FK4") == 0)
		printf (" UJ number    RA1950       Dec1950       Mag  Plate Arcsec\n");
	    else
		printf (" UJ number    RA2000       Dec2000       Mag  Plate Arcsec\n");
	else
	    if (nbg == 0)
		printf ("No Stars Found\n");
	    else if (strcmp (coorout,"FK4") == 0)
		printf (" Number    RA1950       Dec1950       Mag   Peak    Arcsec\n");
	    else
		printf (" Number    RA2000       Dec2000       Mag   Peak    Arcsec\n");
	}

    for (i = 0; i < nbg; i++) {
	if (gx[i] > 0.0 && gy[i] > 0.0) {
	    ra2str (rastr, gra[i], 3);
	    dec2str (decstr, gdec[i], 2);
	    if (refcat == UAC)
		sprintf (headline, "%13.8f	%s	%s	%.1f	%.1f	%.1f	%.1f	%.1f	%d",
		 gnum[i], rastr, decstr, gmb[i], gm[i], gx[i], gy[i], gc[i]);
	    else if (refcat == UJC)
	        sprintf (headline, "%12.7f	%s	%s	%.2f	%d	%.2f",
		 gnum[i], rastr, decstr, gm[i], gc[i], 3600.0*gx[i]);
	    else
	        sprintf (headline, "%9.4f	%s	%s	%.2f	%d	%.2f",
		 gnum[i], rastr, decstr, gm[i], gc[i], 3600.0*gx[i]);
	    if (wfile)
		fprintf (fd, "%s\n", headline);
	    else if (tabout)
		printf ("%s\n", headline);
	    else if (refcat == UAC)
		printf ("%13.8f %s %s %5.1f %5.1f %4d %8.2f\n",
			gnum[i],rastr,decstr,gmb[i],gm[i],gc[i], 3600.0*gx[i]);
	    else if (refcat == UJC)
		printf ("%12.7f %s %s %6.2f %4d %7.2f\n",
			gnum[i], rastr, decstr, gm[i],gc[i], 3600.0*gx[i]);
	    else if (refcat == GSC)
		printf ("%9.4f %s %s %6.2f %2d %7.2f\n",
			gnum[i], rastr, decstr, gm[i],gc[i], 3600.0*gx[i]);
	    else
		printf ("%9.4f %s %s %6.2f %7d %7.2f\n",
			gnum[i], rastr, decstr, gm[i],gc[i], 3600.0*gx[i]);
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

    return (nbg);
}

/* Oct 18 1996	New program based on imtab
 * Nov 13 1996	Set maximum nstar from command line if greater than default
 * Nov 14 1996	Set limits from subroutine
 * Nov 19 1996	Fix usage
 * Dec 12 1996	Allow bright as well as faint magnitude limit
 * Dec 12 1996	Fix header for UAC
 * Dec 12 1996	Fix header for UAC magnitudes
 * Dec 13 1996	Write plate into header if selected
 * Dec 18 1996	Allow WCS sky coordinate format as input argument for center
 * Dec 18 1996	Add option to print entries for specified catalog numbers
 * Dec 30 1996	Clean up closest star message
 * Dec 30 1996	Print message instead of heading if no stars are found
 * Jan 10 1997	Fix bug in RASort Stars which did not sort magnitudes
 */
