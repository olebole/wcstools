/*** File libwcs/tmcread.c
 *** June 27, 2001
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "fitsfile.h"
#include "wcs.h"
#include "wcscat.h"

#define MAXREG 100
#define ABS(a) ((a) < 0 ? (-(a)) : (a))

/* pathname of 2MASS point source catalog root directory
   or catalog search engine URL */
char tmccd[64]="/data/mc4/2MASS";
static double *gdist;	/* Array of distances to stars */
static int ndist = 0;

static int tmcreg();
static int tmcregn();
static int tmczone();
static int tmcsize();
struct StarCat *tmcopen();
void tmcclose();
static int tmcstar();
static int tmcsize();

/* TMCREAD -- Read 2MASS point source catalog stars from CDROM */

int
tmcread (cra,cdec,dra,ddec,drad,distsort,sysout,eqout,epout,mag1,mag2,nstarmax,
	 gnum,gra,gdec,gmag,gmagb,gtype,nlog)

double	cra;		/* Search center J2000 right ascension in degrees */
double	cdec;		/* Search center J2000 declination in degrees */
double	dra;		/* Search half width in right ascension in degrees */
double	ddec;		/* Search half-width in declination in degrees */
double	drad;		/* Limiting separation in degrees (ignore if 0) */
int	distsort;	/* 1 to sort stars by distance from center */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
double	mag1,mag2;	/* Limiting magnitudes (none if equal) */
int	nstarmax;	/* Maximum number of stars to be returned */
double	*gnum;		/* Array of Guide Star numbers (returned) */
double	*gra;		/* Array of right ascensions (returned) */
double	*gdec;		/* Array of declinations (returned) */
double	*gmag;		/* Array of visual magnitudes (returned) */
double	*gmagb;		/* Array of blue magnitudes (returned) */
int	*gtype;		/* Array of object types (returned) */
int	nlog;		/* 1 for diagnostics */
{
    double ra1,ra2;	/* Limiting right ascensions of region in degrees */
    double dec1,dec2;	/* Limiting declinations of region in degrees */
    double dist = 0.0;  /* Distance from search center in degrees */
    double faintmag=0.0; /* Faintest magnitude */
    double maxdist=0.0; /* Largest distance */
    int	faintstar=0;	/* Faintest star */
    int	farstar=0;	/* Most distant star */
    int nreg;		/* Number of 2MASS point source regions in search */
    int rlist[MAXREG];	/* List of regions */
    char inpath[128];	/* Pathname for input region file */
    int sysref=WCS_J2000;	/* Catalog coordinate system */
    double eqref=2000.0;	/* Catalog equinox */
    double epref=2000.0;	/* Catalog epoch */
    struct StarCat *starcat;
    struct Star *star;
    int verbose;
    int wrap;
    int ireg;
    int jstar, iw;
    int zone;
    int nrmax = MAXREG;
    int nstar,i, ntot;
    int istar, istar1, istar2, isp;
    double num, ra, dec, rapm, decpm, mag, magb;
    double rra1, rra2, rra2a, rdec1, rdec2;
    char cstr[32], rastr[32], decstr[32], numstr[32];
    char *str;

    ntot = 0;
    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    /* If pathname is a URL, search and return */
    if ((str = getenv("TMC_PATH")) != NULL ) {
	if (!strncmp (str, "http:",5)) {
	    return (webread (str,"tmcsc",distsort,cra,cdec,dra,ddec,drad,
			     sysout,eqout,epout,mag1,mag2,nstarmax,
			     gnum,gra,gdec,NULL,NULL,gmag,gmagb,gtype,nlog));
	    }
	}
    if (!strncmp (tmccd, "http:",5)) {
	return (webread (tmccd,"tmcsc",distsort,cra,cdec,dra,ddec,drad,
			 sysout,eqout,epout,mag1,mag2,nstarmax,
			 gnum,gra,gdec,NULL,NULL,gmag,gmagb,gtype,nlog));
	}

    wcscstr (cstr, sysout, eqout, epout);

    SearchLim (cra,cdec,dra,ddec,sysout,&ra1,&ra2,&dec1,&dec2,verbose);

/* If RA range includes zero, split it in two */
    wrap = 0;
    if (ra1 > ra2)
	wrap = 1;
    else
	wrap = 0;

/* make mag1 always the smallest magnitude */
    if (mag2 < mag1) {
	mag = mag2;
	mag2 = mag1;
	mag1 = mag;
	}

    /* Allocate table for distances of stars from search center */
    if (nstarmax > ndist) {
	if (ndist > 0)
	    free ((void *)gdist);
	gdist = (double *) malloc (nstarmax * sizeof (double));
	if (gdist == NULL) {
	    fprintf (stderr,"TMCREAD:  cannot allocate separation array\n");
	    return (0);
	    }
	ndist = nstarmax;
	}

    /* Allocate catalog entry buffer */
    star = (struct Star *) calloc (1, sizeof (struct Star));
    star->num = 0.0;

    nstar = 0;
    jstar = 0;

    rra1 = ra1;
    rra2 = ra2;
    rdec1 = dec1;
    rdec2 = dec2;
    RefLim (cra, cdec, dra, ddec, sysout, sysref, eqout, eqref, epout,
	    &rra1, &rra2, &rdec1, &rdec2, verbose);
    if (rra1 > rra2) {
	rra2a = rra2;
	rra2 = 360.0;
	if (!wrap) wrap = 1;
	}

    /* Write header if printing star entries as found */
    if (nstarmax < 1) {
	printf ("catalog	2MASS Point Source Catalog\n");
	ra2str (rastr, 31, cra, 3);
	printf ("ra	%s\n", rastr);
	dec2str (decstr, 31, cdec, 2);
	printf ("dec	%s\n", decstr);
	if (drad != 0.0)
	    printf ("radsec	%.1f\n", drad*3600.0);
	else {
	    printf ("drasec	%.1f\n", dra*3600.0* cos(degrad(cdec)));
	    printf ("ddecsec	%.1f\n", ddec*3600.0);
	    }
	printf ("radecsys	%s\n", cstr);
	printf ("equinox	%.3f\n", eqout);
	printf ("epoch	%.3f\n", epout);
	printf ("program	stmc 2.9.4, 26 June 2001, Doug Mink SAO\n");
	printf ("2mass_id  	ra          	dec         	");
	printf ("magj 	magh 	magk 	arcmin\n");
	printf ("----------	------------	------------	");
	printf ("-----	-----	-----	------\n");
	}

    /* If searching through RA = 0:00, split search in two */
    for (iw = 0; iw <= wrap; iw++) {

	/* Find 2MASS Point Source Catalog regions in which to search */
	nreg = tmcreg (rra1,rra2,rdec1,rdec2,nrmax,rlist,verbose);
	if (nreg <= 0) {
	    fprintf (stderr,"TMCREAD:  no 2MASS regions found\n");
	    return (0);
	    }

	/* Loop through region list */
	for (ireg = 0; ireg < nreg; ireg++) {

	    /* Open file for this region of 2MASS point source catalog */
	    zone = rlist[ireg];
	    starcat = tmcopen (zone);
	    if (starcat == NULL) {
		fprintf (stderr,"TMCREAD: File %s not found\n",inpath);
		return (0);
		}

	    /* Find first and last stars in this region */
	    istar1 = tmcsdec (starcat, star, zone, rdec1);
	    istar2 = tmcsdec (starcat, star, zone, rdec2);
	    if (verbose)
		fprintf (stderr,"TMCREAD: Searching stars %d through %d\n",
			istar1, istar2-1);

	    /* Loop through catalog for this region */
	    for (istar = istar1; istar < istar2; istar++) {
		if (tmcstar (starcat, star, zone, istar)) {
		    fprintf (stderr,"TMCREAD: Cannot read star %d\n", istar);
		    break;
		    }

		/* ID number */
		num = star->num;

		/* Get position in output coordinate system, equinox, epoch */
		rapm = star->rapm;
		decpm = star->decpm;
		ra = star->ra;
		dec = star->dec;
		wcsconp (sysref, sysout, eqref, eqout, epref, epout,
		 	 &ra, &dec, &rapm, &decpm);

		/* Magnitude */
		mag = star->xmag[0];
		magb = star->xmag[1];
		isp = (int) ((star->xmag[2] * 1000.0) + 0.5);

		/* Compute distance from search center */
		if (drad > 0 || distsort)
		    dist = wcsdist (cra,cdec,ra,dec);
		else
		    dist = 0.0;

		/* Check magnitude and position limits */
		if ((mag1 == mag2 || (mag >= mag1 && mag <= mag2)) &&
		    ((wrap && (ra <= ra1 || ra >= ra2)) ||
		    (!wrap && (ra >= ra1 && ra <= ra2))) &&
		    ((drad > 0.0 && dist <= drad) ||
     		    (drad == 0.0 && dec >= dec1 && dec <= dec2))) {

		    /* Write star position and magnitudes to stdout */
		    if (nstarmax < 1) {
			CatNum (TMPSC, -10, 0, num, numstr);
			ra2str (rastr, 31, ra, 3);
			dec2str (decstr, 31, dec, 2);
			dist = wcsdist (cra,cdec,ra,dec) * 60.0;
                        printf ("%s	%s	%s", numstr,rastr,decstr);
			printf ("	%.3f	%.3f	%.3f	%.2f\n",
				star->xmag[0],star->xmag[1],star->xmag[2],
				dist);
			}

		    /* Save star position and magnitudes in table */
		    else if (nstar < nstarmax) {
			gnum[nstar] = num;
			gra[nstar] = ra;
			gdec[nstar] = dec;
			gmag[nstar] = mag;
			gmagb[nstar] = magb;
			gtype[nstar] = isp;
			gdist[nstar] = dist;
			if (dist > maxdist) {
			    maxdist = dist;
			    farstar = nstar;
			    }
			if (mag > faintmag) {
			    faintmag = mag;
			    faintstar = nstar;
			    }
			}

		    /* If too many stars and distance sorting,
		       replace farthest star */
		    else if (distsort) {
			if (dist < maxdist) {
			    gnum[farstar] = num;
			    gra[farstar] = ra;
			    gdec[farstar] = dec;
			    gmag[farstar] = mag;
			    gmagb[farstar] = magb;
			    gtype[farstar] = isp;
			    gdist[farstar] = dist;

			    /* Find new farthest star */
			    maxdist = 0.0;
			    for (i = 0; i < nstarmax; i++) {
				if (gdist[i] > maxdist) {
				    maxdist = gdist[i];
				    farstar = i;
				    }
				}
			    }
			}

		    /* Else if too many stars, replace faintest star */
		    else if (mag < faintmag) {
			gnum[faintstar] = num;
			gra[faintstar] = ra;
			gdec[faintstar] = dec;
			gmag[faintstar] = mag;
			gmagb[faintstar] = magb;
			gtype[faintstar] = isp;
			gdist[faintstar] = dist;
			faintmag = 0.0;

			/* Find new faintest star */
			for (i = 0; i < nstarmax; i++) {
			    if (gmag[i] > faintmag) {
				faintmag = gmag[i];
				faintstar = i;
				}
			    }
			}

		    nstar++;
		    if (nlog == 1)
			fprintf (stderr,"TMCREAD: %11.6f: %9.5f %9.5f %5.2f %5.2f\n",
				 num,ra,dec,magb,mag);

		    /* End of accepted star processing */
		    }

		/* Log operation */
		jstar++;
		if (nlog > 0 && istar%nlog == 0)
		    fprintf (stderr,"TMCREAD: %5d / %5d / %5d sources\r",
			     nstar,jstar,starcat->nstars);

		/* End of star loop */
		}

	    ntot = ntot + starcat->nstars;
	    if (nlog > 0)
		fprintf (stderr,"TMCREAD: %4d / %4d: %5d / %5d  / %5d sources from region %4d    \n",
		 	 ireg+1,nreg,nstar,jstar,starcat->nstars,zone);

	    /* Close region input file */
	    tmcclose (starcat);
	    }
	rra1 = 0.0;
	rra2 = rra2a;
	}

/* close output file and summarize transfer */
    if (nlog > 0) {
	if (nreg > 1)
	    fprintf (stderr,"TMCREAD: %d regions: %d / %d found\n",nreg,nstar,ntot);
	else
	    fprintf (stderr,"TMCREAD: 1 region: %d / %d found\n",nstar,ntot);
	if (nstar > nstarmax)
	    fprintf (stderr,"TMCREAD: %d stars found; only %d returned\n",
		     nstar,nstarmax);
	}
    return (nstar);
}

/* TMCRNUM -- Read HST Guide Star Catalog stars from CDROM */

int
tmcrnum (nstars,sysout,eqout,epout,
	 gnum,gra,gdec,gmag,gmagb,gtype,nlog)

int	nstars;		/* Number of stars to find */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
double	*gnum;		/* Array of Guide Star numbers (returned) */
double	*gra;		/* Array of right ascensions (returned) */
double	*gdec;		/* Array of declinations (returned) */
double	*gmag;		/* Array of V magnitudes (returned) */
double	*gmagb;		/* Array of B magnitudes (returned) */
int	*gtype;		/* Array of object types (returned) */
int	nlog;		/* 1 for diagnostics */
{
    char inpath[128];	/* Pathname for input region file */
    int sysref=WCS_J2000;	/* Catalog coordinate system */
    double eqref=2000.0;	/* Catalog equinox */
    double epref=2000.0;	/* Catalog epoch */
    struct StarCat *starcat;
    struct Star *star;
    char *str;

    int verbose;
    int rnum;
    int jstar;
    int istar, istar1, istar2, nstar, isp;
    double num, ra, dec, rapm, decpm, mag, magb, dstar;

    if (nlog == 1)
	verbose = 1;
    else
	verbose = 0;

    /* If pathname is a URL, search and return */
    if ((str = getenv("TMC_PATH")) != NULL ) {
	if (!strncmp (str, "http:",5)) {
	    return (webrnum (str,"tycho2",nstars,sysout,eqout,epout,
			     gnum,gra,gdec,NULL,NULL,gmag,gmagb,gtype,nlog));
	    }
	}
    if (!strncmp (tmccd, "http:",5)) {
	return (webrnum (tmccd,"tycho2",nstars,sysout,eqout,epout,
			 gnum,gra,gdec,NULL,NULL,gmag,gmagb,gtype,nlog));
	}

    /* Allocate catalog entry buffer */
    star = (struct Star *) calloc (1, sizeof (struct Star));
    star->num = 0.0;
    nstar = 0;

/* Loop through star list */
    for (jstar = 0; jstar < nstars; jstar++) {
	rnum = (int) (gnum[jstar] + 0.0000000001);
	starcat = tmcopen (rnum);
    	if (starcat == NULL) {
	    fprintf (stderr,"TMCRNUM: File %s not found\n",inpath);
	    return (0);
	    }
	dstar = (gnum[jstar] - (double)rnum) * 10000000.0;
	istar = (int) (dstar + 0.5);
	if (tmcstar (starcat, star, rnum, istar)) {
	    fprintf (stderr,"TMCRNUM: Cannot read star %d\n", istar);
	    gra[jstar] = 0.0;
	    gdec[jstar] = 0.0;
	    gmag[jstar] = 0.0;
	    gmagb[jstar] = 0.0;
	    gtype[jstar] = 0;
	    continue;
	    }

	/* If star has been found in catalog */

	/* ID number */
	num = star->num;

	/* Position in degrees at designated epoch */
	ra = star->ra;
	dec = star->dec;
	rapm = star->rapm;
	decpm = star->decpm;
	wcsconp (sysref, sysout, eqref, eqout, epref, epout,
		     &ra, &dec, &rapm, &decpm);

	/* Magnitude */
	mag = star->xmag[0];
	magb = star->xmag[1];
	isp = (int) ((star->xmag[2] * 1000.0) + 0.5);

	/* Save star position and magnitude in table */
	gnum[jstar] = num;
	gra[jstar] = ra;
	gdec[jstar] = dec;
	gmag[jstar] = mag;
	gmagb[jstar] = magb;
	gtype[jstar] = isp;
	if (nlog == 1)
	    fprintf (stderr,"TMCRNUM: %11.6f: %9.5f %9.5f %5.2f %5.2f %s  \n",
		     num, ra, dec, magb, mag, star->isp);

	/* End of star loop */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"TMCRNUM: %d / %d found\n",nstar,starcat->nstars);

    tmcclose (starcat);
    return (nstars);
}

char rdir[50][4]={"0", "1", "2", "3", "4", "5a", "5b", "6a", "6b", "6c",
	"6d", "7a", "7b", "7c", "7d", "8a", "8b", "9", "10", "11", "12",
	"13", "14", "15", "16a", "16b", "17a", "17b", "17c", "17d", "17e",
	"17f", "17g", "17h", "18a", "18b", "18c", "18d", "19a", "19b",
	"19c", "19d", "20a", "20b", "20c", "20d", "21", "22", "23", ""};
double zmax[50]={00.000, 01.000, 02.000, 03.000, 04.000, 05.000, 05.500,
		 06.000, 06.250, 06.500, 06.750, 07.000, 07.250, 07.500,
		 07.750, 08.000, 08.500, 09.000, 10.000, 11.000, 12.000,
		 13.000, 14.000, 15.000, 16.000, 16.500, 17.000, 17.125,
		 17.250, 17.375, 17.500, 17.625, 17.750, 17.875, 18.000,
		 18.250, 18.500, 18.750, 19.000, 19.250, 19.500, 19.750,
		 20.000, 20.250, 20.500, 20.750, 21.000, 22.000, 23.000,
		 24.000};


/* TMCREG -- find the regions contained by the given RA/Dec limits
 * Build lists containing the first star and number of stars for each range.
 */

static int
tmcreg (ra1, ra2, dec1, dec2, nrmax, regions, verbose)

double	ra1, ra2;	/* Right ascension limits in degrees */
double	dec1, dec2; 	/* Declination limits in degrees */
int	nrmax;		/* Maximum number of regions to find */
int	*regions;	/* Region numbers to search (returned)*/
int	verbose;	/* 1 for diagnostics */

{
    int nrgn;		/* Number of regions found (returned) */
    char *tabpath;	/* Pathname for regions table */
    char *line;
    char *str;
    int nwrap;		/* 1 if 0h included in RA span*/
    int ir;
    int iz1,iz2,ir1,ir2,jr1,jr2,i;
    int nsrch;
    double rah1,rah2,ralow, rahi;
    double declow, dechi, decmin, decmax;

    nrgn = 0;

    /* Set path to 2MASS Catalog */
    if ((str = getenv("TMC_PATH")) != NULL ) {
	tabpath = (char *) malloc (strlen (str) + 16);
	strcpy (tabpath, str);
	}
    else {
	tabpath = (char *) malloc (strlen (tmccd) + 16);
	strcpy (tabpath, tmccd);
	}

    /* Find region range to search based on right ascension */
    rah1 = ra1 / 15.0;
    for (i = 1; i < 50; i++) {
	if (rah1 < zmax[i]) {
	    iz1 = i - 1;
	    break;
	    }
	}
    rah2 = ra2 / 15.0;
    for (i = 1; i < 50; i++) {
	if (rah2 < zmax[i]) {
	    iz2 = i - 1;
	    break;
	    }
	}
    if (iz2 >= iz1) {
	ir1 = iz1;
	ir2 = iz2;
	jr1 = 0;
	jr2 = 0;
	nwrap = 1;
	nsrch = iz2 - iz1 + 1;
	}
    else {
	ir1 = iz1;
	ir2 = 48;
	jr1 = 0;
	jr2 = iz2;
	nwrap = 2;
	nsrch = 48 - iz1 + 1 + iz2 + 1;
	}

    /* Search region northern hemisphere or only one region */
    if (verbose) {
	fprintf(stderr,"TMCREG: RA: %.5f - %.5f, Dec: %.5f - %.5f\n",
		ra1,ra2,dec1,dec2);
	if (nsrch == 1)
	    fprintf (stderr,"TMCREG: searching region %d", ir1);
	else
	    fprintf (stderr,"TMCREG: searching %d regions: %d - %d",
		 nsrch, ir1, ir2);
	if (jr1 > 0 && jr2 > 0)
	    fprintf (stderr,", %d - %d", jr1, jr2);
	fprintf (stderr,"\n");
	}

    /* Loop through first section of sky */
    nrgn = 0;
    for (ir = ir1; ir <= ir2; ir++) {
	if (verbose)
	    fprintf (stderr,"TMCREG: Region %d (%s) added to search\n",
		    ir, rdir[ir]);

	/* Add this region to list, if there is space */
	if (nrgn < nrmax) {
	    regions[nrgn] = ir;
	    nrgn++;
	    }
	}
    for (ir = jr1; ir < jr2; ir++) {
	if (verbose)
	    fprintf (stderr,"TMCREG: Region %d %s) added to search\n",
		     ir, rdir[ir]);

	/* Add this region to list, if there is space */
	if (nrgn < nrmax) {
	    regions[nrgn] = ir;
	    nrgn++;
	    }
	}

    return (nrgn);
}


/* TMCOPEN -- Open 2MASS point source catalog file, returning catalog structure */

struct StarCat *
tmcopen (zone)

int	zone;	/* RA zone (hours) to read */

{
    FILE *fcat;
    struct StarCat *sc;
    int lfile, lpath;
    int lread, lskip, nr;
    char *str;
    char *tmcfile;
    char *tmcpath;	/* Full pathname for catalog file */

    /* Set path to 2MASS Point Source Catalog zone */
    if ((str = getenv("TMC_PATH")) != NULL ) {
	lpath = strlen(str) + 18;
	tmcpath = (char *) malloc (lpath);
	sprintf (tmcpath, "%s/idr2psc%s.tbl", str, rdir[zone]);
	}
    else {
	lpath = strlen (tmccd) + 18;
	tmcpath = (char *) malloc (lpath);
	sprintf (tmcpath, "%s/idr2psc%s.tbl", tmccd, rdir[zone]);
	}

    /* Find length of 2MASS catalog file */
    lfile = tmcsize (tmcpath);

    /* Check for existence of catalog */
    if (lfile < 2) {
	fprintf (stderr,"TMCOPEN: Binary catalog %s has no entries\n",tmcpath);
	free (tmcpath);
	return (NULL);
	}

    /* Open 2MASS point source catalog zone file */
    if (!(fcat = fopen (tmcpath, "r"))) {
	fprintf (stderr,"TMCOPEN: 2MASS PSC file %s cannot be read\n",tmcpath);
	free (tmcpath);
	return (0);
	}

    /* Set 2MASS PSC catalog header information */
    sc = (struct StarCat *) calloc (1, sizeof (struct StarCat));
    sc->byteswapped = 0;

    sc->nbent = 302;
    sc->nstars = lfile / sc->nbent;

    /* Separate filename from pathname and save in structure */
    tmcfile = strrchr (tmcpath,'/');
    if (tmcfile)
	tmcfile = tmcfile + 1;
    else
	tmcfile = tmcpath;
    if (strlen (tmcfile) < 24)
	strcpy (sc->isfil, tmcfile);
    else
	strncpy (sc->isfil, tmcfile, 23);

    /* Set other catalog information in structure */
    sc->inform = 'J';
    sc->coorsys = WCS_J2000;
    sc->epoch = 2000.0;
    sc->equinox = 2000.0;
    sc->ifcat = fcat;
    sc->sptype = 2;

    /* 2MASS stars are Dec-sorted within regions */
    sc->rasorted = 0;

    return (sc);
}


void
tmcclose (sc)
struct StarCat *sc;	/* Star catalog descriptor */
{
    fclose (sc->ifcat);
    if (sc->catdata != NULL)
	free (sc->catdata);
    free (sc);
    return;
}


/* TMCSDEC -- Find 2MASS star closest to specified declination */

static int
tmcsdec (starcat, star, zone, decx0)

struct StarCat *starcat; /* Star catalog descriptor */
struct Star *star;	/* Current star entry */
int	zone;		/* RA zone in which search is occuring */
double	decx0;		/* Declination in degrees for which to search */
{
    int istar, istar1, istar2, nrep;
    double decx, dec1, dec, rdiff, rdiff1, rdiff2, sdiff;
    char decstrx[16];
    int debug = 0;

    decx = decx0;
    if (debug)
	dec2str (decstrx, 16, decx, 3);
    istar1 = 1;
    if (tmcstar (starcat, star, zone, istar1))
	return (0);
    dec1 = star->dec;
    istar = starcat->nstars;
    nrep = 0;
    while (istar != istar1 && nrep < 20) {
	if (tmcstar (starcat, star, zone, istar))
	    break;
	else {
	    dec = star->dec;
	    if (dec == dec1)
		break;
	    if (debug) {
		char decstr[16];
		dec2str (decstr, 16, dec, 3);
		fprintf (stderr,"UACSRA %d %d: %s (%s)\n",
			 nrep,istar,decstr,decstrx);
		}
	    rdiff = dec1 - dec;
	    rdiff1 = dec1 - decx;
	    rdiff2 = dec - decx;
	    if (nrep > 20 && ABS(rdiff2) > ABS(rdiff1)) {
		istar = istar1;
		break;
		}
	    nrep++;
	    sdiff = (double)(istar - istar1) * rdiff1 / rdiff;
	    istar2 = istar1 + (int) (sdiff + 0.5);
	    dec1 = dec;
	    istar1 = istar;
	    istar = istar2;
	    if (debug) {
		fprintf (stderr," dec1=    %.5f dec=     %.5f decx=    %.5f\n",
			 dec1,dec,decx);
		fprintf (stderr," rdiff=  %.5f rdiff1= %.5f rdiff2= %.5f\n",
			 rdiff,rdiff1,rdiff2);
		fprintf (stderr," istar1= %d istar= %d istar1= %d\n",
			 istar1,istar,istar2);
		}
	    if (istar < 1)
		istar = 1;
	    if (istar > starcat->nstars)
		istar = starcat->nstars;
	    if (istar == istar1)
		break;
	    }
	}
    return (istar);
}


/* TMCSTAR -- Get 2MASS Point Source Catalog entry for one star;
              return 0 if successful */

static int
tmcstar (sc, st, zone, istar)

struct StarCat *sc;	/* Star catalog descriptor */
struct Star *st;	/* Current star entry */
int	zone;		/* Zone catalog number (1-49) */
int	istar;		/* Star sequence in 2MASS zone file */
{
    char line[500];
    double regnum, starnum, multnum;
    int nbskip, nbr;

    /* Drop out if catalog pointer is not set */
    if (sc == NULL)
	return (1);

    /* Drop out if catalog is not open */
    if (sc->ifcat == NULL)
	return (2);

    /* Drop out if star number is too large */
    if (istar > sc->nstars) {
	fprintf (stderr, "TMCSTAR:  %d  > %d is not in catalog\n",
		 istar, sc->nstars);
	return (3);
	}

    /* Read entry for one star */
    nbskip = sc->nbent * (istar - 1);
    if (fseek (sc->ifcat,nbskip,SEEK_SET))
	return (-1);
    nbr = fread (line, sc->nbent, 1, sc->ifcat) * sc->nbent;
    if (nbr < sc->nbent) {
	fprintf (stderr, "tmcstar %d / %d bytes read\n",nbr, sc->nbent);
	return (-2);
	}

    /* Make up source number from zone number and star number */
    st->num = zone + (0.0000001 * (double) istar);

    /* Read position in degrees */
    st->ra = atof (line);
    st->dec = atof (line+11);

    /* No proper motion */
    st->rapm = 0.0;
    st->decpm = 0.0;

    /* Set J magnitude */
    st->xmag[0] = atof (line+53);

    /* Set H magnitude */
    st->xmag[1] = atof (line+72);

    /* Set K magnitude */
    st->xmag[2] = atof (line+91);

    return (0);
}

/* TMCSIZE -- return size of 2MASS point source catalog file in bytes */

static int
tmcsize (filename)

char	*filename;	/* Name of file for which to find size */
{
    FILE *diskfile;
    long filesize;

    /* Open file */
    if ((diskfile = fopen (filename, "r")) == NULL)
	return (-1);

    /* Move to end of the file */
    if (fseek (diskfile, 0, 2) == 0)

	/* Position is the size of the file */
	filesize = ftell (diskfile);

    else
	filesize = -1;

    fclose (diskfile);

    return (filesize);
}


/* May 29 2001	New program, based on ty2read.c and uacread.c
 * May 30 2001	Round K magnitude to nearest 1000th
 * Jun 13 2001	Round star number up to avoid truncation problem
 * Jun 27 2001	Add code to print one entry at time if nstars < 1
 * Jun 27 2001	Allocate gdist only if larger array is needed
 */