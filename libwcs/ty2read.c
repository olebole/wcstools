/*** File libwcs/ty2read.c
 *** December 11, 2000
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
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

/* pathname of Tycho 2 CDROM or catalog search engine URL */
char ty2cd[64]="/data/catalogs/tycho2";

static int ty2reg();
static int ty2regn();
static int ty2zone();
static int ty2size();
struct StarCat *ty2open();
void ty2close();
static int ty2star();
static int ty2size();

/* TY2READ -- Read Tycho 2 Star Catalog stars from CDROM */

int
ty2read (cra,cdec,dra,ddec,drad,distsort,sysout,eqout,epout,mag1,mag2,nstarmax,
	 gnum,gra,gdec,gpra,gpdec,gmag,gmagb,gtype,nlog)

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
double  *gpra;          /* Array of right ascension proper motions (returned) */
double  *gpdec;         /* Array of declination proper motions (returned) */
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
    double *gdist;	/* Array of distances to stars */
    int nreg;		/* Number of Tycho 2 regions in search */
    int rlist[MAXREG];	/* List of first stars in regions */
    int nlist[MAXREG];	/* List of number of stars per region */
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
    int nrmax = MAXREG;
    int nstar,i, ntot;
    int istar, istar1, istar2, isp;
    double num, ra, dec, rapm, decpm, mag, magb;
    double rra1, rra2, rra2a, rdec1, rdec2;
    char cstr[32];
    char *str;

    ntot = 0;
    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    /* If pathname is a URL, search and return */
    if ((str = getenv("TY2_PATH")) != NULL ) {
	if (!strncmp (str, "http:",5)) {
	    return (webread (str,"tycho2",distsort,cra,cdec,dra,ddec,drad,
			     sysout,eqout,epout,mag1,mag2,nstarmax,
			     gnum,gra,gdec,gpra,gpdec,gmag,gmagb,gtype,nlog));
	    }
	}
    if (!strncmp (ty2cd, "http:",5)) {
	return (webread (ty2cd,"tycho2",distsort,cra,cdec,dra,ddec,drad,
			 sysout,eqout,epout,mag1,mag2,nstarmax,
			 gnum,gra,gdec,gpra,gpdec,gmag,gmagb,gtype,nlog));
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
    gdist = (double *) malloc (nstarmax * sizeof (double));

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

    /* If searching through RA = 0:00, split search in two */
    for (iw = 0; iw <= wrap; iw++) {

	/* Find Tycho 2 Star Catalog regions in which to search */
	nreg = ty2reg (rra1,rra2,rdec1,rdec2,nrmax,rlist,nlist,verbose);
	if (nreg <= 0) {
	    fprintf (stderr,"TY2READ:  no Tycho 2 regions found\n");
	    return (0);
	    }

	/* Loop through region list */
	for (ireg = 0; ireg < nreg; ireg++) {

	    /* Open catalog file for this region */
	    istar1 = rlist[ireg];
	    istar2 = istar1 + nlist[ireg];
	    if (verbose)
		printf ("TY2READ: Searching stars %d through %d\n",
			istar1, istar2-1);

	    /* Open file for this region of Tycho 2 catalog */
	    starcat = ty2open (rlist[ireg], nlist[ireg]);
	    if (starcat == NULL) {
		fprintf (stderr,"TY2READ: File %s not found\n",inpath);
		return (0);
		}

	    /* Loop through catalog for this region */
	    for (istar = istar1; istar < istar2; istar++) {
		if (ty2star (starcat, star, istar)) {
		    fprintf (stderr,"TY2READ: Cannot read star %d\n", istar);
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

		/* Spectral Type */
		isp = (1000 * (int) star->isp[0]) + (int)star->isp[1];

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

		    /* Save star position and magnitude in table */
		    if (nstar < nstarmax) {
			gnum[nstar] = num;
			gra[nstar] = ra;
			gdec[nstar] = dec;
			gpra[nstar] = rapm;
			gpdec[nstar] = decpm;
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
			    gpra[farstar] = rapm;
			    gpdec[farstar] = decpm;
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
			gpra[faintstar] = rapm;
			gpdec[faintstar] = decpm;
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
			fprintf (stderr,"TY2READ: %11.6f: %9.5f %9.5f %5.2f %5.2f\n",
				 num,ra,dec,magb,mag);

		    /* End of accepted star processing */
		    }

		/* Log operation */
		jstar++;
		if (nlog > 0 && istar%nlog == 0)
		    fprintf (stderr,"TY2READ: %5d / %5d / %5d sources\r",
			     nstar,jstar,starcat->nstars);

		/* End of star loop */
		}

	    ntot = ntot + starcat->nstars;
	    if (nlog > 0)
		fprintf (stderr,"TY2READ: %4d / %4d: %5d / %5d  / %5d sources from region %4d    \n",
		 	 ireg+1,nreg,nstar,jstar,starcat->nstars,rlist[ireg]);

	    /* Close region input file */
	    ty2close (starcat);
	    }
	rra1 = 0.0;
	rra2 = rra2a;
	}

/* close output file and summarize transfer */
    if (nlog > 0) {
	if (nreg > 1)
	    fprintf (stderr,"TY2READ: %d regions: %d / %d found\n",nreg,nstar,ntot);
	else
	    fprintf (stderr,"TY2READ: 1 region: %d / %d found\n",nstar,ntot);
	if (nstar > nstarmax)
	    fprintf (stderr,"TY2READ: %d stars found; only %d returned\n",
		     nstar,nstarmax);
	}
    free ((char *)gdist);
    return (nstar);
}

/* TY2RNUM -- Read HST Guide Star Catalog stars from CDROM */

int
ty2rnum (nstars,sysout,eqout,epout,
	 gnum,gra,gdec,gpra,gpdec,gmag,gmagb,gtype,nlog)

int	nstars;		/* Number of stars to find */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
double	*gnum;		/* Array of Guide Star numbers (returned) */
double	*gra;		/* Array of right ascensions (returned) */
double	*gdec;		/* Array of declinations (returned) */
double  *gpra;          /* Array of right ascension proper motions (returned) */
double  *gpdec;         /* Array of declination proper motions (returned) */
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
    double num, ra, dec, rapm, decpm, mag, magb;

    if (nlog == 1)
	verbose = 1;
    else
	verbose = 0;

    /* If pathname is a URL, search and return */
    if ((str = getenv("TY2_PATH")) != NULL ) {
	if (!strncmp (str, "http:",5)) {
	    return (webrnum (str,"tycho2",nstars,sysout,eqout,epout,
			     gnum,gra,gdec,gpra,gpdec,gmag,gmagb,gtype,nlog));
	    }
	}
    if (!strncmp (ty2cd, "http:",5)) {
	return (webrnum (ty2cd,"tycho2",nstars,sysout,eqout,epout,
			 gnum,gra,gdec,gpra,gpdec,gmag,gmagb,gtype,nlog));
	}

    /* Allocate catalog entry buffer */
    star = (struct Star *) calloc (1, sizeof (struct Star));
    star->num = 0.0;
    nstar = 0;

/* Loop through star list */
    for (jstar = 0; jstar < nstars; jstar++) {
	rnum = (int) (gnum[jstar] + 0.0000001);
	if (gnum[jstar]-(double)rnum > 0.0000001) {
	    ty2regn (rnum, &istar1, &istar2, verbose);
	    nstar = istar2 - istar1 + 1;
	    starcat = ty2open (istar1, nstar);
    	    if (starcat == NULL) {
		fprintf (stderr,"TY2RNUM: File %s not found\n",inpath);
		return (0);
		}
	    for (istar = istar1; istar < istar2; istar++) {
		if (ty2star (starcat, star, istar)) {
		    fprintf (stderr,"TY2RNUM: Cannot read star %d\n", istar);
		    gra[jstar] = 0.0;
		    gdec[jstar] = 0.0;
		    gmag[jstar] = 0.0;
		    gmagb[jstar] = 0.0;
		    gtype[jstar] = 0;
		    }
		else {
		    if (fabs (gnum[jstar] - star->num) < 0.0000005)
			break;
		    }
		}
	    ty2close (starcat);
	    }
	/* Find star directly in catalog */
	else {
	    istar = (int) (gnum[jstar] + 0.01);
	    starcat = ty2open (istar, 10);
    	    if (starcat == NULL) {
		fprintf (stderr,"TY2RNUM: File %s not found\n",inpath);
		return (0);
		}
	    if (ty2star (starcat, star, istar)) {
		fprintf (stderr,"TY2RNUM: Cannot read star %d\n", istar);
		gra[jstar] = 0.0;
		gdec[jstar] = 0.0;
		gmag[jstar] = 0.0;
		gmagb[jstar] = 0.0;
		gtype[jstar] = 0;
		continue;
		}
	    ty2close (starcat);
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

	/* Spectral Type */
	isp = (1000 * (int) star->isp[0]) + (int)star->isp[1];

	/* Save star position and magnitude in table */
	gnum[jstar] = num;
	gra[jstar] = ra;
	gdec[jstar] = dec;
	gpra[jstar] = rapm;
	gpdec[jstar] = decpm;
	gmag[jstar] = mag;
	gmagb[jstar] = magb;
	gtype[jstar] = isp;
	if (nlog == 1)
	    fprintf (stderr,"TY2RNUM: %11.6f: %9.5f %9.5f %5.2f %5.2f %s  \n",
		     num, ra, dec, magb, mag, star->isp);

	/* End of star loop */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"TY2RNUM: %d / %d found\n",nstar,starcat->nstars);

    return (nstars);
}


/* Tycho 2 region index for ty2regn() and ty2reg() */

/* First region in each declination zone */
int treg1[25]={1,594,1178,1729,2259,2781,3246,3652,4014,4294,4492,4615,4663,
               5260,5838,6412,6989,7523,8022,8464,8840,9134,9346,9490,9538};
 
/* Last region in each declination zone */
int treg2[24]={593,1177,1728,2258,2780,3245,3651,4013,4293,4491,4614,4662,
               5259,5837,6411,6988,7522,8021,8463,8839,9133,9345,9489,9537};

/* TY2REGN -- read the range of stars in a region from the Tycho 2 Catalog
 * index table.
 */

static int
ty2regn (region, star1, star2, verbose)

int	region;		/* Region to find */
int	*star1;		/* First star number in region (returned)*/
int	*star2;		/* Last star number in region (returned)*/
int	verbose;	/* 1 for diagnostics */

{
    char *tabpath;	/* Pathname for regions table */
    char *buffer;	/* Buffer to hold index table */
    char *line;
    char *str;
    int deczone;
    int nchar=44;	/* Number of characters per line in table */
    int lpath;

    *star1 = 0;
    *star2 = 0;

/* Find declination zone in which this region exists */
    for (deczone = 0; deczone < 24; deczone++) {
	if (region >= treg1[deczone] && region < treg1[deczone+1])
	    break;
	}
    if (deczone > 24)
	return (0);

/* Set path to Tycho 2 Catalog CDROM */
    if ((str = getenv("TY2_PATH")) != NULL ) {
	lpath = strlen (str) + 16;
	tabpath = (char *) malloc (lpath);
	strcpy (tabpath, str);
	}
    else {
	lpath = strlen (ty2cd) + 16;
	tabpath = (char *) malloc (lpath);
	strcpy (tabpath, ty2cd);
	}

/* Set pathname for index table file */
    strcat (tabpath,"/data/index.dat");

/* Read the index table */
    if ((buffer = getfilebuff (tabpath)) == NULL) {
	fprintf (stderr,"TY2REG:  error reading region table %s\n",tabpath);
	return (0);
	}

/* Read first star from regionth line of region table */
    line = buffer + ((region - 1) * nchar);
    *star1 = atoi (line);

/* Read last star + 1 from region+1th line of region table */
    *star2 = atoi (line+nchar);
    free (buffer);
    free (tabpath);
    return (1);
}


/* TY2REG -- search the Tycho 2 Catalog index table for fields
 * in the specified range of coordinates and magnitudes.
 * Build lists containing the first star and number of stars for each range.
 */

static int
ty2reg (ra1, ra2, dec1, dec2, nrmax, rstar, nstar, verbose)

double	ra1, ra2;	/* Right ascension limits in degrees */
double	dec1, dec2; 	/* Declination limits in degrees */
int	nrmax;		/* Maximum number of regions to find */
int	*rstar;		/* Region first star numbers (returned)*/
int	*nstar;		/* Region numbers of stars (returned)*/
int	verbose;	/* 1 for diagnostics */

{
    int nrgn;		/* Number of regions found (returned) */
    char *tabpath;	/* Pathname for regions table */
    char *buffer;	/* Buffer to hold index table */
    char *line;
    char *str;
    int nchar=44;	/* Number of characters per line in table */
    int nwrap;		/* 1 if 0h included in RA span*/
    int iwrap;
    int num1, num2;
    int irow,iz1,iz2,ir1,ir2,jr1,jr2,i;
    int nsrch,nsrch1;
    double ralow, rahi;
    double declow, dechi, decmin, decmax;

    for (i = 0; i < nrmax; i++) {
	rstar[i] = 0;
	nstar[i] = 0;
	}
    nrgn = 0;

/* Set path to Tycho 2 Catalog CDROM */
    if ((str = getenv("TY2_PATH")) != NULL ) {
	tabpath = (char *) malloc (strlen (str) + 16);
	strcpy (tabpath, str);
	}
    else {
	tabpath = (char *) malloc (strlen (ty2cd) + 16);
	strcpy (tabpath, ty2cd);
	}

/* Set pathname for index table file */
    strcat (tabpath,"/data/index.dat");

/* Read the index table */
    if ((buffer = getfilebuff (tabpath)) == NULL) {
	fprintf (stderr,"TY2REG:  error reading region table %s\n",tabpath);
	return (0);
	}

/* Find region range to search based on declination */
    iz1 = ty2zone (dec1);
    iz2 = ty2zone (dec2);
    jr1 = 0;
    jr2 = 0;
    nwrap = 1;

/* Search region northern hemisphere or only one region */
    if (iz2 >= iz1) {
	ir1 = treg1[iz1];
	ir2 = treg2[iz2];
	}

/* Search region in southern hemisphere with multiple regions */
    else if (dec1 < 0 && dec2 < 0) {
	ir1 = treg1[iz2];
	ir2 = treg2[iz1];
	}

/* Search region spans equator */
    else if (dec1 < 0 && dec2 >= 0) {
	ir1 = treg1[12];
	ir2 = treg2[iz1];
	jr1 = treg1[0];
	jr2 = treg2[iz2];
	nwrap = 2;
	}

    nsrch = ir2 - ir1 + 1;
    if (verbose)
	fprintf (stderr,"TY2REG: searching %d regions: %d - %d\n",nsrch,ir1,ir2);
    if (jr1 > 0) {
	nsrch1 = jr2 - jr1 + 1;
	if (verbose)
	    fprintf (stderr,"TY2REG: searching %d regions: %d - %d\n",nsrch1,jr1,jr2);
	}
    if (verbose)
	fprintf(stderr,"TY2REG: RA: %.5f - %.5f, Dec: %.5f - %.5f\n",ra1,ra2,dec1,dec2);

    nrgn = 0;

    for (iwrap = 0; iwrap < nwrap; iwrap++) {

	for (irow = ir1 - 1; irow < ir2; irow++) {

	/* Read next line of region table */
	    line = buffer + ((irow - 1) * nchar);

	/* Declination range of the gs region */
	/* note:  southern dechi and declow are reversed */
	    num1 = atoi (line);
	    num2 = atoi (line+nchar);
	    dechi = atof (line + 29);
	    declow = atof (line + 36);
	    if (dechi > declow) {
		decmin = declow;
		decmax = dechi;
		}
	    else {
		decmax = declow;
		decmin = dechi;
		}

	    if (decmax >= dec1 && decmin <= dec2) {

	    /* Right ascension range of the Guide Star Catalog region */
		ralow = atof (line + 15);
		rahi = atof (line + 22);
		if (rahi <= 0.0) rahi = 360.0;

	    /* Check RA if 0h RA not between region RA limits */
		if (ra1 < ra2) {

		    /* Add this region to list, if there is space */
		    if (ralow <= ra2 && rahi >= ra1) {
			if (verbose)
			    fprintf (stderr,"TY2REG: Region %d added to search\n",irow);
			if (nrgn < nrmax) {
			    rstar[nrgn] = num1;
			    nstar[nrgn] = num2 - num1;
			    nrgn = nrgn + 1;
			    }
			}
		    }

	    /* Check RA if 0h RA is between region RA limits */
		else {
		    if (ralow > rahi) rahi = rahi + 360.0;
		    if (ralow <= ra2 || rahi >= ra1) {

		    /* Add this region to list, if there is space */
			if (verbose)
			    fprintf (stderr,"GSCREG: Region %d added to search\n", irow);

			if (nrgn < nrmax) {
			    rstar[nrgn] = num1;
			    nstar[nrgn] = num2 - num1;
			    nrgn = nrgn + 1;
			    }
			}
		    }
		}
	    }

/* Handle wrap-around through the equator */
	ir1 = jr1;
	ir2 = jr2;
	jr1 = 0;
	jr2 = 0;
	}

    free (buffer);
    return (nrgn);
}

 
 
/*  TY2ZONE -- find the zone number where a declination can be found */
 
static int
ty2zone (dec)
 
double dec;		/* declination in degrees */
 
{
int zone;		/* gsc zone (returned) */
double  zonesize;
int ndeczones = 12;	/* number of declination zones per hemisphere */
 
/* width of declination zones */
    zonesize = 90.0 / ndeczones;
 
    zone = ((int) (dec / zonesize));
    if (dec < 0)
	zone = ndeczones - zone;
 
    return (zone);
}



/* TY2OPEN -- Open Tycho 2 catalog file, returning number of entries */

struct StarCat *
ty2open (nstar, nread)

int	nstar;	/* Number of first star to read */
int	nread;	/* Number of star entries to read */

{
    FILE *fcat;
    struct StarCat *sc;
    int lfile, lpath;
    int lread, lskip, nr;
    char *str;
    char *ty2file;
    char *ty2path;	/* Full pathname for catalog file */

    /* Set path to Tycho 2 Catalog CDROM */
    if ((str = getenv("TY2_PATH")) != NULL ) {
	lpath = strlen(str) + 18;
	ty2path = (char *) malloc (lpath);
	strcpy (ty2path, str);
	}
    else {
	lpath = strlen(ty2cd) + 18;
	ty2path = (char *) malloc (lpath);
	strcpy (ty2path, ty2cd);
	}

    /* Set pathname for catalog file */
    strcat (ty2path, "/data/catalog.dat");

    /* Find length of Tycho 2 catalog file */
    lfile = ty2size (ty2path);

    /* Check for existence of catalog */
    if (lfile < 2) {
	fprintf (stderr,"TY2OPEN: Binary catalog %s has no entries\n",ty2path);
	free (ty2path);
	return (NULL);
	}

    /* Open Tycho 2 file */
    if (!(fcat = fopen (ty2path, "r"))) {
	fprintf (stderr,"TY2OPEN: Tycho 2 file %s cannot be read\n",ty2path);
	free (ty2path);
	return (0);
	}

    /* Set Tycho 2 catalog header information */
    sc = (struct StarCat *) calloc (1, sizeof (struct StarCat));
    sc->byteswapped = 0;

    sc->nbent = 208;
    sc->nstars = lfile / sc->nbent;

    /* Separate filename from pathname and save in structure */
    ty2file = strrchr (ty2path,'/');
    if (ty2file)
	ty2file = ty2file + 1;
    else
	ty2file = ty2path;
    if (strlen (ty2file) < 24)
	strcpy (sc->isfil, ty2file);
    else
	strncpy (sc->isfil, ty2file, 23);

    /* Set other catalog information in structure */
    sc->inform = 'J';
    sc->coorsys = WCS_J2000;
    sc->epoch = 2000.0;
    sc->equinox = 2000.0;
    sc->ifcat = fcat;
    sc->sptype = 2;

    /* Tycho 2 stars are not RA-sorted within regions */
    sc->rasorted = 0;

    /* Read part of catalog into a buffer */
    lread = nread * sc->nbent;
    lskip = (nstar - 1) * sc->nbent;
    sc->catdata = NULL;
    if ((sc->catdata = calloc (1, lread+1)) != NULL) {
	fseek (fcat, lskip, 0);
	nr = fread (sc->catdata, 1, lread, fcat);
	if (nr < lread) {
	    fprintf (stderr,"TY2OPEN: Read %d / %d bytes\n", nr, lread);
            ty2close (sc);
	    free (ty2path);
            return (NULL);
            }
	sc->catlast = sc->catdata + lread;
	}
    else {
	fprintf (stderr,"TY2OPEN: Cannot allocate %d-byte buffer.\n", lread);
        ty2close (sc);
	free (ty2path);
	return (NULL);
	}
    sc->istar = nstar;
    free (ty2path);
    return (sc);
}


void
ty2close (sc)
struct StarCat *sc;	/* Star catalog descriptor */
{
    fclose (sc->ifcat);
    if (sc->catdata != NULL)
	free (sc->catdata);
    free (sc);
    return;
}


/* TY2STAR -- Get Tycho 2 catalog entry for one star;
              return 0 if successful */

static int
ty2star (sc, st, istar)

struct StarCat *sc;	/* Star catalog descriptor */
struct Star *st;	/* Current star entry */
int istar;	/* Star sequence number in Tycho 2 catalog region file */
{
    char *line;
    double regnum, starnum, multnum;

    /* Drop out if catalog pointer is not set */
    if (sc == NULL)
	return (1);

    /* Drop out if catalog is not open */
    if (sc->ifcat == NULL)
	return (2);

    /* Drop out if star number is too large */
    if (istar > sc->nstars) {
	fprintf (stderr, "TY2STAR:  %d  > %d is not in catalog\n",
		 istar, sc->nstars);
	return (3);
	}

    /* Move buffer pointer to start of correct star entry */
    if (istar > 0) {
	line = sc->catdata + ((istar - sc->istar) * sc->nbent);
	if (line >= sc->catlast) {
	    fprintf (stderr, "TY2STAR:  star %d past buffer\n", istar);
	    return (4);
	    }
	}

    /* Read catalog entry */
    if (sc->nbent > sc->catlast-line) {
	fprintf (stderr, "TY2STAR:  %d / %d bytes read\n",
		 sc->catlast - line, sc->nbent);
	return (5);
	}

    regnum = atof (line);
    starnum = atof (line+5);
    multnum = atof (line+11);
    st->num = regnum + (0.0001 * starnum) + (0.00001 * multnum);

    /* Read position in degrees */
    st->ra = atof (line+15);
    st->dec = atof (line+28);

    /* Read proper motion and convert it to to degrees/year */
    st->rapm = (atof (line+41) / 3600000.0) / cosdeg (st->dec);
    st->decpm = atof (line+49) / 3600000.0;

    /* Set B magnitude */
    st->xmag[1] = atof (line+110);

    /* Set V magnitude */
    st->xmag[0] = atof (line+123);

    /* Set main sequence spectral type */
    bv2sp (NULL, st->xmag[1], st->xmag[0], st->isp);

    return (0);
}

/* TY2SIZE -- return size of Tycho 2 catalog file in bytes */

static int
ty2size (filename)

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


/* Jun  2 2000	New program, based on actread.c and gscread.c
 * Jun 13 2000	Correctly order magnitudes: 0=V, 1=B
 * Jun 26 2000	Add coordinate system to SearchLim() arguments
 * Sep 25 2000	Set sc->sptype to 2 to indicate presence of spectral type
 * Nov 29 2000	Add option to read catalog using HTTP
 * Dec 11 2000	Accept catalog search engine URL in ty2cd[]
 */
