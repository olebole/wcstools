/*** File libwcs/ucacread.c
 *** April 14, 2003
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 *** Copyright (C) 2003
 *** Smithsonian Astrophysical Observatory, Cambridge, MA, USA

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Correspondence concerning WCSTools should be addressed as follows:
           Internet email: dmink@cfa.harvard.edu
           Postal address: Doug Mink
                           Smithsonian Astrophysical Observatory
                           60 Garden St.
                           Cambridge, MA 02138 USA
 */

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "fitsfile.h"
#include "wcs.h"
#include "wcscat.h"

#define MAXREG 100

/* pathname of UCAC1 decompressed data files or search engine URL */
char ucacd[64]="/data/astrocat/ucac1";

static double *gdist;	/* Array of distances to stars */
static int ndist = 0;

static int ucacreg();
static int ucacregn();
static int ucaczone();
static int ucacsize();
struct StarCat *ucacopen();
void ucacclose();
static int ucacstar();
static int ucacsize();

/* UCACREAD -- Read UCAC Star Catalog stars */

int
ucacread (cra,cdec,dra,ddec,drad,distsort,sysout,eqout,epout,mag1,mag2,sortmag,
	 nstarmax,gnum,gra,gdec,gpra,gpdec,gmag,gtype,nlog)

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
int	sortmag;	/* Magnitude by which to sort (1 or 2) */
int	nstarmax;	/* Maximum number of stars to be returned */
double	*gnum;		/* Array of Guide Star numbers (returned) */
double	*gra;		/* Array of right ascensions (returned) */
double	*gdec;		/* Array of declinations (returned) */
double  *gpra;          /* Array of right ascension proper motions (returned) */
double  *gpdec;         /* Array of declination proper motions (returned) */
double	**gmag;		/* Array of b and v magnitudes (returned) */
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
    int nreg;		/* Number of UCAC regions in search */
    int regnum[MAXREG];	/* List of region numbers */
    int rlist[MAXREG];	/* List of first stars in regions */
    int nlist[MAXREG];	/* List of number of stars per region */
    char inpath[128];	/* Pathname for input region file */
    int sysref = WCS_J2000;	/* Catalog coordinate system */
    double eqref = 2000.0;	/* Catalog equinox */
    double epref = 2000.0;	/* Catalog epoch */
    struct StarCat *starcat;
    struct Star *star;
    int verbose;
    int wrap;
    int ireg;
    int magsort, magsort1;
    int jstar, iw;
    int nrmax = MAXREG;
    int nstar,i, ntot;
    int istar, istar1, istar2;
/*    int isp; */
    int pass;
    double num, ra, dec, rapm, decpm, mag, magb, magv;
    double rra1, rra2, rra2a, rdec1, rdec2, dmarg;
    double rdist, ddist;
    char cstr[32], rastr[32], decstr[32];
    char *str;

    ntot = 0;
    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    /* If pathname is a URL, search and return */
    if ((str = getenv("UCAC_PATH")) != NULL ) {
	if (!strncmp (str, "http:",5)) {
	    return (webread (str,"ucac1",distsort,cra,cdec,dra,ddec,drad,
			     sysout,eqout,epout,mag1,mag2,sortmag,nstarmax,
			     gnum,gra,gdec,gpra,gpdec,gmag,gtype,nlog));
	    }
	}
    if (!strncmp (ucaccd, "http:",5)) {
	return (webread (ucaccd,"ucac1",distsort,cra,cdec,dra,ddec,drad,
			 sysout,eqout,epout,mag1,mag2,sortmag,nstarmax,
			 gnum,gra,gdec,gpra,gpdec,gmag,gtype,nlog));
	}

    wcscstr (cstr, sysout, eqout, epout);

    SearchLim (cra,cdec,dra,ddec,sysout,&ra1,&ra2,&dec1,&dec2,verbose);

    /* Make mag1 always the smallest magnitude */
    if (mag2 < mag1) {
	mag = mag2;
	mag2 = mag1;
	mag1 = mag;
	}

   if (sortmag == 2) {
	magsort = 0;
	magsort1 = 1;
	}
    else {
	magsort = 1;
	magsort1 = 0;
	}

    /* Allocate table for distances of stars from search center */
    if (nstarmax > ndist) {
	if (ndist > 0)
	    free ((void *)gdist);
	gdist = (double *) malloc (nstarmax * sizeof (double));
	if (gdist == NULL) {
	    fprintf (stderr,"UCACREAD:  cannot allocate separation array\n");
	    return (0);
	    }
	ndist = nstarmax;
	}

    /* Allocate catalog entry buffer */
    star = (struct Star *) calloc (1, sizeof (struct Star));
    star->num = 0.0;

    nstar = 0;
    jstar = 0;

    /* Get RA and Dec limits in catalog (J2000) coordinates */
    rra1 = ra1;
    rra2 = ra2;
    rdec1 = dec1;
    rdec2 = dec2;
    RefLim (cra, cdec, dra, ddec, sysout, sysref, eqout, eqref, epout,
	    &rra1, &rra2, &rdec1, &rdec2, verbose);

    /* Add 60 arcsec/century margins to region to get most stars which move */
    if (epout == 0.0)
	dmarg = 0.0;
    else
	dmarg = (60.0 / 3600.0) * fabs (epout - epref);
    rdec1 = rdec1 - dmarg;
    rdec2 = rdec2 + dmarg;
    rra1 = rra1 - (dmarg / cos (degrad(cdec)));
    rra2 = rra2 + (dmarg / cos (degrad(cdec)));

    /* Deal with search passing through 0:00 RA */
    if (rra1 > rra2) {
	rra2a = rra2;
	rra2 = 360.0;
	wrap = 1;
	}
    else
	wrap = 0;

    /* Write header if printing star entries as found */
    if (nstarmax < 1) {
	char *revmessage;
	revmessage = getrevmsg();
	printf ("catalog	UCAC1\n");
	ra2str (rastr, 31, cra, 3);
	printf ("ra	%s\n", rastr);
	dec2str (decstr, 31, cdec, 2);
	printf ("dec	%s\n", decstr);
	printf ("rpmunit	tsec/century\n");
	printf ("dpmunit	arcsec/century\n");
	if (drad != 0.0)
	    printf ("radmin	%.1f\n", drad*60.0);
	else {
	    printf ("dramin	%.1f\n", dra*60.0* cosdeg (cdec));
	    printf ("ddecmin	%.1f\n", ddec*60.0);
	    }
	printf ("radecsys	%s\n", cstr);
	printf ("equinox	%.3f\n", eqout);
	printf ("epoch	%.3f\n", epout);
	printf ("program	scat %s\n", revmessage);
	printf ("ucac_id	ra          	dec         	");
	printf ("magb 	magv 	ura   	udec  	arcmin\n");
	printf ("----------	------------	------------    ");
	printf ("-----	-----	------	------	------\n");
	}

    /* If searching through RA = 0:00, split search in two */
    for (iw = 0; iw <= wrap; iw++) {

	/* Find UCAC Star Catalog regions in which to search */
	nreg = ucacreg (rra1,rra2,rdec1,rdec2,nrmax,regnum,rlist,nlist,verbose);
	if (nreg <= 0) {
	    fprintf (stderr,"UCACREAD:  no UCAC region for %.2f-%.2f %.2f %.2f\n",
		     rra1, rra2, rdec1, rdec2);
	    rra1 = 0.0;
	    rra2 = rra2a;
	    continue;
	    }

	/* Loop through region list */
	for (ireg = 0; ireg < nreg; ireg++) {

	    /* Open catalog file for this region */
	    istar1 = rlist[ireg];
	    istar2 = istar1 + nlist[ireg];
	    if (verbose)
		fprintf (stderr,"UCACREAD: Searching stars %d through %d\n",
			istar1, istar2-1);

	    /* Open file for this region of UCAC catalog */
	    starcat = ucacopen (rlist[ireg], nlist[ireg]);
	    if (starcat == NULL) {
		fprintf (stderr,"UCACREAD: File %s not found\n",inpath);
		return (0);
		}

	    /* Loop through catalog for this region */
	    for (istar = istar1; istar < istar2; istar++) {
		if (ucacstar (starcat, star, istar)) {
		    fprintf (stderr,"UCACREAD: Cannot read star %d\n", istar);
		    break;
		    }

		/* ID number */
		num = star->num;

		/* Magnitude */
		magv = star->xmag[0];
		magb = star->xmag[1];
		mag = star->xmag[magsort];

		/* Check magnitude limits */
		pass = 1;
		if (mag1 != mag2 && (mag < mag1 || mag > mag2))
		    pass = 0;

		/* Check position limits */
		if (pass) {

		    /* Get position in output coordinate system */
		    rapm = star->rapm;
		    decpm = star->decpm;
		    ra = star->ra;
		    dec = star->dec;
		    wcsconp (sysref, sysout, eqref, eqout, epref, epout,
		 	     &ra, &dec, &rapm, &decpm);

		    /* Compute distance from search center */
		    if (drad > 0 || distsort)
			dist = wcsdist (cra,cdec,ra,dec);
		    else
			dist = 0.0;

		    /* Check radial distance to search center */
		    if (drad > 0) {
			if (dist > drad)
			    pass = 0;
			}

		    /* Check distance along RA and Dec axes */
		    else {
			ddist = wcsdist (cra,cdec,cra,dec);
			if (ddist > ddec)
			    pass = 0;
			rdist = wcsdist (cra,dec,ra,dec);
		        if (rdist > dra)
			   pass = 0;
			}
		    }

		if (pass) {

		/* Spectral Type
		isp = (1000 * (int) star->isp[0]) + (int)star->isp[1]; */

		/* Write star position and magnitudes to stdout */
		    if (nstarmax < 1) {
			ra2str (rastr, 31, ra, 3);
			dec2str (decstr, 31, dec, 2);
			dist = wcsdist (cra,cdec,ra,dec) * 60.0;
			printf ("%010.5f	%s	%s", num,rastr,decstr);
			printf ("	%5.2f	%5.2f	%6.3f	%6.2f	%.2f\n",
				magb, magv, rapm*240000.0, decpm*3600000.0,
				dist / 60.0);
			}

		    /* Save star position and magnitude in table */
		    if (nstar < nstarmax) {
			gnum[nstar] = num;
			gra[nstar] = ra;
			gdec[nstar] = dec;
			gpra[nstar] = rapm;
			gpdec[nstar] = decpm;
			gmag[0][nstar] = magb;
			gmag[1][nstar] = magv;
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
			    gmag[0][farstar] = magb;
			    gmag[1][farstar] = magv;
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
			gmag[0][faintstar] = magb;
			gmag[1][faintstar] = magv;
			gdist[faintstar] = dist;
			faintmag = 0.0;

			/* Find new faintest star */
			for (i = 0; i < nstarmax; i++) {
			    if (gmag[magsort1][i] > faintmag) {
				faintmag = gmag[magsort1][i];
				faintstar = i;
				}
			    }
			}

		    nstar++;
		    if (nlog == 1)
			fprintf (stderr,"UCACREAD: %11.6f: %9.5f %9.5f %5.2f %5.2f\n",
				 num,ra,dec,magb,magv);

		    /* End of accepted star processing */
		    }

		/* Log operation */
		jstar++;
		if (nlog > 0 && istar%nlog == 0)
		    fprintf (stderr,"UCACREAD: %5d / %5d / %5d sources\r",
			     nstar,jstar,starcat->nstars);

		/* End of star loop */
		}

	    ntot = ntot + starcat->nstars;
	    if (nlog > 0)
		fprintf (stderr,"UCACREAD: %4d / %4d: %5d / %5d  / %5d sources from region %4d    \n",
		 	 ireg+1,nreg,nstar,jstar,starcat->nstars,rlist[ireg]);

	    /* Close region input file */
	    ucacclose (starcat);
	    }
	rra1 = 0.0;
	rra2 = rra2a;
	}

/* close output file and summarize transfer */
    if (nlog > 0) {
	if (nreg > 1)
	    fprintf (stderr,"UCACREAD: %d regions: %d / %d found\n",nreg,nstar,ntot);
	else
	    fprintf (stderr,"UCACREAD: 1 region: %d / %d found\n",nstar,ntot);
	if (nstar > nstarmax)
	    fprintf (stderr,"UCACREAD: %d stars found; only %d returned\n",
		     nstar,nstarmax);
	}
    return (nstar);
}

/* UCACRNUM -- Read HST Guide Star Catalog stars from CDROM */

int
ucacrnum (nstars,sysout,eqout,epout,
	 gnum,gra,gdec,gpra,gpdec,gmag,gtype,nlog)

int	nstars;		/* Number of stars to find */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
double	*gnum;		/* Array of Guide Star numbers (returned) */
double	*gra;		/* Array of right ascensions (returned) */
double	*gdec;		/* Array of declinations (returned) */
double  *gpra;          /* Array of right ascension proper motions (returned) */
double  *gpdec;         /* Array of declination proper motions (returned) */
double	**gmag;		/* Array of B and V magnitudes (returned) */
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
    int istar, istar1, istar2, nstar;
/*    int isp; */
    double num, ra, dec, rapm, decpm, magb, magv;

    if (nlog == 1)
	verbose = 1;
    else
	verbose = 0;

    /* If pathname is a URL, search and return */
    if ((str = getenv("UCAC_PATH")) != NULL ) {
	if (!strncmp (str, "http:",5)) {
	    return (webrnum (str,"ucac1",nstars,sysout,eqout,epout,
			     gnum,gra,gdec,gpra,gpdec,gmag,gtype,nlog));
	    }
	}
    if (!strncmp (ucaccd, "http:",5)) {
	return (webrnum (ucaccd,"ucac1",nstars,sysout,eqout,epout,
			 gnum,gra,gdec,gpra,gpdec,gmag,gtype,nlog));
	}

    /* Allocate catalog entry buffer */
    star = (struct Star *) calloc (1, sizeof (struct Star));
    star->num = 0.0;
    nstar = 0;

/* Loop through star list */
    for (jstar = 0; jstar < nstars; jstar++) {
	rnum = (int) (gnum[jstar] + 0.0000001);
	
	/* Find numbered stars (rrr.nnnnnn) */
	if (gnum[jstar]-(double)rnum > 0.0000001) {
	    ucacregn (rnum, &istar1, &istar2, verbose);
	    nstar = istar2 - istar1 + 1;
	    starcat = ucacopen (istar1, nstar);
    	    if (starcat == NULL) {
		fprintf (stderr,"UCACRNUM: File %s not found\n",inpath);
		return (0);
		}
	    for (istar = istar1; istar < istar2; istar++) {
		if (ucacstar (starcat, star, istar)) {
		    fprintf (stderr,"UCACRNUM: Cannot read star %d\n", istar);
		    gra[jstar] = 0.0;
		    gdec[jstar] = 0.0;
		    gmag[0][jstar] = 0.0;
		    gmag[1][jstar] = 0.0;
		    gtype[jstar] = 0;
		    }
		else {
		    if (fabs (gnum[jstar] - star->num) < 0.0000005)
			break;
		    }
		}
	    ucacclose (starcat);
	    }
	/* Find nth sequential stars in catalog (not rrrr.nnnnn) */
	else {
	    istar = (int) (gnum[jstar] + 0.01);
	    starcat = ucacopen (istar, 10);
    	    if (starcat == NULL) {
		fprintf (stderr,"UCACRNUM: File %s not found\n",inpath);
		return (0);
		}
	    if (ucacstar (starcat, star, istar)) {
		fprintf (stderr,"UCACRNUM: Cannot read star %d\n", istar);
		gra[jstar] = 0.0;
		gdec[jstar] = 0.0;
		gmag[0][jstar] = 0.0;
		gmag[1][jstar] = 0.0;
		gtype[jstar] = 0;
		continue;
		}
	    ucacclose (starcat);
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
	magv = star->xmag[0];
	magb = star->xmag[1];

	/* Spectral Type
	isp = (1000 * (int) star->isp[0]) + (int)star->isp[1]; */

	/* Save star position and magnitude in table */
	gnum[jstar] = num;
	gra[jstar] = ra;
	gdec[jstar] = dec;
	gpra[jstar] = rapm;
	gpdec[jstar] = decpm;
	gmag[0][jstar] = magb;
	gmag[1][jstar] = magv;
	/* gtype[jstar] = isp; */
	if (nlog == 1)
	    fprintf (stderr,"UCACRNUM: %11.6f: %9.5f %9.5f %5.2f %5.2f %s  \n",
		     num, ra, dec, magb, magv, star->isp);

	/* End of star loop */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"UCACRNUM: %d / %d found\n",nstar,starcat->nstars);

    return (nstars);
}


/* UCAC region index for ucacregn() and ucacreg() */

/* First region in each declination zone */
int treg1[25]={1,594,1178,1729,2259,2781,3246,3652,4014,4294,4492,4615,4663,
               5260,5838,6412,6989,7523,8022,8464,8840,9134,9346,9490,9538};
 
/* Last region in each declination zone */
int treg2[24]={593,1177,1728,2258,2780,3245,3651,4013,4293,4491,4614,4662,
               5259,5837,6411,6988,7522,8021,8463,8839,9133,9345,9489,9537};

/* UCACREGN -- read the range of stars in a region from the UCAC Catalog
 * index table.
 */

static int
ucaczone (region, star1, star2, verbose)

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

/* Set path to UCAC Catalog */
    if ((str = getenv("UCAC_PATH")) != NULL ) {
	lpath = strlen (str) + 16;
	tabpath = (char *) malloc (lpath);
	strcpy (tabpath, str);
	}
    else {
	lpath = strlen (ucaccd) + 16;
	tabpath = (char *) malloc (lpath);
	strcpy (tabpath, ucaccd);
	}

/* Set pathname for zone file */
    sprintf (temp,"/u1/z%03d", deczone);
    strcat (tabpath, temp);
    nb = getfilesize (tabpath);
    nstars = nb / 67;

    *star1 = 1;
    *star2 = nstars;
    free (tabpath);
    return (1);
}


/* UCACREG -- search the UCAC Catalog index table for fields
 * in the specified range of coordinates and magnitudes.
 * Build lists containing the first star and number of stars for each range.
 */

static int
ucacreg (ra1, ra2, dec1, dec2, nrmax, regnum, rstar, nstar, verbose)

double	ra1, ra2;	/* Right ascension limits in degrees */
double	dec1, dec2; 	/* Declination limits in degrees */
int	nrmax;		/* Maximum number of regions to find */
int	*regnum;	/* Region numbers (returned)*/
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

/* Set path to UCAC Catalog CDROM */
    if ((str = getenv("UCAC_PATH")) != NULL ) {
	tabpath = (char *) malloc (strlen (str) + 16);
	strcpy (tabpath, str);
	}
    else {
	tabpath = (char *) malloc (strlen (ucaccd) + 16);
	strcpy (tabpath, ucaccd);
	}

/* Set pathname for index table file */
    strcat (tabpath,"/data/index.dat");

/* Read the index table */
    if ((buffer = getfilebuff (tabpath)) == NULL) {
	fprintf (stderr,"UCACREG:  error reading region table %s\n",tabpath);
	return (0);
	}

/* Find region range to search based on declination */
    iz1 = ucaczone (dec1);
    iz2 = ucaczone (dec2);
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
	fprintf (stderr,"UCACREG: searching %d regions: %d - %d\n",nsrch,ir1,ir2);
    if (jr1 > 0) {
	nsrch1 = jr2 - jr1 + 1;
	if (verbose)
	    fprintf (stderr,"UCACREG: searching %d regions: %d - %d\n",nsrch1,jr1,jr2);
	}
    if (verbose)
	fprintf(stderr,"UCACREG: RA: %.5f - %.5f, Dec: %.5f - %.5f\n",ra1,ra2,dec1,dec2);

    nrgn = 0;

    for (iwrap = 0; iwrap < nwrap; iwrap++) {

	for (irow = ir1 - 1; irow < ir2; irow++) {

	/* Read next line of region table */
	    line = buffer + (irow * nchar);

	/* Declination range of the gs region */
	/* note:  southern dechi and declow are reversed */
	    num1 = atoi (line);
	    num2 = atoi (line+nchar);
	    dechi = atof (line + 29);
	    declow = atof (line + 36);
	    if (dechi > declow) {
		decmin = declow - 0.1;
		decmax = dechi + 0.1;
		}
	    else {
		decmax = declow + 0.1;
		decmin = dechi - 0.1;
		}

	    if (decmax >= dec1 && decmin <= dec2) {

	    /* Right ascension range of the Guide Star Catalog region */
		ralow = atof (line + 15) - 0.1;
		if (ralow <= 0.0) ralow = 0.0;
		rahi = atof (line + 22) + 0.1;
		if (rahi > 360.0) rahi = 360.0;
		if (rahi <= 0.0) rahi = 360.0;

	    /* Check RA if 0h RA not between region RA limits */
		if (ra1 < ra2) {

		    /* Add this region to list, if there is space */
		    if (ralow <= ra2 && rahi >= ra1) {
			if (verbose)
			    fprintf (stderr,"UCACREG: Region %d added to search\n",irow);
			if (nrgn < nrmax) {
			    regnum[nrgn] = irow;
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
			    regnum[nrgn] = irow;
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

 
 
/*  UCACZONE -- find the zone number where a declination can be found */
 
static int
ucaczone (dec)
 
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



/* UCACOPEN -- Open UCAC catalog file, returning number of entries */

struct StarCat *
ucacopen (nstar, nread)

int	nstar;	/* Number of first star to read */
int	nread;	/* Number of star entries to read */

{
    FILE *fcat;
    struct StarCat *sc;
    int lfile, lpath;
    int lread, lskip, nr;
    char *str;
    char *ucacfile;
    char *ucacpath;	/* Full pathname for catalog file */

    /* Set path to UCAC Catalog CDROM */
    if ((str = getenv("UCAC_PATH")) != NULL ) {
	lpath = strlen(str) + 18;
	ucacpath = (char *) malloc (lpath);
	strcpy (ucacpath, str);
	}
    else {
	lpath = strlen(ucaccd) + 18;
	ucacpath = (char *) malloc (lpath);
	strcpy (ucacpath, ucaccd);
	}

    /* Set pathname for catalog file */
    strcat (ucacpath, "/data/catalog.dat");

    /* Find length of UCAC catalog file */
    lfile = ucacsize (ucacpath);

    /* Check for existence of catalog */
    if (lfile < 2) {
	fprintf (stderr,"UCACOPEN: Binary catalog %s has no entries\n",ucacpath);
	free (ucacpath);
	return (NULL);
	}

    /* Open UCAC file */
    if (!(fcat = fopen (ucacpath, "r"))) {
	fprintf (stderr,"UCACOPEN: UCAC file %s cannot be read\n",ucacpath);
	free (ucacpath);
	return (0);
	}

    /* Set UCAC catalog header information */
    sc = (struct StarCat *) calloc (1, sizeof (struct StarCat));
    sc->byteswapped = 0;

    sc->nbent = 208;
    sc->nstars = lfile / sc->nbent;

    /* Separate filename from pathname and save in structure */
    ucacfile = strrchr (ucacpath,'/');
    if (ucacfile)
	ucacfile = ucacfile + 1;
    else
	ucacfile = ucacpath;
    if (strlen (ucacfile) < 24)
	strcpy (sc->isfil, ucacfile);
    else
	strncpy (sc->isfil, ucacfile, 23);

    /* Set other catalog information in structure */
    sc->inform = 'J';
    sc->coorsys = WCS_J2000;
    sc->epoch = 2000.0;
    sc->equinox = 2000.0;
    sc->ifcat = fcat;
    sc->sptype = 2;

    /* UCAC stars are not RA-sorted within regions */
    sc->rasorted = 0;

    /* Read part of catalog into a buffer */
    lread = nread * sc->nbent;
    lskip = (nstar - 1) * sc->nbent;
    sc->catdata = NULL;
    if ((sc->catdata = calloc (1, lread+1)) != NULL) {
	fseek (fcat, lskip, 0);
	nr = fread (sc->catdata, 1, lread, fcat);
	if (nr < lread) {
	    fprintf (stderr,"UCACOPEN: Read %d / %d bytes\n", nr, lread);
            ucacclose (sc);
	    free (ucacpath);
            return (NULL);
            }
	sc->catlast = sc->catdata + lread;
	}
    else {
	fprintf (stderr,"UCACOPEN: Cannot allocate %d-byte buffer.\n", lread);
        ucacclose (sc);
	free (ucacpath);
	return (NULL);
	}
    sc->istar = nstar;
    free (ucacpath);
    return (sc);
}


void
ucacclose (sc)
struct StarCat *sc;	/* Star catalog descriptor */
{
    fclose (sc->ifcat);
    if (sc->catdata != NULL)
	free (sc->catdata);
    free (sc);
    return;
}


/* UCACSTAR -- Get UCAC catalog entry for one star;
              return 0 if successful */

static int
ucacstar (sc, st, istar)

struct StarCat *sc;	/* Star catalog descriptor */
struct Star *st;	/* Current star entry */
int istar;	/* Star sequence number in UCAC catalog region file */
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
	fprintf (stderr, "UCACSTAR:  %d  > %d is not in catalog\n",
		 istar, sc->nstars);
	return (3);
	}

    /* Move buffer pointer to start of correct star entry */
    if (istar > 0) {
	line = sc->catdata + ((istar - sc->istar) * sc->nbent);
	if (line >= sc->catlast) {
	    fprintf (stderr, "UCACSTAR:  star %d past buffer\n", istar);
	    return (4);
	    }
	}

    /* Read catalog entry */
    if (sc->nbent > sc->catlast-line) {
	fprintf (stderr, "UCACSTAR:  %d / %d bytes read\n",
		 sc->catlast - line, sc->nbent);
	return (5);
	}

    st->num = regnum + (0.000001 * (double) istar);

    /* Read position in degrees */
    st->ra = atof (line) / 3600000.0;
    st->dec = (atof (line+10) / 3600000.0) - 90.0;

    /* Read proper motion and convert it to to degrees/year */
    st->rapm = (atof (line+41) / 3600000.0) / cosdeg (st->dec);
    st->decpm = atof (line+48) / 3600000.0;

    /* Set V magnitude */
    st->xmag[0] = atof (line+20) * 0.01;

    return (0);
}

/* UCACSIZE -- return size of UCAC catalog file in bytes */

static int
ucacsize (filename)

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


/* Apr 10 2003	New subroutines, based on ty2read.c
 * Apr 14 2003	Explicitly get revision date if nstarmax < 1
 */
