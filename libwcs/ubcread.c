/*** File libwcs/ubcread.c
 *** November 25, 2002
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 *** Copyright (C) 2002
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

 * Subroutines to read from the USNO-B1.0 catalog
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "wcs.h"
#include "wcscat.h"

static int ucat=UB1;

/* USNO B-1.0 directory pathname; replaced by UB1_PATH environment variable
 * Use this if CDROMs have been transferred to a single hard disk
 * Otherwise set to null string ("") and use cdroot
 * This may also be a URL to a catalog search engine */
static char ub1path[64]="/data/ub1";

/* Uncomment following line to use ESO USNO-B server for UB1
static char ub1path[64]="http://archive.eso.org/skycat/servers/usnob-server";
 */

static char *ubpath;

typedef struct {
    int rasec, decsec, pm, pmerr, poserr, mag[5], magerr[5], index[5];
} UBCstar;

static int nstars;	/* Number of stars in catalog */
static int cswap = 0;	/* Byte reverse catalog to Intel/DEC order if 1 */
static double *udist;	/* Array of distances to stars */
static int ndist = 0;

static FILE *fcat;
#define ABS(a) ((a) < 0 ? (-(a)) : (a))
#define NZONES 24

static double ubcra();
static double ubcdec();
static double ubcmagr();
static double ubcmagb();
static int ubcmagerr();
static int ubcgsc();
static int ubcplate();

static int ubczones();
static int ubczone();
static int ubcsra();
static int ubcopen();
static int ubcpath();
static int ubcstar();
static void ubcswap();
static int nbent = 80;

static int xplate = 0;	/* If nonzero, use objects only from this plate */
void setuplate (xplate0)
int xplate0;
{ xplate = xplate0; return; }
int getuplate ()
{ return (xplate); }

/* UBCREAD -- Read USNO B Catalog stars */

int
ubcread (refcatname,distsort,cra,cdec,dra,ddec,drad,sysout,eqout,epout,
	 mag1,mag2,sortmag,nstarmax,unum,ura,udec,umag,uplate,nlog)

char	*refcatname;	/* Name of catalog (UBC, USAC, UBC2, USAC2) */
int	distsort;	/* 1 to sort stars by distance from center */
double	cra;		/* Search center J2000 right ascension in degrees */
double	cdec;		/* Search center J2000 declination in degrees */
double	dra;		/* Search half width in right ascension in degrees */
double	ddec;		/* Search half-width in declination in degrees */
double	drad;		/* Limiting separation in degrees (ignore if 0) */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
double	mag1,mag2;	/* Limiting magnitudes (none if equal) */
int	sortmag;	/* Magnitude by which to sort (1 or 2) */
int	nstarmax;	/* Maximum number of stars to be returned */
double	*unum;		/* Array of UB numbers (returned) */
double	*ura;		/* Array of right ascensions (returned) */
double	*udec;		/* Array of declinations (returned) */
double	**umag;		/* Array of red and blue magnitudes (returned) */
int	*uplate;	/* Array of plate numbers (returned) */
int	nlog;		/* Logging interval */
{
    double ra1,ra2;	/* Limiting right ascensions of region in degrees */
    double dec1,dec2;	/* Limiting declinations of region in degrees */
    int nz;		/* Number of input UB zone files */
    int zlist[NZONES];	/* List of input UB zones */
    UBCstar star;	/* UB catalog entry for one star */
    double dist = 0.0;	/* Distance from search center in degrees */
    double faintmag=0.0; /* Faintest magnitude */
    double maxdist=0.0; /* Largest distance */
    int	faintstar=0;	/* Faintest star */
    int	farstar=0;	/* Most distant star */
    int sysref=WCS_J2000;	/* Catalog coordinate system */
    double eqref=2000.0;	/* Catalog equinox */
    double epref=2000.0;	/* Catalog epoch */

    double rra1, rra2, rdec1, rdec2;
    double num;		/* UB numbers */
    int wrap, iwrap, nnfld;
    int verbose;
    int znum, itot,iz, i;
    int itable, jtable,jstar;
    int nstar, nread;
    int uara1, uara2, uadec1, uadec2;
    double ra,dec;
    double mag, magb, magr;
    int istar, istar1, istar2, plate;
    int nzmax = NZONES;	/* Maximum number of declination zones */
    int isp;
    int magsort;
    char ispc[2];
    char *str;
    char cstr[32], rastr[32], numstr[32], decstr[32], catid[32];
    char *title;

    itot = 0;
    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    /* Set catalog code and path to catalog */
    if (strncmp (refcatname,"us",2)==0 ||
        strncmp (refcatname,"US",2)==0) {
	if (strchr (refcatname, '2') != NULL) {
	    if ((str = getenv("USA2_PATH")) != NULL)
		strcpy (usa2path,str);
	    ucat = USA2;
	    ubpath = usa2path;
	    }
	else {
	    if ((str = getenv("USA1_PATH")) != NULL)
		strcpy (usa1path,str);
	    ucat = USA1;
	    ubpath = usa1path;
	    }
	}
    else if (strncmp (refcatname,"ua",2)==0 ||
        strncmp (refcatname,"UB",2)==0) {
	if (strchr (refcatname, '2') != NULL) {
	    if ((str = getenv("UB2_PATH")) != NULL)
		strcpy (ua2path,str);
	    else if ((str = getenv("UB2_ROOT")) != NULL) {
		ua2path[0] = 0;
		strcpy (cdroot,str);
		}
	    ucat = UB2;
	    ubpath = ua2path;
	    }
	else {
	    if ((str = getenv("UB1_PATH")) != NULL)
		strcpy (ua1path,str);
	    else if ((str = getenv("UB1_ROOT")) != NULL) {
		ua1path[0] = 0;
		strcpy (cdroot,str);
		}
	    ucat = UB1;
	    ubpath = ua1path;
	    }
	}
    else {
	fprintf (stderr, "UBCREAD:  %s not a USNO catalog\n", refcatname);
	return (0);
	}

    /* If root pathname is a URL, search and return */
    if (!strncmp (ubpath, "http:",5)) {
	return (webread (ubpath,refcatname,distsort,cra,cdec,dra,ddec,drad,
			 sysout,eqout,epout,mag1,mag2,sortmag,nstarmax,
			 unum,ura,udec,NULL,NULL,umag,uplate,nlog));
	}

    wcscstr (cstr, sysout, eqout, epout);

    SearchLim (cra,cdec,dra,ddec,sysout,&ra1,&ra2,&dec1,&dec2,verbose);

    /* mag1 is always the smallest magnitude */
    if (mag2 < mag1) {
	mag = mag2;
	mag2 = mag1;
	mag1 = mag;
	}
    if (sortmag == 1)
	magsort = 0;
    else
	magsort = 1;

    /* Find UB Star Catalog regions in which to search */
    rra1 = ra1;
    rra2 = ra2;
    rdec1 = dec1;
    rdec2 = dec2;
    RefLim (cra, cdec, dra, ddec, sysout, sysref, eqout, eqref, epout,
	    &rra1, &rra2, &rdec1, &rdec2, verbose);
    nz = ubczones (rra1, rra2, rdec1, rdec2, nzmax, zlist, verbose);
    if (nz <= 0) {
	fprintf (stderr, "UBCREAD:  no USNO A zones found\n");
	return (0);
	}

    /* Write header if printing star entries as found */
    if (nstarmax < 1) {
	title = CatName (ucat, refcatname);
	printf ("catalog        %s\n", title);
	free ((char *)title);
	ra2str (rastr, 31, cra, 3);
	printf ("ra     %s\n", rastr);
	dec2str (decstr, 31, cdec, 2);
	printf ("dec    %s\n", decstr);
	if (drad != 0.0)
	    printf ("radmin     %.1f\n", drad*60.0);
	else {
	    printf ("dramin     %.1f\n", dra*60.0* cosdeg (cdec));
	    printf ("ddecmin    %.1f\n", ddec*60.0);
	    }
	printf ("radecsys       %s\n", cstr);
	printf ("equinox        %.3f\n", eqout);
	printf ("epoch  %.3f\n", epout);
	CatID (catid, ucat);
	printf ("%s	ra          	dec         	", catid);
	printf ("magb 	magr  	arcmin\n");
	printf ("-------------	------------	------------    ");
	printf ("-----	-----	------\n");
	}

    /* If RA range includes zero, set a flat */
    wrap = 0;
    if (rra1 > rra2)
	wrap = 1;
    else
	wrap = 0;

    uara1 = (int) (rra1 * 360000.0 + 0.5);
    uara2 = (int) (rra2 * 360000.0 + 0.5);
    uadec1 = (int) ((rdec1 * 360000.0) + 32400000.5);
    uadec2 = (int) ((rdec2 * 360000.0) + 32400000.5);
    
    if (nstarmax > ndist) {
	if (ndist > 0)
	    free ((void *)udist);
	udist = (double *) malloc (nstarmax * sizeof (double));
	if (udist == NULL) {
	    fprintf (stderr,"UBCREAD:  cannot allocate separation array\n");
	    return (0);
	    }
	ndist = nstarmax;
	}

    /* Loop through region list */
    nstar = 0;
    for (iz = 0; iz < nz; iz++) {

    /* Get path to zone catalog */
	znum = zlist[iz];
	if ((nstars = ubcopen (znum)) != 0) {

	    jstar = 0;
	    jtable = 0;
	    for (iwrap = 0; iwrap <= wrap; iwrap++) {

	    /* Find first star based on RA */
		if (iwrap == 0 || wrap == 0)
		    istar1 = ubcsra (rra1);
		else
		    istar1 = 1;

	    /* Find last star based on RA */
		if (iwrap == 1 || wrap == 0)
		    istar2 = ubcsra (rra2);
		else
		    istar2 = nstars;

		if (istar1 == 0 || istar2 == 0)
		    break;

		nread = istar2 - istar1 + 1;
		itable = 0;

	    /* Loop through zone catalog for this region */
		for (istar = istar1; istar <= istar2; istar++) {
		    itable ++;
		    jtable ++;

		    if (ubcstar (istar, &star)) {
			fprintf (stderr,"UBCREAD: Cannot read star %d\n", istar);
			break;
			}

		/* Extract selected fields */
		    else {

		    /* Check position limits */
     			if ((star.decsec >= uadec1 && star.decsec <= uadec2) &&
			    ((wrap && (star.rasec>=uara1 || star.rasec<=uara2)) ||
			     (!wrap && (star.rasec>=uara1 && star.rasec<=uara2))
			    )){

			/* Check magnitude, distance, and plate number */
			    magb = ubcmagb (star.magetc);
			    magr = ubcmagr (star.magetc);
			    if (magsort == 1)
				mag = magr;
			    else
				mag = magb;
			    plate = ubcplate (star.magetc);
			    ra = ubcra (star.rasec);
			    dec = ubcdec (star.decsec);
			    wcscon (sysref,sysout,eqref,eqout,&ra,&dec,epout);
			    if (distsort || drad > 0)
				dist = wcsdist (cra,cdec,ra,dec);
			    else
				dist = 0.0;
			    if ((mag1==mag2 || (mag>=mag1 && mag<=mag2)) &&
				(drad == 0.0 || dist < drad) &&
				(xplate == 0 || plate == xplate)) {

				/* br2sp (NULL, magb, mag, ispc);
				isp = (1000 * (int)ispc[0]) + (int)ispc[1]; */
				num = (double) znum +
				      (0.00000001 * (double)istar);

			    /* Write star position and magnitudes to stdout */
				if (nstarmax < 1) {
				    CatNum (ucat, -13, 0, num, numstr);
				    ra2str (rastr, 31, ra, 3);
				    dec2str (decstr, 31, dec, 2);
				    dist = wcsdist (cra,cdec,ra,dec) * 60.0;
				    printf ("%s	%s	%s", numstr,rastr,decstr);
				    printf ("	%.2f	%.2f	%.2f\n",
					magb, magr, dist);
				    }

			    /* Save star position and magnitude in table */
				else if (nstar < nstarmax) {
				    unum[nstar] = num;
				    ura[nstar] = ra;
				    udec[nstar] = dec;
				    umag[0][nstar] = magb;
				    umag[1][nstar] = magr;
				    uplate[nstar] = plate;
				    udist[nstar] = dist;
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
				replace furthest star */
				else if (distsort) {
				    if (dist < maxdist) {
					unum[farstar] = num;
					ura[farstar] = ra;
					udec[farstar] = dec;
					umag[0][farstar] = magb;
					umag[1][farstar] = magr;
					uplate[farstar] = plate;
					udist[farstar] = dist;

				    /* Find new farthest star */
					maxdist = 0.0;
					for (i = 0; i < nstarmax; i++) {
					    if (udist[i] > maxdist) {
						maxdist = udist[i];
						farstar = i;
						}
					    }
					}
				    }

			    /* If too many stars, replace faintest star */
				else if (mag < faintmag) {
				    unum[faintstar] = num;
				    ura[faintstar] = ra;
				    udec[faintstar] = dec;
				    umag[0][faintstar] = magb;
				    umag[1][faintstar] = magr;
				    uplate[faintstar] = plate;
				    udist[faintstar] = dist;

			    /* Find new faintest star */
				    faintmag = 0.0;
				    for (i = 0; i < nstarmax; i++) {
					if (umag[magsort][i] > faintmag) {
					    faintmag = umag[magsort][i];
					    faintstar = i;
					    }
					}
				    }
				nstar++;
				jstar++;
				if (nlog == 1)
				    fprintf (stderr,"UBCREAD: %04d.%08d: %9.5f %9.5f %s %5.2f %5.2f\n",
					znum,istar,ra,dec,cstr,magb,magr);

			    /* End of accepted star processing */
				}
			    }

		    /* End of individual star processing */
			}

		/* Log operation */
		    if (nlog > 0 && itable%nlog == 0)
			fprintf (stderr,"UBCREAD: zone %d (%2d / %2d) %8d / %8d / %8d sources\r",
				znum, iz+1, nz, jstar, itable, nread);

		/* End of star loop */
		    }

		/* End of wrap loop */
		}

	/* Close zone input file */
	    (void) fclose (fcat);
	    itot = itot + itable;
	    if (nlog > 0)
		fprintf (stderr,"UBCREAD: zone %d (%2d / %2d) %8d / %8d / %8d sources      \n",
			znum, iz+1, nz, jstar, jtable, nstars);

	/* End of zone processing */
	    }

    /* End of zone loop */
	}

/* Summarize search */
    if (nlog > 0) {
	if (nz > 1)
	    fprintf (stderr,"UBCREAD: %d zones: %d / %d found\n",nz,nstar,itot);
	else
	    fprintf (stderr,"UBCREAD: 1 zone: %d / %d found\n",nstar,itable);
	if (nstar > nstarmax)
	    fprintf (stderr,"UBCREAD: %d stars found; only %d returned\n",
		     nstar,nstarmax);
	}
    return (nstar);
}


int
ubcrnum (refcatname,nnum,sysout,eqout,epout,unum,ura,udec,umag,uplate,nlog)

char	*refcatname;	/* Name of catalog (UBC, USAC, UBC2, USAC2) */
int	nnum;		/* Number of stars to find */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
double	*unum;		/* Array of UB numbers to find */
double	*ura;		/* Array of right ascensions (returned) */
double	*udec;		/* Array of declinations (returned) */
double	**umag;		/* Array of blue and red magnitudes (returned) */
int	*uplate;	/* Array of plate numbers (returned) */
int	nlog;		/* Logging interval */
{
    UBCstar star;	/* UB catalog entry for one star */
    int sysref=WCS_J2000;	/* Catalog coordinate system */
    double eqref=2000.0;	/* Catalog equinox */
    int isp;
    char ispc[2];

    int znum;
    int jnum;
    int nzone;
    int nfound = 0;
    double ra,dec;
    double magr, magb;
    double dstar;
    int istar, plate;
    char *str;

    /* Set catalog code and path to catalog */
    if (strncmp (refcatname,"us",2)==0 ||
        strncmp (refcatname,"US",2)==0) {
	if (strchr (refcatname, '2') != NULL) {
	    if ((str = getenv("USA2_PATH")) != NULL)
		strcpy (usa2path,str);
	    ucat = USA2;
	    ubpath = usa2path;
	    }
	else {
	    if ((str = getenv("USA1_PATH")) != NULL)
		strcpy (usa1path,str);
	    ucat = USA1;
	    ubpath = usa1path;
	    }
	}
    else if (strncmp (refcatname,"ua",2)==0 ||
        strncmp (refcatname,"UB",2)==0) {
	if (strchr (refcatname, '2') != NULL) {
	    if ((str = getenv("UB2_PATH")) != NULL)
		strcpy (ua2path,str);
	    else if ((str = getenv("UB2_ROOT")) != NULL) {
		ua2path[0] = 0;
		strcpy (cdroot,str);
		}
	    ucat = UB2;
	    ubpath = ua2path;
	    }
	else {
	    if ((str = getenv("UB1_PATH")) != NULL)
		strcpy (ua1path,str);
	    else if ((str = getenv("UB1_ROOT")) != NULL) {
		ua1path[0] = 0;
		strcpy (cdroot,str);
		}
	    ucat = UB1;
	    ubpath = ua1path;
	    }
	}
    else {
	fprintf (stderr, "UBCREAD:  %s not a USNO catalog\n", refcatname);
	return (0);
	}

    /* If root pathname is a URL, search and return */
    if (!strncmp (ubpath, "http:",5)) {
	return (webrnum (ubpath,refcatname,nnum,sysout,eqout,epout,
			 unum,ura,udec,NULL,NULL,umag,uplate,nlog));
	}


/* Loop through star list */
    for (jnum = 0; jnum < nnum; jnum++) {

    /* Get path to zone catalog */
	znum = (int) unum[jnum];
	if ((nzone = ubcopen (znum)) != 0) {
	    dstar = (unum[jnum] - znum) * 100000000.0;
	    istar = (int) (dstar + 0.5);
	    if (istar > nzone) {
		fprintf (stderr,"UBCRNUM: Star %d > max. in zone %d\n",
			 istar,nzone);
		break;
		}

	    if (ubcstar (istar, &star)) {
		fprintf (stderr,"UBCRNUM: Cannot read star %d\n", istar);
		break;
		}

	    /* Extract selected fields */
	    else {
		ra = ubcra (star.rasec); /* Right ascension in degrees */
		dec = ubcdec (star.decsec); /* Declination in degrees */
		magb = ubcmagb (star.magetc); /* Blue magnitude */
		magr = ubcmagr (star.magetc); /* Red magnitude */
		plate = ubcplate (star.magetc);	/* Plate number */
		wcscon (sysref, sysout, eqref, eqout, &ra, &dec, epout);

		/* Save star position and magnitude in table */
		ura[nfound] = ra;
		udec[nfound] = dec;
		umag[0][nfound] = magb;
		umag[1][nfound] = magr;
		/* br2sp (NULL, magb, magr, ispc);
		isp = (1000 * (int)ispc[0]) + (int)ispc[1]; */
		uplate[nfound] = plate;

		nfound++;
		if (nlog == 1)
		    fprintf (stderr,"UBCRNUM: %04d.%08d: %9.5f %9.5f %5.2f %5.2f\n",
			     znum,istar,ra,dec,magb,magr);

		/* Log operation */
		if (nlog > 0 && jnum%nlog == 0)
		    fprintf (stderr,"UBCRNUM: %4d.%8d  %8d / %8d sources\r",
			     znum, istar, jnum, nnum);

		(void) fclose (fcat);
		/* End of star processing */
		}

	    /* End of star */
	    }

	/* End of star loop */
	}

    /* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"UBCRNUM:  %d / %d found\n",nfound,nnum);

    return (nfound);
}


/* Declination zone numbers */
int azone[NZONES]={0,75,150,225,300,375,450,525,600,675,750,825,900,
	      975,1050,1125,1200,1275,1350,1425,1500,1575,1650,1725};

/* UBCZONES -- figure out which UB zones will need to be searched */

static int
ubczones (ra1, ra2, dec1, dec2, nzmax, zones, verbose)

double	ra1, ra2;	/* Right ascension limits in degrees */
double	dec1, dec2; 	/* Declination limits in degrees */
int	nzmax;		/* Maximum number of zones to find */
int	*zones;		/* Region numbers (returned)*/
int	verbose;	/* 1 for diagnostics */

{
    int nrgn;		/* Number of zones found (returned) */
    int iz,iz1,iz2,i;

    for (i = 0; i < nzmax; i++)
	zones[i] = 0;

    nrgn = 0;

/* Find zone range to search based on declination */
    iz1 = ubczone (dec1);
    iz2 = ubczone (dec2);

/* Tabulate zones to search */
    i = 0;
    if (iz2 >= iz1) {
	for (iz = iz1; iz <= iz2; iz++)
	    zones[i++] = azone[iz];
	}
    else {
	for (iz = iz2; iz <= iz1; iz++)
	    zones[i++] = azone[iz];
	}

    nrgn = i;
    if (verbose) {
	fprintf(stderr,"UBCZONES:  %d zones: %d - %d\n",nrgn,zones[0],zones[i-1]);
	fprintf(stderr,"UBCZONES: RA: %.5f - %.5f, Dec: %.5f - %.5f\n",ra1,ra2,dec1,dec2);
	}

    return (nrgn);
}


/* UBCRA -- returns right ascension in degrees from the UB star structure */

static double
ubcra (rasec)

int rasec;	/* RA in 100ths of arcseconds from UB catalog entry */
{
    return ((double) (rasec) / 360000.0);
}


/* UBCDEC -- returns the declination in degrees from the UB star structure */

static double
ubcdec (decsec)

int decsec;	/* Declination in 100ths of arcseconds from UB catalog entry */
{
    return ((double) (decsec - 32400000) / 360000.0);
}


/* UBCMAGR -- returns the red magnitude from the UB star structure */

static double
ubcmagr (magetc)

int magetc;	/* Quality, plate, and magnitude from UB catalog entry */
{
    if (magetc < 0)
	return ((double) (-magetc % 1000) * 0.1);
    else
	return ((double) (magetc % 1000) * 0.1);
}


/* UBCMAGB -- returns the blue magnitude from the UB star structure */

static double
ubcmagb (magetc)

int magetc;	/* Quality, plate, and magnitude from UB catalog entry */
{
    if (magetc < 0)
	return ((double) ((-magetc / 1000) % 1000) * 0.1);
    else
	return ((double) ((magetc / 1000) % 1000) * 0.1);
}


/* UBCMAGERR -- returns 1 if magnitude is uncertain from UB star structure */

static int
ubcmagerr (magetc)

int magetc;	/* Quality, plate, and magnitude from UB catalog entry */
{
    if (magetc < 0)
	return ((-magetc / 1000000000) % 10);
    else
	return ((magetc / 1000000000) % 10);
}


/* UBCTY2 -- returns 1 if UB star is in the Tycho-2 Catalog */

static int
ubcty2 (star)

int *star;	/* Quality, plate, and magnitude from UB catalog entry */
{
    int num = star->pmerr;
    if (num < 0)
	return ((-num / 100000000) % 1);
    else
	return ((num / 100000000) % 1);
}


/* UBCPLATE -- returns the plate number from the UB star structure */

static int
ubcplate (magetc)

int magetc;	/* Quality, plate, and magnitude from UB catalog entry */
{
    if (magetc < 0)
	return ((-magetc / 100000000) % 1);
    else
	return ((magetc / 100000000) % 1);
}


/* UBCZONE -- find the UB zone number where a declination can be found */

static int
ubczone (dec)

double dec;	/* declination in degrees */
{
    double zonesize = 7.5;	/* number of degrees per declination zone */
    int zone;

    zone = (int) ((dec + 90.0) / zonesize);
    if (zone > 23)
	zone = 23;
    else if (zone < 0)
	zone = 0;
    return (zone);
}


/* UBCSRA -- Find UB star closest to specified right ascension */

static int
ubcsra (rax0)

double	rax0;		/* Right ascension in degrees for which to search */
{
    int istar, istar1, istar2, nrep;
    double rax, ra1, ra, rdiff, rdiff1, rdiff2, sdiff;
    UBCstar star;	/* UB catalog entry for one star */
    char rastrx[16];
    int debug = 0;

    rax = rax0;
    if (debug)
	ra2str (rastrx, 16, rax, 3);
    istar1 = 1;
    if (ubcstar (istar1, &star))
	return (0);
    ra1 = ubcra (star.rasec);
    istar = nstars;
    nrep = 0;
    while (istar != istar1 && nrep < 20) {
	if (ubcstar (istar, &star))
	    break;
	else {
	    ra = ubcra (star.rasec);
	    if (ra == ra1)
		break;
	    if (debug) {
		char rastr[16];
		ra2str (rastr, 16, ra, 3);
		fprintf (stderr,"UBCSRA %d %d: %s (%s)\n",
			 nrep,istar,rastr,rastrx);
		}
	    rdiff = ra1 - ra;
	    rdiff1 = ra1 - rax;
	    rdiff2 = ra - rax;
	    if (nrep > 20 && ABS(rdiff2) > ABS(rdiff1)) {
		istar = istar1;
		break;
		}
	    nrep++;
	    sdiff = (double)(istar - istar1) * rdiff1 / rdiff;
	    istar2 = istar1 + (int) (sdiff + 0.5);
	    ra1 = ra;
	    istar1 = istar;
	    istar = istar2;
	    if (debug) {
		fprintf (stderr," ra1=    %.5f ra=     %.5f rax=    %.5f\n",
			 ra1,ra,rax);
		fprintf (stderr," rdiff=  %.5f rdiff1= %.5f rdiff2= %.5f\n",
			 rdiff,rdiff1,rdiff2);
		fprintf (stderr," istar1= %d istar= %d istar1= %d\n",
			 istar1,istar,istar2);
		}
	    if (istar < 1)
		istar = 1;
	    if (istar > nstars)
		istar = nstars;
	    if (istar == istar1)
		break;
	    }
	}
    return (istar);
}

/* UBCOPEN -- Open UB Catalog zone catalog, returning number of entries */

static int
ubcopen (znum)

int znum;	/* UB Catalog zone */
{
    char zonepath[64];	/* Pathname for input UB zone file */
    UBCstar star;	/* UB catalog entry for one star */
    struct stat statbuff;
    
/* Get path to zone catalog */
    if (ubcpath (znum, zonepath)) {
	fprintf (stderr, "UBCOPEN: Cannot find zone catalog for %d\n", znum);
	return (0);
	}

/* Find number of stars in zone catalog by its length */
    if (stat (zonepath, &statbuff)) {
	fprintf (stderr,"UB zone catalog %s has no entries\n",zonepath);
	return (0);
	}
    else
	nstars = (int) statbuff.st_size / nbent;

/* Open zone catalog */
    if (!(fcat = fopen (zonepath, "r"))) {
	fprintf (stderr,"UB zone catalog %s cannot be read\n",zonepath);
	return (0);
	}

/* Check to see if byte-swapping is necessary */
    cswap = 0;
    if (ubcstar (1, &star)) {
	fprintf (stderr,"UBCOPEN: cannot read star 1 from UB zone catalog %s\n",
		 zonepath);
	return (0);
	}
    else {
	if (star.rasec > 360 * 360000 || star.rasec < 0) {
	    cswap = 1;
	    /* fprintf (stderr,"UBCOPEN: swapping bytes in UB zone catalog %s\n",
		     zonepath); */
	    }
	else if (star.decsec > 180 * 360000 || star.decsec < 0) {
	    cswap = 1;
	    /* fprintf (stderr,"UBCOPEN: swapping bytes in UB zone catalog %s\n",
		     zonepath); */
	    }
	else
	    cswap = 0;
	}

    return (nstars);
}


/* UBCPATH -- Get UB Catalog region file pathname */

static int
ubcpath (zn, path)

int zn;		/* UB zone number */
char *path;	/* Pathname of UB zone file */

{
    int iz;		/* Zone index (0000 = 0, 0075 = 1, ...) */

    /* Return error code and null path if zone is out of range */
    if (zn < 0 || zn > 1725) {
	fprintf (stderr, "UBCPATH: zone %d out of range 0-1725\n",zn);
	path[0] = 0;
	return (-1);
	}

    /* Set path for USNO SA zone catalog */
    if (ucat == USA1 || ucat == USA2)
	sprintf (path,"%s/zone%04d.cat", ubpath, zn);

    /* Set zone catalog path when USNO A is in a single directory */
    else if (strlen (ubpath) > 0)
	sprintf (path,"%s/zone%04d.cat", ubpath, zn);

    else {
	return (-1);
	}

    return (0);
}


/* UBCSTAR -- Get UB catalog entry for one star; return 0 if successful */

static int
ubcstar (istar, star)

int istar;	/* Star sequence number in UB zone catalog */
UBCstar *star;	/* UB catalog entry for one star */
{
    int nbs, nbr, nbskip;

    if (istar < 1 || istar > nstars) {
	fprintf (stderr, "UBCstar %d is not in catalog\n",istar);
	return (-1);
	}
    nbskip = nbent * (istar - 1);
    if (fseek (fcat,nbskip,SEEK_SET))
	return (-1);
    nbs = sizeof (UBCstar);
    nbr = fread (star, nbs, 1, fcat) * nbs;
    if (nbr < nbs) {
	fprintf (stderr, "UBCstar %d / %d bytes read\n",nbr, nbs);
	return (-2);
	}
    if (cswap)
	ubcswap ((char *)star);
    return (0);
}


/* UBCSWAP -- Reverse bytes of UB Catalog entry */

static void
ubcswap (string)

char *string;	/* Start of vector of 4-byte ints */

{
char *sbyte, *slast;
char temp0, temp1, temp2, temp3;
int nbytes = nbent; /* Number of bytes to reverse */

    slast = string + nbytes;
    sbyte = string;
    while (sbyte < slast) {
	temp3 = sbyte[0];
	temp2 = sbyte[1];
	temp1 = sbyte[2];
	temp0 = sbyte[3];
	sbyte[0] = temp0;
	sbyte[1] = temp1;
	sbyte[2] = temp2;
	sbyte[3] = temp3;
	sbyte = sbyte + 4;
	}
    return;
}

/* Nov 25 2002	New subroutine based on ubcread.c
 */
