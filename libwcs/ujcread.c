/*** File libwcs/ujcread.c
 *** June 24, 1998
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "fitshead.h"
#include "wcs.h"

static char cdu[64]="/data/ujcat/catalog"; /* pathname of UJ 1.0 CDROM */

typedef struct {
    int rasec, decsec, magetc;
} UJCstar;

static int nstars;	/* Number of stars in catalog */
static int cswap = 0;	/* Byte reverse catalog to Intel/DEC order if 1 */
static FILE *fcat;
#define ABS(a) ((a) < 0 ? (-(a)) : (a))
#define NZONES 24

static double ujcra();
static double ujcdec();
static double ujcmag();
static int ujcplate();
static int ujczones();
static int ujczone();
static int ujcsra();
static int ujcopen();
static int ujcpath();
static int ujcstar();
static void ujcswap();

/* UJCREAD -- Read USNO J Catalog stars from CDROM */

int
ujcread (cra,cdec,dra,ddec,drad,mag1,mag2,xplate,nstarmax,unum,ura,udec,umag,
	 uplate,verbose)

double	cra;		/* Search center J2000 right ascension in degrees */
double	cdec;		/* Search center J2000 declination in degrees */
double	dra;		/* Search half width in right ascension in degrees */
double	ddec;		/* Search half-width in declination in degrees */
double	drad;		/* Limiting separation in degrees (ignore if 0) */
double	mag1,mag2;	/* Limiting magnitudes (none if equal) */
int	xplate;		/* If nonzero, use objects only from this plate */
int	nstarmax;	/* Maximum number of stars to be returned */
double	*unum;		/* Array of UJ numbers (returned) */
double	*ura;		/* Array of right ascensions (returned) */
double	*udec;		/* Array of declinations (returned) */
double	*umag;		/* Array of magnitudes (returned) */
int	*uplate;	/* Array of plate numbers (returned) */
int	verbose;	/* 1 for diagnostics */
{
    double ra1,ra2;	/* Limiting right ascensions of region in degrees */
    double dec1,dec2;	/* Limiting declinations of region in degrees */
    double dist = 0.0;  /* Distance from search center in degrees */
    double faintmag=0.0; /* Faintest magnitude */
    double maxdist=0.0; /* Largest distance */
    int	faintstar=0;	/* Faintest star */
    int	farstar=0;	/* Most distant star */
    double *udist;	/* Array of distances to stars */
    int nz;		/* Number of input UJ zone files */
    int zlist[NZONES];	/* List of input UJ zones */
    UJCstar star;	/* UJ catalog entry for one star */

    int wrap, iwrap;
    int znum, itot,iz;
    int nlog,itable,jstar;
    int nstar, i;
    double ra,dec;
    double mag;
    int istar, istar1, istar2, plate;
    int nzmax = NZONES;	/* Maximum number of declination zones */
    char *str;

    itot = 0;

    /* Set path to USNO A Catalog */
    if ((str = getenv("UJ_PATH")) != NULL )
	strcpy (cdu,str);

    /* Set right ascension limits for search */
    ra1 = cra - dra;
    ra2 = cra + dra;

    /* Keep right ascension between 0 and 360 degrees */
    if (ra1 < 0.0)
	ra1 = ra1 + 360.0;
    if (ra2 > 360.0)
	ra2 = ra2 - 360.0;

    /* Set declination limits for search */
    dec1 = cdec - ddec;
    dec2 = cdec + ddec;

    /* dec1 is always the smallest declination */
    if (dec1 > dec2) {
	dec = dec1;
	dec1 = dec2;
	dec2 = dec;
	}

    /* Search zones which include the poles cover 360 degrees in RA */
    if (dec1 < -90.0) {
	dec1 = -90.0;
	ra1 = 0.0;
	ra2 = 359.99999;
	}
    if (dec2 > 90.0) {
	dec2 = 90.0;
	ra1 = 0.0;
	ra2 = 359.99999;
	}
    if (verbose) {
	char rstr1[16],rstr2[16],dstr1[16],dstr2[16];
	ra2str (rstr1, 16, ra1, 3);
        dec2str (dstr1, 16, dec1, 2);
	ra2str (rstr2, 16, ra2, 3);
        dec2str (dstr2, 16, dec2, 2);
	fprintf (stderr,"UJCREAD: RA: %s - %s  Dec: %s - %s\n",
		 rstr1,rstr2,dstr1,dstr2);
	}

/* If RA range includes zero, split it in two */
    wrap = 0;
    if (ra1 > ra2)
	wrap = 1;
    else
	wrap = 0;

/* mag1 is always the smallest magnitude */
    if (mag2 < mag1) {
	mag = mag2;
	mag2 = mag1;
	mag1 = mag;
	}

/* Find UJ Star Catalog regions in which to search */
    nz = ujczones (ra1, ra2, dec1, dec2, nzmax, zlist, verbose);
    if (nz <= 0) {
	fprintf (stderr,"UJCREAD:  no UJ zones found\n");
	return (0);
	}

    udist = (double *) malloc (nstarmax * sizeof (double));

/* Logging interval */
    nlog = 0;
    nstar = 0;

/* Loop through region list */
    for (iz = 0; iz < nz; iz++) {

    /* Get path to zone catalog */
	znum = zlist[iz];
	if ((nstars = ujcopen (znum)) != 0) {

	    jstar = 0;
	    itable = 0;
	    for (iwrap = 0; iwrap <= wrap; iwrap++) {

	    /* Find first star based on RA */
		if (iwrap == 0 || wrap == 0)
		    istar1 = ujcsra (ra1);
		else
		    istar1 = 1;

	    /* Find last star based on RA */
		if (iwrap == 1 || wrap == 0)
		    istar2 = ujcsra (ra2);
		else
		    istar2 = nstars;

		if (istar1 == 0 || istar2 == 0)
		    break;

	    /* Loop through zone catalog for this region */
		for (istar = istar1; istar <= istar2; istar++) {
		    itable ++;

		    if (ujcstar (istar, &star)) {
			fprintf (stderr,"UJCREAD: Cannot read star %d\n", istar);
			break;
			}

		/* Extract selected fields if not probable duplicate */
		    else if (star.magetc > 0) {
			ra = ujcra (star.rasec); /* Right ascension in degrees */
			dec = ujcdec (star.decsec); /* Declination in degrees */
			mag = ujcmag (star.magetc);	/* Magnitude */
			plate = ujcplate (star.magetc);	/* Plate number */
			if (drad > 0)
			    dist = wcsdist (cra,cdec,ra,dec);
			else
			    dist = 0.0;

		    /* Check magnitude and position limits */
			if ((mag1 == mag2 || (mag >= mag1 && mag <= mag2)) &&
     			    (dec >= dec1 && dec <= dec2) &&
			    ((wrap && (ra <= ra1 || ra >= ra2)) ||
			    (!wrap && (ra >= ra1 && ra <= ra2))) &&
			    (drad == 0.0 || dist < drad) &&
			    (xplate == 0 || plate == xplate)) {

			/* Save star position and magnitude in table */
			    if (nstar <= nstarmax) {
				unum[nstar] = (double) znum +
					      (0.0000001 * (double)istar);
				ura[nstar] = ra;
				udec[nstar] = dec;
				umag[nstar] = mag;
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

			/* If too many stars and radial search,
			   replace furthest star */
			    else if (drad > 0 && dist < maxdist) {
				ura[farstar] = ra;
				udec[farstar] = dec;
				umag[farstar] = mag;
				uplate[farstar] = plate;
				udist[farstar] = dist;
				maxdist = 0.0;

			    /* Find new farthest star */
				for (i = 0; i < nstarmax; i++) {
				    if (udist[i] > maxdist) {
					maxdist = udist[i];
					farstar = i;
					}
				    }
				}

			/* If too many stars, replace faintest star */
			    else if (mag < faintmag) {
				ura[faintstar] = ra;
				udec[faintstar] = dec;
				umag[faintstar] = mag;
				uplate[faintstar] = plate;
				udist[faintstar] = dist;
				faintmag = 0.0;

			    /* Find new faintest star */
				for (i = 0; i < nstarmax; i++) {
				    if (umag[i] > faintmag) {
					faintmag = umag[i];
					faintstar = i;
					}
				    }
				}
			    nstar++;
			    jstar++;
			    if (nlog == 1)
				fprintf (stderr,"UJCREAD: %04d.%04d: %9.5f %9.5f %5.2f\n",
				    znum,istar,ra,dec,mag);

			/* End of accepted star processing */
			    }

		    /* End of individual star processing */
			}

		/* Log operation */
		    if (nlog > 0 && itable%nlog == 0)
			fprintf (stderr,"UJCREAD: zone %d (%4d / %4d) %6d / %6d sources\r",
				znum, iz, nz, jstar, itable);

		/* End of star loop */
		    }

		/* End of wrap loop */
		}

	/* Close zone input file */
	    (void) fclose (fcat);
	    itot = itot + itable;
	    if (nlog > 0)
		fprintf (stderr,"UJCREAD: zone %d (%4d / %4d) %6d / %6d / %8d sources\n",
			znum, iz+1, nz, jstar, itable, nstars);

	/* End of zone processing */
	    }

    /* End of zone loop */
	}

/* Summarize search */
    if (nlog > 0) {
	if (nz > 1)
	    fprintf (stderr,"UJCREAD: %d zone: %d / %d found\n",nz,nstar,itot);
	else
	    fprintf (stderr,"UJCREAD: 1 zone: %d / %d found\n",nstar,itable);
	}
    if (nstar > nstarmax) {
	fprintf (stderr,"UJCREAD: %d stars found; only %d returned\n",
		 nstar,nstarmax);
	nstar = nstarmax;
	}
    free ((char *)udist);
    return (nstar);
}


int
ujcrnum (nnum,unum,ura,udec,umag,uplate,nlog)

int	nnum;		/* Number of stars to find */
double	*unum;		/* Array of UA numbers to find */
double	*ura;		/* Array of right ascensions (returned) */
double	*udec;		/* Array of declinations (returned) */
double	*umag;		/* Array of red magnitudes (returned) */
int	*uplate;	/* Array of plate numbers (returned) */
int	nlog;		/* Logging interval */
{
    UJCstar star;	/* UJ catalog entry for one star */

    int znum;
    int jnum;
    int nzone;
    int nfound = 0;
    double ra,dec;
    double mag;
    int istar, plate;
    char *str;

    /* Set path to USNO J Catalog */
    if ((str = getenv("UJ_PATH")) != NULL )
	strcpy (cdu,str);

/* Loop through star list */
    for (jnum = 0; jnum < nnum; jnum++) {

    /* Get path to zone catalog */
	znum = (int) unum[jnum];
	if ((nzone = ujcopen (znum)) != 0) {

	    istar = (int) (((unum[jnum] - znum) * 100000000.0) + 0.5);

	/* Check to make sure star can be in this zone */
	    if (istar > nzone) {
		fprintf (stderr,"UACRNUM: Star %d > zone max. %d\n",
			 istar, nzone);
		break;
		}

	/* Read star entry from catalog */
	    if (ujcstar (istar, &star)) {
		fprintf (stderr,"UACRNUM: Cannot read star %d\n", istar);
		break;
		}

	    /* Extract selected fields if not probable duplicate */
	    else if (star.magetc > 0) {
		ra = ujcra (star.rasec); /* Right ascension in degrees */
		dec = ujcdec (star.decsec); /* Declination in degrees */
		mag = ujcmag (star.magetc); /* Magnitude */
		plate = ujcplate (star.magetc);	/* Plate number */

		/* Save star position and magnitude in table */
		ura[nfound] = ra;
		udec[nfound] = dec;
		umag[nfound] = mag;
		uplate[nfound] = plate;

		nfound++;
		if (nlog == 1)
		    fprintf (stderr,"UJCRNUM: %04d.%08d: %9.5f %9.5f %5.2f\n",
			     znum,istar,ra,dec,mag);

		/* Log operation */
		if (nlog > 0 && jnum%nlog == 0)
		    fprintf (stderr,"UJCRNUM: %04d.%08d  %8d / %8d sources\r",
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
	fprintf (stderr,"UJCRNUM:  %d / %d found\n",nfound,nnum);

    return (nfound);
}


/* Declination zone numbers */
int zone[NZONES]={0,75,150,225,300,375,450,525,600,675,750,825,900,
	      975,1050,1125,1200,1275,1350,1425,1500,1575,1650,1725};

/* UJCZONES -- figure out which UJ zones will need to be searched */

static int
ujczones (ra1, ra2, dec1, dec2, nzmax, zones, verbose)

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
    iz1 = ujczone (dec1);
    iz2 = ujczone (dec2);

/* Tabulate zones to search */
    i = 0;
    if (iz2 >= iz1) {
	for (iz = iz1; iz <= iz2; iz++)
	    zones[i++] = zone[iz];
	}
    else {
	for (iz = iz2; iz <= iz1; iz++)
	    zones[i++] = zone[iz];
	}

    nrgn = i;
    if (verbose)
	printf ("UJCREG:  %d zones: %d - %d\n",nrgn,zones[0],zones[i-1]);
    if (verbose)
	printf("UJCREG: RA: %.5f - %.5f, Dec: %.5f - %.5f\n",ra1,ra2,dec1,dec2);

    return (nrgn);
}


/* UJCRA -- returns right ascension in degrees from the UJ star structure */

static double
ujcra (rasec)

int rasec;	/* RA in 100ths of arcseconds from UJ catalog entry */
{
    return ((double) (rasec) / 360000.0);
}


/* UJCDEC -- returns the declination in degrees from the UJ star structure */

static double
ujcdec (decsec)

int decsec;	/* Declination in 100ths of arcseconds from UJ catalog entry */
{
    return ((double) (decsec - 32400000) / 360000.0);
}


/* UJCMAG -- returns the magnitude from the UJ star structure */

static double
ujcmag (magetc)

int magetc;	/* Quality, plate, and magnitude from UJ catalog entry */
{
    if (magetc < 0)
	return ((double) (-magetc % 10000) * 0.01);
    else
	return ((double) (magetc % 10000) * 0.01);
}


/* UJCPLATE -- returns the plate number from the UJ star structure */

static int
ujcplate (magetc)

int magetc;	/* Quality, plate, and magnitude from UJ catalog entry */
{
    if (magetc < 0)
	return ((double) (-magetc % 10000000) / 10000);
    else
	return ((double) (magetc % 10000000) / 10000);
}


/* UJCZONE -- find the UJ zone number where a declination can be found */

static int
ujczone (dec)

double dec;	/* declination in degrees */
{
    double zonesize = 7.5;	/* number of degrees per declination zone */

    return ((int) ((dec + 90.0) / zonesize));
}


/* UJCSRA -- Find UJ star closest to specified right ascension */

static int
ujcsra (rax0)

double	rax0;		/* Right ascension for which to search */
{
    int istar, istar1, istar2, nrep;
    double rax, ra1, ra, rdiff, rdiff1, rdiff2, sdiff;
    UJCstar star;	/* UJ catalog entry for one star */
    char rastrx[16];
    int debug = 0;

    rax = rax0;
    ra2str (rastrx, 16, rax, 3);
    istar1 = 1;
    if (ujcstar (istar1, &star))
	return (0);
    ra1 = ujcra (star.rasec);
    istar = nstars;
    nrep = 0;
    while (istar != istar1 && nrep < 20) {
	if (ujcstar (istar, &star))
	    break;
	else {
	    ra = ujcra (star.rasec);
	    if (ra == ra1)
		break;
	    if (debug) {
		char rastr[16];
		ra2str (rastr, 16, ra, 3);
		fprintf (stderr,"UJCSRA %d %d: %s (%s)\n",
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
		fprintf (stderr,"UJCSRA: ra1=    %.5f ra=     %.5f rax=    %.5f\n",
			 ra1,ra,rax);
		fprintf (stderr,"UJCSRA: rdiff=  %.5f rdiff1= %.5f rdiff2= %.5f\n",
			 rdiff,rdiff1,rdiff2);
		fprintf (stderr,"UJCSRA: istar1= %d istar= %d istar1= %d\n",
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

/* UJCOPEN -- Open UJ Catalog zone catalog, returning number of entries */

static int
ujcopen (znum)

int znum;	/* UJ Catalog zone */
{
    char zonepath[64];	/* Pathname for input UJ zone file */
    UJCstar star;	/* UJ catalog entry for one star */
    struct stat statbuff;
    
/* Get path to zone catalog */
    if (ujcpath (znum, zonepath)) {
	fprintf (stderr, "UJCOPEN: Cannot find zone catalog for %d\n", znum);
	return (0);
	}

/* Find number of stars in zone catalog by its length */
    if (stat (zonepath, &statbuff)) {
	fprintf (stderr,"UJCOPEN: Zone catalog %s has no entries\n",zonepath);
	return (0);
	}
    else
	nstars = (int) statbuff.st_size / 12;

/* Open zone catalog */
    if (!(fcat = fopen (zonepath, "r"))) {
	fprintf (stderr,"UJCOPEN: Zone catalog %s cannot be read\n",zonepath);
	return (0);
	}

/* Check to see if byte-swapping is necessary */
    cswap = 0;
    if (ujcstar (1, &star)) {
	fprintf (stderr,"UJCOPEN: cannot read star 1 from UJ zone catalog %s\n",
		 zonepath);
	return (0);
	}
    else {
	if (star.rasec > 100000 || star.rasec < 0) {
	    cswap = 1;
	    fprintf (stderr,"UJCOPEN: swapping bytes in UJ zone catalog %s\n",
		     zonepath);
	    }
	else
	    cswap = 0;
	}

    return (nstars);
}


/* UJCPATH -- Get UJ Catalog region file pathname */

static int
ujcpath (zn, path)

int zn;		/* UJ zone number */
char *path;	/* Pathname of UJ zone file */

{
    if (zn < 0 || zn > 1725) {
	fprintf (stderr, "UJCPATH: zone %d out of range 0-1725\n",zn);
	return (-1);
	}
    else if (strchr (cdu,'C'))
	sprintf (path,"%s/ZONE%04d.CAT", cdu, zn);
    else
	sprintf (path,"%s/zone%04d.cat", cdu, zn);

    return (0);
}


/* UJCSTAR -- Get UJ catalog entry for one star; return 0 if successful */

static int
ujcstar (istar, star)

int istar;	/* Star sequence number in UJ zone catalog */
UJCstar *star;	/* UJ catalog entry for one star */
{
    int nbs, nbr, nbskip;

    if (istar < 1 || istar > nstars) {
	fprintf (stderr, "UJCstar %d is not in catalog\n",istar);
	return (-1);
	}
    nbskip = 12 * (istar - 1);
    if (fseek (fcat,nbskip,SEEK_SET))
	return (-1);
    nbs = sizeof (UJCstar);
    nbr = fread (star, nbs, 1, fcat) * nbs;
    if (nbr < nbs) {
	fprintf (stderr, "UJCstar %d / %d bytes read\n",nbr, nbs);
	return (-2);
	}
    if (cswap)
	ujcswap ((char *)star);
    return (0);
}


/* UJCSWAP -- Reverse bytes of UJ Catalog entry */

static void
ujcswap (string)

char *string;	/* Start of vector of 4-byte ints */

{
char *sbyte, *slast;
char temp0, temp1, temp2, temp3;
int nbytes = 12; /* Number of bytes to reverse */

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

/* May 20 1996	New subroutine
 * May 23 1996	Add optional plate number check
 * Aug  6 1996	Remove unused variables after lint
 * Oct 15 1996	Add comparison when testing an assignment
 * Oct 17 1996	Fix after lint: cast argument to UJCSWAP
 * Nov 13 1996	Return no more than maximum star number
 * Nov 13 1996	Write all error messages to stderr with subroutine names
 * Nov 15 1996  Implement search radius; change input arguments
 * Dec 12 1996	Improve logging
 * Dec 17 1996	Add code to keep brightest or closest stars if too many found
 * Dec 18 1996	Add code to read a specific star
 * Mar 20 1997	Clean up UJCRNUM after lint
 *
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Jun 24 1998	Initialize byte-swapping flag in UJCOPEN()
 */
