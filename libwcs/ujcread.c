/*** File libwcs/ujcread.c
 *** August 6, 1996
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "fitshead.h"

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
ujcread (ra1,ra2,dec1,dec2,mag1,mag2,xplate,nstarmax,unum,ura,udec,umag,
	 uplate,verbose)

double	ra1,ra2;	/* Limiting right ascensions of region in degrees */
double	dec1,dec2;	/* Limiting declinations of region in degrees */
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
    int nz;		/* Number of input UJ zone files */
    int zlist[NZONES];	/* List of input UJ zones */
    UJCstar star;	/* UJ catalog entry for one star */

    int wrap, iwrap;
    int znum, itot,iz;
    int nlog,itable,jstar;
    int nstar;
    double ra,dec;
    double mag;
    int istar, istar1, istar2, plate;
    int nzmax = NZONES;	/* Maximum number of declination zones */

    itot = 0;

/* If RA range includes zero, split it in two */
    wrap = 0;
    if (ra1 > ra2)
	wrap = 1;
    else
	wrap = 0;

/* dec1 is always the smallest declination */
    if (dec1 > dec2) {
	dec = dec1;
	dec1 = dec2;
	dec2 = dec;
	}

/* mag1 is always the smallest magnitude */
    if (mag2 < mag1) {
	mag = mag2;
	mag2 = mag1;
	mag1 = mag;
	}

/* Find UJ Star Catalog regions in which to search */
    nz = ujczones (ra1, ra2, dec1, dec2, nzmax, zlist, verbose);
    if (nz <= 0) {
	printf ("UJREAD:  no UJ zones found\n");
	return (0);
	}

/* Logging interval */
    nlog = 0;
    nstar = 0;

/* Loop through region list */
    for (iz = 0; iz < nz; iz++) {

    /* Get path to zone catalog */
	znum = zlist[iz];
	if ((nstars = ujcopen (znum))) {

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

		jstar = 0;
		itable = 0;

	    /* Loop through zone catalog for this region */
		for (istar = istar1; istar <= istar2; istar++) {
		    itable ++;

		    if (ujcstar (istar, &star)) {
			printf ("UJREAD: Cannot read star %d\n", istar);
			break;
			}

		/* Extract selected fields if not probable duplicate */
		    else if (star.magetc > 0) {
			ra = ujcra (star.rasec); /* Right ascension in degrees */
			dec = ujcdec (star.decsec); /* Declination in degrees */
			mag = ujcmag (star.magetc);	/* Magnitude */
			plate = ujcplate (star.magetc);	/* Plate number */

		    /* Check magnitude amd position limits */
			if ((mag1 == mag2 || (mag >= mag1 && mag <= mag2)) &&
			    ((wrap && (ra <= ra1 || ra >= ra2)) ||
			    (!wrap && (ra >= ra1 && ra <= ra2))) &&
			    (xplate == 0 || plate == xplate) &&
     			    (dec >= dec1 && dec <= dec2)) {

			/* Save star position and magnitude in table */
			    if (nstar <= nstarmax) {
				unum[nstar] = (double) znum +
					      (0.0000001 * (double)istar);
				ura[nstar] = ra;
				udec[nstar] = dec;
				umag[nstar] = mag;
				uplate[nstar] = plate;
				}
			    nstar = nstar + 1;
			    jstar = jstar + 1;
			    if (nlog == 1)
				printf ("%04d.%04d: %9.5f %9.5f %5.2f\n",
				    znum,istar,ra,dec,mag);

			/* End of accepted star processing */
			    }

		    /* End of individual star processing */
			}

		/* Log operation */
		    if (nlog > 0 && itable%nlog == 0)
			printf ("%4d / %4d: %6d / %6d sources zone %d\r",
				iz,nz,jstar,itable,znum);

		/* End of star loop */
		    }

		/* End of wrap loop */
		}

	/* Close zone input file */
	    (void) fclose (fcat);
	    itot = itot + itable;
	    if (nlog > 0)
		printf ("%4d / %4d: %6d / %6d sources zone %d\n",
			iz+1, nz, jstar, itable, znum);

	/* End of zone processing */
	    }

    /* End of zone loop */
	}

/* Summarize search */
    if (nlog > 0) {
	if (nz > 1)
	    printf ("%d zone: %d / %d found\n",nz,nstar,itot);
	else
	    printf ("1 zone: %d / %d found\n",nstar,itable);
	}
    return (nstar);
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
    ra2str (rastrx, rax, 3);
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
		ra2str (rastr,ra,3);
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
	fprintf (stderr,"UJ zone catalog %s has no entries\n",zonepath);
	return (0);
	}
    else
	nstars = (int) statbuff.st_size / 12;

/* Open zone catalog */
    if (!(fcat = fopen (zonepath, "r"))) {
	fprintf (stderr,"UJ zone catalog %s cannot be read\n",zonepath);
	return (0);
	}

/* Check to see if byte-swapping is necessary */
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
	ujcswap (star);
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
 */
