/*** File libwcs/uacread.c
 *** December 16, 1996
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

static char cdu[64]="/data/ua"; /* pathname of USNO A 1.0 directory */
				/* Replaced by UA_PATH environment variable */

typedef struct {
    int rasec, decsec, magetc;
} UACstar;

static int nstars;	/* Number of stars in catalog */
static int cswap = 0;	/* Byte reverse catalog to Intel/DEC order if 1 */
static FILE *fcat;
#define ABS(a) ((a) < 0 ? (-(a)) : (a))
#define NZONES 24

static double uacra();
static double uacdec();
static double uacmagr();
static double uacmagb();
static int uacmagerr();
static int uacgsc();
static int uacplate();

static int uaczones();
static int uaczone();
static int uacsra();
static int uacopen();
static int uacpath();
static int uacstar();
static void uacswap();

/* UACREAD -- Read USNO J Catalog stars from CDROM */

int
uacread (cra,cdec,dra,ddec,drad,mag1,mag2,xplate,nstarmax,
	 unum,ura,udec,umag,umagb,uplate,nlog)

double	cra;		/* Search center J2000 right ascension in degrees */
double	cdec;		/* Search center J2000 declination in degrees */
double	dra;		/* Search half width in right ascension in degrees */
double	ddec;		/* Search half-width in declination in degrees */
double	drad;		/* Limiting separation in degrees (ignore if 0) */
double	mag1,mag2;	/* Limiting magnitudes (none if equal) */
int	xplate;		/* If nonzero, use objects only from this plate */
int	nstarmax;	/* Maximum number of stars to be returned */
double	*unum;		/* Array of UA numbers (returned) */
double	*ura;		/* Array of right ascensions (returned) */
double	*udec;		/* Array of declinations (returned) */
double	*umag;		/* Array of red magnitudes (returned) */
double	*umagb;		/* Array of blue magnitudes (returned) */
int	*uplate;	/* Array of plate numbers (returned) */
int	nlog;		/* Logging interval */
{
    double ra1,ra2;	/* Limiting right ascensions of region in degrees */
    double dec1,dec2;	/* Limiting declinations of region in degrees */
    int nz;		/* Number of input UA zone files */
    int zlist[NZONES];	/* List of input UA zones */
    UACstar star;	/* UA catalog entry for one star */
    double faintmag=0.0; /* Faintest magnitude */
    double dist = 0.0;	/* Distance from search center in degrees */
    int	faintstar=0;	/* Faintest star */

    int wrap, iwrap;
    int verbose;
    int znum, itot,iz, i;
    int itable, jtable,jstar;
    int nstar, nread;
    double ra,dec;
    double mag, magb;
    int istar, istar1, istar2, plate;
    int nzmax = NZONES;	/* Maximum number of declination zones */
    char *str;

    itot = 0;
    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    /* Set path to USNO A Catalog */
    if ((str = getenv("UA_PATH")) != NULL )
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
	ra2str (rstr1, ra1, 3);
        dec2str (dstr1, dec1, 2);
	ra2str (rstr2, ra2, 3);
        dec2str (dstr2, dec2, 2);
	fprintf (stderr,"UACREAD: RA: %s - %s  Dec: %s - %s\n",
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

/* Find UA Star Catalog regions in which to search */
    nz = uaczones (ra1, ra2, dec1, dec2, nzmax, zlist, verbose);
    if (nz <= 0) {
	fprintf (stderr, "UACREAD:  no UJ zones found\n");
	return (0);
	}

/* Loop through region list */
    nstar = 0;
    for (iz = 0; iz < nz; iz++) {

    /* Get path to zone catalog */
	znum = zlist[iz];
	if ((nstars = uacopen (znum)) != 0) {

	    jstar = 0;
	    jtable = 0;
	    for (iwrap = 0; iwrap <= wrap; iwrap++) {

	    /* Find first star based on RA */
		if (iwrap == 0 || wrap == 0)
		    istar1 = uacsra (ra1);
		else
		    istar1 = 1;

	    /* Find last star based on RA */
		if (iwrap == 1 || wrap == 0)
		    istar2 = uacsra (ra2);
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

		    if (uacstar (istar, &star)) {
			fprintf (stderr,"UACREAD: Cannot read star %d\n", istar);
			break;
			}

		/* Extract selected fields if not probable duplicate */
		    else if (star.magetc > 0) {
			ra = uacra (star.rasec); /* Right ascension in degrees */
			dec = uacdec (star.decsec); /* Declination in degrees */
			magb = uacmagb (star.magetc); /* Blue magnitude */
			mag = uacmagr (star.magetc); /* Red magnitude */
			plate = uacplate (star.magetc);	/* Plate number */
			if (drad > 0)
			    dist = wcsdist (cra,cdec,ra,dec);

		    /* Check magnitude amd position limits and plate number */
			if ((mag1 == mag2 || (mag >= mag1 && mag <= mag2)) &&
     			    (dec >= dec1 && dec <= dec2) &&
			    ((wrap && (ra <= ra1 || ra >= ra2)) ||
			    (!wrap && (ra >= ra1 && ra <= ra2))) &&
			    (drad == 0.0 || dist < drad) &&
			    (xplate == 0 || plate == xplate)) {

			/* Save star position and magnitude in table */
			    if (nstar <= nstarmax) {
				unum[nstar] = (double) znum +
					      (0.00000001 * (double)istar);
				ura[nstar] = ra;
				udec[nstar] = dec;
				umag[nstar] = mag;
				umagb[nstar] = magb;
				uplate[nstar] = plate;
				if (mag > faintmag) {
				    faintmag = mag;
				    faintstar = nstar;
				    }
				}

			/* If too many stars, replace faintest star */
			    else if (mag < faintmag) {
				ura[faintstar] = ra;
				udec[faintstar] = dec;
				umag[faintstar] = mag;
				umagb[faintstar] = magb;
				uplate[faintstar] = plate;
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
				fprintf (stderr,"UACREAD: %04d.%08d: %9.5f %9.5f %5.2f\n",
				    znum,istar,ra,dec,mag);

			/* End of accepted star processing */
			    }

		    /* End of individual star processing */
			}

		/* Log operation */
		    if (nlog > 0 && itable%nlog == 0)
			fprintf (stderr,"UACREAD: zone %d (%2d / %2d) %8d / %8d / %8d sources\r",
				znum, iz+1, nz, jstar, itable, nread);

		/* End of star loop */
		    }

		/* End of wrap loop */
		}

	/* Close zone input file */
	    (void) fclose (fcat);
	    itot = itot + itable;
	    if (nlog > 0)
		fprintf (stderr,"UACREAD: zone %d (%2d / %2d) %8d / %8d / %8d sources      \n",
			znum, iz+1, nz, jstar, jtable, nstars);

	/* End of zone processing */
	    }

    /* End of zone loop */
	}

/* Summarize search */
    if (nlog > 0) {
	if (nz > 1)
	    fprintf (stderr,"UACREAD: %d zones: %d / %d found\n",nz,nstar,itot);
	else
	    fprintf (stderr,"UACREAD: 1 zone: %d / %d found\n",nstar,itable);
	}
    if (nstar > nstarmax) {
	fprintf (stderr,"UACREAD: %d stars found; only %d returned\n",
		 nstar,nstarmax);
	nstar = nstarmax;
	}
    return (nstar);
}


int
uacrnum (nstars,unum,ura,udec,umag,umagb,uplate,nlog)

int	nstars;		/* Number of stars to find */
double	*unum;		/* Array of UA numbers to find */
double	*ura;		/* Array of right ascensions (returned) */
double	*udec;		/* Array of declinations (returned) */
double	*umag;		/* Array of red magnitudes (returned) */
double	*umagb;		/* Array of blue magnitudes (returned) */
int	*uplate;	/* Array of plate numbers (returned) */
int	nlog;		/* Logging interval */
{
    UACstar star;	/* UA catalog entry for one star */

    int verbose;
    int znum, itot,iz, i;
    int itable, jtable,jstar;
    int nstar, nread, nzone;
    int nfound = 0;
    double ra,dec;
    double mag, magb;
    int istar, istar1, istar2, plate;
    int nzmax = NZONES;	/* Maximum number of declination zones */
    char *str;

    itot = 0;
    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    /* Set path to USNO A Catalog */
    if ((str = getenv("UA_PATH")) != NULL )
	strcpy (cdu,str);

/* Loop through star list */
    for (jstar = 0; jstar < nstars; jstar++) {

    /* Get path to zone catalog */
	znum = (int) unum[jstar];
	if ((nzone = uacopen (znum)) != 0) {

	    istar = (int) (((unum[jstar] - znum) * 100000000.0) + 0.5);

	    if (uacstar (istar, &star)) {
		fprintf (stderr,"UACRNUM: Cannot read star %d\n", istar);
		break;
		}

	    /* Extract selected fields if not probable duplicate */
	    else if (star.magetc > 0) {
		ra = uacra (star.rasec); /* Right ascension in degrees */
		dec = uacdec (star.decsec); /* Declination in degrees */
		magb = uacmagb (star.magetc); /* Blue magnitude */
		mag = uacmagr (star.magetc); /* Red magnitude */
		plate = uacplate (star.magetc);	/* Plate number */

		/* Save star position and magnitude in table */
		ura[nfound] = ra;
		udec[nfound] = dec;
		umag[nfound] = mag;
		umagb[nfound] = magb;
		uplate[nfound] = plate;

		nfound++;
		if (nlog == 1)
		    fprintf (stderr,"UACRNUM: %04d.%08d: %9.5f %9.5f %5.2f\n",
			     znum,istar,ra,dec,mag);

		/* Log operation */
		if (nlog > 0 && jstar%nlog == 0)
		    fprintf (stderr,"UACRNUM: %4d.%8d  %8d / %8d sources\r",
			     znum, istar, jstar, nstar);

		(void) fclose (fcat);
		/* End of star processing */
		}

	    /* End of star */
	    }

	/* End of star loop */
	}

    /* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"UACRNUM:  %d / %d found\n",nfound,nstars);

    return (nfound);
}


/* Declination zone numbers */
int azone[NZONES]={0,75,150,225,300,375,450,525,600,675,750,825,900,
	      975,1050,1125,1200,1275,1350,1425,1500,1575,1650,1725};

/* UACZONES -- figure out which UA zones will need to be searched */

static int
uaczones (ra1, ra2, dec1, dec2, nzmax, zones, verbose)

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
    iz1 = uaczone (dec1);
    iz2 = uaczone (dec2);

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
	fprintf(stderr,"UACZONES:  %d zones: %d - %d\n",nrgn,zones[0],zones[i-1]);
	fprintf(stderr,"UACZONES: RA: %.5f - %.5f, Dec: %.5f - %.5f\n",ra1,ra2,dec1,dec2);
	}

    return (nrgn);
}


/* UACRA -- returns right ascension in degrees from the UA star structure */

static double
uacra (rasec)

int rasec;	/* RA in 100ths of arcseconds from UA catalog entry */
{
    return ((double) (rasec) / 360000.0);
}


/* UACDEC -- returns the declination in degrees from the UA star structure */

static double
uacdec (decsec)

int decsec;	/* Declination in 100ths of arcseconds from UA catalog entry */
{
    return ((double) (decsec - 32400000) / 360000.0);
}


/* UACMAGR -- returns the red magnitude from the UA star structure */

static double
uacmagr (magetc)

int magetc;	/* Quality, plate, and magnitude from UA catalog entry */
{
    if (magetc < 0)
	return ((double) (-magetc % 1000) * 0.1);
    else
	return ((double) (magetc % 1000) * 0.1);
}


/* UACMAGB -- returns the blue magnitude from the UA star structure */

static double
uacmagb (magetc)

int magetc;	/* Quality, plate, and magnitude from UA catalog entry */
{
    if (magetc < 0)
	return ((double) ((-magetc / 1000) % 1000) * 0.1);
    else
	return ((double) ((magetc / 1000) % 1000) * 0.1);
}


/* UACMAGERR -- returns 1 if magnitude is uncertain from UA star structure */

static int
uacmagerr (magetc)

int magetc;	/* Quality, plate, and magnitude from UA catalog entry */
{
    if (magetc < 0)
	return ((-magetc / 1000000000) % 10);
    else
	return ((magetc / 1000000000) % 10);
}


/* UACGSC -- returns 1 if UA star is in the HST Guide Star Catalog */

static int
uacgsc (magetc)

int magetc;	/* Quality, plate, and magnitude from UA catalog entry */
{
    if (magetc < 0)
	return (1);
    else
	return (0);
}


/* UACPLATE -- returns the plate number from the UA star structure */

static int
uacplate (magetc)

int magetc;	/* Quality, plate, and magnitude from UA catalog entry */
{
    if (magetc < 0)
	return ((-magetc / 1000000) % 1000);
    else
	return ((magetc / 1000000) % 1000);
}


/* UACZONE -- find the UA zone number where a declination can be found */

static int
uaczone (dec)

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


/* UACSRA -- Find UA star closest to specified right ascension */

static int
uacsra (rax0)

double	rax0;		/* Right ascension for which to search */
{
    int istar, istar1, istar2, nrep;
    double rax, ra1, ra, rdiff, rdiff1, rdiff2, sdiff;
    UACstar star;	/* UA catalog entry for one star */
    char rastrx[16];
    int debug = 0;

    rax = rax0;
    ra2str (rastrx, rax, 3);
    istar1 = 1;
    if (uacstar (istar1, &star))
	return (0);
    ra1 = uacra (star.rasec);
    istar = nstars;
    nrep = 0;
    while (istar != istar1 && nrep < 20) {
	if (uacstar (istar, &star))
	    break;
	else {
	    ra = uacra (star.rasec);
	    if (ra == ra1)
		break;
	    if (debug) {
		char rastr[16];
		ra2str (rastr,ra,3);
		fprintf (stderr,"UACSRA %d %d: %s (%s)\n",
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

/* UACOPEN -- Open UA Catalog zone catalog, returning number of entries */

static int
uacopen (znum)

int znum;	/* UA Catalog zone */
{
    char zonepath[64];	/* Pathname for input UA zone file */
    UACstar star;	/* UA catalog entry for one star */
    struct stat statbuff;
    
/* Get path to zone catalog */
    if (uacpath (znum, zonepath)) {
	fprintf (stderr, "UACOPEN: Cannot find zone catalog for %d\n", znum);
	return (0);
	}

/* Find number of stars in zone catalog by its length */
    if (stat (zonepath, &statbuff)) {
	fprintf (stderr,"UA zone catalog %s has no entries\n",zonepath);
	return (0);
	}
    else
	nstars = (int) statbuff.st_size / 12;

/* Open zone catalog */
    if (!(fcat = fopen (zonepath, "r"))) {
	fprintf (stderr,"UA zone catalog %s cannot be read\n",zonepath);
	return (0);
	}

/* Check to see if byte-swapping is necessary */
    if (uacstar (1, &star)) {
	fprintf (stderr,"UACOPEN: cannot read star 1 from UA zone catalog %s\n",
		 zonepath);
	return (0);
	}
    else {
	if (star.rasec > 100000 || star.rasec < 0) {
	    cswap = 1;
	    fprintf (stderr,"UACOPEN: swapping bytes in UA zone catalog %s\n",
		     zonepath);
	    }
	else
	    cswap = 0;
	}

    return (nstars);
}


/* UACPATH -- Get UA Catalog region file pathname */

static int
uacpath (zn, path)

int zn;		/* UA zone number */
char *path;	/* Pathname of UA zone file */

{
    if (zn < 0 || zn > 1725) {
	fprintf (stderr, "UACPATH: zone %d out of range 0-1725\n",zn);
	return (-1);
	}
    else if (strchr (cdu,'C'))
	sprintf (path,"%s/ZONE%04d.CAT", cdu, zn);
    else
	sprintf (path,"%s/zone%04d.cat", cdu, zn);

    return (0);
}


/* UACSTAR -- Get UA catalog entry for one star; return 0 if successful */

static int
uacstar (istar, star)

int istar;	/* Star sequence number in UA zone catalog */
UACstar *star;	/* UA catalog entry for one star */
{
    int nbs, nbr, nbskip;

    if (istar < 1 || istar > nstars) {
	fprintf (stderr, "UACstar %d is not in catalog\n",istar);
	return (-1);
	}
    nbskip = 12 * (istar - 1);
    if (fseek (fcat,nbskip,SEEK_SET))
	return (-1);
    nbs = sizeof (UACstar);
    nbr = fread (star, nbs, 1, fcat) * nbs;
    if (nbr < nbs) {
	fprintf (stderr, "UACstar %d / %d bytes read\n",nbr, nbs);
	return (-2);
	}
    if (cswap)
	uacswap ((char *)star);
    return (0);
}


/* UACSWAP -- Reverse bytes of UA Catalog entry */

static void
uacswap (string)

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

/* Nov 15 1996	New subroutine
 * Dec 11 1996	Set ra<0 to ra+360 and ra>360 to ra-360
 * Dec 16 1996	Add code to read a specific star
 */
