/* File libwcs/catutil.c
 * September 21, 1999
 * By Doug Mink
 */

/* int RefCat (refcatname,title,syscat,eqcat,epcat)
 *	Return catalog type code, title, coord. system
 * void CatNum (refcat, nndec, dnum, numstr)
 *	Return formatted source number
 * void SearchLim (cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, verbose)
 *	Comput limiting RA and Dec from center and half-widths
 * int CatNumLen (refcat, nndec)
 *	Return length of source numbers
 * void SearchLim (cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, verbose)
 *	Compute limiting RA and Dec from center and half-widths
 * void RefLim (cra,cdec,dra,ddec,sysc,sysr,eqc,eqr,epc,ramin,ramax,decmin,decmax,verbose)
 *	Compute limiting RA and Dec in new system from center and half-widths
 * struct Range *RangeInit (string, ndef)
 *	Return structure containing ranges of numbers
 * int isrange (string)
 *	Return 1 if string is a range, else 0
 * int rstart (range)
 *	Restart at beginning of range
 * int rgetn (range)
 *	Return number of values from range structure
 * int rgeti4 (range)
 *	Return next number from range structure as 4-byte integer
 * int rgetr8 (range)
 *	Return next number from range structure as 8-byte floating point number
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "wcs.h"
#include "wcscat.h"

/* Return code for reference catalog or its type */

int
RefCat (refcatname, title, syscat, eqcat, epcat)

char	*refcatname;	/* Name of reference catalog */
char	*title;		/* Description of catalog (returned) */
int	*syscat;	/* Catalog coordinate system (returned) */
double	*eqcat;		/* Equinox of catalog (returned) */
double	*epcat;		/* Epoch of catalog (returned) */
{
    struct StarCat *starcat;

    if (strncmp(refcatname,"gs",2)==0 ||
	strncmp (refcatname,"GS",2)== 0) {
	strcpy (title, "HST Guide Stars");
	*syscat = WCS_J2000;
	*eqcat = 2000.0;
	*epcat = 2000.0;
	return (GSC);
	}
    else if (strncmp(refcatname,"us",2)==0 ||
	strncmp(refcatname,"US",2)==0) {
	*syscat = WCS_J2000;
	*eqcat = 2000.0;
	*epcat = 2000.0;
	if (index (refcatname, '1') != NULL) {
	    strcpy (title, "USNO SA-1.0 Catalog Stars");
	    return (USA1);
	    }
	else if (index (refcatname, '2') != NULL) {
	    strcpy (title, "USNO SA-2.0 Catalog Stars");
	    return (USA2);
	    }
	else {
	    strcpy (title, "USNO SA Catalog Stars");
	    return (USAC);
	    }
	}
    else if (strncmp(refcatname,"ua",2)==0 ||
	strncmp(refcatname,"UA",2)==0) {
	*syscat = WCS_J2000;
	*eqcat = 2000.0;
	*epcat = 2000.0;
	if (index (refcatname, '1') != NULL) {
	    strcpy (title, "USNO A-1.0 Catalog Stars");
	    return (UA1);
	    }
	else if (index (refcatname, '2') != NULL) {
	    strcpy (title, "USNO A-2.0 Catalog Stars");
	    return (UA2);
	    }
	else {
	    strcpy (title, "USNO A Catalog Stars");
	    return (UAC);
	    }
	}
    else if (strncmp(refcatname,"uj",2)==0 ||
	strncmp(refcatname,"UJ",2)==0) {
	strcpy (title, "USNO J Catalog Stars");
	*syscat = WCS_J2000;
	*eqcat = 2000.0;
	*epcat = 2000.0;
	return (UJC);
	}
    else if (strncmp(refcatname,"sao",3)==0 ||
	strncmp(refcatname,"SAO",3)==0) {
	strcpy (title, "SAO Catalog Stars");
	if ((starcat = binopen ("SAO")) == NULL)
	    return (0);
	*syscat = starcat->coorsys;
	*eqcat = starcat->equinox;
	*epcat = starcat->epoch;
	binclose (starcat);
	return (SAO);
	}
    else if (strncmp(refcatname,"ppm",3)==0 ||
	strncmp(refcatname,"PPM",3)==0) {
	strcpy (title, "PPM Catalog Stars");
	if ((starcat = binopen ("PPM")) == NULL)
	    return (0);
	*syscat = starcat->coorsys;
	*eqcat = starcat->equinox;
	*epcat = starcat->epoch;
	binclose (starcat);
	return (PPM);
	}
    else if (strncmp(refcatname,"iras",4)==0 ||
	strncmp(refcatname,"IRAS",4)==0) {
	strcpy (title, "IRAS Point Sources");
	if ((starcat = binopen ("IRAS")) == NULL)
	    return (0);
	*syscat = starcat->coorsys;
	*eqcat = starcat->equinox;
	*epcat = starcat->epoch;
	binclose (starcat);
	return (IRAS);
	}
    else if (strncmp(refcatname,"tyc",3)==0 ||
	strncmp(refcatname,"TYC",3)==0) {
	strcpy (title, "Tycho Catalog Stars");
	if ((starcat = binopen ("tycho")) == NULL)
	    return (0);
	*syscat = starcat->coorsys;
	*eqcat = starcat->equinox;
	*epcat = starcat->epoch;
	binclose (starcat);
	return (TYCHO);
	}
    else if (strncmp(refcatname,"hip",3)==0 ||
	strncmp(refcatname,"HIP",3)==0) {
	strcpy (title, "Hipparcos Catalog Stars");
	if ((starcat = binopen ("hipparcos")) == NULL)
	    return (0);
	*syscat = starcat->coorsys;
	*eqcat = starcat->equinox;
	*epcat = starcat->epoch;
	binclose (starcat);
	return (HIP);
	}
    else if (strncmp(refcatname,"act",3)==0 ||
	strncmp(refcatname,"ACT",3)==0) {
	strcpy (title, "ACT Catalog Stars");
	*syscat = WCS_J2000;
	*eqcat = 2000.0;
	*epcat = 2000.0;
	return (ACT);
	}
    else if (strncmp(refcatname,"bsc",3)==0 ||
	strncmp(refcatname,"BSC",3)==0) {
	strcpy (title, "Bright Star Catalog Stars");
	if ((starcat = binopen ("BSC5")) == NULL)
	    return (0);
	*syscat = starcat->coorsys;
	*eqcat = starcat->equinox;
	*epcat = starcat->epoch;
	binclose (starcat);
	return (BSC);
	}
    else if (isbin (refcatname)) {
	strcpy (title, refcatname);
	strcat (title, " Catalog Sources");
	if ((starcat = binopen (refcatname)) == NULL)
	    return (0);
	*syscat = starcat->coorsys;
	*eqcat = starcat->equinox;
	*epcat = starcat->epoch;
	binclose (starcat);
	return (BINCAT);
	}
    else if (istab (refcatname)) {
	strcpy (title, refcatname);
	strcat (title, " Catalog Sources");
	*syscat = WCS_J2000;
	*eqcat = 2000.0;
	*epcat = 2000.0;
	return (TABCAT);
	}
    else {
	strcpy (title, refcatname);
	strcat (title, " Catalog Sources");
	if ((starcat = catopen (refcatname)) == NULL)
	    return (0);
	*syscat = starcat->coorsys;
	*eqcat = starcat->equinox;
	*epcat = starcat->epoch;
	catclose (starcat);
	return (TXTCAT);
	}
}


void
CatNum (refcat, nndec, dnum, numstr)

int	refcat;		/* Catalog code */
int	nndec;		/* Number of decimal places ( >= 0) */
double	dnum;		/* Catalog number of source */
char	*numstr;	/* Formatted number (returned) */

{
    char nform[16];	/* Format for star number */

    /* USNO A1.0, A2.0, SA1.0, or SA2.0 Catalogs */
    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
	refcat == UAC  || refcat == UA1  || refcat == UA2)
	sprintf (numstr, "%13.8f", dnum);

    /* USNO UJ 1.0 Catalog */
    else if (refcat == UJC)
	sprintf (numstr, "%12.7f", dnum);

    /* HST Guide Star Catalog */
    else if (refcat == GSC)
	sprintf (numstr, "%9.4f", dnum);

    /* SAO, PPM, or IRAS Point Source Catalogs (TDC binary format) */
    else if (refcat==SAO || refcat==PPM || refcat==IRAS || refcat==BSC)
	sprintf (numstr, "%6d", (int)(dnum+0.5));

    /* Tycho, Hipparcos, or ACT catalogs */
    else if (refcat == TYCHO || refcat == HIP || refcat == ACT)
	sprintf (numstr, "%10.5f", dnum);

    /* Starbase tab-separated catalogs */
    else if (refcat == TABCAT)
	sprintf (numstr, "%9.4f", dnum);

    /* TDC binary or ASCII catalogs */
    else if (nndec > 0) {
	sprintf (nform,"%%%d.%df", nndec+5, nndec);
	sprintf (numstr, nform, dnum);
	}
    else
	sprintf (numstr, "%6d", (int)(dnum+0.5));

    return;
}


int
CatNumLen (refcat, nndec)

int	refcat;		/* Catalog code */
int	nndec;		/* Number of decimal places ( >= 0) */

{

    /* USNO A1.0, A2.0, SA1.0, or SA2.0 Catalogs */
    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
	refcat == UAC  || refcat == UA1  || refcat == UA2)
	return (13);

    /* USNO UJ 1.0 Catalog */
    else if (refcat == UJC)
	return (12);

    /* HST Guide Star Catalog */
    else if (refcat == GSC)
	return (9);

    /* SAO, PPM, or IRAS Point Source Catalogs (TDC binary format) */
    else if (refcat==SAO || refcat==PPM || refcat==IRAS || refcat==BSC)
	return (6);

    /* Tycho, Hipparcos, or ACT catalogs */
    else if (refcat == TYCHO || refcat == HIP || refcat == ACT)
	return (10);

    /* Starbase tab-separated catalogs */
    else if (refcat == TABCAT)
	return (9);

    /* TDC binary or ASCII catalogs */
    else if (nndec > 0) {
	return (nndec + 5);
	}
    else
	return (6);
}


void
SearchLim (cra, cdec, dra, ddec, ra1, ra2, dec1, dec2, verbose)

double	cra, cdec;	/* Center of search area  in degrees */
double	dra, ddec;	/* Horizontal and verticla half-widths of area */
double	*ra1, *ra2;	/* Right ascension limits in degrees */
double	*dec1, *dec2;	/* Declination limits in degrees */
int	verbose;	/* 1 to print limits, else 0 */

{
    double dec;

    /* Set right ascension limits for search */
    *ra1 = cra - dra;
    *ra2 = cra + dra;

    /* Keep right ascension between 0 and 360 degrees */
    if (*ra1 < 0.0)
	*ra1 = *ra1 + 360.0;
    if (*ra2 > 360.0)
	*ra2 = *ra2 - 360.0;

    /* Set declination limits for search */
    *dec1 = cdec - ddec;
    *dec2 = cdec + ddec;

    /* dec1 is always the smallest declination */
    if (*dec1 > *dec2) {
	dec = *dec1;
	*dec1 = *dec2;
	*dec2 = dec;
	}

    /* Search zones which include the poles cover 360 degrees in RA */
    if (*dec1 < -90.0) {
	*dec1 = -90.0;
	*ra1 = 0.0;
	*ra2 = 359.99999;
	}
    if (*dec2 > 90.0) {
	*dec2 = 90.0;
	*ra1 = 0.0;
	*ra2 = 359.99999;
	}
    if (verbose) {
	char rstr1[16],rstr2[16],dstr1[16],dstr2[16];
	ra2str (rstr1, 16, *ra1, 3);
        dec2str (dstr1, 16, *dec1, 2);
	ra2str (rstr2, 16, *ra2, 3);
        dec2str (dstr2, 16, *dec2, 2);
	fprintf (stderr,"SearchLim: RA: %s - %s  Dec: %s - %s\n",
		 rstr1,rstr2,dstr1,dstr2);
	}
    return;
}


void
RefLim (cra, cdec, dra, ddec, sysc, sysr, eqc, eqr, epc,
	ramin, ramax, decmin, decmax, verbose)

double	cra, cdec;	/* Center of search area  in degrees */
double	dra, ddec;	/* Horizontal and verticla half-widths of area */
int	sysc, sysr;	/* System of search, catalog coordinates */
double	eqc, eqr;	/* Equinox of search, catalog coordinates in years */
double	epc;		/* Epoch of search coordinates in years
			   (catalog is assumed to be at epoch eqr) */
double	*ramin,*ramax;	/* Right ascension search limits in degrees (returned)*/
double	*decmin,*decmax; /* Declination search limits in degrees (returned) */
int	verbose;	/* 1 to print limits, else 0 */

{
    double ra1, ra2, ra3, ra4, dec1, dec2, dec3, dec4;
    double dec;

    /* Set right ascension limits for search */
    ra1 = cra - dra;
    ra2 = cra + dra;

    /* Keep right ascension between 0 and 360 degrees */
    if (ra1 < 0.0)
	ra1 = ra1 + 360.0;
    if (ra2 > 360.0)
	ra2 = ra2 - 360.0;
    ra3 = ra1;
    ra4 = ra2;

    /* Set declination limits for search */
    dec1 = cdec - ddec;
    dec2 = cdec + ddec;

    /* dec1 is always the smallest declination */
    if (dec1 > dec2) {
	dec = dec1;
	dec1 = dec2;
	dec2 = dec;
	}
    dec3 = dec2;
    dec4 = dec1;

    /* Convert search corners to catalog coordinate system and equinox */
    wcscon (sysc, sysr, eqc, eqr, &ra1, &dec1, epc);
    wcscon (sysc, sysr, eqc, eqr, &ra2, &dec2, epc);
    wcscon (sysc, sysr, eqc, eqr, &ra3, &dec3, epc);
    wcscon (sysc, sysr, eqc, eqr, &ra4, &dec4, epc);

    *ramin = ra1;
    if (ra3 < *ramin)
	*ramin = ra3;
    *ramax = ra2;
    if (ra4 > *ramax)
	*ramax = ra4;
    *decmin = dec1;
    if (dec4 < *decmin)
	*decmin = dec4;
    *decmax = dec2;
    if (dec3 > *decmax)
	*decmax = dec3;

    /* Search zones which include the poles cover 360 degrees in RA */
    if (*decmin < -90.0) {
	*decmin = -90.0;
	*ramin = 0.0;
	*ramax = 359.99999;
	}
    if (*decmax > 90.0) {
	*decmax = 90.0;
	*ramin = 0.0;
	*ramax = 359.99999;
	}
    if (verbose) {
	char rstr1[16],rstr2[16],dstr1[16],dstr2[16];
	ra2str (rstr1, 16, *ramin, 3);
        dec2str (dstr1, 16, *decmin, 2);
	ra2str (rstr2, 16, *ramax, 3);
        dec2str (dstr2, 16, *decmax, 2);
	fprintf (stderr,"RefLim: RA: %s - %s  Dec: %s - %s\n",
		 rstr1,rstr2,dstr1,dstr2);
	}
    return;
}

struct Range *
RangeInit (string, ndef)

char	*string;	/* String containing numbers separated by , and - */
int	ndef;		/* Maximum allowable range value */

{
    struct Range *range;
    int ip, nvalues, irange;
    char *slast;
    double first, last, step;

    if (!isrange (string) && !isnum (string))
	return (NULL);
    ip = 0;
    range = (struct Range *)calloc (1, sizeof (struct Range));
    range->irange = -1;
    range->nvalues = 0;
    range->nranges = 0;

    for (irange = 0; irange < MAXRANGE; irange++) {

	/* Default to entire list */
	first = 1.0;
	last = ndef;
	step = 1.0;

	/* Skip delimiters to start of range */
	while (string[ip] == ' ' || string[ip] == '	' ||
	       string[ip] == ',')
	    ip++;

	/* Get first limit
	 * Must be a number, '-', 'x', or EOS.  If not return ERR */
	if (string[ip] == (char)0) {	/* end of list */
	    if (irange == 0) {

		/* Null string defaults */
		range->ranges[0] = first;
		if (first < 1)
		    range->ranges[1] = first;
		else
		    range->ranges[1] = last;
		range->ranges[2] = step;
		range->nvalues = range->nvalues + 1 +
			  ((range->ranges[1]-range->ranges[0])/step);
		range->nranges++;
		return (range);
		}
	    else
		return (range);
	    }
	else if (string[ip] == '-' || string[ip] == ':')
	    ;
	else if (string[ip] == 'x')
	    ;
	else if (string[ip] > (char)47 && string[ip] < 58) {
	    first = strtod (string+ip, &slast);
	    ip = slast - string;
	    }
	else {
	    free (range);
	    return (NULL);
	    }

	/* Skip delimiters */
	while (string[ip] == ' ' || string[ip] == '	' ||
	       string[ip] == ',')
	    ip++;

	/* Get last limit
	* Must be '-', or 'x' otherwise last = first */
	if (string[ip] == 'x')
	    ;
	else if (string[ip] == '-' || string[ip] == ':') {
	    ip++;
	    while (string[ip] == ' ' || string[ip] == '	' ||
	   	   string[ip] == ',')
		ip++;
	    if (string[ip] == (char)0)
		last = first + ndef;
	    else if (string[ip] > (char)47 && string[ip] < 58) {
		last = strtod (string+ip, &slast);
		ip = slast - string;
		}
	    else if (string[ip] == 'x')
		;
	    else
		last = first + ndef;
	    }
	else
	    last = first;

	/* Skip delimiters */
	while (string[ip] == ' ' || string[ip] == '	' ||
	       string[ip] == ',')
	    ip++;

	/* Get step
	 * Must be 'x' or assume default step. */
	if (string[ip] == 'x') {
	    ip++;
	    while (string[ip] == ' ' || string[ip] == '	' ||
	   	   string[ip] == ',')
		ip++;
	    if (string[ip] == (char)0)
		step = 1.0;
	    else if (string[ip] > (char)47 && string[ip] < 58) {
		step = strtod (string+ip, &slast);
		ip = slast - string;
		}
	    else if (string[ip] == '-' || string[ip] == ':')
		;
	    else
		step = 1.0;
            }

	/* Output the range triple */
	range->ranges[irange*3] = first;
	range->ranges[irange*3 + 1] = last;
	range->ranges[irange*3 + 2] = step;
	range->nvalues = range->nvalues + ((last-first+(0.1*step)) / step + 1);
	range->nranges++;
	}

    return (range);
}


/* Return 1 if string is a range, else 0 */

int
isrange (string)

char *string;		/* String which might be a range of numbers */

{
    int i, lstr;

    /* If range separators present, check to make sure string is range */
    if (strchr (string+1, '-') || strchr (string+1, ',')) {
	lstr = strlen (string);
	for (i = 0; i < lstr; i++) {
	    if (strchr ("0123456789-,.x", (int)string[i]) == NULL)
		return (0);
	    }
	return (1);
	}
    else
	return (0);
}


/* Restart at beginning of range */

void
rstart (range)

struct Range *range;	/* Range structure */

{
    range->irange = -1;
    return;
}


/*  Return number of values from range structure */

int
rgetn (range)

struct Range *range;	/* Range structure */

{
    return (range->nvalues);
}


/*  Return next number from range structure as 8-byte floating point number */

double
rgetr8 (range)

struct Range *range;	/* Range structure */

{
    int i;

    if (range == NULL)
	return (0.0);
    else if (range->irange < 0) {
	range->irange = 0;
	range->first = range->ranges[0];
	range->last = range->ranges[1];
	range->step = range->ranges[2];
	range->value = range->first;
	}
    else {
	range->value = range->value + range->step;
	if (range->value > range->last) {
	    range->irange++;
	    if (range->irange < range->nranges) {
		i = range->irange * 3;
		range->first = range->ranges[i];
		range->last = range->ranges[i+1];
		range->step = range->ranges[i+2];
		range->value = range->first;
		}
	    else
		range->value = 0.0;
	    }
	}
    return (range->value);
}


/*  Return next number from range structure as 4-byte integer */

int
rgeti4 (range)

struct Range *range;	/* Range structure */

{
    double value;

    value = rgetr8 (range);
    return ((int) (value + 0.000000001));
}

/* Mar  2 1998	Make number and second magnitude optional
 * Oct 21 1998	Add RefCat() to set reference catalog code
 * Oct 26 1998	Include object names in star catalog entry structure
 * Oct 29 1998	Return coordinate system and title from RefCat
 * Nov 20 1998	Add USNO A-2.0 catalog and return different code
 * Dec  9 1998	Add Hipparcos and Tycho catalogs
 *
 * Jan 26 1999	Add subroutines to deal with ranges of numbers
 * Feb  8 1999	Fix bug initializing ACT catalog
 * Feb 11 1999	Change starcat.insys to starcat.coorsys
 * May 19 1999	Separate catalog subroutines into separate file
 * May 19 1999	Add CatNum() to return properly formatted catalog number
 * May 20 1999	Add date/time conversion subroutines translated from Fortran
 * May 28 1999	Fix bug in CatNum() which omitted GSC
 * Jun  3 1999	Add return to CatNum()
 * Jun  3 1999	Add CatNumLen()
 * Jun 16 1999	Add SearchLim(), used by all catalog search subroutines
 * Jun 30 1999	Add isrange() to check to see whether a string is a range
 * Jul  1 1999	Move date and time utilities to dateutil.c
 * Jul 23 1999	Add Bright Star Catalog
 * Aug 16 1999	Add RefLim() to set catalog search limits
 * Sep 21 1999	In isrange(), check for x
 */
