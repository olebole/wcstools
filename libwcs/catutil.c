/* File libwcs/catutil.c
 * December 15, 2000
 * By Doug Mink, dmink@cfa.harvard.edu
 */

/* int RefCat (refcatname,title,syscat,eqcat,epcat)
 *	Return catalog type code, title, coord. system
 * char *ProgCat (progname)
 *	Return catalog name from program name, NULL if none there
 * char *ProgName (progpath0)
 *	Return program name from pathname by which program is invoked
 * void CatNum (refcat, nndec, dnum, numstr)
 *	Return formatted source number
 * void SearchLim (cra, cdec, dra, ddec, sys, ra1, ra2, dec1, dec2, verbose)
 *	Compute limiting RA and Dec from center and half-widths
 * int CatNumLen (refcat, nndec)
 *	Return length of source number
 * int CatNdec (refcat)
 *	Return number of decimal places in source number, if known
 * int StrNdec (string)
 *	Returns number of decimal places in a numeric string (-1=not number)
 * int CatMaxField (refcat, maxnum, nndec)
 *	Return Maximum field size for source numbers
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
 * int setoken (tokens, string, cwhite)
 *	Tokenize a string for easy decoding
 * int nextoken (tokens, token)
 *	Get next token from tokenized string
 * int getoken (tokens, itok, token)
 *	Get specified token from tokenized string
 * int isfile (filename)
 *	Return 1 if file is a readable file, else 0
 * int agets (string, keyword, lval, value)
 *	Read value from a file where keyword=value, anywhere on a line
 * void bv2sp (bv, b, v, isp)
 *	approximate spectral type given B - V or B and V magnitudes
 * void br2sp (br, b, r, isp)
 *	approximate spectral type given B - R or B and R magnitudes
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "wcs.h"
#include "wcscat.h"

static int catprop = 0;		/* 1 if catalog has proper motion, else 0 */

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
    int refcat;

    catprop = 0;

    if (strncmp(refcatname,"gs",2)==0 ||
	strncmp (refcatname,"GS",2)== 0) {
	strcpy (title, "HST Guide Stars");
	*syscat = WCS_J2000;
	*eqcat = 2000.0;
	*epcat = 2000.0;
	catprop = 0;
	refcat = GSC;
	}
    else if (strncmp(refcatname,"us",2)==0 ||
	strncmp(refcatname,"US",2)==0) {
	*syscat = WCS_J2000;
	*eqcat = 2000.0;
	*epcat = 2000.0;
	catprop = 0;
	if (strncmp (refcatname, ".usno", 5) == 0) {
	    sprintf (title, "USNO %s Stars", refcatname);
	    refcat = USNO;
	    }
	else if (strchr (refcatname, '1') != NULL) {
	    strcpy (title, "USNO SA-1.0 Catalog Stars");
	    refcat = USA1;
	    }
	else if (strchr (refcatname, '2') != NULL) {
	    strcpy (title, "USNO SA-2.0 Catalog Stars");
	    refcat = USA2;
	    }
	else {
	    strcpy (title, "USNO SA Catalog Stars");
	    refcat = USAC;
	    }
	}
    else if (strncmp(refcatname,"ua",2)==0 ||
	strncmp(refcatname,"UA",2)==0) {
	*syscat = WCS_J2000;
	*eqcat = 2000.0;
	*epcat = 2000.0;
	catprop = 0;
	if (strchr (refcatname, '1') != NULL) {
	    strcpy (title, "USNO A-1.0 Catalog Stars");
	    refcat = UA1;
	    }
	else if (strchr (refcatname, '2') != NULL) {
	    strcpy (title, "USNO A-2.0 Catalog Stars");
	    refcat = UA2;
	    }
	else {
	    strcpy (title, "USNO A Catalog Stars");
	    refcat = UAC;
	    }
	}
    else if (strncmp(refcatname,"uj",2)==0 ||
	strncmp(refcatname,"UJ",2)==0) {
	strcpy (title, "USNO J Catalog Stars");
	*syscat = WCS_J2000;
	*eqcat = 2000.0;
	*epcat = 2000.0;
	catprop = 0;
	refcat = UJC;
	}
    else if (strncmp(refcatname,"sao",3)==0 ||
	strncmp(refcatname,"SAO",3)==0) {
	strcpy (title, "SAO Catalog Stars");
	starcat = binopen ("SAO");
	if (starcat == NULL)
	    starcat = binopen ("SAOra");
	if (starcat) {
	    *syscat = starcat->coorsys;
	    *eqcat = starcat->equinox;
	    *epcat = starcat->epoch;
	    catprop = starcat->mprop;
	    binclose (starcat);
	    refcat = SAO;
	    }
	}
    else if (strncmp(refcatname,"ppm",3)==0 ||
	strncmp(refcatname,"PPM",3)==0) {
	strcpy (title, "PPM Catalog Stars");
	starcat = binopen ("PPM");
	if (starcat == NULL)
	    starcat = binopen ("SAOra");
	if (starcat) {
	    *syscat = starcat->coorsys;
	    *eqcat = starcat->equinox;
	    *epcat = starcat->epoch;
	    catprop = starcat->mprop;
	    binclose (starcat);
	    refcat = PPM;
	    }
	}
    else if (strncmp(refcatname,"iras",4)==0 ||
	strncmp(refcatname,"IRAS",4)==0) {
	strcpy (title, "IRAS Point Sources");
	if ((starcat = binopen ("IRAS"))) {
	    *syscat = starcat->coorsys;
	    *eqcat = starcat->equinox;
	    *epcat = starcat->epoch;
	    catprop = starcat->mprop;
	    binclose (starcat);
	    refcat = IRAS;
	    }
	}
    else if (strncmp(refcatname,"ty",2)==0 ||
	strncmp(refcatname,"TY",2)==0) {
	if (strsrch (refcatname, "2") != NULL) {
	    strcpy (title, "Tycho 2 Catalog Stars");
	    *syscat = WCS_J2000;
	    *eqcat = 2000.0;
	    *epcat = 2000.0;
	    catprop = 1;
	    refcat = TYCHO2;
	    }
	else {
	    strcpy (title, "Tycho Catalog Stars");
	    if ((starcat = binopen ("tycho"))) {
		*syscat = starcat->coorsys;
		*eqcat = starcat->equinox;
		*epcat = starcat->epoch;
		catprop = 1;
		binclose (starcat);
		refcat = TYCHO;
		}
	    }
	}
    else if (strncmp(refcatname,"hip",3)==0 ||
	strncmp(refcatname,"HIP",3)==0) {
	strcpy (title, "Hipparcos Catalog Stars");
	if ((starcat = binopen ("hipparcos"))) {
	    *syscat = starcat->coorsys;
	    *eqcat = starcat->equinox;
	    *epcat = starcat->epoch;
	    catprop = starcat->mprop;
	    binclose (starcat);
	    refcat = HIP;
	    }
	}
    else if (strncmp(refcatname,"act",3)==0 ||
	strncmp(refcatname,"ACT",3)==0) {
	strcpy (title, "ACT Catalog Stars");
	*syscat = WCS_J2000;
	*eqcat = 2000.0;
	*epcat = 2000.0;
	catprop = 1;
	refcat = ACT;
	}
    else if (strncmp(refcatname,"bsc",3)==0 ||
	strncmp(refcatname,"BSC",3)==0) {
	strcpy (title, "Bright Star Catalog Stars");
	if ((starcat = binopen ("BSC5"))) {
	    *syscat = starcat->coorsys;
	    *eqcat = starcat->equinox;
	    *epcat = starcat->epoch;
	    catprop = starcat->mprop;
	    binclose (starcat);
	    refcat = BSC;
	    }
	}
    else if (strsrch (refcatname, ".usno")) {
	*syscat = WCS_J2000;
	*eqcat = 2000.0;
	*epcat = 2000.0;
	catprop = 0;
	sprintf (title, "USNO %s Stars", refcatname);
	refcat = USNO;
	}
    else if (isbin (refcatname)) {
	strcpy (title, refcatname);
	strcat (title, " Catalog Sources");
	if ((starcat = binopen (refcatname))) {
	    *syscat = starcat->coorsys;
	    *eqcat = starcat->equinox;
	    *epcat = starcat->epoch;
	    catprop = starcat->mprop;
	    binclose (starcat);
	    refcat = BINCAT;
	    }
	}
    else if (istab (refcatname)) {
	strcpy (title, refcatname);
	strcat (title, " Catalog Sources");
	if ((starcat = ctgopen (refcatname, TABCAT))) {
	    *syscat = starcat->coorsys;
	    *eqcat = starcat->equinox;
	    *epcat = starcat->epoch;
	    catprop = starcat->mprop;
	    ctgclose (starcat);
	    refcat = TABCAT;
	    }
	}
    else {
	strcpy (title, refcatname);
	strcat (title, " Catalog Sources");
	if ((starcat = ctgopen (refcatname, TXTCAT))) {
	    *syscat = starcat->coorsys;
	    *eqcat = starcat->equinox;
	    *epcat = starcat->epoch;
	    catprop = starcat->mprop;
	    ctgclose (starcat);
	    refcat = TXTCAT;
	    }
	}
    return refcat;
}


char *
ProgName (progpath0)

char *progpath0;	/* Pathname by which program is invoked */
{
    char *progpath, *progname;
    int i, lpath;

    lpath = (strlen (progpath0) + 2) / 8;
    lpath = (lpath + 1) * 8;;
    progpath = (char *) calloc (lpath, 1);
    strcpy (progpath, progpath0);
    progname = progpath;
    for (i = strlen (progpath); i > -1; i--) {
        if (progpath[i] > 63 && progpath[i] < 90)
            progpath[i] = progpath[i] + 32;
        if (progpath[i] == '/') {
            progname = progpath + i + 1;
            break;
            }
	}
    return (progname);
}


char *
ProgCat (progname)

char *progname;	/* Program name which might contain catalog code */
{
    char *refcatname;
    refcatname = NULL;

    if (strsrch (progname,"gsc") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "gsc");
	}
    else if (strsrch (progname,"uac") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "uac");
	}
    else if (strsrch (progname,"ua1") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "ua1");
	}
    else if (strsrch (progname,"ua2") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "ua2");
	}
    else if (strsrch (progname,"usac") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "usac");
	}
    else if (strsrch (progname,"usa1") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "usa1");
	}
    else if (strsrch (progname,"usa2") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "usa2");
	}
    else if (strsrch (progname,"ujc") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "ujc");
	}
    else if (strsrch (progname,"sao") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "sao");
	}
    else if (strsrch (progname,"ppm") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "ppm");
	}
    else if (strsrch (progname,"ira") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "iras");
	}
    else if (strsrch (progname,"ty") != NULL) {
	refcatname = (char *) calloc (1,8);
	if (strsrch (progname, "2") != NULL)
	    strcpy (refcatname, "tycho2");
	else
	    strcpy (refcatname, "tycho");
	}
    else if (strsrch (progname,"hip") != NULL) {
	refcatname = (char *) calloc (1,16);
	strcpy (refcatname, "hipparcos");
	}
    else if (strsrch (progname,"act") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "act");
	}
    else if (strsrch (progname,"bsc") != NULL) {
	refcatname = (char *) calloc (1,8);
	strcpy (refcatname, "bsc");
	}

    return (refcatname);
}


void
CatNum (refcat, nnfld, nndec, dnum, numstr)

int	refcat;		/* Catalog code */
int	nnfld;		/* Number of characters in number (from CatNumLen) */
int	nndec;		/* Number of decimal places ( >= 0) */
double	dnum;		/* Catalog number of source */
char	*numstr;	/* Formatted number (returned) */

{
    char nform[16];	/* Format for star number */

    /* USNO A1.0, A2.0, SA1.0, or SA2.0 Catalogs */
    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
	refcat == UAC  || refcat == UA1  || refcat == UA2) {
	if (nnfld < 0)
	    sprintf (numstr, "%013.8f", dnum);
	else
	    sprintf (numstr, "%13.8f", dnum);
	}

    /* USNO Plate Catalog */
    else if (refcat == USNO) {
	if (nnfld < 0)
	    sprintf (numstr, "%07d", (int)(dnum+0.5));
	else
	    sprintf (numstr, "%7d", (int)(dnum+0.5));
	}

    /* USNO UJ 1.0 Catalog */
    else if (refcat == UJC) {
	if (nnfld < 0)
	    sprintf (numstr, "%012.7f", dnum);
	else
	    sprintf (numstr, "%12.7f", dnum);
	}

    /* HST Guide Star Catalog */
    else if (refcat == GSC) {
	if (nnfld < 0)
	    sprintf (numstr, "%09.4f", dnum);
	else
	    sprintf (numstr, "%9.4f", dnum);
	}

    /* SAO, PPM, or IRAS Point Source Catalogs (TDC binary format) */
    else if (refcat==SAO || refcat==PPM || refcat==IRAS || refcat==BSC) {
	if (nnfld < 0)
	    sprintf (numstr, "%06d", (int)(dnum+0.5));
	else
	    sprintf (numstr, "%6d", (int)(dnum+0.5));
	}

    /* Tycho, Hipparcos, or ACT catalogs */
    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==HIP || refcat==ACT) {
	if (nnfld < 0)
	    sprintf (numstr, "%010.5f", dnum);
	else
	    sprintf (numstr, "%10.5f", dnum);
	}

    /* Starbase tab-separated, TDC binary, or TDC ASCII catalogs */
    else if (nndec > 0) {
	if (nnfld > 0)
	    sprintf (nform,"%%%d.%df", nnfld, nndec);
	else if (nnfld < 0)
	    sprintf (nform,"%%0%d.%df", -nnfld, nndec);
	else
	    sprintf (nform,"%%%d.%df", nndec+5, nndec);
	sprintf (numstr, nform, dnum);
	}
    else if (nnfld > 10) {
	sprintf (nform,"%%%d.0f", nnfld);
	sprintf (numstr, nform, dnum+0.49);
	}
    else if (nnfld > 0) {
	sprintf (nform,"%%%dd", nnfld);
	sprintf (numstr, nform, (int)(dnum+0.49));
	}
    else if (nnfld < 0) {
	sprintf (nform,"%%0%dd", -nnfld);
	sprintf (numstr, nform, (int)(dnum+0.49));
	}
    else
	sprintf (numstr, "%6d", (int)(dnum+0.49));

    return;
}


int
CatNumLen (refcat, maxnum, nndec)

int	refcat;		/* Catalog code */
double	maxnum;		/* Maximum ID number */
			/* (Ignored for standard catalogs) */
int	nndec;		/* Number of decimal places ( >= 0) */

{
    int ndp;		/* Number of characters for decimal point */

    /* USNO A1.0, A2.0, SA1.0, or SA2.0 Catalogs */
    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
	refcat == UAC  || refcat == UA1  || refcat == UA2)
	return (13);

    /* USNO Plate Catalogs */
    else if (refcat == USNO)
	return (7);

    /* USNO UJ 1.0 Catalog */
    else if (refcat == UJC)
	return (12);

    /* HST Guide Star Catalog */
    else if (refcat == GSC)
	return (9);

    /* SAO, PPM, or IRAS Point Source Catalogs (TDC binary format) */
    else if (refcat==SAO || refcat==PPM || refcat==IRAS || refcat==BSC)
	return (6);

    /* Tycho, Tycho2, Hipparcos, or ACT catalogs */
    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==HIP || refcat==ACT)
	return (10);

    /* Starbase tab-separated, TDC binary, or TDC ASCII catalogs */
    else {
	if (nndec > 0)
	    ndp = 1;
	else
	    ndp = 0;
	if (maxnum < 10.0)
	    return (1 + nndec + ndp);
	else if (maxnum < 100.0)
	    return (2 + nndec + ndp);
	else if (maxnum < 1000.0)
	    return (3 + nndec + ndp);
	else if (maxnum < 10000.0)
	    return (4 + nndec + ndp);
	else if (maxnum < 100000.0)
	    return (5 + nndec + ndp);
	else if (maxnum < 1000000.0)
	    return (6 + nndec + ndp);
	else if (maxnum < 10000000.0)
	    return (7 + nndec + ndp);
	else if (maxnum < 100000000.0)
	    return (8 + nndec + ndp);
	else if (maxnum < 1000000000.0)
	    return (9 + nndec + ndp);
	else if (maxnum < 10000000000.0)
	    return (10 + nndec + ndp);
	else if (maxnum < 100000000000.0)
	    return (11 + nndec + ndp);
	else if (maxnum < 1000000000000.0)
	    return (12 + nndec + ndp);
	else if (maxnum < 10000000000000.0)
	    return (13 + nndec + ndp);
	else
	    return (14 + nndec + ndp);
	}
}


/* Return number of decimal places in catalogued numbers, if known */
int
CatNdec (refcat)

int	refcat;		/* Catalog code */

{
    /* USNO A1.0, A2.0, SA1.0, or SA2.0 Catalogs */
    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
	refcat == UAC  || refcat == UA1  || refcat == UA2)
	return (8);

    /* USNO Plate Catalogs */
    else if (refcat == USNO)
	return (0);

    /* USNO UJ 1.0 Catalog */
    else if (refcat == UJC)
	return (7);

    /* HST Guide Star Catalog */
    else if (refcat == GSC)
	return (4);

    /* SAO, PPM, or IRAS Point Source Catalogs (TDC binary format) */
    else if (refcat==SAO || refcat==PPM || refcat==IRAS || refcat==BSC)
	return (0);

    /* Tycho, Tycho2, Hipparcos, or ACT catalogs */
    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==HIP || refcat==ACT)
	return (5);

    /* Starbase tab-separated, TDC binary, or TDC ASCII catalogs */
    else
	return (-1);
}


/* Return number of decimal places in numeric string (-1 if not number) */

int
StrNdec (string)

char *string;	/* Numeric string */
{
    char *cdot;
    int lstr;

    if (notnum (string))
	return (-1);
    else {
	lstr = strlen (string);
	if ((cdot = strchr (string, '.')) == NULL)
	    return (0);
	else
	    return (lstr - (cdot - string));
	}
}


void
SearchLim (cra, cdec, dra, ddec, syscoor, ra1, ra2, dec1, dec2, verbose)

double	cra, cdec;	/* Center of search area  in degrees */
double	dra, ddec;	/* Horizontal and verticla half-widths of area */
int	syscoor;	/* Coordinate system */
double	*ra1, *ra2;	/* Right ascension limits in degrees */
double	*dec1, *dec2;	/* Declination limits in degrees */
int	verbose;	/* 1 to print limits, else 0 */

{
    double dec;

    /* Set right ascension limits for search */
    *ra1 = cra - dra;
    *ra2 = cra + dra;

    /* Keep right ascension between 0 and 360 degrees */
    if (syscoor != WCS_XY) {
	if (*ra1 < 0.0)
	    *ra1 = *ra1 + 360.0;
	if (*ra2 > 360.0)
	    *ra2 = *ra2 - 360.0;
	}

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
    if (syscoor != WCS_XY) {
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
	}

    if (verbose) {
	char rstr1[16],rstr2[16],dstr1[16],dstr2[16];
	if (syscoor == WCS_XY) {
	    num2str (rstr1, *ra1, 10, 5);
            num2str (dstr1, *dec1, 10, 5);
	    num2str (rstr2, *ra2, 10, 5);
            num2str (dstr2, *dec2, 10, 5);
	    }
	else {
	    ra2str (rstr1, 16, *ra1, 3);
            dec2str (dstr1, 16, *dec1, 2);
	    ra2str (rstr2, 16, *ra2, 3);
            dec2str (dstr2, 16, *dec2, 2);
	    }
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


/* RANGEINIT -- Initialize range structure from string */

struct Range *
RangeInit (string, ndef)

char	*string;	/* String containing numbers separated by , and - */
int	ndef;		/* Maximum allowable range value */

{
    struct Range *range;
    int ip, irange;
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
	else if (string[ip] > (char)47 && string[ip] < 58) {
	    first = strtod (string+ip, &slast);
	    ip = slast - string;
	    }
	else if (strchr ("-:x", string[ip]) == NULL) {
	    free (range);
	    return (NULL);
	    }

	/* Skip delimiters */
	while (string[ip] == ' ' || string[ip] == '	' ||
	       string[ip] == ',')
	    ip++;

	/* Get last limit
	* Must be '-', or 'x' otherwise last = first */
	if (string[ip] == '-' || string[ip] == ':') {
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
	    else if (string[ip] != 'x')
		last = first + ndef;
	    }
	else if (string[ip] != 'x')
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
	    else if (string[ip] != '-' && string[ip] != ':')
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


/* ISRANGE -- Return 1 if string is a range, else 0 */

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


/* RSTART -- Restart at beginning of range */

void
rstart (range)

struct Range *range;	/* Range structure */

{
    range->irange = -1;
    return;
}


/* RGETN -- Return number of values from range structure */

int
rgetn (range)

struct Range *range;	/* Range structure */

{
    return (range->nvalues);
}


/*  RGETR8 -- Return next number from range structure as 8-byte f.p. number */

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


/*  RGETI4 -- Return next number from range structure as 4-byte integer */

int
rgeti4 (range)

struct Range *range;	/* Range structure */

{
    double value;

    value = rgetr8 (range);
    return ((int) (value + 0.000000001));
}


/* -- SETOKEN -- tokenize a string for easy decoding */

int
setoken (tokens, string, cwhite)

struct Tokens *tokens;	/* Token structure returned */
char	*string;	/* character string to tokenize */
char	*cwhite;	/* additional whitespace characters
			 * if = tab, disallow spaces and commas */
{
    char squote,dquote,jch, newline;
    char *iq, *stri, *wtype, *str0, *inew;
    int i,j,naddw;

    newline = (char) 10;
    squote = (char) 39;
    dquote = (char) 34;
    if (string == NULL)
	return (0);

    /* Line is terminated by newline or NULL */
    inew = strchr (string, newline);
    if (inew != NULL)
	tokens->lline = inew - string - 1;
    else
	tokens->lline = strlen (string);

    /* Save current line in structure */
    tokens->line = string;

    /* Add extra whitespace characters */
    if (cwhite == NULL)
	naddw = 0;
    else
	naddw = strlen (cwhite);

    /* if character is tab, allow only tabs and nulls as separators */
    if (naddw > 0 && !strncmp (cwhite, "tab", 3)) {
	tokens->white[0] = (char) 9;
	tokens->white[1] = (char) 0;
	tokens->nwhite = 2;
	}

    /* otherwise, allow spaces, tabs, commas, nulls, and cwhite */
    else {
	tokens->nwhite = 3 + naddw;;
	tokens->white[0] = ' ';
	tokens->white[1] = (char) 9;
	tokens->white[2] = ',';
	tokens->white[3] = (char) 0;
	if (tokens->nwhite > 20)
	    tokens->nwhite = 20;
	if (naddw > 0) {
	    i = 0;
	    for (j = 3; j < tokens->nwhite; j++) {
		tokens->white[j] = cwhite[i];
		i++;
		}
	    }
	}
    tokens->white[tokens->nwhite] = (char) 0;

    tokens->ntok = 0;
    tokens->itok = 0;
    iq = string - 1;
    for (i = 0; i < MAXTOKENS; i++) {
	tokens->tok1[i] = NULL;
	tokens->ltok[i] = 0;
	}

    /* Process string one character at a time */
    stri = string;
    str0 = string;
    while (stri < string+tokens->lline) {

	/* Keep stuff between quotes in one token */
	if (stri <= iq)
	    continue;
	jch = *stri;

	/* Handle quoted strings */
	if (jch == squote)
	    iq = strchr (stri+1, squote);
	else if (jch == dquote)
	    iq = strchr (stri+1, dquote);
	else
	    iq = stri;
	if (iq > stri) {
	    tokens->ntok = tokens->ntok + 1;
	    if (tokens->ntok < MAXTOKENS) return (MAXTOKENS);
	    tokens->tok1[tokens->ntok] = stri + 1;
	    tokens->ltok[tokens->ntok] = (iq - stri) - 1;
	    stri = iq + 1;
	    str0 = iq + 1;
	    continue;
	    }

	/* Search for unquoted tokens */
	wtype = strchr (tokens->white, jch);

	/* If this is one of the additional whitespace characters,
	 * pass as a separate token */
	if (wtype > tokens->white + 2) {

	    /* Terminate token before whitespace */
	    if (stri > str0) {
		tokens->ntok = tokens->ntok + 1;
		if (tokens->ntok > MAXTOKENS) return (MAXTOKENS);
		tokens->tok1[tokens->ntok] = str0;
		tokens->ltok[tokens->ntok] = stri - str0;
		}

	    /* Make whitespace character next token; start new one */
	    tokens->ntok = tokens->ntok + 1;
	    if (tokens->ntok < MAXTOKENS) return (MAXTOKENS);
	    tokens->tok1[tokens->ntok] = stri;
	    tokens->ltok[tokens->ntok] = 1;
	    stri++;
	    str0 = stri;
	    }

	/* Pass previous token if regular whitespace or NULL */
	else if (wtype != NULL || jch == (char) 0) {

	    /* Ignore leading whitespace */
	    if (stri == str0) {
		stri++;
		str0 = stri;
		}

	    /* terminate token before whitespace; start new one */
	    else {
		tokens->ntok = tokens->ntok + 1;
		if (tokens->ntok > MAXTOKENS) return (MAXTOKENS);
		tokens->tok1[tokens->ntok] = str0;
		tokens->ltok[tokens->ntok] = stri - str0;
		stri++;
		str0 = stri;
		}
	    }

	/* Keep going if not whitespace */
	else
	    stri++;
	}

    /* Add token terminated by end of line */
    if (str0 < stri) {
	tokens->ntok = tokens->ntok + 1;
	if (tokens->ntok > MAXTOKENS)
	    return (MAXTOKENS);
	tokens->tok1[tokens->ntok] = str0;
	tokens->ltok[tokens->ntok] = stri - str0 + 1;
	}

    tokens->itok = 0;

    return (tokens->ntok);
}


/* NEXTOKEN -- get next token from tokenized string */

int
nextoken (tokens, token)
 
struct Tokens *tokens;	/* Token structure returned */
char	*token;		/* token (returned) */
{
    int it, ltok;

    tokens->itok = tokens->itok + 1;
    it = tokens->itok;
    if (it > tokens->ntok)
	it = tokens->ntok;
    else if (it < 1)
	it = 1;
    ltok = tokens->ltok[it];
    strncpy (token, tokens->tok1[it], ltok);
    token[ltok] = (char) 0;
    return (ltok);
}


/* GETOKEN -- get specified token from tokenized string */

int
getoken (tokens, itok, token)

struct Tokens *tokens;	/* Token structure returned */
int itok;		/* token sequence number of token
			 * if <0, get whole string after token -itok
			 * if =0, get whole string */
char *token;		/* token (returned) */
{
    int ltok;		/* length of token string (returned) */
    int it;

    it = itok;
    if (it > 0 ) {
	if (it > tokens->ntok)
	    it = tokens->ntok;
	ltok = tokens->ltok[it];
	strncpy (token, tokens->tok1[it], ltok);
	}
    else if (it < 0) {
	if (it < -tokens->ntok)
	    it  = -tokens->ntok;
	ltok = tokens->line + tokens->lline - tokens->tok1[-it];
	strncpy (token, tokens->tok1[-it], ltok);
	}
    else {
	ltok = tokens->lline;
	strncpy (token, tokens->tok1[1], ltok);
	}
    token[ltok] = (char) 0;

    return (ltok);
}


/* ISFILE -- Return 1 if file is a readable file, else 0 */

int
isfile (filename)

char    *filename;      /* Name of file for which to find size */
{
    if (access (filename, R_OK))
	return (0);
    else
	return (1);
}


/* AGETS -- Get keyword value from ASCII string with keyword=value anywhere */

int
agets (string, keyword0, lval, value)

char *string;  /* character string containing <keyword>= <value> info */
char *keyword0;  /* character string containing the name of the keyword
                   the value of which is returned.  hget searches for a
                   line beginning with this string.  if "[n]" or ",n" is
		   present, the n'th token in the value is returned. */
int lval;       /* Size of str in characters */
char *value;      /* String (returned) */
{
    char skey[16];
    char keyword[81];
    char *pval;
    char squot[2], dquot[2], lbracket[2], rbracket[2], comma[2];
    char *lastval, *rval, *brack1, *brack2, *lastring;
    int ipar, i;

    squot[0] = (char) 39;
    squot[1] = (char) 0;
    dquot[0] = (char) 34;
    dquot[1] = (char) 0;
    lbracket[0] = (char) 91;
    lbracket[1] = (char) 0;
    comma[0] = (char) 44;
    comma[1] = (char) 0;
    rbracket[0] = (char) 93;
    rbracket[1] = (char) 0;
    lastring = string + strlen (string);

    /* Find length of variable name */
    strncpy (keyword,keyword0, sizeof(keyword)-1);
    brack1 = strsrch (keyword,lbracket);
    if (brack1 == NULL)
	brack1 = strsrch (keyword,comma);
    if (brack1 != NULL) {
	*brack1 = '\0';
	brack1++;
	}

    /* First look for keyword= value */
    skey[0] = ' ';
    skey[1] = (char) 0;
    strcat (skey, keyword);
    strcat (skey, "=");
    pval = strsrch (string, skey);

    /* Look for keyword = value */
    if (pval == NULL) {
	skey[0] = ' ';
	skey[1] = (char) 0;
	strcat (skey, keyword);
	strcat (skey, " =");
	pval = strsrch (string, skey);
	}

    /* Look for keyword:  value */
    if (pval == NULL) {
	skey[0] = ' ';
	skey[1] = (char) 0;
	strcat (skey, keyword);
	strcat (skey, ": ");
	pval = strsrch (string, skey);
	}

    /* If keyword has not been found, return 0 */
    if (pval == NULL)
	return (0);

    /* Otherwise, bump pointer past keyword */
    else
	pval = pval + strlen (skey);

    /* Drop leading spaces */
    while (*pval == ' ') pval++;

    /* If keyword has brackets, figure out which token to extract */
    if (brack1 != NULL) {
        brack2 = strsrch (brack1,rbracket);
        if (brack2 != NULL)
            *brack2 = '\0';
        ipar = atoi (brack1);
	}
    else
	ipar = 1;

    /* Move to appropriate token */
    for (i = 1; i < ipar; i++) {
	while (*pval != ' ' && *pval != '/' && pval < lastring)
	    pval++;

	/* Drop leading spaces  or / */
	while (*pval == ' ' || *pval == '/')
	    pval++;
	}

    /* Transfer token value to returned string */
    rval = value;
    lastval = value + lval - 1;
    while (*pval != ' ' && *pval != '\n' && *pval != '/' &&
	   pval < lastring && rval < lastval)
	*rval++ = *pval++;
    if (rval < lastval)
	*rval = (char) 0;
    else
	*lastval = 0;

    return (1);
}

char sptbv[468]={"O5O8B0B0B0B1B1B1B2B2B2B3B3B3B4B5B5B6B6B6B7B7B8B8B8B9B9B9B9A0A0A0A0A0A0A0A0A0A2A2A2A2A2A2A2A2A5A5A5A5A6A7A7A7A7A7A7A7A7A7A7F0F0F0F0F0F0F0F2F2F2F2F2F2F2F5F5F5F5F5F5F5F5F5F8F8F8F8F8F8G0G5G5G2G2G2G3G3G4G4G5G5G5G6G6G6G6G6K6K6K6K6K7K7K7K7K7K7K7K7K7K7K7K7K7K7K8K8K8K8K8K8K8K8K8K8K8K8K8K8K8K8K8K8K8K5K5K5K5K5K6K6K6K6K6K6K6K7K7K7K7K7K7K7K8K8K8K8K9K9K9M0M0M0M0M0M0M1M1M1M1M1M2M2M2M2M3M3M4M4M5M5M5M2M2M2M3M3M4M4M5M5M5M6M6M6M6M6M6M6M6M6M7M7M7M7M7M7M7M7M7M7M7M7M7M7M8M8M8M8M8M8M8"};

void
bv2sp (bv, b, v, isp)

double	*bv;	/* B-V Magnitude */
double	b;	/* B Magnitude used if bv is NULL */
double	v;	/* V Magnitude used if bv is NULL */
char	*isp;	/* Spectral type */
{
    double bmv;	/* B - V magnitude */
    int im;

    if (bv == NULL)
	bmv = b - v;
    else
	bmv = *bv;

    if (bmv < -0.32) {
	isp[0] = '_';
	isp[1] = '_';
	}
    else if (bmv > 2.00) {
	isp[0] = '_';
	isp[1] = '_';
	}
    else if (bmv < 0) {
	im = 2 * (32 + (int)(bmv * 100.0 - 0.5));
	isp[0] = sptbv[im];
	isp[1] = sptbv[im+1];
	}
    else {
	im = 2 * (32 + (int)(bmv * 100.0 + 0.5));
	isp[0] = sptbv[im];
	isp[1] = sptbv[im+1];
	}
    return;
}

char sptbr1[96]={"O5O8O9O9B0B0B0B0B0B1B1B1B2B2B2B2B2B3B3B3B3B3B3B5B5B5B5B6B6B6B7B7B7B7B8B8B8B8B8B9B9B9B9B9A0A0A0"};

char sptbr2[904]={"A0A0A0A0A0A0A0A0A2A2A2A2A2A2A2A2A2A2A2A2A2A2A2A5A5A5A5A5A5A5A5A5A5A5A7A7A7A7A7A7A7A7A7A7A7A7A7A7A7A7F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F2F2F2F2F2F2F2F2F2F2F2F5F5F5F5F5F5F5F5F5F5F5F5F5F5F8F8F8F8F8F8F8F8F8F8F8F8F8F8G0G0G0G0G0G0G0G0G2G2G2G2G2G5G5G5G5G5G5G5G5G8G8G8G8G8G8G8G8G8G8G8G8G8G8K0K0K0K0K0K0K0K0K0K0K0K0K0K0K0K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K2K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K5K7K7K7K7K7K7K7K7K7K7K7K7K7K7K7K7K7K7K7K7K7K7K7K7K7M0M0M0M0M0M0M0M0M0M0M0M0M0M0M0M0M0M0M0M0M0M0M0M0M1M1M1M1M1M1M1M1M1M1M1M1M1M1M1M2M2M2M2M2M2M2M2M2M2M2M2M2M2M2M3M3M3M3M3M3M3M3M3M3M3M4M4M4M4M4M4M4M4M4M4M4M4M4M4M5M5M5M5M5M5M5M5M5M5M5M5M5M5M5M5M5M5M5M5M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M6M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M7M8M8M8M8M8M8M8M8M8M8M8M8M8M8M8M8M8M8M8M8M8M8M8M8"};

void
br2sp (br, b, r, isp)

double	*br;	/* B-R Magnitude */
double	b;	/* B Magnitude used if br is NULL */
double	r;	/* R Magnitude used if br is NULL */
char	*isp;	/* Spectral type */
{
    double bmr;	/* B - R magnitude */
    int im;

    if (br == NULL)
	bmr = b - r;
    else
	bmr = *br;

    if (b == 0.0 && r > 2.0) {
	isp[0] = '_';
	isp[1] = '_';
	}
    else if (bmr < -0.47) {
	isp[0] = '_';
	isp[1] = '_';
	}
    else if (bmr > 4.50) {
	isp[0] = '_';
	isp[1] = '_';
	}
    else if (bmr < 0) {
	im = 2 * (47 + (int)(bmr * 100.0 - 0.5));
	isp[0] = sptbr1[im];
	isp[1] = sptbr1[im+1];
	}
    else {
	im = 2 * ((int)(bmr * 100.0 + 0.49));
	isp[0] = sptbr2[im];
	isp[1] = sptbr2[im+1];
	}
    return;
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
 * Jul 15 1999	Add getfilebuff()
 * Jul 23 1999	Add Bright Star Catalog
 * Aug 16 1999	Add RefLim() to set catalog search limits
 * Sep 21 1999	In isrange(), check for x
 * Oct  5 1999	Add setoken(), nextoken(), and getoken()
 * Oct 15 1999	Fix format eror in error message
 * Oct 20 1999	Use strchr() in range decoding
 * Oct 21 1999	Fix declarations after lint
 * Oct 21 1999	Fix arguments to catopen() and catclose() after lint
 * Nov  3 1999	Fix bug which lost last character on a line in getoken
 * Dec  9 1999	Add next_token(); set pointer to next token in first_token
 *
 * Jan 11 2000	Use nndec for Starbase files, too
 * Feb 10 2000	Read coordinate system, epoch, and equinox from Starbase files
 * Mar  1 2000	Add isfile() to tell whether string is name of readable file
 * Mar  1 2000	Add agets() to return value from keyword = value in string
 * Mar  1 2000	Add isfile() to tell if a string is the name of a readable file
 * Mar  1 2000	Add agets() to read a parameter from a comment line of a file
 * Mar  8 2000	Add ProgCat() to return catalog flag from program name
 * Mar 13 2000	Add PropCat() to return whether catalog has proper motions
 * Mar 27 2000	Clean up code after lint
 * May 22 2000	Add bv2sp() to approximate main sequence spectral type from B-V
 * May 25 2000	Add Tycho 2 catalog
 * May 26 2000	Add field size argument to CatNum() and CatNumLen()
 * Jun  2 2000	Set proper motion for all catalog types in RefCat()
 * Jun 26 2000	Add XY image coordinate system
 * Jul 26 2000	Include math.h to get strtod() on SunOS machines
 * Aug  2 2000	Allow up to 14 digits in catalog IDs
 * Sep  1 2000	Add option in CatNum to print leading zeroes if nnfld > 0
 * Sep 22 2000	Add br2sp() to approximate main sequence spectral type from B-R
 * Oct 24 2000	Add USNO option to RefCat()
 * Nov 21 2000	Clean up logic in RefCat()
 * Nov 28 2000	Try PPMra and SAOra in RefCat() as well as PPM and SAO
 * Dec 13 2000	Add StrNdec() to get number of decimal places in star numbers
 */
