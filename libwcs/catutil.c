/* File libwcs/catutil.c
 * June 3, 1999
 * By Doug Mink
 */

/* int RefCat()			Return catalog type code, title, coord. system
 * struct Range RangeInit()	Return structure containing ranges of numbers
 * int rgetn();			Return number of values from range structure
 * int rgeti4();		Return next number from range structure
 * int rgetr8();		Return next number from range structure
 * void CatNum (refcat, nndec, dnum, numstr) returns formatted source number
 * int CatNumLen (refcat, nndec) returns length of source numbers
 * dt2jd (date,time)
 *	convert date as yyyy.ddmm and time as hh.mmsss to Julian date
 * jd2dt (dj,date,time)
 *	convert Julian date to date as yyyy.mmdd and time as hh.mmssss
 * dt2ts (date,time)
 *	convert date as yyyy.ddmm and time as hh.mmsss to seconds since 1950.0
 * ts2dt (tsec,date,time)
 *	convert seconds since 1950.0 to date as yyyy.ddmm and time as hh.mmsss
 * dt2ep (date, time)
 *	convert date as yyyy.ddmm and time as hh.mmsss to fractional year
 * ep2dt (epoch,date, time)
 *	convert fractional year to date as yyyy.ddmm and time as hh.mmsss
 * jd2ts (dj)
 *	convert Julian day to seconds since 1950.0
 * ts2jd (tsec)
 *	convert seconds since 1950.0 to Julian day
 * jd2ep (dj)
 *	convert Julian date to fractional year as used in epoch
 * ep2jd (epoch)
 *	convert fractional year as used in epoch to Julian date
 * ts2i (tsec,iyr,imon,iday,ihr,imn,sec)
 *	convert sec since 1950.0 to year month day hours minutes seconds .1ms
 * dint (dnum)
 *	returns integer part of floating point number
 * dmod (dnum)
 *	returns Mod of floating point number
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "wcs.h"
#include "wcscat.h"

double dt2jd();
double dt2ts();
void ep2dt();
void ts2dt();
void ts2i();
double dint();
double dmod();

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

struct Range *
RangeInit (string, ndef)

char	*string;	/* String containing numbers separated by , and - */
int	ndef;		/* Maximum allowable range value */

{
    struct Range *range;
    int ip, nvalues, irange;
    char *slast;
    double first, last, step;

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

int
rgetn (range)

struct Range *range;	/* Range structure */

{
    return (range->nvalues);
}


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

int
rgeti4 (range)

struct Range *range;	/* Range structure */

{
    double value;

    value = rgetr8 (range);
    return ((int) (value + 0.000000001));
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
    else if (refcat==SAO || refcat==PPM || refcat==IRAS )
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
    else if (refcat==SAO || refcat==PPM || refcat==IRAS )
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


/* DT2JD-- convert from date as yyyy.mmdd and time as hh.mmsss to Julian Date */

double
dt2jd (date,time)

double	date;		/* Date as yyyy.mmdd
			    yyyy = calendar year (e.g. 1973)
			    mm = calendar month (e.g. 04 = april)
			    dd = calendar day (e.g. 15) */
double	time;		/* Time as hh.mmssxxxx
			    *if time<0, it is time as -(fraction of a day)
			    hh = hour of day (0 .le. hh .le. 23)
			    nn = minutes (0 .le. nn .le. 59)
			    ss = seconds (0 .le. ss .le. 59)
			  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999)
			*/
{
    double dj;		/* Julian date (returned) */
    double tsec;	/* seconds since 1950.0 */

    if (date < 0.0301)
	return (0.0);
    tsec = dt2ts (date,time);
    dj = 2433282.50 + (tsec / 86400.0);

    return (dj);
}


/* jd2dt-- convert Julian date to date as yyyy.mmdd and time as hh.mmssss */

void
jd2dt (dj,date,time)

double	dj;		/* Julian date */
double	*date;		/* Date as yyyy.mmdd (returned)
			    yyyy = calendar year (e.g. 1973)
			    mm = calendar month (e.g. 04 = april)
			    dd = calendar day (e.g. 15) */
double	*time;		/* Time as hh.mmssxxxx (returned)
			    *if time<0, it is time as -(fraction of a day)
			    hh = hour of day (0 .le. hh .le. 23)
			    nn = minutes (0 .le. nn .le. 59)
			    ss = seconds (0 .le. ss .le. 59)
			  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999)
			*/
{
    double tsec;

    tsec = (dj - 2433282.5) * 86400.0;
    ts2dt (tsec, date, time);

    return;
}


/* JD2EP-- convert Julian date to fractional year as used in epoch */

double
jd2ep (dj)

double	dj;		/* Julian date */

{
    double epoch;	/* Date as fractional year (returned) */
    double dj0, dj1, date0, time0, date1, date, time;

    jd2dt (dj, &date, &time);
    time0 = 12.0;
    date0 = dint (date) + 0.0101;
    date1 = dint (date) + 1.0101;
    dj0 = dt2jd (date0, time0);
    dj1 = dt2jd (date1, time0);
    epoch = dint (date) + ((dj - dj0) / (dj1 - dj0));
    return (epoch);
}


/* EP2JD-- convert fractional year as used in epoch to Julian date */

double
ep2jd (epoch)

double	epoch;		/* Date as fractional year */

{
    double dj;		/* Julian date (returned)*/
    double dj0, dj1, date0, time0, date1, date, time;

    ep2dt (epoch, &date, &time);
    dj = dt2jd (date, time);
    return (dj);
}


/* JD2TS-- convert Julian date to seconds since 1950.0 */

double
jd2ts (dj)

double	dj;		/* Julian date */
{
    double tsec;	/* seconds since 1950.0 (returned) */

    tsec = (dj - 2433282.5) * 86400.0;
    return (tsec);
}


/* TS2JD-- convert seconds since 1950.0 to Julian date */

double
ts2jd (tsec)

double	tsec;		/* seconds since 1950.0 */
{
    double dj;		/* Julian date (returned) */

    dj = 2433282.5 + (tsec / 8.6400);
    return (dj);
}


/* DT2EP-- convert from date, time as yyyy.mmdd hh.mmsss to fractional year */

double
dt2ep (date, time)

double	date;		/* Date as yyyy.mmdd
			    yyyy = calendar year (e.g. 1973)
			    mm = calendar month (e.g. 04 = april)
			    dd = calendar day (e.g. 15) */
double	time;		/* Time as hh.mmssxxxx
			    *if time<0, it is time as -(fraction of a day)
			    hh = hour of day (0 .le. hh .le. 23)
			    nn = minutes (0 .le. nn .le. 59)
			    ss = seconds (0 .le. ss .le. 59)
			  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999)
			*/
{
    double epoch;	/* Date as fractional year (returned) */
    double dj, dj0, dj1, date0, time0, date1;

    dj = dt2jd (date, time);
    time0 = 12.0;
    date0 = dint (date) + 0.0101;
    date1 = dint (date) + 1.0101;
    dj0 = dt2jd (date0, time0);
    dj1 = dt2jd (date1, time0);
    epoch = dint (date) + ((dj - dj0) / (dj1 - dj0));
    return (epoch);
}


/* EP2DT-- convert from fractional year to date, time as yyyy.mmdd hh.mmsss */

void
ep2dt (epoch, date, time)

double epoch;		/* Date as fractional year */
double	*date;		/* Date as yyyy.mmdd (returned)
			    yyyy = calendar year (e.g. 1973)
			    mm = calendar month (e.g. 04 = april)
			    dd = calendar day (e.g. 15) */
double	*time;		/* Time as hh.mmssxxxx (returned)
			    *if time<0, it is time as -(fraction of a day)
			    hh = hour of day (0 .le. hh .le. 23)
			    nn = minutes (0 .le. nn .le. 59)
			    ss = seconds (0 .le. ss .le. 59)
			  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999)
			*/
{
    double dj, dj0, dj1, date0, time0, date1, epochi, epochf;

    time0 = 12.0;
    epochi = dint (epoch);
    epochf = epoch - epochi;
    date0 = epochi + 0.0101;
    date1 = epochi + 1.0101;
    dj0 = dt2jd (date0, time0);
    dj1 = dt2jd (date1, time0);
    dj = dj0 + epochf * (dj1 - dj0);
    jd2dt (dj, date, time);
    return;
}


/* DT2TS-- convert from date, time as yyyy.mmdd hh.mmsss to sec since 1950.0 */

double
dt2ts (date,time)

double	date;		/* Date as yyyy.mmdd
			    yyyy = calendar year (e.g. 1973)
			    mm = calendar month (e.g. 04 = april)
			    dd = calendar day (e.g. 15) */
double	time;		/* Time as hh.mmssxxxx
			    *if time<0, it is time as -(fraction of a day)
			    hh = hour of day (0 .le. hh .le. 23)
			    nn = minutes (0 .le. nn .le. 59)
			    ss = seconds (0 .le. ss .le. 59)
			  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999)
			*/
{
    double tsec;	/* Seconds past 1950.0 (returned) */

    double dh,dm,dd;
    int iy,im,id;

    if (date < 0.0301)
	return (0.0);

/* Calculate the number of full years, months, and days already
 * elapsed since 0h, march 1, -1 (up to most recent midnight). */

    /* convert time of day to elapsed seconds */

    if (time < 0.0)
	tsec = time * -86400.0;
    else {
	dh = (int) (time + 0.0000000001);
	dm = (int) (((time - dh) * 100.0) + 0.0000000001);
	tsec = (time * 10000.0) - (dh * 10000.0) - (dm * 100.0);
	tsec = (int) (tsec * 100000.0 + 0.0001) / 100000.0;
	tsec = tsec + (dm * 60.0) + (dh * 3600.0);
	}

    /* Calculate the number of full months elapsed since
     * the current or most recent March */
    iy = (int) (date + 0.0000000001);
    im = (int) (((date - (double) (iy)) * 10000.0) + 0.00000001);
    id = im % 100;
    im = (im / 100) + 9;
    if (im < 12) iy = iy - 1;
    im = im % 12;
    id = id - 1;

    /* starting with March as month 0 and ending with the following
     * February as month 11, the calculation of the number of days
     * per month reduces to a simple formula. the following statement
     * determines the number of whole days elapsed since 3/1/-1 and then
     * subtracts the 712163 days between then and 1/1/1950.  it converts
     * the result to seconds and adds the accumulated seconds above. */
    id = id + ((im+1+im/6+im/11)/2 * 31) + ((im-im/6-im/11)/2 * 30) +
	 (iy / 4) - (iy / 100) + (iy / 400);
    dd = (double) id + (365.0 * (double) iy) - 712163.0;
    tsec = tsec + (dd * 86400.0);

    return (tsec);
}


/* TS2DT-- convert seconds since 1950.0 to date, time as yyyy.mmdd hh.mmssss */

void
ts2dt (tsec,date,time)

double	tsec;		/* Seconds past 1950.0 */
double	*date;		/* Date as yyyy.mmdd (returned)
			    yyyy = calendar year (e.g. 1973)
			    mm = calendar month (e.g. 04 = april)
			    dd = calendar day (e.g. 15) */
double	*time;		/* Time as hh.mmssxxxx (returned)
			    *if time<0, it is time as -(fraction of a day)
			    hh = hour of day (0 .le. hh .le. 23)
			    nn = minutes (0 .le. nn .le. 59)
			    ss = seconds (0 .le. ss .le. 59)
			  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999)
			*/
{
    int iyr,imon,iday,ihr,imn;
    double sec;

    ts2i (tsec,&iyr,&imon,&iday,&ihr,&imn,&sec);

    /* Convert date to yyyy.mmdd */
    *date = (double) iyr + 0.01 * (double) imon + 0.0001 * (double) iday;

    /* Convert time to hh.mmssssss */
    *time = (double) ihr + 0.01 * (double) imn + 0.0001 * sec;

    return;
}


/* TS2I-- convert sec since 1950.0 to year month day hours minutes seconds .1ms */

void
ts2i (tsec,iyr,imon,iday,ihr,imn,sec)

double	tsec;		/* seconds since 1/1/1950 0:00 */
int	*iyr;		/* year (returned) */
int	*imon;		/* month (returned) */
int	*iday;		/* day (returned) */
int	*ihr;		/* hours (returned) */
int	*imn;		/* minutes (returned) */
double	*sec;		/* seconds (returned) */

{
    double t,days;
    int isec,ihms,nc,nc4,nly,ny,m,im;

    /* Time of day (hours, minutes, seconds, .1 msec) */
    t = dint ((tsec + 61530883200.0) * 10000.0 + 0.5);
    *ihr = (int) (dmod (t/36000000.0, 24.0));
    *imn = (int) (dmod (t/60000.0, 60.0));
    if (tsec >= 0) {
	ihms = (int) (dmod (tsec+0.000001, 1.0) * 10000.0);
	isec = (int) (dmod (tsec+0.000001, 60.0));
	}
    else {
	ihms = (int) (dmod (tsec-0.000001, 1.0) * 10000.0);
	isec = (int) (dmod (tsec-0.000001, 60.0));
	}

    /* Seconds */
    *sec = (double) isec + 0.0001 * (double) ihms;

    /* Number of days since 0 hr 0/0/0000 */
    days = dint ((t / 864000000.0) + 0.000001);

    /* Number of leap centuries (400 years) */
    nc4 = (int) ((days / 146097.0) + 0.00001);

    /* Number of centuries since last /400 */
    days = days - (146097.0 * (double) (nc4));
    nc = (int) ((days / 36524.0) + 0.0001);
    if (nc > 3) nc = 3;

    /* Number of leap years since last century */
    days = days - (36524.0 * nc);
    nly = (int) ((days / 1461.0) + 0.0000000001);

    /* Number of years since last leap year */
    days = days - (1461.0 * (double) nly);
    ny = (int) ((days / 365.0) + 0.00000001);
    if (ny > 3) ny = 3;

    /* Day of month */
    days = days - (365.0 * (double) ny);
    *iday = (int) (days + 0.00000001) + 1;
    for (m = 1; m <= 12; m++) {
	im = (m + ((m - 1) / 5)) % 2;
	if (*iday-1 < im+30) break;
	*iday = *iday - im - 30;
	}

    /* Month */
    *imon = ((m+1) % 12) + 1;

    /* Year */
    *iyr = nc4*400 + nc*100 + nly*4 + ny + m/11;

    return;
}

double
dint (dnum)

double	dnum;

{
    double dn;

    if (dnum < 0.0)
	dn = -floor (-dnum);
    else
	dn = floor (dnum);
    return (dn);
}


double
dmod (dnum, dm)

double	dnum, dm;
{
    double dnumx, dnumi, dnumf;
    if (dnum < 0.0)
	dnumx = -dnum;
    else
	dnumx = dnum;
    dnumi = dint (dnumx / dm);
    if (dnum < 0.0)
	dnumf = dnum + (dnumi * dm);
    else if (dnum > 0.0)
	dnumf = dnum - (dnumi * dm);
    else
	dnumf = 0.0;
    return (dnumf);
}

/* May 19 1999	New program, based on iolib/jcon.f and iolib/vcon.f
 */

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
 */
