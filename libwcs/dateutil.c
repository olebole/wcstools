/* File libwcs/dateutil.c
 * December 6, 1999
 * By Doug Mink
 */

/* Date and time conversion routines using the following conventions:
   dt = 2 floating point numbers: yyyy.mmdd, hh.mmssssss
   jd = Julian Date
   ts = seconds since 1950.0 (used for ephemeris computations
   ep = fractional year, often epoch of a position including proper motion
   fd = FITS date string (dd/mm/yy before 2000, then yyyy-mm-dd[Thh:mm:ss.ss])

 * dt2ep (date, time)
 *	convert date as yyyy.ddmm and time as hh.mmsss to fractional year
 * dt2fd (date, time)
 *	convert date as yyyy.ddmm and time as hh.mmsss to FITS date string
 * dt2jd (date,time)
 *	convert date as yyyy.ddmm and time as hh.mmsss to Julian date
 * dt2ts (date,time)
 *	convert date as yyyy.ddmm and time as hh.mmsss to seconds since 1950.0
 * ep2dt (epoch,date, time)
 *	convert fractional year to date as yyyy.ddmm and time as hh.mmsss
 * ep2fd (epoch, string)
 *	convert epoch to FITS ISO date string
 * ep2jd (epoch)
 *	convert fractional year as used in epoch to Julian date
 * ep2ts (epoch)
 *	convert fractional year to seconds since 1950.0
 * fd2jd (string)
 *	convert FITS standard date string to Julian date
 * fd2ep (string)
 *	convert FITS date string to fractional year
 * jd2dt (dj,date,time)
 *	convert Julian date to date as yyyy.mmdd and time as hh.mmssss
 * jd2ep (dj)
 *	convert Julian date to fractional year as used in epoch
 * jd2fd (epoch, string)
 *	convert Julian date to FITS ISO date string
 * jd2ts (dj)
 *	convert Julian day to seconds since 1950.0
 * ts2dt (tsec,date,time)
 *	convert seconds since 1950.0 to date as yyyy.ddmm and time as hh.mmsss
 * ts2fd (tsec)
 *	convert seconds since 1950.0 to FITS standard date string
 * ts2i (tsec,iyr,imon,iday,ihr,imn,sec, ndsec)
 *	convert sec since 1950.0 to year month day hours minutes seconds .1ms
 * ts2jd (tsec)
 *	convert seconds since 1950.0 to Julian day
 * dint (dnum)
 *	returns integer part of floating point number
 * dmod (dnum)
 *	returns Mod of floating point number
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsfile.h"

double dint();
double dmod();

/* DT2JD-- convert from date as yyyy.mmdd and time as hh.mmsss to Julian Date */

double
dt2jd (date,time)

double	date;	/* Date as yyyy.mmdd
		    yyyy = calendar year (e.g. 1973)
		    mm = calendar month (e.g. 04 = april)
		    dd = calendar day (e.g. 15) */
double	time;	/* Time as hh.mmssxxxx
		    *if time<0, it is time as -(fraction of a day)
		    hh = hour of day (0 .le. hh .le. 23)
		    nn = minutes (0 .le. nn .le. 59)
		    ss = seconds (0 .le. ss .le. 59)
		  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999) */
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

double	dj;	/* Julian date */
double	*date;	/* Date as yyyy.mmdd (returned)
		    yyyy = calendar year (e.g. 1973)
		    mm = calendar month (e.g. 04 = april)
		    dd = calendar day (e.g. 15) */
double	*time;	/* Time as hh.mmssxxxx (returned)
		    *if time<0, it is time as -(fraction of a day)
		    hh = hour of day (0 .le. hh .le. 23)
		    nn = minutes (0 .le. nn .le. 59)
		    ss = seconds (0 .le. ss .le. 59)
		  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999) */
{
    double tsec;

    tsec = (dj - 2433282.5) * 86400.0;
    ts2dt (tsec, date, time);

    return;
}


/* JD2EP-- convert Julian date to fractional year as used in epoch */

double
jd2ep (dj)

double	dj;	/* Julian date */

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


/* EP2FD-- convert fractional year to FITS date, yyyy-mm-ddThh:mm:ss.ss */

char *
ep2fd (epoch)

double	epoch;	/* Date as fractional year */
{
    double tsec; /* seconds since 1950.0 (returned) */

    tsec = ep2ts (epoch);

    return (ts2fd (tsec));
}


/* EP2TS-- convert fractional year to seconds since 1950.0 */

double
ep2ts (epoch)

double	epoch;	/* Date as fractional year */
{
    double tsec; /* seconds since 1950.0 (returned) */
    double dj;

    dj = ep2jd (epoch);

    tsec = (dj - 2433282.5) * 86400.0;
    return (tsec);
}



/* EP2JD-- convert fractional year as used in epoch to Julian date */

double
ep2jd (epoch)

double	epoch;	/* Date as fractional year */

{
    double dj;	/* Julian date (returned)*/
    double date, time;

    ep2dt (epoch, &date, &time);
    dj = dt2jd (date, time);
    return (dj);
}


/* JD2FD-- convert Julian date to FITS date, yyyy-mm-ddThh:mm:ss.ss */

char *
jd2fd (dj)

double	dj;	/* Julian date */
{
    double tsec; /* seconds since 1950.0 (returned) */

    tsec = (dj - 2433282.5) * 86400.0;

    return (ts2fd (tsec));
}


/* JD2TS-- convert Julian date to seconds since 1950.0 */

double
jd2ts (dj)

double	dj;	/* Julian date */
{
    double tsec; /* seconds since 1950.0 (returned) */

    tsec = (dj - 2433282.5) * 86400.0;
    return (tsec);
}


/* TS2JD-- convert seconds since 1950.0 to Julian date */

double
ts2jd (tsec)

double	tsec;	/* seconds since 1950.0 */
{
    double dj;	/* Julian date (returned) */

    dj = 2433282.5 + (tsec / 86400.0);
    return (dj);
}


/* DT2EP-- convert from date, time as yyyy.mmdd hh.mmsss to fractional year */

double
dt2ep (date, time)

double	date;	/* Date as yyyy.mmdd
		    yyyy = calendar year (e.g. 1973)
		    mm = calendar month (e.g. 04 = april)
		    dd = calendar day (e.g. 15) */
double	time;	/* Time as hh.mmssxxxx
		    *if time<0, it is time as -(fraction of a day)
		    hh = hour of day (0 .le. hh .le. 23)
		    nn = minutes (0 .le. nn .le. 59)
		    ss = seconds (0 .le. ss .le. 59)
		  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999) */
{
    double epoch; /* Date as fractional year (returned) */
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

double epoch;	/* Date as fractional year */
double	*date;	/* Date as yyyy.mmdd (returned)
		    yyyy = calendar year (e.g. 1973)
		    mm = calendar month (e.g. 04 = april)
		    dd = calendar day (e.g. 15) */
double	*time;	/* Time as hh.mmssxxxx (returned)
		    *if time<0, it is time as -(fraction of a day)
		    hh = hour of day (0 .le. hh .le. 23)
		    nn = minutes (0 .le. nn .le. 59)
		    ss = seconds (0 .le. ss .le. 59)
		  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999) */
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


/* FD2JD-- convert FITS standard date to Julian date */

double
fd2jd (string)

char *string;	/* FITS date string, which may be:
			fractional year
			dd/mm/yy (FITS standard before 2000)
			dd-mm-yy (nonstandard use before 2000)
			yyyy-mm-dd (FITS standard after 1999)
			yyyy-mm-ddThh:mm:ss.ss (FITS standard after 1999) */


{
    double epoch; /* Date as fractional year */
    double dj;	/* Julian date (returned)*/
    double date, time;

    epoch = fd2ep (string);
    ep2dt (epoch, &date, &time);
    dj = dt2jd (date, time);
    return (dj);
}



/* FD2EP -- convert from FITS standard date to fractional year */

double
fd2ep (string)

char *string;	/* FITS date string, which may be:
			yyyy.ffff (fractional year)
			dd/mm/yy (FITS standard before 2000)
			dd-mm-yy (nonstandard FITS use before 2000)
			yyyy-mm-dd (FITS standard after 1999)
			yyyy-mm-ddThh:mm:ss.ss (FITS standard after 1999) */

{
    double yeardays, seconds, fday, epoch;
    char *sstr, *dstr, *tstr, *cstr, *nval;
    int year, month, day, yday, i, hours, minutes;
    static int mday[12] = {31,28,31,30,31,30,31,31,30,31,30,31};

    /* Translate string from ASCII to binary */
    if (string != NULL) {
	sstr = strchr (string,'/');
	dstr = strchr (string,'-');

	/* Original FITS date format: dd/mm/yy */
	if (sstr > string) {
	    *sstr = '\0';
	    day = (int) atof (string);
	    nval = sstr + 1;
	    sstr = strchr (nval,'/');
	    if (sstr == NULL)
		sstr = strchr (nval,'-');
	    if (sstr > string) {
		*sstr = '\0';
		month = (int) atof (nval);
		nval = sstr + 1;
		year = (int) atof (nval);
		if (year >= 0 && year <= 49)
		    year = year + 2000;
		else if (year < 100)
		    year = year + 1900;
		if ((year % 4) == 0)
		    mday[1] = 29;
		else
		    mday[1] = 28;
		if ((year % 100) == 0 && (year % 400) != 0)
		    mday[1] = 28;
		if (day > mday[month-1])
		    day = mday[month-1];
		else if (day < 1)
		    day = 1;
		if (mday[1] == 28)
		    yeardays = 365.0;
		else
		    yeardays = 366.0;
		yday = day - 1;
		for (i = 0; i < month-1; i++)
		    yday = yday + mday[i];
		epoch = (double) year + ((double)yday / yeardays);
		return (epoch);
		}
	    else
		return (0.0);
	    }

	/* New FITS date format: yyyy-mm-ddThh:mm:ss[.sss] */
	else if (dstr > string) {
	    *dstr = '\0';
	    year = (int) atof (string);
	    nval = dstr + 1;
	    dstr = strchr (nval,'-');
	    month = 1;
	    day = 1;
	    tstr = NULL;
	    if (dstr > string) {
		*dstr = '\0';
		month = (int) atof (nval);
		nval = dstr + 1;
		tstr = strchr (nval,'T');
		if (tstr > string)
		    *tstr = '\0';
		day = (int) atof (nval);
		}

	    /* If year is < 32, it is really day of month in old format */
	    if (year < 32) {
		i = year;
		year = day + 1900;
		day = i;
		}

	    if ((year % 4) == 0)
		mday[1] = 29;
	    else
		mday[1] = 28;
	    if ((year % 100) == 0 && (year % 400) != 0)
		mday[1] = 28;
	    if (day > mday[month-1])
		day = mday[month-1];
	    else if (day < 1)
		day = 1;
	    if (mday[1] == 28)
		yeardays = 365.0;
	    else
		yeardays = 366.0;
	    yday = day - 1;
	    for (i = 0; i < month-1; i++)
		yday = yday + mday[i];
	    epoch = (double) year + ((double)yday / yeardays);

	    /* Extract time, if it is present */
	    if (tstr > string) {
		nval = tstr + 1;
		hours = 0.0;
		minutes = 0.0;
		seconds = 0.0;
		cstr = strchr (nval,':');
		if (cstr > string) {
		    *cstr = '\0';
		    hours = (int) atof (nval);
		    nval = cstr + 1;
		    cstr = strchr (nval,':');
		    if (cstr > string) {
			minutes = (int) atof (nval);
			nval = cstr + 1;
			cstr = strchr (nval,':');
			if (cstr > string)
			    seconds = atof (nval);
			}
		    }
		fday = ((3.6e3 * (double)hours) + (6.e1 * (double)minutes) +
		       seconds) / 8.64e4;
		epoch = epoch + (fday / yeardays);
		}
	    return (epoch);
	    }
	else if (isnum (string))
	    return (atof (string));
	else
	    return (0.0);
	}
    else
	return (0.0);
}


/* DT2TS-- convert from date, time as yyyy.mmdd hh.mmsss to sec since 1950.0 */

double
dt2ts (date,time)

double	date;	/* Date as yyyy.mmdd
		    yyyy = calendar year (e.g. 1973)
		    mm = calendar month (e.g. 04 = april)
		    dd = calendar day (e.g. 15) */
double	time;	/* Time as hh.mmssxxxx
		    *if time<0, it is time as -(fraction of a day)
		    hh = hour of day (0 .le. hh .le. 23)
		    nn = minutes (0 .le. nn .le. 59)
		    ss = seconds (0 .le. ss .le. 59)
		  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999) */
{
    double tsec; /* Seconds past 1950.0 (returned) */

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

double	tsec;	/* Seconds past 1950.0 */
double	*date;	/* Date as yyyy.mmdd (returned)
		    yyyy = calendar year (e.g. 1973)
		    mm = calendar month (e.g. 04 = april)
		    dd = calendar day (e.g. 15) */
double	*time;	/* Time as hh.mmssxxxx (returned)
		    *if time<0, it is time as -(fraction of a day)
		    hh = hour of day (0 .le. hh .le. 23)
		    nn = minutes (0 .le. nn .le. 59)
		    ss = seconds (0 .le. ss .le. 59)
		  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999) */
{
    int iyr,imon,iday,ihr,imn;
    double sec;

    ts2i (tsec,&iyr,&imon,&iday,&ihr,&imn,&sec, 4);

    /* Convert date to yyyy.mmdd */
    *date = (double) iyr + 0.01 * (double) imon + 0.0001 * (double) iday;

    /* Convert time to hh.mmssssss */
    *time = (double) ihr + 0.01 * (double) imn + 0.0001 * sec;

    return;
}


/* TS2FD-- convert seconds since 1950.0 to FITS date, yyyy-mm-ddThh:mm:ss.ss */

char *
ts2fd (tsec)

double	tsec;	/* Seconds past 1950.0 */
{
    int iyr,imon,iday,ihr,imn;
    double sec;
    char *string;

    ts2i (tsec,&iyr,&imon,&iday,&ihr,&imn,&sec, 3);

    /* Convert to ISO date format */
    string = (char *) calloc (1, 32);
    sprintf (string, "%4d-%02d-%02dT%02d:%02d:%06.3f",
	     iyr, imon, iday, ihr, imn, sec);

    return (string);
}


/* TS2I-- convert sec since 1950.0 to year month day hours minutes seconds */

void
ts2i (tsec,iyr,imon,iday,ihr,imn,sec, ndsec)

double	tsec;	/* seconds since 1/1/1950 0:00 */
int	*iyr;	/* year (returned) */
int	*imon;	/* month (returned) */
int	*iday;	/* day (returned) */
int	*ihr;	/* hours (returned) */
int	*imn;	/* minutes (returned) */
double	*sec;	/* seconds (returned) */
int	ndsec;	/* Number of decimal places in seconds (0=int) */

{
    double t,days;
    int isec,ihms,nc,nc4,nly,ny,m,im;

    /* Round seconds to 0 - 4 decimal places */
    if (ndsec < 1)
	t = dint (tsec + 61530883200.5) * 10000.0;
    else if (ndsec < 2)
	t = dint ((tsec + 61530883200.0) * 10.0 + 0.5) * 1000.0;
    else if (ndsec < 3)
	t = dint ((tsec + 61530883200.0) * 100.0 + 0.5) * 100.0;
    else if (ndsec < 3)
	t = dint ((tsec + 61530883200.0) * 1000.0 + 0.5) * 10.0;
    else
	t = dint ((tsec + 61530883200.0) * 10000.0 + 0.5);

    /* Time of day (hours, minutes, seconds, .1 msec) */
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

/* Jul  1 1999	New file, based on iolib/jcon.f and iolib/vcon.f and hgetdate()
 * Oct 21 1999	Fix declarations after lint
 * Oct 27 1999	Fix bug to return epoch if fractional year input
 * Dec  6 1999	Fix bug in ts2jd() found by Pete Ratzlaff (SAO)
 */
