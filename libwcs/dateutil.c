/* File libwcs/dateutil.c
 * December 20, 1999
 * By Doug Mink
 */

/* Date and time conversion routines using the following conventions:
   dt = 2 floating point numbers: yyyy.mmdd, hh.mmssssss
   jd = Julian Date
   ts = seconds since 1950.0 (used for ephemeris computations
   ep = fractional year, often epoch of a position including proper motion
   fd = FITS date string (dd/mm/yy before 2000, then yyyy-mm-dd[Thh:mm:ss.ss])
 time = use fd2* with no date to convert time as hh:mm:ss.ss to sec, day, year

 * dt2ep (date, time)
 *	Convert date as yyyy.ddmm and time as hh.mmsss to fractional year
 * dt2fd (date, time)
 *	Convert date as yyyy.ddmm and time as hh.mmsss to FITS date string
 * dt2i (date,time,iyr,imon,iday,ihr,imn,sec, ndsec)
 *	Convert yyyy.mmdd hh.mmssss to year month day hours minutes seconds
 * dt2jd (date,time)
 *	Convert date as yyyy.ddmm and time as hh.mmsss to Julian date
 * dt2ts (date,time)
 *	Convert date as yyyy.ddmm and time as hh.mmsss to seconds since 1950.0
 * ep2dt (epoch,date, time)
 *	Convert fractional year to date as yyyy.ddmm and time as hh.mmsss
 * ep2fd (epoch)
 *	Convert epoch to FITS ISO date string
 * ep2i (epoch,iyr,imon,iday,ihr,imn,sec, ndsec)
 *	Convert fractional year to year month day hours minutes seconds
 * ep2jd (epoch)
 *	Convert fractional year as used in epoch to Julian date
 * ep2ts (epoch)
 *	Convert fractional year to seconds since 1950.0
 * fd2ep (string)
 *	Convert FITS date string to fractional year
 *	Convert time alone to fraction of Besselian year
 * fd2dt (string, date, time)
 *	Convert FITS date string to date as yyyy.ddmm and time as hh.mmsss
 *	Convert time alone to hh.mmssss with date set to 0.0
 * fd2i (string,iyr,imon,iday,ihr,imn,sec, ndsec)
 *	Convert FITS standard date string to year month day hours min sec
 *	Convert time alone to hours min sec, year month day are zero
 * fd2jd (string)
 *	Convert FITS standard date string to Julian date
 *	Convert time alone to fraction of day
 * fd2ts (string)
 *	Convert FITS standard date string to seconds since 1950.0
 *	Convert time alone to seconds of day
 * fd2fd (string)
 *	Convert FITS standard date string to ISO FITS date string
 * jd2dt (dj,date,time)
 *	Convert Julian date to date as yyyy.mmdd and time as hh.mmssss
 * jd2ep (dj)
 *	Convert Julian date to fractional year as used in epoch
 * jd2fd (dj)
 *	Convert Julian date to FITS ISO date string
 * jd2i (dj,iyr,imon,iday,ihr,imn,sec, ndsec)
 *	Convert Julian date to year month day hours min sec
 * jd2ts (dj)
 *	Convert Julian day to seconds since 1950.0
 * ts2dt (tsec,date,time)
 *	Convert seconds since 1950.0 to date as yyyy.ddmm and time as hh.mmsss
 * ts2ep (tsec)
 *	Convert seconds since 1950.0 to fractional year
 * ts2fd (tsec)
 *	Convert seconds since 1950.0 to FITS standard date string
 * ts2i (tsec,iyr,imon,iday,ihr,imn,sec, ndsec)
 *	Convert sec since 1950.0 to year month day hours minutes seconds
 * ts2jd (tsec)
 *	Convert seconds since 1950.0 to Julian date
 * isdate (string)
 *	Return 1 if string is a FITS date (old or ISO)
 *
 * Internally-used subroutines
 *
 * fixdate (iyr, imon, iday, ihr, imn, sec, ndsec)
 *	Round seconds and make sure date and time numbers are within limits
 * caldays (year, month)
 *	Calculate days in month 1-12 given year (Gregorian calendar only
 * dint (dnum)
 *	Return integer part of floating point number
 * dmod (dnum)
 *	Return Mod of floating point number
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "fitsfile.h"

static void fixdate();
static int caldays();
static double dint();
static double dmod();

/* DT2FD-- convert vigesimal date and time to FITS date, yyyy-mm-ddThh:mm:ss.ss */

char *
dt2fd (date, time)

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
    int iyr,imon,iday,ihr,imn;
    double sec;
    char *string;

    dt2i (date, time, &iyr,&imon,&iday,&ihr,&imn,&sec, 3);

    /* Convert to ISO date format */
    string = (char *) calloc (1, 32);
    if (date == 0.0)
	sprintf (string, "%02d:%02d:%06.3f", ihr, imn, sec);
    else if (time == 0.0)
	sprintf (string, "%4d-%02d-%02d", iyr, imon, iday);
    else
	sprintf (string, "%4d-%02d-%02dT%02d:%02d:%06.3f",
		 iyr, imon, iday, ihr, imn, sec);

    return (string);
}

/* DT2JD-- convert from date as yyyy.mmdd and time as hh.mmsss to Julian Date
 *	   Return fractional days if date is zero */

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

    tsec = dt2ts (date,time);
    if (date == 0.0)
	dj = tsec / 86400.0;
    else
	dj = ts2jd (tsec);

    return (dj);
}


/* JD2DT-- convert Julian date to date as yyyy.mmdd and time as hh.mmssss */

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

    tsec = jd2ts (dj);
    ts2dt (tsec, date, time);

    return;
}


/* JD2I-- convert Julian date to date as yyyy.mmdd and time as hh.mmssss */

void
jd2i (dj, iyr, imon, iday, ihr, imn, sec, ndsec)

double	dj;	/* Julian date */
int	*iyr;	/* year (returned) */
int	*imon;	/* month (returned) */
int	*iday;	/* day (returned) */
int	*ihr;	/* hours (returned) */
int	*imn;	/* minutes (returned) */
double	*sec;	/* seconds (returned) */
int	ndsec;	/* Number of decimal places in seconds (0=int) */

{
    double tsec;

    tsec = jd2ts (dj);
    ts2i (tsec, iyr, imon, iday, ihr, imn, sec, ndsec);
    return;
}


/* JD2EP-- convert Julian date to fractional year as used in epoch */

double
jd2ep (dj)

double	dj;	/* Julian date */

{
    double date, time;

    jd2dt (dj, &date, &time);
    return (dt2ep (date, time));
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


/* EP2I-- convert fractional year to year month day hours min sec */

void
ep2i (epoch, iyr, imon, iday, ihr, imn, sec, ndsec)

double	epoch;	/* Date as fractional year */
int	*iyr;	/* year (returned) */
int	*imon;	/* month (returned) */
int	*iday;	/* day (returned) */
int	*ihr;	/* hours (returned) */
int	*imn;	/* minutes (returned) */
double	*sec;	/* seconds (returned) */
int	ndsec;	/* Number of decimal places in seconds (0=int) */
{
    double date, time;

    ep2dt (epoch, &date, &time);
    dt2i (date, time, iyr,imon,iday,ihr,imn,sec, ndsec);
    return;
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


/* TS2EP-- convert seconds since 1950.0 to fractional year as used in epoch */

double
ts2ep (tsec)

double	tsec;	/* Seconds since 1950.0 */

{
    double date, time;

    ts2dt (tsec, &date, &time);
    return (dt2ep (date, time));
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
    if (date == 0.0) {
	epoch = dj / 365.2522;
	return (epoch);
	}
    else {
	time0 = 12.0;
	date0 = dint (date) + 0.0101;
	date1 = dint (date) + 1.0101;
	dj0 = dt2jd (date0, time0);
	dj1 = dt2jd (date1, time0);
	epoch = dint (date) + ((dj - dj0) / (dj1 - dj0));
	return (epoch);
	}
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
    double date, time;

    fd2dt (string, &date, &time);
    return (dt2jd (date, time));
}


/* FD2TS-- convert FITS standard date to seconds since 1950 */

double
fd2ts (string)

char *string;	/* FITS date string, which may be:
			fractional year
			dd/mm/yy (FITS standard before 2000)
			dd-mm-yy (nonstandard use before 2000)
			yyyy-mm-dd (FITS standard after 1999)
			yyyy-mm-ddThh:mm:ss.ss (FITS standard after 1999) */


{
    double date, time;

    fd2dt (string, &date, &time);
    return (dt2ts (date, time));
}


/* FD2FD-- convert any FITS standard date to ISO FITS standard date */

char *
fd2fd (string)

char *string;	/* FITS date string, which may be:
			fractional year
			dd/mm/yy (FITS standard before 2000)
			dd-mm-yy (nonstandard use before 2000)
			yyyy-mm-dd (FITS standard after 1999)
			yyyy-mm-ddThh:mm:ss.ss (FITS standard after 1999) */


{
    double date, time;

    fd2dt (string, &date, &time);
    return (dt2fd (date, time));
}


/* FD2DT-- convert FITS standard date to date, time as yyyy.mmdd hh.mmsss */

void
fd2dt (string, date, time)

char *string;	/* FITS date string, which may be:
		    fractional year
		    dd/mm/yy (FITS standard before 2000)
		    dd-mm-yy (nonstandard use before 2000)
		    yyyy-mm-dd (FITS standard after 1999)
		    yyyy-mm-ddThh:mm:ss.ss (FITS standard after 1999) */
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

    fd2i (string,&iyr,&imon,&iday,&ihr,&imn,&sec, 4);

    /* Convert date to yyyy.mmdd */
    *date = (double) iyr + 0.01 * (double) imon + 0.0001 * (double) iday;

    /* Convert time to hh.mmssssss */
    *time = (double) ihr + 0.01 * (double) imn + 0.0001 * sec;

    return;
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
    double date, time;

    fd2dt (string, &date, &time);
    return (dt2ep (date, time));
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

/* Calculate the number of full years, months, and days already
 * elapsed since 0h, March 1, -1 (up to most recent midnight). */

    /* convert time of day to elapsed seconds */

    /* If time is < 0, it is assumed to be a fractional day */
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
    if (date >= 0.0301) {
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
	}

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
    double date, time;
    char *string;

    ts2dt (tsec, &date, &time);
    return (dt2fd (date, time));
}


/* DT2I-- convert vigesimal date and time to year month day hours min sec */

void
dt2i (date, time, iyr, imon, iday, ihr, imn, sec, ndsec)

double	date;	/* Date as yyyy.mmdd (returned)
		    yyyy = calendar year (e.g. 1973)
		    mm = calendar month (e.g. 04 = april)
		    dd = calendar day (e.g. 15) */
double	time;	/* Time as hh.mmssxxxx (returned)
		    *if time<0, it is time as -(fraction of a day)
		    hh = hour of day (0 .le. hh .le. 23)
		    nn = minutes (0 .le. nn .le. 59)
		    ss = seconds (0 .le. ss .le. 59)
		  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999) */
int	*iyr;	/* year (returned) */
int	*imon;	/* month (returned) */
int	*iday;	/* day (returned) */
int	*ihr;	/* hours (returned) */
int	*imn;	/* minutes (returned) */
double	*sec;	/* seconds (returned) */
int	ndsec;	/* Number of decimal places in seconds (0=int) */

{
    double t,d;
    int days;

    t = time;
    d = date;

    /* Extract components of time */
    *ihr = dint (t + 0.000000001);
    t = 100.0 * (t - (double) *ihr);
    *imn = dint (t + 0.0000001);
    *sec = 100.0 * (t - (double) *imn);

    /* Extract components of date */
    *iyr = dint (d + 0.00001);
    d = 100.0 * (d - (double) *iyr);
    *imon = dint (d + 0.001);
    d = 100.0 * (d - (double) *imon);
    *iday = dint (d + 0.1);

   /* Make sure date and time are legal */
    fixdate (iyr, imon, iday, ihr, imn, sec, ndsec);

    return;
}


/* FD2I -- convert from FITS standard date to year, mon, day, hours, min, sec */

void
fd2i (string, iyr, imon, iday, ihr, imn, sec, ndsec)

char	*string; /* FITS date string, which may be:
			yyyy.ffff (fractional year)
			dd/mm/yy (FITS standard before 2000)
			dd-mm-yy (nonstandard FITS use before 2000)
			yyyy-mm-dd (FITS standard after 1999)
			yyyy-mm-ddThh:mm:ss.ss (FITS standard after 1999) */
int	*iyr;	/* year (returned) */
int	*imon;	/* month (returned) */
int	*iday;	/* day (returned) */
int	*ihr;	/* hours (returned) */
int	*imn;	/* minutes (returned) */
double	*sec;	/* seconds (returned) */
int	ndsec;	/* Number of decimal places in seconds (0=int) */

{
    double tsec;
    int days, i;
    char *sstr, *dstr, *tstr, *cstr, *nval;

    /* Initialize all returned data to zero */
    *iyr = 0;
    *imon = 0;
    *iday = 0;
    *ihr = 0;
    *imn = 0;
    *sec = 0,0;

    /* Return if no input string */
    if (string == NULL)
	return;

    /* Check for various non-numeric characters */
    sstr = strchr (string,'/');
    dstr = strchr (string,'-');
    tstr = strchr (string,'T');
    cstr = strchr (string,':');

    /* Original FITS date format: dd/mm/yy */
    if (sstr > string) {
	*sstr = '\0';
	*iday = (int) atof (string);
	*sstr = '/';
	nval = sstr + 1;
	sstr = strchr (nval,'/');
	if (sstr == NULL)
	    sstr = strchr (nval,'-');
	if (sstr > string) {
	    *sstr = '\0';
	    *imon = (int) atof (nval);
	    *sstr = '/';
	    nval = sstr + 1;
	    *iyr = (int) atof (nval);
	    if (*iyr >= 0 && *iyr <= 49)
		*iyr = *iyr + 2000;
	    else if (*iyr < 100)
		*iyr = *iyr + 1900;
	    }
	else
	    return;
	}

    /* New FITS date format: yyyy-mm-ddThh:mm:ss[.sss] */
    else if (dstr > string) {
	*dstr = '\0';
	*iyr = (int) atof (string);
	*dstr = '-';
	nval = dstr + 1;
	dstr = strchr (nval,'-');
	*imon = 1;
	*iday = 1;

	/* Decode year, month, and day */
	if (dstr > string) {
	    *dstr = '\0';
	    *imon = (int) atof (nval);
	    *dstr = '-';
	    nval = dstr + 1;
	    if (tstr > string)
		*tstr = '\0';
	    *iday = (int) atof (nval);
	    if (tstr > string)
		*tstr = 'T';
	    }

	/* If year is < 32, it is really day of month in old format */
	if (*iyr < 32 || *iday > 31) {
	    i = *iyr;
	    if (*iday < 100)
		*iyr = *iday + 1900;
	    else
		*iyr = *iday;
	    *iday = i;
	    }
	}

    /* In rare cases, a FITS time is entered as an epoch */
    else if (tstr == NULL && cstr == NULL && isnum (string)) {
	tsec = ep2ts (atof (string));
	ts2i (tsec,iyr,imon,iday,ihr,imn,sec, ndsec);
	return;
	}

    /* Extract time, if it is present */
    if (tstr > string || cstr > string) {
	if (tstr > string)
	    nval = tstr + 1;
	else
	    nval = string;
	cstr = strchr (nval,':');
	if (cstr > string) {
	    *cstr = '\0';
	    *ihr = (int) atof (nval);
	    *cstr = ':';
	    nval = cstr + 1;
	    cstr = strchr (nval,':');
	    if (cstr > string) {
		*cstr = '\0';
		*imn = (int) atof (nval);
		*cstr = ':';
		nval = cstr + 1;
		*sec = atof (nval);
		}
	    else
		*imn = (int) atof (nval);
	    }
	else
	    *ihr = (int) atof (nval);
	}
    else
	ndsec = -1;

   /* Make sure date and time are legal */
    fixdate (iyr, imon, iday, ihr, imn, sec, ndsec);

    return;

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

   /* Make sure date and time are legal */
    fixdate (iyr, imon, iday, ihr, imn, sec, ndsec);

    return;
}


/* ISDATE - Return 1 if string is an old or ISO FITS standard date */

int
isdate (string)

char	*string; /* Possible FITS date string, which may be:
			dd/mm/yy (FITS standard before 2000)
			dd-mm-yy (nonstandard FITS use before 2000)
			yyyy-mm-dd (FITS standard after 1999)
			yyyy-mm-ddThh:mm:ss.ss (FITS standard after 1999) */

{
    int iyr;	/* year (returned) */
    int imon;	/* month (returned) */
    int iday;	/* day (returned) */
    int i;
    char *sstr, *dstr, *tstr, *nval;

    /* Translate string from ASCII to binary */
    if (string == NULL) 
	return (0);

    sstr = strchr (string,'/');
    dstr = strchr (string,'-');
    tstr = strchr (string,'T');

    /* Original FITS date format: dd/mm/yy */
    if (sstr > string) {
	*sstr = '\0';
	iday = (int) atof (string);
	*sstr = '/';
	nval = sstr + 1;
	sstr = strchr (nval,'/');
	if (sstr == NULL)
	    sstr = strchr (nval,'-');
	if (sstr > string) {
	    *sstr = '\0';
	    imon = (int) atof (nval);
	    *sstr = '/';
	    nval = sstr + 1;
	    iyr = (int) atof (nval);
	    if (iyr >= 0 && iyr <= 49)
		iyr = iyr + 2000;
	    else if (iyr < 100)
		iyr = iyr + 1900;
	    }
	if (imon > 0 && iday > 0)
	    return (1);
	else
	    return (0);
	}

    /* New FITS date format: yyyy-mm-ddThh:mm:ss[.sss] */
    else if (dstr > string) {
	*dstr = '\0';
	iyr = (int) atof (string);
	nval = dstr + 1;
	*dstr = '-';
	dstr = strchr (nval,'-');
	imon = 0;
	iday = 0;

	/* Decode year, month, and day */
	if (dstr > string) {
	    *dstr = '\0';
	    imon = (int) atof (nval);
	    *dstr = '-';
	    nval = dstr + 1;
	    if (tstr > string)
		*tstr = '\0';
	    iday = (int) atof (nval);
	    if (tstr > string)
		*tstr = 'T';
	    }

	/* If year is < 32, it is really day of month in old format */
	if (iyr < 32 || iday > 31) {
	    i = iyr;
	    if (iday < 100)
		iyr = iday + 1900;
	    else
		iyr = iday;
	    iday = i;
	    }
	if (imon > 0 && iday > 0)
	    return (1);
	else
	    return (0);
	}

    /* If FITS date is entered as an epoch, return 0 anyway */
    else
	return (0);
}


/* Round seconds and make sure date and time numbers are within limits */

static void
fixdate (iyr, imon, iday, ihr, imn, sec, ndsec)

int	*iyr;	/* year (returned) */
int	*imon;	/* month (returned) */
int	*iday;	/* day (returned) */
int	*ihr;	/* hours (returned) */
int	*imn;	/* minutes (returned) */
double	*sec;	/* seconds (returned) */
int	ndsec;	/* Number of decimal places in seconds (0=int) */
{
    double days;

    /* Round seconds to 0 - 4 decimal places (no rounding if <0, >4) */
    if (ndsec == 0)
	*sec = dint (*sec + 0.5);
    else if (ndsec < 2)
	*sec = dint ((*sec * 10.0 + 0.5) / 10.0);
    else if (ndsec < 3)
	*sec = dint ((*sec * 100.0 + 0.5) / 100.0);
    else if (ndsec < 4)
	*sec = dint ((*sec * 1000.0 + 0.5) / 1000.0);
    else if (ndsec < 5)
	*sec = dint ((*sec * 10000.0 + 0.5) / 10000.0);

    /* Adjust minutes and hours */
    if (*sec > 60.0) {
	*sec = *sec - 60.0;
	*imn = *imn + 1;
	}
    if (*imn > 60) {
	*imn = *imn - 60;
	*ihr = *ihr + 1;
	}

    /* Return if no date */
    if (*iyr == 0 && *imon == 0 && *iday == 0)
	return;

   /* Adjust date */
    if (*ihr > 23) {
	*ihr = *ihr - 24;
	*iday = *iday + 1;
	}
    days = caldays (*iyr, *imon);
    if (*iday > days) {
	*iday = *iday - days;
	*imon = *imon + 1;
	}
    if (*iday < 1) {
	*imon = *imon - 1;
	if (*imon < 1) {
	    *imon = *imon + 12;
	    *iyr = *iyr - 1;
	    }
	days = caldays (*iyr, *imon);
	*iday = *iday + days;
	}
    if (*imon < 1) {
	*imon = *imon + 12;
	*iyr = *iyr - 1;
	days = caldays (*iyr, *imon);
	if (*iday > days) {
	    *iday = *iday - days;
	    *imon = *imon + 1;
	    }
	}
    if (*imon > 12) {
	*imon = *imon - 12;
	*iyr = *iyr + 1;
	}
    return;
}


/* Calculate days in month 1-12 given year (Gregorian calendar only */

static int
caldays (year, month)

int	year;	/* 4-digit year */
int	month;	/* Month (1=January, 2=February, etc.) */
{
    if (month < 1) {
	month = month + 12;
	year = year + 1;
	}
    if (month > 12) {
	month = month - 12;
	year = year + 1;
	}
    switch (month) {
	case 1:
	    return (31);
	case 2:
	    if (year%400 == 0)
		return (29);
	    else if (year%100 == 0)
		return (28);
	    else if (year%4 == 0)
		return (29);
	    else
		return (28);
	case 3:
	    return (31);
	case 4:
	    return (30);
	case 5:
	    return (31);
	case 6:
	    return (30);
	case 7:
	    return (31);
	case 8:
	    return (31);
	case 9:
	    return (30);
	case 10:
	    return (31);
	case 11:
	    return (30);
	case 12:
	    return (31);
	default:
	    return (0);
	}
}


static double
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


static double
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
 * Dec  9 1999	Fix bug in ts2jd() found by Pete Ratzlaff (SAO)
 * Dec 17 1999	Add all unimplemented conversions
 * Dec 20 1999	Add isdate(); leave date, time strings unchanged in fd2i()
 * Dec 20 1999	Make all fd2*() subroutines deal with time alone
 */
