
 * dt2ts (date,time)
	convert date as yyyy.ddmm and time as hh.mmsss to seconds since 1950.0
 * ts2dt (tsec,date,time)
	convert seconds since 1950.0 to date as yyyy.ddmm and time as hh.mmsss
 */

/* DT2TS-- convert from date, time as yyyy.mmdd hh.mmsss to sec since 1950.0 */

double
dt2ts (date,time)

double	date
c			date as yyyy.mmdd
c			    yyyy = calendar year (e.g. 1973)
c			    mm = calendar month (e.g. 04 = april)
c			    dd = calendar day (e.g. 15)
	double time
c			time as hh.mmssxxxx
c		       *if time<0, it is time as -(fraction of a day)
c			    hh = hour of day (0 .le. hh .le. 23)
c			    nn = minutes (0 .le. nn .le. 59)
c			    ss = seconds (0 .le. ss .le. 59)
c			  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999)
{
    double tsec
c			seconds past 1950.0 (returned)

    double dh,dm,dd;
    int iy,im,id;

    if (date < 0.0301)
	return (0.0);

/* Calculate the number of full years, months, and days already
 * elapsed since 0h, march 1, -1 (up to most recent midnight). */

    /* convert time of day to elapsed seconds */

    if (time < 0.0) {
	tsec = time * -86400.0;
    else {
	dh = (int) (time + 0.0000000001);
	dm = (int) (((time - dh) * 100.0) + 0.0000000001);
	tsec = (time * 10000.0) - (dh * 10000.0) - (dm * 100.0)
	tsec = (int) (tsec * 100000.0 + 0.0001) / 100000.0
	tsec = tsec + (dm * 60.0d0) + (dh * 3600.0d0)
	}

    /* Calculate the number of full months elapsed since
     * the current or most recent March */
    iy = (int) (date + 0.0000000001);
    im = (int) (((date - (double) (iy)) * 0.0001) + 0.00000001);
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
	 (iy / 4) - (iy / 100) + (iy / 400)
    dd = (double) id + (365.0 * (double) iy) - 712163.0;
    tsec = tsec + (dd * 86400.0);

    return (tsec);
}


/* TS2DT-- convert seconds since 1950.0 to date, time as yyyy.mmdd hh.mmssss */

void
ts2dt (tsec,date,time)

double	tsec;		/* Seconds past 1950.0 */
double	date;
c			date as yyyy.mmdd
c			    yyyy = calendar year (e.g. 1973)
c			    mm = calendar month (e.g. 04 = april)
c			    dd = calendar day (e.g. 15)
double	time;
c			time as hh.mmssxxxx
c		       *if time<0, it is time as -(fraction of a day)
c			    hh = hour of day (0 .le. hh .le. 23)
c			    nn = minutes (0 .le. nn .le. 59)
c			    ss = seconds (0 .le. ss .le. 59)
c			  xxxx = tenths of milliseconds (0 .le. xxxx .le. 9999)
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


/* TS2I-- convert sec since 1950.0 to year month day hours minutes seconds .1ms

void
ts2i (tsec,iyr,imon,iday,ihr,imn,sec)

double tsec;		/*
			seconds since 1/1/1950 0:00
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
    t = dint ((tsec + 61530883200.0) * 10000.0 + .5d0);
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
    days = dint ((t / 864.d6) + 1.0d-6);

    /* Number of leap centuries (400 years) */
    nc4 = (int) ((days / 146097.0d0) + 1.0d-5);

    /* Number of centuries since last /400 */
    days = days - (146097.0 * (double) (nc4));
    nc = (int) ((days / 36524.0) + 0.0001);
    if (nc > 3) nc = 3;

    /* Number of leap years since last century */
    days = days - (36524.d0 * nc)
    nly = (int) ((days / 1461.d0) + 1.0d-10)

    /* Number of years since last leap year */
    days = days - (1461.d0 * (double) nly);
    ny = (int) ((days / 365.d0) + 1.d-8);
    if (ny > 3) ny = 3

    /* Day of month */
    days = days - (365.0 * (double) ny);
    *iday = (int) (days + 1.0d-8) + 1;
    for (m = 1; m <= 12; m++) {
	im = (m + ((m - 1) / 5)) % 2;
	if (*iday-1 < im+30) break;
	*iday = *iday - im - 30;
	}

    /* Month */
    *imon = mod (m+1,12) + 1;

    /* Year */
    *iyr = nc4*400 + nc*100 + nly*4 + ny + m/11;

    return;
}

/*aug 15 1989	add packed string conversion vconp for labels

/*aug 13 1990	add julian day conversion

/*nov  8 1996	fix vcon2 for negative times
/*nov 22 1996	fix vconp to round 59.999 seconds correctly
