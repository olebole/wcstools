/* File getdate.c
 * August 30, 2002
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "libwcs/wcs.h"
#include "libwcs/fitsfile.h"

#define DTVIG	1	/* Vigesimal (yyyy.mmdd hh.mmssss) */
#define DTFITS	2	/* FITS date (new format) */
#define DTJD	3	/* Julian date */
#define DTMJD	4	/* Modified Julian date */
#define DT1950	5	/* Seconds since 1950.0 */
#define DTEP	6	/* Bessellian Epoch */
#define DTEPB	7	/* Bessellian Epoch */
#define DTEPJ	8	/* Julian Epoch */
#define DTOF	9	/* Old FITS date and time */
#define DTOFD	10	/* Old FITS date only */
#define DTOFT	11	/* Old FITS time only */
#define DTLT	12	/* Current local time */
#define DTUT	13	/* Current Universal Time */
#define DTIRAF	14	/* IRAF seconds since 1980-01-01 0:00 */
#define DTUNIX	15	/* Unix seconds since 1970-01-01 0:00 */
#define DTDOY	16	/* Year and day of year (including fraction) */

static void usage();
static void ConvertDate();

static int verbose = 0;		/* Verbose/debugging flag */
static int wfile = 0;		/* True to print output file */
static int debug = 0;		/* True for extra information */
static int sumcol = 0;		/* True to sum columns */
static int meancol = 0;		/* True to compute mean of columns */
static int countcol = 0;	/* True to count entries in columns */
static char *objname = NULL;	/* Object name for output */
static char *keyword = NULL;	/* Column to add to tab table output */
static char progpath[128];
static int ncat = 0;
static int et = 0;	/* If 1, convert output to ET */
static int version = 0;	/* If 1, print only program name and version */
static int nread = 0;	/* Number of lines to read (0=all) */
static int nskip = 0;	/* Number of lines to skip */
static int tabout = 0;	/* If 1, separate output fields with tabs */
static int ndec = 5;	/* Number of decimal places in output */
static int dateonly = 0; /* If 1, print date without time */
static char *outform0;
static int nftok = 2;	/* Number of tokens for FITS date format */
extern void sdatedec();

main (ac, av)
int ac;
char **av;
{
    char *str;
    char *datestring;
    char *timestring;
    int intype = 0;	/* Input date type */
    int outtype = 0;	/* Output date type */
    int appdate = 0;	/* Append date to input file */
    int typeset = 0;
    int lline;
    char line[82];
    FILE *fd;

    datestring = NULL;
    timestring = NULL;
    outform0 = NULL;

    if (ac == 1)
        usage ();

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage ();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage ();
	}

    /* crack arguments */
    for (av++; --ac > 0; av++) {

	/* Set input date format */
	if (!strncmp (*av, "dt2", 3))
	    intype = DTVIG;
	else if (!strncmp (*av, "ep2", 3))
	    intype = DTEP;
	else if (!strncmp (*av, "epb2", 3))
	    intype = DTEPB;
	else if (!strncmp (*av, "epj2", 3))
	    intype = DTEPJ;
	else if (!strncmp (*av, "nfd2", 3)) {
	    intype = DTFITS;
	    nftok = 1;
	    }
	else if (!strncmp (*av, "fd2", 3))
	    intype = DTFITS;
	else if (!strncmp (*av, "jd2", 3))
	    intype = DTJD;
	else if (!strncmp (*av, "mjd2", 4))
	    intype = DTMJD;
	else if (!strncmp (*av, "ts2", 3))
	    intype = DT1950;
	else if (!strncmp (*av, "tsi2", 3))
	    intype = DTIRAF;
	else if (!strncmp (*av, "tsu2", 3))
	    intype = DTUNIX;
	else if (!strncmp (*av, "lt2", 3))
	    intype = DTLT;
	else if (!strncmp (*av, "ut2", 3))
	    intype = DTUT;
	else if (!strncmp (*av, "doy", 3))
	    intype = DTDOY;

	/* Set output date format */
	if (strsrch (*av, "2dt"))
	    outtype = DTVIG;
	else if (strsrch (*av, "2epb"))
	    outtype = DTEPB;
	else if (strsrch (*av, "2epj"))
	    outtype = DTEPJ;
	else if (strsrch (*av, "2ep"))
	    outtype = DTEP;
	else if (strsrch (*av, "2fd"))
	    outtype = DTFITS;
	else if (strsrch (*av, "2jd"))
	    outtype = DTJD;
	else if (strsrch (*av, "2mjd"))
	    outtype = DTMJD;
	else if (strsrch (*av, "2of"))
	    outtype = DTOF;
	else if (strsrch (*av, "2ofd"))
	    outtype = DTOFD;
	else if (strsrch (*av, "2oft"))
	    outtype = DTOFT;
	else if (strsrch (*av, "2tsi"))
	    outtype = DTIRAF;
	else if (strsrch (*av, "2tsu"))
	    outtype = DTUNIX;
	else if (strsrch (*av, "2ts"))
	    outtype = DT1950;
	else if (strsrch (*av, "2doy"))
	    outtype = DTDOY;

	if (typeset == 0 && (intype > 0 || outtype > 0))
	    typeset = 1;
	else if (*(str = *av) == '-' && !isdate(*av) && !isnum(*av)) {
	    char c;
	    while (c = *++str)
	    switch (c) {

		case 'a':	/* Append date to input file */
		    appdate++;
		    break;

		case 'd':	/* Print date without time */
		    dateonly++;
		    break;

		case 'e':	/* Print Ephemeris Time (ET, TDT, TT) */
		    et++;
		    break;

		case 'f':	/* Output format */
		    if (ac < 2)
			usage();
		    outform0 = *++av;
		    ac--;
		    break;

		case 'n':	/* Number of decimal places in output */
		    if (ac < 2)
			usage();
		    ndec = (int) (atof (*++av));
		    setdatedec (ndec);
		    ac--;
		    break;

		case 'v':	/* More verbosity */
		    verbose++;
		    break;

		default:
		    usage ();
		    break;
		}
	    }
	else if (*(str = *av) == '@') {
	    if (strlen (str) > 1) {
		if (!strcmp (str+1, "stdin"))
		    fd = stdin;
		else if (!(fd = fopen (str+1, "r"))) {
		    if (verbose)
			fprintf (stderr, "*** Cannot read file %s\n", str+1);
		    return (-1);
		    }
		datestring = (char *) calloc (80,1);
		timestring = (char *) calloc (80,1);
		while (fgets (line, 80, fd)) {
		    lline = strlen (line);
		    if ((int)line[lline-1] < 32)
			line[lline-1] = (char) 0;
		    if (appdate)
			printf ("%s ", line);
		    if (intype != DTVIG && intype != DTDOY && nftok != 1)
			sscanf (line, "%s %s", datestring, timestring);
		    else
			sscanf (line, "%s", datestring);
		    ConvertDate (intype, outtype, datestring, timestring);
		    }
		free (datestring);
		free (timestring);
		}
	    }
		
	else if (datestring == NULL) {
	    datestring = *av;
	    if ((intype != DTVIG && intype != DTDOY) || ac == 1) {
		ConvertDate (intype, outtype, datestring, timestring);
		datestring = NULL;
		}
	    }
	else if (intype == DTVIG || intype == DTDOY) {
	    timestring = *av;
	    ConvertDate (intype, outtype, datestring, timestring);
	    timestring = NULL;
	    }
	}
    if (intype == DTLT || intype == DTUT)
	ConvertDate (intype, outtype, datestring, timestring);

    return (0);
}

static void
usage ()

{
    if (version)
	exit (-1);
    fprintf (stderr,"Convert date and time between various formats\n");
    fprintf (stderr,"Usage: [-dv][-n dec][-f format] itype2otype [date and/or time]\n");
    fprintf (stderr,"       [-dv][-n dec][-f format] itype2otype @file\n");
    fprintf(stderr,"  itype: nfd=ISOFITS fd=FITS, dt=yyyy.mmdd\n");
    fprintf(stderr,"         jd=Julian Date, mjd=Modified Julian Date\n");
    fprintf(stderr,"         ep=epoch, epj=Julian epoch, epb=Besselian epoch\n");
    fprintf(stderr,"         lt=local time, ut=UT, ts=seconds since 1950-01-01\n");
    fprintf(stderr,"  otype: fd=FITS, dt=yyyy.mmdd, jd=Julian Date, mjd=Modified Julian Date\n");
    fprintf(stderr,"         ep=epoch, epj=Julian epoch, epb=Besselian epoch\n");
    fprintf(stderr,"         ts=seconds since 1950-01-01, tsu=Unix sec, tsi=IRAF sec\n");
    fprintf(stderr,"  @file: First one or two columns are in itype format\n");
    fprintf(stderr,"     -a: Append date to input file, if there is one\n");
    fprintf(stderr,"     -d: Print date without time\n");
    fprintf(stderr,"     -e: Print output as ET/TDT/TT converting from UT\n");
    fprintf(stderr,"     -f: Format for output number (C printf)\n");
    fprintf(stderr,"     -n: Number of decimal places in sec, epoch, JD\n");
    fprintf(stderr,"     -v: Verbose\n");
    fprintf(stderr,"    now: Output current date and time\n");
    exit (1);
}

static void
ConvertDate (intype, outtype, datestring, timestring)

int	intype;		/* Type of input date */
int	outtype;	/* Type of output date */
char	*datestring;	/* Input date string */
char	*timestring;	/* Input time string */

{
    double vdate, vtime, epoch, jd, ts, time1, jd1, ts1, epoch1;
    double doy, vdoy;
    int year, vyear;
    int lfd, oldfits;
    char outform[16];
    char *fitsdate, *newfdate;
    char temp[64];
    char ts0[8];
    char *tchar;
    time_t ut0, utc;	/* UTC seconds since 1970-01-01T0:00 */
    struct tm *ltm;	/* Local time structure */
    int lt;		/* Local time return */
    int its, its1;
    time_t lts;

    strcpy (ts0, "00:00:00");
    if (outform0 == NULL) {
	if (outtype == DTUNIX || outtype == DTIRAF)
	    strcpy (outform,"%d\n");
	else
	    sprintf (outform,"%%.%df\n",ndec);
	}
    else {
	strncpy (outform,outform0, 16);
	strcat (outform, "\n");
	}
    if (outtype < 1)
	outtype = DTFITS;

    if (datestring != NULL) {
	if (!strcmp (datestring, "now")) {
	    intype = DTVIG;

	    /* Get current UTC */
	    ut0 = (time_t) 0;
	    utc = time (&ut0);

	    /* Get local time and convert to vigesimal */
	    ltm = localtime (&utc);
	    vtime = (double) ltm->tm_hour + (0.01 * (double) ltm->tm_min) +
		(0.0001 * (double) ltm->tm_sec);
	    vdate = 1900.0 + (double) ltm->tm_year +
		(0.01 * (double) (1 + ltm->tm_mon)) +
		(0.0001 * (double) ltm->tm_mday);
	    timestring = ts0;
	    }
	}

    switch (intype) {
	case DTVIG:
	    if (datestring != NULL) {
		if (strcmp (datestring, "now")) {
		    vdate = atof (datestring);
		    if (timestring != NULL) 
			vtime = atof (timestring);
		    else {
			timestring = ts0;
			vtime = 0.0;
			}
		    }
		if (verbose)
		    printf ("%s %s -> ", datestring, timestring);

		switch (outtype) {
		    case DTEP:
			epoch = dt2ep (vdate, vtime);
			printf (outform, epoch);
			break;
		    case DTEPB:
			epoch = dt2epb (vdate, vtime);
			printf (outform, epoch);
			break;
		    case DTEPJ:
			epoch = dt2epj (vdate, vtime);
			printf (outform, epoch);
			break;
		    case DTFITS:
			fitsdate = dt2fd (vdate, vtime);
			fitsdate = fd2et (fitsdate);
			if (dateonly) {
			    tchar = strchr (fitsdate, 'T');
			    if (tchar != NULL)
				*tchar = (char) 0;
			    }
			printf ("%s\n", fitsdate);
			break;
		    case DTJD:
			jd = dt2jd (vdate, vtime);
			if (et) jd = jd2jed (jd);
			printf (outform, jd);
			break;
		    case DTMJD:
			jd = dt2mjd (vdate, vtime);
			printf (outform, jd);
			break;
		    case DT1950:
			ts = dt2ts (vdate, vtime);
			if (et) ts = ts2ets (ts);
			printf (outform, ts);
			break;
		    case DTIRAF:
			its = dt2tsi (vdate, vtime);
			printf (outform, its);
			break;
		    case DTUNIX:
			lts = dt2tsu (vdate, vtime);
			printf (outform, lts);
			break;
		    case DTDOY:
			dt2doy (vdate, vtime, &year, &doy);
			sprintf (temp, outform, doy);
			printf ("%04d %s\n", year, temp);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
		}
	    break;
	case DTFITS:
	    if (datestring != NULL) {
		if (strchr (datestring,'/'))
		    oldfits = 1;
		else
		    oldfits = 0;
		if (timestring != NULL) {
		    if (oldfits) {
			time1 = str2dec (timestring);
			lfd = strlen (datestring) + strlen (timestring) + 2;
			fitsdate = (char *) calloc (1, lfd);
			strcpy (fitsdate, datestring);
			}
		    else {
			lfd = strlen (datestring) + strlen (timestring) + 2;
			fitsdate = (char *) calloc (1, lfd);
			strcpy (fitsdate, datestring);
			strcat (fitsdate, "T");
			strcat (fitsdate, timestring);
			}
		    }
		else
		    fitsdate = datestring;
		if (verbose)
		    printf ("%s -> ", fitsdate);
		switch (outtype) {
		    case DTEP:
			epoch = fd2ep (fitsdate);
			if (oldfits && timestring) {
			    epoch1 = fd2ep (timestring);
			    epoch = epoch + epoch1;
			    }
			printf (outform, epoch);
			break;
		    case DTEPB:
			epoch = fd2epb (fitsdate);
			if (oldfits && timestring) {
			    epoch1 = fd2epb (timestring);
			    epoch = epoch + epoch1;
			    }
			printf (outform, epoch);
			break;
		    case DTEPJ:
			epoch = fd2epj (fitsdate);
			if (oldfits && timestring) {
			    epoch1 = fd2epj (timestring);
			    epoch = epoch + epoch1;
			    }
			printf (outform, epoch);
			break;
		    case DTVIG:
			fd2dt (fitsdate, &vdate, &vtime);
			if (oldfits) {
			    if (timestring == NULL)
				vtime = 0.0;
			    else
				vtime = time1;
			    }
			if (dateonly)
			    printf ("%9.4f\n", vdate);
			else
			    printf ("%9.4f %10.7f\n", vdate, vtime);
			break;
		    case DTFITS:
			newfdate = fd2fd (fitsdate);
			if (oldfits && timestring) {
			    strcat (newfdate, "T");
			    strcat (newfdate, timestring);
			    }
			fitsdate = fd2et (fitsdate);
			if (dateonly) {
			    tchar = strchr (newfdate, 'T');
			    if (tchar != NULL)
				*tchar = (char) 0;
			    }
			printf ("%s\n", newfdate);
			break;
		    case DTJD:
			jd = fd2jd (fitsdate);
			if (oldfits && timestring) {
			    jd1 = fd2jd (timestring);
			    jd = jd + jd1;
			    }
			if (et) jd = jd2jed (jd);
			printf (outform, jd);
			break;
		    case DTMJD:
			jd = fd2mjd (fitsdate);
			if (oldfits && timestring) {
			    jd1 = fd2jd (timestring);
			    jd = jd + jd1;
			    }
			printf (outform, jd);
			break;
		    case DTOF:
			newfdate = fd2of (fitsdate);
			if (oldfits && timestring) {
			    free (newfdate);
			    newfdate = fd2ofd (fitsdate);
			    strcat (newfdate, " ");
			    strcat (newfdate, timestring);
			    }
			printf ("%s\n", newfdate);
			break;
		    case DTOFD:
			newfdate = fd2ofd (fitsdate);
			printf ("%s\n", newfdate);
			break;
		    case DTOFT:
			newfdate = fd2oft (fitsdate);
			if (oldfits && timestring) {
			    strcpy (newfdate, timestring);
			    }
			printf ("%s\n", newfdate);
			break;
		    case DT1950:
			ts = fd2ts (fitsdate);
			if (oldfits && timestring) {
			    ts1 = fd2ts (timestring);
			    ts = ts + ts1;
			    }
			if (et) ts = ts2ets (ts);
			printf (outform, ts);
			break;
		    case DTIRAF:
			its = fd2tsi (fitsdate);
			if (oldfits && timestring) {
			    its1 = (int) fd2ts (timestring);
			    its = its + its1;
			    }
			printf (outform, its);
			break;
		    case DTUNIX:
			lts = fd2tsu (fitsdate);
			if (oldfits && timestring) {
			    its1 = (int) fd2ts (timestring);
			    lts = lts + (time_t) its1;
			    }
			printf (outform, lts);
			break;
		    case DTDOY:
			fd2doy (fitsdate, &year, &doy);
			if (oldfits && timestring) {
			    ts1 = (int) fd2ts (timestring);
			    doy = doy + (ts1 / 86400.0);
			    }
			sprintf (temp, outform, doy);
			printf ("%04d %s\n", year, temp);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
		}
	    break;
	case DTDOY:
	    if (datestring != NULL) {
		if (strcmp (datestring, "now")) {
		    vyear = atoi (datestring);
		    if (timestring != NULL) 
			vdoy = atof (timestring);
		    else
			vdoy = 1.0;
		    }
		if (verbose) {
		    sprintf (temp, outform, vdoy);
		    printf ("%04d %s -> ", vyear, temp);
		    }

		switch (outtype) {
		    case DTEP:
			epoch = doy2ep (vyear, vdoy);
			printf (outform, epoch);
			break;
		    case DTEPB:
			epoch = doy2epb (vyear, vdoy);
			printf (outform, epoch);
			break;
		    case DTEPJ:
			epoch = doy2epj (vyear, vdoy);
			printf (outform, epoch);
			break;
		    case DTFITS:
			fitsdate = doy2fd (vyear, vdoy);
			fitsdate = fd2et (fitsdate);
			if (dateonly) {
			    tchar = strchr (fitsdate, 'T');
			    if (tchar != NULL)
				*tchar = (char) 0;
			    }
			printf ("%s\n", fitsdate);
			break;
		    case DTJD:
			jd = doy2jd (vyear, vdoy);
			if (et) jd = jd2jed (jd);
			printf (outform, jd);
			break;
		    case DTMJD:
			jd = doy2mjd (vyear, vdoy);
			printf (outform, jd);
			break;
		    case DT1950:
			ts = doy2ts (vyear, vdoy);
			if (et) ts = ts2ets (ts);
			printf (outform, ts);
			break;
		    case DTIRAF:
			its = doy2tsi (vyear, vdoy);
			printf (outform, its);
			break;
		    case DTUNIX:
			lts = doy2tsu (vyear, vdoy);
			printf (outform, lts);
			break;
		    case DTVIG:
			doy2dt (vyear, vdoy, &vdate, &vtime);
			if (dateonly)
			    printf ("%9.4f\n", vdate);
			else
			    printf ("%9.4f %10.7f\n", vdate, vtime);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
		}
	    break;
	case DTJD:
	    if (datestring != NULL) {
		jd = atof (datestring);
		if (timestring != NULL) {
		    vtime = str2dec (timestring);
		    jd = jd + vtime / 24.0;
		    }
		if (verbose)
		    printf ("%.5f -> ", jd);
		switch (outtype) {
		    case DTEP:
			epoch = jd2ep (jd);
			printf (outform, epoch);
			break;
		    case DTEPB:
			epoch = jd2epb (jd);
			printf (outform, epoch);
			break;
		    case DTEPJ:
			epoch = jd2epj (jd);
			printf (outform, epoch);
			break;
		    case DTMJD:
			jd = jd2mjd (jd);
			printf (outform, jd);
			break;
		    case DTVIG:
			jd2dt (jd, &vdate, &vtime);
			if (dateonly)
			    printf ("%9.4f\n", vdate);
			else
			    printf ("%9.4f %10.7f\n", vdate, vtime);
			break;
		    case DTFITS:
			fitsdate = jd2fd (jd);
			fitsdate = fd2et (fitsdate);
			if (dateonly) {
			    tchar = strchr (fitsdate, 'T');
			    if (tchar != NULL)
				*tchar = (char) 0;
			    }
			printf ("%s\n", fitsdate);
			break;
		    case DT1950:
			ts = jd2ts (jd);
			if (et) ts = ts2ets (ts);
			printf (outform, ts);
			break;
		    case DTIRAF:
			its = jd2tsi (jd);
			printf (outform, its);
			break;
		    case DTUNIX:
			lts = jd2tsu (jd);
			printf (outform, lts);
			break;
		    case DTDOY:
			jd2doy (jd, &year, &doy);
			sprintf (temp, outform, doy);
			printf ("%04d %s\n", year, temp);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
		}
	    break;
	case DTMJD:
	    if (datestring != NULL) {
		jd = atof (datestring);
		if (timestring != NULL) {
		    vtime = str2dec (timestring);
		    jd = jd + vtime / 24.0;
		    }
		if (verbose)
		    printf ("%.5f -> ", jd);
		switch (outtype) {
		    case DTEP:
			epoch = mjd2ep (jd);
			printf (outform, epoch);
			break;
		    case DTEPB:
			epoch = mjd2epb (jd);
			printf (outform, epoch);
			break;
		    case DTEPJ:
			epoch = mjd2epj (jd);
			printf (outform, epoch);
			break;
		    case DTJD:
			jd = mjd2jd (jd);
			if (et) jd = jd2jed (jd);
			printf (outform, jd);
			break;
		    case DTVIG:
			mjd2dt (jd, &vdate, &vtime);
			if (dateonly)
			    printf ("%9.4f\n", vdate);
			else
			    printf ("%9.4f %10.7f\n", vdate, vtime);
			break;
		    case DTFITS:
			fitsdate = mjd2fd (jd);
			fitsdate = fd2et (fitsdate);
			if (dateonly) {
			    tchar = strchr (fitsdate, 'T');
			    if (tchar != NULL)
				*tchar = (char) 0;
			    }
			printf ("%s\n", fitsdate);
			break;
		    case DT1950:
			ts = mjd2ts (jd);
			if (et) ts = ts2ets (ts);
			printf (outform, ts);
			break;
		    case DTDOY:
			mjd2doy (jd, &year, &doy);
			sprintf (temp, outform, doy);
			printf ("%04d %s\n", year, temp);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
		}
	    break;
	case DT1950:
	    if (datestring != NULL) {
    		if (strcmp (datestring, "now"))
		    ts = atof (datestring);
		if (verbose)
		    printf ("%.3f -> ", ts);
		switch (outtype) {
		    case DTEP:
			epoch = ts2ep (ts);
			printf (outform, epoch);
			break;
		    case DTEPB:
			epoch = ts2epb (ts);
			printf (outform, epoch);
			break;
		    case DTEPJ:
			epoch = ts2epj (ts);
			printf (outform, epoch);
			break;
		    case DTVIG:
			ts2dt (ts, &vdate, &vtime);
			if (dateonly)
			    printf ("%9.4f\n", vdate);
			else
			    printf ("%9.4f %10.7f\n", vdate, vtime);
			break;
		    case DTFITS:
			fitsdate = ts2fd (ts);
			fitsdate = fd2et (fitsdate);
			if (dateonly) {
			    tchar = strchr (fitsdate, 'T');
			    if (tchar != NULL)
				*tchar = (char) 0;
			    }
			printf ("%s\n", fitsdate);
			break;
		    case DTJD:
			jd = ts2jd (ts);
			if (et) jd = jd2jed (jd);
			printf (outform, jd);
			break;
		    case DTMJD:
			jd = ts2mjd (ts);
			printf (outform, jd);
			break;
		    case DT1950:
			printf (outform, ts);
			if (et) ts = ts2ets (ts);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
		}
	    break;
	case DTIRAF:
	    if (datestring != NULL) {
    		if (strcmp (datestring, "now"))
		    its = atoi (datestring);
		if (verbose)
		    printf ("%d -> ", its);
		switch (outtype) {
		    case DTVIG:
			tsi2dt (its, &vdate, &vtime);
			if (dateonly)
			    printf ("%9.4f\n", vdate);
			else
			    printf ("%9.4f %10.7f\n", vdate, vtime);
			break;
		    case DTFITS:
			fitsdate = tsi2fd (its);
			fitsdate = fd2et (fitsdate);
			if (dateonly) {
			    tchar = strchr (fitsdate, 'T');
			    if (tchar != NULL)
				*tchar = (char) 0;
			    }
			printf ("%s\n", fitsdate);
			break;
		    case DT1950:
			ts = tsi2ts (its);
			if (et) ts = ts2ets (ts);
			printf (outform, ts);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
		}
	    break;
	case DTUNIX:
	    if (datestring != NULL) {
    		if (strcmp (datestring, "now"))
		    lts = (time_t) (atof (datestring) + 0.5);
		if (verbose)
		    printf ("%d -> ", lts);
		switch (outtype) {
		    case DTVIG:
			tsu2dt (lts, &vdate, &vtime);
			if (dateonly)
			    printf ("%9.4f\n", vdate);
			else
			    printf ("%9.4f %10.7f\n", vdate, vtime);
			break;
		    case DTFITS:
			fitsdate = tsu2fd (lts);
			fitsdate = fd2et (fitsdate);
			if (dateonly) {
			    tchar = strchr (fitsdate, 'T');
			    if (tchar != NULL)
				*tchar = (char) 0;
			    }
			printf ("%s\n", fitsdate);
			break;
		    case DTIRAF:
			its = tsu2tsi (lts);
			printf (outform, its);
			break;
		    case DT1950:
			ts = tsu2ts (lts);
			if (et) ts = ts2ets (ts);
			printf (outform, ts);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
		}
	    break;
	case DTEP:
	    if (datestring != NULL) {
		epoch = atof (datestring);
		if (verbose)
		    printf ("%.5f -> ", epoch);
		switch (outtype) {
		    case DTEPB:
			epoch = ep2epb (epoch);
			printf (outform, epoch);
			break;
		    case DTEPJ:
			epoch = ep2epj (epoch);
			printf (outform, epoch);
			break;
		    case DTVIG:
			ep2dt (epoch, &vdate, &vtime);
			if (dateonly)
			    printf ("%9.4f\n", vdate);
			else
			    printf ("%9.4f %10.7f\n", vdate, vtime);
			break;
		    case DTFITS:
			fitsdate = ep2fd (epoch);
			fitsdate = fd2et (fitsdate);
			if (dateonly) {
			    tchar = strchr (fitsdate, 'T');
			    if (tchar != NULL)
				*tchar = (char) 0;
			    }
			printf ("%s\n", fitsdate);
			break;
		    case DTJD:
			jd = ep2jd (epoch);
			if (et) jd = jd2jed (jd);
			printf (outform, jd);
			break;
		    case DTMJD:
			jd = ep2mjd (epoch);
			printf (outform, jd);
			break;
		    case DT1950:
			ts = ep2ts (epoch);
			if (et) ts = ts2ets (ts);
			printf (outform, ts);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
		}
	    break;
	case DTEPB:
	    if (datestring != NULL) {
		epoch = atof (datestring);
		if (verbose)
		    printf ("%.5f -> ", epoch);
		switch (outtype) {
		    case DTEP:
			epoch = epb2ep (epoch);
			printf (outform, epoch);
			break;
		    case DTEPJ:
			epoch = epb2epj (epoch);
			printf (outform, epoch);
			break;
		    case DTVIG:
			epb2dt (epoch, &vdate, &vtime);
			if (dateonly)
			    printf ("%9.4f\n", vdate);
			else
			    printf ("%9.4f %10.7f\n", vdate, vtime);
			break;
		    case DTFITS:
			fitsdate = epb2fd (epoch);
			fitsdate = fd2et (fitsdate);
			if (dateonly) {
			    tchar = strchr (fitsdate, 'T');
			    if (tchar != NULL)
				*tchar = (char) 0;
			    }
			printf ("%s\n", fitsdate);
			break;
		    case DTJD:
			jd = epb2jd (epoch);
			if (et) jd = jd2jed (jd);
			printf (outform, jd);
			break;
		    case DTMJD:
			jd = epb2mjd (epoch);
			printf (outform, jd);
			break;
		    case DT1950:
			ts = epb2ts (epoch);
			if (et) ts = ts2ets (ts);
			printf (outform, ts);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
		}
	    break;
	case DTEPJ:
	    if (datestring != NULL) {
		epoch = atof (datestring);
		if (verbose)
		    printf ("%.5f -> ", epoch);
		switch (outtype) {
		    case DTEP:
			epoch = epj2ep (epoch);
			printf (outform, epoch);
			break;
		    case DTEPB:
			epoch = epj2epb (epoch);
			printf (outform, epoch);
			break;
		    case DTVIG:
			epj2dt (epoch, &vdate, &vtime);
			if (dateonly)
			    printf ("%9.4f\n", vdate);
			else
			    printf ("%9.4f %10.7f\n", vdate, vtime);
			break;
		    case DTFITS:
			fitsdate = epj2fd (epoch);
			fitsdate = fd2et (fitsdate);
			if (dateonly) {
			    tchar = strchr (fitsdate, 'T');
			    if (tchar != NULL)
				*tchar = (char) 0;
			    }
			printf ("%s\n", fitsdate);
			break;
		    case DTJD:
			jd = epj2jd (epoch);
			if (et) jd = jd2jed (jd);
			printf (outform, jd);
			break;
		    case DTMJD:
			jd = epj2mjd (epoch);
			printf (outform, jd);
			break;
		    case DT1950:
			ts = epj2ts (epoch);
			if (et) ts = ts2ets (ts);
			printf (outform, ts);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
		}
	    break;
	case DTUT:
	    switch (outtype) {
		case DTEP:
		    epoch = ut2ep ();
		    printf (outform, epoch);
		    break;
		case DTEPB:
		    epoch = ut2epb ();
		    printf (outform, epoch);
		    break;
		case DTEPJ:
		    epoch = ut2epj ();
		    printf (outform, epoch);
		    break;
		case DTVIG:
		    ut2dt (&vdate, &vtime);
		    if (dateonly)
			printf ("%9.4f\n", vdate);
		    else
			printf ("%9.4f %10.7f\n", vdate, vtime);
		    break;
		case DTFITS:
		    newfdate = ut2fd ();
		    newfdate = fd2et (newfdate);
		    if (dateonly) {
			tchar = strchr (newfdate, 'T');
			if (tchar != NULL)
			    *tchar = (char) 0;
			}
		    printf ("%s\n", newfdate);
		    break;
		case DTJD:
		    jd = ut2jd ();
		    if (et) jd = jd2jed (jd);
		    printf (outform, jd);
		    break;
		case DTMJD:
		    jd = ut2mjd ();
		    printf (outform, jd);
		    break;
		case DT1950:
		    ts = ut2ts ();
			if (et) ts = ts2ets (ts);
		    printf (outform, ts);
		    break;
		case DTIRAF:
		    its = ut2tsi ();
		    printf (outform, its);
		    break;
		case DTUNIX:
		    lts = ut2tsu ();
		    printf (outform, lts);
		    break;
		case DTDOY:
		    ut2doy (&year, &doy);
		    sprintf (temp, outform, doy);
		    printf ("%04d %s\n", year, temp);
		    break;
		default:
		    printf ("*** Unknown output type %d\n", outtype);
		}
	    break;
	case DTLT:
	    switch (outtype) {
		case DTVIG:
		    lt2dt (&vdate, &vtime);
		    if (dateonly)
			printf ("%9.4f\n", vdate);
		    else
			printf ("%9.4f %10.7f\n", vdate, vtime);
		    break;
		case DTFITS:
		    newfdate = lt2fd ();
		    newfdate = fd2et (newfdate);
		    if (dateonly) {
			tchar = strchr (newfdate, 'T');
			if (tchar != NULL)
			    *tchar = (char) 0;
			}
		    printf ("%s\n", newfdate);
		    break;
		case DT1950:
		    ts = lt2ts ();
		    if (et) ts = ts2ets (ts);
		    printf (outform, ts);
		    break;
		case DTIRAF:
		    its = lt2tsi ();
		    printf (outform, its);
		    break;
		case DTUNIX:
		    lts = lt2tsu ();
		    printf (outform, lts);
		    break;
		default:
		    printf ("*** Unknown output type %d\n", outtype);
		}
	    break;
	default:
	    printf ("*** Unknown input type %d\n", intype);
	}

    return;
}

/* Dec 20 1999	New program
 *
 * Jan 20 2000	Add functions to convert Julian and Besselian epochs
 * Jan 21 2000	Add functions to convert to old FITS date and time
 * Jan 26 2000	Add modified Julian date conversions
 * Mar  2 2000	Add "today" to use current date and time
 * Mar  2 2000	Change date and time variables to vdate and vtime
 * Mar  3 2000	Add -d option to print only date
 * Mar 14 2000	Allow multiple conversions on one command line
 * Mar 24 2000	Add current local time, current UT, and IRAF and Unix seconds
 * May 31 2000	Fix bug which failed to convert to Julian or Bessellian epochs
 *
 * May 25 2001	Add year,day-of-year conversions
 *
 * Apr  8 2002	Change all long declarations to time_t
 * Apr  9 2002	Fix declaration of ConvertDate(); fix other bugs
 * Jul  8 2002	Check for negative dates if argument starts with -
 * Aug 27 2002	Add option to input file of dates to convert; fix output format
 * Aug 29 2002	Add option to append output to input line
 * Aug 30 2002	Add option to convert output to Ephemeris Time (ET/TDT/TT)
 */
