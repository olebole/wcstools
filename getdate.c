/* File getdate.c
 * January 26, 2000
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

static void usage();

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
static int version = 0;	/* If 1, print only program name and version */
static int nread = 0;	/* Number of lines to read (0=all) */
static int nskip = 0;	/* Number of lines to skip */
static int tabout = 0;	/* If 1, separate output fields with tabs */
static int ndec = 5;	/* Number of decimal places in output */

main (ac, av)
int ac;
char **av;
{
    char *str;
    char *temp;
    char *datestring;
    char *timestring;
    int intype= 0;		/* Input date type */
    int outtype= 0;		/* Output date type */

    datestring = NULL;
    timestring = NULL;

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

	/* Read command */
	if (!strncmp (*av, "dt2", 3)) {
	    intype = DTVIG;
	    if (!strcmp (*av, "dt2ep"))
		outtype = DTEP;
	    else if (!strcmp (*av, "dt2epb"))
		outtype = DTEPB;
	    else if (!strcmp (*av, "dt2epj"))
		outtype = DTEPJ;
	    else if (!strcmp (*av, "dt2fd"))
		outtype = DTFITS;
	    else if (!strcmp (*av, "dt2jd"))
		outtype = DTJD;
	    else if (!strcmp (*av, "dt2mjd"))
		outtype = DTMJD;
	    else if (!strcmp (*av, "dt2ts"))
		outtype = DT1950;
	    else {
		printf ("Unknown conversion %s\n", *av);
		exit (1);
		}
	    }
	else if (!strncmp (*av, "ep2", 3)) {
	    intype = DTEP;
	    if (!strcmp (*av, "ep2epb"))
		outtype = DTEPB;
	    else if (!strcmp (*av, "ep2epj"))
		outtype = DTEPJ;
	    else if (!strcmp (*av, "ep2dt"))
		outtype = DTVIG;
	    else if (!strcmp (*av, "ep2fd"))
		outtype = DTFITS;
	    else if (!strcmp (*av, "ep2jd"))
		outtype = DTJD;
	    else if (!strcmp (*av, "ep2mjd"))
		outtype = DTMJD;
	    else if (!strcmp (*av, "ep2ts"))
		outtype = DT1950;
	    else {
		printf ("Unknown conversion %s\n", *av);
		exit (1);
		}
	    }
	else if (!strncmp (*av, "epb2", 3)) {
	    intype = DTEPB;
	    if (!strcmp (*av, "epb2ep"))
		outtype = DTEP;
	    else if (!strcmp (*av, "epb2epj"))
		outtype = DTEPJ;
	    else if (!strcmp (*av, "epb2dt"))
		outtype = DTVIG;
	    else if (!strcmp (*av, "epb2fd"))
		outtype = DTFITS;
	    else if (!strcmp (*av, "epb2jd"))
		outtype = DTJD;
	    else if (!strcmp (*av, "epb2mjd"))
		outtype = DTMJD;
	    else if (!strcmp (*av, "epb2ts"))
		outtype = DT1950;
	    else {
		printf ("Unknown conversion %s\n", *av);
		exit (1);
		}
	    }
	else if (!strncmp (*av, "epj2", 3)) {
	    intype = DTEPJ;
	    if (!strcmp (*av, "epj2ep"))
		outtype = DTEP;
	    else if (!strcmp (*av, "epj2epb"))
		outtype = DTEPB;
	    else if (!strcmp (*av, "epj2dt"))
		outtype = DTVIG;
	    else if (!strcmp (*av, "epj2fd"))
		outtype = DTFITS;
	    else if (!strcmp (*av, "epj2jd"))
		outtype = DTJD;
	    else if (!strcmp (*av, "epj2mjd"))
		outtype = DTMJD;
	    else if (!strcmp (*av, "epj2ts"))
		outtype = DT1950;
	    else {
		printf ("Unknown conversion %s\n", *av);
		exit (1);
		}
	    }
	else if (!strncmp (*av, "fd2", 3)) {
	    intype = DTFITS;
	    if (!strcmp (*av, "fd2dt"))
		outtype = DTVIG;
	    else if (!strcmp (*av, "fd2ep"))
		outtype = DTEP;
	    else if (!strcmp (*av, "fd2epb"))
		outtype = DTEPB;
	    else if (!strcmp (*av, "fd2epj"))
		outtype = DTEPJ;
	    else if (!strcmp (*av, "fd2fd"))
		outtype = DTFITS;
	    else if (!strcmp (*av, "fd2jd"))
		outtype = DTJD;
	    else if (!strcmp (*av, "fd2mjd"))
		outtype = DTMJD;
	    else if (!strcmp (*av, "fd2of"))
		outtype = DTOF;
	    else if (!strcmp (*av, "fd2ofd"))
		outtype = DTOFD;
	    else if (!strcmp (*av, "fd2oft"))
		outtype = DTOFT;
	    else if (!strcmp (*av, "fd2ts"))
		outtype = DT1950;
	    else {
		printf ("Unknown conversion %s\n", *av);
		exit (1);
		}
	    }
	else if (!strncmp (*av, "jd2", 3)) {
	    intype = DTJD;
	    if (!strcmp (*av, "jd2dt"))
		outtype = DTVIG;
	    else if (!strcmp (*av, "jd2ep"))
		outtype = DTEP;
	    else if (!strcmp (*av, "jd2epb"))
		outtype = DTEPB;
	    else if (!strcmp (*av, "jd2epj"))
		outtype = DTEPJ;
	    else if (!strcmp (*av, "jd2fd"))
		outtype = DTFITS;
	    else if (!strcmp (*av, "jd2mjd"))
		outtype = DTMJD;
	    else if (!strcmp (*av, "jd2ts"))
		outtype = DT1950;
	    else {
		printf ("Unknown conversion %s\n", *av);
		exit (1);
		}
	    }
	else if (!strncmp (*av, "mjd2", 3)) {
	    intype = DTJD;
	    if (!strcmp (*av, "mjd2dt"))
		outtype = DTVIG;
	    else if (!strcmp (*av, "mjd2ep"))
		outtype = DTEP;
	    else if (!strcmp (*av, "mjd2epb"))
		outtype = DTEPB;
	    else if (!strcmp (*av, "mjd2epj"))
		outtype = DTEPJ;
	    else if (!strcmp (*av, "mjd2fd"))
		outtype = DTFITS;
	    else if (!strcmp (*av, "mjd2jd"))
		outtype = DTJD;
	    else if (!strcmp (*av, "mjd2ts"))
		outtype = DT1950;
	    else {
		printf ("Unknown conversion %s\n", *av);
		exit (1);
		}
	    }
	else if (!strncmp (*av, "ts2", 3)) {
	    intype = DT1950;
	    if (!strcmp (*av, "ts2dt"))
		outtype = DTVIG;
	    else if (!strcmp (*av, "ts2ep"))
		outtype = DTEP;
	    else if (!strcmp (*av, "ts2epb"))
		outtype = DTEPB;
	    else if (!strcmp (*av, "ts2epj"))
		outtype = DTEPJ;
	    else if (!strcmp (*av, "ts2fd"))
		outtype = DTFITS;
	    else if (!strcmp (*av, "ts2jd"))
		outtype = DTJD;
	    else if (!strcmp (*av, "ts2mjd"))
		outtype = DTMJD;
	    else {
		printf ("Unknown conversion %s\n", *av);
		exit (1);
		}
	    }
	else if (*(str = *av) == '-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

		case 'v':	/* More verbosity */
		    verbose++;
		    break;

		case 'n':	/* Number of decimal places in output */
		    if (ac < 2)
			usage();
		    ndec = (int) (atof (*++av));
		    ac--;
		    break;

		default:
		    usage ();
		    break;
		}
	    }
	else if (datestring == NULL)
	    datestring = *av;
	else
	    timestring = *av;
	}
    ConvertDate (intype, outtype, datestring, timestring);

    return (0);
}

static void
usage ()

{
    if (version)
	exit (-1);
    fprintf (stderr,"Convert date and time between various formats\n");
    fprintf (stderr,"Usage: [-v][-n dec] itype2otype [date and/or time]\n");
    fprintf(stderr,"  itype: fd=FITS, dt=yyyy.mmdd, ep=epoch, jd=Julian Date\n");
    fprintf(stderr,"  otype: fd=FITS, dt=yyyy.mmdd, ep=epoch, jd=Julian Date\n");
    fprintf(stderr,"     -n: Number of decimal places in sec, epoch, JD\n");
    fprintf(stderr,"     -v: Verbose\n");
    exit (1);
}

static int
ConvertDate (intype, outtype, datestring, timestring)

int	intype;		/* Type of input date */
int	outtype;	/* Type of output date */
char	*datestring;	/* Input date string */
char	*timestring;	/* Input time string */

{
    double date, time, epoch, jd, ts, time1, jd1, ts1, epoch1;
    int lfd, oldfits;
    char *fitsdate, *newfdate;
    char ts0[8];
    char outform[16];

    strcpy (ts0, "00:00:00");
    sprintf (outform,"%%.%df\n",ndec);

    switch (intype) {
	case DTVIG:
	    if (datestring != NULL) {
		date = atof (datestring);
		if (timestring != NULL) 
		    time = atof (timestring);
		else {
		    timestring = ts0;
		    time = 0.0;
		    }
		if (verbose)
		    printf ("%s %s -> ", datestring, timestring);
		switch (outtype) {
		    case DTEP:
			epoch = dt2ep (date, time);
			printf (outform, epoch);
			break;
		    case DTEPB:
			epoch = dt2epb (date, time);
			printf (outform, epoch);
			break;
		    case DTEPJ:
			epoch = dt2epj (date, time);
			printf (outform, epoch);
			break;
		    case DTFITS:
			fitsdate = dt2fd (date, time);
			printf ("%s\n", fitsdate);
			break;
		    case DTJD:
			jd = dt2jd (date, time);
			printf (outform, jd);
			break;
		    case DTMJD:
			jd = dt2mjd (date, time);
			printf (outform, jd);
			break;
		    case DT1950:
			ts = dt2ts (date, time);
			printf (outform, ts);
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
			fd2dt (fitsdate, &date, &time);
			if (oldfits) {
			    if (timestring == NULL)
				time = 0.0;
			    else
				time = time1;
			    }
			printf ("%9.4f %10.7f\n", date, time);
			break;
		    case DTFITS:
			newfdate = fd2fd (fitsdate);
			if (oldfits && timestring) {
			    strcat (newfdate, "T");
			    strcat (newfdate, timestring);
			    }
			printf ("%s\n", newfdate);
			break;
		    case DTJD:
			jd = fd2jd (fitsdate);
			if (oldfits && timestring) {
			    jd1 = fd2jd (timestring);
			    jd = jd + jd1;
			    }
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
			printf (outform, ts);
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
		    time = str2dec (timestring);
		    jd = jd + time / 24.0;
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
			jd = jd2mjd (fitsdate);
			printf (outform, jd);
			break;
		    case DTVIG:
			jd2dt (jd, &date, &time);
			printf ("%.4f %.7f\n", date, time);
			break;
		    case DTFITS:
			fitsdate = jd2fd (jd);
			printf ("%s\n", fitsdate);
			break;
		    case DT1950:
			ts = jd2ts (jd);
			printf (outform, ts);
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
		    time = str2dec (timestring);
		    jd = jd + time / 24.0;
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
			jd = mjd2jd (fitsdate);
			printf (outform, jd);
			break;
		    case DTVIG:
			mjd2dt (jd, &date, &time);
			printf ("%.4f %.7f\n", date, time);
			break;
		    case DTFITS:
			fitsdate = mjd2fd (jd);
			printf ("%s\n", fitsdate);
			break;
		    case DT1950:
			ts = mjd2ts (jd);
			printf (outform, ts);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
		}
	    break;
	case DT1950:
	    if (datestring != NULL) {
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
			ts2dt (ts, &date, &time);
			printf ("%.4f %.7f\n", date, time);
			break;
		    case DTFITS:
			fitsdate = ts2fd (ts);
			printf ("%s\n", fitsdate);
			break;
		    case DTJD:
			jd = ts2jd (ts);
			printf (outform, jd);
			break;
		    case DTMJD:
			jd = ts2mjd (ts);
			printf (outform, jd);
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
			ep2dt (epoch, &date, &time);
			printf ("%.4f %.7f\n", date, time);
			break;
		    case DTFITS:
			fitsdate = ep2fd (epoch);
			printf ("%s\n", fitsdate);
			break;
		    case DTJD:
			jd = ep2jd (epoch);
			printf (outform, jd);
			break;
		    case DTMJD:
			jd = ep2mjd (epoch);
			printf (outform, jd);
			break;
		    case DT1950:
			ts = ep2ts (epoch);
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
			epb2dt (epoch, &date, &time);
			printf ("%.4f %.7f\n", date, time);
			break;
		    case DTFITS:
			fitsdate = epb2fd (epoch);
			printf ("%s\n", fitsdate);
			break;
		    case DTJD:
			jd = epb2jd (epoch);
			printf (outform, jd);
			break;
		    case DTMJD:
			jd = epb2mjd (epoch);
			printf (outform, jd);
			break;
		    case DT1950:
			ts = epb2ts (epoch);
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
			epj2dt (epoch, &date, &time);
			printf ("%.4f %.7f\n", date, time);
			break;
		    case DTFITS:
			fitsdate = epj2fd (epoch);
			printf ("%s\n", fitsdate);
			break;
		    case DTJD:
			jd = epj2jd (epoch);
			printf (outform, jd);
			break;
		    case DTMJD:
			jd = epj2mjd (epoch);
			printf (outform, jd);
			break;
		    case DT1950:
			ts = epj2ts (epoch);
			printf (outform, ts);
			break;
		    default:
			printf ("*** Unknown output type %d\n", outtype);
		    }
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
 * Jan 26 200	Addmodified Julian date conversions
 */
