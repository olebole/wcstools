/*** File libwcs/txtread.c
 *** October 16, 1997
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

/* int txtread()	Read ASCII catalog stars in specified region
 * int txtrnum()	Read ASCII catalog stars with specified numbers
 * int txtopen()	Open ASCII catalog catalog, returning number of entries
 * char *txtstar()	Get ASCII catalog catalog entry for one star
 * double txtgetra()	Return double right ascension in degrees
 * double txtgetdec()	Return double declination in degrees
 * double txtgetmag()	Return magnitude from ASCII catalog line
 * int txtgetobj()	Return object name from ASCII catalog line
 * int txtgetc()	Return n'th character entry from ASCII catalog line
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "fitshead.h"
#include "wcs.h"

static int nstars;	/* Number of stars in catalog */
#define ABS(a) ((a) < 0 ? (-(a)) : (a))

int txtopen();
int txtparse();
int txtcol();
int txtgetk();
char *txtstar();
static int txtgeti4();
static double txtgetra();
static double txtgetdec();
static double txtgetmag();
static int txtgetc();

static char newline = 10;
static char txt = 9;
static char *txthead;
static char *txtdata;
static char *txtbuff;
static int nentry = 0;
static int entid, entra, entdec, entmag, entpeak;

#define MAXCOL	50
static char *colname[MAXCOL];	/* Column headers */
static int lcol[MAXCOL];

/* TXTREAD -- Read txt table stars in specified region */

int
txtread (txtcat,cra,cdec,dra,ddec,drad,mag1,mag2,nstarmax,
	 tnum,tra,tdec,tmag,tpeak,nlog)

char	*txtcat;	/* Name of reference star catalog file */
double	cra;		/* Search center J2000 right ascension in degrees */
double	cdec;		/* Search center J2000 declination in degrees */
double	dra;		/* Search half width in right ascension in degrees */
double	ddec;		/* Search half-width in declination in degrees */
double	drad;		/* Limiting separation in degrees (ignore if 0) */
double	mag1,mag2;	/* Limiting magnitudes (none if equal) */
int	nstarmax;	/* Maximum number of stars to be returned */
double	*tnum;		/* Array of UJ numbers (returned) */
double	*tra;		/* Array of right ascensions (returned) */
double	*tdec;		/* Array of declinations (returned) */
double	*tmag;		/* Array of magnitudes (returned) */
int	*tpeak;		/* Array of peak counts (returned) */
int	nlog;
{
    double ra1,ra2;	/* Limiting right ascensions of region in degrees */
    double dec1,dec2;	/* Limiting declinations of region in degrees */
    double dist = 0.0;  /* Distance from search center in degrees */
    double faintmag=0.0; /* Faintest magnitude */
    double maxdist=0.0; /* Largest distance */
    int faintstar=0;    /* Faintest star */
    int farstar=0;      /* Most distant star */
    double *tdist;      /* Array of distances to stars */

    int wrap;
    int jstar;
    int nstar;
    double ra,dec;
    double mag;
    double num;
    int peak, i;
    int istar;
    int verbose;
    char *line;

    line = 0;
    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

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
	fprintf (stderr,"TXTREAD: RA: %s - %s  Dec: %s - %s\n",
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

/* Logging interval */
    nstar = 0;
    tdist = (double *) malloc (nstarmax * sizeof (double));

    if ((nstars = txtopen (txtcat)) > 0) {
	jstar = 0;
	line = txtdata;

    /* Loop through catalog */
	for (istar = 1; istar <= nstars; istar++) {
	    line = txtstar (istar, line);
	    if (line == NULL) {
		fprintf (stderr,"TXTREAD: Cannot read star %d\n", istar);
		break;
		}

	/* Extract selected fields  */
	    num = txtgetr8 (line, entid);	/* ID number */
	    ra = txtgetra (line, entra);	/* Right ascension in degrees */
	    dec = txtgetdec (line, entdec);	/* Declination in degrees */
	    mag = txtgetr8 (line, entmag);	/* Magnitude */
	    peak = txtgeti4 (line, entpeak);	/* Peak counts */
	    if (drad > 0)
		dist = wcsdist (cra,cdec,ra,dec);
	    else
		dist = 0.0;

	/* Check magnitude amd position limits */
	    if ((mag1 == mag2 || (mag >= mag1 && mag <= mag2)) &&
		((wrap && (ra <= ra1 || ra >= ra2)) ||
		(!wrap && (ra >= ra1 && ra <= ra2))) &&
		(drad == 0.0 || dist < drad) &&
     		(dec >= dec1 && dec <= dec2)) {

	    /* Save star position and magnitude in table */
		if (nstar < nstarmax) {
		    tnum[nstar] = num;
		    tra[nstar] = ra;
		    tdec[nstar] = dec;
		    tmag[nstar] = mag;
		    tpeak[nstar] = peak;
		    tdist[nstar] = dist;
		    if (dist > maxdist) {
			maxdist = dist;
			farstar = nstar;
			}
		    if (mag > faintmag) {
			faintmag = mag;
			faintstar = nstar;
			}
		    }

		/* If too many stars and radial search,
		   replace furthest star */
		else if (drad > 0 && dist < maxdist) {
		    tra[farstar] = ra;
		    tdec[farstar] = dec;
		    tmag[farstar] = mag;
		    tpeak[farstar] = peak;
		    tdist[farstar] = dist;
		    maxdist = 0.0;

		    /* Find new farthest star */
		    for (i = 0; i < nstarmax; i++) {
			if (tdist[i] > maxdist) {
			    maxdist = tdist[i];
			    farstar = i;
			    }
			}
		    }

		/* If too many stars, replace faintest star */
		else if (mag < faintmag) {
		    tra[faintstar] = ra;
		    tdec[faintstar] = dec;
		    tmag[faintstar] = mag;
		    tpeak[faintstar] = peak;
		    tdist[faintstar] = dist;
		    faintmag = 0.0;

		    /* Find new faintest star */
		    for (i = 0; i < nstarmax; i++) {
			if (tmag[i] > faintmag) {
			    faintmag = tmag[i];
			    faintstar = i;
			    }
			}
		    }
		
		nstar++;
		jstar++;
		if (nlog == 1)
		    fprintf (stderr,"TXTREAD: %11.6f: %9.5f %9.5f %5.2f %d    \n",
			   num,ra,dec,mag,peak);

	    /* End of accepted star processing */
		}

	/* Log operation */
	    if (nlog > 0 && istar%nlog == 0)
		fprintf (stderr,"TXTREAD: %5d / %5d / %5d sources catalog %s\r",
			jstar,istar,nstars,txtcat);

	/* End of star loop */
	    }

	/* End of open catalog file */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"TXTREAD: Catalog %s : %d / %d / %d found\n",txtcat,jstar,istar,nstars);

    free (txtbuff);

    if (nstar > nstarmax) {
	fprintf (stderr,"TXTREAD: %d stars found; only %d returned\n",
		 nstar,nstarmax);
	nstar = nstarmax;
	}

    free ((char *)tdist);
    return (nstar);
}


/* TXTRNUM -- Read txt table stars with specified numbers */

int
txtrnum (txtcat, nnum, tnum,tra,tdec,tmag,tpeak,nlog)

char	*txtcat;	/* Name of reference star catalog file */
int	nnum;		/* Number of stars to look for */
double	*tnum;		/* Array of star numbers to look for */
double	*tra;		/* Array of right ascensions (returned) */
double	*tdec;		/* Array of declinations (returned) */
double	*tmag;		/* Array of magnitudes (returned) */
int	*tpeak;		/* Array of peak counts (returned) */
int	nlog;
{
    int jnum;
    int nstar;
    double ra,dec;
    double mag;
    double num;
    int peak;
    int istar;
    char *line;

    line = 0;

    nstar = 0;
    if ((nstars = txtopen (txtcat)) <= 0) {
	fprintf (stderr,"TXTRNUM: Cannot read catalog %s\n", txtcat);
	return (0);
	}

    /* Loop through star list */
    for (jnum = 0; jnum < nnum; jnum++) {
	line = txtdata;

	/* Loop through catalog to star */
	for (istar = 1; istar <= nstars; istar++) {
	    line = txtstar (istar, line);
	    if (line == NULL) {
		fprintf (stderr,"TXTRNUM: Cannot read star %d\n", istar);
		break;
		}
	    num = txtgetr8 (line, entid);	/* ID number */
	    if (num == tnum[jnum])
		break;
	    }

	/* If star has been found in table */
	if (num == tnum[jnum]) {

	    /* Extract selected fields  */
	    ra = txtgetra (line, entra);	/* Right ascension in degrees */
	    dec = txtgetdec (line, entdec);	/* Declination in degrees */
	    mag = txtgetr8 (line, entmag);	/* Magnitude */
	    peak = txtgeti4 (line, entpeak);	/* Peak counts */

	    /* Save star position and magnitude in table */
	    tra[jnum] = ra;
	    tdec[jnum] = dec;
	    tmag[jnum] = mag;
	    tpeak[jnum] = peak;
	    nstar++;
	    if (nlog == 1)
		fprintf (stderr,"TXTRNUM: %11.6f: %9.5f %9.5f %5.2f %d    \n",
			 num,ra,dec,mag,peak);

	    /* End of accepted star processing */
	    }

	/* Log operation */
	if (nlog > 0 && jnum%nlog == 0)
	    fprintf (stderr,"TXTRNUM: %5d / %5d / %5d sources catalog %s\r",
		     nstar,jnum,nstars,txtcat);

	/* End of star loop */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"TXTRNUM: Catalog %s : %d / %d found\n",
		 txtcat,nstar,nstars);

    free (txtbuff);
    return (nstars);
}


/* TXTOPEN -- Open txt table catalog, returning number of entries */

int
txtopen (txtfile)

char *txtfile;	/* Tab table catalog file name */
{
    struct stat statbuff;
    FILE *fcat;
    char *headbuff;
    int nr, lfile, ientry;
    char *headlast, *txtnew, *txtline, *lastline;
    char headend[4];

    nentry = 0;
    
/* Find length of txt table catalog */
    if (stat (txtfile, &statbuff)) {
	fprintf (stderr,"TXTOPEN: Tab table catalog %s has no entries\n",txtfile);
	return (0);
	}
    else
	lfile = (int) statbuff.st_size;

/* Open txt table catalog */
    if (!(fcat = fopen (txtfile, "r"))) {
	fprintf (stderr,"TXTOPEN: Tab table catalog %s cannot be read\n",txtfile);
	return (0);
	}

/* Allocate buffer to hold entire catalog and read it */
    if ((txtbuff = malloc (lfile)) != NULL) {
	nr = fread (txtbuff, 1, lfile, fcat);
	if (nr < lfile) {
	    fprintf (stderr,"TXTOPEN: read only %d / %d bytes of file %s\n",
		     nr, lfile, txtfile);
	    (void) fclose (fcat);
	    return (0);
	    }
	txthead = txtbuff;
	headend[0] = '-';
	headend[1] = '-';
	headend[2] = 0;
	txtline = txtbuff;
	while (strncmp (txtline,headend, 2)) {
	    lastline = txtline;
	    txtline = strchr (txtline,newline) + 1;
	    }
	txthead = lastline;
	txtdata = strchr (txtline, newline) + 1;

    /* Extract positions of keywords we will want to use */
	ientry = 0;
	headbuff = txthead;
	if (txtparse (headbuff)) {
	    if (!(entid = txtcol ("ID")))
		entid = txtcol ("id");
	    if (!(entid = txtcol ("RA")))
		entid = txtcol ("ra");
	    if (!(entid = txtcol ("DEC")))
		entid = txtcol ("dec");
	    if (!(entid = txtcol ("MAG")))
		entid = txtcol ("mag");
	    if (!(entid = txtcol ("PEAK")))
		entid = txtcol ("peak");
	    }
	else {
	    fprintf (stderr,"TXTOPEN: No columns in txt table %s\n",txtfile);
	    return (0);
	    }

    /* Enumerate entries in txt table catalog by counting newlines */
	txtnew = strchr (txtdata, newline) + 1;
	txtnew = txtdata;
	nstars = 0;
	while ((txtnew = strchr (txtnew, newline)) != NULL) {
	    txtnew = txtnew + 1;
	    nstars = nstars + 1;
	    }
	}

    (void) fclose (fcat);
    return (nstars);
}


/* TXTSTAR -- Get txt table catalog entry for one star; return 0 if successful */

char *
txtstar (istar, line)

int istar;	/* Star sequence number in txt table catalog */
char *line;	/* Pointer to istar'th entry (returned updated) */
{
    int iline;
    char *nextline;

    if (istar > nstars) {
	fprintf (stderr, "TXTSTAR:  %d is not in catalog\n",istar);
	return (NULL);
	}
    else if (istar < 1 && line) {
	nextline = strchr (line, newline) + 1;
	}
    else {
	nextline = txtdata;
	iline = 1;
	while (iline < istar) {
	    nextline = strchr (line, newline) + 1;
	    iline ++;
	    }
	}

    return (nextline);
}


/* TXTGETRA -- returns double right ascension in degrees */

static double
txtgetra (line, ientry)

char	*line;	/* txt table line */
int	ientry;	/* sequence of entry on line */
{
    char str[24];

    if (txtgetc (line, ientry, str, 24))
	return (0.0);
    else
	return (str2ra (str));
}


/* TXTGETDEC -- returns double declination in degrees */

static double
txtgetdec (line, ientry)

char	*line;	/* txt table line */
int	ientry;	/* sequence of entry on line */
{
    char str[24];

    if (txtgetc (line, ientry, str, 24))
	return (0.0);
    else
	return (str2dec (str));
}



/* TXTGETR8 -- returns 8-byte floating point number from txt table line */

static double
txtgetr8 (line, ientry)

char	*line;	/* txt table line */
int	ientry;	/* sequence of entry on line */
{
    char str[24];

    if (txtgetc (line, ientry, str, 24))
	return (0.0);
    else
        return (atof (str));
}


/* TXTGETI4 -- returns a 4-byte integer from txt table line */

static int
txtgeti4 (line, ientry)

char	*line;	/* txt table line */
int	ientry;	/* sequence of entry on line */
{
    char str[24];

    if (txtgetc (line, ientry, str, 24))
	return (0);
    else
        return ((int) atof (str));
}


/* TXTGETK -- returns a character entry from txt table line for named column */

int
txtgetk (line, keyword, string, maxchar)

char	*line;		/* txt table line */
char	*keyword;	/* column header of desired value */
char	*string;	/* character string (returned) */
int	maxchar;	/* Maximum number of characters in returned string */
{
    int ientry = txtcol (keyword);

    return (txtgetc (line, ientry, string, maxchar));
}


/* TXTGETC -- returns n'th entry from txt table line as character string */

static int
txtgetc (line, ientry, string, maxchar)

char	*line;		/* txt table line */
int	ientry;		/* sequence of entry on line */
char	*string;	/* character string (returned) */
int	maxchar;	/* Maximum number of characters in returned string */
{
    char *entry, *nextxt;
    int ient, ncstr;

    ient = 1;
    entry = line;
    if (ientry > nentry || ientry < 1)
	return (-1);
    for (ient  = 1; ient <= ientry; ient ++) {

    /* End ient'th entry with txt, newline, or end of string */
	if (ient < nentry) 
	    nextxt = strchr (entry, txt);
	else {
	    nextxt = strchr (entry, newline);
	    if (!nextxt)
		nextxt = strchr (entry, 0);
	    }
	if (!nextxt)
	    return (-1);
	if (ient < ientry)
	    entry = nextxt + 1;
	}
    ncstr = nextxt - entry;
    if (ncstr > maxchar - 1)
	ncstr = maxchar - 1;
    strncpy (string, entry, ncstr);
    string[ncstr] = 0;

    return (0);
}


/* TXTHGETR8 -- returns 8-byte floating point number from txt table header */

double
txthgetr8 (keyword)

char	*keyword;	/* sequence of entry on line */
{
    char *str0, *str1, *line, *head, str[24];
    int ncstr;

    head = txthead;
    str0 = 0;
    while (head < txtdata) {
	line = strsrch (head, keyword);
	if (line == txthead || line[-1] == newline) {
	    str0 = strchr (line, txt) + 1;
	    str1 = strchr (line, newline);
	    break;
	    }
	else
	    head = line;
	}
    if (str0) {
	ncstr = str1 - str0 - 1;
	strncpy (str, str0, ncstr);
	return (atof (str));
	}
    else
	return (0.0);
}


/* TXTHGETI4 -- returns a 4-byte integer from txt table header */

int
txthgeti4 (keyword)

char	*keyword;	/* sequence of entry on line */
{
    char *str0, *str1, *line, *head, str[24];
    int ncstr;

    head = txthead;
    str0 = 0;
    while (head < txtdata) {
	line = strsrch (head, keyword);
	if (line == txthead || line[-1] == newline) {
	    str0 = strchr (line, txt) + 1;
	    str1 = strchr (line, newline);
	    break;
	    }
	else
	    head = line;
	}
    if (str0) {
	ncstr = str1 - str0 - 1;
	strncpy (str, str0, ncstr);
	return ((int) atof (str));
	}
    else
	return (0);
}

/* Oct 16 1997	New file
 */
