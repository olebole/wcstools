/*** File libwcs/tabcread.c
 *** December 18, 1996
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

static int nstars;	/* Number of stars in catalog */
#define ABS(a) ((a) < 0 ? (-(a)) : (a))

static int tabopen();
static char *tabstar();
static int tabgeti4();
static double tabgetra();
static double tabgetdec();
static double tabgetr8();
static int tabgetc();

static char newline = 10;
static char tab = 9;
static char *tabhead;
static char *tabdata;
static char *tabbuff;
static int nentry = 5;
static int entid, entra, entdec, entmag, entpeak;

/* TABREAD -- Read tab table stars in specified region */

int
tabread (tabcat,cra,cdec,dra,ddec,drad,mag1,mag2,nstarmax,
	 tnum,tra,tdec,tmag,tpeak,nlog)

char	*tabcat;	/* Name of reference star catalog file */
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
	fprintf (stderr,"TABREAD: RA: %s - %s  Dec: %s - %s\n",
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

    if ((nstars = tabopen (tabcat)) > 0) {
	jstar = 0;
	line = tabdata;

    /* Loop through catalog */
	for (istar = 1; istar <= nstars; istar++) {
	    line = tabstar (istar, line);
	    if (line == NULL) {
		fprintf (stderr,"TABREAD: Cannot read star %d\n", istar);
		break;
		}

	/* Extract selected fields  */
	    num = tabgetr8 (line, entid);	/* ID number */
	    ra = tabgetra (line, entra);	/* Right ascension in degrees */
	    dec = tabgetdec (line, entdec);	/* Declination in degrees */
	    mag = tabgetr8 (line, entmag);	/* Magnitude */
	    peak = tabgeti4 (line, entpeak);	/* Peak counts */
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
		    fprintf (stderr,"TABREAD: %11.6f: %9.5f %9.5f %5.2f %d    \n",
			   num,ra,dec,mag,peak);

	    /* End of accepted star processing */
		}

	/* Log operation */
	    if (nlog > 0 && istar%nlog == 0)
		fprintf (stderr,"TABREAD: %5d / %5d / %5d sources catalog %s\r",
			jstar,istar,nstars,tabcat);

	/* End of star loop */
	    }

	/* End of open catalog file */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"TABREAD: Catalog %s : %d / %d / %d found\n",tabcat,jstar,istar,nstars);

    free (tabbuff);

    if (nstar > nstarmax) {
	fprintf (stderr,"TABREAD: %d stars found; only %d returned\n",
		 nstar,nstarmax);
	nstar = nstarmax;
	}

    free ((char *)tdist);
    return (nstar);
}


/* TABRNUM -- Read tab table stars with specified numbers */

int
tabrnum (tabcat, nstars, tnum,tra,tdec,tmag,tpeak,nlog)

char	*tabcat;	/* Name of reference star catalog file */
int	nstars;		/* Number of stars to look for */
double	*tnum;		/* Array of UJ numbers (returned) */
double	*tra;		/* Array of right ascensions (returned) */
double	*tdec;		/* Array of declinations (returned) */
double	*tmag;		/* Array of magnitudes (returned) */
int	*tpeak;		/* Array of peak counts (returned) */
int	nlog;
{
    int jstar;
    int nstar;
    double ra,dec;
    double mag;
    double num;
    int peak;
    int ncat;
    int istar;
    int verbose;
    char *line;

    line = 0;
    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    nstar = 0;

    /* Loop through star list */
    for (jstar = 0; jstar < nstars; jstar++) {

	/* Open catalog */
        if ((ncat = tabopen (tabcat)) > 0) {
	    jstar = 0;
	    line = tabdata;

	    /* Loop through catalog to star */
	    for (istar = 1; istar <= nstars; istar++) {
		line = tabstar (istar, line);
		if (line == NULL) {
		    fprintf (stderr,"TABRNUM: Cannot read star %d\n", istar);
		    break;
		    }
		num = tabgetr8 (line, entid);	/* ID number */
		if (num == tnum[jstar])
		    break;
		}

	    if (num == tnum[jstar]) {

		/* Extract selected fields  */
		ra = tabgetra (line, entra);	/* Right ascension in degrees */
		dec = tabgetdec (line, entdec);	/* Declination in degrees */
		mag = tabgetr8 (line, entmag);	/* Magnitude */
		peak = tabgeti4 (line, entpeak);	/* Peak counts */

	    /* Save star position and magnitude in table */
		tra[jstar] = ra;
		tdec[jstar] = dec;
		tmag[jstar] = mag;
		tpeak[jstar] = peak;
		nstar++;
		if (nlog == 1)
		    fprintf (stderr,"TABRNUM: %11.6f: %9.5f %9.5f %5.2f %d    \n",
			   num,ra,dec,mag,peak);

	    /* End of accepted star processing */
		}

	/* Log operation */
	    if (nlog > 0 && jstar%nlog == 0)
		fprintf (stderr,"TABRNUM: %5d / %5d / %5d sources catalog %s\r",
			nstar,jstar,nstars,tabcat);

	/* End of star loop */
	    }

	/* End of open catalog file */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"TABRNUM: Catalog %s : %d / %d found\n",
		 tabcat,nstar,nstars);

    free (tabbuff);

    return (nstars);
}


/* TABOPEN -- Open tab table catalog, returning number of entries */

static int
tabopen (tabfile)

char *tabfile;	/* Tab table catalog file name */
{
    struct stat statbuff;
    FILE *fcat;
    char *headbuff;
    int nr, lfile, ientry;
    char *headlast, *tabnew, *tabline, *lastline;
    char headend[4];
    
/* Find length of tab table catalog */
    if (stat (tabfile, &statbuff)) {
	fprintf (stderr,"TABOPEN: Tab table catalog %s has no entries\n",tabfile);
	return (0);
	}
    else
	lfile = (int) statbuff.st_size;

/* Open tab table catalog */
    if (!(fcat = fopen (tabfile, "r"))) {
	fprintf (stderr,"TABOPEN: Tab table catalog %s cannot be read\n",tabfile);
	return (0);
	}

/* Allocate buffer to hold entire catalog and read it */
    if ((tabbuff = malloc (lfile)) != NULL) {
	nr = fread (tabbuff, 1, lfile, fcat);
	if (nr < lfile) {
	    fprintf (stderr,"TABOPEN: read only %d / %d bytes of file %s\n",
		     nr, lfile, tabfile);
	    (void) fclose (fcat);
	    return (0);
	    }
	tabhead = tabbuff;
	headend[0] = '-';
	headend[1] = '-';
	headend[2] = 0;
	tabline = tabbuff;
	while (strncmp (tabline,headend, 2)) {
	    lastline = tabline;
	    tabline = strchr (tabline,newline) + 1;
	    }
	tabhead = lastline;
	tabdata = strchr (tabline, newline) + 1;

    /* Extract positions of keywords we will want to use */
	ientry = 0;
	nentry = 0;
	headbuff = tabhead;
	headlast = strchr (headbuff, newline);
	while (headbuff) {
	    ientry++;
	    if (!strncmp (headbuff,"ID",2))
		entid = ientry;
	    else if (!strncmp (headbuff,"id",2))
		entid = ientry;
	    else if (!strncmp (headbuff,"RA",2))
		entra = ientry;
	    else if (!strncmp (headbuff,"ra",2))
		entra = ientry;
	    else if (!strncmp (headbuff,"DEC",3))
		entdec = ientry;
	    else if (!strncmp (headbuff,"dec",3))
		entdec = ientry;
	    else if (!strncmp (headbuff,"MAG",3))
		entmag = ientry;
	    else if (!strncmp (headbuff,"mag",3))
		entmag = ientry;
	    else if (!strncmp (headbuff,"PEAK",4))
		entpeak = ientry;
	    else if (!strncmp (headbuff,"peak",4))
		entpeak = ientry;
	    nentry++;
	    headbuff = strchr (headbuff, tab) + 1;
	    if (headbuff > headlast)
		break;
	    }

    /* Enumerate entries in tab table catalog by counting newlines */
	tabnew = strchr (tabdata, newline) + 1;
	tabnew = tabdata;
	nstars = 0;
	while ((tabnew = strchr (tabnew, newline)) != NULL) {
	    tabnew = tabnew + 1;
	    nstars = nstars + 1;
	    }
	}

    (void) fclose (fcat);
    return (nstars);
}


/* TABSTAR -- Get tab table catalog entry for one star; return 0 if successful */

static char *
tabstar (istar, line)

int istar;	/* Star sequence number in tab table catalog */
char *line;	/* Pointer to istar'th entry (returned updated) */
{
    int iline;
    char *nextline;

    if (istar > nstars) {
	fprintf (stderr, "TABSTAR:  %d is not in catalog\n",istar);
	return (NULL);
	}
    else if (istar < 1 && line) {
	nextline = strchr (line, newline) + 1;
	}
    else {
	nextline = tabdata;
	iline = 1;
	while (iline < istar) {
	    nextline = strchr (line, newline) + 1;
	    iline ++;
	    }
	}

    return (nextline);
}


/* TABGETRA -- returns double right ascension in degrees */

static double
tabgetra (line, ientry)

char	*line;	/* tab table line */
int	ientry;	/* sequence of entry on line */
{
    char str[24];

    if (tabgetc (line, ientry, str, 24))
	return (0.0);
    else
	return (str2ra (str));
}


/* TABGETDEC -- returns double declination in degrees */

static double
tabgetdec (line, ientry)

char	*line;	/* tab table line */
int	ientry;	/* sequence of entry on line */
{
    char str[24];

    if (tabgetc (line, ientry, str, 24))
	return (0.0);
    else
	return (str2dec (str));
}



/* TABGETR8 -- returns 8-byte floating point number from tab table line */

static double
tabgetr8 (line, ientry)

char	*line;	/* tab table line */
int	ientry;	/* sequence of entry on line */
{
    char str[24];

    if (tabgetc (line, ientry, str, 24))
	return (0.0);
    else
        return (atof (str));
}


/* TABGETI4 -- returns a 4-byte integer from tab table line */

static int
tabgeti4 (line, ientry)

char	*line;	/* tab table line */
int	ientry;	/* sequence of entry on line */
{
    char str[24];

    if (tabgetc (line, ientry, str, 24))
	return (0);
    else
        return ((int) atof (str));
}


/* TABGETC -- Returns n'th entry from tab table line as character string */

static int
tabgetc (line, ientry, string, maxchar)

char	*line;		/* tab table line */
int	ientry;		/* sequence of entry on line */
char	*string;	/* character string (returned) */
int	maxchar;	/* Maximum number of characters in returned string */
{
    char *entry, *nextab;
    int ient, ncstr;

    ient = 1;
    entry = line;
    if (ientry > nentry || ientry < 0)
	return (-1);
    for (ient  = 1; ient <= ientry; ient ++) {

    /* End ient'th entry with tab, newline, or end of string */
	if (ient < nentry) 
	    nextab = strchr (entry, tab);
	else {
	    nextab = strchr (entry, newline);
	    if (!nextab)
		nextab = strchr (entry, 0);
	    }
	if (!nextab)
	    return (-1);
	if (ient < ientry)
	    entry = nextab + 1;
	}
    ncstr = nextab - entry;
    if (ncstr > maxchar - 1)
	ncstr = maxchar - 1;
    strncpy (string, entry, ncstr);
    string[ncstr] = 0;

    return (0);
}


/* TABHGETR8 -- returns 8-byte floating point number from tab table header */

double
tabhgetr8 (keyword)

char	*keyword;	/* sequence of entry on line */
{
    char *str0, *str1, *line, *head, str[24];
    int ncstr;

    head = tabhead;
    str0 = 0;
    while (head < tabdata) {
	line = strsrch (head, keyword);
	if (line == tabhead || line[-1] == newline) {
	    str0 = strchr (line, tab) + 1;
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


/* TABHGETI4 -- returns a 4-byte integer from tab table header */

int
tabhgeti4 (keyword)

char	*keyword;	/* sequence of entry on line */
{
    char *str0, *str1, *line, *head, str[24];
    int ncstr;

    head = tabhead;
    str0 = 0;
    while (head < tabdata) {
	line = strsrch (head, keyword);
	if (line == tabhead || line[-1] == newline) {
	    str0 = strchr (line, tab) + 1;
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

/* Jul 18 1996	New subroutines
 * Aug  6 1996	Remove unused variables after lint
 * Aug  8 1996	Fix bugs in entry reading and logging
 * Oct 15 1996  Add comparison when testing an assignment
 * Nov  5 1996	Drop unnecessary static declarations
 * Nov 13 1996	Return no more than maximum star number
 * Nov 13 1996	Write all error messages to stderr with subroutine names
 * Nov 15 1996  Implement search radius; change input arguments
 * Nov 19 1996	Allow lower case column headings
 * Dec 18 1996	Add UJCRNUM to read specified catalog entries
 * Dec 18 1996  Keep closest stars, not brightest, if searching within radius
 */
