/*** File libwcs/tabcread.c
 *** August 8, 1996
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "fitshead.h"

static int nstars;	/* Number of stars in catalog */
#define ABS(a) ((a) < 0 ? (-(a)) : (a))

static int tabopen();
static char *tabstar();
static int tabgeti4();
static double tabgetra();
static double tabgetdec();
static double tabgetr8();
static int tabhgeti4();
static double tabhgetr8();
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
tabread (tabcat,ra1,ra2,dec1,dec2,mag1,mag2,nstarmax,tnum,tra,tdec,tmag,tpeak,nlog)

char	*tabcat;	/* Name of reference star catalog file */
double	ra1,ra2;	/* Limiting right ascensions of region in degrees */
double	dec1,dec2;	/* Limiting declinations of region in degrees */
double	mag1,mag2;	/* Limiting magnitudes (none if equal) */
int	nstarmax;	/* Maximum number of stars to be returned */
double	*tnum;		/* Array of UJ numbers (returned) */
double	*tra;		/* Array of right ascensions (returned) */
double	*tdec;		/* Array of declinations (returned) */
double	*tmag;		/* Array of magnitudes (returned) */
int	*tpeak;		/* Array of peak counts (returned) */
int	nlog;
{
    int wrap;
    int jstar;
    int nstar;
    double ra,dec;
    double mag;
    double num;
    int peak;
    int istar;
    char *line;

    line = 0;

/* If RA range includes zero, split it in two */
    wrap = 0;
    if (ra1 > ra2)
	wrap = 1;
    else
	wrap = 0;

/* dec1 is always the smallest declination */
    if (dec1 > dec2) {
	dec = dec1;
	dec1 = dec2;
	dec2 = dec;
	}

/* mag1 is always the smallest magnitude */
    if (mag2 < mag1) {
	mag = mag2;
	mag2 = mag1;
	mag1 = mag;
	}

/* Logging interval */
    nstar = 0;

    if ((nstars = tabopen (tabcat))) {
	jstar = 0;
	line = tabdata;

    /* Loop through catalog */
	for (istar = 1; istar <= nstars; istar++) {
	    line = tabstar (istar, line);
	    if (line == NULL) {
		printf ("TABREAD: Cannot read star %d\n", istar);
		break;
		}

	/* Extract selected fields  */
	    num = tabgetr8 (line, entid);	/* ID number */
	    ra = tabgetra (line, entra);	/* Right ascension in degrees */
	    dec = tabgetdec (line, entdec);	/* Declination in degrees */
	    mag = tabgetr8 (line, entmag);	/* Magnitude */
	    peak = tabgeti4 (line, entpeak);	/* Peak counts */

	/* Check magnitude amd position limits */
	    if ((mag1 == mag2 || (mag >= mag1 && mag <= mag2)) &&
		((wrap && (ra <= ra1 || ra >= ra2)) ||
		(!wrap && (ra >= ra1 && ra <= ra2))) &&
     		(dec >= dec1 && dec <= dec2)) {

	    /* Save star position and magnitude in table */
		if (nstar <= nstarmax) {
		    tnum[nstar] = num;
		    tra[nstar] = ra;
		    tdec[nstar] = dec;
		    tmag[nstar] = mag;
		    tpeak[nstar] = peak;
		    }
		nstar++;
		jstar++;
		if (nlog == 1)
		    printf ("%11.6f: %9.5f %9.5f %5.2f %d    \n",
			   num,ra,dec,mag,peak);

	    /* End of accepted star processing */
		}

	/* Log operation */
	    if (nlog > 0 && istar%nlog == 0)
		printf ("%5d / %5d / %5d sources catalog %s\r",
			jstar,istar,nstars,tabcat);

	/* End of star loop */
	    }

	/* End of open catalog file */
	}

/* Summarize search */
    if (nlog > 0)
	printf ("Catalog %s : %d / %d / %d found\n",tabcat,jstar,istar,nstars);

    free (tabbuff);

    return (nstar);
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
	fprintf (stderr,"Tab table catalog %s has no entries\n",tabfile);
	return (0);
	}
    else
	lfile = (int) statbuff.st_size;

/* Open tab table catalog */
    if (!(fcat = fopen (tabfile, "r"))) {
	fprintf (stderr,"Tab table catalog %s cannot be read\n",tabfile);
	return (0);
	}

/* Allocate buffer to hold entire catalog and read it */
    if ((tabbuff = malloc (lfile))) {
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
	    else if (!strncmp (headbuff,"RA",2))
		entra = ientry;
	    else if (!strncmp (headbuff,"DEC",3))
		entdec = ientry;
	    else if (!strncmp (headbuff,"MAG",3))
		entmag = ientry;
	    else if (!strncmp (headbuff,"PEAK",4))
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
	while (tabnew = strchr (tabnew, newline)) {
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

static double
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

static int
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
 */
