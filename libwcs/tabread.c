/*** File libwcs/tabcread.c
 *** December 8, 1998
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

/* int tabread()	Read tab table stars in specified region
 * int tabrnum()	Read tab table stars with specified numbers
 * int tabrkey()	Read single keyword from tab table stars with specified numbers
 * int tabopen()	Open tab table catalog, returning number of entries
 * char *tabstar()	Get tab table catalog entry for one star
 * double tabgetra()	Return double right ascension in degrees
 * double tabgetdec()	Return double declination in degrees
 * double tabgetr8()	Return 8-byte floating point number from tab table line
 * int tabgeti4()	Return 4-byte integer from tab table line
 * int tabgetk()	Return character entry from tab table line for column
 * int tabgetc()	Return n'th character entry from tab table line
 * double tabhgetr8()	Return 8-byte floating point from tab table header
 * int tabhgeti4()	Return 4-byte integer from tab table header
 * int tabparse()	Make a table of column headings
 * int tabcol()		Search a table of column headings for a particlar entry
 * int tabsize()	Return length of file in bytes
 * int istab()		Return 1 if first line of file contains a tab, else 0
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include "fitshead.h"
#include "wcs.h"

static int nstars;	/* Number of stars in catalog */
#define ABS(a) ((a) < 0 ? (-(a)) : (a))

int tabopen();
int tabparse();
int tabcol();
int tabgetk();
char *tabstar();
static int tabgeti4();
static double tabgetra();
static double tabgetdec();
static double tabgetr8();
static int tabgetc();
void tabclose();

static int entid, entra, entdec, entmag, entpeak;
static int iline = 0;	/* Current line of data */

static char *tabdata;


/* TABREAD -- Read tab table stars in specified region */

int
tabread (tabcat,cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,mag2,nstarmax,
	 tnum,tra,tdec,tmag,tpeak,nlog)

char	*tabcat;	/* Name of reference star catalog file */
double	cra;		/* Search center J2000 right ascension in degrees */
double	cdec;		/* Search center J2000 declination in degrees */
double	dra;		/* Search half width in right ascension in degrees */
double	ddec;		/* Search half-width in declination in degrees */
double	drad;		/* Limiting separation in degrees (ignore if 0) */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
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
    int sysref=WCS_J2000;	/* Catalog coordinate system */
    double eqref=2000.0;	/* Catalog equinox */
    double epref=2000.0;	/* Catalog epoch */
    char cstr[32];

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
    wcscstr (cstr, sysout, eqout, epout);

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
	ra2str (rstr1, 16, ra1, 3);
        dec2str (dstr1, 16, dec1, 2);
	ra2str (rstr2, 16, ra2, 3);
        dec2str (dstr2, 16, dec2, 2);
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
	    if (num == 0.0)
		num = (double)istar;
	    ra = tabgetra (line, entra);	/* Right ascension in degrees */
	    dec = tabgetdec (line, entdec);	/* Declination in degrees */
	    wcscon (sysref, sysout, eqref, eqout, &ra, &dec, epout);
	    mag = tabgetr8 (line, entmag);	/* Magnitude */
	    peak = tabgeti4 (line, entpeak);	/* Peak counts */
	    if (drad > 0)
		dist = wcsdist (cra,cdec,ra,dec);
	    else
		dist = 0.0;

	/* Check magnitude and position limits */
	    if ((mag1 == mag2 || (mag >= mag1 && mag <= mag2)) &&
		((wrap && (ra <= ra1 || ra >= ra2)) ||
		(!wrap && (ra >= ra1 && ra <= ra2))) &&
		((drad > 0.0 && dist < drad) ||
     		(drad == 0.0 && dec >= dec1 && dec <= dec2))) {

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
		    tnum[farstar] = num;
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
		    tnum[faintstar] = num;
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
		    fprintf (stderr,"TABREAD: %11.6f: %9.5f %9.5f %s %5.2f %d    \n",
			   num,ra,dec,cstr,mag,peak);

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


    if (nstar > nstarmax) {
	fprintf (stderr,"TABREAD: %d stars found; only %d returned\n",
		 nstar,nstarmax);
	nstar = nstarmax;
	}

    tabclose();
    free ((char *)tdist);
    return (nstar);
}


/* TABRNUM -- Read tab table stars with specified numbers */

int
tabrnum (tabcat, nnum, sysout, eqout, epout, tnum,tra,tdec,tmag,tpeak,nlog)

char	*tabcat;	/* Name of reference star catalog file */
int	nnum;		/* Number of stars to look for */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
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
    int sysref=WCS_J2000;	/* Catalog coordinate system */
    double eqref=2000.0;	/* Catalog equinox */
    double epref=2000.0;	/* Catalog epoch */
    char cstr[32];

    line = 0;

    nstar = 0;
    if ((nstars = tabopen (tabcat)) <= 0) {
	fprintf (stderr,"TABRNUM: Cannot read catalog %s\n", tabcat);
	return (0);
	}

    wcscstr (cstr, sysout, eqout, epout);

    /* Loop through star list */
    for (jnum = 0; jnum < nnum; jnum++) {
	line = tabdata;

	/* Loop through catalog to star */
	for (istar = 1; istar <= nstars; istar++) {
	    line = tabstar (istar, line);
	    if (line == NULL) {
		fprintf (stderr,"TABRNUM: Cannot read star %d\n", istar);
		break;
		}
	    num = tabgetr8 (line, entid);	/* ID number */
	    if (num == 0.0)
		num = (double) istar;
	    if (num == tnum[jnum])
		break;
	    }

	/* If star has been found in table */
	if (num == tnum[jnum]) {

	    /* Extract selected fields  */
	    ra = tabgetra (line, entra);	/* Right ascension in degrees */
	    dec = tabgetdec (line, entdec);	/* Declination in degrees */
	    wcscon (sysref, sysout, eqref, eqout, &ra, &dec, epout);
	    mag = tabgetr8 (line, entmag);	/* Magnitude */
	    peak = tabgeti4 (line, entpeak);	/* Peak counts */

	    /* Save star position and magnitude in table */
	    tra[jnum] = ra;
	    tdec[jnum] = dec;
	    tmag[jnum] = mag;
	    tpeak[jnum] = peak;
	    nstar++;
	    if (nlog == 1)
		fprintf (stderr,"TABRNUM: %11.6f: %9.5f %9.5f %s %5.2f %d    \n",
			 num,ra,dec,cstr,mag,peak);

	    /* End of accepted star processing */
	    }

	/* Log operation */
	if (nlog > 0 && jnum%nlog == 0)
	    fprintf (stderr,"TABRNUM: %5d / %5d / %5d sources catalog %s\r",
		     nstar,jnum,nstars,tabcat);

	/* End of star loop */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"TABRNUM: Catalog %s : %d / %d found\n",
		 tabcat,nstar,nstars);

    tabclose();
    return (nstars);
}


#define TABMAX 64

/* TABRKEY -- Read single keyword from tab table stars with specified numbers */

int
tabrkey (tabcat, nnum, tnum, keyword, tval)

char	*tabcat;	/* Name of reference star catalog file */
int	nnum;		/* Number of stars to look for */
double	*tnum;		/* Array of star numbers to look for */
char	*keyword;	/* Keyword for which to return values */
char	**tval;		/* Returned values for specified keyword */
{
    int jnum, lval;
    int nstar;
    int istar;
    double num;
    char *line;
    char *tvalue;
    char value[TABMAX];

    nstar = 0;
    if ((nstars = tabopen (tabcat)) <= 0) {
	fprintf (stderr,"TABRNUM: Cannot read catalog %s\n", tabcat);
	return (0);
	}

    /* Loop through star list */
    for (jnum = 0; jnum < nnum; jnum++) {
	line = tabdata;

	/* Loop through catalog to star */
	for (istar = 1; istar <= nstars; istar++) {
	    line = tabstar (istar, line);
	    if (line == NULL) {
		fprintf (stderr,"TABRNUM: Cannot read star %d\n", istar);
		break;
		}
	    num = tabgetr8 (line, entid);	/* ID number */
	    if (num == 0.0)
		num = (double) istar;
	    if (num == tnum[jnum])
		break;
	    }

	/* If star has been found in table */
	if (num == tnum[jnum]) {
	    nstar++;

	    /* Extract selected field */
	    (void) tabgetk (line, keyword, value, TABMAX);
	    lval = strlen (value);
	    if (lval > 0) {
		tvalue = (char *) calloc (1, lval+1);
		strcpy (tvalue, value);
		}
	    else
		tvalue = NULL;
	    tval[jnum] = tvalue;
	    }
	}

    tabclose();
    return (nstars);
}

static char newline = 10;
static char tab = 9;
static char *tabhead;
static char *tabbuff;
static int nentry = 0;


/* TABOPEN -- Open tab table catalog, returning number of entries */

int
tabopen (tabfile)

char *tabfile;	/* Tab table catalog file name */
{
    FILE *fcat;
    char *headbuff;
    int nr, lfile, ientry;
    char *headlast, *tabnew, *tabline, *lastline;
    char headend[4];

    nentry = 0;
    
/* Find length of tab table catalog */
    lfile = tabsize (tabfile);
    if (lfile < 1) {
	fprintf (stderr,"TABOPEN: Tab table catalog %s has no entries\n",tabfile);
	return (0);
	}

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
	headbuff = tabhead;
	if (tabparse (headbuff)) {
	    if (!(entid = tabcol ("ID")))
		entid = tabcol ("id");
	    if (!(entra = tabcol ("RA")))
		entra = tabcol ("ra");
	    if (!(entdec = tabcol ("DEC")))
		entdec = tabcol ("dec");
	    if (!(entmag = tabcol ("MAG")))
		entmag = tabcol ("mag");
	    if (!(entpeak = tabcol ("PEAK")))
		entpeak = tabcol ("peak");
	    }
	else {
	    fprintf (stderr,"TABOPEN: No columns in tab table %s\n",tabfile);
	    return (0);
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
    iline = 1;
    return (nstars);
}


void
tabclose()
{
    free (tabbuff);
    return;
}


/* TABSTAR -- Get tab table catalog entry for one star; return 0 if successful */

char *
tabstar (istar, line)

int istar;	/* Star sequence number in tab table catalog */
char *line;	/* Pointer to istar'th entry (updated as returned value) */
{
    char *nextline = line;

    if (istar > nstars) {
	fprintf (stderr, "TABSTAR:  %d is not in catalog\n",istar);
	return (NULL);
	}
    else if (istar < 1 && line) {
	nextline = strchr (nextline, newline) + 1;
	}
    else if (istar < iline || iline < 1) {
	iline = 1;
	nextline = tabdata;
	while (iline < istar) {
	    nextline = strchr (nextline, newline) + 1;
	    line = nextline;
	    iline ++;
	    }
	}
    else {
	while (iline < istar) {
	    nextline = strchr (nextline, newline) + 1;
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


/* TABGETK -- returns a character entry from tab table line for named column */

int
tabgetk (line, keyword, string, maxchar)

char	*line;		/* tab table line */
char	*keyword;	/* column header of desired value */
char	*string;	/* character string (returned) */
int	maxchar;	/* Maximum number of characters in returned string */
{
    int ientry = tabcol (keyword);

    return (tabgetc (line, ientry, string, maxchar));
}


/* TABGETC -- returns n'th entry from tab table line as character string */

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
    if (ientry > nentry || ientry < 1)
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

#define MAXCOL	50
static char *colname[MAXCOL];	/* Column headers */
static int lcol[MAXCOL];

/* Make a table of column headings */
int
tabparse (headbuff)

char *headbuff;		/* Line containing column headings */

{
    char *colhead;	/* Column heading first character */
    char *endcol;	/* Column heading last character */
    char *headlast;
    char *nextab;
    int ientry = 0;

    headlast = strchr (headbuff, newline);
    colhead = headbuff;
    while (colhead) {
	nextab = strchr (colhead, tab);
	if (nextab < headlast)
	    endcol = nextab - 1;
	else
	    endcol = headlast - 1;
	while (*endcol == ' ')
	    endcol = endcol - 1;
	lcol[ientry] = (int) (endcol - colhead) + 1;
	colname[ientry] = colhead;
	ientry++;
	colhead = nextab + 1;
	if (colhead > headlast)
	    break;
	}
    nentry = ientry;
    return (ientry);
}


/* Search table of column headings for a particlar entry (case-dependent) */

int
tabcol (keyword)

char	*keyword;	/* Column heading to find */

{
    int i;

    for (i = 0; i < nentry; i++) {
	if (!strncmp (keyword, colname[i], lcol[i])) {
	    return (i + 1);
	    }
	}
    return (0);
}


/* TABSIZE -- return size of file in bytes */

int
tabsize (filename)

char	*filename;	/* Name of file for which to find size */
{
    FILE *diskfile;
    long filesize;
    long position;

    /* Open file */
    if ((diskfile = fopen (filename, "r")) == NULL)
	return (-1);

    /* Move to end of the file */
    if (fseek (diskfile, 0, 2) == 0)

 	/* Position is the size of the file */
	filesize = ftell (diskfile);

    else
	filesize = -1;

    fclose (diskfile);

    return (filesize);
}


/* ISTAB -- Return 1 if tab table file, else 0 */

int
istab (filename)

char    *filename;      /* Name of file to check */
{
    FILE *diskfile;
    char line[81];
    char *endline;
    int nbr;

    /* First check file extension */
    if (strsrch (filename, ".tab"))
	return (1);

    /* If no .tab file extension, try opening the file */
    else {
	if ((diskfile = fopen (filename, "r")) == NULL)
	    return (0);
	else {
	    nbr = fread (line, 1, 80, diskfile);
	    if (nbr < 8)
		return (0);
	    line[nbr] = (char) 0;
	    if ((endline = strchr (line, newline)) == NULL)
		return (0);
	    else
		*endline = (char) 0 ;
	    fclose (diskfile);
	    if (strchr (line, tab))
		return (1);
	    else
		return (0);
	    }
	}
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
 *
 * Mar 20 1997	Clean up code in TABRNUM
 * May  7 1997	Set entry number to zero if column not found
 * May 29 1997	Add TABPARSE and TABCOL to more easily extract specific columns
 * May 29 1997	Add TABCLOSE to free memory from outside this file
 * Jun  4 1997	Set ID to sequence number in table if no ID/id entry present
 *
 * Jun  2 1998	Fix bug parsing last column of header
 * Jun 24 1998	Add string lengths to ra2str() and dec2str() calls
 * Sep 16 1998	Use limiting radius correctly
 * Sep 22 1998	Convert to output coordinate system
 * Oct 15 1998	Add tabsize() and istab()
 * Oct 21 1998	Add tabrkey() to read values of keyword for list of stars
 * Oct 29 1998	Correctly assign numbers when too many stars are found
 * Oct 30 1998	Fix istab() to check only first line of file
 * Dec  8 1998	Do not declare tabsize() static
 */
