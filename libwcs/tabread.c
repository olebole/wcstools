/*** File libwcs/tabread.c
 *** August 25, 1999
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

/* int tabread()	Read tab table stars in specified region
 * int tabrnum()	Read tab table stars with specified numbers
 * int tabrkey()	Read single keyword from specified tab table stars
 * struct StarCat tabcatopen()	Open tab table catalog, return number of entries
 * int tabopen()	Open tab table, returning number of entries
 * char *tabline()	Get tab table entry for one line
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
#include "wcscat.h"

static int nstars;	/* Number of stars in catalog */
#define ABS(a) ((a) < 0 ? (-(a)) : (a))

int tabparse();
static int tabgeti4();
static double tabgetra();
static double tabgetdec();
static double tabgetr8();

static int iline = 0;	/* Current line of data */

static char *tabdata;


/* TABREAD -- Read tab table stars in specified region */

int
tabread (tabcatname,cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,mag2,
	 nstarmax,tnum,tra,tdec,tmag,tpeak,nlog)

char	*tabcatname;	/* Name of reference star catalog file */
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
    int sysref;		/* Catalog coordinate system */
    double eqref;	/* Catalog equinox */
    double epref;	/* Catalog epoch */
    char cstr[32];
    struct TabTable *startab;
    struct StarCat *starcat;
    struct Star *star;

    int wrap;
    int jstar;
    int nstar;
    double ra,dec, rapm, decpm;
    double mag;
    double num;
    int peak, i;
    int istar, nstars;
    int verbose;
    char *line;

    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    SearchLim (cra, cdec, dra, ddec, &ra1, &ra2, &dec1, &dec2, verbose);

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

    starcat = tabcatopen (tabcatname);
    if (starcat == NULL || starcat->nstars <= 0) {
	fprintf (stderr,"TABRNUM: Cannot read catalog %s\n", tabcatname);
	return (0);
	}

    nstars = starcat->nstars;
    jstar = 0;

    /* Set catalog coordinate system */
    if (starcat->equinox != 0.0)
	eqref = starcat->equinox;
    else
	eqref = eqout;
    if (starcat->epoch != 0.0)
	epref = starcat->epoch;
    else
	epref = epout;
    if (starcat->coorsys)
	sysref = starcat->coorsys;
    else
	sysref = sysout;
    wcscstr (cstr, sysout, eqout, epout);
    star = (struct Star *) calloc (1, sizeof (struct Star));

    /* Loop through catalog */
    for (istar = 1; istar <= nstars; istar++) {

	/* Read position of next star */
	if (tabstar (istar, starcat, star)) {
	    fprintf (stderr,"TABREAD: Cannot read star %d\n", istar);
	    break;
	    }

	/* Set coordinate system for this star */
	sysref = star->coorsys;
	eqref = star->equinox;
	epref = star->epoch;

	/* Extract selected fields  */
	num = star->num;
	ra = star->ra;
	dec = star->dec;
	rapm = star->rapm;
	decpm = star->decpm;
	if (starcat->mprop)
	    wcsconp (sysref, sysout, eqref, eqout, epref, epout,
		     &ra, &dec, &rapm, &decpm);
	else
	    wcscon (sysref, sysout, eqref, eqout, &ra, &dec, epout);
	mag = star->xmag[0];
	peak = star->peak;
	if (drad > 0)
	    dist = wcsdist (cra,cdec,ra,dec);
	else
	    dist = 0.0;

	/* Check magnitude and position limits */
	if ((mag1 == mag2 || (mag >= mag1 && mag <= mag2)) &&
	    ((wrap && (ra >= ra1 || ra <= ra2)) ||
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

	    /* If radial search & too many stars, replace furthest star */
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

	    /* Otherwise if too many stars, replace faintest star */
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
			jstar,istar,nstars,tabcatname);

	/* End of star loop */
	}

    /* Summarize search */
    if (nlog > 0) {
	fprintf (stderr,"TABREAD: Catalog %s : %d / %d / %d found\n",tabcatname,
		 jstar,istar,nstars);
	if (nstar > nstarmax)
	    fprintf (stderr,"TABREAD: %d stars found; only %d returned\n",
		     nstar,nstarmax);
	}

    tabcatclose (starcat);
    free ((char *)tdist);
    return (nstar);
}


/* TABRNUM -- Read tab table stars with specified numbers */

int
tabrnum (tabcatname, nnum, sysout, eqout, epout, tnum,tra,tdec,tmag,tpeak,nlog)

char	*tabcatname;	/* Name of reference star catalog file */
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
    int istar, nstars;
    char *line;
    int sysref;		/* Catalog coordinate system */
    double eqref;	/* Catalog equinox */
    double epref;	/* Catalog epoch */
    char cstr[32];
    struct TabTable *startab;
    struct StarCat *starcat;
    struct Star *star;

    line = 0;

    nstar = 0;
    starcat = tabcatopen (tabcatname);
    if (starcat == NULL || starcat->nstars <= 0) {
	fprintf (stderr,"TABRNUM: Cannot read catalog %s\n", tabcatname);
	return (0);
	}
    startab = starcat->startab;
    nstars = starcat->nstars;

    /* Set catalog coordinate system */
    if (starcat->equinox != 0.0)
	eqref = starcat->equinox;
    else
	eqref = eqout;
    if (starcat->epoch != 0.0)
	epref = starcat->epoch;
    else
	epref = epout;
    if (starcat->coorsys)
	sysref = starcat->coorsys;
    else
	sysref = sysout;

    wcscstr (cstr, sysout, eqout, epout);

    /* Allocate catalog entry buffer */
    star = (struct Star *) calloc (1, sizeof (struct Star));
    star->num = 0.0;

    /* Loop through star list */
    line = startab->tabdata;
    for (jnum = 0; jnum < nnum; jnum++) {

	/* Loop through catalog to star */
	for (istar = 1; istar <= nstars; istar++) {
	    if ((line = tabline (startab, istar)) == NULL) {
		fprintf (stderr,"TABRNUM: Cannot read star %d\n", istar);
		break;
		}

	    /* Check ID number first */
	    if ((num = tabgetr8 (startab,line,starcat->entid)) == 0.0)
		num = (double) istar;
	    if (num == tnum[jnum])
		break;
	    }

	/* If star has been found in table, read rest of entry */
	if (num == tnum[jnum]) {
	    starcat->istar = startab->iline;
	    tabstar (istar, starcat, star);

	    /* Set coordinate system for this star */
	    sysref = star->coorsys;
	    eqref = star->equinox;
	    epref = star->epoch;

	    /* Extract selected fields  */
	    ra = star->ra;
	    dec = star->dec;
	    wcscon (sysref, sysout, eqref, eqout, &ra, &dec, epout);
	    mag = star->xmag[0];
	    peak = star->peak;

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
		     nstar,jnum,nstars,tabcatname);

	/* End of star loop */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"TABRNUM: Catalog %s : %d / %d found\n",
		 tabcatname,nstar,nstars);

    tabcatclose (starcat);
    return (nstar);
}


#define TABMAX 64

/* TABRKEY -- Read single keyword from tab table stars with specified numbers */

int
tabrkey (tabcatname, nnum, tnum, keyword, tval)

char	*tabcatname;	/* Name of reference star catalog file */
int	nnum;		/* Number of stars to look for */
double	*tnum;		/* Array of star numbers to look for */
char	*keyword;	/* Keyword for which to return values */
char	**tval;		/* Returned values for specified keyword */
{
    int jnum, lval;
    int nstar;
    int istar, nstars;
    double num;
    char *line;
    char *tvalue;
    char value[TABMAX];
    struct TabTable *startab;
    struct StarCat *starcat;

    nstar = 0;
    starcat = tabcatopen (tabcatname);
    startab = starcat->startab;
    if (startab == NULL || startab->nlines <= 0) {
	fprintf (stderr,"TABRKEY: Cannot read catalog %s\n", tabcatname);
	return (0);
	}

    /* Loop through star list */
    nstars = startab->nlines;
    for (jnum = 0; jnum < nnum; jnum++) {

	/* Loop through catalog to star */
	for (istar = 1; istar <= nstars; istar++) {
	    if ((line = tabline (startab, istar)) == NULL) {
		fprintf (stderr,"TABRKEY: Cannot read star %d\n", istar);
		num = 0.0;
		break;
		}

	    /* Check ID number */
	    if ((num = tabgetr8 (startab,line,starcat->entid)) == 0.0)
		num = (double) istar;
	    if (num == tnum[jnum])
		break;
	    }

	/* If star has been found in table */
	if (num == tnum[jnum]) {
	    nstar++;

	    /* Extract selected field */
	    (void) tabgetk (startab, line, keyword, value, TABMAX);
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

    tabclose (startab);
    return (nstars);
}

static char newline = 10;
static char tab = 9;


/* TABCATOPEN -- Open tab table catalog, returning number of entries */

struct StarCat *
tabcatopen (tabfile)

char *tabfile;	/* Tab table catalog file name */
{
    char cstr[32];
    struct TabTable *startab;
    struct StarCat *sc;

    /* Allocate catalog data structure */
    sc = (struct StarCat *) calloc (1, sizeof (struct StarCat));

    /* Open the tab table file */
    if ((startab = tabopen (tabfile)) == NULL)
	return (NULL);
    sc->startab = startab;

    /* Extract positions of keywords we will want to use */
    if (!(sc->entid = tabcol (startab, "ID")))
	sc->entid = tabcol (startab, "id");
    if (!(sc->entra = tabcol (startab, "RA")))
	sc->entra = tabcol (startab, "ra");
    if (!(sc->entdec = tabcol (startab, "DEC")))
	sc->entdec = tabcol (startab, "dec");
    if (!(sc->entmag = tabcol (startab, "MAG")))
	sc->entmag = tabcol (startab, "mag");
    if (!(sc->entpeak = tabcol (startab, "PEAK")))
	sc->entpeak = tabcol (startab, "peak");

    /* Set catalog coordinate system */
    sc->coorsys = 0;
    sc->equinox = 0.0;
    sc->epoch = 0.0;
    if (tabhgetc (startab,"RADECSYS", cstr)) {
	sc->coorsys = wcscsys (cstr);
	if (!tabhgetr8 (startab,"EQUINOX", &sc->equinox))
	    sc->equinox = wcsceq (cstr);
	if (!tabhgetr8 (startab,"EPOCH",&sc->epoch))
	    sc->epoch = sc->equinox;
	}
    else if (tabhgetr8 (startab,"EQUINOX", &sc->equinox)) {
	if (!tabhgetr8 (startab,"EPOCH",&sc->epoch))
	    sc->epoch = sc->equinox;
	if (sc->equinox = 1950.0)
	    sc->coorsys = WCS_B1950;
	else
	    sc->coorsys = WCS_J2000;
	}
    else if (tabhgetr8 (startab,"EPOCH", &sc->epoch)) {
	sc->equinox = sc->epoch;
	if (sc->equinox = 1950.0)
	    sc->coorsys = WCS_B1950;
	else
	    sc->coorsys = WCS_J2000;
	}

    /* Set other stuff */
    sc->rasorted = 0;
    sc->nstars = startab->nlines;
    if (sc->entmag)
	sc->nmag = 1;
    else
	sc->nmag = 0;
    sc->mprop = 0;
    if (sc->entid)
	sc->stnum = 1;
    else
	sc->stnum = 0;

    iline = 1;
    return (sc);
}


/* TABCATCLOSE -- Close tab table catalog and free associated data structures */

void
tabcatclose (sc)
    struct StarCat *sc;
{
    tabclose (sc->startab);
    free (sc);
    return;
}


/* TABSTAR -- Get tab table catalog entry for one star;
   return NULL if unsuccessful */

int
tabstar (istar, sc, st)

int istar;      /* Star sequence number in tab table catalog */
struct StarCat *sc; /* Star catalog data structure */
struct Star *st; /* Star data structure, updated on return */
{
    struct TabTable *startab = sc->startab;
    char *line;

    if ((line = tabline (startab, istar)) == NULL) {
	fprintf (stderr,"TABSTAR: Cannot read star %d\n", istar);
	return (-1);
	}

    /* Extract selected fields  */
    if (sc->entid)
	st->num = tabgetr8 (startab,line,sc->entid);
    else
	st->num = (double) istar;
    st->ra = tabgetra (startab, line, sc->entra);       /* Right ascension */
    st->dec = tabgetdec (startab, line, sc->entdec);    /* Declination */
    st->xmag[0] = tabgetr8 (startab, line, sc->entmag);     /* Magnitude */
    st->peak = tabgeti4 (startab, line, sc->entpeak);   /* Peak counts */

    st->rapm = 0.0;
    st->decpm = 0.0;
    st->coorsys = sc->coorsys;
    st->equinox = sc->equinox;
    st->epoch = sc->epoch;
    return (0);
}


/* TABOPEN -- Open tab table file, returning number of entries */

struct TabTable *
tabopen (tabfile)

char *tabfile;	/* Tab table catalog file name */
{
    FILE *fcat;
    char *headbuff;
    int nr, lfile, ientry, lfname;
    char *headlast, *tabnew, *tabline, *lastline;
    char headend[4];
    struct TabTable *tabtable;

    if (!strcmp (tabfile, "stdin")) {
	lfile = 100000;
	fcat = stdin;
	}
    else {

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
	}

    /* Allocate tab table structure */
    if ((tabtable=(struct TabTable *) malloc(sizeof(struct TabTable))) == NULL){
	fprintf (stderr,"TABOPEN: cannot allocate Tab Table structure\n");
	return (0);
	}

    /* Allocate space in structure for filename and save it */
    lfname = strlen (tabfile) + 1;
    if ((tabtable->filename = malloc (lfname)) == NULL) {
	fprintf (stderr,"TABOPEN: cannot allocate filename in sructure\n");
	return (0);
	}
    strncpy (tabtable->filename, tabfile, lfname);

    /* Allocate buffer to hold entire catalog and read it */
    if ((tabtable->tabbuff = malloc (lfile)) != NULL) {
	nr = fread (tabtable->tabbuff, 1, lfile, fcat);
	if (fcat != stdin && nr < lfile) {
	    fprintf (stderr,"TABOPEN: read only %d / %d bytes of file %s\n",
		     nr, lfile, tabfile);
	    (void) fclose (fcat);
	    return (0);
	    }
	tabtable->tabhead = tabtable->tabbuff;
	headend[0] = '-';
	headend[1] = '-';
	headend[2] = 0;
	tabline = tabtable->tabbuff;
	while (strncmp (tabline,headend, 2)) {
	    lastline = tabline;
	    tabline = strchr (tabline,newline) + 1;
	    }
	tabtable->tabhead = lastline;
	tabtable->tabdata = strchr (tabline, newline) + 1;

	/* Extract positions of keywords we will want to use */
	ientry = 0;
	if (!tabparse (tabtable)) {
	    fprintf (stderr,"TABOPEN: No columns in tab table %s\n",tabfile);
	    return (0);
	    }

    /* Enumerate entries in tab table catalog by counting newlines */
	tabnew = strchr (tabtable->tabdata, newline) + 1;
	tabnew = tabtable->tabdata;
	tabtable->nlines = 0;
	while ((tabnew = strchr (tabnew, newline)) != NULL) {
	    tabnew = tabnew + 1;
	    tabtable->nlines = tabtable->nlines + 1;
	    }
	}

    (void) fclose (fcat);
    tabtable->tabline = tabtable->tabdata;
    tabtable->iline = 1;
    return (tabtable);
}


void
tabclose (tabtable)

    struct TabTable *tabtable;
{
    if (tabtable->filename != NULL) free (tabtable->filename);
    if (tabtable->tabbuff != NULL) free (tabtable->tabbuff);
    if (tabtable->colname != NULL) free (tabtable->colname);
    if (tabtable->lcol != NULL) free (tabtable->lcol);
    if (tabtable->lcfld != NULL) free (tabtable->lcfld);
    free (tabtable);
    return;
}


/* TABLINE -- Get tab table entry for one line;
	      return NULL if unsuccessful */

char *
tabline (tabtable, iline)

struct TabTable *tabtable;	/* Tab table structure */
int iline;	/* Line sequence number in tab table */
{
    char *nextline = tabtable->tabline;

    /* Return NULL if tab table has not been opened */
    if (tabtable == NULL)
	return (NULL);

    /* Return NULL if trying to read past last line */
    if (iline > tabtable->nlines) {
	fprintf (stderr, "TABLINE:  line %d is not in table\n",iline);
	return (NULL);
	}

    /* If iline is 0 or less, just read next line from table */
    else if (iline < 1 && nextline) {
	tabtable->iline++;
	if (tabtable->iline > tabtable->nlines) {
	    fprintf (stderr, "TABLINE:  line %d is not in table\n",iline);
	    return (NULL);
	    }
	nextline = strchr (nextline, newline) + 1;
	}

    /* If iline is before current line, read from start of file */
    else if (iline < tabtable->iline) {
	tabtable->iline = 1;
	tabtable->tabline = tabtable->tabdata;
	while (tabtable->iline < iline) {
	    tabtable->tabline = strchr (tabtable->tabline, newline) + 1;
	    tabtable->iline ++;
	    }
	}
    /* If iline is after current line, read forward */
    else if (iline > tabtable->iline) {
	while (tabtable->iline < iline) {
	    tabtable->tabline = strchr (tabtable->tabline, newline) + 1;
	    tabtable->iline ++;
	    }
	}

    return (tabtable->tabline);
}


/* TABGETRA -- returns double right ascension in degrees */

static double
tabgetra (startab, line, ientry)

struct TabTable *startab;	/* Tab table structure */
char	*line;	/* tab table line */
int	ientry;	/* sequence of entry on line */
{
    char str[24];

    if (tabgetc (startab, line, ientry, str, 24))
	return (0.0);
    else
	return (str2ra (str));
}


/* TABGETDEC -- returns double declination in degrees */

static double
tabgetdec (startab, line, ientry)

struct TabTable *startab;	/* Tab table structure */
char	*line;			/* tab table line */
int	ientry;			/* sequence of entry on line */
{
    char str[24];

    if (tabgetc (startab, line, ientry, str, 24))
	return (0.0);
    else
	return (str2dec (str));
}

/* TABGETR8 -- returns 8-byte floating point number from tab table line */

static double
tabgetr8 (tabtable, line, ientry)

struct TabTable *tabtable;	/* Tab table structure */
char	*line;	/* tab table line */
int	ientry;	/* sequence of entry on line */
{
    char str[24];

    if (tabgetc (tabtable, line, ientry, str, 24))
	return (0.0);
    else
        return (atof (str));
}


/* TABGETI4 -- returns a 4-byte integer from tab table line */

static int
tabgeti4 (tabtable, line, ientry)

struct TabTable *tabtable;	/* Tab table structure */
char	*line;	/* tab table line */
int	ientry;	/* sequence of entry on line */
{
    char str[24];

    if (tabgetc (tabtable, line, ientry, str, 24))
	return (0);
    else
        return ((int) atof (str));
}


/* TABGETK -- returns a character entry from tab table line for named column */

int
tabgetk (tabtable, line, keyword, string, maxchar)

struct TabTable *tabtable;	/* Tab table structure */
char	*line;			/* tab table line */
char	*keyword;		/* column header of desired value */
char	*string;		/* character string (returned) */
int	maxchar;	/* Maximum number of characters in returned string */
{
    int ientry = tabcol (tabtable, keyword);

    return (tabgetc (tabtable, line, ientry, string, maxchar));
}


/* TABGETC -- returns n'th entry from tab table line as character string */

int
tabgetc (tabtable, line, ientry, string, maxchar)

struct TabTable *tabtable;	/* Tab table structure */
char	*line;			/* tab table line */
int	ientry;			/* sequence of entry on line */
char	*string;		/* character string (returned) */
int	maxchar;	/* Maximum number of characters in returned string */
{
    char *entry, *nextab;
    int ient, ncstr;

    ient = 1;
    entry = line;
    if (ientry > tabtable->ncols || ientry < 1)
	return (-1);
    for (ient  = 1; ient <= ientry; ient ++) {

    /* End ient'th entry with tab, newline, or end of string */
	if (ient < tabtable->ncols) 
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


/* TABHGETR8 -- read an 8-byte floating point number from a tab table header */

int
tabhgetr8 (tabtable, keyword, result)

struct TabTable *tabtable;	/* Tab table structure */
char	*keyword;		/* sequence of entry on line */
double	*result;
{
    char value[24];

    if (tabhgetc (tabtable, keyword, value)) {
	*result = atof (value);
	return (1);
	}
    else
	return (0);
}


/* TABHGETI4 -- read a 4-byte integer from a tab table header */

int
tabhgeti4 (tabtable, keyword, result)

struct TabTable *tabtable;	/* Tab table structure */
char	*keyword;		/* sequence of entry on line */
int	*result;
{
    char value[24];

    if (tabhgetc (tabtable, keyword, value)) {
	*result  = (int) atof (value);
	return (1);
	}
    else
	return (0);
}


/* TABHGETC -- read a string from a tab table header */

int
tabhgetc (tabtable, keyword, result)

struct TabTable *tabtable;	/* Tab table structure */
char	*keyword;		/* sequence of entry on line */
char	*result;
{
    char *str0, *str1, *line, *head, str[24], keylow[24], keyup[24];
    int ncstr, lkey, i;

    head = tabtable->tabbuff;
    str0 = 0;

    /* Make all-upper-case and all-lower-case versions of keyword */
    lkey = strlen (keyword);
    for (i = 0; i < lkey; i++) {
	if (keyword[i] > 96 && keyword[i] < 123)
	    keyup[i] = keyword[i] - 32;
	else
	    keyup[i] = keyword[i];
	if (keyword[i] > 64 && keyword[i] < 91)
	    keylow[i] = keyword[i] + 32;
	else
	    keylow[i] = keyword[i];
	}

    /* Find keyword or all-upper-case or all-lower-case version in header */
    while (head < tabtable->tabhead) {
	line = strsrch (head, keyword);
	if (line == NULL)
	    line = strsrch (head, keylow);
	if (line == NULL)
	    line = strsrch (head, keyup);
	if (line == NULL)
	    break;
	if (line == tabtable->tabhead || line[-1] == newline) {
	    str0 = strchr (line, tab) + 1;
	    str1 = strchr (line, newline);
	    break;
	    }
	else
	    head = line;
	}

    /* Return value as a character string and 1 if found */
    if (str0) {
	ncstr = str1 - str0 - 1;
	strncpy (result, str0, ncstr);
	return (1);
	}
    else
	return (0);
}


/* TABPARSE -- Make a table of column headings */

int
tabparse (tabtable)

struct TabTable *tabtable;	/* Tab table structure */
{
    char *colhead;	/* Column heading first character */
    char *endcol;	/* Column heading last character */
    char *headlast;
    char *hyphens;
    char *hyphlast;
    char *nextab;
    int i, ientry = 0;
    int nbytes;

    /* Return if no column names in header */
    headlast = strchr (tabtable->tabhead, newline);
    if (headlast == tabtable->tabhead)
	return (0);

    /* Count columns in table header */
    tabtable->ncols = 1;
    for (colhead = tabtable->tabhead; colhead < headlast; colhead++) {
	if (*colhead == tab)
	    tabtable->ncols++;
	}

    /* Tabulate column names */
    nbytes = tabtable->ncols * sizeof (char *);
    tabtable->colname = (char **)malloc (nbytes);
    nbytes = tabtable->ncols * sizeof (int);
    tabtable->lcol = malloc (nbytes);
    colhead = tabtable->tabhead;
    while (colhead) {
	nextab = strchr (colhead, tab);
	if (nextab < headlast)
	    endcol = nextab - 1;
	else
	    endcol = headlast - 1;
	while (*endcol == ' ')
	    endcol = endcol - 1;
	tabtable->lcol[ientry] = (int) (endcol - colhead) + 1;
	tabtable->colname[ientry] = colhead;
	ientry++;
	colhead = nextab + 1;
	if (colhead > headlast)
	    break;
	}

    /* Tabulate field widths */
    hyphens = headlast + 1;
    hyphlast = strchr (hyphens, newline);
    if (hyphlast == hyphens)
	return (0);
    tabtable->lcfld = malloc (nbytes);
    colhead = hyphens;
    i = 0;
    while (colhead) {
	nextab = strchr (colhead, tab);
	if (nextab < hyphlast)
	    endcol = nextab - 1;
	else
	    endcol = headlast - 1;
	tabtable->lcfld[i] = (int) (endcol - colhead) + 1;
	i++;
	colhead = nextab + 1;
	if (colhead > hyphlast)
	    break;
	}

    return (ientry);
}


/* Search table of column headings for a particlar entry (case-dependent) */

int
tabcol (tabtable, keyword)

struct TabTable *tabtable;	/* Tab table structure */
char	*keyword;		/* Column heading to find */

{
    int i;

    for (i = 0; i < tabtable->ncols; i++) {
	if (!strncmp (keyword, tabtable->colname[i], tabtable->lcol[i])) {
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
    char line[1025];
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
	    nbr = fread (line, 1, 1024, diskfile);
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

 * Jan 20 1999	Add tabcatopen() and keep table info in structure, not global
 * Jan 25 1999	Add lcfld to structure to keep track of field widths
 * Jan 29 1999	Add current line number and pointer to table structure
 * Feb  2 1999	Allow for equinox other than 2000.0 in tab table header
 * Feb  2 1999	Add tabhgetc() to read char string values from tab table header
 * Feb 17 1999	Increase maximum line in istab() from 80 to 1024
 * Mar  2 1999	Fix bugs calling tabhgetx()
 * Mar  2 1999	Rewrite tabhgetx() to use tabhgetc() for all header reading
 * May 28 1999	Add tabcatopen() and tabstar() and use them
 * Jun  3 1999	Fix bug so header parameters are read correctly
 * Jun 16 1999	Use SearchLim()
 * Aug 16 1999	Fix bug to fix failure to search across 0:00 RA
 * Aug 25 1999  Return real number of stars from tabread()
 */
