/*** File libwcs/catread.c
 *** August 25, 1999
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

/* int catread()	Read ASCII catalog stars in specified region
 * int catrnum()	Read ASCII catalog stars with specified numbers
 * int catopen()	Open ASCII catalog, returning number of entries
 * char *catstar()	Get ASCII catalog entry for one star
 * int catsize()	Return length of file in bytes
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include "fitshead.h"
#include "wcs.h"
#include "wcscat.h"

#define MAXTOKENS 20	/* Maximum number of tokens to parse */
#define MAXWHITE 20	/* Maximum number of whitespace characters */
struct Tokens {
    char *line;		/* Line which has been parsed */
    int lline;		/* Number of characters in line */
    int ntok;		/* Number of tokens on line */
    int nwhite;		/* Number of whitespace characters */
    char white[MAXWHITE];	/* Whitespace (separator) characters */
    char *tok1[MAXTOKENS];	/* Pointers to start of tokens */
    int ltok[MAXTOKENS];	/* Lengths of tokens */
    int itok;		/* Current token number */
};

/* default pathname for catalog,  used if catalog file not found in current
   working directory, but overridden by WCS_CATDIR environment variable */
char catdir[64]="/data/catalogs";

int catstar();
static double cat2ra();
static double cat2dec();
double dt2ep();		/* Julian Date to epoch (fractional year) */


/* CATREAD -- Read ASCII stars in specified region */

int
catread (catfile, cra, cdec, dra, ddec, drad, sysout, eqout, epout,
	 mag1, mag2, nstarmax, tnum, tra, tdec, tmag, tobj, nlog)

char	*catfile;	/* Name of reference star catalog file */
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
char	**tobj;		/* Array of object names (returned) */
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
    struct StarCat *starcat;
    struct Star *star;
    char *objname;
    int lname;
    int wrap;
    int jstar;
    int nstar;
    double ra,dec,rapm,decpm;
    double mag;
    double num;
    int peak, i;
    int istar;
    int verbose;

    star = NULL;
    starcat = NULL;

    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    wcscstr (cstr, sysout, eqout, epout);

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

    /* Allocate catalog entry buffer */
    star = (struct Star *) calloc (1, sizeof (struct Star));
    star->num = 0.0;

/* Logging interval */
    nstar = 0;
    tdist = (double *) malloc (nstarmax * sizeof (double));

    /* Open catalog file */
    starcat = catopen (catfile);
    if (starcat == NULL)
	return (0);
    if (starcat->nstars <= 0) {
	free (starcat);
	if (star != NULL)
	    free (star);
	return (0);
	}

    jstar = 0;

    /* Loop through catalog */
    for (istar = 1; istar <= starcat->nstars; istar++) {
	if (catstar (istar, starcat, star)) {
	    fprintf (stderr,"CATREAD: Cannot read star %d\n", istar);
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
	peak = 0;
	if (drad > 0)
	    dist = wcsdist (cra,cdec,ra,dec);
	else
	    dist = 0.0;

	/* If catalog is RA-sorted, stop reading if past highest RA */
	if (starcat->rasorted && !wrap && ra > ra2)
	    break;

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
		tdist[nstar] = dist;
		if (tobj != NULL) {
		    lname = strlen (star->objname) + 1;
		    objname = (char *)calloc (lname, 1);
		    strcpy (objname, star->objname);
		    tobj[nstar] = objname;
		    }
		if (dist > maxdist) {
			maxdist = dist;
			farstar = nstar;
			}
		if (mag > faintmag) {
			faintmag = mag;
			faintstar = nstar;
			}
		}

	    /* If too many stars and radial search, replace furthest star */
	    else if (drad > 0 && dist < maxdist) {
		tnum[farstar] = num;
		tra[farstar] = ra;
		tdec[farstar] = dec;
		tmag[farstar] = mag;
		tdist[farstar] = dist;
		if (tobj != NULL) {
		    free (tobj[farstar]);
		    lname = strlen (star->objname) + 1;
		    objname = (char *)calloc (lname, 1);
		    strcpy (objname, star->objname);
		    tobj[farstar] = objname;
		    }
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
		tdist[faintstar] = dist;
		if (tobj != NULL) {
		    free (tobj[faintstar]);
		    lname = strlen (star->objname) + 1;
		    objname = (char *)calloc (lname, 1);
		    strcpy (objname, star->objname);
		    tobj[faintstar] = objname;
		    }
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
		fprintf (stderr,"CATREAD: %11.6f: %9.5f %9.5f %s %5.2f %d    \n",
			 num,ra,dec,cstr,mag,peak);

	    /* End of accepted star processing */
	    }

	/* Log operation */
	if (nlog > 0 && istar%nlog == 0)
	    fprintf (stderr,"CATREAD: %5d / %5d / %5d sources catalog %s\r",
		     jstar,istar,starcat->nstars,catfile);

	/* End of star loop */
	}

    /* Summarize search */
    if (nlog > 0) {
	fprintf (stderr,"CATREAD: Catalog %s : %d / %d / %d found\n",
		 catfile,jstar,istar,starcat->nstars);
	if (nstar > nstarmax)
	    fprintf (stderr,"CATREAD: %d stars found; only %d returned\n",
		     nstar,nstarmax);
	}

    catclose (starcat);
    free ((char *)tdist);
    free (star);
    return (nstar);
}


/* CATRNUM -- Read ASCII stars with specified numbers */

int
catrnum (catfile,nnum,sysout,eqout,epout,match,tnum,tra,tdec,tmag,tobj,nlog)

char	*catfile;	/* Name of reference star catalog file */
int	nnum;		/* Number of stars to look for */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
int	match;		/* 1 to match star number exactly, else sequence num.*/
double	*tnum;		/* Array of star numbers to look for */
double	*tra;		/* Array of right ascensions (returned) */
double	*tdec;		/* Array of declinations (returned) */
double	*tmag;		/* Array of magnitudes (returned) */
char	**tobj;		/* Array of object names (returned) */
int	nlog;
{
    int jnum;
    int nstar;
    double ra,dec;
    double rapm, decpm;
    double mag;
    double num;
    int peak;
    int istar;
    int sysref;		/* Catalog coordinate system */
    double eqref;	/* Catalog equinox */
    double epref;	/* Catalog epoch */
    char cstr[32];
    struct StarCat *starcat;
    struct Star *star;
    char *objname;
    int lname;
    int starfound;

    nstar = 0;
    if ((starcat = catopen (catfile)) == NULL) {
	fprintf (stderr,"CATRNUM: Cannot read catalog %s\n", catfile);
	return (0);
	}
    sysref = starcat->coorsys;
    eqref = starcat->equinox;
    epref = starcat->epoch;
    if (!sysout)
	sysout = sysref;
    if (!eqout)
	eqout = eqref;
    if (!epout)
	epout = epref;

    /* Allocate catalog entry buffer */
    star = (struct Star *) calloc (1, sizeof (struct Star));
    star->num = 0.0;

    /* Loop through star list */
    for (jnum = 0; jnum < nnum; jnum++) {

	/* Loop through catalog to star */
	starfound = 0;
	if (match && starcat->stnum > 0 && starcat->stnum < 5) {
	    for (istar = 1; istar <= starcat->nstars; istar++) {
		if (catstar (istar, starcat, star)) {
		    fprintf (stderr,"CATRNUM: Cannot read star %d\n", istar);
		    break;
		    }
		if (star->num == tnum[jnum]) {
		    starfound = 1;
		    break;
		    }
		}
	    }
	else {
	    istar = (int) (tnum[jnum] + 0.5);
	    if (catstar (istar, starcat, star)) {
		fprintf (stderr,"CATRNUM: Cannot read star %d\n", istar);
		continue;
		}
	    starfound = 1;
	    }

	/* If star has been found in catalog */
	if (starfound) {

	    /* Extract selected fields  */
	    ra = star->ra;
	    dec = star->dec;
	    rapm = star->rapm;
	    decpm = star->decpm;

	    /* Set coordinate system for this star */
	    sysref = star->coorsys;
	    eqref = star->equinox;
	    epref = star->epoch;
    
	    if (starcat->mprop)
		wcsconp (sysref, sysout, eqref, eqout, epref, epout,
			 &ra, &dec, &rapm, &decpm);
	    else
		wcscon (sysref, sysout, eqref, eqout, &ra, &dec, epout);
	    mag = star->xmag[0];
	    peak = 0;

	    /* Save star position and magnitude in table */
	    tnum[jnum] = star->num;
	    tra[jnum] = ra;
	    tdec[jnum] = dec;
	    tmag[jnum] = mag;
	    if (tobj != NULL) {
		lname = strlen (star->objname) + 1;
		objname = (char *)calloc (lname, 1);
		strcpy (objname, star->objname);
		tobj[nstar] = objname;
		}
	    nstar++;
	    if (nlog == 1)
		fprintf (stderr,"CATRNUM: %11.6f: %9.5f %9.5f %s %5.2f %d    \n",
			 num,ra,dec,cstr,mag,peak);

	    /* End of accepted star processing */
	    }

	/* Log operation */
	if (nlog > 0 && jnum%nlog == 0)
	    fprintf (stderr,"CATRNUM: %5d / %5d / %5d sources catalog %s\r",
		     nstar,jnum,starcat->nstars,catfile);

	/* End of star loop */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"CATRNUM: Catalog %s : %d / %d found\n",
		 catfile,nstar,starcat->nstars);

    catclose(starcat);
    free (star);
    return (nstar);
}


static char newline = 10;
static char cat = 9;
static char *cathead;

/* CATOPEN -- Open ASCII catalog, returning number of entries */

struct StarCat *
catopen (catfile)

char *catfile;	/* ASCII catalog file name */
{
    struct StarCat *sc;
    struct Tokens tokens;
    FILE *fcat;
    char *headbuff;
    char header[80];
    char catpath[128];	/* Full pathname for catalog file */
    char *str;
    int nr, lfile, ientry, lhead, ldesc;
    char *catnew, *catline, *lastline, *catdesc;
    char ctemp, *line, *linend, *cdot;
    char token[80];
    int ntok, itok;
    int ltok;

/* Find length of ASCII catalog */
    lfile = catsize (catfile);

    /* Check for existence of catalog */
    if (lfile < 2) {

	/* Prepend directory name file not in working directory */
	if ((str = getenv("WCS_CATDIR")) != NULL )
	    strcpy (catdir, str);
	strcpy (catpath, catdir);
	strcat (catpath, "/");
	strcat (catpath, catfile);
	lfile = catsize (catpath);
	if (lfile < 2) {
	    fprintf (stderr,"CATOPEN: ASCII catalog %s has no entries\n",catfile);
	    return (0);
	    }
	}
    else
	strcpy (catpath, catfile);

    /* Allocate catalog data structure */
    sc = (struct StarCat *) calloc (1, sizeof (struct StarCat));

    /* Open ASCII catalog */
    if (!(fcat = fopen (catpath, "r"))) {
	fprintf (stderr,"CATOPEN: ASCII catalog %s cannot be read\n",catfile);
	return (0);
	}

    /* Allocate buffer to hold entire catalog and read it */
    if ((sc->catbuff = malloc (lfile+2)) == NULL) {
	fprintf (stderr,"CATOPEN: Cannot allocate memory for ASCII catalog %s\n",
		 catfile);
	fclose (fcat);
	return (0);
	}
    sc->catbuff[lfile] = (char) 0;
    sc->catbuff[lfile+1] = (char) 0;

    /* Read entire catalog into memory at once */
    nr = fread (sc->catbuff, 1, lfile, fcat);
    if (nr < lfile) {
	fprintf (stderr,"CATOPEN: read only %d / %d bytes of file %s\n",
		 nr, lfile, catfile);
	(void) fclose (fcat);
	return (0);
	}
    cathead = sc->catbuff;

    /* Extract catalog information from first line */
    sc->inform = 'H';
    sc->coorsys = WCS_B1950;
    sc->epoch = 1950.0;
    sc->equinox = 1950.0;
    sc->nmag = 1;
    sc->mprop = 0;
    sc->rasorted = 0;
    sc->stnum = 1;
    sc->entepoch = 0;

    catdesc = strchr (sc->catbuff, newline) + 1;
    lhead = catdesc - sc->catbuff;
    if (lhead < 80) {
	strncpy (header, sc->catbuff, lhead);
	header[lhead] = (char) 0;
	}
    else {
	strncpy (header, sc->catbuff, 79);
	header[79] = (char) 0;
	}

    /* Catalog positions are in B1950 (FK4) coordinates */
    if (strsrch (header, "/b") || strsrch (header, "/B")) {
	sc->coorsys = WCS_B1950;
	sc->epoch = 1950.0;
	sc->equinox = 1950.0;
	}

    /* Catalog positions are in ecliptic coordinates */
    if (strsrch (header, "/e") || strsrch (header, "/E")) {
	sc->coorsys = WCS_ECLIPTIC;
	sc->inform = 'D';
	sc->epoch = 2000.0;
	sc->equinox = 2000.0;
	}

    /* Catalog positions are galactic coordinates */
    if (strsrch (header, "/g") || strsrch (header, "/G")) {
	sc->coorsys = WCS_GALACTIC;
	sc->inform = 'D';
	sc->epoch = 2000.0;
	sc->equinox = 2000.0;
	}

    /* Catalog positions are in hh.mmsssss dd.mmssss format */
    if (strsrch (header, "/h") || strsrch (header, "/H"))
	sc->inform = 'H';

    /* Catalog positions are J2000 (FK5) coordinates */
    if (strsrch (header, "/j") || strsrch (header, "/J")) {
	sc->coorsys = WCS_J2000;
	sc->epoch = 2000.0;
	sc->equinox = 2000.0;
	}
    if (strsrch (header, "/q") || strsrch (header, "/Q")) {
	sc->coorsys = 0;
	sc->epoch = 0.0;
	sc->equinox = 0.0;
	}

    /* No number in first column, RA or object name first */
    if (strsrch (header, "/n") || strsrch (header, "/N"))
	sc->stnum = 0;

    /* Object name instead of number in first column */
    if (strsrch (header, "/o") || strsrch (header, "/O"))
	sc->stnum = 5;

    /* No magnitude */
    if (strsrch (header, "/m") || strsrch (header, "/M"))
	sc->nmag = 0;
    else
	sc->nmag = 1;

    /* Proper motion */
    if (strsrch (header, "/p") || strsrch (header, "/P"))
	sc->mprop = 1;
    if (strsrch (header, "/r") || strsrch (header, "/R"))
	sc->rasorted = 1;
    if (strsrch (header, "/t") || strsrch (header, "/T"))
	sc->inform = 'T';
    if (strsrch (header, "/y") || strsrch (header, "/Y"))
	sc->nepoch = 1;
    else
	sc->nepoch = 0;

    /* Second line is description */
    sc->catdata = strchr (catdesc, newline) + 1;
    ldesc = sc->catdata - catdesc;
    if (ldesc < 64)
	strncpy (sc->isname, catdesc, ldesc);
    else
	strncpy (sc->isname, catdesc, 63);

    /* Enumerate entries in ASCII catalog by counting newlines */
    catnew = sc->catdata;
    sc->nstars = 0;
    while ((catnew = strchr (catnew, newline)) != NULL) {
	catnew = catnew + 1;
	sc->nstars = sc->nstars + 1;
	}
    sc->catline = sc->catdata;
    sc->catlast = sc->catdata + lfile;
    sc->istar = 1;

    /* Check number of decimal places in star number, if present */
    if (sc->stnum == 1) {

	/* Temporarily terminate line with 0 */
	line = sc->catline;
	linend = strchr (sc->catline, newline);
	if (linend == NULL)
	    linend = sc->catlast;
	ctemp = *linend;
	*linend = (char) 0;

	/* Extract information from line of catalog */
	ntok = setoken (&tokens, line, NULL);
	ltok = nextoken (&tokens, token);
	if (ltok > 0) {
	    sc->nndec = 0;
	    if ((cdot = strchr (token,'.')) != NULL) {
		while (*(++cdot) != (char)0)
		    sc->nndec++;
		}
	    }
	else
	    sc->nndec = 0;
	*linend = ctemp;
	}

    (void) fclose (fcat);
    return (sc);
}


/* CATCLOSE -- Close ASCII catalog and free associated data structures */

void
catclose (sc)
    struct StarCat *sc;
{
    free (sc->catbuff);
    free (sc);
    return;
}



/* CATSTAR -- Get ASCII catalog entry for one star; return NULL if unsuccessful */

int
catstar (istar, sc, st)

int istar;	/* Star sequence number in ASCII catalog */
struct StarCat *sc; /* Star catalog data structure */
struct Star *st; /* Star data structure, updated on return */
{
    struct Tokens tokens;
    double ydate, dj;
    char *line;
    char *nextline;
    char token[80];
    char ctemp, *linend;
    int ntok, itok;
    int ltok;
    int nmag, imag;
    double dno;

    /* Return error if requested number beyond catalog */
    if (istar > sc->nstars) {
	fprintf (stderr, "CATSTAR:  %d is not in catalog\n",istar);
	return (-1);
	}

    /* If star is 0, read next star in catalog */
    else if (istar < 1) {
	line = strchr (sc->catline, newline) + 1;
	if (line == NULL)
	    return (-1);
	else
	    sc->catline = line;
	}

    /* If star is before current star, read from start of catalog */
    else if (istar < sc->istar) {
	sc->istar = 1;
	sc->catline = sc->catdata;
	while (sc->istar < istar) {
	    nextline = strchr (sc->catline, newline) + 1;
	    if (nextline == NULL)
		return (-1);
	    sc->catline = nextline;
	    sc->istar++;
	    }
	}

    /* If star is after current star, read forward to it */
    else if (istar > sc->istar) {
	while (sc->istar < istar) {
	    nextline = strchr (sc->catline, newline) + 1;
	    if (nextline == NULL)
		return (-1);
	    sc->catline = nextline;
	    sc->istar++;
	    }
	}

    /* Temporarily terminate line with 0 */
    line = sc->catline;
    linend = strchr (sc->catline, newline);
    if (linend == NULL)
	linend = sc->catlast;
    ctemp = *linend;
    *linend = (char) 0;

    /* Extract information from line of catalog */
    ntok = setoken (&tokens, line, NULL);

    /* Source number */
    if (sc->stnum > 0 && sc->stnum < 5) {
	ltok = nextoken (&tokens, token);
	if (ltok > 0) {
	    st->num = atof (token);
	    if (dno - ((double) ((int) dno)) < 0.000001)
		sc->stnum = 1;
	    else
		sc->stnum = 2;
	    }
	}
    else
	st->num = (double) sc->istar;

    /* Object name, if at start of line */
    if (sc->stnum == 5) {
	ltok = nextoken (&tokens, token);
	if (ltok > 31) {
	    strncpy (st->objname, token, 31);
	    st->objname[31] = 0;
	    }
	else if (ltok > 0)
	    strcpy (st->objname, token);
	}

    /* Right ascension or longitude */
    ltok = nextoken (&tokens, token);
    if (ltok < 1)
	return (-1);

    /* Translate 3-token RA (hh mm ss.ss) */
    if (sc->inform == 'T') {
	int hr, mn;
	double sec;

	hr = (int) (atof (token) + 0.5);
	ltok = nextoken (&tokens, token);
	if (ltok < 1)
	    return (-1);
	mn = (int) (atof (token) + 0.5);
	ltok = nextoken (&tokens, token);
	if (ltok < 1)
	    return (-1);
	sec = atof (token);
	st->ra = 15.0 * ((double) hr + ((double) mn / 60.0) + (sec / 3600.0));
	}

    /* Translate single-token RA (hh:mm:ss.ss or hh.mmssss) */
    else
	st->ra = cat2ra (token);

    /* Declination or latitude */
    ltok = nextoken (&tokens, token);
    if (ltok < 1)
	return (-1);

    /* Translate 3-token Dec (hh mm ss.ss) */
    if (sc->inform == 'T') {
	int deg, min;
	double sec;

	deg = (int) (atof (token) + 0.5);
	ltok = nextoken (&tokens, token);
	if (ltok < 1)
	    return (-1);
	min = (int) (atof (token) + 0.5);
	ltok = nextoken (&tokens, token);
	if (ltok < 1)
	    return (-1);
	sec = atof (token);
	st->dec = (double) deg + ((double) min / 60.0) + (sec / 3600.0);
	}

    /* Translate single-token Dec (dd:mm:ss.ss or dd.mmssss) */
    else
	st->dec = cat2dec (token);

    /* Equinox, if not set by header flag */
    if (sc->coorsys == 0) {
	ltok = nextoken (&tokens, token);
	if (ltok < 1)
	    return (-1);
	st->coorsys = wcscsys (token);
	st->equinox = wcsceq (token);
	st->epoch = sc->epoch;
	if (st->epoch == 0.0)
	    st->epoch = st->equinox;
	}

    else {
	st->coorsys = sc->coorsys;
	st->equinox = sc->equinox;
	st->epoch = sc->epoch;
	}

    /* Magnitude, if present */
    if (sc->nmag > 0) {
	for (imag = 0; imag < sc->nmag; imag++) {
	    ltok = nextoken (&tokens, token);
	    if (ltok > 0)
		st->xmag[imag] = atof (token);
	    }
	}

    /* Proper motion, if present */
    if (sc->mprop) {
	ltok = nextoken (&tokens, token);
	if (ltok > 1)
	    st->rapm = atof (token) / 3600.0;
	ltok = nextoken (&tokens, token);
	if (ltok > 1)
	    st->decpm = atof (token) / 3600.0;
	}

    /* Epoch, if present */
    if (sc->nepoch) {
	ltok = nextoken (&tokens, token);
	ydate = atof (token);
	st->epoch = dt2ep (ydate, 12.0);
	}

    /* Object name */
    itok = tokens.itok;
    if (sc->stnum < 5 && itok < ntok) {
	itok = -(itok+1);
	ltok = getoken (&tokens, itok, token);
	if (ltok > 31) {
	    strncpy (st->objname, token, 31);
	    st->objname[31] = 0;
	    }
	else if (ltok > 0)
	    strcpy (st->objname, token);
	}

    *linend = ctemp;
    return (0);
}


/* CATSIZE -- return size of ASCII catalog file in bytes */

int
catsize (filename)

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


/* -- SETOKEN -- tokenize a string for easy decoding */

int
setoken (tokens, string, cwhite)

struct Tokens *tokens;	/* Token structure returned */
char	*string;	/* character string to tokenize */
char	*cwhite;	/* additional whitespace characters
			 * if = tab, disallow spaces and commas */
{
    int ntok;		/* number of tokens found (returned) */
    char squote,dquote,jch;
    char token[32];
    char *iq, *stri, *wtype, *str0;
    int i0,i,j,naddw;

    squote = (char) 39;
    dquote = (char) 34;
    if (string == NULL)
	return (0);

    /* Line is terminated by newline or NULL */
    if (strchr (string, newline))
	tokens->lline = strchr (string, newline) - string - 1;
    else
	tokens->lline = strlen (string);

    /* Save current line in structure */
    tokens->line = string;

    /* Add extra whitespace characters */
    if (cwhite == NULL)
	naddw = 0;
    else
	naddw = strlen (cwhite);

    /* if character is tab, allow only tabs and nulls as separators */
    if (naddw > 0 && !strncmp (cwhite, "tab", 3)) {
	tokens->white[0] = (char) 9;
	tokens->white[1] = (char) 0;
	tokens->nwhite = 2;
	}

    /* otherwise, allow spaces, tabs, commas, nulls, and cwhite */
    else {
	tokens->nwhite = 3 + naddw;;
	tokens->white[0] = ' ';
	tokens->white[1] = (char) 9;
	tokens->white[2] = ',';
	tokens->white[3] = (char) 0;
	if (tokens->nwhite > 20)
	    tokens->nwhite = 20;
	if (naddw > 0) {
	    i = 0;
	    for (j = 3; j < tokens->nwhite; j++) {
		tokens->white[j] = cwhite[i];
		i++;
		}
	    }
	}
    tokens->white[tokens->nwhite] = (char) 0;

    tokens->ntok = 0;
    tokens->itok = 0;
    iq = string - 1;
    for (i = 0; i < MAXTOKENS; i++) {
	tokens->tok1[i] = NULL;
	tokens->ltok[i] = 0;
	}

    /* Process string one character at a time */
    stri = string;
    str0 = string;
    while (stri < string+tokens->lline) {

	/* Keep stuff between quotes in one token */
	if (stri <= iq)
	    continue;
	jch = *stri;

	/* Handle quoted strings */
	if (jch == squote)
	    iq = strchr (stri+1, squote);
	else if (jch == dquote)
	    iq = strchr (stri+1, dquote);
	else
	    iq = stri;
	if (iq > stri) {
	    tokens->ntok = tokens->ntok + 1;
	    if (tokens->ntok < MAXTOKENS) return (MAXTOKENS);
	    tokens->tok1[tokens->ntok] = stri + 1;
	    tokens->ltok[tokens->ntok] = (iq - stri) - 1;
	    stri = iq + 1;
	    str0 = iq + 1;
	    continue;
	    }

	/* Search for unquoted tokens */
	wtype = strchr (tokens->white, jch);

	/* If this is one of the additional whitespace characters,
	 * pass as a separate token */
	if (wtype > tokens->white + 2) {

	    /* Terminate token before whitespace */
	    if (stri > str0) {
		tokens->ntok = tokens->ntok + 1;
		if (tokens->ntok > MAXTOKENS) return (MAXTOKENS);
		tokens->tok1[tokens->ntok] = str0;
		tokens->ltok[tokens->ntok] = stri - str0;
		}

	    /* Make whitespace character next token; start new one */
	    tokens->ntok = tokens->ntok + 1;
	    if (tokens->ntok < MAXTOKENS) return (MAXTOKENS);
	    tokens->tok1[tokens->ntok] = stri;
	    tokens->ltok[tokens->ntok] = 1;
	    stri++;
	    str0 = stri;
	    }

	/* Pass previous token if regular whitespace or NULL */
	else if (wtype != NULL || jch == (char) 0) {

	    /* Ignore leading whitespace */
	    if (stri == str0) {
		stri++;
		str0 = stri;
		}

	    /* terminate token before whitespace; start new one */
	    else {
		tokens->ntok = tokens->ntok + 1;
		if (tokens->ntok > MAXTOKENS) return (MAXTOKENS);
		tokens->tok1[tokens->ntok] = str0;
		tokens->ltok[tokens->ntok] = stri - str0;
		stri++;
		str0 = stri;
		}
	    }

	/* Keep going if not whitespace */
	else
	    stri++;
	}

    /* Add token terminated by end of line */
    if (str0 < stri) {
	tokens->ntok = tokens->ntok + 1;
	if (tokens->ntok > MAXTOKENS)
	    return (MAXTOKENS);
	tokens->tok1[tokens->ntok] = str0;
	tokens->ltok[tokens->ntok] = stri - str0;
	}

    tokens->itok = 0;

    return (tokens->ntok);
}


/* NEXTOKEN -- get next token from tokenized string */

int
nextoken (tokens, token)
 
struct Tokens *tokens;	/* Token structure returned */
char	*token;		/* token (returned) */
{
    int it, ltok;

    tokens->itok = tokens->itok + 1;
    it = tokens->itok;
    if (it > tokens->ntok)
	it = tokens->ntok;
    else if (it < 1)
	it = 1;
    ltok = tokens->ltok[it];
    strncpy (token, tokens->tok1[it], ltok);
    token[ltok] = (char) 0;
    return (ltok);
}


/* GETOKEN -- get specified token from tokenized string */

int
getoken (tokens, itok, token)

struct Tokens *tokens;	/* Token structure returned */
int itok;		/* token sequence number of token
			 * if <0, get whole string after token -itok
			 * if =0, get whole string */
char *token;		/* token (returned) */
{
    int ltok;		/* length of token string (returned) */
    int it;

    it = itok;
    if (it > 0 ) {
	if (it > tokens->ntok)
	    it = tokens->ntok;
	ltok = tokens->ltok[it];
	strncpy (token, tokens->tok1[it], ltok);
	}
    else if (it < 0) {
	if (it < -tokens->ntok)
	    it  = -tokens->ntok;
	ltok = tokens->line + tokens->lline - tokens->tok1[-it];
	strncpy (token, tokens->tok1[-it], ltok);
	}
    else {
	ltok = tokens->lline;
	strncpy (token, tokens->tok1[1], ltok);
	}
    token[ltok] = (char) 0;

    return (ltok);
}


/* Read the right ascension, ra, in sexagesimal hours from in[] */

static double
cat2ra (in)

char	*in;	/* Character string */

{
    double ra;	/* Right ascension in degrees (returned) */

    ra = cat2dec (in);
    ra = ra * 15.0;

    return (ra);
}


/* Read the declination, dec, in sexagesimal degrees from in[] */

static double
cat2dec (in)

char	*in;	/* Character string */

{
    double dec;		/* Declination in degrees (returned) */
    double deg, min, sec, sign;
    char *value, *c1;

    dec = 0.0;

    /* Translate value from ASCII colon-delimited string to binary */
    if (!in[0])
	return (dec);
    else
	value = in;

    /* Set sign */
    if (!strchr (value,'-'))
	sign = 1.0;
    else {
	sign = -1.0;
	value = strchr (value,'-') + 1;
	}

    /* Translate value from ASCII colon-delimited string to binary */
    if ((c1 = strchr (value,':')) != NULL) {
	*c1 = 0;
	deg = (double) atoi (value);
	*c1 = ':';
	value = c1 + 1;
	if ((c1 = strchr (value,':')) != NULL) {
	    *c1 = 0;
	    min = (double) atoi (value);
	    *c1 = ':';
	    value = c1 + 1;
	    sec = atof (value);
	    }
	else {
	    sec = 0.0;
	    if ((c1 = strchr (value,'.')) != NULL)
		min = atof (value);
	    if (strlen (value) > 0)
		min = (double) atoi (value);
	    }
	dec = sign * (deg + (min / 60.0) + (sec / 3600.0));
	}
    
    /* Translate value from dd.mmssss */
    else if ((c1 = strchr (value,'.')) != NULL) {
	double xnum, deg, min, sec;
	xnum = atof (value);
	deg = (double)((int) (xnum + 0.000000001));
	xnum = (xnum - deg) * 100.0;
	min = (double)((int) (xnum + 0.000000001));
	sec = (xnum - deg) * 100.0;
	dec = sign * (deg + (min / 60.0) + (sec / 3600.0));
	}

    /* Translate integer */
    else 
	dec = sign * (double) atoi (value);

    return (dec);
}

/* Oct 16 1998	New subroutines
 * Oct 20 1998	Clean up error messages
 * Oct 26 1998	Return object names in catread() and catrnum()
 * Oct 29 1998	Correctly assign numbers when too many stars are found
 * Oct 30 1998	Fix epoch and equinox for J2000
 * Nov  9 1998	Drop out of star loop if rasorted catalog and past max RA
 * Dec  8 1998	Do not declare catsize() static
 * Dec 21 1998	Fix parsing so first character of line is not dropped

 * Jan 20 1999	Use strchr() instead of strsrch() for single char searches
 * Jan 29 1999	Default to star id number present
 * Feb  1 1999	Rewrite tokenizing subroutines for clarity
 * Feb  1 1999	Add match argument to catrnum()
 * Feb  2 1999	Add code to count decimal places in numbers
 * Feb  2 1999	Set sysout, eqout, and epout in catrnum() if not set
 * Feb 10 1999	Implement per star coordinate system
 * Feb 11 1999	Change starcat.insys to starcat.coorsys for consistency
 * Feb 17 1999	Fix per star coordinate system bugs
 * May 20 1999	Add option to read epoch of coordinates
 * Jun 16 1999	Use SearchLim()
 * Aug 16 1999	Fix bug to fix failure to search across 0:00 RA
 * Aug 25 1999	Return real number of stars from catread()
 */
