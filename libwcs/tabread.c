/*** File libwcs/tabread.c
 *** September 27, 2000
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

/* int tabread()	Read tab table stars in specified region
 * int tabrnum()	Read tab table stars with specified numbers
 * int tabxyread()	Read x, y, and magnitude from tab table star list
 * int tabrkey()	Read single keyword from specified tab table stars
 * struct StarCat tabcatopen()	Open tab table catalog, return number of entries
 * struct TabTable *tabopen()	Open tab table, returning number of entries
 * char *tabline()	Get tab table entry for one line
 * double tabgetra()	Return double right ascension in degrees
 * double tabgetdec()	Return double declination in degrees
 * double tabgetr8()	Return 8-byte floating point number from tab table line
 * int tabgeti4()	Return 4-byte integer from tab table line
 * int tabgetk()	Return character entry from tab table line for column
 * int tabgetc()	Return n'th character entry from tab table line
 * int tabhgetr8()	Return 8-byte floating point keyword value from header
 * int tabhgeti4()	Return 4-byte integer keyword value from header
 * int tabhgetc()	Return character keyword value from header
 * int tabparse()	Make a table of column headings
 * int tabcol()		Search a table of column headings for a particlar entry
 * int tabsize()	Return length of file in bytes
 * int istab()		Return 1 if first line of file contains a tab, else 0
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include "fitshead.h"
#include "wcs.h"
#include "wcscat.h"

#define ABS(a) ((a) < 0 ? (-(a)) : (a))

static int tabparse();
static int tabhgetr8();
static int tabhgeti4();
static int tabhgetc();
static int tabcont();
static int tabsize();
static double lastnum = 0.0;
static int nndec0 = 0;
static int nndec = 0;
static char *taberr;

char *gettaberr ()
{ return (taberr); }

int gettabndec()
{ return (nndec); }

static char *kwo = NULL;	/* Keyword returned by tabread(), tabrnum() */
void settabkey (keyword0)
char *keyword0;
{ kwo = keyword0; return; }


/* TABREAD -- Read tab table stars in specified region */

int
tabread (tabcatname,distsort,cra,cdec,dra,ddec,drad,
	 sysout,eqout,epout,mag1,mag2,nstarmax,starcat,
	 tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tpeak,tkey,nlog)

char	*tabcatname;	/* Name of reference star catalog file */
double	cra;		/* Search center J2000 right ascension in degrees */
double	cdec;		/* Search center J2000 declination in degrees */
double	dra;		/* Search half width in right ascension in degrees */
double	ddec;		/* Search half-width in declination in degrees */
double	drad;		/* Limiting separation in degrees (ignore if 0) */
int	distsort;	/* 1 to sort stars by distance from center */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
double	mag1,mag2;	/* Limiting magnitudes (none if equal) */
int	nstarmax;	/* Maximum number of stars to be returned */
struct StarCat **starcat; /* Star catalog data structure */
double	*tnum;		/* Array of UJ numbers (returned) */
double	*tra;		/* Array of right ascensions (returned) */
double	*tdec;		/* Array of declinations (returned) */
double	*tpra;		/* Array of right ascension proper motions (returned) */
double	*tpdec;		/* Array of declination proper motions (returned) */
double	*tmag;		/* Array of magnitudes (returned) */
double	*tmagb;		/* Array of second magnitudes (returned) */
int	*tpeak;		/* Array of peak counts (returned) */
char	**tkey;		/* Array of additional keyword values */
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
    struct Star *star;
    struct StarCat *sc;	/* Star catalog data structure */

    int wrap;
    int jstar;
    int nstar;
    char *objname;
    int lname;
    double ra,dec, rapm, decpm;
    double mag, magb;
    double num;
    int peak, i;
    int istar, nstars;
    int verbose;

    sc = *starcat;

    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    SearchLim (cra,cdec,dra,ddec,sysout,&ra1,&ra2,&dec1,&dec2,verbose);

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
    tdist = (double *) calloc (nstarmax, sizeof (double));

    star = (struct Star *) calloc (1, sizeof (struct Star));
    if (sc == NULL)
	sc = tabcatopen (tabcatname);
    *starcat = sc;
    if (sc == NULL || sc->nstars <= 0) {
	if (taberr != NULL)
	    fprintf (stderr,"%s\n", taberr);
	fprintf (stderr,"TABREAD: Cannot read catalog %s\n", tabcatname);
	free (star);
	sc = NULL;
	return (0);
	}

    nstars = sc->nstars;
    jstar = 0;

    /* Set catalog coordinate system */
    if (sc->equinox != 0.0)
	eqref = sc->equinox;
    else
	eqref = eqout;
    if (sc->epoch != 0.0)
	epref = sc->epoch;
    else
	epref = epout;
    if (sc->coorsys)
	sysref = sc->coorsys;
    else
	sysref = sysout;
    wcscstr (cstr, sysout, eqout, epout);

    /* Loop through catalog */
    for (istar = 1; istar <= nstars; istar++) {

	/* Read position of next star */
	if (tabstar (istar, sc, star)) {
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
	if (sc->mprop)
	    wcsconp (sysref, sysout, eqref, eqout, epref, epout,
		     &ra, &dec, &rapm, &decpm);
	else
	    wcscon (sysref, sysout, eqref, eqout, &ra, &dec, epout);
	mag = star->xmag[0];
	magb = star->xmag[1];
	peak = star->peak;
	if (drad > 0 || distsort)
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
		if (sc->mprop) {
		    tpra[nstar] = rapm;
		    tpdec[nstar] = decpm;
		    }
		tmag[nstar] = mag;
		tmagb[nstar] = magb;
		tpeak[nstar] = peak;
		tdist[nstar] = dist;
		if (kwo != NULL) {
		    lname = strlen (star->objname) + 1;
		    objname = (char *)calloc (lname, 1);
		    strcpy (objname, star->objname);
		    tkey[nstar] = objname;
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

	    /* If radial search & too many stars, replace furthest star */
	    else if (distsort) {
		if (dist < maxdist) {
		    tnum[farstar] = num;
		    tra[farstar] = ra;
		    tdec[farstar] = dec;
		    if (sc->mprop) {
			tpra[farstar] = rapm;
			tpdec[farstar] = decpm;
			}
		    tmag[farstar] = mag;
		    tmagb[farstar] = magb;
		    tpeak[farstar] = peak;
		    tdist[farstar] = dist;
		    if (kwo != NULL) {
			lname = strlen (star->objname) + 1;
			objname = (char *)calloc (lname, 1);
			strcpy (objname, star->objname);
			tkey[farstar] = objname;
			}

		    /* Find new farthest star */
		    maxdist = 0.0;
		    for (i = 0; i < nstarmax; i++) {
			if (tdist[i] > maxdist) {
			    maxdist = tdist[i];
			    farstar = i;
			    }
			}
		    }
		}

	    /* Otherwise if too many stars, replace faintest star */
	    else if (mag < faintmag) {
		tnum[faintstar] = num;
		tra[faintstar] = ra;
		tdec[faintstar] = dec;
		if (sc->mprop) {
		    tpra[faintstar] = rapm;
		    tpdec[faintstar] = decpm;
		    }
		tmag[faintstar] = mag;
		tmagb[faintstar] = magb;
		tpeak[faintstar] = peak;
		tdist[faintstar] = dist;
		if (kwo != NULL) {
		    lname = strlen (star->objname) + 1;
		    objname = (char *)calloc (lname, 1);
		    strcpy (objname, star->objname);
		    tkey[faintstar] = objname;
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
		fprintf (stderr,"TABREAD: %11.6f: %9.5f %9.5f %s %5.2f %5.2f %d    \n",
			 num,ra,dec,cstr,mag,magb,peak);

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

    free ((char *)tdist);
    return (nstar);
}


/* TABRNUM -- Read tab table stars with specified numbers */

int
tabrnum (tabcatname, nnum, sysout, eqout, epout,
	 tnum, tra, tdec, tpra, tpdec, tmag, tmagb, tpeak, tkey, nlog)

char	*tabcatname;	/* Name of reference star catalog file */
int	nnum;		/* Number of stars to look for */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
double	*tnum;		/* Array of star numbers to look for */
double	*tra;		/* Array of right ascensions (returned) */
double	*tdec;		/* Array of declinations (returned) */
double	*tpra;		/* Array of right ascension proper motions (returned) */
double	*tpdec;		/* Array of declination proper motions (returned) */
double	*tmag;		/* Array of magnitudes (returned) */
double	*tmagb;		/* Array of second magnitudes (returned) */
int	*tpeak;		/* Array of peak counts (returned) */
char	**tkey;		/* Array of additional keyword values */
int	nlog;
{
    int jnum;
    int nstar;
    double ra,dec, rapm, decpm;
    double mag, magb;
    double num;
    int peak;
    int istar, nstars;
    char *line;
    int sysref;		/* Catalog coordinate system */
    double eqref;	/* Catalog equinox */
    double epref;	/* Catalog epoch */
    char cstr[32];
    char *objname;
    int lname;
    struct TabTable *startab;
    struct StarCat *starcat;
    struct Star *star;

    line = 0;

    nstar = 0;
    nndec = 0;
    nndec0 = 0;

    /* Allocate catalog entry buffer */
    star = (struct Star *) calloc (1, sizeof (struct Star));

    /* Open star catalog */
    starcat = tabcatopen (tabcatname);
    if (starcat == NULL || starcat->nstars <= 0) {
	if (taberr != NULL)
	    fprintf (stderr,"%s\n", taberr);
	fprintf (stderr,"TABRNUM: Cannot read catalog %s\n", tabcatname);
	free (star);
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
	    if (tabstar (istar, starcat, star))
		fprintf (stderr,"TABRNUM: Cannot read star %d\n", istar);

	    /* If star entry has been read successfully */
	    else {

		/* Set coordinate system for this star */
		sysref = star->coorsys;
		eqref = star->equinox;

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
		magb = star->xmag[1];
		peak = star->peak;

		/* Save star position and magnitude in table */
		tnum[jnum] = num;
		tra[jnum] = ra;
		tdec[jnum] = dec;
		if (starcat->mprop) {
		    tpra[jnum] = rapm;
		    tpdec[jnum] = decpm;
		    }
		tmag[jnum] = mag;
		tmagb[jnum] = magb;
		tpeak[jnum] = peak;
		if (kwo != NULL) {
		    lname = strlen (star->objname) + 1;
		    objname = (char *)calloc (lname, 1);
		    strcpy (objname, star->objname);
		    tkey[nstar] = objname;
		    }
		nstar++;
		if (nlog == 1)
		    fprintf (stderr,"TABRNUM: %11.6f: %9.5f %9.5f %s %5.2f %5.2f %d    \n",
			     num,ra,dec,cstr,mag,magb,peak);
		/* End of accepted star processing */
		}
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


/* TABXYREAD -- Read X, Y, and magnitude of tab table stars */

int
tabxyread (tabcatname, xa, ya, ba, pa, nlog)

char	*tabcatname;	/* Name of reference star catalog file */
double	**xa;		/* Array of x coordinates (returned) */
double	**ya;		/* Array of y coordinates (returned) */
double	**ba;		/* Array of fluxes (returned) */
int	**pa;		/* Array of magnitudes*100 (returned) */
int	nlog;
{
    double xi, yi, magi, flux;
    char *line;
    int istar, nstars;
    int verbose;
    struct TabTable *startab;
    int entx, enty, entmag;

    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    /* Open tab table file */
    nndec = 0;
    nndec0 = 0;
    startab = tabopen (tabcatname);
    if (startab == NULL || startab->nlines <= 0) {
	fprintf (stderr,"TABXYREAD: Cannot read catalog %s\n", tabcatname);
	return (0);
	}

    /* Find columns for X, Y, and magnitude */
    if (!(entx = tabcol (startab, "X")))
        entx = tabcol (startab, "x");
    if (!(enty = tabcol (startab, "Y")))
        enty = tabcol (startab, "y");
    if (!(entmag = tabcol (startab, "MAG")))
        entmag = tabcol (startab, "mag");

    /* Allocate vectors for x, y, magnitude, and flux */
    nstars = startab->nlines;
    *xa = (double *) realloc(*xa, nstars*sizeof(double));
    if (*xa == NULL) {
	fprintf (stderr,"TABXYREAD: Cannot allocate memory for x\n");
	return (0);
	}
    *ya = (double *) realloc(*ya, nstars*sizeof(double));
    if (*ya == NULL) {
	fprintf (stderr,"TABXYREAD: Cannot allocate memory for y\n");
	return (0);
	}
    *ba = (double *) realloc(*ba, nstars*sizeof(double));
    if (*ba == NULL) {
	fprintf (stderr,"TABXYREAD: Cannot allocate memory for mag\n");
	return (0);
	}
    *pa = (int *) realloc(*pa, nstars*sizeof(int));
    if (*pa == NULL) {
	fprintf (stderr,"TABXYREAD: Cannot allocate memory for flux\n");
	return (0);
	}

    /* Loop through catalog */
    for (istar = 0; istar < nstars; istar++) {

	/* Read line for next star */
	if ((line = tabline (startab, istar)) == NULL) {
	    fprintf (stderr,"TABXYREAD: Cannot read star %d\n", istar);
	    break;
	    }

	/* Extract x, y, and magnitude */
	xi = tabgetr8 (startab, line, entx);
	yi = tabgetr8 (startab, line, enty);
	magi = tabgetr8 (startab, line, entmag);

	(*xa)[istar] = xi;
	(*ya)[istar] = yi;
	flux = 10000.0 * pow (10.0, (-magi / 2.5));
	(*ba)[istar] = flux;
	(*pa)[istar] = (int)(magi * 100.0);

	if (nlog == 1)
	    fprintf (stderr,"DAOREAD: %6d/%6d: %9.5f %9.5f %15.2f %6.2f\n",
		     istar,nstars,xi,yi,flux,magi);

	/* Log operation */
	if (nlog > 1 && istar%nlog == 0)
		fprintf (stderr,"TABXYREAD: %5d / %5d sources catalog %s\r",
			istar,nstars,tabcatname);

	/* End of star loop */
	}

    /* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"TABXYREAD: Catalog %s : %d / %d found\n",tabcatname,
		 istar,nstars);

    /* Free table */
    tabclose (startab);
    if (istar < nstars-1)
	return (istar + 1);
    else
	return (nstars);
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
    if (starcat == NULL) {
	if (taberr != NULL)
	    fprintf (stderr,"%s\n", taberr);
	fprintf (stderr,"%s\n", taberr);
	return (0);
	}
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
tabcatopen (tabpath)

char *tabpath;		/* Tab table catalog file pathname */
{
    char cstr[32];
    char *tabname;
    struct TabTable *startab;
    struct StarCat *sc;
    int i, lnum, ndec, istar;
    char *line;
    double dnum;

    /* Allocate catalog data structure */
    sc = (struct StarCat *) calloc (1, sizeof (struct StarCat));

    /* Open the tab table file */
    if ((startab = tabopen (tabpath)) == NULL)
	return (NULL);
    sc->startab = startab;

    /* Save name of catalog */
    tabname = strrchr (tabpath, '/');
    if (tabname)
	tabname = tabname + 1;
    else
	tabname = tabpath;
    if (strlen (tabname) < 24)
	strcpy (sc->isfil, tabname);
    else {
	strncpy (sc->isfil, tabname, 23);
	sc->isfil[23] = (char) 0;
	}

    /* Find column and name of object identifier */
    sc->entid = -1;
    sc->keyid[0] = (char) 0;
    if ((sc->entid = tabcol (startab, "ID")))
	strcpy (sc->keyid, "ID");
    else if ((sc->entid = tabcol (startab, "id")))
	strcpy (sc->keyid, "id");
    else if ((sc->entid = tabcont (startab, "_id"))) {
	i = sc->entid - 1;
	strncpy (sc->keyid, startab->colname[i], startab->lcol[i]);
	}
    else if ((sc->entid = tabcont (startab, "ID"))) {
	i = sc->entid - 1;
	strncpy (sc->keyid, startab->colname[i], startab->lcol[i]);
	}
    else if ((sc->entid = tabcont (startab, "name"))) {
	i = sc->entid - 1;
	strncpy (sc->keyid, startab->colname[i], startab->lcol[i]);
	}
    sc->nndec = nndec;

    /* Find column and name of object right ascension */
    sc->entra = -1;
    sc->keyra[0] = (char) 0;
    if ((sc->entra = tabcol (startab, "RA")))
	strcpy (sc->keyra, "RA");
    else if ((sc->entra = tabcol (startab, "ra")))
	strcpy (sc->keyra, "ra");
    else if ((sc->entra = tabcol (startab, "Ra")))
	strcpy (sc->keyra, "ra");
    else if ((sc->entra = tabcont (startab, "ra"))) {
	i = sc->entra - 1;
	strncpy (sc->keyra, startab->colname[i], startab->lcol[i]);
	}

    /* Find column and name of object declination */
    sc->entdec = -1;
    sc->keydec[0] = (char) 0;
    if ((sc->entdec = tabcol (startab, "DEC")))
	strcpy (sc->keydec, "DEC");
    else if ((sc->entdec = tabcol (startab, "dec")))
	strcpy (sc->keydec, "dec");
    else if ((sc->entdec = tabcol (startab, "Dec")))
	strcpy (sc->keydec, "dec");
    else if ((sc->entdec = tabcont (startab, "dec"))) {
	i = sc->entdec;
	strncpy (sc->keydec, startab->colname[i], startab->lcol[i]);
	}

    /* Find column and name of object first magnitude */
    sc->entmag1 = -1;
    sc->keymag1[0] = (char) 0;
    if ((sc->entmag1 = tabcol (startab, "MAG")))
	strcpy (sc->keymag1, "MAG");
    else if ((sc->entmag1 = tabcol (startab, "mag")))
	strcpy (sc->keymag1, "MAG");
    else if ((sc->entmag1 = tabcol (startab, "magr")))
	strcpy (sc->keymag1, "magr");
    else if ((sc->entmag1 = tabcont (startab, "mag"))) {
	i = sc->entmag1 - 1;
	strncpy (sc->keymag1, startab->colname[i], startab->lcol[i]);
	}

    /* Find column and name of object second magnitude */
    sc->entmag2 = -1;
    sc->keymag2[0] = (char) 0;
    if ((sc->entmag2 = tabcol (startab, "magb")))
	strcpy (sc->keymag2, "magb");

    /* Find column and name of object right ascension proper motion */
    sc->entrpm = -1;
    sc->keyrpm[0] = (char) 0;
    if ((sc->entrpm = tabcol (startab, "URA")))
	strcpy (sc->keyrpm, "URA");
    else if ((sc->entrpm = tabcol (startab, "ura")))
	strcpy (sc->keyrpm, "ura");
    else if ((sc->entrpm = tabcol (startab, "Ura")))
	strcpy (sc->keyrpm, "Ura");
    else if ((sc->entrpm = tabcol (startab, "Ux")))
	strcpy (sc->keyrpm, "Ux");

    /* Find column and name of object declination proper motion */
    sc->entdpm = -1;
    sc->keydpm[0] = (char) 0;
    if ((sc->entdpm = tabcol (startab, "UDEC")))
	strcpy (sc->keydpm, "UDEC");
    else if ((sc->entdpm = tabcol (startab, "udec")))
	strcpy (sc->keydpm, "udec");
    else if ((sc->entdpm = tabcol (startab, "Udec")))
	strcpy (sc->keydpm, "Udec");
    else if ((sc->entdpm = tabcol (startab, "Uy")))
	strcpy (sc->keydpm, "Uy");

    /* Find units for RA proper motion */
    sc->mprop = 0;
    cstr[0] = 0;
    if (!tabhgetc (startab,"RPMUNIT", cstr)) {
	if (!tabhgetc (startab,"rpmunit", cstr)) {
	    if (!tabhgetc (startab,"pmunit", cstr))
		tabhgetc (startab,"pmunit", cstr);
	    }
	}
    if (strlen (cstr) > 0) {
	sc->mprop = 1;
	if (!strcmp (cstr, "mas/yr") || !strcmp (cstr, "mas/year"))
	    sc->rpmunit = PM_MASYR;
	else if (!strcmp (cstr, "arcsec/yr") || !strcmp (cstr, "arcsec/year"))
	    sc->rpmunit = PM_ARCSECYR;
	else if (!strcmp (cstr, "rad/yr") || !strcmp (cstr, "rad/year"))
	    sc->rpmunit = PM_RADYR;
	else if (!strcmp (cstr, "sec/yr") || !strcmp (cstr, "sec/year"))
	    sc->rpmunit = PM_TSECYR;
	else if (!strcmp (cstr, "tsec/yr") || !strcmp (cstr, "tsec/year"))
	    sc->rpmunit = PM_TSECYR;
	else
	    sc->rpmunit = PM_DEGYR;
	}

    /* Find units for Dec proper motion */
    cstr[0] = 0;
    if (!tabhgetc (startab,"DPMUNIT", cstr)) {
	if (!tabhgetc (startab,"dpmunit", cstr)) {
	    if (!tabhgetc (startab,"pmunit", cstr))
		tabhgetc (startab,"pmunit", cstr);
	    }
	}
    if (strlen (cstr) > 0) {
	sc->mprop = 1;
	if (!strcmp (cstr, "mas/yr") || !strcmp (cstr, "mas/year"))
	    sc->dpmunit = PM_MASYR;
	else if (!strcmp (cstr, "sec/yr") || !strcmp (cstr, "sec/year"))
	    sc->dpmunit = PM_ARCSECYR;
	else if (!strcmp (cstr, "tsec/yr") || !strcmp (cstr, "tsec/year"))
	    sc->dpmunit = PM_ARCSECYR;
	else if (!strcmp (cstr, "arcsec/yr") || !strcmp (cstr, "arcsec/year"))
	    sc->dpmunit = PM_ARCSECYR;
	else if (!strcmp (cstr, "rad/yr") || !strcmp (cstr, "rad/year"))
	    sc->dpmunit = PM_RADYR;
	else
	    sc->dpmunit = PM_DEGYR;
	}

    /* Find column and name of object peak or plate number */
    sc->entpeak = -1;
    sc->keypeak[0] = (char) 0;
    if ((sc->entpeak = tabcol (startab, "PEAK")))
	strcpy (sc->keypeak, "PEAK");
    else if ((sc->entpeak = tabcol (startab, "peak")))
	strcpy (sc->keypeak, "peak");
    else if ((sc->entpeak = tabcol (startab, "plate")))
	strcpy (sc->keypeak, "plate");

    sc->entadd = -1;
    sc->keyadd[0] = (char) 0;
    if (kwo != NULL)
	sc->entadd = tabcol (startab, kwo);

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
    if (tabhgetc (startab,"radecsys", cstr)) {
	sc->coorsys = wcscsys (cstr);
	if (!tabhgetr8 (startab,"equinox", &sc->equinox))
	    sc->equinox = wcsceq (cstr);
	if (!tabhgetr8 (startab,"epoch",&sc->epoch))
	    sc->epoch = sc->equinox;
	}
    else if (tabhgetr8 (startab,"EQUINOX", &sc->equinox)) {
	if (!tabhgetr8 (startab,"EPOCH",&sc->epoch))
	    sc->epoch = sc->equinox;
	if (sc->equinox == 1950.0)
	    sc->coorsys = WCS_B1950;
	else
	    sc->coorsys = WCS_J2000;
	}
    else if (tabhgetr8 (startab,"equinox", &sc->equinox)) {
	if (!tabhgetr8 (startab,"epoch",&sc->epoch))
	    sc->epoch = sc->equinox;
	if (sc->equinox == 1950.0)
	    sc->coorsys = WCS_B1950;
	else
	    sc->coorsys = WCS_J2000;
	}
    else if (tabhgetr8 (startab,"EPOCH", &sc->epoch)) {
	sc->equinox = sc->epoch;
	if (sc->equinox == 1950.0)
	    sc->coorsys = WCS_B1950;
	else
	    sc->coorsys = WCS_J2000;
	}
    else if (tabhgetr8 (startab,"epoch", &sc->epoch)) {
	sc->equinox = sc->epoch;
	if (sc->equinox == 1950.0)
	    sc->coorsys = WCS_B1950;
	else
	    sc->coorsys = WCS_J2000;
	}

    /* Set other stuff */
    sc->rasorted = 0;
    sc->nstars = startab->nlines;
    if (sc->entmag1 && sc->entmag2)
	sc->nmag = 2;
    else if (sc->entmag1)
	sc->nmag = 1;
    else
	sc->nmag = 0;
    if (sc->entid)
	sc->stnum = 1;
    else
	sc->stnum = 0;

    /* Find out if ID is number, and if so, how many decimal places it has */
    if (sc->entid) {

	istar = 1;
	if ((line = tabline (startab, istar)) == NULL) {
	    fprintf (stderr,"TABCATOPEN: Cannot read first star\n");
	    tabcatclose (sc);
	    return (NULL);
	    }
	tabgetc (startab, line, sc->entid, cstr, 32);

	/* Find number of decimal places in identifier */
	if (tabhgeti4 (startab, "ndec", &nndec)) {
	    sc->nndec = nndec;
	    if (!isnum (cstr))
		sc->stnum = 5;
	    }
	else if (isnum (cstr)) {
	    dnum = tabgetr8 (startab,line,sc->entid);
	    sprintf (cstr,"%.0f", (dnum * 100000000.0) + 0.1);
	    lnum = strlen (cstr);
	    for (i = 0; i < 8; i++) {
		if (cstr[lnum-i-1] != '0') {
		    ndec = 8 - i;
		    if (ndec > nndec) {
			nndec = ndec;
			sc->nndec = nndec;
			}
		    break;
		    }
		}
	    sc->nndec = nndec;
	    }
	else {
	    sc->stnum = 5;
	    sc->nndec = nndec;
	    }
	}
    sc->refcat = TABCAT;

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
   return 0 if successful, else -1 */

int
tabstar (istar, sc, st)

int istar;      /* Star sequence number in tab table catalog */
struct StarCat *sc; /* Star catalog data structure */
struct Star *st; /* Star data structure, updated on return */
{
    struct TabTable *startab = sc->startab;
    char *line;
    char cnum[32];
    int ndec, i;
    int lnum;

    if ((line = tabline (startab, istar)) == NULL) {
	fprintf (stderr,"TABSTAR: Cannot read star %d\n", istar);
	return (-1);
	}

    /* Extract selected fields  */
    if (sc->entid) {
	tabgetc (startab, line, sc->entid, cnum, 32);
	if (isnum (cnum) || isnum (cnum+1)) {
	    if (isnum(cnum))
		st->num = atof (cnum);
	    else
		st->num = atof (cnum+1);
	    lastnum = st->num;

	    /* Find number of decimal places in identifier */
	    if (strchr (cnum,'.') == NULL) {
		nndec = 0;
		sc->nndec = nndec;
		}
	    else {
		sprintf (cnum,"%.0f", (st->num * 100000000.0) + 0.1);
		lnum = strlen (cnum);
		for (i = 0; i < 8; i++) {
		    if (cnum[lnum-i-1] != '0') {
			ndec = 8 - i;
			if (ndec > nndec) {
			    nndec = ndec;
			    sc->nndec = nndec;
			    }
			break;
			}
		    }
		}
	    nndec0 = nndec;
	    }
	else {
	    strcpy (st->objname, cnum);
	    st->num = st->num + 1.0;
	    }
	}
    else {
	st->num = (double) istar;
	nndec = 0;
	sc->nndec = nndec;
	}

    /* Right ascension */
    st->ra = tabgetra (startab, line, sc->entra);

    /* Declination */
    st->dec = tabgetdec (startab, line, sc->entdec);

    /* Magnitudes */
    st->xmag[0] = tabgetr8 (startab, line, sc->entmag1);
    if (sc->entmag2)
	st->xmag[1] = tabgetr8 (startab, line, sc->entmag2);
    else
	st->xmag[1] = 0.0;

    /* Right ascension proper motion */
    st->rapm = tabgetpm (startab, line, sc->entrpm);
    if (sc->rpmunit == PM_MASYR)
	st->rapm = (st->rapm / 3600000.0) / cosdeg (st->dec);
    else if (sc->rpmunit == PM_ARCSECYR)
	st->rapm = (st->rapm / 3600.0) / cosdeg (st->dec);
    else if (sc->rpmunit == PM_TSECYR)
	st->rapm = 15.0 * st->rapm / 3600.0;
    else if (sc->rpmunit == PM_RADYR)
	st->rapm = raddeg (st->rapm);

    /* Declination proper motion */
    st->decpm = tabgetpm (startab, line, sc->entdpm);
    if (sc->dpmunit == PM_MASYR)
	st->decpm = st->decpm / 3600000.0;
    else if (sc->dpmunit == PM_ARCSECYR)
	st->decpm = st->decpm / 3600.0;
    else if (sc->dpmunit == PM_RADYR)
	st->decpm = raddeg (st->decpm);

    /* Peak counts */
    st->peak = tabgeti4 (startab, line, sc->entpeak);

    /* Extract selected field */
    if (kwo != NULL)
	(void) tabgetk (startab, line, kwo, st->objname, 32);
    else
	st->objname[0] = (char) 0;

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
    int nr, lfile, lfname, lfa, lname;
    char *tabnew, *tabline, *lastline;
    char *tabcomma;
    char *thisname, *tabname;
    int thistab, itab, nchar;
    int formfeed = (char) 12;
    struct TabTable *tabtable;

    tabcomma = NULL;
    if (taberr != NULL) {
	free (taberr);
	taberr = NULL;
	}

    tabname = NULL;
    if (!strcmp (tabfile, "stdin")) {
	lfile = 100000;
	fcat = stdin;
	}
    else {

	/* Separate table name from file name, if necessary */
	if ((tabcomma = strchr (tabfile, ',')) != NULL) {
	    tabname = (char *) calloc (1,64);
	    strcpy (tabname, tabcomma+1);
	    *tabcomma = (char) 0;
	    lname = strlen (tabname);
	    }

	/* Find length of tab table catalog */
	lfile = tabsize (tabfile);
	if (lfile < 1) {
	    taberr = (char *) calloc (64 + strlen (tabfile), 1);
	    sprintf (taberr,"TABOPEN: Tab table file %s has no entries",
		     tabfile);
	    if (tabcomma != NULL) *tabcomma = ',';
	    return (NULL);
	    }

	/* Open tab table catalog */
	if (!(fcat = fopen (tabfile, "r"))) {
	    taberr = (char *) calloc (64 + strlen (tabfile), 1);
	    sprintf (taberr,"TABOPEN: Tab table file %s cannot be read",
		     tabfile);
	    if (tabcomma != NULL) *tabcomma = ',';
	    return (NULL);
	    }
	}

    /* Allocate tab table structure */
    if ((tabtable=(struct TabTable *) calloc(1,sizeof(struct TabTable)))==NULL){
	taberr = (char *) calloc (64 + strlen (tabfile), 1);
	sprintf (taberr,"TABOPEN: cannot allocate Tab Table structure for %s",
		 tabfile);
	if (tabcomma != NULL) *tabcomma = ',';
	return (NULL);
	}

    tabtable->tabname = tabname;

    /* Allocate space in structure for filename and save it */
    lfname = strlen (tabfile) + 1;
    lfa = lfname / 64;
    if (lfa*64 < lfname)
	lfa = 64 * (lfa + 1);
    else
	lfa = 64 * lfa;
    if ((tabtable->filename = (char *)calloc (1, lfa)) == NULL) {
	taberr = (char *) calloc (64 + strlen (tabfile), 1);
	sprintf (taberr,"TABOPEN: cannot allocate filename %s in structure",
		 tabfile);
	(void) fclose (fcat);
	tabclose (tabtable);
	if (tabcomma != NULL) *tabcomma = ',';
	return (NULL);
	}
    strncpy (tabtable->filename, tabfile, lfname);

    /* Allocate buffer to hold entire catalog and read it */
    if ((tabtable->tabbuff = (char *) calloc (1, lfile+2)) == NULL) {
	taberr = (char *) calloc (64 + strlen (tabfile), 1);
	sprintf (taberr,"TABOPEN: cannot allocate buffer for tab table %s",
		 tabfile);
	(void) fclose (fcat);
	tabclose (tabtable);
	if (tabcomma != NULL) *tabcomma = ',';
	return (NULL);
	}
    else {
	nr = fread (tabtable->tabbuff, 1, lfile, fcat);
	if (fcat != stdin && nr < lfile) {
	    fprintf (stderr,"TABOPEN: read only %d / %d bytes of file %s\n",
		     nr, lfile, tabfile);
	    (void) fclose (fcat);
	    tabclose (tabtable);
	    if (tabcomma != NULL) *tabcomma = ',';
	    return (NULL);
	    }

	/* Check for named table within a file */
	if (tabname != NULL) {
	    if (isnum (tabname)) {
		itab = atoi (tabname);
		thisname = tabtable->tabbuff;
		thistab = 1;
		if (itab > 1) {
		    while (thistab < itab && thisname != NULL) {
			thisname = strchr (thisname, formfeed);
			if (thisname != NULL)
			    thisname++;
			thistab++;
			}
		    }
		if (thisname == NULL) {
		    fprintf (stderr, "GETTAB:  There are < %d tables in %s\n",
			itab, tabfile);
		    return (NULL);
		    }
		while (*thisname==' ' || *thisname==newline ||
		       *thisname==formfeed || *thisname==(char)13)
		    thisname++;
		tabline = strchr (thisname, newline);
		if (tabline != NULL) {
		    nchar = tabline - thisname;
		    if (strchr (thisname, tab) > tabline)
			strncpy (tabtable->tabname, thisname, nchar);
		    }
		}
	    else {
		thisname = tabtable->tabbuff;
		while (*thisname != NULL) {
		    while (*thisname==' ' || *thisname==newline ||
		   	   *thisname==formfeed || *thisname==(char)13)
			thisname++;
		    if (!strncmp (tabname, thisname, lname))
			break;
		    else
			thisname = strchr (thisname, formfeed);
		    }
		}
	    if (thisname == NULL) {
		fprintf (stderr, "TABOPEN: table %s in file %s not found\n",
			 tabname, tabfile);
		if (tabcomma != NULL) *tabcomma = ',';
		return (NULL);
		}
	    else
		tabtable->tabheader = strchr (thisname, newline) + 1;
	    }
	else
	    tabtable->tabheader = tabtable->tabbuff;

	tabline = tabtable->tabheader;
	while (*tabline != '-' && tabline < tabtable->tabbuff+lfile) {
	    lastline = tabline;
	    tabline = strchr (tabline,newline) + 1;
	    }
	if (*tabline != '-') {
	    taberr = (char *) calloc (64 + strlen (tabfile), 1);
	    sprintf (taberr,"TABOPEN: No - line in tab table %s",tabfile);
	    (void) fclose (fcat);
	    tabclose (tabtable);
	    if (tabcomma != NULL) *tabcomma = ',';
	    return (NULL);
	    }
	tabtable->tabhead = lastline;
	tabtable->tabdata = strchr (tabline, newline) + 1;

	/* Extract positions of keywords we will want to use */
	if (!tabparse (tabtable)) {
	    fprintf (stderr,"TABOPEN: No columns in tab table %s\n",tabfile);
	    (void) fclose (fcat);
	    tabclose (tabtable);
	    if (tabcomma != NULL) *tabcomma = ',';
	    return (NULL);
	    }

    /* Enumerate entries in tab table catalog by counting newlines */
	tabnew = strchr (tabtable->tabdata, newline) + 1;
	tabnew = tabtable->tabdata;
	tabtable->nlines = 0;
	while ((tabnew = strchr (tabnew, newline)) != NULL) {
	    tabnew = tabnew + 1;
	    tabtable->nlines = tabtable->nlines + 1;
	    if (*tabnew == formfeed)
		break;
	    }
	}

    (void) fclose (fcat);
    tabtable->tabline = tabtable->tabdata;
    tabtable->iline = 1;
    if (tabcomma != NULL) *tabcomma = ',';
    return (tabtable);
}


void
tabclose (tabtable)

    struct TabTable *tabtable;
{
    if (tabtable != NULL) {
	if (tabtable->filename != NULL) free (tabtable->filename);
	if (tabtable->tabname != NULL) free (tabtable->tabname);
	if (tabtable->tabbuff != NULL) free (tabtable->tabbuff);
	if (tabtable->colname != NULL) free (tabtable->colname);
	if (tabtable->lcol != NULL) free (tabtable->lcol);
	if (tabtable->lcfld != NULL) free (tabtable->lcfld);
	free (tabtable);
	}
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

double
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

double
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

/* TABGETPM -- returns double RA or Dec proper motion in degrees/year */

double
tabgetpm (startab, line, ientry)

struct TabTable *startab;	/* Tab table structure */
char	*line;			/* tab table line */
int	ientry;			/* sequence of entry on line */
{
    char str[24];

    if (tabgetc (startab, line, ientry, str, 24))
	return (0.0);
    else
	return (atof (str));
}

/* TABGETR8 -- returns 8-byte floating point number from tab table line */

double
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

int
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

static int
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

static int
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

static int
tabhgetc (tabtable, keyword, result)

struct TabTable *tabtable;	/* Tab table structure */
char	*keyword;		/* sequence of entry on line */
char	*result;
{
    char *str0, *str1, *line, *head, keylow[24], keyup[24];
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
	if (line == tabtable->tabbuff || line[-1] == newline) {
	    str0 = strchr (line, tab) + 1;
	    str1 = strchr (line, newline);
	    break;
	    }
	else
	    head = line + 1;
	}

    /* Return value as a character string and 1 if found */
    if (str0) {
	ncstr = str1 - str0;
	strncpy (result, str0, ncstr);
	result[ncstr] = (char)0;
	return (1);
	}
    else
	return (0);
}


/* TABPARSE -- Make a table of column headings */

static int
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
    int nbytes, nba;

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
    nba = nbytes / 64;
    if (nbytes > nba * 64)
	nba = (nba + 1) * 64;
    else
	nba = nba * 64;
    tabtable->colname = (char **)calloc (nba, 1);
    nbytes = tabtable->ncols * sizeof (int);
    nba = nbytes / 64;
    if (nbytes > nba * 64)
	nba = (nba + 1) * 64;
    else
	nba = nba * 64;
    tabtable->lcol = (int *) calloc (nba, 1);
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
    nbytes = tabtable->ncols * sizeof (int);
    nba = nbytes / 64;
    if (nbytes > nba * 64)
	nba = (nba + 1) * 64;
    else
	nba = nba * 64;
    tabtable->lcfld = (int *) calloc (nba, 1);
    colhead = hyphens;
    i = 0;
    while (colhead) {
	nextab = strchr (colhead, tab);
	if (nextab < hyphlast)
	    endcol = nextab - 1;
	else
	    endcol = hyphlast - 1;
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


/* Search table of column headings for first with string (case-dependent) */

static int
tabcont (tabtable, keyword)

struct TabTable *tabtable;	/* Tab table structure */
char	*keyword;		/* Part of column heading to find */

{
    int i;

    for (i = 0; i < tabtable->ncols; i++) {
	if (strnsrch (tabtable->colname[i], keyword, tabtable->lcol[i])) {
	    return (i + 1);
	    }
	}
    return (0);
}


/* TABSIZE -- return size of file in bytes */

static int
tabsize (filename)

char	*filename;	/* Name of file for which to find size */
{
    FILE *diskfile;
    long filesize;

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
    struct TabTable *tabtable;

    /* First check file extension */
    if (strsrch (filename, ".tab"))
	return (1);

    /* If no .tab file extension, try opening the file */
    else {
	if ((tabtable = tabopen (filename)) != NULL) {
	    tabclose (tabtable);
	    return (1);
	    }
	else
	    return (0);
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
 * Sep 10 1999	Fix bug setting equinox and coordinate system in tabcatopen()
 * Sep 10 1999	Set additional keyword selection with subroutine
 * Sep 13 1999	Fix comment for tabstar()
 * Sep 16 1999	Fix bug which didn't always return closest stars
 * Sep 16 1999	Add distsort argument so brightest stars in circle works, too
 * Oct 21 1999	Clean up code after lint
 * Oct 25 1999	Fix subroutine declaration inconsistency
 * Oct 25 1999	Replace malloc() calls with calloc()
 * Oct 29 1999	Add tabxyread() for image catalogs
 * Nov 23 1999	Improve error checking on Starbase tables; istab() opens file
 * Nov 30 1999	Fix bugs found when compiling under SunOS 4.1.3
 *
 * Jan  4 2000	Always close file and print error message on tabopen() failures
 * Jan  6 2000	If "id" not found, try heading with "_id" to catch scat output
 * Jan 10 2000	Add second magnitude; save column headers in catalog structure
 * Feb 10 2000	Implement proper motion in source catalogs
 * Feb 10 2000	Accept first mag-containing column as first magnitude
 * Feb 10 2000	Clean up id reading: no. decimals, non-numeric id
 * Feb 14 2000	Save table opening errors in string
 * Feb 16 2000	Lengthen short calloc() lengths
 * Feb 16 2000	Pad tabbuff with 2 nulls so end can be found
 * Mar 10 2000	Return proper motions from tabread() and tabrnum()
 * Mar 13 2000	Do not free tabtable structure if it is null
 * Mar 27 2000	Clean up code after lint
 * May 26 2000	Add ability to read named tables in a multi-table file
 * Jun 26 2000	Add coordinate system to SearchLim() arguments
 * Jul  5 2000	Check for additional column heading variations
 * Jul 10 2000	Deal with number of decimal places and name/number in tabcatopen()
 * Jul 12 2000	Add star catalog data structure to tabread() argument list
 * Jul 13 2000	Use nndec header parameter to optionally set decimal places
 * Jul 25 2000	Pass star catalog address of data structure address
 * Aug  3 2000	Skip first character of ID if rest is number
 * Aug  3 2000	If no decimal point in numeric ID, set ndec to zero
 * Sep 27 2000	Use first column with name if no id column is found
 */
