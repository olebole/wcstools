/*** File libwcs/ctgread.c
 *** May 20, 2000
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

/* int ctgread()	Read ASCII catalog stars in specified region
 * int ctgrnum()	Read ASCII catalog stars with specified numbers
 * int ctgopen()	Open ASCII catalog, returning number of entries
 * int ctgstar()	Get ASCII catalog entry for one star
 * int ctgsize()	Return length of file in bytes
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include "fitshead.h"
#include "wcs.h"
#include "wcscat.h"

/* default pathname for catalog,  used if catalog file not found in current
   working directory, but overridden by WCS_CATDIR environment variable */
char catdir[64]="/data/catalogs";

static double ctg2ra();
static double ctg2dec();
static int ctgsize();
double dt2ep();		/* Julian Date to epoch (fractional year) */

static char newline = 10;


/* CTGREAD -- Read ASCII stars in specified region */

int
ctgread (catfile, refcat, distsort, cra, cdec, dra, ddec, drad,
	 sysout, eqout, epout, mag1, mag2, nsmax,
	 tnum, tra, tdec, tpra, tpdec, tmag, tmagb, tc, tobj, nlog)

char	*catfile;	/* Name of reference star catalog file */
int	refcat;		/* Catalog code from wcctg.h */
int	distsort;	/* 1 to sort stars by distance from center */
double	cra;		/* Search center J2000 right ascension in degrees */
double	cdec;		/* Search center J2000 declination in degrees */
double	dra;		/* Search half width in right ascension in degrees */
double	ddec;		/* Search half-width in declination in degrees */
double	drad;		/* Limiting separation in degrees (ignore if 0) */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
double	mag1,mag2;	/* Limiting magnitudes (none if equal) */
int	nsmax;		/* Maximum number of stars to be returned */
double	*tnum;		/* Array of UJ numbers (returned) */
double	*tra;		/* Array of right ascensions (returned) */
double	*tdec;		/* Array of declinations (returned) */
double	*tpra;		/* Array of right ascension proper motions (returned) */
double	*tpdec;		/* Array of declination proper motions (returned) */
double	*tmag;		/* Array of magnitudes (returned) */
double	*tmagb;		/* Array of second magnitudes (returned) */
int	*tc;		/* Array of fluxes (returned) */
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
    int nameobj;	/* Save object name if 1, else do not */
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

    nstar = 0;

    /* Call the appropriate search program if not TDC ASCII catalog */
    if (refcat != TXTCAT) {
        if (refcat == GSC)
            nstar = gscread (cra,cdec,dra,ddec,drad,distsort,
			     sysout,eqout,epout,mag1,mag2,nsmax,
			     tnum,tra,tdec,tmag,tc,nlog);
        else if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
                 refcat == UAC  || refcat == UA1  || refcat == UA2)
            nstar = uacread (catfile,distsort,
                          cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,
                          mag2,nsmax,tnum,tra,tdec,tmag,tmagb,tc,nlog);
        else if (refcat == UJC)
            nstar = ujcread (cra,cdec,dra,ddec,drad,distsort,
			     sysout,eqout,epout,mag1,mag2,nsmax,
			     tnum,tra,tdec,tmag,tc,nlog);
        else if (refcat == ACT)
            nstar = actread (cra,cdec,dra,ddec,drad,distsort,
			     sysout,eqout,epout,mag1, mag2,nsmax,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,nlog);
        else if (refcat == TYCHO2)
            nstar = ty2read (cra,cdec,dra,ddec,drad,distsort,
			     sysout,eqout,epout,mag1,mag2,nsmax,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
        else if (refcat == SAO)
            nstar = binread ("SAOra", distsort,cra,cdec,dra,ddec,drad,
			     sysout,eqout,epout,mag1,mag2,nsmax,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
        else if (refcat == PPM)
            nstar = binread ("PPMra",distsort,cra,cdec,dra,ddec,drad,
			     sysout,eqout,epout,mag1,mag2,nsmax,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
        else if (refcat == IRAS)
            nstar = binread ("IRAS", distsort, cra,cdec,dra,ddec,drad,
			     sysout,eqout,epout,mag1,mag2,nsmax,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
        else if (refcat == TYCHO)
            nstar = binread ("tychora", distsort,cra,cdec,dra,ddec,drad,
			     sysout,eqout,epout,mag1,mag2,nsmax,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
        else if (refcat == HIP)
            nstar = binread ("hipparcosra", distsort, cra,cdec,dra,ddec,drad,
			     sysout,eqout,epout,mag1,mag2,nsmax,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
        else if (refcat == BSC)
            nstar = binread ("BSC5ra", distsort, cra,cdec,dra,ddec,drad,
			     sysout,eqout,epout,mag1,mag2,nsmax,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
        else if (refcat == BINCAT)
            nstar = binread (catfile, distsort, cra,cdec,dra,ddec,drad,
			     sysout,eqout,epout,mag1,mag2,nsmax,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,tobj,nlog);
        else if (refcat == TABCAT)
            nstar = tabread (catfile, distsort,cra,cdec,dra,ddec,drad,
			     sysout,eqout,epout,mag1,mag2,nsmax,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,tobj,nlog);
	return (nstar);
	}

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
    tdist = (double *) malloc (nsmax * sizeof (double));

    /* Open catalog file */
    starcat = ctgopen (catfile, refcat);
    if (starcat == NULL)
	return (0);
    if (starcat->nstars <= 0) {
	free (starcat);
	if (star != NULL)
	    free (star);
	return (0);
	}

    jstar = 0;
    if (tobj == NULL || starcat->ignore)
	nameobj = 0;
    else
	nameobj = 1;

    /* Loop through catalog */
    for (istar = 1; istar <= starcat->nstars; istar++) {
	if (ctgstar (istar, starcat, star)) {
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
	if (drad > 0 || distsort)
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
	    if (nstar < nsmax) {
		tnum[nstar] = num;
		tra[nstar] = ra;
		tdec[nstar] = dec;
		tmag[nstar] = mag;
		tdist[nstar] = dist;
		if (nameobj) {
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

	    /* If too many stars and distance sorting, replace furthest star */
	    else if (distsort) {
		if (dist < maxdist) {
		    tnum[farstar] = num;
		    tra[farstar] = ra;
		    tdec[farstar] = dec;
		    tmag[farstar] = mag;
		    tdist[farstar] = dist;
		    if (nameobj) {
			free (tobj[farstar]);
			lname = strlen (star->objname) + 1;
			objname = (char *)calloc (lname, 1);
			strcpy (objname, star->objname);
			tobj[farstar] = objname;
			}

		    /* Find new farthest star */
		    maxdist = 0.0;
		    for (i = 0; i < nsmax; i++) {
			if (tdist[i] > maxdist) {
			    maxdist = tdist[i];
			    farstar = i;
			    }
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
		if (nameobj) {
		    free (tobj[faintstar]);
		    lname = strlen (star->objname) + 1;
		    objname = (char *)calloc (lname, 1);
		    strcpy (objname, star->objname);
		    tobj[faintstar] = objname;
		    }
		faintmag = 0.0;

		/* Find new faintest star */
		for (i = 0; i < nsmax; i++) {
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
	if (nstar > nsmax)
	    fprintf (stderr,"CATREAD: %d stars found; only %d returned\n",
		     nstar,nsmax);
	}

    ctgclose (starcat, TXTCAT);
    free ((char *)tdist);
    free (star);
    return (nstar);
}


/* CTGRNUM -- Read ASCII stars with specified numbers */

int
ctgrnum (catfile,refcat, nnum,sysout,eqout,epout,match,
	 tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,tobj,nlog)

char	*catfile;	/* Name of reference star catalog file */
int	refcat;		/* Catalog code from wcctg.h */
int	nnum;		/* Number of stars to look for */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
int	match;		/* 1 to match star number exactly, else sequence num.*/
double	*tnum;		/* Array of star numbers to look for */
double	*tra;		/* Array of right ascensions (returned) */
double	*tdec;		/* Array of declinations (returned) */
double	*tpra;		/* Array of right ascension proper motions (returned) */
double	*tpdec;		/* Array of declination proper motions (returned) */
double	*tmag;		/* Array of magnitudes (returned) */
double	*tmagb;		/* Array of second magnitudes (returned) */
int	*tc;		/* Array of fluxes (returned) */
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
    int nameobj;

    nstar = 0;

    /* Call the appropriate search program if not TDC ASCII catalog */
    if (refcat != TXTCAT) {
	if (refcat == GSC)
	    nstar = gscrnum (nnum,sysout,eqout,epout,
			     tnum,tra,tdec,tmag,tc,nlog);
	else if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
	         refcat == UAC  || refcat == UA1  || refcat == UA2)
	    nstar = uacrnum (catfile,nnum,sysout,eqout,epout,
			     tnum,tra,tdec,tmag,tmagb,tc,nlog);
	else if (refcat == UJC)
	    nstar = ujcrnum (nnum,sysout,eqout,epout,
			     tnum,tra,tdec,tmag,tc,nlog);
	else if (refcat == SAO)
	    nstar = binrnum ("SAO",nnum,sysout,eqout,epout,match,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
	else if (refcat == PPM)
	    nstar = binrnum ("PPM",nnum,sysout,eqout,epout,match,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
	else if (refcat == IRAS)
	    nstar = binrnum ("IRAS",nnum,sysout,eqout,epout,match,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
	else if (refcat == TYCHO)
	    nstar = binrnum ("tycho",nnum,sysout,eqout,epout,match,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
	else if (refcat == HIP)
	    nstar = binrnum ("hipparcos",nnum,sysout,eqout,epout,match,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
	else if (refcat == BSC)
	    nstar = binrnum ("BSC5",nnum,sysout,eqout,epout,match,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,NULL,nlog);
	else if (refcat == ACT)
	    nstar = actrnum (nnum,sysout,eqout,epout,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,nlog);
	else if (refcat == TYCHO2)
	    nstar = ty2rnum (nnum,sysout,eqout,epout,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,nlog);
	else if (refcat == TABCAT)
	    nstar = tabrnum (catfile,nnum,sysout,eqout,epout,
			     tnum,tra,tdec,tpra,tpdec,tmag,tmagb,tc,tobj,nlog);
	else if (refcat == BINCAT)
	    nstar = binrnum (catfile,nnum,sysout,eqout,epout,match,
			     tnum,tra,tdec,tmag,tmagb,tc,tobj,nlog);
	return (nstar);
	}

    if ((starcat = ctgopen (catfile, refcat)) == NULL) {
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
    if (tobj == NULL || starcat->ignore)
	nameobj = 0;
    else
	nameobj = 1;

    /* Loop through star list */
    for (jnum = 0; jnum < nnum; jnum++) {

	/* Loop through catalog to star */
	starfound = 0;
	if (match && starcat->stnum > 0 && starcat->stnum < 5) {
	    for (istar = 1; istar <= starcat->nstars; istar++) {
		if (ctgstar (istar, starcat, star)) {
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
	    if (ctgstar (istar, starcat, star)) {
		fprintf (stderr,"CATRNUM: Cannot read star %d\n", istar);
		continue;
		}
	    starfound = 1;
	    }

	/* If star has been found in catalog */
	if (starfound) {

	    /* Extract selected fields  */
	    num = star->num;
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
	    if (nameobj) {
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

    ctgclose(starcat, TXTCAT);
    free (star);
    return (nstar);
}

/* CTGOPEN -- Open ASCII catalog, returning number of entries */

struct StarCat *
ctgopen (catfile, refcat)

char	*catfile;	/* ASCII catalog file name */
int	refcat;		/* Catalog code from wcctg.h (TXTCAT,BINCAT,TABCAT) */
{
    struct StarCat *sc;
    struct Tokens tokens;
    FILE *fcat;
    char header[80];
    char catpath[128];	/* Full pathname for catalog file */
    char *str;
    int nr, lfile, lhead, ldesc;
    char *catnew, *catdesc;
    char ctemp, *line, *linend, *cdot;
    char token[80];
    int ntok;
    int ltok;

    if (refcat != TXTCAT) {
	if (refcat == BINCAT)
	    sc = binopen (catfile);
	else if (refcat == TABCAT)
	    sc = tabcatopen (catfile);
	else
	    sc = NULL;
	return (sc);
	}

/* Find length of ASCII catalog */
    lfile = ctgsize (catfile);

    /* If catalog is not in current directory, look elsewhere */
    if (lfile < 2) {

	/* Prepend directory name file not in working directory */
	if ((str = getenv("WCS_CATDIR")) != NULL )
	    strcpy (catpath, str);
	else
	    strcpy (catpath, catdir);
	strcat (catpath, "/");
	strcat (catpath, catfile);
	lfile = ctgsize (catpath);
	if (lfile < 2) {
	    fprintf (stderr,"CATOPEN: ASCII catalog %s has no entries\n",catfile);
	    return (NULL);
	    }
	}
    else
	strcpy (catpath, catfile);

    /* Allocate catalog data structure */
    sc = (struct StarCat *) calloc (1, sizeof (struct StarCat));

    /* Open ASCII catalog */
    if (!(fcat = fopen (catpath, "r"))) {
	fprintf (stderr,"CATOPEN: ASCII catalog %s cannot be read\n",catfile);
	return (NULL);
	}

    /* Allocate buffer to hold entire catalog and read it */
    if ((sc->catbuff = malloc (lfile+2)) == NULL) {
	fprintf (stderr,"CATOPEN: Cannot allocate memory for ASCII catalog %s\n",
		 catfile);
	fclose (fcat);
	return (NULL);
	}
    sc->catbuff[lfile] = (char) 0;
    sc->catbuff[lfile+1] = (char) 0;

    /* Read entire catalog into memory at once */
    nr = fread (sc->catbuff, 1, lfile, fcat);
    if (nr < lfile) {
	fprintf (stderr,"CATOPEN: read only %d / %d bytes of file %s\n",
		 nr, lfile, catfile);
	(void) fclose (fcat);
	return (NULL);
	}

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
    if (lhead > 79)
	lhead = 79;
    if (sc->catbuff[0] == '#' && sc->catbuff[1] == ' ') {
	strncpy (header, sc->catbuff+2, lhead-2);
	header[lhead-2] = (char) 0;
	}
    else if (sc->catbuff[0] == '#') {
	strncpy (header, sc->catbuff+1, lhead-1);
	header[lhead-1] = (char) 0;
	}
    else {
	strncpy (header, sc->catbuff, lhead);
	header[lhead] = (char) 0;
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

    /* Catalog positions are in degrees */
    if (strsrch (header, "/d") || strsrch (header, "/D"))
	sc->inform = 'D';

    /* Ignore information after position and magnitude */
    if (strsrch (header, "/i") || strsrch (header, "/I"))
	sc->ignore = 1;
    else
	sc->ignore = 0;

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

    /* RA-sorted catalog */
    if (strsrch (header, "/r") || strsrch (header, "/R"))
	sc->rasorted = 1;

    /* Table format (hh mm ss dd mm ss) */
    if (strsrch (header, "/t") || strsrch (header, "/T"))
	sc->inform = 'T';

    /* Epoch per catalog object as yyyy.mmdd */
    if (strsrch (header, "/y") || strsrch (header, "/Y"))
	sc->nepoch = 1;
    else
	sc->nepoch = 0;

    /* Second line is description */
    sc->catdata = strchr (catdesc, newline) + 1;
    ldesc = sc->catdata - catdesc;
    if (ldesc > 63)
	ldesc = 63;
    if (catdesc[0] == '#' && catdesc[1] == ' ') {
	strncpy (sc->isname, catdesc+2, ldesc-2);
	sc->isname[ldesc-2] = (char) 0;
	}
    else if (catdesc[0] == '#') {
	strncpy (sc->isname, catdesc+1, ldesc-1);
	sc->isname[ldesc-1] = (char) 0;
	}
    else {
	strncpy (sc->isname, catdesc, ldesc);
	sc->isname[ldesc] = (char) 0;
	}

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


/* CTGCLOSE -- Close ASCII catalog and free associated data structures */

void
ctgclose (sc, refcat)

struct	StarCat *sc;
int	refcat;		/* Catalog code from wcctg.h (TXTCAT,BINCAT,TABCAT) */
{
    if (refcat == BINCAT)
	binclose (sc);
    else if (refcat == TABCAT)
	tabcatclose (sc);
    else if (refcat == TXTCAT) {
	free (sc->catbuff);
	free (sc);
	}
    return;
}



/* CTGSTAR -- Get ASCII catalog entry for one star; return 0 if successful */

int
ctgstar (istar, sc, st)

int istar;	/* Star sequence number in ASCII catalog */
struct StarCat *sc; /* Star catalog data structure */
struct Star *st; /* Star data structure, updated on return */
{
    struct Tokens tokens;
    double ydate;
    char *line;
    char *nextline;
    char token[80];
    char ctemp, *linend;
    int ntok, itok;
    int ltok;
    int imag;

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
    st->objname[0] = (char) 0;

    /* Extract information from line of catalog */
    ntok = setoken (&tokens, line, NULL);

    /* Source number */
    if (sc->stnum > 0 && sc->stnum < 5) {
	ltok = nextoken (&tokens, token);
	if (ltok > 0) {
	    st->num = atof (token);
	    if (st->num - ((double) ((int) st->num)) < 0.000001)
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

    /* Translate 3-token right ascension (hh mm ss.ss) */
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

    /* Translate single-token right ascension as degrees */
    else if (sc->inform == 'D')
	st->ra = atof (token);

    /* Translate single-token right ascension (hh:mm:ss.ss or hh.mmssss) */
    else
	st->ra = ctg2ra (token);

    /* Declination or latitude */
    ltok = nextoken (&tokens, token);
    if (ltok < 1)
	return (-1);

    /* Translate 3-token declination (sdd mm ss.ss) */
    if (sc->inform == 'T') {
	int deg, min;
	int decsgn = 0;
	double sec;

	deg = (int) (atof (token) + 0.5);
	if (strchr (token, '-') != NULL)
	    decsgn = 1;
	ltok = nextoken (&tokens, token);
	if (ltok < 1)
	    return (-1);
	min = (int) (atof (token) + 0.5);
	ltok = nextoken (&tokens, token);
	if (ltok < 1)
	    return (-1);
	sec = atof (token);
	st->dec = (double) deg + ((double) min / 60.0) + (sec / 3600.0);
	if (deg > 0 && decsgn)
	    st->dec = -st->dec;
	}

    /* Translate single-token declination as degrees */
    else if (sc->inform == 'D')
	st->dec = atof (token);

    /* Translate single-token declination (dd:mm:ss.ss or dd.mmssss) */
    else
	st->dec = ctg2dec (token);

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
    if (sc->stnum < 5 && itok < ntok && !sc->ignore) {
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


/* CTGSIZE -- return size of ASCII catalog file in bytes */

static int
ctgsize (filename)

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


/* CTG2RA -- Read the right ascension, ra, in sexagesimal hours from in[] */

static double
ctg2ra (in)

char	*in;	/* Character string */

{
    double ra;	/* Right ascension in degrees (returned) */

    ra = ctg2dec (in);
    ra = ra * 15.0;

    return (ra);
}



/* CTG2DEC -- Read the declination, dec, in sexagesimal degrees from in[] */

static double
ctg2dec (in)

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
 * Sep 13 1999	Use these subroutines for general catalog access
 * Sep 16 1999	Fix bug which didn't always return closest stars
 * Sep 16 1999	Add distsort argument so brightest stars in circle works, too
 * Oct  5 1999	Move token subroutines to catutil.c
 * Oct 15 1999	Fix calls to catclose(); eliminate dno in catstar()
 * Oct 22 1999	Fix declarations after lint
 * Oct 22 1999	Rename subroutines from cat* to ctg* to avoid system conflict
 * Nov 16 1999	Transfer dec degree sign if table format
 *
 * Jan 10 2000	Add second magnitude to tabread() and tabrnum()
 * Mar 10 2000	Add proper motions to ctgread(), ctgrnum(), tabread(), tabrnum()
 * Mar 15 2000	Add proper motions to binread(), binrnum(), actread(), actrnum()
 * Apr  3 2000	Implement fractional degrees (/d) for positions
 * Apr  3 2000	Add /i option to ignore stuff at end of line
 * Apr  3 2000	Ignore leading # on first two lines of ASCII catalog file
 * May 20 2000	Add Tycho 2 catalog support
 */
