

/* CTGRDATE -- Read ASCII stars with specified date range */

int
ctgrdate (catfile,refcat,sysout,eqout,epout,match,starcat,date1,date2,
	 nmax,tnum,tra,tdec,tpra,tpdec,tmag,tc,tobj,nlog)

char	*catfile;	/* Name of reference star catalog file */
int	refcat;		/* Catalog code from wcscat.h */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
int	match;		/* 1 to match star number exactly, else sequence num.*/
struct StarCat **starcat; /* Catalog data structure */
double	date1;		/* Start time as Modified Julian Date or Julian Date */
double	date2;		/* End time as Modified Julian Date or Julian Date */
int	nmax;		/* Maximum number of stars to look for */
double	*tnum;		/* Array of star numbers to look for */
double	*tra;		/* Array of right ascensions (returned) */
double	*tdec;		/* Array of declinations (returned) */
double	*tpra;		/* Array of right ascension proper motions (returned) */
double	*tpdec;		/* Array of declination proper motions (returned) */
double	**tmag;		/* 2-D Array of magnitudes (returned) */
int	*tc;		/* Array of fluxes (returned) */
char	**tobj;		/* Array of object names (returned) */
int	nlog;
{
    int jnum;
    int nstar;
    double ra,dec;
    double rapm, decpm;
    int istar;
    int sysref;		/* Catalog coordinate system */
    double eqref;	/* Catalog equinox */
    double epref;	/* Catalog epoch */
    double epoch1, epoch2;
    char cstr[32];
    struct StarCat *sc;
    struct Star *star;
    char *objname;
    int imag;
    int lname;
    int starfound;
    int nameobj;

    nstar = 0;

    /* Call the appropriate search program if not TDC ASCII catalog */
    if (refcat != TXTCAT) {
	if (refcat == TABCAT || refcat == WEBCAT)
	    nstar = tabrdate (catfile,sysout,eqout,epout,starcat,date1,date2,
			     nmax,tnum,tra,tdec,tpra,tpdec,tmag,tc,tobj,nlog);
	else if (refcat == BINCAT)
	    nstar = binrdate (catfile,nnum,sysout,eqout,epout,date1,date2,
			     nmax,tnum,tra,tdec,tpra,tpdec,tmag,tc,tobj,nlog);
	return (nstar);
	}

    sc = *starcat;
    if (sc == NULL) {
	if ((sc = ctgopen (catfile, refcat)) == NULL) {
	    fprintf (stderr,"CTGRDATE: Cannot read catalog %s\n", catfile);
	    *starcat = sc;
	    return (0);
	    }
	}
    *starcat = sc;
    sysref = sc->coorsys;
    eqref = sc->equinox;
    epref = sc->epoch;
    if (!sysout)
	sysout = sysref;
    if (!eqout)
	eqout = eqref;
    if (!epout)
	epout = epref;

    /* Convert date limits to fractional year for easy testing */
    if (date1 < 3000.0 && date1 > 0.0)
	epoch1 = dt2ep (date1, 12.0);
    else if (date1 < 100000.0)
	epoch1 = mjd2ep (date1);
    else
        epoch1 = jd2ep (date1);
    if (date2 < 3000.0 && date2 > 0.0)
	epoch2 = dt2ep (date2, 12.0);
    else if (date2 < 100000.0)
	epoch2 = mjd2ep (date2);
    else
        epoch2 = jd2ep (date2);

    /* Allocate catalog entry buffer */
    star = (struct Star *) calloc (1, sizeof (struct Star));
    star->num = 0.0;
    if (tobj == NULL || sc->ignore)
	nameobj = 0;
    else
	nameobj = 1;

    /* Loop through catalog to star from first date */
    starfound = 0;
    for (istar = 1; istar <= sc->nstars; istar++) {
	if (ctgstar (istar, sc, star)) {
	    fprintf (stderr,"CTGRDATE: Cannot read star %d\n", istar);
	    break;
	    }

	/* If before start date, skip to next star */
	if (star->epoch < epoch1) {
	    continue;
	    }

	/* If past starting date, drop out of reading loop */
	else if (star->epoch > epoch2) {
	    break;
	    }

	/* Extract selected fields  */
	ra = star->ra;
	dec = star->dec;
	rapm = star->rapm;
	decpm = star->decpm;

	/* Set coordinate system for this star */
	sysref = star->coorsys;
	eqref = star->equinox;
	epref = star->epoch;
    
	if (sc->inform != 'X') {
	    if (sc->mprop == 1)
		wcsconp (sysref, sysout, eqref, eqout, epref, epout,
			 &ra, &dec, &rapm, &decpm);
	    else
		wcscon (sysref, sysout, eqref, eqout, &ra, &dec, epout);
	    }

	/* Save star position and magnitude in table */
	tnum[jnum] = star->num;
	tra[jnum] = ra;
	tdec[jnum] = dec;
	for (imag = 0; imag < sc->nmag; imag++) {
	    if (tmag[imag] != NULL)
		tmag[imag][nstar] = star->xmag[imag];
	    }

	/* Spectral type */
	if (sc->sptype)
	    tc[jnum] = (1000 * (int) star->isp[0]) + (int)star->isp[1];

	if (nameobj) {
	    lname = strlen (star->objname) + 1;
	    if (lname > 1) {
		objname = (char *)calloc (lname, 1);
		strcpy (objname, star->objname);
		tobj[nstar] = objname;
		}
	    else
		tobj[nstar] = NULL;
	    }
	nstar++;
	if (nlog == 1)
	    fprintf (stderr,"CTGRDATE: %11.6f: %9.5f %9.5f %s %5.2f    \n",
		     star->num,ra,dec,cstr,star->xmag[0]);

	/* Log operation */
	else if (nlog > 0 && jnum%nlog == 0)
	    fprintf (stderr,"CTGRDATE: %5d / %5d / %5d sources catalog %s\r",
		     nstar,jnum,sc->nstars,catfile);

	/* End of star loop */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"CTGRDATE: Catalog %s : %d / %d found\n",
		 catfile,nstar,sc->nstars);

    free (star);
    return (nstar);
}
