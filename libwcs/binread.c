/*** File libwcs/binread.c
 *** October 21, 1999
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

/* int binread()	Read binary catalog sources + names in specified region
 * int binrnum()	Read binary catalog sources + names with specified numbers
 * int binopen()	Open binary catalog, returning number of entries
 * int binstar()	Get binary catalog entry for one source
 * void binclose()	Close binary catalog
 */

/* default pathname for catalog,  used if catalog file not found in current
   working directory, but overridden by WCS_BINDIR environment variable */
char bindir[64]="/data/catalogs";

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include "fitshead.h"
#include "wcs.h"
#include "wcscat.h"

static int binsra();
void binclose();
static void binswap8();
static void binswap4();
static void binswap2();


/* BINREAD -- Read binary catalog sources + names in specified region */

int
binread (bincat,distsort,cra,cdec,dra,ddec,drad,sysout,eqout,epout,mag1,mag2,
	 nstarmax,tnum,tra,tdec,tmag,tmagb,tpeak,tobj,nlog)

char	*bincat;	/* Name of reference star catalog file */
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
int	nstarmax;	/* Maximum number of sources to be returned */
double	*tnum;		/* Array of catalog numbers (returned) */
double	*tra;		/* Array of right ascensions (returned) */
double	*tdec;		/* Array of declinations (returned) */
double	*tmag;		/* Array of V magnitudes (returned) */
double	*tmagb;		/* Array of B magnitudes (returned) */
int	*tpeak;		/* Array of peak counts (returned) */
char	**tobj;		/* Array of object names (returned) */
int	nlog;
{
    double rra1,rra2;	/* Limiting catalog right ascensions of region */
    double rdec1,rdec2;	/* Limiting catalog declinations of region */
    double ra1,ra2;	/* Limiting output right ascensions of region */
    double dec1,dec2;	/* Limiting output declinations of region */
    double dist = 0.0;  /* Distance from search center in degrees */
    double faintmag=0.0; /* Faintest magnitude */
    double maxdist=0.0; /* Largest distance */
    int faintstar=0;    /* Faintest star */
    int farstar=0;      /* Most distant star */
    double *tdist;      /* Array of distances to sources from search center */
    int sysref;		/* Catalog coordinate system */
    double eqref;	/* Catalog coordinate equinox */
    double epref;	/* Catalog position epoch */
    double ra, dec, rapm, decpm;
    double rra1a, rra2a;
    struct StarCat *starcat;
    struct Star *star;
    int rwrap, wrap, iwrap, istar1,istar2;
    char *objname;
    int lname;
    int jstar;
    int nstar;
    double mag, magb;
    double num;
    int i;
    int istar;
    int isp;
    int verbose;
    char cstr[16];

    star = NULL;
    starcat = NULL;

    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    /* Keep mag1 the smallest magnitude */
    if (mag2 < mag1) {
	mag = mag2;
	mag2 = mag1;
	mag1 = mag;
	}

    /* Logging interval */
    nstar = 0;
    tdist = (double *) malloc (nstarmax * sizeof (double));

    /* Open catalog file */
    starcat = binopen (bincat);
    if (starcat == NULL)
	return (0);
    if (starcat->nstars <= 0) {
	free (starcat);
	if (star != NULL)
	    free (star);
	return (0);
	}

    SearchLim (cra, cdec, dra, ddec, &ra1, &ra2, &dec1, &dec2, verbose);
  
    sysref = starcat->coorsys;
    eqref = starcat->equinox;
    epref = starcat->epoch;
    RefLim (cra, cdec, dra, ddec, sysout, sysref, eqout, eqref, epout,
	    &rra1, &rra2, &rdec1, &rdec2, verbose);
    if (verbose) {
	char rstr1[16],rstr2[16],dstr1[16],dstr2[16];
	ra2str (rstr1, 16, rra1, 3);
        dec2str (dstr1, 16, rdec1, 2);
	ra2str (rstr2, 16, rra2, 3);
        dec2str (dstr2, 16, rdec2, 2);

	wcscstr (cstr, sysref,eqref,epref);
	fprintf (stderr,"BINREAD: RA: %s - %s  Dec: %s - %s %s\n",
		 rstr1, rstr2, dstr1, dstr2, cstr);
	}

    /* If catalog RA range includes zero, split search in two */
    rwrap = 0;
    if (rra1 > rra2) {
	rwrap = 1;
	rra1a = 0.0;
	rra2a = rra2;
	rra2 = 360.0;
	}
    else
	rwrap = 0;

    /* Make sure first declination is always the smallest one */
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

    /* If output RA range includes zero, set flag */
    if (ra1 > ra2)
	wrap = 1;
    else
	wrap = 0;

    /* Allocate catalog entry buffer */
    star = (struct Star *) calloc (1, sizeof (struct Star));
    star->num = 0.0;

    jstar = 0;

    /* Loop through wraps (do not cross 360 degrees in search */
    for (iwrap = 0; iwrap <= rwrap; iwrap++) {

	/* Set first and last stars to check */
	if (starcat->rasorted) {
	    istar1 = binsra (starcat, star, rra1, 0);
	    istar2 = binsra (starcat, star, rra2, 1);
	    }
	else {
	    istar1 = starcat->star1;
	    istar2 = starcat->star0 + starcat->nstars;
	    }
	if (verbose)
	    printf ("BINREAD: Searching stars %d through %d\n",istar1,istar2);

	/* Loop through catalog */
	for (istar = istar1; istar <= istar2; istar++) {
	    if (binstar (starcat, star, istar)) {
		fprintf (stderr,"BINREAD: Cannot read star %d\n", istar);
		break;
		}

	    /* ID number */
	    num = star->num;

	    /* Get position in output coordinate system, equinox, and epoch */
	    rapm = star->rapm;
	    decpm = star->decpm;
	    ra = star->ra;
	    dec = star->dec;
	    wcsconp (sysref, sysout, eqref, eqout, epref, epout,
		     &ra, &dec, &rapm, &decpm);

	    /* Magnitude */
	    mag = 0.01 * (double) star->mag[0];
	    if (starcat->nmag > 1) {
		magb = mag;
		mag = 0.01 * (double) star->mag[1];
		}
	    else
		magb = 20.0;

	    /* Spectral Type */
	    isp = (1000 * (int) star->isp[0]) + (int)star->isp[1];

	    /* Compute distance from search center */
	    if (drad > 0 || distsort)
		dist = wcsdist (cra,cdec,ra,dec);
	    else
		dist = 0.0;

	    /* Check magnitude and position limits */
	    if ((mag1 == mag2 || (mag >= mag1 && mag <= mag2)) &&
		((wrap && (ra <= ra1 || ra >= ra2)) ||
		(!wrap && (ra >= ra1 && ra <= ra2))) &&
		((drad > 0.0 && dist <= drad) ||
     		(drad == 0.0 && dec >= dec1 && dec <= dec2))) {

		/* Save star position and magnitude in table */
		if (nstar < nstarmax) {
		    tnum[nstar] = num;
		    tra[nstar] = ra;
		    tdec[nstar] = dec;
		    tmag[nstar] = mag;
		    tmagb[nstar] = magb;
		    tpeak[nstar] = isp;
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

		/* If too many stars and distance sorting,
		   replace furthest star */
		else if (distsort) {
		    if (dist < maxdist) {
			tnum[farstar] = num;
			tra[farstar] = ra;
			tdec[farstar] = dec;
			tmag[farstar] = mag;
			tmagb[farstar] = magb;
			tpeak[farstar] = isp;
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
		    }

		/* Else if too many stars, replace faintest star */
		else if (mag < faintmag) {
		    tnum[faintstar] = num;
		    tra[faintstar] = ra;
		    tdec[faintstar] = dec;
		    tmag[faintstar] = mag;
		    tmagb[faintstar] = magb;
		    tpeak[faintstar] = isp;
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
		    fprintf (stderr,"BINREAD: %11.6f: %9.5f %9.5f %5.2f %5.2f\n",
			   num,ra,dec,magb,mag);

	    /* End of accepted star processing */
		}

	/* Log operation */
	    if (nlog > 0 && istar%nlog == 0)
		fprintf (stderr,"BINREAD: %5d / %5d / %5d sources catalog %s\r",
			jstar,istar,starcat->nstars,bincat);

	/* End of star loop */
	    }

	/* Set second set of RA limits if passing through 0h */
	rra1 = rra1a;
	rra2 = rra2a;
	}

    /* Summarize search */
    if (nlog > 0) {
	fprintf (stderr,"BINREAD: Catalog %s : %d / %d / %d found\n",
		 bincat,jstar,istar,starcat->nstars);
	if (nstar > nstarmax)
	    fprintf (stderr,"BINREAD: %d stars found; only %d returned\n",
		     nstar,nstarmax);
	}

    binclose(starcat);
    free (star);
    free ((char *)tdist);
    return (nstar);
}


/* BINRNUM -- Read binary catalog stars with specified numbers */

int
binrnum (bincat, nnum, sysout, eqout, epout, match,
	 tnum,tra,tdec,tmag,tmagb,tpeak,tobj,nlog)

char	*bincat;	/* Name of reference star catalog file */
int	nnum;		/* Number of stars to look for */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
int	match;		/* If 1, match number exactly, else number is sequence*/
double	*tnum;		/* Array of star numbers to look for */
double	*tra;		/* Array of right ascensions (returned) */
double	*tdec;		/* Array of declinations (returned) */
double	*tmag;		/* Array of V magnitudes (returned) */
double	*tmagb;		/* Array of B magnitudes (returned) */
int	*tpeak;		/* Array of peak counts (returned) */
char	**tobj;		/* Array of object names (returned) */
int	nlog;
{
    int sysref;		/* Catalog coordinate system */
    double eqref;	/* Catalog coordinate equinox */
    double epref;	/* Catalog position epoch */
    int jnum;
    int nstar;
    double ra, dec, rapm, decpm;
    double mag, magb;
    double num;
    int istar;
    int isp;
    int lname;
    char *objname;
    struct StarCat *starcat;
    struct Star *star;

    nstar = 0;
    starcat = binopen (bincat);
    if (starcat == NULL)
	return (0);
    if (starcat->nstars <= 0) {
	free (starcat);
	if (star != NULL)
	    free (star);
	fprintf (stderr,"BINRNUM: Cannot read catalog %s\n", bincat);
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

	/* Find star in catalog */
	istar = (int) tnum[jnum];
	if (match) {
	    istar = 1;
	    while (istar <= starcat->nstars) {
		if (binstar (starcat, star, istar)) {
		    fprintf (stderr,"BINRNUM: Cannot read star %d\n", istar);
		    tra[jnum] = 0.0;
		    tdec[jnum] = 0.0;
		    tmag[jnum] = 0.0;
		    tmagb[jnum] = 0.0;
		    tpeak[jnum] = 0;
		    continue;
		    }
		if (star->num == tnum[jnum])
		    break;
		istar++;
		}
	    if (star->num != tnum[jnum])
		continue;
	    }

	else if (binstar (starcat, star, istar)) {
	    fprintf (stderr,"BINRNUM: Cannot read star %d\n", istar);
	    tra[jnum] = 0.0;
	    tdec[jnum] = 0.0;
	    tmag[jnum] = 0.0;
	    tmagb[jnum] = 0.0;
	    tpeak[jnum] = 0;
	    continue;
	    }

	/* If star has been found in catalog */

	/* ID number */
	num = star->num;

	/* Position in degrees at designated epoch */
	ra = star->ra;
	dec = star->dec;
	rapm = star->rapm;
	decpm = star->decpm;
	wcsconp (sysref, sysout, eqref, eqout, epref, epout,
		     &ra, &dec, &rapm, &decpm);

	/* Magnitude */
	mag = 0.01 * (double) star->mag[0];
	if (starcat->nmag > 1) {
	    magb = mag;
	    mag = 0.01 * (double) star->mag[1];
	    }
	else
	    magb = 20.0;

	/* Spectral Type */
	isp = (1000 * (int) star->isp[0]) + (int)star->isp[1];

	/* Save star position and magnitude in table */
	tnum[jnum] = num;
	tra[jnum] = ra;
	tdec[jnum] = dec;
	tmag[jnum] = mag;
	tpeak[jnum] = isp;
	if (tobj != NULL) {
	    lname = strlen (star->objname) + 1;
	    objname = (char *)calloc (lname, 1);
	    strcpy (objname, star->objname);
	    tobj[nstar] = objname;
	    }
	nstar++;
	if (nlog == 1)
	    fprintf (stderr,"BINRNUM: %11.6f: %9.5f %9.5f %5.2f %5.2f %s  \n",
		     num, ra, dec, magb, mag, star->isp);

	/* End of star loop */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"BINRNUM: Catalog %s : %d / %d found\n",
		 bincat,nstar,starcat->nstars);

    binclose (starcat);
    return (nstar);
}


/* BINOPEN -- Open binary catalog, returning number of entries */

struct StarCat *
binopen (bincat)

char *bincat;	/* Binary catalog file name */
{
    FILE *fcat;
    struct StarCat *sc;
    int nr, lfile;
    char *binfile;
    char binpath[128];	/* Full pathname for catalog file */
    char *str;
    int lf;
    static int binsize();

    /* Find length of binary catalog */
    lfile = binsize (bincat);

    /* Check for existence of catalog */
    if (lfile < 2) {

	/* Prepend directory name file not in working directory */
	if ((str = getenv("WCS_BINDIR")) != NULL )
	    strcpy (bindir, str);
	strcpy (binpath, bindir);
	strcat (binpath, "/");
	strcat (binpath, bincat);
	lfile = binsize (binpath);
	if (lfile < 2) {
	    fprintf (stderr,"BINOPEN: Binary catalog %s has no entries\n",bincat);
	    return (0);
	    }
	}
    else
	strcpy (binpath, bincat);

    /* Open binary catalog */
    if (!(fcat = fopen (binpath, "r"))) {
	fprintf (stderr,"BINOPEN: Binary catalog %s cannot be read\n",binpath);
	return (0);
	}

    /* Read binary catalog header information */
    sc = (struct StarCat *) calloc (1, sizeof (struct StarCat));
    nr = fread (sc, 1, 28, fcat);
    if (nr < 28) {
	fprintf (stderr,"BINOPEN: read only %d / %d bytes of file %s\n",
		 nr, lfile, binpath);
	(void) fclose (fcat);
	return (0);
	}

    /* Check for byte reversal */
    if (sc->nmag > 10 || sc->nmag < -10) {
	sc->byteswapped = 1;
	binswap4 (&sc->star0);
	binswap4 (&sc->star1);
	binswap4 (&sc->nstars);
	binswap4 (&sc->stnum);
	binswap4 (&sc->mprop);
	binswap4 (&sc->nmag);
	binswap4 (&sc->nbent);
	}
    else
	sc->byteswapped = 0;

    if (sc->stnum < 0)
	sc->ncobj = -sc->stnum;
    else
	sc->ncobj = 0;

    /* Set number of decimal places in star numbers */
    if (sc->stnum == 2)
	sc->nndec = 4;
    else if (sc->stnum == 3)
	sc->nndec = 5;
    else
	sc->nndec = 0;

    strcpy (sc->incdir, bindir);
    strcpy (sc->incfile, bincat);

    /* Separate filename from pathname and save in structure */
    binfile = strrchr (binpath,'/');
    if (binfile)
	binfile = binfile + 1;
    else
	binfile = binpath;
    if (strlen (binfile) < 24)
	strcpy (sc->isfil, binfile);
    else
	strncpy (sc->isfil, binfile, 23);

    /* Set other catalog information in structure */
    if (sc->nmag < 0) {
	sc->inform = 'J';
	sc->coorsys = WCS_J2000;
	sc->epoch = 2000.0;
	sc->equinox = 2000.0;
	sc->nmag = -sc->nmag;
	}
    else if (sc->nstars < 0) {
	sc->inform = 'J';
	sc->coorsys = WCS_J2000;
	sc->epoch = 2000.0;
	sc->equinox = 2000.0;
	sc->nstars = -sc->nstars;
	}
    else {
	sc->inform = 'B';
	sc->coorsys = WCS_B1950;
	sc->epoch = 1950.0;
	sc->equinox = 1950.0;
	}
    sc->ifcat = fcat;

    /* Check name to see if file is RA-sorted */
    lf = strlen (binfile);
    sc->rasorted = 0;
    if (binfile[lf-2] == 'r' && binfile[lf-1] == 'a')
	sc->rasorted = 1;
    if (!strncmp (binfile, "IRAS", 4))
	sc->rasorted = 1;
    return (sc);
}


void
binclose (sc)
struct StarCat *sc;	/* Star catalog descriptor */
{
    fclose (sc->ifcat);
    free (sc);
    return;
}


/* BINSRA -- Find star closest to given RA in RA-sorted catalog */

static int
binsra (sc, st, dra, end)

struct StarCat *sc;	/* Star catalog descriptor */
struct Star *st;	/* Current star entry */
double	dra;		/* Right ascension in degrees */
int	end;		/* 0: start of range, 1: end of range */

{
    char rastr[16], raxstr[16], ramins[16], ramaxs[16];
    int istar0, istarx, nrep, ismin, ismax;
    double rax, ramin, ramax;
    int verbose = 0;

    /* Keep RA between 0 and 360 degrees */
    if (dra > 360.0)
	rax = dra - 360.0;
    else
	rax = dra;

    ismin = 1;
    if (binstar (sc, st, ismin)) {
	fprintf (stderr,"BINSRA: Cannot read star %d\n", ismin);
	return (0);
	}
    else
	ramin = st->ra;

    ismax = sc->nstars;
    if (binstar (sc, st, ismax)) {
	fprintf (stderr,"BINSRA: Cannot read star %d\n", ismax);
	return (0);
	}
    else
	ramax = st->ra;

    istarx = sc->nstars / 2;

    for (nrep = 0; nrep < 32; nrep++) {
	if (binstar (sc, st, istarx)) {
	    fprintf (stderr,"BINSRA: Cannot read star %d\n", istarx);
            return (0);
	    }

	/* Find next catalog number to read */
	if (st->ra < rax) {
	    ismin = istarx;
	    ramin = st->ra;
	    istar0 = istarx;
	    if (ismax - istarx > 1)
		istarx = istarx + (ismax - istarx) / 2;
	    else if (ismax - istarx > 0)
		istarx = istarx + 1;
	    }
	else if (st->ra > rax) {
	    ismax = istarx;
	    ramax = st->ra;
	    istar0 = istarx;
	    if (istarx - ismin > 1)
		istarx = istarx - ((istarx - ismin) / 2);
	    else if (istarx - ismin > 0)
		istarx = istarx - 1;
	    }
	else
	    break;

	if (verbose) {
	    ra2str (rastr, 16, st->ra, 3);
	    ra2str (raxstr, 16, rax, 3);
	    ra2str (ramins, 16, ramin, 3);
	    ra2str (ramaxs, 16, ramax, 3);
	    printf ("%9d: %s -> %s  %9d: %s  %9d: %s\n",
		    istarx, rastr, raxstr, ismin,ramins,ismax,ramaxs);
	    }
	if (istarx == istar0)
	    break;
	}

    /* Make sure final star is real */
    if (binstar (sc, st, istarx)) {
	fprintf (stderr,"BINSRA: Cannot read star %d\n", istarx);
        return (0);
	}
    else
	return (istarx);
}


/* BINSTAR -- Get binary catalog entry for one star;
              return 0 if successful */

int
binstar (sc, st, istar)

struct StarCat *sc;	/* Star catalog descriptor */
struct Star *st;	/* Current star entry */
int istar;	/* Star sequence number in binary catalog */
{
    int ino, nbent;
    long offset;
    float pm[2];

    /* Drop out if catalog pointer is not set */
    if (sc == NULL)
	return (1);

    /* Drop out if catalog is not open */
    if (sc->ifcat == NULL)
	return (2);

    /* Drop out if star number is too large */
    if (istar > sc->nstars) {
	fprintf (stderr, "BINSTAR:  %d  > %d is not in catalog\n",
		 istar, sc->nstars);
	return (3);
	}

    /* Move file pointer to start of correct star entry */
    if (istar > 0) {
	offset = 28 + (istar - sc->star1) * sc->nbent;
	if (fseek (sc->ifcat, offset, SEEK_SET))
	    return (0);
	}

    /* Read catalog entry */
    if (sc->mprop)
	nbent = sc->nbent - 8;
    else
	nbent = sc->nbent;
    if (sc->stnum < 0)
	nbent = nbent + sc->stnum;
    sc->ncobj = 0;
    if (sc->stnum <= 0) {
	sc->ncobj = -sc->stnum;
	if (fread (&st->ra, nbent, 1, sc->ifcat) < 1)
	    return (4);
	if (sc->stnum < 0) {
	    if (fread (st->objname, sc->ncobj, 1, sc->ifcat) < 1)
		return (4);
	    }
	}
    else {
	if (fread (&st->xno, nbent, 1, sc->ifcat) < 1)
	    return (4);
	if (sc->byteswapped)
	    binswap4 (&st->xno);
	}
    if (sc->byteswapped) {
	binswap8 (&st->ra);
	binswap8 (&st->dec);
	binswap2 (st->mag, sc->nmag*2);
	}

    /* Convert position to degrees */
    st->ra = raddeg (st->ra);
    st->dec = raddeg (st->dec);

    /* Read proper motion, if it is present, and convert it to degrees */
    if (sc->mprop) {
	if (fread (pm, 8, 1, sc->ifcat) < 1)
	    return (5);
	if (sc->byteswapped) {
	    binswap4 (&pm[0]);
	    binswap4 (&pm[1]);
	    }
	st->rapm = raddeg ((double) pm[0]);
	st->decpm = raddeg ((double) pm[1]);
	}

    /* Interpret catalog number */
    switch (sc->stnum) {
	case 1:
	    st->num = (double) st->xno;
	    break;
	case 2:
	    bcopy ((char *)&st->xno, (char *) &ino, 4);
	    st->num = 0.0001 * (double) ino;
	    break;
	case 3:
	    bcopy ((char *)&st->xno, (char *) &ino, 4);
	    st->num = 0.00001 * (double) ino;
	    break;
	case 4:
	    bcopy ((char *)&st->xno, (char *) &ino, 4);
	    st->num = (double) ino;
	    break;
	default:
	    if (istar > 0)
		st->num = (double) istar;
	    else
		st->num = st->num + 1.0;
	    break;
	}

    return (0);
}


/* BINSWAP2 -- Swap bytes in string in place */

static void
binswap2 (string,nbytes)


char *string;	/* Address of starting point of bytes to swap */
int nbytes;	/* Number of bytes to swap */

{
    char *sbyte, temp, *slast;

    slast = string + nbytes;
    sbyte = string;
    while (sbyte < slast) {
	temp = sbyte[0];
	sbyte[0] = sbyte[1];
	sbyte[1] = temp;
	sbyte= sbyte + 2;
	}
    return;
}


/* BINSWAP4 -- Reverse bytes of Integer*4 or Real*4 number in place */

static void
binswap4 (string)

char *string;	/* Address of Integer*4 or Real*4 vector */

{
    char temp0, temp1, temp2, temp3;

    temp3 = string[0];
    temp2 = string[1];
    temp1 = string[2];
    temp0 = string[3];
    string[0] = temp0;
    string[1] = temp1;
    string[2] = temp2;
    string[3] = temp3;

    return;
}


/* BINSWAP8 -- Reverse bytes of Real*8 vector in place */

static void
binswap8 (string)

char *string;	/* Address of Real*8 vector */

{
    char temp[8];

    temp[7] = string[0];
    temp[6] = string[1];
    temp[5] = string[2];
    temp[4] = string[3];
    temp[3] = string[4];
    temp[2] = string[5];
    temp[1] = string[6];
    temp[0] = string[7];
    string[0] = temp[0];
    string[1] = temp[1];
    string[2] = temp[2];
    string[3] = temp[3];
    string[4] = temp[4];
    string[5] = temp[5];
    string[6] = temp[6];
    string[7] = temp[7];
    return;
}


/* BINSIZE -- return size of binary catalog file in bytes */

static int
binsize (filename)

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


/* ISBIN -- Return 1 if TDC binary catalog file, else 0 */

int
isbin (filename)

char    *filename;      /* Name of file to check */
{
    FILE *diskfile;
    char line[8];
    int nbr;

    if ((diskfile = fopen (filename, "r")) == NULL)
	return (0);
    else {
	nbr = fread (line, 1, 4, diskfile);
	fclose (diskfile);
	if (nbr < 4)
	    return (0);
	else if (line[0] == 0 || line[1] == 0 || line[2] == 0 || line[3] == 0)
	    return (1);
	else
	    return (0);
	}
}

/* Sep 10 1998	New subroutines
 * Sep 15 1998	Add byte swapping
 * Sep 16 1998	Use limiting radius correctly; use arbitrary search system
 * Sep 21 1998	Return converted coordinates
 * Sep 23 1998	Set search limits in catalog system, but test in output system
 * Sep 24 1998	Return second magnitude if more than one in catalog
 * Oct 15 1998	Add binsize() to compute file length in bytes
 * Oct 20 1998	Add binrobj() to read object names for specified stars
 * Oct 21 1998	Fix binsra() to get out of 3-in-a-row identical RAs
 * Oct 21 1998	Add isbin() to check whether file could be TDC binary catalog
 * Oct 26 1998	Add object retrieval to binread() and binrnum()
 * Oct 29 1998	Correctly assign numbers when too many stars are found
 * Oct 30 1998	Be more careful with coordinate system of catalog
 * Oct 30 1998	Fix convergence at edges of catalog
 * Nov  9 1998	Add flag to catalog structure for RA-sorted catalogs
 * Dec  8 1998	Have binstar() return 0 instead of null
 *
 * Feb  1 1999	Add match argument to binrnum() 
 * Feb  2 1999	Set number of decimal places in star number
 * Feb 11 1999	Change starcat.insys to starcat.coorsys
 * Jun 16 1999	Use SearchLim()
 * Aug 16 1999	Add RefLim() to get converted search coordinates right
 * Aug 24 1999	Fix declination limit bug which broke search 
 * Aug 25 1999	Return real number of stars from binread()
 * Sep 16 1999	Fix bug which didn't always return closest stars
 * Sep 16 1999	Add distsort argument so brightest stars in circle works, too
 * Oct 21 1999	Fix declarations after lint
 */
