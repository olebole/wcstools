/*** File libwcs/sortstar.c
 *** April 8, 2002
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 */

/* void FluxSortStars()		Sort star list based on brightness
 * int StarFluxSort()		Return brightest of two stars based on flux
 * void MagSortStars()		Sort stars list based on magnitude
 * int StarMagSort()		Return brightest of two stars based on mag.
 * void RASortStars()		Sort stars based on right ascension
 * int StarRASort()		Return star with lowest right ascension
 * void DecSortStars()		Sort stars based on declination
 * int StarDecSort()		Return star with lowest declination
 * void XSortStars()		Sort stars based on X coordinate in image
 * int StarXSort()		Return star with lowest X coordinate
 * void YSortStars()		Sort stars based on Y coordinate in image
 * int StarYSort()		Return star with lowest Y coordinate
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "wcs.h"
#include "wcscat.h"

/* structure for star lists needed for sorting */
typedef struct {
    double n;		/* Identifying number */
    double ra;		/* Right Ascension */
    double dec;		/* Declination */
    double pra;		/* Right Ascension proper motion */
    double pdec;	/* Declination proper motion */
    double m[11];	/* Magnitude */
    double b;		/* flux */
    double x;		/* Image X coordinate */
    double y;		/* Image Y coordinate */
    int    c;		/* Other 4-byte information */
    char   *obj;	/* Object name */
} StarInfo;

/* Sort image stars by decreasing flux */

void
FluxSortStars (sx, sy, sb, sc, ns)

double	*sx;
double	*sy;
double	*sb;
int	*sc;
int	ns;

{
    StarInfo *stars;
    int StarFluxSort ();
    int i;

    stars = (StarInfo *) calloc ((unsigned int)ns, sizeof(StarInfo));

    for (i = 0; i < ns; i++) {
	stars[i].x = sx[i];
	stars[i].y = sy[i];
	stars[i].b = sb[i];
	stars[i].c = sc[i];
	}

    qsort ((char *)stars, ns, sizeof(StarInfo), StarFluxSort);

    for (i = 0; i < ns; i++) {
	sx[i] = stars[i].x;
	sy[i] = stars[i].y;
	sb[i] = stars[i].b;
	sc[i] = stars[i].c;
	}

    free ((char *)stars);
    return;
}


/* StarFluxSort -- Order stars in decreasing flux called by qsort */

int
StarFluxSort (ssp1, ssp2)

void *ssp1, *ssp2;

{
    double b1 = ((StarInfo *)ssp1)->b;
    double b2 = ((StarInfo *)ssp2)->b;

    if (b2 > b1)
	return (1);
    else if (b2 < b1)
	return (-1);
    else
	return (0);
}


/* MagSortStars -- Sort image stars by increasing magnitude */
static int magsort = 0;

void
MagSortStars (sn, sra, sdec, spra, spdec, sx, sy, sm, sc, sobj, ns, nm, ms)

double *sn;		/* Identifying number */
double *sra;		/* Right Ascension */
double *sdec;		/* Declination */
double *spra;		/* Right Ascension proper motion */
double *spdec;		/* Declination proper motion */
double *sx;		/* Image X coordinate */
double *sy;		/* Image Y coordinate */
double **sm;		/* Magnitudes */
int    *sc;		/* Other 4-byte information */
char   **sobj;		/* Object name */
int	ns;		/* Number of stars to sort */
int	nm;		/* Number of magnitudes per star */
int	ms;		/* Magnitude by which to sort (1 to nmag) */

{
    StarInfo *stars;
    int i, j, hasnum, haspos, hasobj, haspm, hasxy;
    int StarMagSort();

    stars = (StarInfo *) calloc ((unsigned int)ns, sizeof(StarInfo));

    if (ms > 0 && ms <= nm)
	magsort = ms - 1;

    if (sn == NULL)
	hasnum = 0;
    else
	hasnum = 1;
    if (sra != NULL && sdec != NULL)
	haspos = 1;
    else
	haspos = 0;
    if (spra != NULL && spdec != NULL)
	haspm = 1;
    else
	haspm = 0;
    if (sx != NULL && sy != NULL)
	hasxy = 1;
    else
	hasxy = 0;
    if (sobj == NULL)
	hasobj = 0;
    else
	hasobj = 1;

    for (i = 0; i < ns; i++) {
	if (hasnum)
	    stars[i].n = sn[i];
	if (haspos) {
	    stars[i].ra = sra[i];
	    stars[i].dec = sdec[i];
	    }
	if (haspm) {
	    stars[i].pra = spra[i];
	    stars[i].pdec = spdec[i];
	    }
	if (hasxy) {
	    stars[i].x = sx[i];
	    stars[i].y = sy[i];
	    }
	for (j = 0; j < nm; j++)
	    stars[i].m[j] = sm[j][i];
	stars[i].c = sc[i];
	if (hasobj)
	    stars[i].obj = sobj[i];
	}

    qsort ((char *)stars, ns, sizeof(StarInfo), StarMagSort);

    for (i = 0; i < ns; i++) {
	if (hasnum)
	    sn[i] = stars[i].n;
	if (haspos) {
	    sra[i] = stars[i].ra;
	    sdec[i] = stars[i].dec;
	    }
	if (haspm) {
	    spra[i] = stars[i].pra;
	    spdec[i] = stars[i].pdec;
	    }
	if (hasxy) {
	    sx[i] = stars[i].x;
	    sy[i] = stars[i].y;
	    }
	for (j = 0; j < nm; j++)
	    sm[j][i] = stars[i].m[j];
	sc[i] = stars[i].c;
	if (hasobj)
	    sobj[i] = stars[i].obj;
	}

    free ((char *)stars);
    return;
}


/* StarMagSort -- Order stars in decreasing flux called by qsort */

int
StarMagSort (ssp1, ssp2)

void *ssp1, *ssp2;

{
    double b1 = ((StarInfo *)ssp1)->m[magsort];
    double b2 = ((StarInfo *)ssp2)->m[magsort];

    /* If sort magnitude is not set, check the others until one is found */
    if (b1 > 100.0)
	b1 = b1 - 100.0;
    if (b1 == 99.90)
	b1 = ((StarInfo *)ssp1)->m[0];
    if (b1 == 99.90)
	b1 = ((StarInfo *)ssp1)->m[1];
    if (b1 == 99.90)
	b1 = ((StarInfo *)ssp1)->m[2];
    if (b1 == 99.90)
	b1 = ((StarInfo *)ssp1)->m[3];

    /* If sort magnitude is not set, check the others until one is found */
    if (b2 > 100.0)
	b2 = b2 - 100.0;
    if (b2 == 99.90)
	b2 = ((StarInfo *)ssp2)->m[0];
    if (b2 == 99.90)
	b2 = ((StarInfo *)ssp2)->m[1];
    if (b2 == 99.90)
	b2 = ((StarInfo *)ssp2)->m[2];
    if (b2 == 99.90)
	b2 = ((StarInfo *)ssp2)->m[3];

    if (b2 < b1)
	return (1);
    else if (b2 > b1)
	return (-1);
    else
	return (0);
}


/* Sort image stars by increasing right ascension */

void
RASortStars (sn, sra, sdec, spra, spdec, sx, sy, sm, sc, sobj, ns, nm)

double *sn;		/* Identifying number */
double *sra;		/* Right Ascension */
double *sdec;		/* Declination */
double *spra;		/* Right Ascension proper motion */
double *spdec;		/* Declination proper motion */
double *sx;		/* Image X coordinate */
double *sy;		/* Image Y coordinate */
double **sm;		/* Magnitudes */
int    *sc;		/* Other 4-byte information */
char   **sobj;		/* Object name */
int	ns;		/* Number of stars to sort */
int	nm;		/* Number of magnitudes per star */
{
    StarInfo *stars;
    int i, j, hasnum, hasobj, haspm, hasxy;
    int StarRASort();

    stars = (StarInfo *) calloc ((unsigned int)ns, sizeof(StarInfo));

    if (sn == NULL)
	hasnum = 0;
    else
	hasnum = 1;
    if (spra != NULL && spdec != NULL)
	haspm = 1;
    if (sx != NULL && sy != NULL)
	hasxy = 1;
    else
	hasxy = 0;
    if (sobj == NULL)
	hasobj = 0;
    else
	hasobj = 1;

    for (i = 0; i < ns; i++) {
	if (hasnum)
	    stars[i].n = sn[i];
	stars[i].ra = sra[i];
	stars[i].dec = sdec[i];
	if (haspm) {
	    stars[i].pra = spra[i];
	    stars[i].pdec = spdec[i];
	    }
	if (hasxy) {
	    stars[i].x = sx[i];
	    stars[i].y = sy[i];
	    }
	for (j = 0; j < nm; j++)
	    stars[i].m[j] = sm[j][i];
	stars[i].c = sc[i];
	if (hasobj)
	    stars[i].obj = sobj[i];
	}

    qsort ((char *)stars, ns, sizeof(StarInfo), StarRASort);

    for (i = 0; i < ns; i++) {
	if (hasnum)
	    sn[i] = stars[i].n;
	sra[i] = stars[i].ra;
	sdec[i] = stars[i].dec;
	if (haspm) {
	    spra[i] = stars[i].pra;
	    spdec[i] = stars[i].pdec;
	    }
	if (hasxy) {
	    sx[i] = stars[i].x;
	    sy[i] = stars[i].y;
	    }
	for (j = 0; j < nm; j++)
	    sm[j][i] = stars[i].m[j];
	sc[i] = stars[i].c;
	if (hasobj)
	    sobj[i] = stars[i].obj;
	}

    free ((char *)stars);
    return;
}


/* Order stars in increasing right ascension (called by qsort) */

int
StarRASort (ssp1, ssp2)

void *ssp1, *ssp2;

{
    double ra1 = ((StarInfo *)ssp1)->ra;
    double ra2 = ((StarInfo *)ssp2)->ra;

    if (ra2 > ra1)
	return (-1);
    else if (ra2 < ra1)
	return (1);
    else
	return (0);
}


/* Sort image stars by increasing declination */

void
DecSortStars (sn, sra, sdec, spra, spdec, sx, sy, sm, sc, sobj, ns, nm)

double *sn;		/* Identifying number */
double *sra;		/* Right Ascension */
double *sdec;		/* Declination */
double *spra;		/* Right Ascension proper motion */
double *spdec;		/* Declination proper motion */
double *sx;		/* Image X coordinate */
double *sy;		/* Image Y coordinate */
double **sm;		/* Magnitudes */
int    *sc;		/* Other 4-byte information */
char   **sobj;		/* Object name */
int	ns;		/* Number of stars to sort */
int	nm;		/* Number of magnitudes per star */
{
    StarInfo *stars;
    int i, j, hasnum, hasobj, haspm, hasxy;
    int StarDecSort ();

    stars = (StarInfo *) calloc ((unsigned int)ns, sizeof(StarInfo));

    if (sn == NULL)
	hasnum = 0;
    else
	hasnum = 1;
    if (spra != NULL && spdec != NULL)
	haspm = 1;
    else
	haspm = 0;
    if (sx != NULL && sy != NULL)
	hasxy = 1;
    else
	hasxy = 0;
    if (sobj == NULL)
	hasobj = 0;
    else
	hasobj = 1;

    for (i = 0; i < ns; i++) {
	if (hasnum)
	    stars[i].n = sn[i];
	stars[i].ra = sra[i];
	stars[i].dec = sdec[i];
	if (haspm) {
	    stars[i].pra = spra[i];
	    stars[i].pdec = spdec[i];
	    }
	if (hasxy) {
	    stars[i].x = sx[i];
	    stars[i].y = sy[i];
	    }
	for (j = 0; j < nm; j++)
	    stars[i].m[j] = sm[j][i];
	stars[i].c = sc[i];
	if (hasobj)
	    stars[i].obj = sobj[i];
	}

    qsort ((char *)stars, ns, sizeof(StarInfo), StarDecSort);

    for (i = 0; i < ns; i++) {
	if (hasnum)
	    sn[i] = stars[i].n;
	sra[i] = stars[i].ra;
	sdec[i] = stars[i].dec;
	if (haspm) {
	    spra[i] = stars[i].pra;
	    spdec[i] = stars[i].pdec;
	    }
	if (hasxy) {
	    sx[i] = stars[i].x;
	    sy[i] = stars[i].y;
	    }
	for (j = 0; j < nm; j++)
	    sm[j][i] = stars[i].m[j];
	sc[i] = stars[i].c;
	if (hasobj)
	    sobj[i] = stars[i].obj;
	}

    free ((char *)stars);
    return;
}


/* Order stars in increasing declination (called by qsort) */

int
StarDecSort (ssp1, ssp2)

void *ssp1, *ssp2;

{
    double dec1 = ((StarInfo *)ssp1)->dec;
    double dec2 = ((StarInfo *)ssp2)->dec;

    if (dec2 > dec1)
	return (-1);
    else if (dec2 < dec1)
	return (1);
    else
	return (0);
}


/* XSortStars -- Sort image stars by increasing X value */

void
XSortStars (sn, sra, sdec, spra, spdec, sx, sy, sm, sc, sobj, ns, nm)

double *sn;		/* Identifying number */
double *sra;		/* Right Ascension */
double *sdec;		/* Declination */
double *spra;		/* Right Ascension proper motion */
double *spdec;		/* Declination proper motion */
double *sx;		/* Image X coordinate */
double *sy;		/* Image Y coordinate */
double **sm;		/* Magnitudes */
int    *sc;		/* Other 4-byte information */
char   **sobj;		/* Object name */
int	ns;		/* Number of stars to sort */
int	nm;		/* Number of magnitudes per star */
{
    StarInfo *stars;
    int i, j, hasnum, hasobj, haspos, haspm;
    int StarXSort ();

    stars = (StarInfo *) calloc ((unsigned int)ns, sizeof(StarInfo));
    if (sn == NULL)
	hasnum = 0;
    else
	hasnum = 1;
    if (sra != NULL && sdec != NULL)
	haspos = 1;
    else
	haspos = 0;
    if (spra != NULL && spdec != NULL)
	haspm = 1;
    else
	haspm = 0;
    if (sobj == NULL)
	hasobj = 0;
    else
	hasobj = 1;

    for (i = 0; i < ns; i++) {
	if (hasnum)
	    stars[i].n = sn[i];
	if (haspos) {
	    stars[i].ra = sra[i];
	    stars[i].dec = sdec[i];
	    }
	if (haspm) {
	    stars[i].pra = spra[i];
	    stars[i].pdec = spdec[i];
	    }
	stars[i].x = sx[i];
	stars[i].y = sy[i];
	for (j = 0; j < nm; j++)
	    stars[i].m[j] = sm[j][i];
	stars[i].c = sc[i];
	if (hasobj)
	    stars[i].obj = sobj[i];
	}

    qsort ((char *)stars, ns, sizeof(StarInfo), StarXSort);

    for (i = 0; i < ns; i++) {
	if (hasnum)
	    sn[i] = stars[i].n;
	if (haspos) {
	    sra[i] = stars[i].ra;
	    sdec[i] = stars[i].dec;
	    }
	if (haspm) {
	    spra[i] = stars[i].pra;
	    spdec[i] = stars[i].pdec;
	    }
	sx[i] = stars[i].x;
	sy[i] = stars[i].y;
	for (j = 0; j < nm; j++)
	    sm[j][i] = stars[i].m[j];
	sc[i] = stars[i].c;
	if (hasobj)
	    sobj[i] = stars[i].obj;
	}

    free ((char *)stars);
    return;
}


/* StarXSort -- Order stars in decreasing X value called by qsort */

int
StarXSort (ssp1, ssp2)

void *ssp1, *ssp2;

{
    double x1 = ((StarInfo *)ssp1)->x;
    double x2 = ((StarInfo *)ssp2)->x;

    if (x2 < x1)
	return (1);
    else if (x2 > x1)
	return (-1);
    else
	return (0);
}


/* YSortStars -- Sort image stars by increasing Y value */

void
YSortStars (sn, sra, sdec, spra, spdec, sx, sy, sm, sc, sobj, ns, nm)

double *sn;		/* Identifying number */
double *sra;		/* Right Ascension */
double *sdec;		/* Declination */
double *spra;		/* Right Ascension proper motion */
double *spdec;		/* Declination proper motion */
double *sx;		/* Image X coordinate */
double *sy;		/* Image Y coordinate */
double **sm;		/* Magnitudes */
int    *sc;		/* Other 4-byte information */
char   **sobj;		/* Object name */
int	ns;		/* Number of stars to sort */
int	nm;		/* Number of magnitudes per star */
{
    StarInfo *stars;
    int i, j, hasnum, hasobj, haspm, haspos;
    int StarYSort ();

    stars = (StarInfo *) calloc ((unsigned int)ns, sizeof(StarInfo));
    if (sn == NULL)
	hasnum = 0;
    else
	hasnum = 1;
    if (sra != NULL && sdec != NULL)
	haspos = 1;
    else
	haspos = 0;
    if (spra != NULL && spdec != NULL)
	haspm = 1;
    else
	haspm = 0;
    if (sobj == NULL)
	hasobj = 0;
    else
	hasobj = 1;

    for (i = 0; i < ns; i++) {
	if (hasnum)
	    stars[i].n = sn[i];
	if (haspos) {
	    stars[i].ra = sra[i];
	    stars[i].dec = sdec[i];
	    }
	if (haspm) {
	    stars[i].pra = spra[i];
	    stars[i].pdec = spdec[i];
	    }
	stars[i].x = sx[i];
	stars[i].y = sy[i];
	for (j = 0; j < nm; j++)
	    stars[i].m[j] = sm[j][i];
	stars[i].c = sc[i];
	if (hasobj)
	    stars[i].obj = sobj[i];
	}

    qsort ((char *)stars, ns, sizeof(StarInfo), StarYSort);

    for (i = 0; i < ns; i++) {
	if (hasnum)
	    sn[i] = stars[i].n;
	if (haspos) {
	    sra[i] = stars[i].ra;
	    sdec[i] = stars[i].dec;
	    }
	if (haspm) {
	    spra[i] = stars[i].pra;
	    spdec[i] = stars[i].pdec;
	    }
	sx[i] = stars[i].x;
	sy[i] = stars[i].y;
	for (j = 0; j < nm; j++)
	    sm[j][i] = stars[i].m[j];
	sc[i] = stars[i].c;
	if (hasobj)
	    sobj[i] = stars[i].obj;
	}

    free ((char *)stars);
    return;
}


/* StarYSort -- Order stars in decreasing Y value called by qsort */

int
StarYSort (ssp1, ssp2)

void *ssp1, *ssp2;

{
    double y1 = ((StarInfo *)ssp1)->y;
    double y2 = ((StarInfo *)ssp2)->y;

    if (y2 < y1)
	return (1);
    else if (y2 > y1)
	return (-1);
    else
	return (0);
}

/* Jun 13 1996	New program
 * Oct 18 1996	Add sorting by X value
 * Nov 13 1996	Add second magnitude
 * Jan 10 1997	Fix bug in RASortStars to return correct red magnitude
 *
 * Mar  2 1998	Make number and second magnitude optional
 * Oct 21 1998	Add RefCat() to set reference catalog code
 * Oct 26 1998	Include object names in star catalog entry structure
 * Oct 29 1998	Return coordinate system and title from RefCat
 * Nov 20 1998	Add USNO A-2.0 catalog and return different code
 * Dec  9 1998	Add Hipparcos and Tycho catalogs
 *
 * Jan 26 1999	Add subroutines to deal with ranges of numbers
 * Feb  8 1999	Fix bug initializing ACT catalog
 * Feb 11 1999	Change starcat.insys to starcat.coorsys
 * May 19 1999	Move catalog subroutines to catutil()
 * Aug 26 1999	Compare pointers to NULL, not 0
 *
 * Mar 14 2000	Add proper motions
 *
 * May 22 2001	Add sort by declination
 * Jun 28 2001	In MagSort, if b mag is 99.9, try r mag
 * Jul 20 2001	In MagSort, allow for absence of ra and dec
 * Sep 12 2001	Allow up to 11 magnitudes; add nm and magsort
 * Sep 13 2001	Add YSortStars() to sort by Y coordinate
 * Sep 18 2001	Subtract 100 in MagSort if magnitude is greater than 100
 * Nov  6 2001	Allow missing x and y in MagSort, RASort, and DecSort
 * Nov  6 2001	Allow missing ra and dec in XSort and YSort
 *
 * Apr  8 2002	Drop static subroutine declarations
 */
