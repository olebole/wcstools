/* File libwcs/sortstar.c
 * March 14, 2000
 * By Doug Mink
 */

/* void FluxSortStars()		Sort star list based on brightness
 * int StarFluxSort()		Return brightest of two stars based on flux
 * void MagSortStars()		Sort stars list based on magnitude
 * int StarMagSort()		Return brightest of two stars based on mag.
 * void RASortStars()		Sort stars based on right ascension
 * int StarRASort()		Return star with lowest right ascension
 * void XSortStars()		Sort stars based on X coordinate in image
 * int StarXSort()		Return star with lowest X coordinate
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
    double b;		/* First magnitude */
    double r;		/* Second magnitude */
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
    static int StarFluxSort ();
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

static int
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

void
MagSortStars (sn, sra, sdec, spra, spdec, sx, sy, sm, sr, sc, sobj, ns)

double *sn;		/* Identifying number */
double *sra;		/* Right Ascension */
double *sdec;		/* Declination */
double *spra;		/* Right Ascension proper motion */
double *spdec;		/* Declination proper motion */
double *sx;		/* Image X coordinate */
double *sy;		/* Image Y coordinate */
double *sm;		/* First magnitude */
double *sr;		/* Second magnitude */
int    *sc;		/* Other 4-byte information */
char   **sobj;		/* Object name */
int	ns;		/* Number of stars to sort */

{
    StarInfo *stars;
    int i, hasnum, hasmagr, hasobj, haspm;
    static int StarMagSort ();

    stars = (StarInfo *) calloc ((unsigned int)ns, sizeof(StarInfo));

    if (sn == NULL)
	hasnum = 0;
    else
	hasnum = 1;
    if (sr == NULL)
	hasmagr = 0;
    else
	hasmagr = 1;
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
	stars[i].ra = sra[i];
	stars[i].dec = sdec[i];
	if (haspm) {
	    stars[i].pra = spra[i];
	    stars[i].pdec = spdec[i];
	    }
	stars[i].x = sx[i];
	stars[i].y = sy[i];
	stars[i].b = sm[i];
	if (hasmagr)
	    stars[i].r = sr[i];
	stars[i].c = sc[i];
	if (hasobj)
	    stars[i].obj = sobj[i];
	}

    qsort ((char *)stars, ns, sizeof(StarInfo), StarMagSort);

    for (i = 0; i < ns; i++) {
	if (hasnum)
	    sn[i] = stars[i].n;
	sra[i] = stars[i].ra;
	sdec[i] = stars[i].dec;
	if (haspm) {
	    spra[i] = stars[i].pra;
	    spdec[i] = stars[i].pdec;
	    }
	sx[i] = stars[i].x;
	sy[i] = stars[i].y;
	sm[i] = stars[i].b;
	if (hasmagr)
	    sr[i] = stars[i].r;
	sc[i] = stars[i].c;
	if (hasobj)
	    sobj[i] = stars[i].obj;
	}

    free ((char *)stars);
    return;
}


/* StarMagSort -- Order stars in decreasing flux called by qsort */

static int
StarMagSort (ssp1, ssp2)

void *ssp1, *ssp2;

{
    double b1 = ((StarInfo *)ssp1)->b;
    double b2 = ((StarInfo *)ssp2)->b;

    if (b2 < b1)
	return (1);
    else if (b2 > b1)
	return (-1);
    else
	return (0);
}


/* Sort image stars by increasing right ascension */

void
RASortStars (sn, sra, sdec, spra, spdec, sx, sy, sm, sm1, sc, sobj, ns)

double *sn;		/* Identifying number */
double *sra;		/* Right Ascension */
double *sdec;		/* Declination */
double *spra;		/* Right Ascension proper motion */
double *spdec;		/* Declination proper motion */
double *sx;		/* Image X coordinate */
double *sy;		/* Image Y coordinate */
double *sm;		/* First magnitude */
double *sm1;		/* Second magnitude */
int    *sc;		/* Other 4-byte information */
char   **sobj;		/* Object name */
int	ns;		/* Number of stars to sort */
{
    StarInfo *stars;
    int i, hasnum, hasmag1, hasobj, haspm;
    static int StarRASort ();

    stars = (StarInfo *) calloc ((unsigned int)ns, sizeof(StarInfo));

    if (sn == NULL)
	hasnum = 0;
    else
	hasnum = 1;
    if (spra != NULL && spdec != NULL)
	haspm = 1;
    else
	haspm = 0;
    if (sm1 == NULL)
	hasmag1 = 0;
    else
	hasmag1 = 1;
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
	stars[i].x = sx[i];
	stars[i].y = sy[i];
	stars[i].b = sm[i];
	if (hasmag1)
	    stars[i].r = sm1[i];
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
	sx[i] = stars[i].x;
	sy[i] = stars[i].y;
	sm[i] = stars[i].b;
	if (hasmag1)
	    sm1[i] = stars[i].r;
	sc[i] = stars[i].c;
	if (hasobj)
	    sobj[i] = stars[i].obj;
	}

    free ((char *)stars);
    return;
}


/* Order stars in increasing right ascension (called by qsort) */

static int
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


/* XSortStars -- Sort image stars by increasing X value */

void
XSortStars (sn, sra, sdec, spra, spdec, sx, sy, sm, sm1, sc, sobj, ns)

double *sn;		/* Identifying number */
double *sra;		/* Right Ascension */
double *sdec;		/* Declination */
double *spra;		/* Right Ascension proper motion */
double *spdec;		/* Declination proper motion */
double *sx;		/* Image X coordinate */
double *sy;		/* Image Y coordinate */
double *sm;		/* First magnitude */
double *sm1;		/* Second magnitude */
int    *sc;		/* Other 4-byte information */
char   **sobj;		/* Object name */
int	ns;		/* Number of stars to sort */
{
    StarInfo *stars;
    int i, hasnum, hasmag1, hasobj, haspm;
    static int StarXSort ();

    stars = (StarInfo *) calloc ((unsigned int)ns, sizeof(StarInfo));
    if (sn == NULL)
	hasnum = 0;
    else
	hasnum = 1;
    if (spra != NULL && spdec != NULL)
	haspm = 1;
    else
	haspm = 0;
    if (sm1 == NULL)
	hasmag1 = 0;
    else
	hasmag1 = 1;
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
	stars[i].x = sx[i];
	stars[i].y = sy[i];
	stars[i].b = sm[i];
	if (hasmag1)
	    stars[i].r = sm1[i];
	stars[i].c = sc[i];
	if (hasobj)
	    stars[i].obj = sobj[i];
	}

    qsort ((char *)stars, ns, sizeof(StarInfo), StarXSort);

    for (i = 0; i < ns; i++) {
	if (hasnum)
	    sn[i] = stars[i].n;
	sra[i] = stars[i].ra;
	sdec[i] = stars[i].dec;
	if (haspm) {
	    spra[i] = stars[i].pra;
	    spdec[i] = stars[i].pdec;
	    }
	sx[i] = stars[i].x;
	sy[i] = stars[i].y;
	sm[i] = stars[i].b;
	if (hasmag1)
	    sm1[i] = stars[i].r;
	sc[i] = stars[i].c;
	if (hasobj)
	    sobj[i] = stars[i].obj;
	}

    free ((char *)stars);
    return;
}


/* StarXSort -- Order stars in decreasing X value called by qsort */

static int
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
 */
