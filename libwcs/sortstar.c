/* File libwcs/sortstar.c
 * June 13, 1996
 * By Doug Mink
 */

#include <stdlib.h>

/* structure for star lists needed for sorting */
typedef struct {
    double n;
    double ra;
    double dec;
    double b;
    double x;
    double y;
    int    c;
} StarInfo;

static int StarFluxSort ();
static int StarMagSort ();
static int StarRASort ();

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
MagSortStars (sn, sra, sdec, sx, sy, sb, sc, ns)

double	*sn;
double	*sra;
double	*sdec;
double	*sx;
double	*sy;
double	*sb;
int	*sc;
int	ns;

{
    StarInfo *stars;
    int i;

    stars = (StarInfo *) calloc ((unsigned int)ns, sizeof(StarInfo));

    for (i = 0; i < ns; i++) {
	stars[i].n = sn[i];
	stars[i].ra = sra[i];
	stars[i].dec = sdec[i];
	stars[i].x = sx[i];
	stars[i].y = sy[i];
	stars[i].b = sb[i];
	stars[i].c = sc[i];
	}

    qsort ((char *)stars, ns, sizeof(StarInfo), StarMagSort);

    for (i = 0; i < ns; i++) {
	sn[i] = stars[i].n;
	sra[i] = stars[i].ra;
	sdec[i] = stars[i].dec;
	sx[i] = stars[i].x;
	sy[i] = stars[i].y;
	sb[i] = stars[i].b;
	sc[i] = stars[i].c;
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
RASortStars (sn, sra, sdec, sx, sy, sb, sc, ns)

double	*sn;
double	*sra;
double	*sdec;
double	*sx;
double	*sy;
double	*sb;
int	*sc;
int	ns;

{
    StarInfo *stars;
    int i;

    stars = (StarInfo *) calloc ((unsigned int)ns, sizeof(StarInfo));

    for (i = 0; i < ns; i++) {
	stars[i].n = sn[i];
	stars[i].ra = sra[i];
	stars[i].dec = sdec[i];
	stars[i].x = sx[i];
	stars[i].y = sy[i];
	stars[i].b = sb[i];
	stars[i].c = sc[i];
	}

    qsort ((char *)stars, ns, sizeof(StarInfo), StarRASort);

    for (i = 0; i < ns; i++) {
	sn[i] = stars[i].n;
	sra[i] = stars[i].ra;
	sdec[i] = stars[i].dec;
	sx[i] = stars[i].x;
	sy[i] = stars[i].y;
	sb[i] = stars[i].b;
	sc[i] = stars[i].c;
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
