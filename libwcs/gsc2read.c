/*** File libwcs/gsc2read.c
 *** April 24, 2003
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 *** Copyright (C) 2001-2003
 *** Smithsonian Astrophysical Observatory, Cambridge, MA, USA

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.
    
    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    Correspondence concerning WCSTools should be addressed as follows:
           Internet email: dmink@cfa.harvard.edu
           Postal address: Doug Mink
                           Smithsonian Astrophysical Observatory
                           60 Garden St.
                           Cambridge, MA 02138 USA
 */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "fitsfile.h"
#include "wcs.h"
#include "wcscat.h"

#define ABS(a) ((a) < 0 ? (-(a)) : (a))
#define LINE    1024

/* GSC II search engine URL */
char gsc2url[64]="http://www-gsss.stsci.edu/cgi-bin/gsc22query.exe";

/* GSC2READ -- Read GSC II catalog stars over the web */

int
gsc2read (cra,cdec,dra,ddec,drad,distsort,sysout,eqout,epout,mag1,mag2,
	 sortmag,nstarmax,gnum,gra,gdec,gmag,gtype,nlog)

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
int	sortmag;	/* Magnitude by which to sort (1 to nmag) */
int	nstarmax;	/* Maximum number of stars to be returned */
double	*gnum;		/* Array of Guide Star numbers (returned) */
double	*gra;		/* Array of right ascensions (returned) */
double	*gdec;		/* Array of declinations (returned) */
double	**gmag;		/* 2-D array of magnitudes (returned) */
int	*gtype;		/* Array (returned) */
int	nlog;		/* 1 for diagnostics */
{
    double *gpra, *gpdec;
    char srchurl[LINE];
    char temp[64];
    struct TabTable *tabtable;
    double dtemp;
    struct StarCat *starcat;
    int nstar;
    double ra, dec, mag;
    char rastr[32], decstr[32];

    if (nstarmax < 1)
	nlog = -1;

    gpra = NULL;
    gpdec = NULL;

/* make mag1 always the smallest magnitude */
    if (mag2 < mag1) {
	mag = mag2;
	mag2 = mag1;
	mag1 = mag;
	}

    /* Set up query for STScI GSC II server */
    ra = cra;
    dec = cdec;
    if (sysout != WCS_J2000)
	wcscon (sysout, WCS_J2000, eqout, 2000.0, &ra, &dec, epout);
    ra2str (rastr, 32, ra, 3);
    dec2str (decstr, 32, dec, 2);
    sprintf (srchurl, "?ra=%s&dec=%s&", rastr, decstr);
    if (drad != 0.0) {
	dtemp = drad * 60.0;
	sprintf (temp, "r1=0&r2=%.3f&",dtemp);
	}
    else {
	if (dra > ddec)
	    sprintf (temp, "r1=0&r2=%.3f&",dra*60.0);
	else
	    sprintf (temp, "r1=0&r2=%.3f&",ddec*60.0);
	}
    strcat (srchurl, temp);
    sprintf (temp, "m1=%.2f&m2=99.9&", mag1);
    strcat (srchurl, temp);
    nstar = 100000;
    sprintf (temp, "n=%d",nstar);
    strcat (srchurl, temp);
    if (nlog > 0)
	fprintf (stderr,"%s%s\n", gsc2url, srchurl);

    /* Run search across the web */
    if ((tabtable = webopen (gsc2url, srchurl, nlog)) == NULL) {
	if (nlog > 0)
	    fprintf (stderr, "WEBREAD: %s failed\n", srchurl);
	return (0);
	}

    /* Return if no data */
    if (tabtable->tabdata == NULL || strlen (tabtable->tabdata) == 0 ||
	!strncasecmp (tabtable->tabdata, "[EOD]", 5)) {
	if (nlog > 0)
	    fprintf (stderr, "WEBRNUM: No data returned\n");
	return (0);
	}

    /* Dump returned file and stop */
    if (nlog < 0) {
	fwrite  (tabtable->tabbuff, tabtable->lbuff, 1, stdout);
	exit (0);
	}

    /* Open returned Starbase table as a catalog */
    if ((starcat = tabcatopen (gsc2url, tabtable)) == NULL) {
	if (nlog > 0)
	    fprintf (stderr, "WEBREAD: Could not open Starbase table as catalog\n");
	return (0);
	}

    /* Set reference frame, epoch, and equinox of catalog */
    starcat->rpmunit = PM_MASYR;
    starcat->dpmunit = PM_MASYR;
    starcat->coorsys = WCS_J2000;
    starcat->epoch = 2000.0;
    starcat->equinox = 2000.0;
    starcat->nmag = 5;

    /* Extract desired sources from catalog  and return them */
    nstar = tabread (gsc2url,distsort,cra,cdec,dra,ddec,drad,
	     sysout,eqout,epout,mag1,mag2,sortmag,nstarmax,&starcat,
	     gnum,gra,gdec,gpra,gpdec,gmag,gtype,NULL,nlog);

    tabcatclose (starcat);
    starcat = NULL;

    return (nstar);
}

/* Jun 22 2001	New program
 * Jun 28 2001	Set proper motion to milliarcseconds/year
 * Jun 29 2001	Always set maximum magnitude to 99.9 to get Tycho-2 stars, too
 * Sep 13 2001	Pass array of magnitudes, not vector
 * Sep 14 2001	Add option to print entire returned file if nlog < 0
 * Sep 20 2001	Make argument starcat, not *starcat in tabcatclose()
 *
 * Apr  8 2002	Fix bugs in null subroutine gsc2rnum()
 * Oct  3 2002	If nstarmax is less than 1, print everything returned
 *
 * Feb  6 2003	Reset nmag to 4 because there is an epoch column
 * Mar 11 2003	Fix URL for search
 * Apr  3 2003	Drop unused variables after lint; drop gsc2rnum()
 * Apr 24 2003	Set nmag to 5 to include epoch, which is not printed
 */
