/*** File libwcs/wcs.c
 *** October 9, 1997
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics

 * Module:	wcs.c (World Coordinate Systems)
 * Purpose:	Convert FITS WCS to pixels and vice versa:
 * Subroutine:	wcsinit (hstring) sets a WCS structure from an image header
 * Subroutine:	wcsninit (hstring,lh) sets a WCS structure from an image header
 * Subroutine:	wcsset (cra,cdec,secpix,xrpix,yrpix,nxpix,nypix,rotate,equinox,epoch,proj)
 *		sets a WCS structure from arguments
 * Subroutine:	iswcs(wcs) returns 1 if WCS structure is filled, else 0
 * Subroutine:	nowcs(wcs) returns 0 if WCS structure is filled, else 1
 * Subroutine:	wcscent (wcs) prints the image center and size in WCS units
 * Subroutine:	wcssize (wcs, cra, cdec, dra, ddec) returns image center and size
 * Subroutine:	wcsfull (wcs, cra, cdec, width, height) returns image center and size
 * Subroutine:	wcsshift (wcs,cra,cdec) resets the center of a WCS structure
 * Subroutine:	wcsdist (x1,y1,x2,y2) compute angular distance between ra/dec or lat/long
 * Subroutine:	wcscominit (wcs,command) sets up a command format for execution by wcscom
 * Subroutine:	wcsoutinit (wcs,coor) sets up the output coordinate system
 * Subroutine:	wcsout(wcs) returns current output coordinate system
 * Subroutine:	wcscom (wcs,file,x,y) executes a command using the current world coordinates
 * Subroutine:	pix2wcs (wcs,xpix,ypix,xpos,ypos) pixel coordinates -> sky coordinates
 * Subroutine:	pix2wcst (wcs,xpix,ypix,wcstring,lstr) pixels -> sky coordinate string
 * Subroutine:	wcs2pix (wcs,xpos,ypos,xpix,ypix,offscl) sky coordinates -> pixel coordinates

 * Copyright:   1996 Smithsonian Astrophysical Observatory
 *              You may do anything you like with this file except remove
 *              this copyright.  The Smithsonian Astrophysical Observatory
 *              makes no representations about the suitability of this
 *              software for any purpose.  It is provided "as is" without
 *              express or implied warranty.

 */

#include <string.h>		/* strstr, NULL */
#include <stdio.h>		/* stderr */
#include <math.h>
#include "wcs.h"
#include "fitshead.h"
#ifndef VMS
#include <stdlib.h>
#endif

static void wcseq();
static char wcserrmsg[80];

/* set up a WCS structure from a FITS image header lhstring bytes long */

struct WorldCoor *
wcsninit (hstring, lhstring)

char	*hstring;	/* character string containing FITS header information
		   	in the format <keyword>= <value> {/ <comment>} */
int	lhstring;	/* Length of FITS header in bytes */
{
    hlength (hstring, lhstring);
    return (wcsinit (hstring));
}

/* set up a WCS structure from a FITS image header */

struct WorldCoor *
wcsinit (hstring)

char *hstring;	/* character string containing FITS header information
		   in the format <keyword>= <value> {/ <comment>} */
{
	struct WorldCoor *wcs;
	char wcstemp[16];
	char *hcoeff;		/* pointer to first coeff's in header */
	char decsign;
	double rah,ram,ras, dsign,decd,decm,decs;
	double dec_deg,ra_hours, secpix, cddet;
	int ieq, i;
	char ctypes[8][5];
	char *str;

#ifdef USE_SAOLIB
	set_saolib((void *)hstring);
#endif

	strcpy (ctypes[0],"SIN");
	strcpy (ctypes[1],"TAN");
	strcpy (ctypes[2],"ARC");
	strcpy (ctypes[3],"NCP");
	strcpy (ctypes[4],"GLS");
	strcpy (ctypes[5],"MER");
	strcpy (ctypes[6],"AIT");
	strcpy (ctypes[7],"STG");

	wcs = (struct WorldCoor *) calloc (1, sizeof(struct WorldCoor));

	/* Plate solution coefficients */
	wcs->plate_fit = 0;
	hgetr8 (hstring,"NAXIS1",&wcs->nxpix);
	hgetr8 (hstring,"NAXIS2",&wcs->nypix);
	if (ksearch (hstring,"PLTRAH") != NULL) {
	    wcs->plate_fit = 1;
	    hcoeff = ksearch (hstring,"PLTRAH");
	    hgetr8 (hcoeff,"PLTRAH",&rah);
	    hgetr8 (hcoeff,"PLTRAM",&ram);
	    hgetr8 (hcoeff,"PLTRAS",&ras);
	    ra_hours = rah + (ram / (double)60.0) + (ras / (double)3600.0);
	    wcs->plate_ra = hrrad (ra_hours);
	    decsign = '+';
	    hgets (hcoeff,"PLTDECSN", 1, &decsign);
	    if (decsign == '-')
		dsign = -1.;
	    else
		dsign = 1.;
	    hgetr8 (hcoeff,"PLTDECD",&decd);
	    hgetr8 (hcoeff,"PLTDECM",&decm);
	    hgetr8 (hcoeff,"PLTDECS",&decs);
	    dec_deg = dsign * (decd+(decm/(double)60.0)+(decs/(double)3600.0));
	    wcs->plate_dec = degrad (dec_deg);
	    hgetr8 (hstring,"EQUINOX",&wcs->equinox);
	    hgeti4 (hstring,"EQUINOX",&ieq);
	    if (ieq == 1950)
		strcpy (wcs->radecsys,"FK4");
	    else
		strcpy (wcs->radecsys,"FK5");
	    wcs->epoch = wcs->equinox;
	    hgetr8 (hstring,"EPOCH",&wcs->epoch);
	    (void)sprintf (wcs->center,"%2.0f:%2.0f:%5.3f %c%2.0f:%2.0f:%5.3f %s",
		    rah,ram,ras,decsign,decd,decm,decs,wcs->radecsys);
	    hgetr8 (hstring,"PLTSCALE",&wcs->plate_scale);
	    hgetr8 (hstring,"XPIXELSZ",&wcs->x_pixel_size);
	    hgetr8 (hstring,"YPIXELSZ",&wcs->y_pixel_size);
	    hgetr8 (hstring,"CNPIX1",&wcs->x_pixel_offset);
	    hgetr8 (hstring,"CNPIX2",&wcs->y_pixel_offset);
	    hcoeff = ksearch (hstring,"PPO1");
	    hgetr8 (hcoeff,"PPO1",&wcs->ppo_coeff[0]);
	    hgetr8 (hcoeff,"PPO2",&wcs->ppo_coeff[1]);
	    hgetr8 (hcoeff,"PPO3",&wcs->ppo_coeff[2]);
	    hgetr8 (hcoeff,"PPO4",&wcs->ppo_coeff[3]);
	    hgetr8 (hcoeff,"PPO5",&wcs->ppo_coeff[4]);
	    hgetr8 (hcoeff,"PPO6",&wcs->ppo_coeff[5]);
	    hcoeff = ksearch (hstring,"AMDX1");
	    hgetr8 (hcoeff,"AMDX1",&wcs->amd_x_coeff[0]);
	    hgetr8 (hcoeff,"AMDX2",&wcs->amd_x_coeff[1]);
	    hgetr8 (hcoeff,"AMDX3",&wcs->amd_x_coeff[2]);
	    hgetr8 (hcoeff,"AMDX4",&wcs->amd_x_coeff[3]);
	    hgetr8 (hcoeff,"AMDX5",&wcs->amd_x_coeff[4]);
	    hgetr8 (hcoeff,"AMDX6",&wcs->amd_x_coeff[5]);
	    hgetr8 (hcoeff,"AMDX7",&wcs->amd_x_coeff[6]);
	    hgetr8 (hcoeff,"AMDX8",&wcs->amd_x_coeff[7]);
	    hgetr8 (hcoeff,"AMDX9",&wcs->amd_x_coeff[8]);
	    hgetr8 (hcoeff,"AMDX10",&wcs->amd_x_coeff[9]);
	    hgetr8 (hcoeff,"AMDX11",&wcs->amd_x_coeff[10]);
	    hgetr8 (hcoeff,"AMDX12",&wcs->amd_x_coeff[11]);
	    hgetr8 (hcoeff,"AMDX13",&wcs->amd_x_coeff[12]);
	    hgetr8 (hcoeff,"AMDX14",&wcs->amd_x_coeff[13]);
	    hgetr8 (hcoeff,"AMDX15",&wcs->amd_x_coeff[14]);
	    hgetr8 (hcoeff,"AMDX16",&wcs->amd_x_coeff[15]);
	    hgetr8 (hcoeff,"AMDX17",&wcs->amd_x_coeff[16]);
	    hgetr8 (hcoeff,"AMDX18",&wcs->amd_x_coeff[17]);
	    hgetr8 (hcoeff,"AMDX19",&wcs->amd_x_coeff[18]);
	    hgetr8 (hcoeff,"AMDX20",&wcs->amd_x_coeff[19]);
	    hcoeff = ksearch (hstring,"AMDY1");
	    hgetr8 (hcoeff,"AMDY1",&wcs->amd_y_coeff[0]);
	    hgetr8 (hcoeff,"AMDY2",&wcs->amd_y_coeff[1]);
	    hgetr8 (hcoeff,"AMDY3",&wcs->amd_y_coeff[2]);
	    hgetr8 (hcoeff,"AMDY4",&wcs->amd_y_coeff[3]);
	    hgetr8 (hcoeff,"AMDY5",&wcs->amd_y_coeff[4]);
	    hgetr8 (hcoeff,"AMDY6",&wcs->amd_y_coeff[5]);
	    hgetr8 (hcoeff,"AMDY7",&wcs->amd_y_coeff[6]);
	    hgetr8 (hcoeff,"AMDY8",&wcs->amd_y_coeff[7]);
	    hgetr8 (hcoeff,"AMDY9",&wcs->amd_y_coeff[8]);
	    hgetr8 (hcoeff,"AMDY10",&wcs->amd_y_coeff[9]);
	    hgetr8 (hcoeff,"AMDY11",&wcs->amd_y_coeff[10]);
	    hgetr8 (hcoeff,"AMDY12",&wcs->amd_y_coeff[11]);
	    hgetr8 (hcoeff,"AMDY13",&wcs->amd_y_coeff[12]);
	    hgetr8 (hcoeff,"AMDY14",&wcs->amd_y_coeff[13]);
	    hgetr8 (hcoeff,"AMDY15",&wcs->amd_y_coeff[14]);
	    hgetr8 (hcoeff,"AMDY16",&wcs->amd_y_coeff[15]);
	    hgetr8 (hcoeff,"AMDY17",&wcs->amd_y_coeff[16]);
	    hgetr8 (hcoeff,"AMDY18",&wcs->amd_y_coeff[17]);
	    hgetr8 (hcoeff,"AMDY19",&wcs->amd_y_coeff[18]);
	    hgetr8 (hcoeff,"AMDY20",&wcs->amd_y_coeff[19]);
	    wcs->wcson = 1;
	    (void)strcpy (wcs->c1type, "RA");
	    (void)strcpy (wcs->c2type, "DEC");
	    (void)strcpy (wcs->ptype, "PLATE");
	    wcs->degout = 0;
	    wcs->ndec = 3;
	    }

	/* World coordinate system reference coordinate information */
	else if (hgets (hstring,"CTYPE1", 16, wcstemp)) {

	/* Deal appropriately with linear coordinates */
	    if (!strncmp (wcstemp,"LINEAR",6)) {
		wcs->pcode = 0;
		strcpy (wcs->c1type, wcstemp);
		strcpy (wcs->ptype, wcstemp);
		}

	/* Deal appropriately with pixel coordinates */
	    else if (!strncmp (wcstemp,"PIXEL",6)) {
		wcs->pcode = -1;
		strcpy (wcs->c1type, wcstemp);
		strcpy (wcs->ptype, wcstemp);
		}

	/* Set up right ascension, declination, latitude, or longitude */
	    else if (wcstemp[0] == 'R' ||
		     wcstemp[0] == 'D' ||
		     wcstemp[0] == 'A' ||
		     wcstemp[1] == 'L') {

		strncpy (wcs->c1type, wcstemp,4);
		if (wcs->c1type[2] == '-')
		    wcs->c1type[2] = 0;
		if (wcs->c1type[3] == '-')
		    wcs->c1type[3] = 0;
		wcs->c1type[4] = 0;

		if (wcstemp[4] == '-')
		    strncpy (wcs->ptype, wcstemp+5,3);
		else
		    strncpy (wcs->ptype, wcstemp+4,4);
		wcs->ptype[4] = 0;

	    /*  Find projection type  */
		celset (wcs->ptype, wcs->cel, wcs->prj);
		wcs->pcode = 0;  /* default type is linear */
		for (i = 0; i < 8; i++)
		    if (!strncmp(wcs->ptype, ctypes[i], 3))
			wcs->pcode = i + 1;
		}

	/* If not linear or sky coordinates, drop out with error message */
	    else {
		(void)sprintf (wcserrmsg,"WCSINIT: CTYPE1 not sky coordinates or LINEAR -> no WCS\n");
		free (wcs);
		return (NULL);
		}

	/* Second coordinate type */
	    if (!hgets (hstring,"CTYPE2", 16, wcstemp)) {
		(void)sprintf (wcserrmsg,"WCSINIT: No CTYPE2 -> no WCS\n");
		free (wcs);
		return (NULL);
		}

	/* Deal appropriately with linear coordinates */
	    if (!strncmp (wcstemp,"LINEAR",6)) {
		wcs->pcode = 0;
		strcpy (wcs->c2type, wcstemp);
		}

	/* Deal appropriately with pixel coordinates */
	    else if (!strncmp (wcstemp,"PIXEL",6)) {
		wcs->pcode = -1;
		strcpy (wcs->c2type, wcstemp);
		}

	/* Set up right ascension, declination, latitude, or longitude */
	    else if (wcstemp[0] == 'R' ||
		     wcstemp[0] == 'D' ||
		     wcstemp[0] == 'A' ||
		     wcstemp[1] == 'L') {
		wcs->c2type[0] = wcstemp[0];
		wcs->c2type[1] = wcstemp[1];
		if (wcstemp[2] == '-')
		    wcs->c2type[2] = 0;
		else
		    wcs->c2type[2] = wcstemp[2];
		if (wcstemp[3] == '-')
		    wcs->c2type[3] = 0;
		else
		    wcs->c2type[3] = wcstemp[3];
		wcs->c2type[4] = 0;

		if (!strncmp (wcs->c1type, "DEC", 3) ||
		    !strncmp (wcs->c1type, "GLAT", 4))
		    wcs->coorflip = 1;
		else
		    wcs->coorflip = 0;
		if (wcstemp[1] == 'L' || wcstemp[0] == 'A') {
		    wcs->degout = 1;
		    wcs->ndec = 5;
		    }
		else {
		    wcs->degout = 0;
		    wcs->ndec = 3;
		    }
		}

	/* If not linear or sky coordinates, drop out with error message */
	    else {
		(void)sprintf (wcserrmsg,"WCSINIT: CTYPE2 not sky coordinates or LINEAR -> no WCS\n");
		free (wcs);
		return (NULL);
		}

	/* Reference pixel coordinates and WCS value */
	    wcs->xrefpix = 1.0;
	    wcs->yrefpix = 1.0;
	    wcs->xref = 1.0;
	    wcs->yref = 1.0;
	    hgetr8 (hstring,"CRPIX1",&wcs->xrefpix);
	    hgetr8 (hstring,"CRPIX2",&wcs->yrefpix);
	    hgetr8 (hstring,"CRVAL1",&wcs->xref);
	    hgetr8 (hstring,"CRVAL2",&wcs->yref);
	    if (hgetr8 (hstring,"CDELT1",&wcs->xinc) != 0) {
		wcs->yinc = wcs->xinc;
		hgetr8 (hstring,"CDELT2",&wcs->yinc);
		wcs->rot = 0.;
		hgetr8 (hstring,"CROTA1",&wcs->rot);
		if (wcs->rot == 0.)
		    hgetr8 (hstring,"CROTA2",&wcs->rot);
		wcs->cd11 = 0.;
		wcs->cd21 = 0.;
		wcs->cd12 = 0.;
		wcs->cd22 = 0.;
		wcs->crot = cos (degrad (wcs->rot));
		wcs->srot = sin (degrad (wcs->rot));
		wcs->rotmat = 0;
		}
	    else if (hgetr8 (hstring,"CD1_1",&wcs->cd11) != 0) {
		wcs->rotmat = 1;
		wcs->cd12 = 0.;
		hgetr8 (hstring,"CD1_2",&wcs->cd12);
		wcs->cd21 = 0.;
		hgetr8 (hstring,"CD2_1",&wcs->cd21);
		wcs->cd22 = wcs->cd11;
		hgetr8 (hstring,"CD2_2",&wcs->cd22);
		cddet = (wcs->cd11 * wcs->cd22) - (wcs->cd12 * wcs->cd21);
		if (cddet != 0.0) {
		    wcs->dc11 = wcs->cd22 / cddet;
		    wcs->dc12 = -wcs->cd21 / cddet;
		    wcs->dc21 = -wcs->cd12 / cddet;
		    wcs->dc22 = wcs->cd11 / cddet;
		    }
		wcs->xinc = 0.;
		wcs->yinc = 0.;
		wcs->rot = 0.;
		wcs->crot = 1.;
		wcs->srot = 0.;
		}
	    else {
		wcs->xinc = 1.0;
		wcs->yinc = 1.0;
		(void)sprintf (wcserrmsg,"WCSINIT: setting CDELT to 1\n");
		}

	/* Coordinate reference frame, equinox, and epoch */
	    if (strncmp (wcs->ptype,"LINEAR",6) &&
		strncmp (wcs->ptype,"PIXEL",5))
		wcseq (hstring,wcs);
	    else
		wcs->degout = -1;

	    wcs->wcson = 1;
	    }

	/* Approximate world coordinate system if plate scale is known */
	else if (ksearch (hstring,"SECPIX") != NULL ||
		 ksearch (hstring,"PIXSCAL1") != NULL ||
		 ksearch (hstring,"SECPIX1") != NULL) {
	    secpix = 0.0;
	    hgetr8 (hstring,"SECPIX",&secpix);
	    if (secpix == 0.0) {
		hgetr8 (hstring,"SECPIX1",&secpix);
		if (secpix != 0.0) {
		    wcs->xinc = -secpix / 3600.0;
		    hgetr8 (hstring,"SECPIX2",&secpix);
		    wcs->yinc = secpix / 3600.0;
		    }
		else {
		    hgetr8 (hstring,"PIXSCAL1",&secpix);
		    wcs->xinc = -secpix / 3600.0;
		    hgetr8 (hstring,"PIXSCAL2",&secpix);
		    wcs->yinc = secpix / 3600.0;
		    }
		}
	    else {
		wcs->yinc = secpix / 3600.0;
		wcs->xinc = -wcs->yinc;
		}
	    wcs->xrefpix = wcs->nxpix * 0.5;
	    wcs->yrefpix = wcs->nypix * 0.5;
	    hgetr8 (hstring,"CRPIX1",&wcs->xrefpix);
	    hgetr8 (hstring,"CRPIX2",&wcs->yrefpix);

	    wcs->xref = 0.0;
	    if (!hgetra (hstring,"RA",&wcs->xref)) {
		(void)sprintf (wcserrmsg,"WCSINIT: No RA with SECPIX, no WCS\n");
		free (wcs);
		return (NULL);
		}
	    wcs->yref = 0.0;
	    if (!hgetdec (hstring,"DEC",&wcs->yref)) {
		(void)sprintf (wcserrmsg,"WCSINIT No DEC with SECPIX, no WCS\n");
		free (wcs);
		return (NULL);
		}
	    strcpy (wcs->c1type,"RA");
	    strcpy (wcs->c2type,"DEC");
	    strcpy (wcs->ptype,"-TAN");
	    wcs->pcode = 1;
	    wcs->coorflip = 0;
	    wcs->rot = 0.;
	    wcs->degout = 0;
	    wcs->ndec = 3;
	    hgetr8 (hstring,"CROTA1",&wcs->rot);
	    if (wcs->rot == 0.)
	    hgetr8 (hstring,"CROTA2",&wcs->rot);
	    wcs->cd11 = 0.;
	    wcs->cd21 = 0.;
	    wcs->cd12 = 0.;
	    wcs->cd22 = 0.;
	    wcs->dc11 = 0.;
	    wcs->dc21 = 0.;
	    wcs->dc12 = 0.;
	    wcs->dc22 = 0.;
	    wcs->crot = cos (degrad (wcs->rot));
	    wcs->srot = sin (degrad (wcs->rot));
	    wcs->rotmat = 0;

	/* Coordinate reference frame and equinox */
	    wcseq (hstring,wcs);

	/* Epoch of image (from observation date, if possible) */
	    if (!hgetdate (hstring,"DATE-OBS",&wcs->epoch)) {
		if (!hgetr8 (hstring,"EPOCH",&wcs->epoch)) {
		    wcs->epoch = wcs->equinox;
		    }
		}
	    wcs->wcson = 1;
	    }

	else {
	    free (wcs);
	    return (NULL);
	    }

	strcpy (wcs->sysout,wcs->radecsys);
	wcs->changesys = 0;
	wcs->printsys = 1;
	wcs->tabsys = 0;
	if ((str = getenv("WCS_COMMAND")) != NULL ) {
	    int icom;
	    int lcom = strlen (str);
	    for (icom = 0; icom < lcom; icom++) {
		if (str[icom] == '_')
		    str[icom] = ' ';
		}
	    strcpy (wcs->search_format, str);
	    }
	else
	    strcpy (wcs->search_format, "rgsc %s");
	return (wcs);
}


static void
wcseq (hstring, wcs)

char	*hstring;
struct WorldCoor *wcs;
{
    int ieq = 0;
    char wcstemp[16];

    /* Set equinox from EQUINOX, EPOCH, or RADECSYS; default to 2000 */
    if (hgeti4 (hstring,"EQUINOX",&ieq))
	hgetr8 (hstring,"EQUINOX",&wcs->equinox);

    else if (hgeti4 (hstring,"EPOCH",&ieq))
	if (ieq == 0) {
	    ieq = 1950;
	    wcs->equinox = 1950.0;
	    }
	else
            hgetr8 (hstring,"EPOCH",&wcs->equinox);

    else if (hgets (hstring,"RADECSYS", 16, wcstemp)) {
	if (!strncmp (wcstemp,"FK4",3)) {
	    wcs->equinox = 1950.0;
	    ieq = 1950;
	    }
	else if (!strncmp (wcstemp,"FK5",3)) {
	    wcs->equinox = 2000.0;
	    ieq = 2000;
	    }
	else if (!strncmp (wcstemp,"GAL",3)) {
	    wcs->equinox = 2000.0;
	    ieq = 2000;
	    }
	else if (!strncmp (wcstemp,"ECL",3)) {
	    wcs->equinox = 2000.0;
	    ieq = 2000;
	    }
	}

    if (ieq == 0) {
	wcs->equinox = 2000.0;
	ieq = 2000;
	}

    /* Epoch of image (from observation date, if possible) */
    if (!hgetdate (hstring,"DATE-OBS",&wcs->epoch)) {
	if (!hgetr8 (hstring,"EPOCH",&wcs->epoch)) {
	    wcs->epoch = wcs->equinox;
	    }
	}
    if (wcs->epoch == 0.0)
	wcs->epoch = wcs->equinox;

    /* Set coordinate system from keyword, if it is present */
    if (hgets (hstring,"RADECSYS", 16, wcstemp)) {
	strcpy (wcs->radecsys,wcstemp);
	if (!strncmp (wcs->radecsys,"FK4",3))
	    wcs->equinox = 1950.0;
	else if (!strncmp (wcs->radecsys,"FK5",3))
	    wcs->equinox = 2000.0;
	else if (!strncmp (wcs->radecsys,"GAL",3) && ieq == 0)
	    wcs->equinox = 2000.0;
	}

    /* Set galactic coordinates if GLON or GLAT are in C1TYPE */
    else if (wcs->c1type[0] == 'G')
	strcpy (wcs->radecsys,"GALACTIC");
    else if (wcs->c1type[0] == 'E')
	strcpy (wcs->radecsys,"ECLIPTIC");
    else if (wcs->c1type[0] == 'S')
	strcpy (wcs->radecsys,"SGALACTC");
    else if (wcs->c1type[0] == 'H')
	strcpy (wcs->radecsys,"HELIOECL");
    else if (wcs->c1type[0] == 'A')
	strcpy (wcs->radecsys,"ALTAZ");

    /* Otherwise set coordinate system from equinox */
    /* Systemless coordinates cannot be translated using b, j, or g commands */
    else {
	if (ieq > 1980)
	    strcpy (wcs->radecsys,"FK5");
	else
	    strcpy (wcs->radecsys,"FK4");
	}

    return;
}


/* set up a WCS structure */

struct WorldCoor *
wcsset (cra,cdec,secpix,xrpix,yrpix,nxpix,nypix,rotate,equinox,epoch,proj)

double	cra;	/* Center right ascension in degrees */
double	cdec;	/* Center declination in degrees */
double	secpix;	/* Number of arcseconds per pixel */
double	xrpix;	/* Reference pixel X coordinate */
double	yrpix;	/* Reference pixel X coordinate */
int	nxpix;	/* Number of pixels along x-axis */
int	nypix;	/* Number of pixels along y-axis */
double	rotate;	/* Rotation angle (clockwise positive) in degrees */
int	equinox; /* Equinox of coordinates, 1950 and 2000 supported */
double	epoch;	/* Epoch of coordinates, used for FK4/FK5 conversion
		 * no effect if 0 */
char	*proj;	/* Projection */

{
	struct WorldCoor *wcs;
	char *str;

	wcs = (struct WorldCoor *) calloc (1, sizeof(struct WorldCoor));

	/* Plate solution coefficients */
	wcs->plate_fit = 0;
	wcs->nxpix = nxpix;
	wcs->nypix = nypix;

/* Approximate world coordinate system from a known plate scale */
	wcs->yinc = secpix / 3600.0;
	wcs->xinc = -wcs->yinc;
	wcs->xrefpix = xrpix;
	wcs->yrefpix = yrpix;

	wcs->xref = cra;
	wcs->yref = cdec;
	strcpy (wcs->c1type,"RA");
	strcpy (wcs->c2type,"DEC");
	strcpy (wcs->ptype,proj);
	wcs->pcode = 1;
	wcs->coorflip = 0;
	wcs->rot = rotate;
	wcs->crot = cos (degrad (wcs->rot));
	wcs->srot = sin (degrad (wcs->rot));
	wcs->rotmat = 0;
	wcs->cd11 = 0.;
	wcs->cd21 = 0.;
	wcs->cd12 = 0.;
	wcs->cd22 = 0.;
	wcs->dc11 = 0.;
	wcs->dc21 = 0.;
	wcs->dc12 = 0.;
	wcs->dc22 = 0.;

	/* Coordinate reference frame and equinox */
	wcs->equinox =  (double) equinox;
	if (equinox > 1980)
	    strcpy (wcs->radecsys,"FK5");
	else
	    strcpy (wcs->radecsys,"FK4");
	if (epoch > 0)
	    wcs->epoch = epoch;
	else
	    wcs->epoch = 0.0;
	wcs->wcson = 1;

	strcpy (wcs->sysout,wcs->radecsys);
	wcs->changesys = 0;
	wcs->printsys = 1;
	wcs->tabsys = 0;
	if ((str = getenv("WCS_COMMAND")) != NULL ) {
	    int icom;
	    int lcom = strlen (str);
	    for (icom = 0; icom < lcom; icom++) {
		if (str[icom] == '_')
		    str[icom] = ' ';
		}
	    strcpy (wcs->search_format, str);
	    }
	else
	    strcpy (wcs->search_format, "rgsc %s");
	return (wcs);
}


/* Return 1 if WCS structure is filled, else 0 */

int
iswcs (wcs)

struct WorldCoor *wcs;		/* World coordinate system structure */

{
    if (wcs == NULL)
	return (0);
    else
	return (wcs->wcson);
}


/* Return 0 if WCS structure is filled, else 1 */

int
nowcs (wcs)

struct WorldCoor *wcs;		/* World coordinate system structure */

{
    if (wcs == NULL)
	return (1);
    else
	return (!wcs->wcson);
}


/* Reset the center of a WCS structure */

void
wcsshift (wcs,cra,cdec,coorsys)

struct WorldCoor *wcs;	/* World coordinate system structure */
double	cra;		/* New center right ascension in degrees */
double	cdec;		/* New center declination in degrees */
char	*coorsys;	/* FK4 or FK5 coordinates (1950 or 2000) */

{
    if (nowcs (wcs))
	return;

/* Approximate world coordinate system from a known plate scale */
    wcs->xrefpix = wcs->nxpix * 0.5;
    wcs->yrefpix = wcs->nypix * 0.5;
    wcs->xref = cra;
    wcs->yref = cdec;

/* Coordinate reference frame */
    strcpy (wcs->radecsys,coorsys);
    if (strcmp (coorsys,"FK4") == 0)
	wcs->equinox = 1950.0;
    else
	wcs->equinox = 2000.0;

    return;
}

/* Print position of WCS center, if WCS is set */

void
wcscent (wcs)

struct WorldCoor *wcs;		/* World coordinate system structure */

{
	double	xpix,ypix, xpos1, xpos2, ypos1, ypos2;
	char wcstring[32];
	double width, height, secpix;
	int lstr = 32;

	if (wcs == NULL)
	    (void)fprintf (stderr,"No WCS information available\n");
	else {
	    if (wcs->plate_fit)
		(void)fprintf (stderr,"WCS plate center  %s\n", wcs->center);
	    xpix = 0.5 * wcs->nxpix;
	    ypix = 0.5 * wcs->nypix;
	    (void) pix2wcst (wcs,xpix,ypix,wcstring, lstr);
	    (void)fprintf (stderr,"WCS center %s %s %s %s at pixel (%.2f,%.2f)\n",
		     wcs->c1type,wcs->c2type,wcstring,wcs->ptype,xpix,ypix);

	/* Image width */
	    (void) pix2wcs (wcs,1.0,ypix,&xpos1,&ypos1);
	    (void) pix2wcs (wcs,wcs->nxpix,ypix,&xpos2,&ypos2);
	    width = wcsdist (xpos1,ypos1,xpos2,ypos2);
	    if (width < 1/60.0)
		(void) fprintf (stderr, "WCS width = %.2f arcsec ",width*3600.0);
	    else if (width < 1.0)
		(void) fprintf (stderr, "WCS width = %.2f arcmin ",width*60.0);
	    else
		(void) fprintf (stderr, "WCS width = %.3f degrees ",width);

	/* Image height */
	    (void) pix2wcs (wcs,xpix,1.0,&xpos1,&ypos1);
	    (void) pix2wcs (wcs,xpix,wcs->nypix,&xpos2,&ypos2);
	    height = wcsdist (xpos1,ypos1,xpos2,ypos2);
	    if (height < 1/60.0)
	       (void) fprintf (stderr, " height = %.2f arcsec",height*3600.0);
	    else if (height < 1.0)
	       (void) fprintf (stderr, " height = %.2f arcmin",height*60.0);
	    else
	       (void) fprintf (stderr, " height = %.3f degrees",height);

	/* Image scale */
	    secpix = 3600.0 * height / (double)wcs->nypix;
	    if (secpix < 100.0)
	       (void) fprintf (stderr, "  %.3f arcsec/pixel\n",secpix);
	    else if (secpix < 3600.0)
	       (void) fprintf (stderr, "  %.3f arcmin/pixel\n",secpix*60.0);
	    else
	       (void) fprintf (stderr, "  %.3f degrees/pixel\n",secpix*3600.0);
	    }
	return;
}

/* Return RA and Dec of image center, plus size in RA and Dec */

void
wcssize (wcs, cra, cdec, dra, ddec)

struct WorldCoor *wcs;	/* World coordinate system structure */
double	*cra;		/* Right ascension of image center (deg) (returned) */
double	*cdec;		/* Declination of image center (deg) (returned) */
double	*dra;		/* Half-width in right ascension (deg) (returned) */
double	*ddec;		/* Half-width in declination (deg) (returned) */

{
	double	xpix,ypix, xpos1, xpos2, ypos1, ypos2;
	double	xcent, ycent;
	double width, height;

	/* Find right ascension and declination of coordinates */
	if (iswcs(wcs)) {
	    xpix = 0.5 * wcs->nxpix;
	    ypix = 0.5 * wcs->nypix;
	    (void) pix2wcs (wcs,xpix,ypix,&xcent, &ycent);
	    *cra = xcent;
	    *cdec = ycent;

	/* Compute image half-width in degrees of right ascension */
	    (void) pix2wcs (wcs,1.0,ypix,&xpos1,&ypos1);
	    (void) pix2wcs (wcs,wcs->nxpix,ypix,&xpos2,&ypos2);
	    if (strncmp (wcs->ptype,"LINEAR",6) &&
		strncmp (wcs->ptype,"PIXEL",5)) {
		width = wcsdist (xpos1,ypos1,xpos2,ypos2);
		*dra = (width * 0.5) / cos (degrad (*cdec));
		}
	    else
		*dra = sqrt (((ypos2-ypos1) * (ypos2-ypos1)) +
				((xpos2-xpos1) * (xpos2-xpos1)));

	/* Compute image half-height in degrees of declination*/
	    (void) pix2wcs (wcs,xpix,1.0,&xpos1,&ypos1);
	    (void) pix2wcs (wcs,xpix,wcs->nypix,&xpos2,&ypos2);
	    if (strncmp (wcs->ptype,"LINEAR",6) &&
		strncmp (wcs->ptype,"PIXEL",5)) {
		height = wcsdist (xpos1,ypos1,xpos2,ypos2);
		*ddec = height * 0.5;
		}
	    else
		*ddec = sqrt (((ypos2-ypos1) * (ypos2-ypos1)) +
				((xpos2-xpos1) * (xpos2-xpos1)));
	    }
	return;
}


/* Return RA and Dec of image center, plus size in degrees */

void
wcsfull (wcs, cra, cdec, width, height)

struct WorldCoor *wcs;	/* World coordinate system structure */
double	*cra;		/* Right ascension of image center (deg) (returned) */
double	*cdec;		/* Declination of image center (deg) (returned) */
double	*width;		/* Width in degrees (returned) */
double	*height;	/* Height in degrees (returned) */

{
	double	xpix,ypix, xpos1, xpos2, ypos1, ypos2;
	double	xcent, ycent;

	/* Find right ascension and declination of coordinates */
	if (iswcs(wcs)) {
	    xpix = 0.5 * wcs->nxpix;
	    ypix = 0.5 * wcs->nypix;
	    (void) pix2wcs (wcs,xpix,ypix,&xcent, &ycent);
	    *cra = xcent;
	    *cdec = ycent;

	/* Compute image width in degrees */
	    (void) pix2wcs (wcs,1.0,ypix,&xpos1,&ypos1);
	    (void) pix2wcs (wcs,wcs->nxpix,ypix,&xpos2,&ypos2);
	    if (strncmp (wcs->ptype,"LINEAR",6) &&
		strncmp (wcs->ptype,"PIXEL",5))
		*width = wcsdist (xpos1,ypos1,xpos2,ypos2);
	    else
		*width = sqrt (((ypos2-ypos1) * (ypos2-ypos1)) +
				((xpos2-xpos1) * (xpos2-xpos1)));

	/* Compute image height in degrees */
	    (void) pix2wcs (wcs,xpix,1.0,&xpos1,&ypos1);
	    (void) pix2wcs (wcs,xpix,wcs->nypix,&xpos2,&ypos2);
	    if (strncmp (wcs->ptype,"LINEAR",6) &&
		strncmp (wcs->ptype,"PIXEL",5))
		*height = wcsdist (xpos1,ypos1,xpos2,ypos2);
	    else
		*height = sqrt (((ypos2-ypos1) * (ypos2-ypos1)) +
				((xpos2-xpos1) * (xpos2-xpos1)));
	    }
	return;
}


/* Compute distance in degrees between two sky coordinates */

double
wcsdist (x1,y1,x2,y2)

double	x1,y1;	/* (RA,Dec) or (Long,Lat) in degrees */
double	x2,y2;	/* (RA,Dec) or (Long,Lat) in degrees */

{
	double xr1, xr2, yr1, yr2;
	double pos1[3], pos2[3], w, diff, cosb;
	int i;

	/* Convert two vectors to direction cosines */
	xr1 = degrad (x1);
	yr1 = degrad (y1);
	cosb = cos (yr1);
	pos1[0] = cos (xr1) * cosb;
	pos1[1] = sin (xr1) * cosb;
	pos1[2] = sin (yr1);

	xr2 = degrad (x2);
	yr2 = degrad (y2);
	cosb = cos (yr2);
	pos2[0] = cos (xr2) * cosb;
	pos2[1] = sin (xr2) * cosb;
	pos2[2] = sin (yr2);

	/* Modulus squared of half the difference vector */
	w = 0.0;
	for (i = 0; i < 3; i++) {
	    w = w + (pos1[i] - pos2[i]) * (pos1[i] - pos2[i]);
	    }
	w = w / 4.0;
	if (w > 1.0) w = 1.0;

	/* Angle beween the vectors */
	diff = 2.0 * atan2 (sqrt (w), sqrt (1.0 - w));
	diff = raddeg (diff);
	return (diff);
}


/* Initialize catalog search command set by -wcscom */

void
wcscominit (wcs, command)

struct WorldCoor *wcs;		/* World coordinate system structure */
char *command;		/* command with %s where coordinates will go */

{
    int lcom,icom;

    if (iswcs(wcs)) {
	lcom = strlen (command);
	if (lcom > 0) {
	    for (icom = 0; icom < lcom; icom++) {
		if (command[icom] == '_')
		    wcs->search_format[icom] = ' ';
		else
		    wcs->search_format[icom] = command[icom];
		}
	    wcs->search_format[lcom] = 0;
	    }
	}
    return;
}


/* Execute Unix command with world coordinates (from x,y) and/or filename */

void
wcscom ( wcs, filename, xfile, yfile )

struct WorldCoor *wcs;		/* World coordinate system structure */
char	*filename;		/* Image file name */
double	xfile,yfile;		/* Image pixel coordinates for WCS command */
{
    char wcstring[32];
    int lstr = 32;
    char command[120];
    char comform[120];
    char *fileform, *posform;
    int ier;

    if (wcs->search_format[0] > 0)
	strcpy (comform, wcs->search_format);
    else
	strcpy (comform, "rgsc %s");

    if (!iswcs(wcs))
	(void)fprintf(stderr,"WCSCOM: no WCS\n");

    else if (comform[0] > 0) {

	/* Get WCS coordinates for this image coordinate */
	(void) pix2wcst (wcs,xfile,yfile,wcstring,lstr);

	/* Create and execute search command */
	if ((fileform = strsrch (comform,"%f")) != NULL) {
	    posform = strsrch (comform,"%s");
	    *(fileform+1) = 's';
	    if (fileform < posform)
		(void)sprintf(command, comform, filename, wcstring);
	    else
		(void)sprintf(command, comform, wcstring, filename);
	    }
	else
	    (void)sprintf(command, comform, wcstring);
	ier = system (command);
	if (ier)
	    (void)fprintf(stderr,"WCSCOM: %s failed %d\n",command,ier);
	}
    return;
}

/* Initialize WCS output coordinate system set by -wcsout */

void
wcsoutinit (wcs, coorsys)

struct WorldCoor *wcs;		/* World coordinate system structure */
char *coorsys;

{
    if (!iswcs(wcs))
	return;

    if (wcs->radecsys[0] == 0)
	wcs->sysout[0] = 0;
    else if (strcmp (coorsys,"fk4") == 0 || strcmp (coorsys,"FK4") == 0 ||
	strncmp (coorsys,"b1",2) == 0 || strncmp (coorsys,"B1",2) == 0)
	strcpy (wcs->sysout,"FK4");
    else if (strcmp (coorsys,"fk5") == 0 || strcmp (coorsys,"FK5") == 0 ||
	     strncmp (coorsys,"j2",2) == 0 || strncmp (coorsys,"J2",2) == 0)
	strcpy (wcs->sysout,"FK5");
    else if (strncmp (coorsys,"g",1) == 0 || strncmp (coorsys,"G",1) == 0)
	strcpy (wcs->sysout,"GALACTIC");
    else {
	strcpy (wcs->sysout,wcs->radecsys);
	wcs->changesys = 0;
	}
    if (wcs->wcson) {
	if (strncmp (wcs->radecsys,"FK4",3) == 0 &&
	    strncmp(wcs->sysout,"FK5",3) == 0)
	    wcs->changesys = 1;
	else if (strncmp (wcs->radecsys,"FK5",3) == 0 &&
	    strncmp(wcs->sysout,"FK4",3) == 0)
	    wcs->changesys = 2;
	else if (strncmp (wcs->radecsys,"FK4",3) == 0 &&
	    strncmp(wcs->sysout,"GAL",3) == 0)
	    wcs->changesys = 3;
	else if (strncmp (wcs->radecsys,"FK5",3) == 0 &&
	    strncmp(wcs->sysout,"GAL",3) == 0)
	    wcs->changesys = 4;
	else
	    wcs->changesys = 0;
	if (strncmp(wcs->sysout,"GAL",3) == 0) {
	    wcs->degout = 1;
	    wcs->ndec = 5;
	    }
	else if (strncmp(wcs->sysout,"ALT",3) == 0) {
	    wcs->degout = 1;
	    wcs->ndec = 5;
	    }
	else {
	    wcs->degout = 0;
	    wcs->ndec = 3;
	    }
	}
    return;
}

/* Return current value of WCS output coordinate system set by -wcsout */

char *
wcsout(wcs)
struct	WorldCoor *wcs; /* World coordinate system structure */
{
    return(wcs->sysout);
}




/* Convert pixel coordinates to World Coordinate string */

int
pix2wcst (wcs, xpix, ypix, wcstring, lstr)

struct	WorldCoor *wcs;	/* World coordinate system structure */
double	xpix,ypix;	/* Image coordinates in pixels */
char	*wcstring;	/* World coordinate string (returned) */
int	lstr;		/* Length of world coordinate string (returned) */
{
	double	xpos,ypos;
	char	rastr[16], decstr[16];
	int	minlength;

	if (!iswcs(wcs)) {
	    if (lstr > 0)
		wcstring[0] = 0;
	    return(0);
	    }

	pix2wcs (wcs,xpix,ypix,&xpos,&ypos);

	/* Keep ra/longitude within range */
	if ((!strncmp (wcs->sysout,"GAL",3) ||
	    !strncmp (wcs->sysout,"ECL",3)) && xpos > 180.0) 
	    xpos = xpos - 360.0;

	else if (!strncmp (wcs->sysout,"FK",2) && xpos > 360.0)
	    xpos = xpos - 360.0;

	/* If point is off scale, set string accordingly */
	if (wcs->offscl) {
	    (void)sprintf (wcstring,"Off map");
	    return (1);
	    }
	else if (wcs->degout == 1) {
	    minlength = 9 + (2 * wcs->ndec);
	    if (lstr > minlength) {
		deg2str (rastr, xpos, wcs->ndec);
		deg2str (decstr, ypos, wcs->ndec);
		if (wcs->tabsys)
		    (void)sprintf (wcstring,"%s	%s", rastr, decstr);
		else
		    (void)sprintf (wcstring,"%s %s", rastr, decstr);
		lstr = lstr - minlength;
		}
	    else {
		strncpy (wcstring,"*******************",lstr);
		lstr = 0;
		}
	    }

	else if (wcs->degout == 0) {
	    minlength = 18 + (2 * wcs->ndec);
	    if (lstr > minlength) {
		ra2str (rastr, xpos, wcs->ndec);
		dec2str (decstr, ypos, wcs->ndec-1);
		if (wcs->tabsys)
		    (void)sprintf (wcstring,"%s	%s", rastr, decstr);
		else
		    (void)sprintf (wcstring,"%s %s", rastr, decstr);
	        lstr = lstr - minlength;
		}
	    else {
		strncpy (wcstring,"**************************",lstr);
		lstr = 0;
		}
	    }

	/* Label galactic coordinates */
	if (!strncmp (wcs->sysout,"GAL",3)) {
	    if (lstr > 9 && wcs->printsys)
		strcat (wcstring," galactic");
	    }

	/* Label ecliptic coordinates */
	else if (!strncmp (wcs->sysout,"ECL",3)) {
	    if (lstr > 9 && wcs->printsys)
		strcat (wcstring," ecliptic");
	    }

	/* Label alt-az coordinates */
	else if (!strncmp (wcs->sysout,"ALT",3)) {
	    if (lstr > 7 && wcs->printsys)
		strcat (wcstring," alt-az");
	    }

	/* Label equatorial coordinates */
	else if (!strncmp (wcs->sysout,"FK",2)) {
	    if (lstr > 6 && wcs->printsys) {
		if (!strncmp (wcs->sysout,"FK5",3))
		    strcat (wcstring," J2000");
		else if (!strncmp (wcs->sysout,"FK4",3))
		    strcat (wcstring," B1950");
		}
	    }

	/* Output linear coordinates */
	else {
	    if (!strncmp (wcs->ptype, "LINEAR",6) &&
		xpos > 180.0)
		xpos = xpos - 360.0;
	    if (lstr > 23)
		(void)sprintf (wcstring,"%11.5f %11.5f", xpos,ypos);
	    else
		strncpy (wcstring,"*******************",lstr);
	    }
	return (1);
}


/* Convert pixel coordinates to World Coordinates */

void
pix2wcs (wcs,xpix,ypix,xpos,ypos)

struct WorldCoor *wcs;		/* World coordinate system structure */
double	xpix,ypix;	/* x and y image coordinates in pixels */
double	*xpos,*ypos;	/* RA and Dec in degrees (returned) */
{
	double	xp,yp;
    extern int platepos(), worldpos();
    extern void fk4prec(),fk5prec(),fk425e(),fk524e(),fk42gal(),fk52gal();

	if (!iswcs(wcs))
	    return;
	wcs->xpix = xpix;
	wcs->ypix = ypix;
	wcs->offscl = 0;

	/* Convert image coordinates to sky coordinates */
	if (wcs->plate_fit) {
	    if (platepos (xpix, ypix, wcs, &xp, &yp)) {
		wcs->offscl = 1;
		}
	    }
	else if (worldpos (xpix, ypix, wcs, &xp, &yp)) {
	    wcs->offscl = 1;
	    }

	if (wcs->pcode > 0) {

	    /* Convert coordinates to FK4 or FK5 */
	    if (strncmp (wcs->radecsys,"FK4",3) == 0) {
		if (wcs->equinox != 1950.0)
		    fk4prec (wcs->equinox, 1950.0, &xp, &yp);
		}
	    else if (strncmp (wcs->radecsys,"FK5",3) == 0) {
		if (wcs->equinox != 2000.0)
		    fk5prec (wcs->equinox, 2000.0, &xp, &yp);
		}

	    /* Convert coordinates to desired output system */
	    if (wcs->changesys == 1)
		fk425e (&xp, &yp, wcs->epoch);
	    else if (wcs->changesys == 2)
		fk524e (&xp, &yp, wcs->epoch);
	    else if (wcs->changesys == 3)
		fk42gal (&xp, &yp);
	    else if (wcs->changesys == 4)
		fk52gal (&xp, &yp);
	    }

	if (!wcs->offscl) {
	    wcs->xpos = xp;
	    wcs->ypos = yp;
	    *xpos = xp;
	    *ypos = yp;
	    }
	return;
}


/* Convert World Coordinates to pixel coordinates */

void
wcs2pix (wcs,xpos,ypos,xpix,ypix,offscl)

struct WorldCoor *wcs;		/* World coordinate system structure */
double	xpos,ypos;	/* World coordinates in degrees */
double	*xpix,*ypix;	/* Image coordinates in pixels */

int	*offscl;
{
    double xp,yp;
    extern int platepix(), worldpix();
    extern void fk4prec(), fk5prec(), fk425e(), fk524e();

	if (!iswcs(wcs))
	    return;

	*offscl = 0;
	xp = xpos;
	yp = ypos;

	/* Convert coordinates to same system as image */
	if (wcs->changesys == 1)
	    fk524e (&xp, &yp, wcs->epoch);
	else if (wcs->changesys == 2)
	    fk425e (&xp, &yp, wcs->epoch);

	/* Convert coordinates from FK4 or FK5 to equinox used */
	if (strncmp (wcs->radecsys,"FK4",3) == 0) {
	    if (wcs->equinox != 1950.0)
		fk4prec (1950.0, wcs->equinox, &xp, &yp);
	    }
	else if (strncmp (wcs->radecsys,"FK5",3) == 0) {
	    if (wcs->equinox != 2000.0)
		fk5prec (2000.0, wcs->equinox, &xp, &yp);
	    }

	/* Convert sky coordinates to image coordinates */
	if (wcs->plate_fit) {
	    if (platepix (xp, yp, wcs, xpix, ypix)) {
		*offscl = 1;
		}
	    }
	else if (worldpix (xp,yp,wcs,xpix,ypix)) {
	    *offscl = 1;
	    }
	if (*xpix < 0 || *ypix < 0)
	    *offscl = 1;
	else if (*xpix > wcs->nxpix + 1 || *ypix > wcs->nypix + 1)
	    *offscl = 1;

	wcs->xpix = *xpix;
	wcs->ypix = *ypix;
	wcs->offscl = *offscl;
	wcs->xpos = xpos;
	wcs->ypos = ypos;
	return;
}

void
wcserr ()

{
    fprintf (stderr, "%s\n",wcserrmsg);
    return;
}
/* Oct 28 1994	new program
 * Dec 21 1994	Implement CD rotation matrix
 * Dec 22 1994	Allow RA and DEC to be either x,y or y,x
 *
 * Mar  6 1995	Add Digital Sky Survey plate fit
 * May  2 1995	Add prototype of PIX2WCST to WCSCOM
 * May 25 1995	Print leading zero for hours and degrees
 * Jun 21 1995	Add WCS2PIX to get pixels from WCS
 * Jun 21 1995	Read plate scale from FITS header for plate solution
 * Jul  6 1995	Pass WCS structure as argument; malloc it in WCSINIT
 * Jul  6 1995	Check string lengths in PIX2WCST
 * Aug 16 1995	Add galactic coordinate conversion to PIX2WCST
 * Aug 17 1995	Return 0 from iswcs if wcs structure is not yet set
 * Sep  8 1995	Do not include malloc.h if VMS
 * Sep  8 1995	Check for legal WCS before trying anything
 * Sep  8 1995	Do not try to set WCS if missing key keywords
 * Oct 18 1995	Add WCSCENT and WCSDIST to print center and size of image
 * Nov  6 1995	Include stdlib.h instead of malloc.h
 * Dec  6 1995	Fix format statement in PIX2WCST
 * Dec 19 1995	Change MALLOC to CALLOC to initialize array to zeroes
 * Dec 19 1995	Explicitly initialize rotation matrix and yinc
 * Dec 22 1995	If SECPIX is set, use approximate WCS
 * Dec 22 1995	Always print coordinate system
 *
 * Jan 12 1996	Use plane-tangent, not linear, projection if SECPIX is set
 * Jan 12 1996  Add WCSSET to set WCS without an image
 * Feb 15 1996	Replace all calls to HGETC with HGETS
 * Feb 20 1996	Add tab table output from PIX2WCST
 * Apr  2 1996	Convert all equinoxes to B1950 or J2000
 * Apr 26 1996	Get and use image epoch for accurate FK4/FK5 conversions
 * May 16 1996	Clean up internal documentation
 * May 17 1996	Return width in right ascension degrees, not sky degrees
 * May 24 1996	Remove extraneous print command from WCSSIZE
 * May 28 1996	Add NOWCS and WCSSHIFT subroutines
 * Jun 11 1996	Drop unused variables after running lint
 * Jun 12 1996	Set equinox as well as system in WCSSHIFT
 * Jun 14 1996	Make DSS keyword searches more robust
 * Jul  1 1996	Allow for SECPIX1 and SECPIX2 keywords
 * Jul  2 1996	Test for CTYPE1 instead of CRVAL1
 * Jul  5 1996	Declare all subroutines in wcs.h
 * Jul 19 1996	Add subroutine WCSFULL to return real image size
 * Aug 12 1996	Allow systemless coordinates which cannot be converted
 * Aug 15 1996	Allow LINEAR WCS to pass numbers through transparently
 * Aug 15 1996	Add WCSERR to print error message under calling program control
 * Aug 16 1996	Add latitude and longitude as image coordinate types
 * Aug 26 1996	Fix arguments to HLENGTH in WCSNINIT
 * Aug 28 1996	Explicitly set OFFSCL in WCS2PIX if coordinates outside image
 * Sep  3 1996	Return computed pixel values even if they are offscale
 * Sep  6 1996	Allow filename to be passed by WCSCOM
 * Oct  8 1996	Default to 2000 for EQUINOX and EPOCH and FK5 for RADECSYS
 * Oct  8 1996	If EPOCH is 0 and EQUINOX is not set, default to 1950 and FK4
 * Oct 15 1996  Add comparison when testing an assignment
 * Oct 16 1996  Allow PIXEL CTYPE which means WCS is same as image coordinates
 * Oct 21 1996	Add WCS_COMMAND environment variable
 * Oct 25 1996	Add image scale to WCSCENT
 * Oct 30 1996	Fix bugs in WCS2PIX
 * Oct 31 1996	Fix CD matrix rotation angle computation
 * Oct 31 1996	Use inline degree <-> radian conversion functions
 * Nov  1 1996	Add option to change number of decimal places in PIX2WCST
 * Nov  5 1996	Set wcs->crot to 1 if rotation matrix is used
 * Dec  2 1996	Add altitide/azimuth coordinates
 * Dec 13 1996	Fix search format setting from environment
 *
 * Jan 22 1997	Add ifdef for Eric Mandel (SAOtng)
 * Feb  5 1997	Add wcsout for Eric Mandel
 * Mar 20 1997	Drop unused variable STR in WCSCOM
 * May 21 1997	Do not make pixel coordinates mod 360 in PIX2WCST
 * May 22 1997	Add PIXEL pcode = -1;
 * Jul 11 1997	Get center pixel x and y from header even if no WCS
 * Aug  7 1997	Add NOAO PIXSCALi keywords for default WCS
 * Oct  3 1997	Add MROT and XMPIX,YMPIX for chip rotation
 * Oct  9 1997	Implement Mark Calabretta's WCSLIB
 */
