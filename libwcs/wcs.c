/*** File saoimage/wcslib/wcs.c
 *** April 2, 1996
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics

 * Module:	wcs.c (World Coordinate Systems)
 * Purpose:	Convert FITS WCS to pixels and vice versa:
 * Subroutine:	wcsinit (hstring) sets a WCS structure from an image header
 * Subroutine:	wcsset (cra,cdec,secpix,nxpix,nypix,rotate,equinox)
		sets a WCS structure from arguments
 * Subroutine:	iswcs(wcs) returns 1 if WCS structure is filled, else 0
 * Subroutine:	wcscent (wcs) prints the image center and size in WCS units
 * Subroutine:	wcsdist (x1,y1,x2,y2) compute angular distance between ra/dec or lat/long
 * Subroutine:	wcscominit (wcs,command) sets up a command format for execution by wcscom
 * Subroutine:	wcsoutinit (wcs,coor) sets up the output coordinate system
 * Subroutine:	wcscom (wcs,event) executes a command using the current world coordinates
 * Subroutine:	pix2wcs (wcs,xpix,ypix,xpos,ypos) pixel coordinates -> sky coordinates
 * Subroutine:	pix2wcst (wcs,xpix,ypix,wcstring,lstr) pixels -> sky coordinate string
 * Subroutine:	wcs2pix (wcs,xpos,ypos,xpix,ypix) sky coordinates -> pixel coordinates

 * Copyright:   1996 Smithsonian Astrophysical Observatory
 *              You may do anything you like with this file except remove
 *              this copyright.  The Smithsonian Astrophysical Observatory
 *              makes no representations about the suitability of this
 *              software for any purpose.  It is provided "as is" without
 *              express or implied warranty.

 */

#include <string.h>		/* strstr, NULL */
#include <stdio.h>		/* stderr */
#include <math.h>		/* stderr */
#include "wcs.h"
#include "fitshead.h"
#ifndef VMS
#include <stdlib.h>
#endif

/* set up a WCS structure from a FITS or IRAF image header */

struct WorldCoor *wcsinit (hstring)

char *hstring;	/* character string containing FITS header information
		   in the format <keyword>= <value> {/ <comment>} */
{
	struct WorldCoor *wcs;
	char wcstemp[16];
	char *hcoeff;		/* pointer to first coeff's in header */
	char decsign;
	double rah,ram,ras, dsign,decd,decm,decs;
	double dec_deg,ra_hours, secpix;
	int ieq, i;
	double cond2r = 1.745329252e-2;
	char ctypes[8][5];

	strcpy (ctypes[0],"-SIN");
	strcpy (ctypes[1],"-TAN");
	strcpy (ctypes[2],"-ARC");
	strcpy (ctypes[3],"-NCP");
	strcpy (ctypes[4],"-GLS");
	strcpy (ctypes[5],"-MER");
	strcpy (ctypes[6],"-AIT");
	strcpy (ctypes[7],"-STG");

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
	    wcs->plate_ra = ra_hours * 15.0 * cond2r;
	    decsign = '+';
	    hgets (hcoeff,"PLTDECSN", 1, decsign);
	    if (decsign == '-')
		dsign = -1.;
	    else
		dsign = 1.;
	    hgetr8 (hcoeff,"PLTDECD",&decd);
	    hgetr8 (hcoeff,"PLTDECM",&decm);
	    hgetr8 (hcoeff,"PLTDECS",&decs);
	    dec_deg = dsign * (decd+(decm/(double)60.0)+(decs/(double)3600.0));
	    wcs->plate_dec = dec_deg * cond2r;
	    hgetr8 (hcoeff,"EQUINOX",&wcs->equinox);
	    hgeti4 (hcoeff,"EQUINOX",&ieq);
	    if (ieq == 1950)
		strcpy (wcs->radecsys,"FK4");
	    else
		strcpy (wcs->radecsys,"FK5");
	    (void)sprintf (wcs->center,"%2.0f:%2.0f:%5.3f %3.0f:%2.0f:%5.3f %s",
		    rah,ram,ras,dsign*decd,decm,decs,wcs->radecsys);
	    hgetr8 (hcoeff,"PLTSCALE",&wcs->plate_scale);
	    hgetr8 (hcoeff,"CNPIX1",&wcs->x_pixel_offset);
	    hgetr8 (hcoeff,"CNPIX2",&wcs->y_pixel_offset);
	    hcoeff = ksearch (hstring,"XPIXELSZ");
	    hgetr8 (hcoeff,"XPIXELSZ",&wcs->x_pixel_size);
	    hgetr8 (hcoeff,"YPIXELSZ",&wcs->y_pixel_size);
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
	    }

	/* World coordinate system reference coordinate information */
	else if (ksearch (hstring,"CRPIX1") != NULL) {
	    hgetr8 (hstring,"CRPIX1",&wcs->xrefpix);
	    hgetr8 (hstring,"CRPIX2",&wcs->yrefpix);
	    hgetr8 (hstring,"CRVAL1",&wcs->xref);
	    hgetr8 (hstring,"CRVAL2",&wcs->yref);
	    if (hgetr8 (hstring,"CDELT1",&wcs->xinc) != 0) {
		wcs->yinc = wcs->xinc;
		hgetr8 (hstring,"CDELT2",&wcs->yinc);
		wcs->rot = 0.;
		hgetr8 (hstring,"CROTA1",&wcs->rot);
		wcs->cd11 = 0.;
		wcs->cd21 = 0.;
		wcs->cd12 = 0.;
		wcs->cd22 = 0.;
		wcs->crot = cos (wcs->rot*cond2r);
		wcs->srot = sin (wcs->rot*cond2r);
		wcs->rotmat = 0;
		}
	    else if (hgetr8 (hstring,"CD1_1",&wcs->cd11) != 0) {
		wcs->cd12 = 0.;
		hgetr8 (hstring,"CD1_2",&wcs->cd12);
		wcs->cd21 = 0.;
		hgetr8 (hstring,"CD2_1",&wcs->cd21);
		wcs->cd22 = wcs->cd11;
		hgetr8 (hstring,"CD2_2",&wcs->cd22);
		wcs->xinc = 0.;
		wcs->yinc = 0.;
		wcs->rot = atan2 (wcs->cd12, wcs->cd22);
		wcs->crot = cos (wcs->rot);
		wcs->srot = sin (wcs->rot);
		wcs->xinc = wcs->cd11 / wcs->crot;
		wcs->yinc = wcs->cd22 / wcs->crot;
		wcs->rotmat = 1;
		}
	    else {
		(void)fprintf (stderr,"WCSINIT no CD or CDELT, so no WCS\n");
		free (wcs);
		return (NULL);
		}

	/* Coordinate reference frame and equinox */
	    ieq = 1950;
	    if (hgetr8 (hstring,"EQUINOX",&wcs->equinox)) {
		hgeti4 (hstring,"EQUINOX",&ieq);
		}
	    else if (hgetr8 (hstring,"EPOCH",&wcs->equinox)) {
		hgeti4 (hstring,"EPOCH",&ieq);
	        }
	    if (ieq == 0) {
		ieq = 1950;
		wcs->equinox = 1950.0;
		}
	    if (hgets (hstring,"RADECSYS", 16, wcstemp))
		strcpy (wcs->radecsys,wcstemp);
	    else {
		if (ieq > 1980)
		    strcpy (wcs->radecsys,"FK5");
		else
		    strcpy (wcs->radecsys,"FK4");
		}

	/* First coordinate type and projection */
	    if (!hgets (hstring,"CTYPE1", 16, wcstemp)) {
		(void)fprintf (stderr,"WCSINIT CTYPE1 is missing, no WCS\n");
		free (wcs);
		return (NULL);
		}
	    wcs->c1type[0] = wcstemp[0];
	    wcs->c1type[1] = wcstemp[1];
	    if (wcstemp[2] == '-')
		wcs->c1type[2] = 0;
	    else
		wcs->c1type[2] = wcstemp[2];
	    if (wcstemp[3] == '-')
		wcs->c1type[3] = 0;
	    else
		wcs->c1type[3] = wcstemp[3];
	    wcs->c1type[4] = 0;
	    wcs->ptype[0] = wcstemp[4];
	    wcs->ptype[1] = wcstemp[5];
	    wcs->ptype[2] = wcstemp[6];
	    wcs->ptype[3] = wcstemp[7];
	    wcs->ptype[4] = 0;

	/*  Find projection type  */
	    wcs->pcode = 0;  /* default type is linear */
	    for (i=0; i<8; i++)
		if (!strncmp(wcs->ptype, ctypes[i], 4))
		    wcs->pcode = i + 1;

	/* Second coordinate type */
	    if (!hgets (hstring,"CTYPE1", 16, wcstemp)) {
		(void)fprintf (stderr,"WCSINIT CTYPE2 is missing, no WCS\n");
		free (wcs);
		return (NULL);
		}
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

	    if (!strncmp (wcs->c1type, "DEC", 3))
		wcs->coorflip = 1;
	    else
		wcs->coorflip = 0;

	    wcs->wcson = 1;
	    }

	/* Approximate world coordinate system if plate scale is known */
	else if (ksearch (hstring,"SECPIX") != NULL) {
	    hgetr8 (hstring,"SECPIX",&secpix);
	    wcs->yinc = secpix / 3600.0;
	    wcs->xinc = -wcs->yinc;
	    wcs->xrefpix = wcs->nxpix * 0.5;
	    wcs->yrefpix = wcs->nypix * 0.5;

	    wcs->xref = 0.0;
	    if (!hgetra (hstring,"RA",&wcs->xref)) {
		(void)fprintf (stderr,"WCSINIT: No RA with SECPIX, no WCS\n");
		free (wcs);
		return (NULL);
		}
	    wcs->yref = 0.0;
	    if (!hgetdec (hstring,"DEC",&wcs->yref)) {
		(void)fprintf (stderr,"WCSINIT No DEC with SECPIX, no WCS\n");
		free (wcs);
		return (NULL);
		}
	    strcpy (wcs->c1type,"RA");
	    strcpy (wcs->c2type,"DEC");
	    strcpy (wcs->ptype,"-TAN");
	    wcs->pcode = 1;
	    wcs->coorflip = 0;
	    wcs->rot = 0.;
	    hgetr8 (hstring,"CROTA1",&wcs->rot);
	    wcs->cd11 = 0.;
	    wcs->cd21 = 0.;
	    wcs->cd12 = 0.;
	    wcs->cd22 = 0.;
	    wcs->crot = cos (wcs->rot*cond2r);
	    wcs->srot = sin (wcs->rot*cond2r);
	    wcs->rotmat = 0;

	/* Coordinate reference frame and equinox */
	    ieq = 1950;
	    if (hgetr8 (hstring,"EQUINOX",&wcs->equinox))
		hgeti4 (hstring,"EQUINOX",&ieq);
	    else if (hgetr8 (hstring,"EPOCH",&wcs->equinox))
		hgeti4 (hstring,"EPOCH",&ieq);
	    if (ieq == 0) {
		ieq = 1950;
		wcs->equinox = 1950.0;
		}
	    if (!hgets (hstring,"RADECSYS", 16, wcstemp)) {
		if (ieq > 1980)
		    strcpy (wcs->radecsys,"FK5");
		else
		    strcpy (wcs->radecsys,"FK4");
		}
	    else
		strcpy (wcs->radecsys,wcstemp);
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
	strcpy (wcs->search_format,"rgsc_%s");
	return (wcs);
}


/* set up a WCS structure */

struct WorldCoor *wcsset (cra,cdec,secpix,nxpix,nypix,rotate,equinox)

double	cra;	/* Center right ascension in degrees */
double	cdec;	/* Center declination in degrees */
double	secpix;	/* Number of arcseconds per pixel */
int	nxpix;	/* Number of pixels along x-axis */
int	nypix;	/* Number of pixels along y-axis */
double	rotate;	/* Rotation angle (clockwise positive) in degrees */
int	equinox; /* Equinox of coordinates, 1950 and 2000 supported */

{
	struct WorldCoor *wcs;
	double cond2r = 1.745329252e-2;

	wcs = (struct WorldCoor *) calloc (1, sizeof(struct WorldCoor));

	/* Plate solution coefficients */
	wcs->plate_fit = 0;
	wcs->nxpix = nxpix;
	wcs->nypix = nypix;

/* Approximate world coordinate system from a known plate scale */
	wcs->yinc = secpix / 3600.0;
	wcs->xinc = -wcs->yinc;
	wcs->xrefpix = wcs->nxpix * 0.5;
	wcs->yrefpix = wcs->nypix * 0.5;

	wcs->xref = cra;
	wcs->yref = cdec;
	strcpy (wcs->c1type,"RA");
	strcpy (wcs->c2type,"DEC");
	strcpy (wcs->ptype,"-TAN");
	wcs->pcode = 1;
	wcs->coorflip = 0;
	wcs->rot = rotate;
	wcs->cd11 = 0.;
	wcs->cd21 = 0.;
	wcs->cd12 = 0.;
	wcs->cd22 = 0.;
	wcs->crot = cos (wcs->rot*cond2r);
	wcs->srot = sin (wcs->rot*cond2r);
	wcs->rotmat = 0;

	/* Coordinate reference frame and equinox */
	wcs->equinox =  (double) equinox;
	if (equinox > 1980)
	    strcpy (wcs->radecsys,"FK5");
	else
	    strcpy (wcs->radecsys,"FK4");
	wcs->wcson = 1;

	strcpy (wcs->sysout,wcs->radecsys);
	wcs->changesys = 0;
	wcs->printsys = 1;
	wcs->tabsys = 0;
	strcpy (wcs->search_format,"rgsc_%s");
	return (wcs);
}


/* Return 1 if WCS structure is filled, else 0 */

int iswcs (wcs)

struct WorldCoor *wcs;		/* World coordinate system structure */

{
	if (wcs == NULL)
	    return (0);
	else
	    return (wcs->wcson);
}

/* Print position of WCS center, if WCS is set */

void wcscent (wcs)

struct WorldCoor *wcs;		/* World coordinate system structure */

{
	double	xpix,ypix, xpos1, xpos2, ypos1, ypos2;
	char wcstring[32];
	double width, height;
	int lstr = 32;
	double wcsdist();
	void pix2wcs();

	if (wcs == NULL)
	    (void)fprintf (stderr,"No WCS info available\n");
	else {
	    if (wcs->plate_fit)
		(void)fprintf (stderr,"WCS plate center  %s\n", wcs->center);
	    xpix = 0.5 * wcs->nxpix;
	    ypix = 0.5 * wcs->nypix;
	    (void) pix2wcst (wcs,xpix,ypix,wcstring, lstr);
	    (void)fprintf (stderr,"WCS center %s %s %s %s at pixel (%.2f,%.2f)\n",
		     wcs->c1type,wcs->c2type,wcstring,wcs->ptype,xpix,ypix);

	/* Compute image width in degrees */
	    (void) pix2wcs (wcs,1.0,ypix,&xpos1,&ypos1);
	    (void) pix2wcs (wcs,wcs->nxpix,ypix,&xpos2,&ypos2);
	    width = wcsdist (xpos1,ypos1,xpos2,ypos2);
	    if (width < 1/60.0)
		(void) fprintf (stderr, "WCS width = %.2f arcsec ",width*3600.0);
	    else if (width < 1.0)
		(void) fprintf (stderr, "WCS width = %.2f arcmin ",width*60.0);
	    else
		(void) fprintf (stderr, "WCS width = %.3f degrees ",width);

	/* Compute image height in degrees */
	    (void) pix2wcs (wcs,xpix,1.0,&xpos1,&ypos1);
	    (void) pix2wcs (wcs,xpix,wcs->nypix,&xpos2,&ypos2);
	    height = wcsdist (xpos1,ypos1,xpos2,ypos2);
	    if (height < 1/60.0)
	       (void) fprintf (stderr, " height = %.2f arcsec\n",height*3600.0);
	    else if (height < 1.0)
	       (void) fprintf (stderr, " height = %.2f arcmin\n",height*60.0);
	    else
	       (void) fprintf (stderr, " height = %.3f degrees\n",height);
	    }
	return;
}

/* Return RA and Dec of image center, plus size in RA and Dec */

void wcssize (wcs, cra, cdec, dra, ddec)

struct WorldCoor *wcs;		/* World coordinate system structure */
double	*cra;			/* Right ascension of image center (deg) (returned) */
double	*cdec;			/* Declination of image center (deg) (returned) */
double	*dra;
double	*ddec;

{
	double	xpix,ypix, xpos1, xpos2, ypos1, ypos2;
	double	xcent, ycent;
	char wcstring[32];
	double width, height;
	int lstr = 32;
	double wcsdist();
	void pix2wcs();

	if (wcs == NULL)
	    (void)fprintf (stderr,"No WCS info available\n");
	else {
	    if (wcs->plate_fit)
		(void)fprintf (stderr,"WCS plate center  %s\n", wcs->center);
	    xpix = 0.5 * wcs->nxpix;
	    ypix = 0.5 * wcs->nypix;
	    (void) pix2wcs (wcs,xpix,ypix,&xcent, &ycent);
	    *cra = xcent;
	    *cdec = ycent;

	/* Compute image width in degrees */
	    (void) pix2wcs (wcs,1.0,ypix,&xpos1,&ypos1);
	    (void) pix2wcs (wcs,wcs->nxpix,ypix,&xpos2,&ypos2);
	    width = wcsdist (xpos1,ypos1,xpos2,ypos2);
	    *dra = width * 0.5;

	/* Compute image height in degrees */
	    (void) pix2wcs (wcs,xpix,1.0,&xpos1,&ypos1);
	    (void) pix2wcs (wcs,xpix,wcs->nypix,&xpos2,&ypos2);
	    height = wcsdist (xpos1,ypos1,xpos2,ypos2);
	    *ddec = height * 0.5;
	    }
	return;
}


/* Compute distance in degrees between two sky coordinates */

double wcsdist (x1,y1,x2,y2)

double	x1,y1;	/* (RA,Dec) or (Long,Lat) in degrees */
double	x2,y2;	/* (RA,Dec) or (Long,Lat) in degrees */

{
	double xr1, xr2, yr1, yr2;
	double pos1[3], pos2[3], w, diff, cosb;
	double cond2r = 1.745329252e-2;
	int i;

	/* Convert two vectors to direction cosines */
	xr1 = x1 * cond2r;
	yr1 = y1 * cond2r;
	cosb = cos (yr1);
	pos1[0] = cos (xr1) * cosb;
	pos1[1] = sin (xr1) * cosb;
	pos1[2] = sin (yr1);

	xr2 = x2 * cond2r;
	yr2 = y2 * cond2r;
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
	diff = diff / cond2r;
	return (diff);
}


/* Initialize catalog search command set by -wcscom */

void wcscominit (wcs, command)

struct WorldCoor *wcs;		/* World coordinate system structure */
char *command;		/* command with %s where coordinates will go */

{
	int lcom,icom;

	if (!iswcs(wcs))
	    return;
	lcom = strlen (command);
	for (icom = 0; icom < lcom; icom++) {
	    if (command[icom] == '_')
		wcs->search_format[icom] = ' ';
	    else
		wcs->search_format[icom] = command[icom];
	    }
	wcs->search_format[lcom] = 0;
}


/* Execute catalog search command set by -wcscom */

void wcscom ( wcs, xfile, yfile )

struct WorldCoor *wcs;		/* World coordinate system structure */
double xfile,yfile;		/* Image pixel coordinates for WCS command */
{
	char wcstring[32];
	int lstr = 32;
	char command[120];
	int ier;

	if (!iswcs(wcs)) {
	    (void)fprintf(stderr,"WCSCOM: no WCS\n");
	    return;
	    }
	if (wcs->search_format[0] > 0) {

	    /* Get WCS coordinates for this image coordinate */
	    (void) pix2wcst (wcs,xfile,yfile,wcstring,lstr);

	    /* Create and execute search command */
	    (void)sprintf(command, wcs->search_format, wcstring);
	    ier = system (command);
	    if (ier)
		(void)fprintf(stderr,"WCSCOM: %s failed %d\n",command,ier);
	    }
}

/* Initialize WCS output coordinate system set by -wcsout */

void wcsoutinit (wcs, coorsys)

struct WorldCoor *wcs;		/* World coordinate system structure */
char *coorsys;

{
	if (!iswcs(wcs))
	    return;

	if (strcmp (coorsys,"fk4") == 0 || strcmp (coorsys,"FK4") == 0 ||
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
	    }
	return;
}


/* Convert pixel coordinates to World Coordinate string */

int pix2wcst (wcs, xpix, ypix, wcstring, lstr)

struct	WorldCoor *wcs;	/* World coordinate system structure */
double	xpix,ypix;	/* Image coordinates in pixels */
char	*wcstring;	/* World coordinate string (returned) */
int	lstr;		/* Length of world coordinate string (returned) */
{
	double	xpos,ypos, xp, yp;
	double	xpos_deg, ypos_deg;
	int	rah,ram,decd,decm;
	double	ras,decs;
	char	decp;
	char	rastr[16], decstr[16];
	void	pix2wcs();

	if (!iswcs(wcs)) {
	    if (lstr > 0)
		wcstring[0] = 0;
	    return(0);
	    }

	pix2wcs (wcs,xpix,ypix,&xpos,&ypos);

	/* Save degrees -- we might need to use them in the output string */
	xpos_deg = xpos;
	ypos_deg = ypos;

	/* If point is off scale, set string accordingly */
	if (wcs->offscl) {
	    (void)sprintf (wcstring,"Off map");
	    }

	/* Output galactic coordinates in degrees */
	else if (!strncmp (wcs->sysout,"GAL",3)) {
	    if (xpos > 180.0)
		xpos = xpos - 360.0;
	    if (lstr > 19)
		(void)sprintf (wcstring,"%9.5f %9.5f", xpos,ypos);
	    else
		strncpy (wcstring,"*******************",lstr);
	    if (lstr > 28 && wcs->printsys)
		strcat (wcstring," galactic");
	    }
	else if (!strncmp (wcs->sysout,"FK",2)) {
	    ra2str (rastr, xpos, 3);
	    dec2str (decstr, ypos, 2);
	    if (ypos < 0) {
		ypos = -ypos;
		decp = '-';
		}
	    else {
		decp = '+';
		}
	    decd = (int) ypos;
	    yp = (double) 60.0 * (ypos - (double) decd);
	    decm = (int) yp;
	    decs = (double) 60.0 * (yp - (double) decm);
	    if (lstr > 26) {
		if (wcs->tabsys)
		    (void)sprintf (wcstring,"%s	%s", rastr, decstr);
		else
		    (void)sprintf (wcstring,"%s %s", rastr, decstr);
		if (lstr > 31 && wcs->printsys) {
		    if (!strncmp (wcs->sysout,"FK5",3))
			strcat (wcstring," J2000");
		    else if (!strncmp (wcs->sysout,"FK4",3))
			strcat (wcstring," B1950");
		    }
		}
	    else if (lstr >19)
		sprintf (wcstring,"%9.5f %9.5f", xpos_deg, ypos_deg);
	    else
		strncpy (wcstring,"**************************",lstr);
	    }
	else {
	    if (xpos > 180.0)
		xpos = xpos - 360.0;
	    if (lstr > 19)
		(void)sprintf (wcstring,"%9.5f %9.5f", xpos,ypos);
	    else
		strncpy (wcstring,"*******************",lstr);
	    }
	return (1);
}


/* Convert pixel coordinates to World Coordinates */

void pix2wcs (wcs,xpix,ypix,xpos,ypos)

struct WorldCoor *wcs;		/* World coordinate system structure */
double	xpix,ypix;	/* x and y image coordinates in pixels */
double	*xpos,*ypos;	/* RA and Dec in degrees (returned) */
{
	double	xp,yp;

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
	    fk425 (&xp, &yp);
	else if (wcs->changesys == 2)
	    fk524 (&xp, &yp);
	else if (wcs->changesys == 3)
	    fk42gal (&xp, &yp);
	else if (wcs->changesys == 4)
	    fk52gal (&xp, &yp);

	if (!wcs->offscl) {
	    wcs->xpos = xp;
	    wcs->ypos = yp;
	    *xpos = xp;
	    *ypos = yp;
	    }
	return;
}


/* Convert World Coordinates to pixel coordinates */

void wcs2pix (wcs,xpos,ypos,xpix,ypix,offscl)

struct WorldCoor *wcs;		/* World coordinate system structure */
double	xpos,ypos;	/* World coordinates in degrees */
double	*xpix,*ypix;	/* Image coordinates in pixels */

int	*offscl;
{
double	xp,yp;

	if (!iswcs(wcs))
	    return;

	/* Convert coordinates to same system as image */
	if (wcs->changesys == 1)
	    fk524 (&xp, &yp);
	else if (wcs->changesys == 2)
	    fk425 (&xp, &yp);

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
	    if (platepix (xpos, ypos, wcs, &xp, &yp)) {
		wcs->offscl = 1;
		}
	    }
	else if (worldpix (xpos,ypos,wcs,&xp,&yp)) {
	    wcs->offscl = 1;
	    }

	if (wcs->offscl)
	    *offscl = 1;
	else {
	    wcs->xpos = xpos;
	    wcs->ypos = ypos;
	    wcs->offscl = 0;
	    wcs->xpix = xp;
	    wcs->ypix = yp;
	    *xpix = xp;
	    *ypix = yp;
	    wcs->offscl = 0;
	    *offscl = 0;
	    }
	return;
}
/* Oct 28 1994	new program
 * Dec 21 1994	Implement CD rotation matrix
 * Dec 22 1994	Allow RA and DEC to be either x,y or y,x

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

 * Jan 12 1996	Use plane-tangent, not linear, projection if SECPIX is set
 * Jan 12 1996  Add WCSSET to set WCS without an image
 * Feb 15 1996	Replace all calls to HGETC with HGETS
 * Feb 20 1996	Add tab table output from PIX2WCST
 * Apr  2 1996	Convert all equinoxes to B1950 or J2000
 */
