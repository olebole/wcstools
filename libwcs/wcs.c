/*** File libwcs/wcs.c
 *** May 15, 1998
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics

 * Module:	wcs.c (World Coordinate Systems)
 * Purpose:	Convert FITS WCS to pixels and vice versa:
 * Subroutine:	wcsinit (hstring) sets a WCS structure from an image header
 * Subroutine:	wcsninit (hstring,lh) sets a WCS structure from an image header
 * Subroutine:	wcsxinit (cra,cdec,secpix,xrpix,yrpix,nxpix,nypix,rotate,equinox,epoch,proj)
 *		sets a WCS structure from arguments
 * Subroutine:	wcsreset (wcs,crpix1,crpix2,crval1,crval2,cdelt1,cdelt2,crota,cd, equinox)
 *		resets an existing WCS structure from arguments
 * Subroutine:	wcseqset (wcs, equinox) resets an existing WCS structure to new equinox
 * Subroutine:	iswcs(wcs) returns 1 if WCS structure is filled, else 0
 * Subroutine:	nowcs(wcs) returns 0 if WCS structure is filled, else 1
 * Subroutine:	wcscent (wcs) prints the image center and size in WCS units
 * Subroutine:	wcssize (wcs, cra, cdec, dra, ddec) returns image center and size
 * Subroutine:	wcsfull (wcs, cra, cdec, width, height) returns image center and size
 * Subroutine:	wcsshift (wcs,cra,cdec) resets the center of a WCS structure
 * Subroutine:	wcsdist (x1,y1,x2,y2) compute angular distance between ra/dec or lat/long
 * Subroutine:	wcscominit (wcs,command) sets up a command format for execution by wcscom
 * Subroutine:	wcsoutinit (wcs,coor) sets up the coordinate system used by pix2wcs
 * Subroutine:	getwcsout (wcs) returns current output coordinate system used by pix2wcs
 * Subroutine:	wcsininit (wcs,coor) sets up the coordinate system used by wcs2pix
 * Subroutine:	getwcsin (wcs) returns current input coordinate system used by wcs2pix
 * Subroutine:	setdegout(wcs, new) sets WCS output in degrees or hh:mm:ss
 * Subroutine:	getradecsys(wcs) returns current coordinate system type
 * Subroutine:	wcscom (wcs,file,x,y) executes a command using the current world coordinates
 * Subroutine:	setlinmode (wcs, mode) sets output string mode for LINEAR
 * Subroutine:	pix2wcst (wcs,xpix,ypix,wcstring,lstr) pixels -> sky coordinate string
 * Subroutine:	pix2wcs (wcs,xpix,ypix,xpos,ypos) pixel coordinates -> sky coordinates
 * Subroutine:	wcsc2pix (wcs,xpos,ypos,coorsys,xpix,ypix,offscl) sky coordinates -> pixel coordinates
 * Subroutine:	wcs2pix (wcs,xpos,ypos,xpix,ypix,offscl) sky coordinates -> pixel coordinates

 * Copyright:   1998 Smithsonian Astrophysical Observatory
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

static int oldwcs0 = 0;

/* set up a WCS structure from a FITS image header lhstring bytes long */

struct WorldCoor *
wcsninit (hstring, lhstring)

char	*hstring;	/* character string containing FITS header information
		   	in the format <keyword>= <value> [/ <comment>] */
int	lhstring;	/* Length of FITS header in bytes */
{
    hlength (hstring, lhstring);
    return (wcsinit (hstring));
}

/* set up a WCS structure from a FITS image header */

struct WorldCoor *
wcsinit (hstring)

char *hstring;	/* character string containing FITS header information
		   in the format <keyword>= <value> [/ <comment>] */
{
    struct WorldCoor *wcs;
    char wcstemp[16];
    char *hcoeff;		/* pointer to first coeff's in header */
    char decsign;
    double rah,ram,ras, dsign,decd,decm,decs;
    double dec_deg,ra_hours, secpix, ra0, ra1, dec0, dec1, radiff;
    char keyword[16];
    int ieq, i, j, naxes, mem;
    /* int ix1, ix2, iy1, iy2, idx1, idx2, idy1, idy2;
    double dxrefpix, dyrefpix;
    char temp[32];
    char *ic; */
    char *str;
    double *pci;
    double mjd;
    double s, srot, crot, sdelt;
    extern int matinv();
    extern int tnxinit();
    extern int platepos();
    int  nctype = 30;
    char ctypes[30][4];

    strcpy (ctypes[0], "LIN");
    strcpy (ctypes[1], "AZP");
    strcpy (ctypes[2], "TAN");
    strcpy (ctypes[3], "SIN");
    strcpy (ctypes[4], "STG");
    strcpy (ctypes[5], "ARC");
    strcpy (ctypes[6], "ZPN");
    strcpy (ctypes[7], "ZEA");
    strcpy (ctypes[8], "AIR");
    strcpy (ctypes[9], "CYP");
    strcpy (ctypes[10], "CAR");
    strcpy (ctypes[11], "MER");
    strcpy (ctypes[12], "CEA");
    strcpy (ctypes[13], "COP");
    strcpy (ctypes[14], "COD");
    strcpy (ctypes[15], "COE");
    strcpy (ctypes[16], "COO");
    strcpy (ctypes[17], "BON");
    strcpy (ctypes[18], "PCO");
    strcpy (ctypes[19], "GLS");
    strcpy (ctypes[20], "PAR");
    strcpy (ctypes[21], "AIT");
    strcpy (ctypes[22], "MOL");
    strcpy (ctypes[23], "CSC");
    strcpy (ctypes[24], "QSC");
    strcpy (ctypes[25], "TSC");
    strcpy (ctypes[26], "NCP");
    strcpy (ctypes[27], "DSS");
    strcpy (ctypes[28], "PLT");
    strcpy (ctypes[29], "TNX");

    wcs = (struct WorldCoor *) calloc (1, sizeof(struct WorldCoor));

    /* Set WCSLIB flags so that structures will be reinitialized */
    wcs->cel.flag = 0;
    wcs->lin.flag = 0;
    wcs->wcsl.flag = 0;

    /* Initialize to no plate fit */
    wcs->ncoeff1 = 0;
    wcs->ncoeff2 = 0;

    /* Initialize to no CD matrix */
    wcs->cd[0] = 0.0;

    /* Header parameters independent of projection */
    hgeti4 (hstring, "NAXIS", &wcs->lin.naxis);
    hgeti4 (hstring, "NAXIS", &wcs->naxes);
    hgeti4 (hstring, "NAXIS", &naxes);
    hgetr8 (hstring, "NAXIS1", &wcs->nxpix);
    hgetr8 (hstring, "NAXIS2", &wcs->nypix);
    hgets (hstring, "INSTRUME", 16, wcs->instrument);
    hgeti4 (hstring, "DETECTOR", &wcs->detector);
    wcs->oldwcs = oldwcs0;
    for (i = 0; i < 16; i++) wcs->pc[i] = 0.0;
    for (i = 0; i < naxes; i++) wcs->pc[(i*naxes)+i] = 1.0;

    /* World coordinate system reference coordinate information */
    if (hgets (hstring,"CTYPE1", 16, wcstemp)) {
	strcpy (wcs->ctype[0], wcstemp);

	/* Deal appropriately with linear coordinates */
	if (!strncmp (wcstemp,"LINEAR",6)) {
	    wcs->prjcode = 0;
	    strcpy (wcs->c1type, wcstemp);
	    if (!hgets (hstring, "CUNIT1", 16, wcs->units[0])) {
		if (!mgets (hstring, "WAT1", "units", 16, wcs->units[0])) {
		    wcs->units[0][0] = 0;
		    }
		}
	    strcpy (wcs->ptype, wcstemp);
	    }

	/* Deal appropriately with pixel coordinates */
	else if (!strncmp (wcstemp,"PIXEL",6)) {
	    wcs->prjcode = WCS_PIX;
	    strcpy (wcs->c1type, wcstemp);
	    strcpy (wcs->ptype, wcstemp);
	    }

	/* Set up right ascension, declination, latitude, or longitude */
	else if (wcstemp[0] == 'R' ||
		 wcstemp[0] == 'D' ||
		 wcstemp[0] == 'A' ||
		 wcstemp[1] == 'L') {
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
	    wcs->ptype[0] = wcstemp[5];
	    wcs->ptype[1] = wcstemp[6];
	    wcs->ptype[2] = wcstemp[7];
	    wcs->ptype[3] = 0;

	    /*  Find projection type  */
	    wcs->prjcode = 0;  /* default type is linear */
	    for (i = 1; i < nctype; i++) {
		if (!strncmp(wcs->ptype, ctypes[i], 3))
		    wcs->prjcode = i;
		}

	    /* Handle obsolete projection */
	    if (wcs->prjcode == WCS_NCP)
		wcs->oldwcs = 1;
	    if (wcs->oldwcs && (
		wcs->prjcode != WCS_STG && wcs->prjcode != WCS_AIT &&
		wcs->prjcode != WCS_MER && wcs->prjcode != WCS_GLS &&
		wcs->prjcode != WCS_ARC && wcs->prjcode != WCS_TAN &&
		wcs->prjcode != WCS_TNX && wcs->prjcode != WCS_SIN &&
		wcs->prjcode != WCS_PIX && wcs->prjcode != WCS_LPR))
		wcs->oldwcs = 0;

	    /* Handle NOAO corrected TNX as TAN if oldwcs is set */
	    if (wcs->oldwcs && wcs->prjcode == WCS_TNX) {
		wcs->ctype[0][6] = 'A';
		wcs->ctype[0][7] = 'N';
		wcs->prjcode = WCS_TAN;
		}
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
	    wcs->prjcode = 0;
	    strcpy (wcs->c2type, wcstemp);
	    if (!hgets (hstring, "CUNIT2", 16, wcs->units[1])) {
		if (!mgets (hstring, "WAT2", "units", 16, wcs->units[1])) {
		    wcs->units[1][0] = 0;
		    }
		}
	    }

	/* Deal appropriately with pixel coordinates */
	else if (!strncmp (wcstemp,"PIXEL",6)) {
	    wcs->prjcode = WCS_PIX;
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
	    strcpy (wcs->ctype[1], wcstemp);

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
	wcs->crpix[0] = 1.0;
	hgetr8 (hstring,"CRPIX1",&wcs->crpix[0]);
	wcs->crpix[1] = 1.0;
	hgetr8 (hstring,"CRPIX2",&wcs->crpix[1]);
	wcs->xrefpix = wcs->crpix[0];
	wcs->yrefpix = wcs->crpix[1];
	wcs->crval[0] = 0.0;
	hgetr8 (hstring,"CRVAL1",&wcs->crval[0]);
	wcs->crval[1] = 0.0;
	hgetr8 (hstring,"CRVAL2",&wcs->crval[1]);
	wcs->xref = wcs->crval[0];
	wcs->yref = wcs->crval[1];
	if (wcs->coorflip) {
	    wcs->cel.ref[0] = wcs->crval[1];
	    wcs->cel.ref[1] = wcs->crval[0];
	    }
	else {
	    wcs->cel.ref[0] = wcs->crval[0];
	    wcs->cel.ref[1] = wcs->crval[1];
	    }
	wcs->longpole = 999.0;
	hgetr8 (hstring,"LONGPOLE",&wcs->longpole);
	wcs->cel.ref[2] = wcs->longpole;
	wcs->latpole = 999.0;
	hgetr8 (hstring,"LATPOLE",&wcs->latpole);
	wcs->cel.ref[3] = wcs->latpole;

	/* Use polynomial fit instead of projection, if present */
	wcs->ncoeff1 = 0;
	wcs->ncoeff2 = 0;
	if (!wcs->oldwcs && (hcoeff = ksearch (hstring,"CO1_1")) != NULL) {
	    wcs->prjcode = WCS_PLT;
	    (void)strcpy (wcs->ptype, "PLATE");
	    for (i = 0; i < 20; i++) {
		sprintf (keyword,"CO1_%d", i+1);
		wcs->x_coeff[i] = 0.0;
		if (hgetr8 (hcoeff, keyword, &wcs->x_coeff[i]))
		    wcs->ncoeff1 = i + 1;
		}
	    hcoeff = ksearch (hstring,"CO2_1");
	    for (i = 0; i < 20; i++) {
		sprintf (keyword,"CO2_%d",i+1);
		wcs->y_coeff[i] = 0.0;
		if (hgetr8 (hcoeff, keyword, &wcs->y_coeff[i]))
		    wcs->ncoeff2 = i + 1;
		}
	    platepos (wcs->crpix[0], wcs->crpix[1], wcs, &ra0, &dec0);
	    platepos (wcs->crpix[0]+100.0, wcs->crpix[1], wcs, &ra1, &dec1);
	    radiff = -(ra1 - ra0) / cos (degrad(0.5*dec0+dec1));
	    wcs->rot = atan2 ((dec1 - dec0), radiff);
	    wcs->cdelt[0] = 0.01 * wcsdist (ra0, dec0, ra1, dec1);
	    wcs->xinc = wcs->cdelt[0];
	    platepos (wcs->crpix[0], wcs->crpix[1]+100.0, wcs, &ra1, &dec1);
	    wcs->cdelt[1] = 0.01 * wcsdist (ra0, dec0, ra1, dec1);
	    wcs->yinc = wcs->cdelt[1];
	    }
	else if (hgetr8 (hstring,"CD1_1",&wcs->cd[0]) != 0) {
	    wcs->rotmat = 1;
	    wcs->cd[1] = 0.;
	    hgetr8 (hstring,"CD1_2",&wcs->cd[1]);
	    wcs->cd[2] = 0.;
	    hgetr8 (hstring,"CD2_1",&wcs->cd[2]);
	    wcs->cd[3] = wcs->cd[0];
	    hgetr8 (hstring,"CD2_2",&wcs->cd[3]);
	    (void) matinv (2, wcs->cd, wcs->dc);
	    wcs->xinc = sqrt (wcs->cd[0]*wcs->cd[0] +
			      wcs->cd[2]*wcs->cd[2]);
	    wcs->yinc = sqrt (wcs->cd[1]*wcs->cd[1] +
			      wcs->cd[3]*wcs->cd[3]);
	    sdelt = (wcs->cd[0] * wcs->cd[3]) - (wcs->cd[1]*wcs->cd[2]);
	    if (sdelt < 0.0) {
		if (wcs->coorflip) {
		    wcs->xinc = -wcs->xinc;
		    wcs->yinc = -wcs->yinc;
		    }
		else
		    wcs->xinc = -wcs->xinc;
		wcs->rot = raddeg (atan2 (-wcs->cd[1], wcs->cd[3]));
		}
	    else
		wcs->rot = raddeg (atan2 (wcs->cd[1], wcs->cd[3]));
	    if (wcs->coorflip)
		wcs->rot = wcs->rot - 90.0;
	    wcs->cdelt[0] = wcs->xinc;
	    wcs->cdelt[1] = wcs->yinc;
	    mem = 4 * sizeof(double);
	    wcs->lin.piximg = (double*)malloc(mem);
	    if (wcs->lin.piximg != NULL) {
		wcs->lin.imgpix = (double*)malloc(mem);
		if (wcs->lin.imgpix != NULL) {
		    wcs->lin.flag = LINSET;
		    for (i = 0; i < 4; i++) {
			wcs->lin.piximg[i] = wcs->cd[i];
			}
		    if (wcs->coorflip) {
			wcs->lin.piximg[1] = -wcs->cd[2];
			wcs->lin.piximg[2] = -wcs->cd[1];
			}
		    (void) matinv (2, wcs->lin.piximg, wcs->lin.imgpix);
		    }
		}
	    }
	else if (hgetr8 (hstring,"CDELT1",&wcs->xinc) != 0) {
	    wcs->yinc = wcs->xinc;
	    hgetr8 (hstring,"CDELT2",&wcs->yinc);
	    wcs->cdelt[0] = wcs->xinc;
	    wcs->cdelt[1] = wcs->yinc;
	    if (naxes > 2)
		hgetr8 (hstring,"CDELT3",&wcs->cdelt[2]);
	    if (naxes > 3)
		hgetr8 (hstring,"CDELT4",&wcs->cdelt[3]);
	    pci = wcs->pc;
	    for (i = 0; i < naxes; i++) {
		for (j = 0; j < naxes; j++) {
		    if (i ==j)
			*pci = 1.0;
		    else
			*pci = 0.0;
		    pci++;
		    }
		}
	    wcs->rot = 0.;
	    wcs->cd[0] = 1.;
	    wcs->cd[2] = 0.;
	    wcs->cd[1] = 0.;
	    wcs->cd[3] = 1.;
	    wcs->rotmat = 0;
	    if (hgetr8 (hstring,"PC001001",&wcs->cd[0]) != 0) {
		wcs->cd[1] = 0.0;
		hgetr8 (hstring,"PC001002",&wcs->cd[1]);
		wcs->cd[2] = 0.0;
		hgetr8 (hstring,"PC002001",&wcs->cd[2]);
		wcs->cd[3] = wcs->cd[1];
		hgetr8 (hstring,"PC002002",&wcs->cd[3]);
		wcs->pc[0] = wcs->cd[0];
		wcs->pc[1] = wcs->cd[1];
		wcs->pc[naxes] = wcs->cd[2];
		wcs->pc[naxes+1] = wcs->cd[3];
		wcs->cd[0] = wcs->cd[0] * wcs->xinc;
		wcs->cd[1] = wcs->cd[1] * wcs->yinc;
		wcs->cd[2] = wcs->cd[2] * wcs->xinc;
		wcs->cd[3] = wcs->cd[3] * wcs->yinc;
		wcs->rotmat = 1;
		if (naxes > 2) {
		    wcs->pc[2] = 0.0;
		    hgetr8 (hstring,"PC001003",&wcs->pc[2]);
		    wcs->pc[naxes+2] = 0.0;
		    hgetr8 (hstring,"PC002003",&wcs->pc[naxes+2]);
		    wcs->pc[2*naxes] = 0.0;
		    hgetr8 (hstring,"PC003001",&wcs->pc[2*naxes]);
		    wcs->pc[(2*naxes)+1] = 0.0;
		    hgetr8 (hstring,"PC003002",&wcs->pc[(2*naxes)+1]);
		    wcs->pc[(2*naxes)+2] = 1.0;
		    hgetr8 (hstring,"PC003003",&wcs->pc[(2*naxes)+2]);
		}
		if (naxes > 3) {
		    wcs->pc[3] = 0.0;
		    hgetr8 (hstring,"PC001004",&wcs->pc[3]);
		    wcs->pc[naxes+2] = 0.0;
		    hgetr8 (hstring,"PC002004",&wcs->pc[naxes+3]);
		    wcs->pc[(2*naxes)+3] = 0.0;
		    hgetr8 (hstring,"PC003004",&wcs->pc[(2*naxes)+3]);
		    wcs->pc[3*naxes] = 0.0;
		    hgetr8 (hstring,"PC004001",&wcs->pc[3*naxes]);
		    wcs->pc[(3*naxes)+1] = 0.0;
		    hgetr8 (hstring,"PC004002",&wcs->pc[(3*naxes)+1]);
		    wcs->pc[(3*naxes)+2] = 0.0;
		    hgetr8 (hstring,"PC004003",&wcs->pc[(3*naxes)+2]);
		    wcs->pc[(3*naxes)+3] = 1.0;
		    hgetr8 (hstring,"PC004004",&wcs->pc[(3*naxes)+3]);
		    }
		}
	    else {
		hgetr8 (hstring,"CROTA1",&wcs->rot);
		if (wcs->rot == 0.)
		    hgetr8 (hstring,"CROTA2",&wcs->rot);
		s = wcs->cdelt[1] / wcs->cdelt[0];
		crot = cos (degrad(wcs->rot));
		srot = sin (degrad(wcs->rot));
		wcs->pc[0] = crot;
		wcs->pc[1] = -srot * s;
		wcs->pc[naxes] = srot / s;
		wcs->pc[naxes+1] = crot;
		wcs->cd[0] = wcs->xinc * crot;
		wcs->cd[1] = -wcs->yinc * srot;
		wcs->cd[2] = wcs->xinc * srot;
		wcs->cd[3] = wcs->yinc * crot;
		if (naxes > 2) {
		    wcs->pc[2] = 0.0;
		    wcs->pc[naxes+2] = 0.0;
		    wcs->pc[2*naxes] = 0.0;
		    wcs->pc[(2*naxes)+1] = 0.0;
		    wcs->pc[(2*naxes)+2] = 1.0;
		    }
		if (naxes > 3) {
		    wcs->pc[3] = 0.0;
		    wcs->pc[naxes+2] = 0.0;
		    wcs->pc[(2*naxes)+3] = 0.0;
		    wcs->pc[3*naxes] = 0.0;
		    wcs->pc[(3*naxes)+1] = 0.0;
		    wcs->pc[(3*naxes)+2] = 0.0;
		    wcs->pc[(3*naxes)+3] = 1.0;
		    }
		}
	    }
	else {
	    wcs->xinc = 1.0;
	    wcs->yinc = 1.0;
	    (void)sprintf (wcserrmsg,"WCSINIT: setting CDELT to 1\n");
	    }

	/* Initialize TNX, defaulting to TAN if there is a problem */
	if (wcs->prjcode == WCS_TNX) {
	    if (tnxinit (hstring, wcs)) {
		wcs->ctype[1][6] = 'A';
		wcs->ctype[1][7] = 'N';
		wcs->prjcode = WCS_TAN;
		}
	    }

	/* Coordinate reference frame, equinox, and epoch */
	if (strncmp (wcs->ptype,"LINEAR",6) &&
	    strncmp (wcs->ptype,"PIXEL",5))
	    wcseq (hstring,wcs);
	else {
	    wcs->degout = -1;
	    wcs->ndec = 5;
	    }

	wcs->wcson = 1;
	}

    /* Plate solution coefficients */
    else if (ksearch (hstring,"PLTRAH") != NULL) {
	wcs->prjcode = WCS_DSS;
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
	for (i = 0; i < 6; i++) {
	    sprintf (keyword,"PPO%d", i+1);
	    wcs->ppo_coeff[i] = 0.0;
	    hgetr8 (hcoeff,keyword,&wcs->ppo_coeff[i]);
	    }
	hcoeff = ksearch (hstring,"AMDX1");
	for (i = 0; i < 20; i++) {
	    sprintf (keyword,"AMDX%d", i+1);
	    wcs->x_coeff[i] = 0.0;
	    hgetr8 (hcoeff, keyword, &wcs->x_coeff[i]);
	    }
	hcoeff = ksearch (hstring,"AMDY1");
	for (i = 0; i < 20; i++) {
	    sprintf (keyword,"AMDY%d",i+1);
	    wcs->y_coeff[i] = 0.0;
	    hgetr8 (hcoeff, keyword, &wcs->y_coeff[i]);
	    }
	wcs->wcson = 1;
	wcs->wcson = 1;
	(void)strcpy (wcs->c1type, "RA");
	(void)strcpy (wcs->c2type, "DEC");
	(void)strcpy (wcs->ptype, "DSS");
	wcs->degout = 0;
	wcs->ndec = 3;
	wcs->crpix[0] = 0.5 * wcs->nxpix;
	wcs->crpix[1] = 0.5 * wcs->nypix;
	wcs->xrefpix = wcs->crpix[0];
	wcs->yrefpix = wcs->crpix[1];
	dsspos (wcs->crpix[0], wcs->crpix[1], wcs, &ra0, &dec0);
	wcs->crval[0] = ra0;
	wcs->crval[1] = dec0;
	wcs->xref = wcs->crval[0];
	wcs->yref = wcs->crval[1];
	dsspos (wcs->crpix[0]+100.0, wcs->crpix[1], wcs, &ra1, &dec1);
	radiff = -(ra1 - ra0) / cos (degrad(0.5*dec0+dec1));
	wcs->rot = atan2 ((dec1 - dec0), radiff);
	wcs->cdelt[0] = 0.01 * wcsdist (ra0, dec0, ra1, dec1);
	wcs->xinc = wcs->cdelt[0];
	dsspos (wcs->crpix[0], wcs->crpix[1]+100.0, wcs, &ra1, &dec1);
	wcs->cdelt[1] = 0.01 * wcsdist (ra0, dec0, ra1, dec1);
	wcs->yinc = wcs->cdelt[1];
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
	wcs->cdelt[0] = wcs->xinc;
	wcs->cdelt[1] = wcs->yinc;

	/* Get rotation angle from the header, if it's there */
	if (ksearch (hstring,"CROTA2") != NULL)
	    hgetr8 (hstring,"CROTA2",&wcs->rot);
	else
	    wcs->rot = 0.0;

	/* Set CD and PC matrices */
	s = wcs->cdelt[1] / wcs->cdelt[0];
	crot = cos (degrad(wcs->rot));
	srot = sin (degrad(wcs->rot));
	wcs->pc[0] = crot;
	wcs->pc[1] = -srot * s;
	wcs->pc[naxes] = srot / s;
	wcs->pc[naxes+1] = crot;
	wcs->cd[0] = wcs->xinc * crot;
	wcs->cd[1] = -wcs->yinc * srot;
	wcs->cd[2] = wcs->xinc * srot;
	wcs->cd[3] = wcs->yinc * crot;
	if (naxes > 2) {
	    wcs->pc[2] = 0.0;
	    wcs->pc[naxes+2] = 0.0;
	    wcs->pc[2*naxes] = 0.0;
	    wcs->pc[(2*naxes)+1] = 0.0;
	    wcs->pc[(2*naxes)+2] = 1.0;
	    }
	if (naxes > 3) {
	    wcs->pc[3] = 0.0;
	    wcs->pc[naxes+2] = 0.0;
	    wcs->pc[(2*naxes)+3] = 0.0;
	    wcs->pc[3*naxes] = 0.0;
	    wcs->pc[(3*naxes)+1] = 0.0;
	    wcs->pc[(3*naxes)+2] = 0.0;
	    wcs->pc[(3*naxes)+3] = 1.0;
	    }

	/* By default, set reference pixel to center of image */
	wcs->crpix[0] = wcs->nxpix * 0.5;
	wcs->crpix[1] = wcs->nypix * 0.5;

	/* Get reference pixel from the header, if it's there */
	if (ksearch (hstring,"CRPIX1") != NULL) {
	    hgetr8 (hstring,"CRPIX1",&wcs->crpix[0]);
	    hgetr8 (hstring,"CRPIX2",&wcs->crpix[1]);
	    }

	/* Use center of detector array as reference pixel
	else if (ksearch (hstring,"DETSIZE") != NULL ||
		 ksearch (hstring,"DETSEC") != NULL) {
	    hgets (hstring, "DETSIZE", 32, temp);
	    ic = strchr (temp, ':');
	    if (ic != NULL)
		*ic = ' ';
	    ic = strchr (temp, ',');
	    if (ic != NULL)
		*ic = ' ';
	    ic = strchr (temp, ':');
	    if (ic != NULL)
		*ic = ' ';
	    ic = strchr (temp, ']');
	    if (ic != NULL)
		*ic = (char) 0;
	    sscanf (temp, "%d %d %d %d", &idx1, &idx2, &idy1, &idy2);
	    dxrefpix = 0.5 * (double) (idx1 + idx2 - 1);
	    dyrefpix = 0.5 * (double) (idy1 + idy2 - 1);
	    hgets (hstring, "DETSEC", 32, temp);
	    ic = strchr (temp, ':');
	    if (ic != NULL)
		*ic = ' ';
	    ic = strchr (temp, ',');
	    if (ic != NULL)
		*ic = ' ';
	    ic = strchr (temp, ':');
	    if (ic != NULL)
		*ic = ' ';
	    ic = strchr (temp, ']');
	    if (ic != NULL)
		*ic = (char) 0;
	    sscanf (temp, "%d %d %d %d", &ix1, &ix2, &iy1, &iy2);
	    wcs->crpix[0] = dxrefpix - (double) (ix1 - 1);
	    wcs->crpix[1] = dyrefpix - (double) (iy1 - 1);
	    } */
	wcs->xrefpix = wcs->crpix[0];
	wcs->yrefpix = wcs->crpix[1];

	/* Get rotation angle from the header, if it's there */
	if (ksearch (hstring,"CROTA2") != NULL) {
	    hgetr8 (hstring,"CROTA2",&wcs->rot);
	    }

	wcs->crval[0] = 0.0;
	if (!hgetra (hstring,"RA",&wcs->crval[0])) {
	    (void)sprintf (wcserrmsg,"WCSINIT: No RA with SECPIX, no WCS\n");
	    free (wcs);
	    return (NULL);
	    }
	wcs->crval[1] = 0.0;
	if (!hgetdec (hstring,"DEC",&wcs->crval[1])) {
	    (void)sprintf (wcserrmsg,"WCSINIT No DEC with SECPIX, no WCS\n");
	    free (wcs);
	    return (NULL);
	    }
	wcs->xref = wcs->crval[0];
	wcs->yref = wcs->crval[1];
	wcs->coorflip = 0;

	wcs->cel.ref[0] = wcs->crval[0];
	wcs->cel.ref[1] = wcs->crval[1];
	wcs->cel.ref[2] = 999.0;
	hgetr8 (hstring,"LONGPOLE",&wcs->cel.ref[2]);
	wcs->cel.ref[3] = 999.0;
	hgetr8 (hstring,"LATPOLE",&wcs->cel.ref[3]);
	    
	strcpy (wcs->ctype[0], "RA---TAN");
	strcpy (wcs->ctype[1], "DEC--TAN");
	strcpy (wcs->c1type,"RA");
	strcpy (wcs->c2type,"DEC");
	strcpy (wcs->ptype,"TAN");
	wcs->prjcode = WCS_TAN;
	wcs->coorflip = 0;
	wcs->rot = 0.;
	wcs->degout = 0;
	wcs->ndec = 3;
	hgetr8 (hstring,"CROTA1",&wcs->rot);
	if (wcs->rot == 0.)
	    hgetr8 (hstring,"CROTA2",&wcs->rot);
	wcs->dc[0] = 0.;
	wcs->dc[1] = 0.;
	wcs->dc[2] = 0.;
	wcs->dc[3] = 0.;
	wcs->rotmat = 0;
	s = wcs->cdelt[1] / wcs->cdelt[0];
	crot = cos (degrad(wcs->rot));
	srot = sin (degrad(wcs->rot));
	wcs->pc[0] = crot;
	wcs->pc[1] = -srot * s;
	wcs->pc[2] = srot / s;
	wcs->pc[3] = crot;

	/* Coordinate reference frame and equinox */
	wcseq (hstring,wcs);

	/* Epoch of image (from observation date, if possible) */
	if (hgetr8 (hstring, "MJD-OBS", &mjd))
	    wcs->epoch = 1950.0 + (mjd / 365.22);
	else if (!hgetdate (hstring,"DATE-OBS",&wcs->epoch)) {
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

    wcs->lin.crpix = wcs->crpix;
    wcs->lin.cdelt = wcs->cdelt;
    wcs->lin.pc = wcs->pc;
    if (strlen (wcs->radecsys) == 0 || wcs->prjcode == WCS_LPR)
	strcpy (wcs->radecsys, "LINEAR");
    wcs->syswcs = wcscsys (wcs->radecsys);
    if (wcs->syswcs == WCS_B1950)
	strcpy (wcs->radecout, "B1950");
    else if (wcs->syswcs == WCS_J2000)
	strcpy (wcs->radecout, "J2000");
    else
	strcpy (wcs->radecout, wcs->radecsys);
    wcs->sysout = wcscsys (wcs->radecout);
    wcs->eqout = 0.0;
    strcpy (wcs->radecin, wcs->radecsys);
    wcs->sysin = wcscsys (wcs->radecin);
    wcs->printsys = 1;
    wcs->tabsys = 0;
    wcs->linmode = 0;
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
    wcstemp[0] = 0;
    hgets (hstring,"EQUINOX",16,wcstemp);
    if (wcstemp[0] == 'J') {
	wcs->equinox = atof (wcstemp+1);
	ieq = atoi (wcstemp+1);
	}
    else if (wcstemp[0] == 'B') {
	wcs->equinox = atof (wcstemp+1);
	ieq = atoi (wcstemp+1);
	}
    else if (hgeti4 (hstring,"EQUINOX",&ieq))
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
    else if (wcs->c1type[0] == 'L')
	strcpy (wcs->radecsys,"LINEAR");

    /* Otherwise set coordinate system from equinox */
    /* Systemless coordinates cannot be translated using b, j, or g commands */
    else {
	if (ieq > 1980)
	    strcpy (wcs->radecsys,"FK5");
	else
	    strcpy (wcs->radecsys,"FK4");
	}
    wcs->syswcs = wcscsys (wcs->radecsys);

    return;
}


/* Set up a WCS structure from subroutine arguments */

struct WorldCoor *
wcsxinit (cra,cdec,secpix,xrpix,yrpix,nxpix,nypix,rotate,equinox,epoch,proj)

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
    double s, srot, crot;

    wcs = (struct WorldCoor *) calloc (1, sizeof(struct WorldCoor));

    /* Set WCSLIB flags so that structures will be reinitialized */
    wcs->cel.flag = 0;
    wcs->lin.flag = 0;
    wcs->wcsl.flag = 0;

    /* Plate solution coefficients */
    wcs->nxpix = nxpix;
    wcs->nypix = nypix;

    wcs->oldwcs = oldwcs0;

    /* Approximate world coordinate system from a known plate scale */
    wcs->yinc = secpix / 3600.0;
    wcs->xinc = -wcs->yinc;
    wcs->cdelt[0] = wcs->xinc;
    wcs->cdelt[1] = wcs->yinc;

    wcs->crpix[0] = xrpix;
    wcs->crpix[1] = yrpix;
    wcs->xrefpix = wcs->crpix[0];
    wcs->yrefpix = wcs->crpix[1];

    wcs->crval[0] = cra;
    wcs->crval[1] = cdec;
    wcs->xref = wcs->crval[0];
    wcs->yref = wcs->crval[1];

    strcpy (wcs->c1type,"RA");
    strcpy (wcs->c2type,"DEC");

/* Allan Brighton: 28.4.98: for backward compat., remove leading "--" */
    while (proj && *proj == '-')
	proj++;
    strcpy (wcs->ptype,proj);
    strcpy (wcs->ctype[0],"RA---");
    strcpy (wcs->ctype[1],"DEC--");
    strcat (wcs->ctype[0],proj);
    strcat (wcs->ctype[1],proj);
    
    wcs->prjcode = WCS_TAN;
    wcs->coorflip = 0;
    wcs->rot = rotate;
    wcs->rotmat = 0;
    s = wcs->cdelt[1] / wcs->cdelt[0];
    crot = cos (degrad(wcs->rot));
    srot = sin (degrad(wcs->rot));
    wcs->pc[0] = crot;
    wcs->pc[1] = -srot * s;
    wcs->pc[2] = srot / s;
    wcs->pc[3] = crot;
    wcs->cd[0] = wcs->xinc * crot;
    wcs->cd[1] = -wcs->yinc * srot;
    wcs->cd[2] = wcs->xinc * srot;
    wcs->cd[3] = wcs->yinc * crot;
    wcs->cel.ref[0] = wcs->crval[0];
    wcs->cel.ref[1] = wcs->crval[1];
    wcs->cel.ref[2] = 999.0;

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

    strcpy (wcs->radecout, wcs->radecsys);
    wcs->syswcs = wcscsys (wcs->radecsys);
    wcsoutinit (wcs, wcs->radecsys);
    wcsininit (wcs, wcs->radecsys);
    wcs->eqout = 0.0;
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


int
wcsreset (wcs, crpix1, crpix2, crval1, crval2, cdelt1, cdelt2, crota, cd)

struct WorldCoor *wcs;		/* World coordinate system data structure */
double crpix1, crpix2;		/* Reference pixel coordinates */
double crval1, crval2;		/* Coordinates at reference pixel in degrees */
double cdelt1, cdelt2;		/* scale in degrees/pixel, ignored if cd is not NULL */
double crota;			/* Rotation angle in degrees, ignored if cd is not NULL */
double *cd;			/* Rotation matrix, used if not NULL */
{
    int i, j, mem;
    double *pci;
    double s, srot, crot;
    extern int matinv();

    if (nowcs (wcs))
	return (-1);

    /* Set WCSLIB flags so that structures will be reinitialized */
    wcs->cel.flag = 0;
    wcs->lin.flag = 0;
    wcs->wcsl.flag = 0;

    /* Reference pixel coordinates and WCS value */
    wcs->crpix[0] = crpix1;
    wcs->crpix[1] = crpix2;
    wcs->xrefpix = wcs->crpix[0];
    wcs->yrefpix = wcs->crpix[1];
    wcs->crval[0] = crval1;
    wcs->crval[1] = crval2;
    wcs->xref = wcs->crval[0];
    wcs->yref = wcs->crval[1];
    if (wcs->coorflip) {
	wcs->cel.ref[1] = wcs->crval[0];
	wcs->cel.ref[0] = wcs->crval[1];
	}
    else {
	wcs->cel.ref[0] = wcs->crval[0];
	wcs->cel.ref[1] = wcs->crval[1];
	}
    /* Keep ref[2] and ref[3] from input */

    if (cd != NULL) {
	wcs->rotmat = 1;
	wcs->cd[0] = cd[0];
	wcs->cd[1] = cd[1];
	wcs->cd[2] = cd[2];
	wcs->cd[3] = cd[3];
	(void) matinv (2, wcs->cd, wcs->dc);
	wcs->xinc = sqrt (wcs->cd[0]*wcs->cd[0] + wcs->cd[2]*wcs->cd[2]);
	wcs->yinc = sqrt (wcs->cd[1]*wcs->cd[1] + wcs->cd[3]*wcs->cd[3]);
	if ((wcs->cd[0]*wcs->cd[0] - wcs->cd[1]*wcs->cd[1]) < 0) {
	    if (!strncmp(wcs->c1type,"RA",2) || !strncmp(wcs->c1type,"GLON",4))
		wcs->xinc = -wcs->xinc;
	    if (!strncmp(wcs->c2type,"RA",2) || !strncmp(wcs->c2type,"GLON",4))
		wcs->yinc = -wcs->yinc;
	    wcs->rot = raddeg (atan2 (-wcs->cd[1], wcs->cd[3]));
	    }
	else
	    wcs->rot = raddeg (atan2 (wcs->cd[1], wcs->cd[3]));
	wcs->cdelt[0] = wcs->xinc;
	wcs->cdelt[1] = wcs->yinc;
	mem = 4 * sizeof(double);
	wcs->lin.piximg = (double*)malloc(mem);
	if (wcs->lin.piximg != NULL) {
	    wcs->lin.imgpix = (double*)malloc(mem);
	    if (wcs->lin.imgpix != NULL) {
		wcs->lin.flag = LINSET;
		for (i = 0; i < 4; i++) {
		    wcs->lin.piximg[i] = wcs->cd[i];
		    }
		if (wcs->coorflip) {
		    wcs->lin.piximg[1] = -wcs->cd[2];
		    wcs->lin.piximg[2] = -wcs->cd[1];
		    }
		(void) matinv (2, wcs->lin.piximg, wcs->lin.imgpix);
		wcs->lin.flag = LINSET;
		}
	    }
	}
    else if (cdelt1 != 0.0) {
	wcs->xinc = cdelt1;
	if (cdelt2 != 0.0)
	    wcs->yinc = cdelt2;
	else
	    wcs->yinc = cdelt1;
	wcs->cdelt[0] = wcs->xinc;
	wcs->cdelt[1] = wcs->yinc;
	pci = wcs->pc;
	for (i = 0; i < wcs->lin.naxis; i++) {
	    for (j = 0; j < wcs->lin.naxis; j++) {
		if (i ==j)
		    *pci = 1.0;
		else
		    *pci = 0.0;
		pci++;
		}
	    }
	wcs->rotmat = 0;
	wcs->rot = crota;
	s = wcs->cdelt[1] / wcs->cdelt[0];
	crot = cos (degrad(wcs->rot));
	srot = sin (degrad(wcs->rot));
	wcs->cd[0] = wcs->xinc * crot;
	wcs->cd[1] = -wcs->yinc * srot;
	wcs->cd[2] = wcs->xinc * srot;
	wcs->cd[3] = wcs->yinc * crot;
	wcs->pc[0] = crot;
	wcs->pc[1] = -srot * s;
	wcs->pc[wcs->lin.naxis] = srot / s;
	wcs->pc[wcs->lin.naxis+1] = crot;
	wcs->lin.flag = 0;
	}
    else {
	wcs->xinc = 1.0;
	wcs->yinc = 1.0;
	(void)sprintf (wcserrmsg,"WCSRESET: setting CDELT to 1\n");
	}

    /* Coordinate reference frame, equinox, and epoch */
    if (!strncmp (wcs->ptype,"LINEAR",6) ||
	!strncmp (wcs->ptype,"PIXEL",5))
	wcs->degout = -1;

    /* Initialize to no plate fit */
    wcs->ncoeff1 = 0;
    wcs->ncoeff2 = 0;

    wcs->wcson = 1;
    return (0);
}

void
wcseqset (wcs, equinox)

struct WorldCoor *wcs;		/* World coordinate system data structure */
double equinox;			/* Desired equinox as fractional year */
{
    extern void fk425e(), fk524e();

    if (nowcs (wcs))
	return;

    /* Leave WCS alone if already at desired equinox */
    if (wcs->equinox == equinox)
	return;

    /* Convert center from B1950 (FK4) to J2000 (FK5) */
    if (equinox == 2000.0 && wcs->equinox == 1950.0) {
	if (wcs->coorflip) { 
	    fk425e (&wcs->crval[1], &wcs->crval[0], wcs->epoch);
	    wcs->cel.ref[1] = wcs->crval[0];
	    wcs->cel.ref[0] = wcs->crval[1];
	    }
	else {
	    fk425e (&wcs->crval[0], &wcs->crval[1], wcs->epoch);
	    wcs->cel.ref[0] = wcs->crval[0];
	    wcs->cel.ref[1] = wcs->crval[1];
	    }
	wcs->xref = wcs->crval[0];
	wcs->yref = wcs->crval[1];
	wcs->equinox = 2000.0;
	strcpy (wcs->radecsys, "FK5");
	wcs->syswcs = WCS_J2000;
	wcs->cel.flag = 0;
	wcs->wcsl.flag = 0;
	}

    /* Convert center from J2000 (FK5) to B1950 (FK4) */
    else if (equinox == 1950.0 && wcs->equinox == 2000.0) {
	if (wcs->coorflip) { 
	    fk524e (&wcs->crval[1], &wcs->crval[0], wcs->epoch);
	    wcs->cel.ref[1] = wcs->crval[0];
	    wcs->cel.ref[0] = wcs->crval[1];
	    }
	else {
	    fk524e (&wcs->crval[0], &wcs->crval[1], wcs->epoch);
	    wcs->cel.ref[0] = wcs->crval[0];
	    wcs->cel.ref[1] = wcs->crval[1];
	    }
	wcs->xref = wcs->crval[0];
	wcs->yref = wcs->crval[1];
	wcs->equinox = 1950.0;
	strcpy (wcs->radecsys, "FK4");
	wcs->syswcs = WCS_B1950;
	wcs->cel.flag = 0;
	wcs->wcsl.flag = 0;
	}
    wcsoutinit (wcs, wcs->radecsys);
    wcsininit (wcs, wcs->radecsys);
    return;
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
wcsshift (wcs,rra,rdec,coorsys)

struct WorldCoor *wcs;	/* World coordinate system structure */
double	rra;		/* Reference pixel right ascension in degrees */
double	rdec;		/* Reference pixel declination in degrees */
char	*coorsys;	/* FK4 or FK5 coordinates (1950 or 2000) */

{
    if (nowcs (wcs))
	return;

/* Approximate world coordinate system from a known plate scale */
    wcs->crval[0] = rra;
    wcs->crval[1] = rdec;
    wcs->xref = wcs->crval[0];
    wcs->yref = wcs->crval[1];


/* Coordinate reference frame */
    strcpy (wcs->radecsys,coorsys);
    wcs->syswcs = wcscsys (coorsys);
    if (wcs->syswcs == WCS_B1950)
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
    double width, height, secpix, secpixh, secpixw;
    int lstr = 32;

    if (nowcs (wcs))
	(void)fprintf (stderr,"No WCS information available\n");
    else {
	if (wcs->prjcode == WCS_DSS)
	    (void)fprintf (stderr,"WCS plate center  %s\n", wcs->center);
	xpix = 0.5 * wcs->nxpix;
	ypix = 0.5 * wcs->nypix;
	(void) pix2wcst (wcs,xpix,ypix,wcstring, lstr);
	(void)fprintf (stderr,"WCS center %s %s %s %s at pixel (%.2f,%.2f)\n",
		     wcs->ctype[0],wcs->ctype[1],wcstring,wcs->ptype,xpix,ypix);

	/* Image width */
	(void) pix2wcs (wcs,1.0,ypix,&xpos1,&ypos1);
	(void) pix2wcs (wcs,wcs->nxpix,ypix,&xpos2,&ypos2);
	if (wcs->syswcs == WCS_LINEAR) {
	    width = xpos2 - xpos1;
	    if (width < 100.0)
	    (void)fprintf (stderr, "WCS width = %.5f %s ",width, wcs->units[0]);
	    else
	    (void)fprintf (stderr, "WCS width = %.3f %s ",width, wcs->units[0]);
	    }
	else {
	    width = wcsdist (xpos1,ypos1,xpos2,ypos2);
	    if (width < 1/60.0)
		(void)fprintf (stderr, "WCS width = %.2f arcsec ",width*3600.0);
	    else if (width < 1.0)
		(void)fprintf (stderr, "WCS width = %.2f arcmin ",width*60.0);
	    else
		(void)fprintf (stderr, "WCS width = %.3f degrees ",width);
	    }
	secpixw = width / (wcs->nxpix - 1.0);

	/* Image height */
	(void) pix2wcs (wcs,xpix,1.0,&xpos1,&ypos1);
	(void) pix2wcs (wcs,xpix,wcs->nypix,&xpos2,&ypos2);
	if (wcs->syswcs == WCS_LINEAR) {
	    height = ypos2 - ypos1;
	    if (width < 100.0)
	    (void)fprintf (stderr, " height = %.5f %s ",width, wcs->units[1]);
	    else
	    (void)fprintf (stderr, " height = %.3f %s ",width, wcs->units[1]);
	    }
	else {
	    height = wcsdist (xpos1,ypos1,xpos2,ypos2);
	    if (height < 1/60.0)
		(void) fprintf (stderr, " height = %.2f arcsec",height*3600.0);
	    else if (height < 1.0)
		(void) fprintf (stderr, " height = %.2f arcmin",height*60.0);
	    else
		(void) fprintf (stderr, " height = %.3f degrees",height);
	    }
	secpixh = height / (wcs->nypix - 1.0);

	/* Image scale */
	if (wcs->syswcs == WCS_LINEAR) {
	    (void) fprintf (stderr,"\n");
	    (void) fprintf (stderr,"WCS  %.5f %s/pixel, %.5f %s/pixel\n",
			    wcs->xinc,wcs->units[0],wcs->yinc,wcs->units[1]);
	    }
	else {
	    if (wcs->xinc != 0.0 && wcs->yinc != 0.0)
		secpix = (fabs(wcs->xinc) + fabs(wcs->yinc)) * 0.5 * 3600.0;
	    else if (secpixh > 0.0 && secpixw > 0.0)
		secpix = (secpixw + secpixh) * 0.5 * 3600.0;
	    else if (wcs->xinc != 0.0 || wcs->yinc != 0.0)
		secpix = (fabs(wcs->xinc) + fabs(wcs->yinc)) * 3600.0;
	    else
		secpix = (secpixw + secpixh) * 3600.0;
	    if (secpix < 100.0)
		(void) fprintf (stderr, "  %.3f arcsec/pixel\n",secpix);
	    else if (secpix < 3600.0)
		(void) fprintf (stderr, "  %.3f arcmin/pixel\n",secpix*60.0);
	    else
		(void) fprintf (stderr, "  %.3f degrees/pixel\n",secpix*3600.0);
	    }
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
    double width, height;

    /* Find right ascension and declination of coordinates */
    if (iswcs(wcs)) {
	wcsfull (wcs, cra, cdec, &width, &height);
	*dra = 0.5 * width;
	*ddec = 0.5 * height;
	}
    else {
	*cra = 0.0;
	*cdec = 0.0;
	*dra = 0.0;
	*ddec = 0.0;
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
    double xpix, ypix, xpos[4], ypos[4], xmin, xmax, ymin, ymax, xi, yi;
    double xcent, ycent, xlast, ylast;
    int i;

    /* Find right ascension and declination of coordinates */
    if (iswcs(wcs)) {
	xpix = 0.5 * wcs->nxpix;
	ypix = 0.5 * wcs->nypix;
	(void) pix2wcs (wcs,xpix,ypix,&xcent, &ycent);
	*cra = xcent;
	*cdec = ycent;

	/* Compute coordinates of all four corners of the image */
	xlast = wcs->nxpix;
	ylast = wcs->nypix;
	(void) pix2wcs (wcs, 1.0, 1.0, &xpos[0], &ypos[0]);
	(void) pix2wcs (wcs, 1.0, ylast, &xpos[1], &ypos[1]);
	(void) pix2wcs (wcs, xlast, 1.0, &xpos[2], &ypos[2]);
	(void) pix2wcs (wcs, xlast, ylast, &xpos[3], &ypos[3]);

	/* Find maximum and minimum ra/theta and dec/phi */
	xmax = xpos[0];
	xmin = xpos[0];
	ymax = ypos[0];
	ymin = ypos[0];
	for (i = 1; i < 4; i++) {
	    xi = xpos[i];
	    yi = ypos[i];
	    if (xi > xmax)
		xmax = xi;
	    if (xi < xmin)
		xmin = xi;
	    if (yi > ymax)
		ymax = yi;
	    if (yi < ymin)
		ymin = yi;
	    }

	/* Compute image width in degrees */
	*width = xmax - xmin;

	/* Compute image height in degrees */
	*height = ymax - ymin;
	}
    else {
	*cra = 0.0;
	*cdec = 0.0;
	*width = 0.0;
	*height = 0.0;
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

    if (nowcs (wcs))
	return;

    if (wcs->search_format[0] > 0)
	strcpy (comform, wcs->search_format);
    else
	strcpy (comform, "rgsc %s");

    if (nowcs (wcs))
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

/* Initialize WCS output coordinate system for use by PIX2WCS() */

void
wcsoutinit (wcs, coorsys)

struct WorldCoor *wcs;	/* World coordinate system structure */
char	*coorsys;	/* Input world coordinate system:
			   FK4, FK5, B1950, J2000, GALACTIC, ECLIPTIC
			   fk4, fk5, b1950, j2000, galactic, ecliptic */
{
    int sysout, i;

    if (nowcs (wcs))
	return;

    /* If argument is null, set to image system and equinox */
    if (coorsys == NULL || strlen (coorsys) < 1 ||
	!strcmp(coorsys,"IMSYS") || !strcmp(coorsys,"imsys")) {
	wcs->sysout = wcs->syswcs;
	strcpy (wcs->radecout, wcs->radecsys);
	wcs->eqout = wcs->equinox;
	if (wcs->sysout == WCS_B1950) {
	    if (wcs->eqout != 1950.0) {
		wcs->radecout[0] = 'B';
		sprintf (wcs->radecout+1,"%.4f", wcs->equinox);
		i = strlen(wcs->radecout) - 1;
		if (wcs->radecout[i] == '0')
		    wcs->radecout[i] = (char)0;
		i = strlen(wcs->radecout) - 1;
		if (wcs->radecout[i] == '0')
		    wcs->radecout[i] = (char)0;
		i = strlen(wcs->radecout) - 1;
		if (wcs->radecout[i] == '0')
		    wcs->radecout[i] = (char)0;
		}
	    else
		strcpy (wcs->radecout, "B1950");
	    }
	else if (wcs->sysout == WCS_J2000) {
	    if (wcs->eqout != 2000.0) {
		wcs->radecout[0] = 'J';
		sprintf (wcs->radecout+1,"%.4f", wcs->equinox);
		i = strlen(wcs->radecout) - 1;
		if (wcs->radecout[i] == '0')
		    wcs->radecout[i] = (char)0;
		i = strlen(wcs->radecout) - 1;
		if (wcs->radecout[i] == '0')
		    wcs->radecout[i] = (char)0;
		i = strlen(wcs->radecout) - 1;
		if (wcs->radecout[i] == '0')
		    wcs->radecout[i] = (char)0;
		}
	    else
		strcpy (wcs->radecout, "J2000");
	    }
	}

    /* Ignore new coordinate system if it is not supported */
    if ((sysout = wcscsys (coorsys)) < 0)
	return;

    /* Do not try to convert linear or alt-az coordinates */
    if (sysout != wcs->syswcs &&
	(wcs->syswcs == WCS_LINEAR || wcs->syswcs == WCS_ALTAZ))
	return;

    strcpy (wcs->radecout, coorsys);
    wcs->sysout = sysout;
    wcs->eqout = wcsceq (coorsys);
    if (wcs->wcson) {

	/* Set output in degrees flag and number of decimal places */
	if (wcs->sysout == WCS_GALACTIC || wcs->sysout == WCS_ECLIPTIC) {
	    wcs->degout = 1;
	    wcs->ndec = 5;
	    }
	else if (wcs->sysout == WCS_ALTAZ) {
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
getwcsout(wcs)
struct	WorldCoor *wcs; /* World coordinate system structure */
{
    return(wcs->radecout);
}


/* Initialize WCS input coordinate system for use by WCS2PIX() */

void
wcsininit (wcs, coorsys)

struct WorldCoor *wcs;	/* World coordinate system structure */
char	*coorsys;	/* Input world coordinate system:
			   FK4, FK5, B1950, J2000, GALACTIC, ECLIPTIC
			   fk4, fk5, b1950, j2000, galactic, ecliptic */
{
    int sysin, i;

    if (nowcs (wcs))
	return;

    /* If argument is null, set to image system and equinox */
    if (coorsys == NULL || strlen (coorsys) < 1) {
	wcs->sysin = wcs->syswcs;
	strcpy (wcs->radecin, wcs->radecsys);
	wcs->eqin = wcs->equinox;
	if (wcs->sysin == WCS_B1950) {
	    if (wcs->eqin != 1950.0) {
		wcs->radecin[0] = 'B';
		sprintf (wcs->radecin+1,"%.4f", wcs->equinox);
		i = strlen(wcs->radecin) - 1;
		if (wcs->radecin[i] == '0')
		    wcs->radecin[i] = (char)0;
		i = strlen(wcs->radecin) - 1;
		if (wcs->radecin[i] == '0')
		    wcs->radecin[i] = (char)0;
		i = strlen(wcs->radecin) - 1;
		if (wcs->radecin[i] == '0')
		    wcs->radecin[i] = (char)0;
		}
	    else
		strcpy (wcs->radecin, "B1950");
	    }
	else if (wcs->sysin == WCS_J2000) {
	    if (wcs->eqin != 2000.0) {
		wcs->radecin[0] = 'J';
		sprintf (wcs->radecin+1,"%.4f", wcs->equinox);
		i = strlen(wcs->radecin) - 1;
		if (wcs->radecin[i] == '0')
		    wcs->radecin[i] = (char)0;
		i = strlen(wcs->radecin) - 1;
		if (wcs->radecin[i] == '0')
		    wcs->radecin[i] = (char)0;
		i = strlen(wcs->radecin) - 1;
		if (wcs->radecin[i] == '0')
		    wcs->radecin[i] = (char)0;
		}
	    else
		strcpy (wcs->radecin, "J2000");
	    }
	}

    /* Ignore new coordinate system if it is not supported */
    if ((sysin = wcscsys (coorsys)) < 0)
	return;

    wcs->sysin = sysin;
    wcs->eqin = wcsceq (coorsys);
    strcpy (wcs->radecin, coorsys);
    return;
}


/* Return current value of WCS input coordinate system set by wcsininit */
char *
getwcsin (wcs)
struct	WorldCoor *wcs; /* World coordinate system structure */
{
    if (nowcs (wcs))
	return (NULL);
    else
	return (wcs->radecin);
}


/* Set WCS output in degrees or hh:mm:ss dd:mm:ss, returning old flag value */
int
setdegout(wcs, new)
     struct     WorldCoor *wcs; /* World coordinate system structure */
     int new;
{
    int old;

    if (nowcs (wcs))
	return (0);
    old = wcs->degout;
    wcs->degout = new;
    return(old);
}


/* Return current value of coordinate system */
char *
getradecsys(wcs)
struct     WorldCoor *wcs; /* World coordinate system structure */
{
    if (nowcs (wcs))
	return (NULL);
    else
	return (wcs->radecsys);
}


/* Set output string mode for LINEAR coordinates */

void
setlinmode (wcs, mode)
struct	WorldCoor *wcs; /* World coordinate system structure */
int	mode;		/* mode = 0: x y linear
			   mode = 1: x units x units
			   mode = 2: x y linear units */
{
    wcs->linmode = mode;
    return;
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
	char	rastr[32], decstr[32];
	int	minlength, lunits, lstring;

	if (nowcs (wcs)) {
	    if (lstr > 0)
		wcstring[0] = 0;
	    return(0);
	    }

	pix2wcs (wcs,xpix,ypix,&xpos,&ypos);

	/* Keep ra/longitude within range
	if (xpos < 0.0)
	    xpos = xpos + 360.0;

	else if (xpos > 360.0)
	    xpos = xpos - 360.0; */

	/* If point is off scale, set string accordingly */
	if (wcs->offscl) {
	    (void)sprintf (wcstring,"Off map");
	    return (1);
	    }

	/* Print coordinates in degrees */
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
		if (wcs->tabsys)
		    strncpy (wcstring,"*********	**********",lstr);
		else
		    strncpy (wcstring,"*******************",lstr);
		lstr = 0;
		}
	    }

	/* print coordinates in sexagesimal notation */
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
		if (wcs->tabsys)
		    strncpy (wcstring,"*************	*************",lstr);
		else
		    strncpy (wcstring,"**************************",lstr);
		lstr = 0;
		}
	    }

	/* Label galactic coordinates */
	if (wcs->sysout == WCS_GALACTIC) {
	    if (lstr > 9 && wcs->printsys)
		strcat (wcstring," galactic");
	    }

	/* Label ecliptic coordinates */
	else if (wcs->sysout == WCS_ECLIPTIC) {
	    if (lstr > 9 && wcs->printsys)
		if (wcs->tabsys)
		    strcat (wcstring,"	ecliptic");
		else
		    strcat (wcstring," ecliptic");
	    }

	/* Label alt-az coordinates */
	else if (wcs->sysout == WCS_ALTAZ) {
	    if (lstr > 7 && wcs->printsys)
		if (wcs->tabsys)
		    strcat (wcstring,"	alt-az");
		else
		    strcat (wcstring," alt-az");
	    }

	/* Label equatorial coordinates */
	else if (wcs->sysout==WCS_B1950 || wcs->sysout==WCS_J2000) {
	    if (lstr > strlen(wcs->radecout)+1 && wcs->printsys) {
		if (wcs->tabsys)
		    strcat (wcstring,"	");
		else
		    strcat (wcstring," ");
		strcat (wcstring, wcs->radecout);
		}
	    }

	/* Output linear coordinates */
	else {
	    num2str (rastr, xpos, 0, wcs->ndec);
	    num2str (decstr, ypos, 0, wcs->ndec);
	    lstring = strlen (rastr) + strlen (decstr) + 1;
	    lunits = strlen (wcs->units[0]) + strlen (wcs->units[1]) + 2;
	    if (wcs->syswcs == WCS_LINEAR && wcs->linmode == 1) {
		if (lstr > lstring + lunits) {
		    if (strlen (wcs->units[0]) > 0) {
			strcat (rastr, " ");
			strcat (rastr, wcs->units[0]);
			}
		    if (strlen (wcs->units[1]) > 0) {
			strcat (decstr, " ");
			strcat (decstr, wcs->units[1]);
			}
		    lstring = lstring + lunits;
		    }
		}
	    if (lstr > lstring) {
		if (wcs->tabsys)
		    (void)sprintf (wcstring,"%s	%s", rastr, decstr);
		else
		    (void)sprintf (wcstring,"%s %s", rastr, decstr);
		}
	    else {
		if (wcs->tabsys)
		    strncpy (wcstring,"**********	*********",lstr);
		else
		    strncpy (wcstring,"*******************",lstr);
		}
	    if (wcs->syswcs == WCS_LINEAR && wcs->linmode != 1 &&
		lstr > lstring + 7)
		strcat (wcstring, " linear");
	    if (wcs->syswcs == WCS_LINEAR && wcs->linmode == 2 &&
		lstr > lstring + lunits + 7) {
		if (strlen (wcs->units[0]) > 0) {
		    strcat (wcstring, " ");
		    strcat (wcstring, wcs->units[0]);
		    }
		if (strlen (wcs->units[1]) > 0) {
		    strcat (wcstring, " ");
		    strcat (wcstring, wcs->units[1]);
		    }
		    
		}
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
    double	eqin, eqout;
    int wcspos();
    extern int dsspos();
    extern int platepos();
    extern int worldpos();
    extern int tnxpos();
    extern void fk4prec(),fk5prec(), wcscon();

    if (nowcs (wcs))
	return;
    wcs->xpix = xpix;
    wcs->ypix = ypix;
    wcs->offscl = 0;

    /* Convert image coordinates to sky coordinates */

/* Use Digitized Sky Survey plate fit */
    if (wcs->prjcode == WCS_DSS) {
	if (dsspos (xpix, ypix, wcs, &xp, &yp))
	    wcs->offscl = 1;
	}

    /* Use SAO plate fit */
    else if (wcs->prjcode == WCS_PLT) {
	if (platepos (xpix, ypix, wcs, &xp, &yp))
	    wcs->offscl = 1;
	}

    /* Use NOAO IRAF corrected plane tangent projection */
    else if (wcs->prjcode == WCS_TNX) {
	if (tnxpos (xpix, ypix, wcs, &xp, &yp))
	    wcs->offscl = 1;
	}

    /* Use Classic AIPS projections */
    else if (wcs->oldwcs == 1 || wcs->prjcode <= 0) {
	if (worldpos (xpix, ypix, wcs, &xp, &yp))
	    wcs->offscl = 1;
	}

    /* Use Mark Calabretta's WCSLIB projections */
    else if (wcspos (xpix, ypix, wcs, &xp, &yp))
	wcs->offscl = 1;
	    	

    /* Do not change coordinates if offscale */
    if (!wcs->offscl) {

	/* Convert coordinates to output system, if not LINEAR */
        if (wcs->prjcode > 0) {

	    /* Convert coordinates to desired output system */
	    eqin = wcs->equinox;
	    eqout = wcs->eqout;
	    wcscon (wcs->syswcs,wcs->sysout,eqin,eqout,&xp,&yp,wcs->epoch);
	    }
	wcs->xpos = xp;
	wcs->ypos = yp;
	*xpos = xp;
	*ypos = yp;
	}
    else {
	*xpos = 0.0;
	*ypos = 0.0;
	}
    return;
}


/* Convert World Coordinates to pixel coordinates */

void
wcs2pix (wcs, xpos, ypos, xpix, ypix, offscl)

struct WorldCoor *wcs;	/* World coordinate system structure */
double	xpos,ypos;	/* World coordinates in degrees */
double	*xpix,*ypix;	/* Image coordinates in pixels */
int	*offscl;	/* 0 if within bounds, else off scale */
{
    wcsc2pix (wcs, xpos, ypos, wcs->radecin, xpix, ypix, offscl);
    return;
}

/* Convert World Coordinates to pixel coordinates */

void
wcsc2pix (wcs, xpos, ypos, coorsys, xpix, ypix, offscl)

struct WorldCoor *wcs;	/* World coordinate system structure */
double	xpos,ypos;	/* World coordinates in degrees */
char	*coorsys;	/* Input world coordinate system:
			   FK4, FK5, B1950, J2000, GALACTIC, ECLIPTIC
			   fk4, fk5, b1950, j2000, galactic, ecliptic */
double	*xpix,*ypix;	/* Image coordinates in pixels */
int	*offscl;	/* 0 if within bounds, else off scale */
{
    double xp,yp;
    double eqin, eqout;
    int sysin;
    int wcspix();
    extern int dsspix();
    extern int platepix();
    extern int worldpix();
    extern int tnxpix();
    extern void fk4prec(), fk5prec(), wcscon();

    if (nowcs (wcs))
	return;

    *offscl = 0;
    xp = xpos;
    yp = ypos;
    sysin = wcscsys (coorsys);
    eqin = wcsceq (coorsys);

    /* Convert coordinates to same system as image */
    eqout = wcs->equinox;
    wcscon (sysin, wcs->syswcs, eqin, eqout, &xp, &yp, wcs->epoch);

    /* Convert sky coordinates to image coordinates */

    /* Use Digitized Sky Survey plate fit */
    if (wcs->prjcode == WCS_DSS) {
	if (dsspix (xp, yp, wcs, xpix, ypix))
	    *offscl = 1;
	}

    /* Use SAO polynomial plate fit */
    else if (wcs->prjcode == WCS_PLT) {
	if (platepix (xp, yp, wcs, xpix, ypix))
	    wcs->offscl = 1;
	}

    /* Use NOAO IRAF corrected plane tangent projection */
    else if (wcs->prjcode == WCS_TNX) {
	if (tnxpix (xp, yp, wcs, xpix, ypix))
	    wcs->offscl = 1;
	}

    /* Use Classic AIPS projections */
    else if (wcs->oldwcs == 1 || wcs->prjcode <= 0) {
	if (worldpix (xp, yp, wcs, xpix, ypix))
	    *offscl = 1;
	}

    /* Use Mark Calabretta's WCSLIB projections */
    else if (wcspix (xp, yp, wcs, xpix, ypix)) {
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


int
wcspos (xpix, ypix, wcs, xpos, ypos)

/* Input: */
double  xpix;          /* x pixel number  (RA or long without rotation) */
double  ypix;          /* y pixel number  (dec or lat without rotation) */
struct WorldCoor *wcs;  /* WCS parameter structure */

/* Output: */
double  *xpos;           /* x (RA) coordinate (deg) */
double  *ypos;           /* y (dec) coordinate (deg) */
{
    int offscl;
    int wcsrev();
    double wcscrd[4], imgcrd[4], pixcrd[4];
    double phi, theta;
    
    *xpos = 0.0;
    *ypos = 0.0;

    pixcrd[0] = xpix;
    pixcrd[1] = ypix;
    pixcrd[2] = 1.0;
    pixcrd[3] = 1.0;
    offscl = wcsrev (wcs->ctype, &wcs->wcsl, pixcrd, &wcs->lin, imgcrd,
		    &wcs->prj, &phi, &theta, wcs->crval, &wcs->cel, wcscrd);
    if (offscl == 0) {
	*xpos = wcscrd[wcs->wcsl.lng];
	*ypos = wcscrd[wcs->wcsl.lat];
	}
    return (offscl);
}

int
wcspix (xpos, ypos, wcs, xpix, ypix)

/* Input: */
double  xpos;           /* x (RA) coordinate (deg) */
double  ypos;           /* y (dec) coordinate (deg) */
struct WorldCoor *wcs;  /* WCS parameter structure */

/* Output: */
double  *xpix;          /* x pixel number  (RA or long without rotation) */
double  *ypix;          /* y pixel number  (dec or lat without rotation) */

{
    int offscl;
    int wcsfwd();
    double wcscrd[4], imgcrd[4], pixcrd[4];
    double phi, theta;

    *xpix = -1.0;
    *ypix = -1.0;
    if (wcs->wcsl.flag != WCSSET) {
	if (wcsset (wcs->lin.naxis, wcs->ctype, &wcs->wcsl) )
	    return (1);
	}
    wcscrd[wcs->wcsl.lng] = xpos;
    wcscrd[wcs->wcsl.lat] = ypos;
    offscl = wcsfwd (wcs->ctype, &wcs->wcsl, wcscrd, wcs->crval, &wcs->cel,
		     &phi, &theta, &wcs->prj, imgcrd, &wcs->lin, pixcrd);
    if (offscl == 0) {
	*xpix = pixcrd[0];
	*ypix = pixcrd[1];
	}
    return (offscl);
}

void
setdefwcs (oldwcs)
int oldwcs;
{ oldwcs0 = oldwcs; return;}

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
 * May 22 1997	Add PIXEL prjcode = -1;
 * Jul 11 1997	Get center pixel x and y from header even if no WCS
 * Aug  7 1997	Add NOAO PIXSCALi keywords for default WCS
 * Oct 15 1997	Do not reset reference pixel in WCSSHIFT
 * Oct 20 1997	Set chip rotation
 * Oct 24 1997	Keep longitudes between 0 and 360, not -180 and +180
 * Nov  5 1997	Do no compute crot and srot; they are now computed in worldpos
 * Dec 16 1997	Set rotation and axis increments from CD matrix
 *
 * Jan  6 1998	Deal with J2000 and B1950 as EQUINOX values (from ST)
 * Jan  7 1998	Read INSTRUME and DETECTOR header keywords
 * Jan  7 1998	Fix tab-separated output
 * Jan  9 1998	Precess coordinates if either FITS projection or *DSS plate*
 * Jan 16 1998	Change PTYPE to not include initial hyphen
 * Jan 16 1998	Change WCSSET to WCSXINIT to avoid conflict with Calabretta
 * Jan 23 1998	Change PCODE to PRJCODE to avoid conflict with Calabretta
 * Jan 27 1998	Add LATPOLE and LONGPOLE for Calabretta projections
 * Feb  5 1998	Make cd and dc into vectors; use matinv() to invert cd
 * Feb  5 1998	In wcsoutinit(), check that corsys is a valid pointer
 * Feb 18 1998	Fix bugs for Calabretta projections
 * Feb 19 1998	Add wcs structure access subroutines from Eric Mandel
 * Feb 19 1998	Add wcsreset() to make sure derived values are reset
 * Feb 19 1998	Always set oldwcs to 1 if NCP projection
 * Feb 19 1998	Add subroutine to set oldwcs default
 * Feb 20 1998	Initialize projection types one at a time for SunOS C
 * Feb 23 1998	Add TNX projection from NOAO; treat it as TAN
 * Feb 23 1998	Compute size based on max and min coordinates, not sides
 * Feb 26 1998	Add code to set center pixel if part of detector array
 * Mar  6 1998	Write 8-character values to RADECSYS
 * Mar  9 1998	Add naxis to WCS structure
 * Mar 16 1998	Use separate subroutine for IRAF TNX projection
 * Mar 20 1998	Set PC matrix if more than two axes and it's not in header
 * Mar 20 1998	Reset lin flag in WCSRESET if CDELTn
 * Mar 20 1998	Set CD matrix with CDELTs and CROTA in wcsinit and wcsreset
 * Mar 20 1998	Allow initialization of rotation angle alone
 * Mar 23 1998	Use dsspos() and dsspix() instead of platepos() and platepix()
 * Mar 24 1998	Set up PLT/PLATE plate polynomial fit using platepos() and platepix()
 * Mar 25 1998	Read plate fit coefficients from header in getwcscoeffs()
 * Mar 27 1998	Check for FITS WCS before DSS WCS
 * Mar 27 1998	Compute scale from edges if xinc and yinc not set in wcscent()
 * Apr  6 1998	Change plate coefficient keywords from PLTij to COi_j
 * Apr  6 1998	Read all coefficients in line instead of with subroutine
 * Apr  7 1998	Change amd_i_coeff to i_coeff
 * Apr  8 1998	Add wcseqset to change equinox after wcs has been set
 * Apr 10 1998	Use separate counters for x and y polynomial coefficients
 * Apr 13 1998	Use CD/CDELT+CDROTA if oldwcs is set
 * Apr 14 1998	Use codes instead of strings for various coordinate systems
 * Apr 14 1998	Separate input coordinate conversion from output conversion
 * Apr 14 1998	Use wcscon() for most coordinate conversion
 * Apr 17 1998	Always compute cdelt[]
 * Apr 17 1998	Deal with reversed axis more correctly
 * Apr 17 1998	Compute rotation angle and approximate CDELTn for polynomial
 * Apr 23 1998	Deprecate xref/yref in favor of crval[]
 * Apr 23 1998	Deprecate xrefpix/yrefpix in favor of crpix[]
 * Apr 23 1998	Add LINEAR to coordinate system types
 * Apr 23 1998	Always use AIPS subroutines for LINEAR or PIXEL
 * Apr 24 1998	Format linear coordinates better
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * Apr 28 1998  Change projection flags to WCS_*
 * Apr 28 1998	Add subroutine wcsc2pix for coordinates to pixels with system
 * Apr 28 1998	Add setlinmode() to set output string mode for LINEAR coordinates
 * Apr 30 1998	Fix bug by setting degree flag for lat and long in wcsinit()
 * Apr 30 1998	Allow leading "-"s in projecting in wcsxinit()
 * May  1 1998	Assign new output coordinate system only if legitimate system
 * May  1 1998	Do not allow oldwcs=1 unless allowed projection
 * May  4 1998	Fix bug in units reading for LINEAR coordinates
 * May  6 1998	Initialize to no CD matrix
 * May  6 1998	Use TAN instead of TNX if oldwcs
 * May 12 1998	Set 3rd and 4th coordinates in wcspos()
 * May 12 1998	Return *xpos and *ypos = 0 in pix2wcs() if offscale
 * May 12 1998	Declare undeclared external subroutines after lint
 * May 13 1998	Add equinox conversion to specified output equinox
 * May 13 1998	Set output or input system to image with null argument
 * May 15 1998	Return reference pixel, cdelts, and rotation for DSS
 */
