/*** File 
 */

int
platefwd (xpos, ypos, wcs, xpix, ypix)

    double tdec,ctan,ccos,traoff, craoff, etar, xir;
    double ra0, dec0, ra, dec;

    /* Convert RA and Dec in radians to standard coordinates on a plate */
    ra = degrad (xpos);
    dec = degrad (ypos);
    tdec = tan (dec);
    ra0 = degrad (wcs->crval[0]);
    dec0 = degrad (wcs->crval[1]);
    ctan = tan (dec0);
    ccos = cos (dec0);
    traoff = tan (ra - ra0);
    craoff = cos (ra - ra0);
    etar = (1.0 - ctan * craoff / tdec) / (ctan + (craoff / tdec));
    xir = traoff * ccos * (1.0 - (etar * ctan));
    xi = raddeg (xir);
    eta = raddeg (etar);

    polyfwd (xi, eta, wcs, xpix, ypix);

    /* Convert from plate pixels to image pixels */
    *xpix = x + wcs->crval[0];
    *ypix = y + wcs->crval[1];

    /* If position is off of the image, return offscale code */
    if (*xpix < 0.5 || *xpix > wcs->nxpix+0.5)
	return -1;
    if (*ypix < 0.5 || *ypix > wcs->nypix+0.5)
	return -1;
