

int
uacrnum (nstars,unum,ura,udec,umag,umagb,uplate,nlog)

int	nstars;		/* Number of stars to find */
double	*unum;		/* Array of UA numbers to find */
double	*ura;		/* Array of right ascensions (returned) */
double	*udec;		/* Array of declinations (returned) */
double	*umag;		/* Array of red magnitudes (returned) */
double	*umagb;		/* Array of blue magnitudes (returned) */
int	*uplate;	/* Array of plate numbers (returned) */
int	nlog;		/* Logging interval */
{
    UACstar star;	/* UA catalog entry for one star */

    int verbose;
    int znum, itot,iz, i;
    int itable, jtable,jstar;
    int nstar, nread, nzone;
    int nfound = 0;
    double ra,dec;
    double mag, magb;
    int istar, istar1, istar2, plate;
    int nzmax = NZONES;	/* Maximum number of declination zones */
    char *str;

    itot = 0;
    if (nlog > 0)
	verbose = 1;
    else
	verbose = 0;

    /* Set path to USNO A Catalog */
    if ((str = getenv("UA_PATH")) != NULL )
	strcpy (cdu,str);

/* Loop through star list */
    for (jstar = 0; jstar < nstars; jstar++) {

    /* Get path to zone catalog */
	znum = (int) unum[jstar];
	if ((nzone = uacopen (znum)) != 0) {

	    istar = (int) (((unum[jstar] - znum) * 100000000.0) + 0.5);

	    if (uacstar (istar, &star)) {
		fprintf (stderr,"UACRNUM: Cannot read star %d\n", istar);
		break;
		}

	    /* Extract selected fields if not probable duplicate */
	    else if (star.magetc > 0) {
		ra = uacra (star.rasec); /* Right ascension in degrees */
		dec = uacdec (star.decsec); /* Declination in degrees */
		magb = uacmagb (star.magetc); /* Blue magnitude */
		mag = uacmagr (star.magetc); /* Red magnitude */
		plate = uacplate (star.magetc);	/* Plate number */

		/* Save star position and magnitude in table */
		ura[nfound] = ra;
		udec[nfound] = dec;
		umag[nfound] = mag;
		umagb[nfound] = magb;
		uplate[nfound] = plate;

		nfound++;
		if (nlog == 1)
		    fprintf (stderr,"UACRNUM: %04d.%08d: %9.5f %9.5f %5.2f\n",
			     znum,istar,ra,dec,mag);

		/* Log operation */
		if (nlog > 0 && jstar%nlog == 0)
		    fprintf (stderr,"UACRNUM: %4d.%8d  %8d / %8d sources\r",
			     znum, istar, jstar, nstar);

		(void) fclose (fcat);
		/* End of star processing */
		}

	    /* End of star */
	    }

	/* End of star loop */
	}

    /* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"UACRNUM:  %d / %d found\n",nfound,nstars);

    return (nfound);
}


