

/* Set parameter values from the command line as keyword=value
 * Return 1 if successful, else 0 */

static int
scatparm (parstring)

char *parstring;
{
    char *parname;
    char *parvalue;
    char *parequal;
    char *temp;
    char *refcatn;
    int lcat, lrange;

    /* Check for request for help */
    if (!strncasecmp (parstring,"help", 4)) {
	PrintWebHelp ();
	}

    /* Check for request for command line help */
    if (!strncasecmp (parstring,"comhelp", 7)) {
	PrintUsage (NULL);
	}

    /* Check for request for GSC bandpasses */
    if (!strncasecmp (parstring, "band", 4) ||
	!strncasecmp (parstring, "filt", 4)) {
	PrintGSCBand ();
	}

    /* Check for request for GSC object classes */
    if (!strncasecmp (parstring, "clas", 4) ||
	!strncasecmp (parstring, "obje", 4)) {
	PrintGSClass ();
	}

    /* Check for scat version request */
    if (!strcasecmp (parstring, "version")) {
	PrintUsage ("version");
	}

    /* Separate parameter name and value */
    parname = parstring;
    if ((parequal = strchr (parname,'=')) == NULL)
	return (0);
    *parequal = (char) 0;
    parvalue = parequal + 1;

    /* Get closest source */
    if (!strcasecmp (parname, "closest")) {
	if (!strncasecmp (parvalue, "y", 1)) {
	    catsort = SORT_DIST;
	    nstars = 1;
	    closest++;
	    }
	}

    /* Set range of source numbers to print */
    else if (!strncasecmp (parname,"num",3)) {
	if (ranges) {
	    temp = ranges;
	    lrange = strlen(ranges) + strlen(parvalue) + 2;
	    ranges = (char *) malloc (lrange);
	    strcpy (ranges, temp);
	    strcat (ranges, ",");
	    strcat (ranges, parvalue);
	    free (temp);
	    }
	else {
	    lrange = strlen (parvalue) + 2;
	    ranges = (char *) malloc (lrange);
	    if (strchr (parvalue,'.'))
		match = 1;
	    strcpy (ranges, parvalue);
	    }
	}

    /* Radius in arcseconds */
    else if (!strncasecmp (parname,"rad",3)) {
	if (strchr (parvalue,':'))
	    rad0 = 3600.0 * str2dec (parvalue);
	else
	    rad0 = atof (parvalue);
	}

    /* Radius in degrees */
    else if (!strncasecmp (parname,"sr",3)) {
	if (strchr (parvalue,':'))
	    rad0 = 3600.0 * str2dec (parvalue);
	else
	    rad0 = 3600.0 * atof (parvalue);
	votab = 1;
	degout0 = 1;
	}

    /* Search center right ascension */
    else if (!strcasecmp (parname,"ra"))
	ra0 = str2ra (parvalue);

    /* Search center declination */
    else if (!strcasecmp (parname,"dec"))
	dec0 = str2dec (parvalue);

    /* Search center coordinate system */
    else if (!strncasecmp (parname,"sys",3)) {
	syscoor = wcscsys (parvalue);
	eqcoor = wcsceq (parvalue);
	}

    /* Output coordinate system */
    else if (!strcasecmp (parname, "outsys")) {

	/* B1950 (FK4) coordinates */
	if (!strcasecmp (parvalue, "B1950") ||
	    !strcasecmp (parvalue, "FK4")) {
	    sysout0 = WCS_B1950;
	    eqout = 1950.0;
    	    }

	/* J2000 (FK5) coordinates */
	else if (!strcasecmp (parvalue, "J2000") ||
	    !strcasecmp (parvalue, "FK5")) {
	    sysout0 = WCS_J2000;
	    eqout = 2000.0;
    	    }

	/* Galactic coordinates */
	else if (!strncasecmp (parvalue, "GAL", 3))
	    sysout0 = WCS_GALACTIC;

	/* Ecliptic coordinates */
	else if (!strncasecmp (parvalue, "ECL", 3))
	    sysout0 = WCS_ECLIPTIC;
	}

    /* Set reference catalog */
    else if (!strcasecmp (parname, "catalog")) {
	lcat = strlen (parvalue) + 2;
	refcatn = (char *) malloc (lcat);
	strcpy (refcatn, parvalue);
	refcatname[ncat] = refcatn;
	ncat = ncat + 1;
	}

    /* Set output coordinate epoch */
    else if (!strcasecmp (parname, "epoch"))
	epoch0 = fd2ep (parvalue);

    /* Output equinox in years */
    else if (!strcasecmp (parname, "equinox"))
    	eqout = fd2ep (parvalue);

    /* Output in degrees instead of sexagesimal */
    else if (!strcasecmp (parname, "cformat")) {
	if (!strncasecmp (parvalue, "deg", 3))
	    degout0 = 1;
	else if (!strncasecmp (parvalue, "rad", 3))
	    degout0 = 2;
	else
	    degout0 = 0;
	}

    /* Number of decimal places in output positions */
    else if (!strcasecmp (parname, "ndec")) {
	if (isnum (parvalue)) {
	    if (degout0) {
		nddeg = atoi (parvalue);
		if (nddeg < 0 || nddeg > 10)
		    nddeg = 7;
		}
	    else {
		ndra = atoi (parvalue);
		if (ndra < 0 || ndra > 10)
		    ndra = 3;
		nddec = ndra - 1;
		}
	    }
	}

    /* Minimum proper motion quality for USNO-B1.0 catalog */
    else if (!strcasecmp (parname, "minpmq")) {
	if (isnum (parvalue)) {
	    minid = atoi (parvalue);
	    setminpmqual (minid);
	    }
	}

    /* Minimum number of plate ID's for USNO-B1.0 catalog */
    else if (!strcasecmp (parname, "minid")) {
	if (isnum (parvalue)) {
	    minid = atoi (parvalue);
	    setminid (minid);
	    }
	}

    /* Output in VOTable XML instead of tab-separated table */
    else if (!strcasecmp (parname, "format")) {
	if (!strncasecmp (parvalue, "vot", 3)) {
	    votab = 1;
	    degout0 = 1;
	    }
	else
	    votab = 0;
	}

    /* Print center and closest star on one line */
    else if (!strcasecmp (parname, "oneline")) {
	if (parvalue[0] == 'y' || parvalue[0] == 'Y') {
	    oneline++;
	    catsort = SORT_DIST;
	    closest++;
	    nstars = 1;
	    }
	}

    /* Magnitude limit */
    else if (!strncasecmp (parname,"mag",3)) {
	maglim2 = atof (parvalue);
	if (MAGLIM1 == MAGLIM2)
	    maglim1 = -2.0;
	}
    else if (!strncasecmp (parname,"max",3))
	maglim2 = atof (parvalue);
    else if (!strncasecmp (parname,"min",3))
	maglim1 = atof (parvalue);

    /* Number of brightest stars to read */
    else if (!strncasecmp (parname,"nstar",5))
	nstars = atoi (parvalue);

    /* Object name */
    else if (!strcmp (parname, "object") || !strcmp (parname, "OBJECT")) {
	lcat = strlen (parvalue) + 2;
	objname = (char *) malloc (lcat);
	strcpy (objname, parvalue);
	}

    /* Magnitude by which to sort */
    else if (!strcasecmp (parname, "sortmag")) {
	if (isnum (parvalue+1))
	    sortmag = atoi (parvalue);
	}

    /* Output sorting */
    else if (!strcasecmp (parname, "sort")) {

	/* Sort by distance from center */
	if (!strncasecmp (parvalue,"di",2))
	    catsort = SORT_DIST;

	/* Sort by RA */
	if (!strncasecmp (parvalue,"r",1))
	    catsort = SORT_RA;

	/* Sort by Dec */
	if (!strncasecmp (parvalue,"de",2))
		catsort = SORT_DEC;

	/* Sort by magnitude */
	if (!strncasecmp (parvalue,"m",1)) {
	    catsort = SORT_MAG;
	    if (strlen (parvalue) > 1) {
		if (isnum (parvalue+1))
		    sortmag = atoi (parvalue+1);
		}
	    }

	/* No sort */
	if (!strncasecmp (parvalue,"n",1))
	    catsort = NOSORT;
	}

    /* Search box half-width in RA */
    else if (!strcasecmp (parname,"dra")) {
	if (strchr (parvalue,':'))
	    dra0 = 3600.0 * str2ra (parvalue);
	else
	    dra0 = atof (parvalue);
	if (ddec0 <= 0.0)
	    ddec0 = dra0;
	/* rad0 = sqrt (dra0*dra0 + ddec0*ddec0); */
	}

    /* Search box half-height in Dec */
    else if (!strcasecmp (parname,"ddec")) {
	if (strchr (parvalue,':'))
	    ddec0 = 3600.0 * str2dec (parvalue);
	else
	    ddec0 = atof (parvalue);
	if (dra0 <= 0.0)
	    dra0 = ddec0;
	rad0 = sqrt (dra0*dra0 + ddec0*ddec0);
	}

    /* Guide Star object class */
    else if (!strncasecmp (parname,"cla",3)) {
	classd = (int) atof (parvalue);
	setgsclass (classd);
	}

    /* Catalog to be searched */
    else if (!strncasecmp (parname,"cat",3)) {
	lcat = strlen (parvalue) + 2;
	refcatn = (char *) malloc (lcat);
	strcpy (refcatn, parvalue);
	refcatname[ncat] = refcatn;
	ncat = ncat + 1;
	}
    else {
	*parequal = '=';
	return (1);
	}
    *parequal = '=';
    return (0);
}
