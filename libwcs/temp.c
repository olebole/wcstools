void
CatTabHead (refcat,sysout,mprop,ranges,keyword,gcset,tabout,printxy,gobj1,fd)

int	refcat;		/* Catalog being searched */
int	sysout;		/* Output coordinate system */
int	mprop;		/* 1 if proper motion in catalog */
char	*ranges;	/* Catalog numbers to print */
char	*keyword;	/* Column to add to tab table output */
int	gcset;		/* 1 if there are any values in gc[] */
int	tabout;		/* 1 if output is tab-delimited */
int	printxy;	/* 1 if X and Y included in output */
char	**gobj1;	/* Pointer to array of object names; NULL if none */
FILE	*fd;		/* Output file descriptor; none if NULL */

{    int typecol;

    /* Set flag for plate, class, type, or 3rd magnitude column */
    if (refcat == BINCAT || refcat == SAO  || refcat == PPM ||
	refcat == ACT  || refcat == TYCHO2 || refcat == BSC)
	typecol = 1;
    else if ((refcat == GSC || refcat == GSCACT) && classd < -1)
	typecol = 3;
    else if (refcat == TMPSC)
	typecol = 4;
    else if (refcat == GSC || refcat == GSCACT ||
	refcat == UJC || refcat == IRAS ||
	refcat == USAC || refcat == USA1   || refcat == USA2 ||
	refcat == UAC  || refcat == UA1    || refcat == UA2 ||
	refcat == BSC  || (refcat == TABCAT&&gcset))
	typecol = 2;
    else
	typecol = 0;


    /* Print column headings */
    if (refcat == ACT)
	strcpy (headline, "act_id       ");
    else if (refcat == BSC)
	strcpy (headline, "bsc_id       ");
    else if (refcat == GSC || refcat == GSCACT)
	strcpy (headline, "gsc_id       ");
    else if (refcat == USAC)
	strcpy (headline,"usac_id       ");
    else if (refcat == USA1)
	strcpy (headline,"usa1_id       ");
    else if (refcat == USA2)
	strcpy (headline,"usa2_id       ");
    else if (refcat == UAC)
	strcpy (headline,"usnoa_id      ");
    else if (refcat == UA1)
	strcpy (headline,"usnoa1_id     ");
    else if (refcat == UA2)
	strcpy (headline,"usnoa2_id     ");
    else if (refcat == UJC)
	strcpy (headline,"usnoj_id      ");
    else if (refcat == TMPSC)
	strcpy (headline,"2mass_id      ");
    else if (refcat == SAO)
	strcpy (headline,"sao_id        ");
    else if (refcat == PPM)
	strcpy (headline,"ppm_id        ");
    else if (refcat == IRAS)
	strcpy (headline,"iras_id       ");
    else if (refcat == TYCHO)
	strcpy (headline,"tycho_id      ");
    else if (refcat == TYCHO2)
	strcpy (headline,"tycho2_id     ");
    else if (refcat == HIP)
	strcpy (headline,"hip_id        ");
    else
	strcpy (headline,"id            ");
    headline[nnfld] = (char) 0;

    if (sysout == WCS_GALACTIC)
	strcat (headline,"	long_gal   	lat_gal  ");
    else if (sysout == WCS_ECLIPTIC)
	strcat (headline,"	long_ecl   	lat_ecl  ");
    else if (sysout == WCS_B1950)
	strcat (headline,"	ra1950      	dec1950  ");
    else
	strcat (headline,"	ra      	dec      ");
    if (refcat == USAC || refcat == USA1 || refcat == USA2 ||
	refcat == UAC  || refcat == UA1  || refcat == UA2)
	strcat (headline,"	magb	magr	plate");
    if (refcat == TMPSC)
	strcat (headline,"	magj	magh	magk");
    else if (refcat==TYCHO || refcat==TYCHO2 || refcat==HIP || refcat==ACT)
	strcat (headline,"	magb	magv");
    else if (refcat == GSC || refcat == GSCACT)
	strcat (headline,"	mag	class	band	N");
    else if (refcat == UJC)
	strcat (headline,"	mag	plate");
    else
	strcat (headline,"	mag");
    if (typecol == 1)
	strcat (headline,"	type");
    if (mprop)
	strcat (headline,"	Ura    	Udec  ");
    if (ranges == NULL)
	strcat (headline,"	arcsec");
    if (refcat == TABCAT && keyword != NULL) {
	strcat (headline,"	");
	strcat (headline, keyword);
	}
    if (gobj1 != NULL)
	strcat (headline,"	object");
    if (printxy)
	strcat (headline, "	X      	Y      ");
    if (tabout) {
	printf ("%s\n", headline);
	if (fd != NULL)
	    fprintf (fd, "%s\n", headline);
	}

    strcpy (headline, "---------------------");
    headline[nnfld] = (char) 0;
    strcat (headline,"	------------	------------");
    if (nmag == 2)
	strcat (headline,"	-----	-----");
    else
	strcat (headline,"	-----");
    if (refcat == GSC || refcat == GSCACT)
	strcat (headline,"	-----	----	-");
    else if (typecol == 1)
	strcat (headline,"	----");
    else if (typecol == 2)
	strcat (headline,"	-----");
    else if (typecol == 4)
	strcat (headline,"	-----");
    if (mprop)
	strcat (headline,"	-------	------");
    if (ranges == NULL)
	strcat (headline, "	------");
    if (refcat == TABCAT && keyword != NULL)
	strcat (headline,"	------");
    if (printxy)
	strcat (headline, "	-------	-------");
    if (tabout) {
	printf ("%s\n", headline);
	if (fd != NULL)
	    fprintf (fd, "%s\n", headline);
	}
    }
