/*** File libwcs/gscread.c
 *** February 16, 1996
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

#include <stdlib.h>
#include <string.h>
#include "fitshead.h"

/* GSCREAD -- Read HST Guide Star Catalog stars from CDROM */

int gscread (ra1,ra2,dec1,dec2,mag1,mag2,classd,nstarmax,gnum,gra,gdec,gmag)

double	ra1,ra2;	/* Limiting right ascensions of region in degrees */
double dec1,dec2;	/* Limiting declinations of region in degrees */
double mag1,mag2;	/* Limiting magnitudes (none if equal) */
int classd;		/* Desired object class (-1=all, 0=stars, 3=nonstars) */
int nstarmax;		/* Maximum number of stars to be returned */
double *gnum;		/* Array of Guide Star numbers (returned) */
double *gra;		/* Array of right ascensions (returned) */
double *gdec;		/* Array of declinations (returned) */
double *gmag;		/* Array of magnitudes (returned) */
{
    char *header;	/* FITS header */
    char *table;	/* FITS table */
    int nreg;		/* Number of input FITS tables files */
    int rlist[100];	/* List of input FITS tables files */
    char inpath[64];	/* Pathname for input FITS table file */
    char entry[100];	/* Buffer for FITS table row */
    int class, class0;	/* Object class (0>star, 3>other) */
    struct Keyword kw[8];	/* Keyword structure */

    int first;
    int wrap;
    char temp[50];
    int rnum, num0, num, itot,jtot,ireg,ltab,lhead;
    int nlog,ik,nk,ier,itable,ntable,jstar;
    int nbline,npos,nbhead;
    int ift,nbr,nrmax,nstar,i;
    double ra,ra0,rasum,dec,dec0,decsum,perr,perr0,perr2,perrsum,msum;
    double mag,mag0,merr,merr0,merr2,merrsum;
    int fitsrtline();

    itot = 0;
    jtot = 0;
    ltab = 10000;
    table = malloc (10000);
    lhead = 14400;
    header = malloc (14400);
    for (i = 0; i < 100; i++)
	entry[i] = 0;

/* If RA range includes zero, split it in two */
    wrap = 0;
    if (ra1 > ra2)
	wrap = 1;
    else
	wrap = 0;

/* Make dec1 always the smallest declination */
    if (dec1 > dec2) {
	dec = dec1;
	dec1 = dec2;
	dec2 = dec;
	}

/* make mag1 always the smallest magnitude */
    if (mag2 < mag1) {
	mag = mag2;
	mag2 = mag1;
	mag1 = mag;
	}

/* Find Guide Star Catalog regions in which to search */
    nrmax = 100;
    nreg = gscreg (ra1,ra2,dec1,dec2, lhead,header, ltab,table, nrmax,rlist);
    if (nreg <= 0) {
	printf ("GSCREAD:  no Guide Star regions found\n");
	free (table);
	free (header);
	return (0);
	}

    classd = -1;

/* Logging interval */
	nlog = 0;

/* Set keyword list */
    nk = 8;
    strcpy (kw[0].kname,"GSC_ID");
    strcpy (kw[1].kname,"RA_DEG");
    strcpy (kw[2].kname,"DEC_DEG");
    strcpy (kw[3].kname,"POS_ERR");
    strcpy (kw[4].kname,"MAG");
    strcpy (kw[5].kname,"MAG_ERR");
    strcpy (kw[6].kname,"MAG_BAND");
    strcpy (kw[7].kname,"CLASS");
    for (ik = 0; ik < nk; ik++) {
	kw[ik].kn = 0;
	kw[ik].kf = 0;
	kw[ik].kl = 0;
	}
    first = 1;
    nstar = 0;

/* Loop through region list */
    for (ireg = 0; ireg < nreg; ireg++) {
	gscpath (rlist[ireg],inpath);

	/* Read keyword info from FITS table header for first region */
	if (first) {
	    ift = fitsrtopen (inpath,lhead,header,nk,kw,&ntable,&nbline,&nbhead);
	    first = 0;
	    }

	/* Open FITS table for this region */
	else {
	    ift =  fitsropen (inpath, lhead,header, &nbhead);
	    if (ift < 0) {
		printf ("GSCREAD: file %s table not found\n",inpath);
		break;
		}
	    }

	rnum = rlist[ireg];
	num0 = 0;
	rasum = 0.0;
	decsum = 0.0;
	msum = 0.0;
	perrsum = 0.0;
	merrsum = 0.0;
	npos = 0;
	num = 0;
	fitsrtlset();
	jstar = 0;

	/* Loop through FITS table for this region */
	for (itable = 0; itable <= ntable; itable++) {

	    if (itable < ntable) {
		nbr = fitsrtline (ift,nbhead,ltab,table,itable,nbline,entry);
		if (nbr < nbline) {
		    printf ("GSCREAD: %d / %d bytes read, line %d / %d, region %d\n",
			      nbr,nbline,itable,ntable,rnum);
		    break;
		    }

	 /* Extract selected fields */

		/* Star number within region */
		strncpy (temp,entry+kw[0].kf,kw[0].kl);
		temp[kw[0].kl] = 0;
		num0 = atoi (temp);

		/* Right ascension in degrees */
		strncpy (temp,entry+kw[1].kf,kw[1].kl);
		temp[kw[1].kl] = 0;
		ra0 = atof (temp);

		/* Declination in degrees */
		strncpy (temp,entry+kw[2].kf,kw[2].kl);
		temp[kw[2].kl] = 0;
		dec0 = atof (temp);

		/* Position error */
		strncpy (temp,entry+kw[3].kf,kw[3].kl);
		temp[kw[3].kl] = 0;
		perr0 = atof (temp);

		/* Magnitude */
		strncpy (temp,entry+kw[4].kf,kw[4].kl);
		temp[kw[4].kl] = 0;
		mag0 = atof (temp);

		/* Magnitude error */
		strncpy (temp,entry+kw[5].kf,kw[5].kl);
		temp[kw[5].kl] = 0;
		merr0 = atof (temp);

		/* Bandpass code
		strncpy (temp,entry+kw[6].kf,kw[6].kl);
		temp[kw[6].kl] = 0;
		band0 = atoi (temp); */

		/* Object class code */
		strncpy (temp,entry+kw[7].kf,kw[7].kl);
		temp[kw[7].kl] = 0;
		class0 = atoi (temp);
		}
	    else
		num0 = 0;

	/* Compute mean position and magnitude for object */
	    if (num != num0 && itable > 0) {
		if (npos <= 0) continue;
		ra = rasum / perrsum;
		dec = decsum / perrsum;
		mag = msum / merrsum;

	/* Check magnitude amd position limits */
		if ((mag1 != mag2 && (mag >= mag1 && mag <= mag2)) &&
		    ((wrap && (ra <= ra1 || ra >= ra2)) ||
		    (!wrap && (ra >= ra1 && ra <= ra2))) &&
     		    (dec >= dec1 && dec <= dec2)) {

	/* Save star position in table */
		    if (nstar <= nstarmax) {
			gnum[nstar] = (double)rnum + (0.0001 * (double) num);
			gra[nstar] = ra;
			gdec[nstar] = dec;
			gmag[nstar] = mag;
			}
		    nstar = nstar + 1;
		    jstar = jstar + 1;
		    if (nlog == 1)
			printf ("%04d.%04d: %9.5f %9.5f %5.2f %d %d\n",
				rnum,num,ra,dec,mag,class,npos);
		    }

	/* Reset star position for averaging */
		rasum = 0.0;
		decsum = 0.0;
		msum = 0.0;
		perrsum = 0.0;
		merrsum = 0.0;
		npos = 0;
		}

	/* Add information from current line to current object */

	    /* Check object class */
     	    if ((classd > 0 && class0 == classd) ||
		(classd >= -1 && class0 != 5)) {
		perr = perr0;
		perr2 = perr * perr;
		if (perr2 <= 0.0) perr2 = 0.01;
		rasum = rasum + (ra0 / perr2);
		decsum = decsum + (dec0 / perr2);
		perrsum = perrsum + (1.0 / perr2);
		if (merr0 <= 0.0) merr0 = 0.01;
		merr = merr0;
		merr2 = merr * merr;
		msum = msum + (mag0 / merr2);
		merrsum = merrsum + (1.0 / merr2);
		num = num0;
		class = class0;
		npos = npos + 1;
		}

/* Log operation */

	    if (nlog > 0 && itable%nlog == 0)
		printf ("%4d / %4d: %6d / %6d sources from %s\r",
			 ireg,nreg,jstar,itable,inpath);

/* End of region */
	    }

/* Close region input file */
	ier = close (ift);
	itot = itot + itable;
	if (nlog > 0)
	    printf ("%4d / %4d: %6d / %6d sources from %s\n",
		     ireg+1,nreg,jstar,itable,inpath);
	}

/* close output file and summarize transfer */
    if (nlog > 0) {
	if (nreg > 1)
	    printf ("%d regions: %d / %d found\n",nreg,nstar,itot);
	else
	    printf ("1 region: %d / %d found\n",nstar,itable);
	}
    free (table);
    free (header);
    return (nstar);
}


char cdn[64]="/data/gsc1";	/* pathname of northern hemisphere GSC CDROM */
char cds[64]="/data/gsc2";	/* pathname of southern hemisphere GSC CDROM */

/* First region in each declination zone */
int zreg1[24]={1,594,1178,1729,2259,2781,3246,3652,4014,4294, 4492,4615,
	      4663,5260,5838,6412,6989,7523,8022,8464,8840,9134,9346,9490};

/* Last region in each declination zone */
int zreg2[24]={593,1177,1728,2258,2780,3245,3651,4013,4293,4491,4614,4662,
	       5259,5837,6411,6988,7522,8021,8463,8839,9133,9345,9489,9537};

/* Directory for each declination zone */
char zdir[24][8]={"n0000","n0730","n1500","n2230","n3000","n3730","n4500",
		"n5230","n6000","n6730","n7500","n8230","s0000","s0730",
		"s1500","s2230","s3000","s3730","s4500","s5230","s6000",
		"s6730","s7500","s8230"};

static struct Keyword rkw[15];
static int nrkw = 13;

/* GSCREG -- search the HST Guide Star Catalog index table for fields
 * in the specified range of coordinates and magnitudes.  Build a
 * list containing the pathnames of the files on the cd-rom.
 */

int gscreg (ra1,ra2,dec1,dec2, lhead,header, ltab,table, nrmax,rgns)

double ra1, ra2;	/* Right ascension limits in degrees */
double dec1, dec2; 	/* Declination limits in degrees */
int lhead;		/* Maximum length of FITS header in bytes */
char *header;		/* Table data fits header */
int ltab;		/* Maximum length of table buffer in bytes */
char *table;		/* Table data buffer */
int nrmax;		/* Maximum number of regions to find */
int *rgns;		/* Region numbers (returned)*/

{
    int nrgn;		/* Number of regions found (returned) */
    char tabpath[64];	/* Pathname for regions table */
    int nrows;		/* Number of entries in table */
    int nchar;		/* Number of characters per line in table */
    int nwrap;		/* 1 if 0h included in RA span*/
    int iwrap;
    double gscra(), gscdec();
    int gsczone();

    int verbose;
    char fitsline[120];
    char temp[8];
    int irow,iz1,iz2,ir1,ir2,jr1,jr2,i;
    int nsrch,nsrch1,nbhead,ift,nbr;
    double ralow, rahi;
    double declow, dechi, decmin, decmax;
    int regnum;
    int fitsrtline();

/* Set up keyword list for table entries to extract */
    strcpy (rkw[0].kname,"REG_NO");
    strcpy (rkw[1].kname,"RA_H_LOW");
    strcpy (rkw[2].kname,"RA_M_LOW");
    strcpy (rkw[3].kname,"RA_S_LOW");
    strcpy (rkw[4].kname,"RA_H_HI");
    strcpy (rkw[5].kname,"RA_M_HI");
    strcpy (rkw[6].kname,"RA_S_HI");
    strcpy (rkw[7].kname,"DECSI_LOW");
    strcpy (rkw[8].kname,"DEC_D_LOW");
    strcpy (rkw[9].kname,"DEC_M_LOW");
    strcpy (rkw[10].kname,"DECSI_HI");
    strcpy (rkw[11].kname,"DEC_D_HI");
    strcpy (rkw[12].kname,"DEC_M_HI");
    rkw[13].kname[0] = 0;
    for (i = 0; i < nrmax; i++)
	rgns[i] = 0;

    nrgn = 0;
    verbose = 0;

/* Set pathnames to guide star catalog cdroms */
    strcpy (tabpath,cdn);

/* Set pathname for region table file */
    strcat (tabpath,"/tables/regions.tbl");

/* Open the index table */
    ift = fitsrtopen (tabpath,lhead,header,nrkw,rkw, &nrows, &nchar, &nbhead);
    if (ift <= 0) {

/* If the northern hemisphere CDROM cannot be read, try the southern */
	strcpy (tabpath,cds);
	strcat (tabpath,"/tables/regions.tbl");
	ift = fitsrtopen (tabpath,lhead,header,nrkw,rkw,&nchar,&nrows,&nbhead);
	if (ift <= 0) {
	    printf ("GSCREG:  error reading region table %s\n",tabpath);
	    return (0);
	    }
	}

/* Find region range to search based on declination */
    iz1 = gsczone (dec1);
    iz2 = gsczone (dec2);
    jr1 = 0;
    jr2 = 0;
    nwrap = 1;

/* Search region northern hemisphere or only one region */
    if (iz2 >= iz1) {
	ir1 = zreg1[iz1];
	ir2 = zreg2[iz2];
	}

/* Search region in southern hemisphere with multiple regions */
    else if (dec1 < 0 && dec2 < 0) {
	ir1 = zreg1[iz2];
	ir2 = zreg2[iz1];
	}

/* Search region spans equator */
    else if (dec1 < 0 && dec2 >= 0) {
	ir1 = zreg1[12];
	ir2 = zreg2[iz1];
	jr1 = zreg1[0];
	jr2 = zreg2[iz2];
	nwrap = 2;
	}

    nsrch = ir2 - ir1 + 1;
    if (verbose)
	printf ("GSCREG: searching %d regions: %d - %d\n",nsrch,ir1,ir2);
    if (jr1 > 0) {
	nsrch1 = jr2 - jr1 + 1;
	if (verbose)
	    printf ("GSCREG: searching %d regions: %d - %d\n",nsrch1,jr1,jr2);
	}
    if (verbose)
	printf("GSCREG: RA: %.5f - %.5f, Dec: %.5f - %.5f\n",ra1,ra2,dec1,dec2);

    nrgn = 0;

    for (iwrap = 0; iwrap < nwrap; iwrap++) {

	for (irow = ir1 - 1; irow < ir2; irow++) {

	/* Read next line of region table */
	    nbr = fitsrtline (ift,nbhead,ltab,table,irow,nchar,fitsline);
	    if (nbr < nchar) {
		printf ("GSREG: %d / %d bytes read for row %d\n",nbr,nchar,irow);
		break;
		}

	/* Declination range of the gs region */
	/* note:  southern dechi and declow are reversed */
	    dechi = gscdec (fitsline, 10, 11, 12);
	    declow = gscdec (fitsline, 7, 8, 9);
	    if (dechi > declow) {
		decmin = declow;
		decmax = dechi;
		}
	    else {
		decmax = declow;
		decmin = dechi;
		}
	    if (decmax >= dec1 && decmin <= dec2) {

	    /* right ascension range of the gs region */
		ralow = gscra (fitsline, 1, 2, 3);
		rahi = gscra (fitsline, 4, 5, 6);
		if (rahi <= 0.0) rahi = 360.0;
		if (ralow > rahi) rahi = rahi + 360.0;

	/* Check RA if 0h RA not between region RA limits */
		if (ra1 < ra2) {
		    if (ralow <= ra1 && rahi >= ra2)  {

			/* Get region number from FITS table */
    			strncpy (temp, fitsline+rkw[0].kf, rkw[0].kl);
			regnum = atoi (temp);
			if (verbose)
			    printf ("GSCREG: Region %d added to search\n",regnum);

			/* Add this region to list, if there is space */
			if (nrgn < nrmax) {
			    rgns[nrgn] = regnum;
			    nrgn = nrgn + 1;
			    }
			}

	/* Check RA if 0h RA is between region RA limits */
		    else if (ralow >= ra2 && rahi <= ra1) {

			/* Get region number from FITS table */
			strncpy (temp, fitsline+rkw[0].kf, rkw[0].kl);
			regnum = atoi (temp);
			if (verbose)
			    printf ("GSCREG: Region %d added to search\n",regnum);

			/* Add this region to list, if there is space */
			if (nrgn < nrmax) {
			    rgns[nrgn] = regnum;
			    nrgn = nrgn + 1;
			    }
			}
		    }
		}
	    }

/* Handle wrap-around through the equator */
	ir1 = jr1;
	ir2 = jr2;
	jr1 = 0;
	jr2 = 0;
	}

    close (ift);
    return (nrgn);
}


/* GSCRA -- returns right ascension in degrees from the GSC index table
 *  This is converted from the hours, minutes, and seconds columns.
 */

double gscra (fitsline, hcol, mcol, scol)

char *fitsline;		/* index table line */
int hcol;		/* column index for hours */
int mcol;		/* column index for minutes */
int scol;		/* column index for seconds */

{
    double ra;		/* right ascension in fractional degrees */
    double hrs;		/* hours of right ascension */
    double min;		/* minutes of right ascension */
    double sec;		/* seconds of right ascension */
    char temp[16];

/*  hours of right ascension */
    strncpy (temp, fitsline+rkw[hcol].kf, rkw[hcol].kl);
    hrs = atof (temp);

/* minutes of right ascension */
    strncpy (temp, fitsline+rkw[mcol].kf, rkw[mcol].kl);
    min = atof (temp);

/* seconds of right ascension */
    strncpy (temp, fitsline+rkw[scol].kf, rkw[scol].kl);
    sec = atof (temp);

/* right ascension in degrees */
    ra = hrs + (min / 60.0) + (sec / 3600.0);
    ra = ra * 15.0;

    return (ra);
}


/*  GSCDEC -- returns the declination in degrees from the GSC index table.
 *  This is converted from the sign, degrees, minutes, and seconds columns.
 */

double gscdec (fitsline, sgncol, dcol, mcol)

char *fitsline;		/* index table line */
int sgncol;		/* column index for sign */
int dcol;		/* column index for degrees */
int mcol;		/* column index for minutes */

{
double	dec;		/* declination in fractional degrees */

char sign[4];		/* sign of declination */
double deg;		/* degrees of declination*/
double min;		/* minutes of declination */
char temp[8];

/* get declination sign from table */
    strncpy (sign, fitsline+rkw[sgncol].kf, rkw[sgncol].kl);

/* get degrees of declination from table */
    strncpy (temp, fitsline+rkw[dcol].kf, rkw[dcol].kl);
    deg = atof (temp);

/* get minutes of declination from table */
    strncpy (temp, fitsline+rkw[mcol].kf, rkw[mcol].kl);
    min = atof (temp);

    dec = deg + (min / 60.0);

/* negative declination */
    if (strchr (sign, '-') != NULL)
	dec = -dec;

    return (dec);
}


/*  GSCZONE -- find the zone number where a declination can be found */

int gsczone (dec)

double dec;	/* declination in degrees */

{
int zone;		/* gsc zone (returned) */
double	zonesize;
int numdir = 24;	/* number of declination zone directories */
int ndeczones = 12;	/* number of declination zones per hemisphere */

/* width of declination zones */
    zonesize = 90.0 / ndeczones;

    zone = ((int) (dec / zonesize));
    if (dec < 0)
	zone = ndeczones - zone;
	
    return (zone);
}


/* GSCPATH -- Get HST Guide Star Catalog region FITS file pathname */

gscpath (regnum,path)

int regnum;	/* Guide Star Catalog region number */
char *path;	/* Pathname of GSC region FITS file */

{
    int zone;		/* Name of Guide Star Catalog zone directory */
    char root[20];	/* Name of Guide Star Catalog region file */
    int i;

/* get zone directory name given region number */
    for (i = 0; i < 24; i++) {
	if (regnum >= zreg1[i] && regnum <= zreg2[i]) {
	    zone = i;
	    break;
	    }
	}

/* Set the pathname using the appropriate GSC CDROM directory */

/* Northern hemisphere disk (volume 1) */
    if (regnum < zreg1[13])
	sprintf (path,"%s/gsc/%s/%04d.gsc", cdn, zdir[zone], regnum);

/* Southern hemisphere disk (volume 2) */
    else
	sprintf (path,"%s/%s/gsc/%04d.gsc", cds, zdir[zone], regnum);

    return;
}
