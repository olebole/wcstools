/*** File libwcs/gscread.c
 *** November 15, 1996
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "fitshead.h"
#include "wcs.h"

char cdn[64]="/data/gsc1";	/* pathname of northern hemisphere GSC CDROM */
char cds[64]="/data/gsc2";	/* pathname of southern hemisphere GSC CDROM */

static void gscpath();
static int gscreg();

/* GSCREAD -- Read HST Guide Star Catalog stars from CDROM */

int
gscread (cra,cdec,dra,ddec,drad,mag1,mag2,classd,nstarmax,
	 gnum,gra,gdec,gmag,gtype,nlog)

double	cra;		/* Search center J2000 right ascension in degrees */
double	cdec;		/* Search center J2000 declination in degrees */
double	dra;		/* Search half width in right ascension in degrees */
double	ddec;		/* Search half-width in declination in degrees */
double	drad;		/* Limiting separation in degrees (ignore if 0) */
double	mag1,mag2;	/* Limiting magnitudes (none if equal) */
int	classd;		/* Desired object class (-1=all, 0=stars, 3=nonstars) */
int	nstarmax;	/* Maximum number of stars to be returned */
double	*gnum;		/* Array of Guide Star numbers (returned) */
double	*gra;		/* Array of right ascensions (returned) */
double	*gdec;		/* Array of declinations (returned) */
double	*gmag;		/* Array of magnitudes (returned) */
int	*gtype;		/* Array of object types (returned) */
int	nlog;		/* 1 for diagnostics */
{
    double ra1,ra2;	/* Limiting right ascensions of region in degrees */
    double dec1,dec2;	/* Limiting declinations of region in degrees */
    double dist = 0.0;  /* Distance from search center in degrees */
    char *table;	/* FITS table */
    int nreg;		/* Number of input FITS tables files */
    int rlist[100];	/* List of input FITS tables files */
    char inpath[64];	/* Pathname for input FITS table file */
    char entry[100];	/* Buffer for FITS table row */
    int class, class0;	/* Object class (0>star, 3>other) */
    struct Keyword kw[8];	/* Keyword structure */
    struct Keyword *kwn;

    int verbose;
    int wrap;
    int rnum, num0, num, itot,ireg,ltab;
    int ik,nk,itable,ntable,jstar;
    int nbline,npos,nbhead;
    int nbr,nrmax,nstar,i;
    FILE *ift;
    double ra,ra0,rasum,dec,dec0,decsum,perr,perr0,perr2,perrsum,msum;
    double mag,mag0,merr,merr0,merr2,merrsum;
    char *str;

    itot = 0;
    if (nlog == 1)
	verbose = 1;
    else
	verbose = 0;
    ltab = 10000;
    table = malloc (10000);
    for (i = 0; i < 100; i++)
	entry[i] = 0;

    /* Set path to Guide Star Catalog */
    if ((str = getenv("GSC_NORTH")) != NULL )
	strcpy (cdn,str);
    if ((str = getenv("GSC_SOUTH")) != NULL )
	strcpy (cds,str);

    /* Set right ascension limits for search */
    ra1 = cra - dra;
    ra2 = cra + dra;

    /* Keep right ascension between 0 and 360 degrees */
    if (ra1 < 0.0)
	ra1 = 360.0 + ra1;
    if (ra2 > 360.0)
	ra2 = ra2 - 360.0;

    /* Set declination limits for search */
    dec1 = cdec - ddec;
    dec2 = cdec + ddec;

    /* dec1 is always the smallest declination */
    if (dec1 > dec2) {
	dec = dec1;
	dec1 = dec2;
	dec2 = dec;
	}

    /* Search zones which include the poles cover 360 degrees in RA */
    if (dec1 < -90.0) {
	dec1 = -90.0;
	ra1 = 0.0;
	ra2 = 359.99999;
	}
    if (dec2 > 90.0) {
	dec2 = 90.0;
	ra1 = 0.0;
	ra2 = 359.99999;
	}
    if (verbose) {
	char rstr1[16],rstr2[16],dstr1[16],dstr2[16];
	ra2str (rstr1, ra1, 3);
        dec2str (dstr1, dec1, 2);
	ra2str (rstr2, ra2, 3);
        dec2str (dstr2, dec2, 2);
	fprintf (stderr,"GSCREAD: RA: %s - %s  Dec: %s - %s\n",
		 rstr1,rstr2,dstr1,dstr2);
	}

/* If RA range includes zero, split it in two */
    wrap = 0;
    if (ra1 > ra2)
	wrap = 1;
    else
	wrap = 0;

/* make mag1 always the smallest magnitude */
    if (mag2 < mag1) {
	mag = mag2;
	mag2 = mag1;
	mag1 = mag;
	}

/* Find Guide Star Catalog regions in which to search */
    nrmax = 100;
    nreg = gscreg (ra1,ra2,dec1,dec2,ltab,table,nrmax,rlist,verbose);
    if (nreg <= 0) {
	fprintf (stderr,"GSCREAD:  no Guide Star regions found\n");
	free (table);
	return (0);
	}

    classd = -1;

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
    nstar = 0;

/* Loop through region list */
    for (ireg = 0; ireg < nreg; ireg++) {
	gscpath (rlist[ireg],inpath);

    /* Read size and keyword info from FITS table header */
	kwn = kw;
	ift = fitsrtopen (inpath,&nk,&kwn,&ntable,&nbline,&nbhead);

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
	class = 0;

	/* Loop through FITS table for this region */
	for (itable = 0; itable <= ntable; itable++) {

	    if (itable < ntable) {
		nbr = fitsrtline (ift,nbhead,ltab,table,itable,nbline,entry);
		if (nbr < nbline) {
		    fprintf (stderr,"GSCREAD: %d / %d bytes read, line %d / %d, region %d\n",
			      nbr,nbline,itable,ntable,rnum);
		    break;
		    }

	 /* Extract selected fields */

		/* Star number within region */
		num0 = ftgeti4 (entry, &kw[0]);

		/* Right ascension in degrees */
		ra0 = ftgetr8 (entry, &kw[1]);

		/* Declination in degrees */
		dec0 = ftgetr8 (entry, &kw[2]);

		/* Position error */
		perr0 = ftgetr8 (entry, &kw[3]);

		/* Magnitude */
		mag0 = ftgetr8 (entry, &kw[4]);

		/* Magnitude error */
		merr0 = ftgetr8 (entry, &kw[5]);

		/* Bandpass code */
		/* band0 = ftgeti4 (entry, &kw[6]); */

		/* Object class code */
		class0 = ftgeti4 (entry, &kw[7]);
		}
	    else
		num0 = 0;

	/* Compute mean position and magnitude for object */
	    if (num != num0 && itable > 0 && npos > 0) {
		ra = rasum / perrsum;
		dec = decsum / perrsum;
		mag = msum / merrsum;
		if (drad > 0)
		    dist = wcsdist (cra,cdec,ra,dec);

	/* Check magnitude amd position limits */
		if (((mag1 != mag2 && (mag >= mag1 && mag <= mag2)) ||
		    (mag1 == mag2)) &&
		    ((wrap && (ra <= ra1 || ra >= ra2)) ||
		    (!wrap && (ra >= ra1 && ra <= ra2))) &&
		    (drad == 0.0 || dist < drad) &&
     		    (dec >= dec1 && dec <= dec2)) {

	/* Save star position in table */
		    if (nstar <= nstarmax) {
			gnum[nstar] = (double)rnum + (0.0001 * (double) num);
			gra[nstar] = ra;
			gdec[nstar] = dec;
			gmag[nstar] = mag;
			gtype[nstar] = class;
			}
		    nstar = nstar + 1;
		    jstar = jstar + 1;
		    if (nlog == 1)
			fprintf (stderr,"GSCREAD: %04d.%04d: %9.5f %9.5f %5.2f %d %d\n",
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
		fprintf (stderr,"GSCREAD: %4d / %4d: %5d / %5d  / %5d sources, region %4d.%04d\n",
			 ireg,nreg,jstar,itable,ntable,rlist[ireg],num0);

/* End of region */
	    }

/* Close region input file */
	(void) fclose (ift);
	itot = itot + itable;
	if (nlog > 0)
	    fprintf (stderr,"GSCREAD: %4d / %4d: %5d / %5d  / %5d sources from region %4d\n",
		     ireg+1,nreg,jstar,itable,ntable,rlist[ireg]);
	}

/* close output file and summarize transfer */
    if (nlog > 0) {
	if (nreg > 1)
	    fprintf (stderr,"GSCREAD: %d regions: %d / %d found\n",nreg,nstar,itot);
	else
	    fprintf (stderr,"GSCREAD: 1 region: %d / %d found\n",nstar,itable);
	}
    free (table);
    if (nstar > nstarmax) {
	fprintf (stderr,"GSCREAD: %d stars found; only %d returned\n",
		 nstar,nstarmax);
	nstar = nstarmax;
	}
    return (nstar);
}


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

static int
gscreg (ra1, ra2, dec1, dec2, ltab, table, nrmax, rgns, verbose)

double	ra1, ra2;	/* Right ascension limits in degrees */
double	dec1, dec2; 	/* Declination limits in degrees */
int	ltab;		/* Maximum length of table buffer in bytes */
char	*table;		/* Table data buffer */
int	nrmax;		/* Maximum number of regions to find */
int	*rgns;		/* Region numbers (returned)*/
int	verbose;	/* 1 for diagnostics */

{
    int nrgn;		/* Number of regions found (returned) */
    char tabpath[64];	/* Pathname for regions table */
    int nrows;		/* Number of entries in table */
    int nchar;		/* Number of characters per line in table */
    int nwrap;		/* 1 if 0h included in RA span*/
    int iwrap;
    struct Keyword *kwn;
    double gscra(), gscdec();
    int gsczone();

    char fitsline[120];
    int irow,iz1,iz2,ir1,ir2,jr1,jr2,i;
    int nsrch,nsrch1,nbhead,nbr;
    FILE *ift;
    double ralow, rahi;
    double declow, dechi, decmin, decmax;
    int regnum;

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

/* Set pathnames to guide star catalog cdroms */
    strcpy (tabpath,cdn);

/* Set pathname for region table file */
    strcat (tabpath,"/tables/regions.tbl");

/* Open the index table */
    kwn = rkw;
    ift = fitsrtopen (tabpath,&nrkw,&kwn, &nrows, &nchar, &nbhead);
    if (!ift) {

/* If the northern hemisphere CDROM cannot be read, try the southern */
	strcpy (tabpath,cds);
	strcat (tabpath,"/tables/regions.tbl");
	ift = fitsrtopen (tabpath,&nrkw,&kwn,&nchar,&nrows,&nbhead);
	if (!ift) {
	    fprintf (stderr,"GSCREG:  error reading region table %s\n",tabpath);
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
	fprintf (stderr,"GSCREG: searching %d regions: %d - %d\n",nsrch,ir1,ir2);
    if (jr1 > 0) {
	nsrch1 = jr2 - jr1 + 1;
	if (verbose)
	    fprintf (stderr,"GSCREG: searching %d regions: %d - %d\n",nsrch1,jr1,jr2);
	}
    if (verbose)
	fprintf(stderr,"GSCREG: RA: %.5f - %.5f, Dec: %.5f - %.5f\n",ra1,ra2,dec1,dec2);

    nrgn = 0;

    for (iwrap = 0; iwrap < nwrap; iwrap++) {

	for (irow = ir1 - 1; irow < ir2; irow++) {

	/* Read next line of region table */
	    nbr = fitsrtline (ift,nbhead,ltab,table,irow,nchar,fitsline);
	    if (nbr < nchar) {
		fprintf (stderr,"GSREG: %d / %d bytes read for row %d\n",nbr,nchar,irow);
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

	    /* Right ascension range of the Guide Star Catalog region */
		ralow = gscra (fitsline, 1, 2, 3);
		rahi = gscra (fitsline, 4, 5, 6);
		if (rahi <= 0.0) rahi = 360.0;

	    /* Check RA if 0h RA not between region RA limits */
		if (ra1 < ra2) {
		    if (ralow <= ra2 && rahi >= ra1) {

			/* Get region number from FITS table */
			regnum = ftgeti4 (fitsline, &rkw[0]);
			if (verbose)
			    fprintf (stderr,"GSCREG: Region %d added to search\n",regnum);

			/* Add this region to list, if there is space */
			if (nrgn < nrmax) {
			    rgns[nrgn] = regnum;
			    nrgn = nrgn + 1;
			    }
			}
		    }

	    /* Check RA if 0h RA is between region RA limits */
		else {
		    if (ralow > rahi) rahi = rahi + 360.0;
		    if (ralow <= ra2 || rahi >= ra1) {

		    /* Get region number from FITS table */
			regnum = ftgeti4 (fitsline, &rkw[0]);
			if (verbose)
			    fprintf (stderr,"GSCREG: Region %d added to search\n",regnum);

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

    (void) fclose (ift);
    return (nrgn);
}


/* GSCRA -- returns right ascension in degrees from the GSC index table
 *  This is converted from the hours, minutes, and seconds columns.
 */

double
gscra (fitsline, hcol, mcol, scol)

char *fitsline;		/* Index table line */
int hcol;		/* Column index for hours */
int mcol;		/* Column index for minutes */
int scol;		/* Column index for seconds */

{
    double ra;		/* Right ascension in fractional degrees */
    double hrs;		/* Hours of right ascension */
    double min;		/* Minutes of right ascension */
    double sec;		/* Seconds of right ascension */

/*  hours of right ascension */
    hrs = ftgetr8 (fitsline, &rkw[hcol]);

/* minutes of right ascension */
    min = ftgetr8 (fitsline, &rkw[mcol]);

/* seconds of right ascension */
    sec = ftgetr8 (fitsline, &rkw[scol]);

/* right ascension in degrees */
    ra = hrs + (min / 60.0) + (sec / 3600.0);
    ra = ra * 15.0;

    return (ra);
}


/*  GSCDEC -- returns the declination in degrees from the GSC index table.
 *  This is converted from the sign, degrees, minutes, and seconds columns.
 */

double
gscdec (fitsline, sgncol, dcol, mcol)

char *fitsline;		/* Index table line */
int sgncol;		/* Column index for sign */
int dcol;		/* Column index for degrees */
int mcol;		/* Column index for minutes */

{
double	dec;		/* Declination in fractional degrees */

char sgn[4];		/* Sign of declination */
double deg;		/* Degrees of declination*/
double min;		/* Minutes of declination */

/* get declination sign from table */
    (void) ftgetc (fitsline, &rkw[sgncol], sgn, 3);

/* get degrees of declination from table */
    deg = ftgetr8 (fitsline, &rkw[dcol]);

/* get minutes of declination from table */
    min = ftgetr8 (fitsline, &rkw[mcol]);

    dec = deg + (min / 60.0);

/* negative declination */
    if (strchr (sgn, '-') != NULL)
	dec = -dec;

    return (dec);
}


/*  GSCZONE -- find the zone number where a declination can be found */

int
gsczone (dec)

double dec;	/* declination in degrees */

{
int zone;		/* gsc zone (returned) */
double	zonesize;
int ndeczones = 12;	/* number of declination zones per hemisphere */

/* width of declination zones */
    zonesize = 90.0 / ndeczones;

    zone = ((int) (dec / zonesize));
    if (dec < 0)
	zone = ndeczones - zone;
	
    return (zone);
}


/* GSCPATH -- Get HST Guide Star Catalog region FITS file pathname */

static void
gscpath (regnum,path)

int regnum;	/* Guide Star Catalog region number */
char *path;	/* Pathname of GSC region FITS file */

{
    int zone;		/* Name of Guide Star Catalog zone directory */
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
	sprintf (path,"%s/gsc/%s/%04d.gsc", cds, zdir[zone], regnum);

    return;
}

/* Feb 16 1996	New program
 * May 14 1996	Change arguments to FITSRTOPEN
 * May 17 1996	Fix bug so equal magnitude limits accepts anything
 * May 17 1996	Use new FITSRTOPEN which internalizes FITS header
 * May 24 1996	Fix string decoding bug and region search
 * May 31 1996	Use stream I/O
 * Jun 10 1996	Remove unused variables after using lint
 * Jul  1 1996	Fix GSC pathname
 * Aug  6 1996	Minor changes after lint
 * Sep 17 1996	Fix bug causing incomplete region access; improve logging
 * Oct 17 1996	Change FITS table entry calls
 * Nov 13 1996	Return no more than maximum star number
 * Nov 13 1996	Write all error messages to stderr with subroutine names
 * Nov 15 1996	Implement search radius; change input arguments
 */
