/* File libwcs/wcscat.h
 * September 13, 1999
 * By Doug Mink, SAO
 */

/* Source catalog flags and programs */

#define GSC		1	/* refcat value for HST Guide Star Catalog */
#define UJC		2	/* refcat value for USNO UJ Star Catalog */
#define UAC		3	/* refcat value for USNO A Star Catalog */
#define USAC		4	/* refcat value for USNO SA Star Catalog */
#define SAO		5	/* refcat value for SAO Star Catalog */
#define IRAS		6	/* refcat value for IRAS Point Source Catalog */
#define PPM		7	/* refcat value for PPM Star Catalog */
#define TYCHO		8	/* refcat value for Tycho Star Catalog */
#define UA1		9	/* refcat value for USNO A-1.0 Star Catalog */
#define UA2		10	/* refcat value for USNO A-2.0 Star Catalog */
#define USA1		11	/* refcat value for USNO SA-1.0 Star Catalog */
#define USA2		12	/* refcat value for USNO SA-2.0 Star Catalog */
#define HIP		13	/* refcat value for Hipparcos Star Catalog */
#define ACT		14	/* refcat value for USNO ACT Star Catalog */
#define BSC		15	/* refcat value for Yale Bright Star Catalog */
#define TABCAT		-1	/* refcat value for StarBase tab table catalog */
#define BINCAT		-2	/* refcat value for TDC binary catalog */
#define TXTCAT		-3	/* refcat value for TDC ASCII catalog */

/* Subroutines for dealing with catalogs */
int RefCat();
void CatNum();
int CatNumLen();
void SearchLim();
void RefLim();

/* Subroutines for extracting sources from catalogs by sky region */
int gscread();		/* Read sources from HST Guide Star Catalog */
int uacread();		/* Read sources from USNO A or SA Catalog */
int ujcread();		/* Read sources from USNO J Catalog */
int tabread();		/* Read sources from tab table catalog */
int binread();		/* Read sources from SAO TDC binary format catalog */
int catread();		/* Read sources from SAO TDC ASCII format catalog */
int actread();		/* Read sources from USNO ACT Catalog */

/* Subroutines for extracting sources from catalogs by ID number */
int gscrnum();		/* Read sources from HST Guide Star Catalog */
int uacrnum();		/* Read sources from USNO A or SA Catalog */
int ujcrnum();		/* Read sources from USNO J Catalog */
int tabrnum();		/* Read sources from tab table catalog */
int binrnum();		/* Read sources from SAO TDC binary format catalog */
int catrnum();		/* Read sources from SAO TDC ASCII format catalog */
int actrnum();		/* Read sources from USNO ACT Catalog */
void setgsclass();	/* Set GSC object class */
void setuplate();	/* Set USNO catalog plate number to search */
int getuplate();	/* Get USNO catalog plate number to search */
void settabkey();	/* Set tab table keyword to read for object */
int catstar();		/* Read one star entry from ASCII catalog, 0 if OK */
int binstar();		/* Read one star entry from binary catalog, 0 if OK */
int tabstar();		/* Read one star entry from tab table catalog, 0 if OK */

/* Subroutines for sorting tables of star positions and magnitudes */
void XSortStars();
void RASortStars();
void MagSortStars();

/* Data structures for SAO TDC ASCII and binary star catalogs */

struct Star {
    float rdum;
    float xno;		/* Catalog number */
    double ra;		/* Right Ascension (degrees) */
    double dec;		/* Declination (degrees) */
    char isp[2];	/* Spectral type or other 2-char identifier */
    short mag[11];	/* Up to 10 Magnitudes * 100 */
    double rapm;	/* RA proper motion (degrees per year) */
    double decpm;	/* Dec proper motion (degrees per year) */
    double xmag[11];	/* Up to 10 Magnitudes */
    double num;		/* Actual star number */
    int coorsys;	/* Coordinate system (WCS_J2000, WCS_B1950,...) */
    double equinox;	/* Equinox of coordinate system as fractional year */
    double epoch;	/* Epoch of position as fractional year */
    char objname[32];	/* Object name */
    int peak;		/* Peak flux per pixel in star image */
};

struct StarCat {
    int star0;		/* Subtract from star number for file sequence number */
    int star1;		/* First star number in file */
    int nstars;		/* Number of stars in file */
    int stnum;		/* Star number format in catalog file:
			  <0: -stnum-character name at end instead of number
			   0:  no star i.d. numbers
			   1: Real*4 star i.d. numbers
			   2: Integer*4 <region><nnnn>
			   3: Integer*4 <region><nnnnn>
			   4: Integer*4 <nnnnnnnnn>
			   5: Character ID instead of number in ASCII files */
    int mprop;		/* 1 if proper motion is included */
    int nmag;		/* Number of magnitudes present
			   Negative for J2000 catalog */
    int nbent;		/* Number of bytes per star entry */
    int	rasorted;	/* 1 if RA-sorted, else 0 */
    FILE *ifcat;	/* File descriptor for catalog file */
    char isfil[24];	/* Star catalog file name */
    char isname[64];	/* Star catalog description */
    int  byteswapped;	/* 1 if catalog is byte-reversed from CPU */
    int  coorsys;	/* Coordinate system
			   B1950 J2000 Galactic Ecliptic */
    double epoch;	/* Epoch of catalog coordinates in years */
    double equinox;	/* Equinox of catalog coordinates in years */
    char inform;	/* Coordinate format
			   (B>inary D>egrees H>MS T>able U>SNO) */
    char incdir[128];	/* Catalog directory pathname */
    char incfile[32];	/* Catalog file name */
    int ncobj;		/* Length of object name in binary star entry */
    int nndec;		/* Number of decimal places in star number */
    int nepoch;		/* 1 if epoch of coordinates is present */
    char *catbuff;	/* Pointer to start of catalog */
    char *catdata;	/* Pointer to first entry in catalog */
    char *catline;	/* Pointer to current entry in catalog */
    char *catlast;	/* Pointer to one past end of last entry in catalog */
    int  istar;		/* Number of current catalog entry */
    struct TabTable *startab;	/* Structure for tab table catalog */
    int entid;		/* Entry number for ID */
    int entra;		/* Entry number for right ascension */
    int entdec;		/* Entry number for declination */
    int entmag;		/* Entry number for magnitude */
    int entpeak;	/* Entry number for peak counts */
    int entepoch;	/* Entry number for epoch of observation */
    int entname;	/* Entry number for object name */
    int entprop;	/* Entry number for proper motion */
    int entkey;		/* Entry number for additional keyword */
};

/* Subroutines for reading headers of TDC binary and ASCII catalogs */
int isbin();
struct StarCat *binopen();
void binclose();
struct StarCat *catopen();
void catclose();

/* Data structure for tab table files */
struct TabTable {
    char *filename;	/* Name of tab table file */
    int nlines;		/* Number of entries in table */
    char *tabbuff;	/* Pointer to start of saved tab table in memory */
    char *tabhead;	/* Pointer to start of line containing table header */
    char *tabdata;	/* Pointer to start of first line of table data */
    int iline;		/* Number of current line (1=first) */
    char *tabline;	/* Pointer to start of current line */
    int ncols;		/* Number of columns per table entry */
    char **colname;	/* Column names */
    int *lcol;		/* Lengths of column header names */
    int *lcfld;		/* Number of columns in field (hyphens) */
};

/* Subroutines for extracting tab table information */
struct TabTable *tabopen();	/* Open tab table file */
struct StarCat *tabcatopen();	/* Open tab table catalog */
void tabcatclose();	/* Close tab table catalog */
char *tabline();	/* Find a specified line in a tab table */
int tabrkey();		/* Keyword values from tab table catalogs */
int tabcol();		/* Find column for name */
int tabgetk();		/* Get tab table entries for named column */
int tabgetc();		/* Get tab table entry for named column */
void tabclose();	/* Free all arrays left open by tab table structure */
int istab();

#define MAXRANGE 20

/* Structure for dealing with ranges */
struct Range {
    double first;	/* Current minimum value */
    double last;	/* Current maximum value */
    double step;	/* Current step in value */
    double value;	/* Current value */
    double ranges[MAXRANGE*3];	/* nranges sets of first, last, step */
    int nvalues;	/* Total number of values in all ranges */
    int nranges;	/* Number of ranges */
    int irange;		/* Index of current range */
};

/* Subroutines for dealing with ranges */
struct Range *RangeInit();	/* Initialize range structure from string */
int isrange();		/* Return 1 if string is a range of numbers, else 0 */
int rgetn();		/* Return number of values in all ranges */
int rgeti4();		/* Return next number in range as integer */
double rgetr8();	/* Return next number in range as double */
void rstart();		/* Restart range */

/* Subroutines for translating dates and times */
double fd2ep();	/* FITS standard date string to fractional year (epoch) */
double fd2jd();	/* FITS standard date string to Julian date */
void jd2dt();	/* Julian date to yyyy.mmdd hh.mmssss */
double dt2jd();	/* yyyy.ddmm and hh.mmsss to Julian date */
void jd2dt();	/* Julian date to yyyy.mmdd hh.mmssss */
double dt2ts();	/* yyyy.ddmm and hh.mmsss to seconds since 1950.0 */ 
void ts2dt();	/* seconds since 1950.0 to yyyy.mmdd hh.mmssss */
double dt2ep();	/* yyyy.ddmm and hh.mmsss to fractional year (epoch) */
void ep2dt();	/* fractional year to yyyy.mmdd hh.mmssss */
double ts2jd();	/* seconds since 1950.0 to yyyy.mmdd hh.mmssss */
double jd2ts();	/* Julian date to seconds since 1950.0 */
void ts2i();	/* seconds since 1950.0 to year, month, day, hours, min, sec */
double ep2jd();	/* fractional year to Julian Date */

/* Shapes for SAOimage region file output */
#define WCS_CIRCLE 1	/* shape for SAOimage plotting */
#define WCS_SQUARE 2	/* shape for SAOimage plotting */
#define WCS_DIAMOND 3	/* shape for SAOimage plotting */
#define WCS_CROSS 4	/* shape for SAOimage plotting */
#define WCS_EX 5	/* shape for SAOimage plotting */
#define WCS_VAR 6	/* shape for HST GSC SAOimage plotting (+ and x)*/

/* Sep 22 1998  New header file (star.h)
 * Oct 16 1998  Add more options for ASCII catalogs
 * Oct 20 1998  Add object name to binary files
 * Oct 21 1998	New file (wcscat.h)
 * Oct 26 1998	Combined wcscat.h and star.h
 * Oct 27 1998	Add SAOimage region shapes
 * Nov  9 1998	Add rasorted flag to catalog structure
 * Nov 20 1998	Add support for USNO A-2.0 and SA-2.0 catalogs
 * Dec  8 1998	Add support for the Hipparcos and ACT catalogs
 *
 * Jan 25 1999	Add declarations for tab table access
 * Jan 25 1999	Add declarations for dealing with ranges of numbers
 * Feb  2 1999	Add number of decimal places in star number to StarCat
 * Feb 11 1999	Add coordinate system info to star structure
 * Feb 11 1999	Change starcat.insys to starcat.coorsys for consistency
 * May 14 1999	Update Star and StarCat structure to cover tab tables
 * May 19 1999	Update StarCat structure to include epoch from catalog
 * June 4 1999	Add CatNumLen()
 * Jun 14 1999	Add SearchLim()
 * Jun 30 1999	Add isrange()
 * Jul  1 1999	Add declarations for date/time conversions in dateutil.c
 * Jul  2 1999	Add rstart()
 * Jul 26 1999	Add Yale Bright Star Catalog
 * Aug 16 1999	Add RefLim() to get converted search coordinates right
 * Aug 25 1999	Add ACT catalog
 * Sep 10 1999	Move special case setting from argument list to subroutines
 * Sep 13 1999	Add subroutines to access data structure for single stars
 */
