/* File libwcs/wcscat.h
 * March 27, 2000
 * By Doug Mink, dmink@cfa.harvard.edu
 */

/* Source catalog flags and subroutines */

/* Source catalog flags returned from RefCat */
#define GSC		1	/* HST Guide Star Catalog */
#define UJC		2	/* USNO UJ Star Catalog */
#define UAC		3	/* USNO A Star Catalog */
#define USAC		4	/* USNO SA Star Catalog */
#define SAO		5	/* SAO Star Catalog */
#define IRAS		6	/* IRAS Point Source Catalog */
#define PPM		7	/* PPM Star Catalog */
#define TYCHO		8	/* Tycho Star Catalog */
#define UA1		9	/* USNO A-1.0 Star Catalog */
#define UA2		10	/* USNO A-2.0 Star Catalog */
#define USA1		11	/* USNO SA-1.0 Star Catalog */
#define USA2		12	/* USNO SA-2.0 Star Catalog */
#define HIP		13	/* Hipparcos Star Catalog */
#define ACT		14	/* USNO ACT Star Catalog */
#define BSC		15	/* Yale Bright Star Catalog */
#define TYCHO2		16	/* Tycho-2 Star Catalog */
#define TABCAT		-1	/* StarBase tab table catalog */
#define BINCAT		-2	/* TDC binary catalog */
#define TXTCAT		-3	/* TDC ASCII catalog */

/* Subroutines for dealing with catalogs */
int RefCat();		/* Return catalog type code, title, coord. system */
char *ProgCat();	/* Return catalog name given program name used */
char *ProgName();	/* Return program name given program path used */
int PropCat();		/* Return 1 if catalog has proper motions, else 0 */
void CatNum();		/* Return formatted source number */
int CatNumLen();	/* Return length of source numbers */
void SearchLim();	/* Compute limiting RA and Dec */
void RefLim();		/* Compute limiting RA and Dec in new system */
int isfile();		/* Return 1 if string is name of readable file */
int agets();		/* Extract value from keyword= value in string */

/* Subroutines for extracting sources from catalogs by sky region */
int gscread();		/* Read sources from HST Guide Star Catalog */
int uacread();		/* Read sources from USNO A or SA Catalog */
int ujcread();		/* Read sources from USNO J Catalog */
int tabread();		/* Read sources from tab table catalog */
int binread();		/* Read sources from SAO TDC binary format catalog */
int ctgread();		/* Read sources from SAO TDC ASCII format catalog */
int actread();		/* Read sources from USNO ACT Catalog */

/* Subroutines for extracting sources from catalogs by ID number */
int gscrnum();		/* Read sources from HST Guide Star Catalog */
int uacrnum();		/* Read sources from USNO A or SA Catalog */
int ujcrnum();		/* Read sources from USNO J Catalog */
int tabrnum();		/* Read sources from tab table catalog */
int binrnum();		/* Read sources from SAO TDC binary format catalog */
int ctgrnum();		/* Read sources from SAO TDC ASCII format catalog */
int actrnum();		/* Read sources from USNO ACT Catalog */
void setgsclass();	/* Set GSC object class */
void setuplate();	/* Set USNO catalog plate number to search */
int getuplate();	/* Get USNO catalog plate number to search */
void settabkey();	/* Set tab table keyword to read for object */
int ctgstar();		/* Read one star entry from ASCII catalog, 0 if OK */
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

/* Catalog proper motion units */
#define PM_MASYR		1	/* milliarcseconds per year */
#define PM_ARCSECYR		2	/* arcseconds per year */
#define PM_DEGYR		3	/* degrees per year */
#define PM_RADYR		4	/* radians per year */
#define PM_TSECYR		5	/* seconds of time (RA) per year */

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
    int entmag1;	/* Entry number for first or only magnitude */
    int entmag2;	/* Entry number for second magnitude, if present */
    int entpeak;	/* Entry number for peak counts */
    int entepoch;	/* Entry number for epoch of observation */
    int entname;	/* Entry number for object name */
    int entadd;		/* Entry number for additional keyword */
    int entrpm;		/* Entry number for proper motion in right ascension */
    int entdpm;		/* Entry number for proper motion in declination */
    int rpmunit;	/* Units for RA proper motion (PM_x) */
    int dpmunit;	/* Units for DEC proper motion (PM_x) */
    char keyid[16];	/* Entry name for ID */
    char keyra[16];	/* Entry name for right ascension */
    char keydec[16];	/* Entry name for declination */
    char keymag1[16];	/* Entry name for first or only magnitude */
    char keymag2[16];	/* Entry name for second magnitude, if present */
    char keyrpm[16];	/* Entry name for right ascension proper motion */
    char keydpm[16];	/* Entry name for declination proper motion */
    char keypeak[16];	/* Entry name for integer code */
    char keyadd[16];	/* Entry name for additional keyword */
};

/* Subroutines for reading headers of TDC binary and ASCII catalogs */
int isbin();
struct StarCat *binopen();
void binclose();
struct StarCat *ctgopen();
void ctgclose();

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
int tabxyread();	/* Read x, y, and magnitude from tab table star list */
char *tabline();	/* Find a specified line in a tab table */
int tabrkey();		/* Keyword values from tab table catalogs */
int tabcol();		/* Find column for name */
int tabgetk();		/* Get tab table entries for named column */
int tabgetc();		/* Get tab table entry for named column */
int tabgeti4();		/* Return 4-byte integer from tab table line */
double tabgetra();	/* Return right ascension in degrees from tab table*/
double tabgetdec();	/* Return declination in degrees from tab table*/
double tabgetpm();	/* Return RA or Dec p.m. in degrees from tab table*/
double tabgetr8();	/* Return double number from tab table line */
void tabclose();	/* Free all arrays left open by tab table structure */
char *tgettaberr();	/* Retrun most recent tab table error message */
int istab();
int gettabndec();	/* Return number of decimal places in tab catalog ids */

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

/* Shapes for SAOimage region file output */
#define WCS_CIRCLE 1	/* shape for SAOimage plotting */
#define WCS_SQUARE 2	/* shape for SAOimage plotting */
#define WCS_DIAMOND 3	/* shape for SAOimage plotting */
#define WCS_CROSS 4	/* shape for SAOimage plotting */
#define WCS_EX 5	/* shape for SAOimage plotting */
#define WCS_VAR 6	/* shape for HST GSC SAOimage plotting (+ and x)*/

/* Structire and subroutines for access to tokens within a string */
#define MAXTOKENS 100    /* Maximum number of tokens to parse */
#define MAXWHITE 20     /* Maximum number of whitespace characters */
struct Tokens {
    char *line;         /* Line which has been parsed */
    int lline;          /* Number of characters in line */
    int ntok;           /* Number of tokens on line */
    int nwhite;         /* Number of whitespace characters */
    char white[MAXWHITE];       /* Whitespace (separator) characters */
    char *tok1[MAXTOKENS];      /* Pointers to start of tokens */
    int ltok[MAXTOKENS];        /* Lengths of tokens */
    int itok;           /* Current token number */
};
int setoken();		/* Tokenize a string for easy decoding */
int nextoken();		/* Get next token from tokenized string */
int getoken();		/* Get specified token from tokenized string */

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
 * Oct  1 1999	Add structure and subroutines for tokenized strings
 * Oct 22 1999	Change cat*() to ctg*() to avoid system conflict
 * Oct 29 1999	Add tabget() subroutines
 * Nov  1 1999	Increase maximum number of tokens on a line from 20 to 100
 * Nov  2 1999	Move date utilities to fitsfile.h
 *
 * Jan 10 2000	Add column names to catalog data structure
 * Jan 11 2000	Add gettabndec()
 * Feb  9 2000	Add proper motion entry information to star data structure
 * Feb 16 2000	Add gettaberr() to return tab table error message
 * Mar  1 2000	Add isfile() and agets() to help with ASCII files
 * Mar  8 2000	Add ProgCat() to return catalog name from program name used
 * Mar  8 2000	Add ProgName() to extract program name from path used
 * Mar 10 2000	Add PropCat() to tell whether a catalog has proper motions
 * Mar 27 2000	Add tabxyread()
 */
