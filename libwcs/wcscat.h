/* File libwcs/wcscat.h
 * September 25, 2001
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
#define USNO		17	/* USNO-format plate catalog */
#define TMPSC		18	/* 2MASS Point Source Catalog */
#define GSCACT		19	/* GSC-ACT revised Guide Star Catalog */
#define GSC2		20	/* GSC II version 2.2 */
#define TABCAT		-1	/* StarBase tab table catalog */
#define BINCAT		-2	/* TDC binary catalog */
#define TXTCAT		-3	/* TDC ASCII catalog */
#define WEBCAT		-4	/* Tab catalog via the web */
#define NUMCAT		19	/* Number of predefined catalogs */

/* Subroutines for dealing with catalogs */
int RefCat();		/* Return catalog type code, title, coord. system */
char *CatName();	/* Return catalog name given catalog type code */
char *ProgCat();	/* Return catalog name given program name used */
char *ProgName();	/* Return program name given program path used */
char *CatName();	/* Return catalog name given catalog type code */
void CatID();		/* Return catalog ID keyword given catalog type code */
void CatNum();		/* Return formatted source number */
int CatNumLen();	/* Return length of source numbers */
int CatNdec();		/* Return number of decimal places in source numbers */
void CatMagName();	/* Return name of specified magnitude */

int StrNdec();		/* Return number of decimal places in numeric string */
void SearchLim();	/* Compute limiting RA and Dec */
void RefLim();		/* Compute limiting RA and Dec in new system */
int isacat();		/* Return 1 if string is name of ASCII catalog file */
int ageti4();		/* Extract int value from keyword= value in string */
int agetr8();		/* Extract double value from keyword= value in string */
int agets();		/* Extract value from keyword= value in string */
void bv2sp();		/* Approximate main sequence spectral type from B - V */

/* Subroutines for extracting sources from catalogs by sky region */
int gscread();		/* Read sources from HST Guide Star Catalog */
int tmcread();		/* Read sources from 2MASS Point Source Catalog */
int uacread();		/* Read sources from USNO A or SA Catalog */
int ujcread();		/* Read sources from USNO J Catalog */
int tabread();		/* Read sources from tab table catalog */
int binread();		/* Read sources from SAO TDC binary format catalog */
int ctgread();		/* Read sources from SAO TDC ASCII format catalog */
int actread();		/* Read sources from USNO ACT Catalog */
int ty2read();		/* Read sources from Tycho 2 Catalog */
int webread();		/* Read sources from catalog on the World Wide Web */


/* Subroutines for extracting sources from catalogs by ID number */
int gscrnum();		/* Read sources from HST Guide Star Catalog */
int tmcrnum();		/* Read sources from 2MASS Point Source Catalog */
int uacrnum();		/* Read sources from USNO A or SA Catalog */
int ujcrnum();		/* Read sources from USNO J Catalog */
int tabrnum();		/* Read sources from tab table catalog */
int binrnum();		/* Read sources from SAO TDC binary format catalog */
int ctgrnum();		/* Read sources from SAO TDC ASCII format catalog */
int actrnum();		/* Read sources from USNO ACT Catalog */
int ty2rnum();		/* Read sources from Tycho 2 Catalog */
int webrnum();		/* Read sources from catalog on the World Wide Web */

void setgsclass();	/* Set GSC object class */
void setuplate();	/* Set USNO catalog plate number to search */
int getuplate();	/* Get USNO catalog plate number to search */
void settabkey();	/* Set tab table keyword to read for object */
int ctgstar();		/* Read one star entry from ASCII catalog, 0 if OK */
int binstar();		/* Read one star entry from binary catalog, 0 if OK */
int tabstar();		/* Read one star entry from tab table catalog, 0 if OK */
struct TabTable *webopen();	/* Open tab table across the web */
char *webbuff();	/* Read URL into buffer across the web */

/* Subroutines for sorting tables of star positions and magnitudes */
#define NOSORT		0	/* Do not sort catalog output */
#define SORT_MAG	1	/* Sort output by magnitude */
#define SORT_DIST	2	/* Sort output by distance from center */
#define SORT_RA		3	/* Sort output by right ascension */
#define SORT_DEC	4	/* Sort output by declination */
#define SORT_X		5	/* Sort output by image X coordinate */
#define SORT_Y		6	/* Sort output by image Y coordinate */
void XSortStars();
void YSortStars();
void RASortStars();
void DecSortStars();
void MagSortStars();

/* Data structure for SAO TDC ASCII and binary star catalog entries */
struct Star {
    float rdum;
    float xno;		/* Catalog number */
    double ra;		/* Right Ascension (degrees) */
    double dec;		/* Declination (degrees) */
    char isp[24];	/* Spectral type or other 2-char identifier */
    short mag[11];	/* Up to 10 Magnitudes * 100 */
    double rapm;	/* RA proper motion (degrees per year) */
    double decpm;	/* Dec proper motion (degrees per year) */
    double xmag[11];	/* Up to 10 Magnitudes */
    double num;		/* Actual star number */
    int coorsys;	/* Coordinate system (WCS_J2000, WCS_B1950,...) */
    double equinox;	/* Equinox of coordinate system as fractional year */
    double epoch;	/* Epoch of position as fractional year */
    double parallax;	/* Parallax in arcseconds */
    double pxerror;	/* Parallax error in arcseconds */
    double radvel;	/* Radial velocity in km/sec, positive away */
    double dist;	/* Distance from search center in arcseconds */
    char *entry;	/* Line copied from input catalog */
    char objname[32];	/* Object name */
    int peak;		/* Peak flux per pixel in star image */
};

/* Catalog proper motion units */
#define PM_MASYR		1	/* milliarcseconds per year */
#define PM_ARCSECYR		2	/* arcseconds per year */
#define PM_DEGYR		3	/* degrees per year */
#define PM_RADYR		4	/* radians per year */
#define PM_TSECYR		5	/* seconds of time (RA) per century */
#define PM_ARCSECCEN		6	/* arcseconds per year */
#define PM_TSECCEN		7	/* seconds of time (RA) per century */
#define PM_MTSYR		8	/* milliseconds of time (RA) per year */

/* Data structure for SAO TDC ASCII and binary star catalogs */
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
			/* 2 if radial velocity is included */
    int nmag;		/* Number of magnitudes present
			   Negative for J2000 catalog */
    int nbent;		/* Number of bytes per star entry */
    int	rasorted;	/* 1 if RA-sorted, else 0 */
    int	ignore;		/* 1 if ignoring info after position and magnitude */
    FILE *ifcat;	/* File descriptor for catalog file */
    char isfil[24];	/* Star catalog file name */
    char isname[64];	/* Star catalog description */
    int  byteswapped;	/* 1 if catalog is byte-reversed from CPU */
    int  refcat;	/* Code for type of catalog (TXTCAT, BINCAT, etc.) */
    int  coorsys;	/* Coordinate system
			   B1950 J2000 Galactic Ecliptic */
    double epoch;	/* Epoch of catalog coordinates in years */
    double equinox;	/* Equinox of catalog coordinates in years */
    char inform;	/* Coordinate format
			   (B>inary D>egrees H>MS T>able U>SNO) */
    char incdir[128];	/* Catalog directory pathname */
    char incfile[32];	/* Catalog file name */
    int ncobj;		/* Length of object name in binary star entry */
    int nnfld;		/* Length of star number  */
    int nndec;		/* Number of decimal places in star number */
    int nepoch;		/* 1 if epoch of coordinates is present */
    int sptype;		/* 1 if spectral type is present in catalog */
    int plate;		/* 1 if plate or field number is present in catalog */
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
    int entmag3;	/* Entry number for third magnitude, if present */
    int entmag4;	/* Entry number for fourth magnitude, if present */
    int entpeak;	/* Entry number for peak counts */
    int entepoch;	/* Entry number for epoch of observation */
    int entname;	/* Entry number for object name */
    int entadd;		/* Entry number for additional keyword */
    int entrpm;		/* Entry number for proper motion in right ascension */
    int entdpm;		/* Entry number for proper motion in declination */
    int entpx;		/* Entry number for parallax */
    int entpxe;		/* Entry number for parallax error */
    int entrv;		/* Entry number for radial velocity */
    int enttype;	/* Entry number for spectral type */
    int rpmunit;	/* Units for RA proper motion (PM_x) */
    int dpmunit;	/* Units for DEC proper motion (PM_x) */
    char *caturl;	/* set if web search, else NULL */
    char keyid[16];	/* Entry name for ID */
    char keyra[16];	/* Entry name for right ascension */
    char keydec[16];	/* Entry name for declination */
    char keymag1[16];	/* Entry name for first or only magnitude */
    char keymag2[16];	/* Entry name for second magnitude, if present */
    char keymag3[16];	/* Entry name for third magnitude, if present */
    char keymag4[16];	/* Entry name for fourth magnitude, if present */
    char keyrpm[16];	/* Entry name for right ascension proper motion */
    char keydpm[16];	/* Entry name for declination proper motion */
    char keypeak[16];	/* Entry name for integer code */
    char keytype[16];	/* Entry name for spectral type */
    char keyrv[16];	/* Entry name for radial velocity */
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
    char *tabname;	/* Name of this table or NULL */
    char *tabbuff;	/* Pointer to start of saved tab table in memory */
    char *tabheader;	/* Pointer to start of line containing table header */
    char *tabhead;	/* Pointer to start of line containing column heading */
    char *tabdata;	/* Pointer to start of first line of table data */
    int lhead;		/* Number of bytes before first data line */
    int iline;		/* Number of current line (1=first) */
    int lline;		/* Length in bytes of line buffer */
    char *tabline;	/* Pointer to start of current line */
    FILE *tcat;		/* File descriptor for tab table file */
    int ncols;		/* Number of columns per table entry */
    char **colname;	/* Column names */
    int *lcol;		/* Lengths of column header names */
    int *lcfld;		/* Number of columns in field (hyphens) */
    int lbuff;		/* Number of bytes in entire tab table */
};

/* Subroutines for extracting tab table information */
struct TabTable *tabopen();	/* Open tab table file */
struct StarCat *tabcatopen();	/* Open tab table catalog */
void tabcatclose();	/* Close tab table catalog */
int tabxyread();	/* Read x, y, and magnitude from tab table star list */
char *gettabline();	/* Find a specified line in a tab table */
int tabrkey();		/* Keyword values from tab table catalogs */
int tabcol();		/* Find column for name */
int tabgetk();		/* Get tab table entries for named column */
int tabgetc();		/* Get tab table entry for named column */
int tabgeti4();		/* Return 4-byte integer from tab table line */
int tabparse();		/* Aeturn column names and positions in tabtable */
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
#define WCS_CIRCLE 1	/* circle shape for SAOimage plotting */
#define WCS_SQUARE 2	/* square shape for SAOimage plotting */
#define WCS_DIAMOND 3	/* diamond shape for SAOimage plotting */
#define WCS_CROSS 4	/* cross shape for SAOimage plotting */
#define WCS_EX 5	/* x shape for SAOimage plotting */
#define WCS_VAR 6	/* variable (+ and x) shape for HSTGSC plotting */
#define WCS_PCIRCLE 11	/* pixel circle shape for SAOimage plotting */
#define WCS_PSQUARE 12	/* pixel square shape for SAOimage plotting */
#define WCS_PDIAMOND 13	/* pixel diamond shape for SAOimage plotting */
#define WCS_PCROSS 14	/* pixel cross shape for SAOimage plotting */
#define WCS_PEX 15	/* pixel ex shape for SAOimage plotting */
#define WCS_PVAR 16	/* pixel variable (+ and x) shape for HSTGSC plotting */

/* Structure and subroutines for access to tokens within a string */
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

/* Subroutines for fitting and evaluating polynomials */
void polfit();		/* Fit polynomial coefficients */
double polcomp();	/* Evaluate polynomial from polfit coefficients */

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
 * Apr  3 2000	Add option in catalog structure to ignore extra info
 * May 22 2000	Add Tycho 2 support, bv2sp()
 * May 26 2000	Add separate pointer to header in tab table structure
 * May 26 2000	Add separate pointer to table name in tab table structure
 * Jul 12 2000	Add catalog type code to ctalog data structure
 * Sep 20 2000	Add isacat() to detect ASCII catalog files
 * Sep 25 2000	Add starcat.sptype to flag spectral type in catalog
 * Oct 23 2000	Add USNO plate catalog to catalog type table
 * Oct 26 2000	Add proper motion flags for seconds and arcseconds per century
 * Oct 31 2000	Add proper motion flags for milliseconds per year
 * Nov  2 2000	Add parallax and radial velocity to star structure
 * Nov 21 2000	Add WEBCAT catalog type for tab ctalogs returned from the Web
 * Nov 22 2000	Add webread() and webrnum()
 * Nov 28 2000	Add tabparse()
 * Nov 30 2000	Add spectral type to catalog header; make star->isp 4 char.
 * Dec 13 2000	Add StrNdec() to get number of decimal places in number strings
 * Dec 15 2000	Add CatNdec() to get number of decimal places in source numbers
 * Dec 18 2000	Drop PropCat(), a cludgy proper motion flag
 *
 * Mar 22 2001	Add web search flag in catalog data structure
 * Mar 27 2001	Add shapes in pixels to SAOimage region options
 * May 14 2001	Add 2MASS Point Source Catalog flags
 * May 22 2001	Add declination sorting
 * May 24 2001	Add 2MASS Point Source Catalog subroutines
 * May 29 2001	Add length of star number to catalog structure
 * May 30 2001	Add third magnitude for tab tables to catalog structure
 * Jun 15 2001	Add CatName() and CatID()
 * Jun 19 2001	Add parallax error to catalog and star structures
 * Jun 20 2001	Add webopen(), GSC2, fourth magnitude to star and starcat
 * Jul 12 2001	Add separate web access subroutine, webbuff()
 * Jul 23 2001	Add ageti4() and agetr8()
 * Jul 24 2001	Add polfit() and polcomp()
 * Aug  8 2001	Add keyrv and option to set mprop to 2 to include rv/cz
 * Sep 10 2001	Add entry line and distance from search center to Star
 * Sep 13 2001	Add YSortStars() and SORT_Y
 * Sep 14 2001	Add lbuff to TabTable structure
 * Sep 20 2001	Add CatMagName()
 * Sep 25 2001	Move isfile() to fitsfile.h
 */
