/* File libwcs/wcscat.h
 * November 20, 1998
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
#define TABCAT		-1	/* refcat value for StarBase tab table catalog */
#define BINCAT		-2	/* refcat value for TDC binary catalog */
#define TXTCAT		-3	/* refcat value for TDC ASCII catalog */

int RefCat();

/* Subroutines for extracting sources from catalogs by sky region */
int gscread();
int uacread();
int usacread();
int ujcread();
int tabread();
int binread();
int catread();

/* Subroutines for extracting sources from catalogs by ID number */
int gscrnum();
int uacrnum();
int usacrnum();
int ujcrnum();
int tabrnum();
int binrnum();
int catrnum();

/* Subroutines for extracting supplemental information by ID number */
int tabrkey();	/* Keyword values from tab table catalogs */

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
    char objname[32];	/* Object name */
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
    int  insys;		/* Coordinate system
			   B1950 J2000 Galactic Ecliptic */
    double epoch;	/* Epoch of catalog coordinates in years */
    double equinox;	/* Equinox of catalog coordinates in years */
    char inform;	/* Coordinate format
			   (B>inary D>egrees H>MS T>able U>SNO) */
    char incdir[128];	/* Catalog directory pathname */
    char incfile[32];	/* Catalog file name */
    int ncobj;		/* Length of object name in binary star entry */
    char *catbuff;	/* Pointer to start of catalog */
    char *catdata;	/* Pointer to first entry in catalog */
    char *catline;	/* Pointer to current entry in catalog */
    char *catlast;	/* Pointer to one past end of last entry in catalog */
    int  istar;		/* Number of current catalog entry */
};

/* Subroutines for reading headers of TDC binary and ASCII catalogs */
int istab();
int isbin();
struct StarCat *binopen();
void binclose();
struct StarCat *catopen();
void catclose();

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
 */
