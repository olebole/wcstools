/* fitsfile.h  FITS and IRAF file access subroutines
 * December 15, 1999
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

#ifndef fitsfile_h_
#define fitsfile_h_
#include "fitshead.h"

/* FITS table keyword structure */
struct Keyword {
    char kname[10];	/* Keyword for table entry */
    int kn;		/* Index of entry on line */
    int kf;		/* Index in line of first character of entry */
    int kl;		/* Length of entry value */
};

#define FITSBLOCK 2880

/* Subroutines in fitsfile.c */

/* FITS file access subroutines in fitsfile.c */
extern int fitsropen();
extern char *fitsrhead();
extern char *fitsrimage();
extern int fitswhead();
extern int fitswimage();
extern int fitscimage();
extern int isfits();

/* FITS table file access subroutines in fitsfile.c */
extern int fitsrtopen();
extern int fitsrthead();
extern void fitsrtlset();
extern int fitsrtline();
extern short ftgeti2();
extern int ftgeti4();
extern float ftgetr4();
extern double ftgetr8();
extern int ftgetc();

/* IRAF file access subroutines in imhfile.c */
extern char *irafrhead();
extern char *irafrimage();
extern int irafwhead();
extern int irafwimage();
extern int isiraf();
extern char *iraf2fits();
extern char *fits2iraf();

/* Image pixel access subroutines in imio.c */
extern double getpix();	/* Read one pixel from any data type 2-D array (0,0)*/
extern double getpix1(); /* Read one pixel from any data type 2-D array (1,1)*/
extern void putpix();	/* Write one pixel to any data type 2-D array (0,0)*/
extern void putpix1();	/* Write one pixel to any data type 2-D array (1,1) */
extern void addpix();	/* Add to one pixel in any data type 2-D array (0,0)*/
extern void addpix1();	/* Add to one pixel in any data type 2-D array (1,1)*/
extern void movepix();	/* Move one pixel value between two 2-D arrays (0,0) */
extern void movepix1();	/* Move one pixel value between two 2-D arrays (1,1) */
extern void getvec();	/* Read vector from 2-D array */
extern void putvec();	/* Write vector into 2-D array */
extern void imswap();	/* Swap alternating bytes in a vector */
extern void imswap2();	/* Swap bytes in a vector of 2-byte (short) integers */
extern void imswap4();	/* Reverse bytes in a vector of 4-byte numbers */
extern void imswap8();	/* Reverse bytes in a vector of 8-byte numbers */
extern int imswapped();	/* Return 1 if machine byte order is not FITS order */

/* File utilities from fileutil.c */
extern int getfilelines();
extern char *getfilebuff();
extern int getfilesize();
extern int isimlist();
extern int first_token();

/* Subroutines for translating dates and times */
double dt2ep();	/* yyyy.ddmm and hh.mmsss to fractional year (epoch) */
char *dt2fd();	/* yyyy.ddmm and hh.mmsss to FITS date string */
void dt2i();	/* yyyy.ddmm and hh.mmsss to year, month, day, hrs, min, sec */
double dt2jd();	/* yyyy.ddmm and hh.mmsss to Julian date */
double dt2ts();	/* yyyy.ddmm and hh.mmsss to seconds since 1950.0 */ 
void ep2dt();	/* Fractional year to yyyy.mmdd hh.mmssss */
char *ep2fd();	/* Fractional year to FITS date string yyyy-mm-ddThh:mm:ss.ss */
void ep2i();	/* Fractional year to year, month, day, hours, min., sec. */
double ep2jd();	/* Fractional year to Julian Date */
double ep2ts();	/* Fractional year to seconds since 1950.0 */
void fd2dt();	/* FITS standard date string to Julian date */
double fd2ep();	/* FITS standard date string to fractional year (epoch) */
char *fd2fd();	/* Any FITS standard date string to ISO FITS date string */
void fd2i();	/* FITS standard date string to year, mon, day, hrs, min, sec */
double fd2jd();	/* FITS standard date string to Julian date */
double fd2ts();	/* FITS standard date to seconds since 1950.0 */
void jd2dt();	/* Julian date to yyyy.mmdd hh.mmssss */
double jd2ep();	/* Julian date to fractional year */
char *jd2fd();	/* Julian date to FITS date string yyyy-mm-ddThh:mm:ss.ss */
void jd2i();	/* Julian date to year, month, day, hours, min., sec. */
double jd2ts();	/* Julian date to seconds since 1950.0 */
void ts2dt();	/* Seconds since 1950.0 to yyyy.mmdd hh.mmssss */
double ts2ep();	/* Seconds since 1950.0 to fractional year */
char *ts2fd();	/* Seconds since 1950.0 to FITS date, yyyy-mm-ddT00:00:00.000 */
void ts2i();	/* Seconds since 1950.0 to year, month, day, hours, min, sec */
double ts2jd();	/* Seconds since 1950.0 to Julian date */
int isdate();	/* Return 1 if string is FITS old or ISO date */

#endif /* fitsfile_h_ */

/* May 31 1996	Use stream I/O for reading as well as writing
 * Jun 12 1996	Add byte-swapping subroutines
 * Jul 10 1996	FITS header now allocated in subroutines
 * Jul 17 1996	Add FITS table column extraction subroutines
 * Aug  6 1996	Add MOVEPIX, HDEL and HCHANGE declarations
 *
 * Oct 10 1997	FITS file opening subroutines now return int instead of FILE *
 *
 * May 27 1998	Split off fitsio and imhio subroutines to fitsio.h
 * Jun  4 1998	Change fits2iraf from int to int *
 * Jul 24 1998	Make IRAF header char instead of int
 * Aug 18 1998	Change name to fitsfile.h from fitsio.h
 * Oct  5 1998	Add isiraf() and isfits()
 * Oct  7 1998	Note separation of imhfile.c into two files
 *
 * Jul 15 1999	Add fileutil.c subroutines
 * Sep 28 1999	Add (1,1)-based image access subroutines
 * Oct 21 1999	Add fitswhead()
 * Nov  2 1999	Add date utilities from wcscat.h
 * Nov 23 1999	Add fitscimage()
 * Dec 15 1999	Fix misdeclaration of *2fd() subroutines, add fd2i(), dt2i()
 * Dec 20 1999	Add isdate()
 */
