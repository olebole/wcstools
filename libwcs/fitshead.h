/* fitshead.h  FITS and IRAF file access subroutines
 * December 12, 1996
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

/* Extract a value from a FITS header for given keyword */
extern int hgeti4 ();	/* int */
extern int hgeti2 ();	/* short */
extern int hgetr4 ();	/* float */
extern int hgetr8 ();	/* double */
extern int hgetra ();	/* Right ascension in degrees from string */
extern int hgetdec ();	/* Declination in degrees from string */
extern int hgetdate (); /* Date in years from FITS date string */
extern int hgetl ();	/* T->1, F->0 from FITS logical entry */
extern int hgets ();	/* Previously allocated string */

/* Implant a value into a FITS header for given keyword */
extern void hputi4 ();	/* int */
extern void hputi2 ();	/* short */
extern void hputr4 ();	/* float */
extern void hputr8 ();	/* double */
extern void hputnr8 ();	/* double with specified number of decimal places */
extern void hputra ();	/* Right ascension in degrees into hh:mm:ss.sss */
extern void hputdec ();	/* Declination in degrees into dd:mm:ss.ss */
extern void hputl ();	/* 0 -> F, else T FITS logical entry */
extern void hputs ();	/* Character string */
extern void hputc ();	/* Character string without quotes */

extern int hdel ();	/* Delete a keyword line from a FITS header */
extern int hchange ();	/* Change a keyword name in a FITS header */

/* Subroutines for conversion between strings and RA and Dec */
extern void ra2str ();
extern void dec2str ();
extern void deg2str ();
extern double str2ra ();
extern double str2dec ();

/* Find given keyword entry in FITS header */
extern char *ksearch ();

/* Check to see whether a string is a number or not */
extern int isnum ();

/* Search for substring s2 within string s1 */
extern char *strsrch ();	/* s1 null-terminated */
extern char *strnsrch ();	/* s1 ls1 characters long */

/* Set length of header which is not null-terminated */
extern int hlength();

/* FITS table keyword structure */
struct Keyword {
    char kname[10];	/* Keyword for table entry */
    int kn;		/* Index of entry on line */
    int kf;		/* Index in line of first character of entry */
    int kl;		/* Length of entry value */
};

#define FITSBLOCK 2880

/* FITS file access subroutines */
extern FILE *fitsropen ();
extern char *fitsrhead ();
extern char *fitsrimage ();
extern int fitswimage ();

/* FITS table file access subroutines */
extern FILE *fitsrtopen ();
extern int fitsrthead ();
extern void fitsrtlset();
extern int fitsrtline ();
extern short ftgeti2 ();
extern int ftgeti4 ();
extern float ftgetr4 ();
extern double ftgetr8 ();
extern int ftgetc ();

/* IRAF file access subroutines */
extern int *irafrhead ();
extern char *irafrimage ();
extern int irafwhead ();
extern int irafwimage ();
extern char *iraf2fits();
extern int fits2iraf();

/* Image pixel access subroutines */
extern double getpix();
extern void putpix();
extern void movepix();
extern void getvec();
extern void putvec();
extern void imswap();
extern void imswap2();
extern void imswap4();
extern void imswap8();
extern int imswapped();

/* Apr 26 1996	Add HGETDATE to get year from date string
 * May 22 1996	Return double from STR2RA and STR2DEC
 * May 31 1996	Use stream I/O for reading as well as writing
 * Jun 12 1996	Add byte-swapping subroutines
 * Jul 10 1996	FITS header now allocated in subroutines
 * Jul 17 1996	Add FITS table column extraction subroutines
 * Jul 19 1996	Add declarations for header implanting subroutines
 * Aug  5 1996	Add HLENGTH for FITS headers which are not null-terminated
 * Aug  5 1996	Add STRNSRCH for FITS headers which are not null-terminated
 * Aug  6 1996	Add HPUTNR8 to save a specified number of decimal places
 * Aug  6 1996	Add MOVEPIX, HDEL and HCHANGE declarations
 * Nov  1 1996	Add DEG2STR
 * Dec 12 1996	Add ISNUM
 */
