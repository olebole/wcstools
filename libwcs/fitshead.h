/* fitshead.h  FITS header value extraction subroutines
 * February 15, 1996
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 */

/* Extract a value from a FITS header for given keyword */
extern int hgeti4 ();
extern int hgeti2 ();
extern int hgetr4 ();
extern int hgetr8 ();

/* Extract RA in degrees from FITS header keyword entry */
extern int hgetra ();

/* Extract Dec in degrees from FITS header keyword entry */
extern int hgetdec ();

/* Extract T->1, F->0 from FITS header logical keyword entry */
extern int hgetl ();

/* Fill previously allocated string from FITS header keyword entry */
extern int hgets ();

/* Subroutines for conversion between strings and RA and Dec */
extern void ra2str ();
extern void dec2str ();
extern void str2ra ();
extern void str2dec ();

/* Find given keyword entry in FITS header */
extern char *ksearch ();

/* Search for substring s2 within string s1 */
extern char *strsrch ();

/* FITS table keyword structure */
struct Keyword {
    char kname[10];	/* Keyword for table entry */
    int kn;		/* Index of entry on line */
    int kf;		/* Index in line of first character of entry */
    int kl;		/* Length of entry value */
};

#define FITSBLOCK 2880

/* FITS file access subroutines */
extern int fitsrtopen ();
extern int fitsropen ();
extern int fitsrtline ();
extern int fitsrhead ();
extern char *fitsrimage ();
extern int fitswimage ();

/* IRAF file access subroutines */
extern int *irafrhead ();
extern char *irafrimage ();
extern int irafwhead ();
extern int irafwimage ();
