/* libwcs/wcs.h
   November 1, 1996
   By Doug Mink, Harvard-Smithsonian Center for Astrophysics */

struct WorldCoor {
  double	xref;		/* x reference coordinate value (deg) */
  double	yref;		/* y reference coordinate value (deg) */
  double	xrefpix;	/* x reference pixel */
  double	yrefpix;	/* y reference pixel */
  double	xinc;		/* x coordinate increment (deg) */
  double	yinc;		/* y coordinate increment (deg) */
  double	rot;		/* rotation (deg)  (from N through E) */
  double	crot,srot;	/* Cosine and sine of rotation angle */
  double	cd11,cd12,cd21,cd22;
				/* rotation matrix */
  double	dc11,dc12,dc21,dc22;
				/* inverse rotation matrix */
  double	equinox;	/* Equinox of coordinates default to 1950.0 */
  double	epoch;		/* Epoch of coordinates default to equinox */
  double	nxpix;		/* Number of pixels in X-dimension of image */
  double	nypix;		/* Number of pixels in Y-dimension of image */
  double	plate_ra;	/* Right ascension of plate center */
  double	plate_dec;	/* Declination of plate center */
  double	plate_scale;	/* Plate scale in arcsec/mm */
  double	x_pixel_offset;	/* X pixel offset of image lower right */
  double	y_pixel_offset;	/* Y pixel offset of image lower right */
  double	x_pixel_size;	/* X pixel_size */
  double	y_pixel_size;	/* Y pixel_size */
  double	ppo_coeff[6];
  double	amd_x_coeff[20]; /* X coefficients for plate model */
  double	amd_y_coeff[20]; /* Y coefficients for plate model */
  double	xpix;		/* x (RA) coordinate (pixels) */
  double	ypix;		/* y (dec) coordinate (pixels) */
  double	xpos;		/* x (RA) coordinate (deg) */
  double	ypos;		/* y (dec) coordinate (deg) */
  int		pcode;		/* projection code (1-8) */
  int		changesys;	/* 1 for FK4->FK5, 2 for FK5->FK4 */
  				/* 3 for FK4->galactic, 4 for FK5->galactic */
  int		printsys;	/* 1 to print coordinate system, else 0 */
  int		ndec;		/* Number of decimal places in PIX2WCST */
  int		degout;		/* 1 to always print degrees in PIX2WCST */
  int		tabsys;		/* 1 to put tab between RA & Dec, else 0 */
  int		rotmat;		/* 0 if CDELT, CROTA; 1 if CD */
  int		coorflip;	/* 0 if x=RA, y=Dec; 1 if x=Dec, y=RA */
  int		offscl;		/* 0 if OK, 1 if offscale */
  int		plate_fit;	/* 1 if plate fit, else 0 */
  int		wcson;		 /* 1 if WCS is set, else 0 */
  char		c1type[8];	/*  1st coordinate type code:
					RA--, GLON, ELON */
  char		c2type[8];	/*  2nd coordinate type code:
					DEC-, GLAT, ELAT */
  char		ptype[8];	/*  projection type code:
				    -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT */
  char		radecsys[16];	/* Reference frame: FK4, FK4-NO-E, FK5, GAPPT*/
  char		sysout[16];	/* Reference frame for output: FK4, FK5 */
  char		center[32];	/* Center coordinates (with frame) */
  char		search_format[120];	/* search command format */
				/* where %s is replaced by WCS coordinates */
};

#ifndef PI
#define PI		3.141592653589793
#endif

/* Conversions among hours of RA, degrees and radians. */
#define degrad(x)	((x)*PI/180.)
#define raddeg(x)	((x)*180./PI)
#define hrdeg(x)	((x)*15.)
#define deghr(x)	((x)/15.)
#define hrrad(x)	degrad(hrdeg(x))
#define radhr(x)	deghr(raddeg(x))

/* WCS subroutines in wcs.c */
struct WorldCoor *wcsinit (); /* set up a WCS structure from a FITS image header */
struct WorldCoor *wcsninit (); /* set up a WCS structure from a FITS image header */
struct WorldCoor *wcsset (); /* set up a WCS structure */
int iswcs ();		/* Return 1 if WCS structure is filled, else 0 */
int nowcs ();		/* Return 0 if WCS structure is filled, else 1 */
void wcsshift ();	/* Reset the center of a WCS structure */
void wcscent ();
void wcssize ();	/* Return RA and Dec of image center, size in RA and Dec */
void wcsfull ();	/* Return RA and Dec of image center, size in degrees */
double wcsdist ();	/* Distance in degrees between two sky coordinates */
void wcscominit ();	/* Initialize catalog search command set by -wcscom */
void wcscom ();		/* Execute catalog search command set by -wcscom */
void wcsoutinit ();	/* Initialize WCS output coordinate system set by -wcsout */
int pix2wcst ();	/* Convert pixel coordinates to World Coordinate string */
void pix2wcs ();	/* Convert pixel coordinates to World Coordinates */
void wcs2pix ();	/* Convert World Coordinates to pixel coordinates */

/* Oct 26 1994	New file
 * Dec 21 1994	Add rotation matrix
 * Dec 22 1994	Add flag for coordinate reversal

 * Mar  6 1995	Add parameters for Digital Sky Survey plate fit
 * Jun  8 1995	Add parameters for coordinate system change
 * Jun 21 1995	Add parameter for plate scale
 * Jul  6 1995	Add parameter to note whether WCS is set
 * Aug  8 1995	Add parameter to note whether to print coordinate system
 * Oct 16 1995	Add parameters to save image dimensions and center coordinates

 * Feb 15 1996	Add coordinate conversion functions
 * Feb 20 1996	Add flag for tab tables
 * Apr 26 1996	Add epoch of positions (actual date of image)
 * Jul  5 1996	Add subroutine declarations
 * Jul 19 1996	Add WCSFULL declaration
 * Aug  5 1996	Add WCSNINIT to initialize WCS for non-terminated header
 * Oct 31 1996	Add DCnn inverse rotation matrix
 * Nov  1 1996	Add NDEC number of decimal places in output
 */
