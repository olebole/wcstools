
struct Star {
    float rdum;
    float xno;			/* Catalog number */
    double sra0;		/* Right Ascension (radians) */
    double sdec0;		/* Declination (radians) */
    char isp[2];		/* Spectral type or other 2-char identifier */
    short mag[11];		/* Up to 10 Magnitudes * 100 */
    float xrpm;			/* RA proper motion (radians per year) */
    float xdpm;			/* Dec proper motion (radians per year) */
    double dcra;		/* Change in RA from catalog position */
    double dcdec;		/* Change in Dec from catalog position */
    int ixno;			/* Integer star number */
    char objname[32];		/* Object name */
    }

struct StarCat {
    int star0;		/* Subtract from star number for file sequence number */
    int star1;		/* First star number in file */
    int starn;		/* Number of stars in file */
    int stnum;		/* Star number format in catalog file:
			   0:  no star i.d. numbers
			   1: Real*4 star i.d. numbers
			   2: Integer*4 <region><nnnn>
			   3: Integer*4 <region><nnnnn>
			   4: Integer*4 <nnnnnnnnn> */
    int mprop;		/* 1 if proper motion is included */
    int nmag;		/* Number of magnitudes present
			   Negative for J2000 catalog */
    int nbent;		/* Number of bytes per star entry */
    int ifcat;		/* File descriptor for catalog file */
    char isfil[24];	/* Star catalog file name */
    char isname[64];	/* Star catalog description */
    char inform;	/* Coordinate format
			   (B>inary D>EGREES H>MS T>ABLE U>SNO) */
    char incoor;	/* Coordinate system
			   B>1950 J>2000 G>alactic E>cliptic */
    char incdir[128];	/* Catalog directory pathname */
    char incfile[32];	/* Catalog file name */
    }
