/* File irafio.c
 * August 13, 1996
 * By Doug Mink, based on Mike VanHilst's readiraf.c

 * Module:      irafio.c (IRAF image file reading and writing)
 * Purpose:     Read and write IRAF image files (and translate headers)
 * Subroutine:  check_immagic (irafheader, teststring )
 *		Verify that file is valid IRAF imhdr or impix
 * Subroutine:  irafrhead (filename, lfhead, fitsheader, lihead)
 *              Read IRAF image header
 * Subroutine:  irafrimage (fitsheader)
 *              Read IRAF image pixels (call after irafrhead)
 * Subroutine:	same_path (pixname, hdrname)
 *		Put filename and header path together
 * Subroutine:	iraf2fits (irafheader, nbiraf, nbfits)
 *		Convert IRAF image header to FITS image header
 * Subroutine:	irafwhead (hdrname, irafheader, fitsheader)
 *		Write IRAF header file
 * Subroutine:	irafwimage (hdrname, irafheader, fitsheader, image )
 *		Write IRAF image and header files
 * Subroutine:	fits2iraf (fitsheader, irafheader)
 *		Convert FITS image header to IRAF image header
 * Subroutine:	irafswap (bitpix,string,nbytes)
 *		Swap bytes in string in place, with FITS bits/pixel code
 * Subroutine:	irafswap2 (string,nbytes)
 *		Swap bytes in string in place
 * Subroutine	irafswap4 (string,nbytes)
 *		Reverse bytes of Integer*4 or Real*4 vector in place
 * Subroutine	irafswap8 (string,nbytes)
 *		Reverse bytes of Real*8 vector in place


 * Copyright:   1996 Smithsonian Astrophysical Observatory
 *              You may do anything you like with this file except remove
 *              this copyright.  The Smithsonian Astrophysical Observatory
 *              makes no representations about the suitability of this
 *              software for any purpose.  It is provided "as is" without
 *              express or implied warranty.
 */

#include <stdio.h>		/* define stderr, FD, and NULL */
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <sys/stat.h>
#include "fitshead.h"


/* Parameters from iraf/lib/imhdr.h */
#define SZ_IMPIXFILE	79
#define LEN_IMHDR	513

/* Offsets into header in sizeof(int) units for various parameters */
#define	IM_PIXTYPE	4		/* datatype of the pixels */
#define	IM_NDIM		5		/* number of dimensions */
#define	IM_LEN		6		/* length (as stored) */
#define	IM_PHYSLEN	13		/* physical length (as stored) */
#define	IM_PIXOFF	22		/* offset of the pixels */
#define	IM_CTRAN	52		/* coordinate transformations */
#define	IM_PIXFILE	103		/* name of pixel storage file */
#define	IM_TITLE	183		/* image name string */

/* The Coordinate Transformation Structure (IM_CTRAN) from iraf/lib/imhdr.h */
/* None of this stuff is supported by IRAF applications, yet */
#define	LEN_CTSTRUCT	50
#define	CT_VALID	0		/* (y/n) is structure valid? */
#define	CT_BSCALE	1		/* pixval scale factor */
#define	CT_BZERO	2		/* pixval offset */
#define	CT_CRVAL	3		/* value at pixel (for each dim) */
#define	CT_CRPIX	10		/* index of pixel (for each dim) */
#define	CT_CDELT	17		/* increment along axis for each dim */
#define	CT_CROTA	24		/* rotation angle (for each dim) */
#define	CT_CTYPE	36		/* coord units string */

/* Codes from iraf/unix/hlib/iraf.h */
#define	TY_CHAR		2
#define	TY_SHORT	3
#define	TY_INT		4
#define	TY_LONG		5
#define	TY_REAL		6
#define	TY_DOUBLE	7
#define	TY_COMPLEX	8

#define LEN_IRAFHDR	25000
#define LEN_PIXHDR	1024
#define LEN_FITSHDR	11520

int check_immagic();
char *irafstr();
static int swapiraf=0;	/* =1 to swap bytes of foreign IRAF file */
static void irafswap();
static void irafswap2();
static void irafswap4();
static void irafswap8();

/* Subroutine:	irafrhead
 * Purpose:	Open and read the iraf .imh file, translating it to FITS, too.
 * Returns:	NULL if failure, else pointer to IRAF .imh image header
 * Notes:	The imhdr format is defined in iraf/lib/imhdr.h, some of
 *		which defines or mimicked, above.
 */

int *
irafrhead (filename, lihead)

char	*filename;	/* Name of IRAF header file */
int	*lihead;	/* Length of IRAF image header in bytes (returned) */
{
    FILE *fd;
    int nbr;
    int *irafheader;
    int nbhead;
    struct stat buff;

    /* Find size of image header file */
    if (stat (filename,&buff)) {
	fprintf (stderr, "IRAFRHEAD:  cannot read file %s\n", filename);
	perror ("IRAFRHEAD");
	return (NULL);
	}
    nbhead = (int) buff.st_size;
    *lihead = 0;

    /* open the image header file */
    fd = fopen (filename, "r");
    if (!fd) {
	fprintf (stderr, "IRAFRHEAD:  cannot read file %s\n", filename);
	perror ("IRAFRHEAD");
	return (NULL);
	}

    /* allocate initial sized buffer */
    *lihead = nbhead + 500;
    irafheader = (int *)malloc ((unsigned int) *lihead);

    /* read in IRAF header */
    nbr = fread (irafheader, 1, nbhead, fd);
    fclose (fd);
    if (nbr < LEN_PIXHDR) {
	(void)fprintf(stderr, "IRAF header file %s: %d / %d bytes read.\n",
		      filename,nbr,LEN_PIXHDR);
	free ((char *)irafheader);
	return (NULL);
	}

    /* check header magic word */
    if (check_immagic ((short *)irafheader, "imhdr") ) {
	free ((char *)irafheader);
	(void)fprintf(stderr, "File %s not valid IRAF image header\n", filename);
	return(NULL);
	}

    /* check number of image dimensions */
    if (irafheader[IM_NDIM]< 2 ) {
	free ((char *)irafheader);
	(void)fprintf(stderr, "File %s does not contain 2d image\n", filename);
	return (NULL);
	}

    return (irafheader);
}


char *
irafrimage (fitsheader)

char	*fitsheader;	/* FITS image header (filled) */
{
    FILE *fd;
    char *bang;
    int naxis1, naxis2, bitpix, bytepix;
    char *pixname;
    char *image;
    int nbr, nbimage;
    int *pixheader;

    /* Convert pixel file name to character string */
    pixname = (char *) malloc (SZ_IMPIXFILE);
    hgets (fitsheader, "PIXFILE", pixname);

    if ((bang = strchr (pixname, '!')) != NULL ) {
	fd = fopen (bang + 1, "r");
	if (!fd) {
	    (void)fprintf(stderr, "IRAFRIMAGE: Cannot open IRAF pixel file %s\n", pixname);
	    perror ("IRAFRIMAGE");
	    return (NULL);
	    }
	}
    else {
	fd = fopen (pixname, "r");
	if (!fd) {
	    (void)fprintf(stderr, "Cannot open IRAF pixel file %s\n", pixname);
	    perror ("IRAFRIMAGE");
	    return (NULL);
	    }
	}

    /* Read pixel header */
    pixheader = (int *)malloc (LEN_PIXHDR);
    nbr = fread (pixheader, 1, LEN_PIXHDR, fd);

    /* Check size of pixel header */
    if (nbr < LEN_PIXHDR) {
	(void)fprintf(stderr, "IRAF pixel file %s: %d / %d bytes read.\n",
		      pixname,nbr,LEN_PIXHDR);
	free (pixheader);
	fclose (fd);
	return (NULL);
	}

    /* check pixel header magic word */
    if (check_immagic ((short *)pixheader, "impix") ) {
	(void)fprintf(stderr, "File %s not valid IRAF pixel file.\n", pixname);
	free (pixheader);
	fclose (fd);
	return(NULL);
	}
    free (pixheader);

    /* Find number of bytes to read */
    hgeti4 (fitsheader,"NAXIS1",&naxis1);
    hgeti4 (fitsheader,"NAXIS2",&naxis2);
    hgeti4 (fitsheader,"BITPIX",&bitpix);
    if (bitpix < 0)
	bytepix = -bitpix / 8;
    else
	bytepix = bitpix / 8;
    nbimage = naxis1 * naxis2 * bytepix;
    image =  (char *) malloc (nbimage);

    /* read in IRAF image */
    nbr = fread (image, 1, nbimage, fd);
    fclose (fd);

    /* Check size of image */
    if (nbr < nbimage) {
	(void)fprintf(stderr, "IRAF pixel file %s: %d / %d bytes read.\n",
		      pixname,nbr,nbimage);
	free (image);
	return (NULL);
	}

	    /* Byte-reverse image, if necessary */
    if (swapiraf)
	irafswap (bitpix, image, nbimage);

    return (image);
}


/* Verify that file is valid IRAF imhdr or impix by checking first 5 chars
 * Returns:	0 on success, 1 on failure */

int
check_immagic (irafheader, teststring )

short	*irafheader;	/* IRAF image header from file */
char	*teststring;	/* "imhdr" for header file or "impix" for pixel file */

{
    char *line;

    line = irafstr (irafheader,5,0);
    if (strncmp(line, teststring, 5) == 0) {
	free (line);
	swapiraf = 0;
	return (0);
	}
    else {
	free (line);
	line = irafstr (irafheader,5,1);
	if (strncmp(line, teststring, 5) == 0) {
	    free (line);
	    swapiraf = 1;
	    return (0);
	    }
	free (line);
	return (1);
	}
}


/* Convert IRAF 2-byte/char string to 1-byte/char string */

char *
irafstr (irafstring, nchar, byteorder)

short	*irafstring;	/* IRAF 2-byte/character string */
int	nchar;		/* Number of characters in string */
int	byteorder;	/* 0=big-endian (Sun, Mac), 1=little-endian (PC,DEC) */
{
    char *string;
    int i, nbytes;

    string = (char *) malloc (nchar+1);

    /* Swap bytes, if requested */
    if (byteorder) {
	nbytes = nchar * 2;
	irafswap2 (irafstring, nbytes);
	}

    /* Convert appropriate byte of input to output character */
    for (i = 0; i < nchar; i++)
	string[i] = (char) irafstring[i];

    return (string);
}


/* Put filename and header path together */

void
same_path (pixname, hdrname)

char	*pixname;	/* IRAF pixel file pathname */
char	*hdrname;	/* IRAF image header file pathname */

{
    int len;
    char temp[SZ_IMPIXFILE];

    if (strncmp(pixname, "HDR$", 4) == 0 ) {

	/* load entire header name string into name buffer */
	(void)strncpy (temp, &pixname[4], SZ_IMPIXFILE);
	(void)strncpy (pixname, hdrname, SZ_IMPIXFILE);

	/* find the end of the pathname */
	len = strlen(pixname);
#ifndef VMS
	while( (len > 0) && (pixname[len-1] != '/') )
#else
	while( (len > 0) && (pixname[len-1] != ']') && (pixname[len-1] != ':') )
#endif
      len--;

	/* add name */
	pixname[len] = '\0';
	(void)strncat(pixname, temp, SZ_IMPIXFILE);
	}
}

/* Convert IRAF image header to FITS image header, returning FITS header */

char *
iraf2fits (hdrname, irafheader, nbiraf, nbfits)

char	*hdrname;	/* IRAF header file name (may be path) */
int	*irafheader;	/* IRAF image header */
int	nbiraf;		/* Number of bytes in IRAF image */
int	*nbfits;	/* Number of bytes in FITS header (returned) */

{
    int lfhead;		/* Actual length of FITS header (returned) */
    char *objname;	/* object name from FITS file */
    int i, j, nax, nbits, nbytes;
    short *irafobj;
    char *pixname;
    short *irafpix;
    char *fitsheader;

    int ncr, nblock, nlines;
    char *fhead, *fhead1, *fp, endline[81];
    short *irafline;
    char fitsline[81];

    (void)strncpy (endline,"END                                     ",40);
    (void)strncpy (endline+40,"                                        ",40);

    /*  Initialize FITS header */
    nlines = 7 + (nbiraf - (4 * LEN_IMHDR) / 162);
    nblock = (nlines * 80) / 2880;
    *nbfits = (nblock + 3) * 2880;
    fitsheader = (char *) malloc ((unsigned int) *nbfits);
    fhead = fitsheader;
    lfhead = 0;
    (void)strncpy (fitsheader, endline, 80);
    hputl (fitsheader, "SIMPLE", 1);
    fhead = fhead + 80;
    lfhead = lfhead + 80;

    /*  Set pixel size in FITS header */
    if (swapiraf)
	irafswap4 ((char *)&irafheader[IM_PIXTYPE],4);
    switch (irafheader[IM_PIXTYPE] ) {
    case TY_CHAR:
	nbits = 8;
	break;
    case TY_SHORT:
	nbits = 16;
	break;
    case TY_INT:
    case TY_LONG:
	nbits = 32;
	break;
    case TY_REAL:
	nbits = -32;
	break;
    case TY_DOUBLE:
	nbits = -64;
	break;
    default:
	(void)fprintf(stderr,"Unsupported data type: %d\n",
		      irafheader[IM_PIXTYPE]);
	return (NULL);
    }
    hputi4 (fitsheader,"BITPIX",nbits);
    fhead = fhead + 80;
    lfhead = lfhead + 80;

    /*  Set image dimensions in FITS header */
   if (swapiraf)
	irafswap4 ((char *)&irafheader[IM_NDIM],4);
    nax = irafheader[IM_NDIM];
    hputi4 (fitsheader,"NAXIS",nax);
    fhead = fhead + 80;
    lfhead = lfhead + 80;

    if (swapiraf)
	irafswap4 ((char *)&irafheader[IM_PHYSLEN],4);
    hputi4 (fitsheader,"NAXIS1",irafheader[IM_PHYSLEN]);
    fhead = fhead + 80;
    lfhead = lfhead + 80;

    if (swapiraf)
	irafswap4 ((char *)&irafheader[IM_PHYSLEN+1],4);
    hputi4 (fitsheader,"NAXIS2",irafheader[IM_PHYSLEN+1]);
    fhead = fhead + 80;
    lfhead = lfhead + 80;

    if (nax > 2) {
	if (swapiraf)
	    irafswap4 ((char *)&irafheader[IM_PHYSLEN+2],4);
	hputi4 (fitsheader,"NAXIS3",irafheader[IM_PHYSLEN+2]);
	fhead = fhead + 80;
	lfhead = lfhead + 80;
	}
    if (nax > 3) {
	hputi4 (fitsheader,"NAXIS4",irafheader[IM_PHYSLEN+3]);
	fhead = fhead + 80;
	lfhead = lfhead + 80;
	}

    /* Set object name in FITS header */
    irafobj = (short *) (irafheader + IM_TITLE);
    objname = irafstr (irafobj, 64, swapiraf);
    hputs (fitsheader,"OBJECT",objname);
    free (objname);
    fhead = fhead + 80;
    lfhead = lfhead + 80;

    /* Save image header filename in header */
    hputs (fitsheader,"IMHFILE",hdrname);
    fhead = fhead + 80;
    lfhead = lfhead + 80;

    /* Save image pixel file pathname in header */
    irafpix = (short *) (irafheader + IM_PIXFILE);
    pixname = irafstr (irafpix, SZ_IMPIXFILE, swapiraf);
    if (strncmp(pixname, "HDR$", 4) == 0 )
	same_path (pixname, hdrname);
    hputs (fitsheader,"PIXFILE",pixname);
    fhead = fhead + 80;
    lfhead = lfhead + 80;

    /* Swap user portion of IRAF header, if necessary */
    if (swapiraf) {
	nbytes = nbiraf - (LEN_IMHDR * 4);
	fhead1 = ((char *) irafheader) + LEN_IMHDR*4;
	irafswap2 (fhead1, nbytes);
	}

    /* Add user portion of IRAF header to FITS header */
    ncr = nbiraf / 2;
    fitsline[80] = 0;
    for (i = LEN_IMHDR*2; i < ncr; i = i + 81) {
	irafline = ((short *) irafheader) + i;
	if (irafline[0] == 0) break;
	for (j = 0; j < 80; j++) {
	    if (irafline[j] < 32)
		irafline[j] = 32;
	    fitsline[j] = (char) irafline[j];
	    }
	(void)strncpy (fhead, fitsline, 80);
	fhead = fhead + 80;
	/* printf ("%80s\n",fitsline); */
	lfhead = lfhead + 80;
	}

    /* Add END to last line */
    (void)strncpy (fhead, endline, 80);
    lfhead = lfhead + 80;

    /* Find end of last 2880-byte block of header */
    nblock = lfhead / 2880;
    if (nblock*2880 < lfhead)
	nblock = nblock + 1;
    fhead1 = ksearch (fitsheader, "END") + 80;

    /* Pad rest of header with spaces */
    strncpy (endline,"   ",3);
    for (fp = fhead+80; fp < fhead1; fp = fp + 80) {
	(void)strncpy (fp, endline,80);
	lfhead = lfhead + 80;
	}

    return (fitsheader);
}


int
irafwhead (hdrname, irafheader, fitsheader)

char	*hdrname;	/* Name of IRAF header file */
int	*irafheader;	/* IRAF header */
char	*fitsheader;	/* FITS image header */

{
    FILE *fd;
    int nbw, nbhead;

    /* Convert FITS header to IRAF header */
    nbhead = fits2iraf (fitsheader, irafheader);

    /* Write IRAF header to disk file */
    fd = fopen (hdrname, "w");
    if (!fd) {
	(void)fprintf(stderr, "IRAFWHEAD: Cannot open IRAF header file %s\n", hdrname);
	perror ("IRAFWHEAD");
	return (0);
	}
    nbw = fwrite (irafheader, 1, nbhead, fd);
    fclose (fd);
    if (nbw < nbhead) {
	(void)fprintf(stderr, "IRAF header file %s: %d / %d bytes written.\n",
		      hdrname, nbw, nbhead);
	return (-1);
	}

    return (nbw);
}


int
irafwimage (hdrname, irafheader, fitsheader, image )

char	*hdrname;	/* Name of IRAF header file */
int	*irafheader;	/* IRAF header */
char	*fitsheader;	/* FITS image header */
char	*image;		/* IRAF image */

{
    FILE *fd;
    char *bang;
    int i, nbw, bytepix, bitpix, naxis1, naxis2, nbhead, nbimage, ioff;
    void same_path();
    char *pixname;
    short *irafpix;

    /* Convert FITS header to IRAF header */
    nbhead = fits2iraf (fitsheader, irafheader);

    /* Write IRAF header to disk file */
    fd = fopen (hdrname, "w");
    if (!fd) {
	(void)fprintf(stderr, "IRAFWIMAGE: Cannot open IRAF header file %s\n", hdrname);
	perror ("IRAFWIMAGE");
	return (0);
	}
    nbw = fwrite (irafheader, 1, nbhead, fd);
    fclose (fd);

    /* Convert pixel file name to character string */
    irafpix = (short *) (irafheader + IM_PIXFILE);
    pixname = (char *) malloc (SZ_IMPIXFILE);
    for (i = 0; i < SZ_IMPIXFILE; i++)
	pixname[i] = (char) irafpix[i];
    if (strncmp(pixname, "HDR$", 4) == 0 )
	same_path (pixname, hdrname);

    /* Find number of bytes to write */
    hgeti4 (fitsheader,"NAXIS1",&naxis1);
    hgeti4 (fitsheader,"NAXIS2",&naxis2);
    hgeti4 (fitsheader,"BITPIX",&bitpix);
    if (bitpix < 0)
	bytepix = -bitpix / 8;
    else
	bytepix = bitpix / 8;
    nbimage = naxis1 * naxis2 * bytepix;

    /* Open pixel file */
    if ((bang = strchr (pixname, '!')) != NULL ) {
	fd = fopen (bang + 1, "w");
	if (!fd) {
	    (void)fprintf(stderr, "Cannot open IRAF pixel file %s\n", pixname);
	    perror ("IRAFWIMAGE");
	    return (0);
	    }
	}
    else {
	fd = fopen (pixname, "w");
	if (!fd) {
	    (void)fprintf(stderr, "Cannot open IRAF pixel file %s\n", pixname);
	    perror ("IRAFWIMAGE");
	    return (0);
	    }
	}

    /* Skip over pixel file header */
    ioff = fseek (fd, (irafheader[IM_PIXOFF]-1) * sizeof(short), SEEK_SET);
    if (ioff) {
	fclose (fd);
	return (0);
	}

    /* Write IRAF pixel file to disk */
    nbw = fwrite (image, 1, nbimage, fd);
    fclose (fd);

    return (nbw);
}

/* Convert FITS image header to IRAF image header, returning IRAF header */

int
fits2iraf (fitsheader, irafheader, nbhead)

char	*fitsheader;	/* FITS image header */
int	*irafheader;	/* IRAF image header (returned updated) */
int	nbhead;		/* Length of IRAF header */

{
    int i;
    char *fitsend, *fitsp, pixfile[64];
    short *irafp, *irafs;
    int	nax,nax4,nbiraf, nlfits, lpixfile;

    /* Delete FITS header keywords not needed by IRAF */
    hgeti4 (fitsheader,"NAXIS",&nax);
    hdel (fitsheader,"SIMPLE");
    hdel (fitsheader,"BITPIX");
    hdel (fitsheader,"NAXIS");
    hgeti4 (fitsheader,"NAXIS1",&irafheader[IM_PHYSLEN]);
    hgeti4 (fitsheader,"NAXIS1",&irafheader[IM_LEN]);
    hdel (fitsheader,"NAXIS1");
    hgeti4 (fitsheader,"NAXIS2",&irafheader[IM_PHYSLEN+1]);
    hgeti4 (fitsheader,"NAXIS2",&irafheader[IM_LEN+1]);
    hdel (fitsheader,"NAXIS2");
    if (nax > 2) {
	hgeti4 (fitsheader,"NAXIS3",&irafheader[IM_PHYSLEN+2]);
	hgeti4 (fitsheader,"NAXIS3",&irafheader[IM_LEN+2]);
	hdel (fitsheader,"NAXIS3");
	}
    if (nax > 3) {
	hgeti4 (fitsheader,"NAXIS4",&irafheader[IM_PHYSLEN+3]);
	hgeti4 (fitsheader,"NAXIS4",&irafheader[IM_LEN+3]);
	hgeti4 (fitsheader,"NAXIS4",&nax4);
	hdel (fitsheader,"NAXIS4");
	}
    hdel (fitsheader,"OBJECT");

    /* Find length of FITS header */
    fitsend = ksearch (fitsheader,"END");
    nlfits = (fitsend - fitsheader) / 80;
    nbiraf = (4 * LEN_IMHDR) + (162 * nlfits);
    if (nbiraf > nbhead)
	irafheader = (int *) realloc (irafheader, nbiraf);

    /* Replace pixel file name, if it is in FITS header */
    if (hgets (fitsheader, "PIXFILE", pixfile)) {
	lpixfile = strlen (pixfile);
	for (i = 0; i < SZ_IMPIXFILE; i++) {
	    if (i < lpixfile)
		irafheader[IM_PIXFILE+i] = (short) pixfile[i];
	    else
		irafheader[IM_PIXFILE+i] = (short) 0;
	    }
	hdel (fitsheader,"PIXFILE");
	}

    /* Remove heder file name from header, if it is there */
    if (ksearch (fitsheader, "IMHFILE"))
	hdel (fitsheader, "IMHFILE");

    /*  Replace user portion of IRAF header with remaining FITS header */
    irafs = (short *)irafheader;
    irafp = irafs + (2 * LEN_IMHDR);
    for (fitsp = fitsheader; fitsp < fitsend; fitsp = fitsp + 80) {
	for (i = 0; i < 80; i++)
	    *irafp++ = (short) fitsp[i];
	*irafp++ = 10;
	}

    /* Return number of bytes in new IRAF header */
    return (2 * (irafp - irafs));
}


/* IRAFSWAP -- Reverse bytes of any type of vector in place */

static void
irafswap (bitpix, string, nbytes)

int	bitpix;		/* Number of bits per pixel */
			/*  16 = short, -16 = unsigned short, 32 = int */
			/* -32 = float, -64 = double */
char	*string;	/* Address of starting point of bytes to swap */
int	nbytes;		/* Number of bytes to swap */

{
    switch (bitpix) {

	case 16:
	    if (nbytes < 2) return;
	    irafswap2 (string,nbytes);
	    break;

	case 32:
	    if (nbytes < 4) return;
	    irafswap4 (string,nbytes);
	    break;

	case -16:
	    if (nbytes < 2) return;
	    irafswap2 (string,nbytes);
	    break;

	case -32:
	    if (nbytes < 4) return;
	    irafswap4 (string,nbytes);
	    break;

	case -64:
	    if (nbytes < 8) return;
	    irafswap8 (string,nbytes);
	    break;

	}
    return;
}


/* IRAFSWAP2 -- Swap bytes in string in place */

static void
irafswap2 (string,nbytes)


char *string;	/* Address of starting point of bytes to swap */
int nbytes;	/* Number of bytes to swap */

{
    char *sbyte, temp, *slast;

    slast = string + nbytes;
    sbyte = string;
    while (sbyte < slast) {
	temp = sbyte[0];
	sbyte[0] = sbyte[1];
	sbyte[1] = temp;
	sbyte= sbyte + 2;
	}
    return;
}


/* IRAFSWAP4 -- Reverse bytes of Integer*4 or Real*4 vector in place */

static void
irafswap4 (string,nbytes)

char *string;	/* Address of Integer*4 or Real*4 vector */
int nbytes;	/* Number of bytes to reverse */

{
    char *sbyte, *slast;
    char temp0, temp1, temp2, temp3;

    slast = string + nbytes;
    sbyte = string;
    while (sbyte < slast) {
	temp3 = sbyte[0];
	temp2 = sbyte[1];
	temp1 = sbyte[2];
	temp0 = sbyte[3];
	sbyte[0] = temp0;
	sbyte[1] = temp1;
	sbyte[2] = temp2;
	sbyte[3] = temp3;
	sbyte = sbyte + 4;
	}

    return;
}


/* IRAFSWAP8 -- Reverse bytes of Real*8 vector in place */

static void
irafswap8 (string,nbytes)

char *string;	/* Address of Real*8 vector */
int nbytes;	/* Number of bytes to reverse */

{
    char *sbyte, *slast;
    char temp[8];

    slast = string + nbytes;
    sbyte = string;
    while (sbyte < slast) {
	temp[7] = sbyte[0];
	temp[6] = sbyte[1];
	temp[5] = sbyte[2];
	temp[4] = sbyte[3];
	temp[3] = sbyte[4];
	temp[2] = sbyte[5];
	temp[1] = sbyte[6];
	temp[0] = sbyte[7];
	sbyte[0] = temp[0];
	sbyte[1] = temp[1];
	sbyte[2] = temp[2];
	sbyte[3] = temp[3];
	sbyte[4] = temp[4];
	sbyte[5] = temp[5];
	sbyte[6] = temp[6];
	sbyte[7] = temp[7];
	sbyte = sbyte + 8;
	}
    return;
}

/*
 * Feb 15 1996	New file
 * Apr 10 1996	Add more documentation
 * Apr 17 1996	Print error message on open failure
 * Jun  5 1996	Add byte swapping (reversal); use streams
 * Jun 10 1996	Make fixes after running lint
 * Jun 12 1996	Use IMSWAP subroutines instead of local ones
 * Jul  3 1996	Go back to using local IRAFSWAP subroutines
 * Jul  3 1996	Write to pixel file from FITS header
 * Jul 10 1996	Allocate all headers
 * Aug 13 1996	Add unistd.h to include list
 */
