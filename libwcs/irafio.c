/* File irafio.c
 * February 15, 1996
 * By Doug Mink, based on Mike VanHilst's readiraf.c
 */

#include <stdio.h>		/* define stderr, FD, and NULL */
#include <fcntl.h>
#include <string.h>
#include "fitshead.h"


/* Parameters from iraf/lib/imhdr.h */
#define SZ_IMPIXFILE	79
#define LEN_IMHDR	513

/* Offsets into header in sizeof(int) units for various parameters */
#define	IM_PIXTYPE	4		/* datatype of the pixels */
#define	IM_NDIM		5		/* number of dimensions */
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

#define LEN_IRAFHDR	15000
#define LEN_PIXHDR	1024
#define LEN_FITSHDR	11520

int check_immagic ();

/*
 * Subroutine:	init_irafimh
 * Purpose:	Open and read the iraf .imh file and set up the pixfile
 		for reading
 * Returns:	-1 if failure, else FD to image file ready for reading data
 * Notes:	The imhdr format is defined in iraf/lib/imhdr.h, some of
 *		which defines or mimiced, above.
 */
int *
irafrhead (filename, lhead, fitsheader)

char	*filename;	/* Name of IRAF header file */
int	lhead;		/* Maximum length of FITS image header in bytes */
char	*fitsheader;	/* FITS image header (filled) */
{
    int fd, nbr;
    int *irafheader;
    int nbfhead, ibpx;

    /* open the image header file */
    if ((fd = open (filename, O_RDONLY, 0)) <= 0 )
	return (NULL);

    /* allocate initial sized buffer */
    irafheader = (int *)malloc (LEN_IRAFHDR, sizeof(short));

    /* read in IRAF header */
    nbr = read (fd, irafheader, LEN_IRAFHDR, 0);
    close (fd);
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

    /* Conver IRAF header to FITS header */
    nbfhead = iraf2fits (irafheader, LEN_IRAFHDR, lhead, fitsheader);

    return (irafheader);
}


char *
irafrimage (hdrname, irafheader, fitsheader)

char	*hdrname;	/* Name of IRAF image header file */
int	*irafheader;	/* IRAF image header */
char	*fitsheader;	/* FITS image header (filled) */
{
    int fd;
    char *bang;
    int i, naxis1, naxis2, bitpix, bytepix;
    void same_path();
    char *pixname;
    char *image;
    int nbr, nbimage;
    int *pixheader;
    short *irafpix;

    /* Convert pixel file name to character string */
    irafpix = (short *) (irafheader + IM_PIXFILE);
    pixname = (char *) malloc (SZ_IMPIXFILE);
    for (i = 0; i < SZ_IMPIXFILE; i++)
	pixname[i] = (char) irafpix[i];

    if (strncmp(pixname, "HDR$", 4) == 0 )
	same_path (pixname, hdrname);

    if ((bang = strchr (pixname, '!')) != NULL ) {
	if ((fd = open (bang + 1, O_RDONLY, 0)) <= 1 ) {
	    (void)fprintf(stderr, "Cannot open IRAF pixel file %s\n", pixname);
	    return (NULL);
	    }
	}
    else {
	if ((fd = open (pixname, O_RDONLY, 0)) <= 1 ) {
	    (void)fprintf(stderr, "Cannot open IRAF pixel file %s\n", pixname);
	    return (NULL);
	    }
	}

    /* Read pixel header */
    pixheader = (int *)malloc (LEN_PIXHDR);
    nbr = read (fd, pixheader, LEN_PIXHDR, 0);

    /* Check size of pixel header */
    if (nbr < LEN_PIXHDR) {
	(void)fprintf(stderr, "IRAF pixel file %s: %d / %d bytes read.\n",
		      pixname,nbr,LEN_PIXHDR);
	free (pixheader);
	return (NULL);
	}

    /* check pixel header magic word */
    if (check_immagic ((short *)pixheader, "impix") ) {
	(void)fprintf(stderr, "File %s not valid IRAF pixel file.\n", pixname);
	free (pixheader);
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
    nbr = read (fd, image, nbimage, 0);
    close (fd);

    /* Check size of image */
    if (nbr < nbimage) {
	(void)fprintf(stderr, "IRAF pixel file %s: %d / %d bytes read.\n",
		      pixname,nbr,nbimage);
	free (image);
	return (NULL);
	}
    return (image);
}


/* Verify that file is valid IRAF imhdr or impix by checking first 5 chars
 * Returns:	0 on success, 1 on failure */

int
check_immagic (irafheader, teststring )

short	*irafheader;	/* IRAF image header from file */
char	*teststring;	/* "imhdr" for header file or "impix" for pixel file */

{
    char line[8];
    int i;

    for (i = 0; i < 5; i++)
	line[i] = (char) irafheader[i];
    if (strncmp(line, teststring, 5) != 0 )
	return(1);
    else
	return(0);
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

/* Convert IRAF image heaader to FITS image header */

int
iraf2fits (irafheader, nbiraf, nbfits, fitsheader)

int	*irafheader;	/* IRAF image header */
int	nbiraf;		/* Number of bytes in IRAF image */
int	nbfits;		/* Maximum length of FITS header */
char	*fitsheader;	/* FITS image header (returned) */

{
    int lfhead;		/* Actual length of FITS header (returned) */
    char objname[64];	/* object name from FITS file */
    int i, j, nax, nbits;
    short *irafobj;

    int ncr, nblock;
    char *fhead, *fhead1, *fp, *fhmax, endline[81];
    short *irafline;
    char fitsline[81];

    (void)strncpy (endline,"END                                     ",40);
    (void)strncpy (endline+40,"                                        ",40);

    /*  Initialize FITS header */
    fhead = fitsheader;
    lfhead = 0;
    fhmax = fhead + nbfits;
    if (fhead+400 > fhmax)
	return 0;
    (void)strncpy (fitsheader, endline, 80);
    hputl (fitsheader, "SIMPLE", 1);
    fhead = fhead + 80;
    lfhead = lfhead + 80;

    /*  Set pixel size in FITS header */
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
	return (0);
    }
    hputi4 (fitsheader,"BITPIX",nbits);
    fhead = fhead + 80;
    lfhead = lfhead + 80;

    /*  Set image dimensions in FITS header */
    nax = irafheader[11];
    hputi4 (fitsheader,"NAXIS",nax);
    fhead = fhead + 80;
    lfhead = lfhead + 80;
    hputi4 (fitsheader,"NAXIS1",irafheader[IM_PHYSLEN]);
    fhead = fhead + 80;
    lfhead = lfhead + 80;
    hputi4 (fitsheader,"NAXIS2",irafheader[IM_PHYSLEN+1]);
    fhead = fhead + 80;
    if (fhead > fhmax) {
	(void)strncpy (fhead-80, endline,80);
	return (lfhead);
	}
    lfhead = lfhead + 80;
    if (nax > 2) {
	hputi4 (fitsheader,"NAXIS3",irafheader[IM_PHYSLEN+2]);
	fhead = fhead + 80;
	if (fhead > fhmax) {
	    (void)strncpy (fhead-80, endline,80);
	    return (lfhead);
	    }
	lfhead = lfhead + 80;
	}
    if (nax > 3) {
	hputi4 (fitsheader,"NAXIS4",irafheader[IM_PHYSLEN+3]);
	fhead = fhead + 80;
	if (fhead > fhmax) {
	    (void)strncpy (fhead-80, endline,80);
	    return (lfhead);
	    }
	lfhead = lfhead + 80;
	}

    /*  Set object name in FITS header */
    irafobj = (short *) (irafheader + IM_TITLE);
    i = 0;
    for (i = 0; i < 64; i++) {
	objname[i] = (char) irafobj[i];
	if (irafobj[i] == 0)
	    break;
	}
    hputs (fitsheader,"OBJECT",objname);
    fhead = fhead + 80;
    if (fhead > fhmax) {
	(void)strncpy (fhead-80, endline,80);
	return (lfhead);
	}
    if (fhead > fhmax) {
	(void)strncpy (fhead-80, endline,80);
	return (lfhead);
	}
    lfhead = lfhead + 80;

    /*  Add user portion of IRAF header to FITS header */
    ncr = nbiraf / 2;
    fitsline[81] = 0;
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
	if (fhead > fhmax) {
	    (void)strncpy (fhead-80, endline,80);
	    return (lfhead);
	    }
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
    fhead1 = fitsheader + (nblock * 2880);

    /* Pad rest of header with spaces */
    strncpy (endline,"   ",3);
    for (fp = fhead+80; fp < fhead1; fp = fp + 80) {
	(void)strncpy (fp, endline,80);
	lfhead = lfhead + 80;
	}

    return (lfhead);
}


int
irafwhead (hdrname, irafheader, fitsheader)

char	*hdrname;	/* Name of IRAF header file */
int	*irafheader;	/* IRAF header */
char	*fitsheader;	/* FITS image header */

{
    int fd;
    char *bang;
    int i, nbw, nbhead;
    short *irafs;

    /* Convert FITS header to IRAF header */
    nbhead = fits2iraf (fitsheader, irafheader);

    /* Write IRAF header to disk file */
    if ((fd = open (hdrname, O_RDWR, 0)) <= 1 ) {
	(void)fprintf(stderr, "Cannot open IRAF header file %s\n", hdrname);
	return (0);
	}
    irafs = (short *)irafheader;
    nbw = write (fd, irafheader, nbhead);
    close (fd);
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
    int fd;
    char *bang;
    int i, nbw, bytepix, bitpix, naxis1, naxis2, nbhead, nbimage;
    void same_path();
    char *pixname;
    short *irafpix;

    /* Convert FITS header to IRAF header */
    nbhead = fits2iraf (fitsheader, irafheader);

    /* Write IRAF header to disk file */
    if ((fd = open (hdrname, O_RDWR, 0)) <= 1 ) {
	(void)fprintf(stderr, "Cannot open IRAF header file %s\n", hdrname);
	return (0);
	}
    nbw = write (fd, irafheader, nbhead);
    close (fd);

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
	if ((fd = open (bang + 1, O_RDWR, 0)) <= 1 ) {
	    (void)fprintf(stderr, "Cannot open IRAF pixel file %s\n", pixname);
	    return (0);
	    }
	}
    else {
	if ((fd = open (pixname, O_RDWR, 0)) <= 1 ) {
	    (void)fprintf(stderr, "Cannot open IRAF pixel file %s\n", pixname);
	    return (0);
	    }
	}

    /* Skip over pixel file header */
    if (lseek (fd, (irafheader[IM_PIXOFF]-1) * sizeof(short), 0) < 0 ) {
	close (fd);
	return (0);
	}

    /* Write IRAF pixel file to disk */
    nbw = write (fd, irafheader, nbimage);
    close (fd);

    return (nbw);
}

/* Convert FITS image heaader to IRAF image header */

int fits2iraf (fitsheader, irafheader)

char	*fitsheader;	/* FITS image header */
int	*irafheader;	/* IRAF image header (returned updated) */

{
    int i;
    char *fitsend, *fitsp;
    short *irafp, *irafs;
    int	naxis,nbiraf, nlfits;

    /* Delete FITS header keywords not needed by IRAF */
    hgeti4 (fitsheader,"NAXIS",&naxis);
    hdel (fitsheader,"SIMPLE");
    hdel (fitsheader,"BITPIX");
    hdel (fitsheader,"NAXIS");
    hdel (fitsheader,"NAXIS1");
    hdel (fitsheader,"NAXIS2");
    if (naxis > 2)
	hdel (fitsheader,"NAXIS3");
    if (naxis > 3)
	hdel (fitsheader,"NAXIS4");
    hdel (fitsheader,"OBJECT");

    /* Find length of FITS header */
    fitsend = ksearch (fitsheader,"END");
    nlfits = (fitsend - fitsheader) / 80;
    nbiraf = (4 * LEN_IMHDR) + (162 * nlfits);
    if (nbiraf > LEN_IRAFHDR)
	irafheader = (int *) realloc (irafheader, nbiraf);

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
