/*** File wcslib/fitsio.c
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 *** February 8, 1996
 */

#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/file.h>
#include <string.h>
#include "fitshead.h"

static int verbose=0;		/* if 1 print diagnostics */

/* FITSRTOPEN -- Open FITS table file and return header and pointers to
 *		 selected keywords, as well as file descriptor
 */

int fitsrtopen (inpath, lhead, header, nk, kw, nrows, nchar, nbhead)

char	*inpath;	/* Pathname for FITS tables file to read */
int	lhead;		/* Maximum length in bytes of FITS header */
char	*header;	/* Header for FITS tables file to read */
int	nk;		/* Number of keywords to use */
struct Keyword	*kw;	/* Structure for desired entries */
int	*nrows;		/* Number of rows in table (returned) */
int	*nchar;		/* Number of characters in one table row (returned) */
int	*nbhead;	/* Number of characters before table starts */

{
    struct Keyword pw[100];	/* Structure for desired entries */
    int lpnam[100];		/* length of name for each field */
    int nfields;
    int nleft,ifield,ipt,lf0,lfl,lfield,ier,ik,ln,lnm;
    char *h0, *h1, *tf1, *tf2;
    char tname[12];
    char temp[16];
    char tform[16];
    int tverb;
    int fd;

/* Allocate space for FITS header */
    h0 = header;

/* Read FITS header from input file */
    fd = fitsropen (inpath, lhead, header, nbhead);
    if (fd <= 0) {
	fprintf (stderr,"FITSRTOPEN:  %s is not a FITS file\n",inpath);
	return (-1);
	}

/* Make sure this file is really a FITS table file */
    hgets (header,"XTENSION",16,temp);
    if (strncmp (temp, "TABLE", 5) != 0) {
	fprintf (stderr, "FITSRTOPEN:  %s is not a FITS table file\n",inpath);
	free (temp);
	return (-1);
	}

/* Get table size from FITS header */
    *nchar = 0;
    hgeti4 (header,"NAXIS1",nchar);
    *nrows = 0;
    hgeti4 (header,"NAXIS2", nrows);
    if (*nrows <= 0 || *nchar <= 0) {
	fprintf (stderr, "FITSRTOPEN: cannot read %d x %d table %s\n",
		 *nrows,*nchar,inpath);
	return (-1);
	}

/* Set up table for access to individual fields */
    hgeti4 (header,"TFIELDS",&nfields);
    if (verbose)
	printf ("FITSRTHEAD: %d fields per table entry\n", nfields);
    nleft = nk;
    tverb = verbose;
    verbose = 0;
    for (ifield = 0; ifield < nfields; ifield++) {
	if (nleft <= 0)
	    break;

    /* First column of field */
	sprintf (tname, "TBCOL%d", ifield+1);
	h1 = ksearch (h0,tname);
	hgeti4 (h0,tname, &pw[ifield].kf);

    /* Length of field */
	sprintf (tname, "TFORM%d", ifield+1);;
	hgets (h0,tname,16,tform);
	tf1 = tform + 1;
	tf2 = strchr (tform,'.');
	if (tf2 != NULL)
	    *tf2 = ' ';
	pw[ifield].kl = atoi (tf1);

    /* Name of field */
	sprintf (tname, "TTYPE%d", ifield+1);;
	hgets (h0,tname,16,temp);
	strcpy (pw[ifield].kname,temp);
	lpnam[ifield] = strlen (pw[ifield].kname);
	h0 = h1 - 1;
	}

/* Set up table for access to desired fields */
    verbose = tverb;
    if (verbose)
	printf ("FITSRTHEAD: keywords from  %s\n", inpath);

/* Find each desired keyword in the header */
    for (ik = 0; ik < nk; ik++) {
	if (kw[ik].kn <= 0) {
	    for (ifield = 0; ifield < nfields; ifield++) {
		ln = lpnam[ifield];
		if (strncmp (pw[ifield].kname, kw[ik].kname, ln) == 0) {
		    break;
		    }
		}
	    }
	else
	    ifield = kw[ik].kn - 1;

/* Set pointer, lentth, and name in returned array of structures */
	kw[ik].kn = ifield + 1;
	kw[ik].kf = pw[ifield].kf - 1;
	kw[ik].kl = pw[ifield].kl;
	strcpy (kw[ik].kname, pw[ifield].kname);
	}

    return (fd);
}


int fitsropen (inpath, lhead, header, nbhead)

char	*inpath;	/* Pathname for FITS tables file to read */
int	lhead;		/* Maximum length in bytes of FITS header */
char	*header;	/* Header for FITS tables file to read */
int	*nbhead;	/* Number of bytes of header information */
			/* This includes the first, simple header block */

{
    int extend;
    int nbytes,ih0,ih1,iend,nbentry,naxis;
    int ntry,nbr,ier,irec,nrec;
    char fitsbuf[2884];
    char *headend;	/* Pointer to last line of header */
    char *headnext;	/* Pointer to next line of header to be added */
    int fd;		/* file descriptor for FITS tables file */

/* Open input file */
    *nbhead = 0;
    fitsbuf[FITSBLOCK] = 0;
    for (ntry = 0; ntry < 3; ntry++) {
	fd = open (inpath, O_RDONLY, 0);
	if (fd > 2)
	    break;
	else if (ntry == 2) {
	    fprintf (stderr, "FITSRHEAD:  cannot read file %s\n", inpath);
	    return (-1);
	    }
	}

    if (verbose)
	printf ("FITSRHEAD:  input file %s opened\n",inpath);

    nbytes = FITSBLOCK;
    ih0 = 1;
    ih1 = FITSBLOCK;
    *nbhead = 0;
    headend = NULL;
    headnext = header;
    nrec = lhead / FITSBLOCK;

/* Read FITS header from input file one FITS block at a time */
    irec = 0;
    while (irec < nrec) {
	nbytes = FITSBLOCK;
	for (ntry = 0; ntry < 10; ntry++) {
	    nbr = read (fd, fitsbuf, nbytes);

/* Short records are allowed only if they contain the last header line */
	    if (nbr < nbytes) {
		headend = ksearch (fitsbuf,"END");
		if (headend == NULL) {
		    if (ntry < 9) {
			if (verbose)
			    printf ("FITSRHEAD: %d / %d bytes read %d\n",
				     nbr,nbytes,ntry);
			}
		    else {
			fprintf(stderr, "FITSRHEAD: '%d / %d bytes of header read from %s\n"
				,nbr,nbytes,inpath);
			close (fd);
			return (-1);
			}
		    }
		else
		    break;
		}
	    else
		break;
	    }

/* Move current FITS record into header string */
	strncpy (headnext, fitsbuf, nbr);
	*nbhead = *nbhead + nbr;
	irec = irec + 1;

/* Check to see if this is the final record in this header */
	headend = ksearch (fitsbuf,"END");
	if (headend == NULL)
	    headnext = headnext + nbr;

/* If this is not the real header, read it, starting with the next record */
	else {
	    naxis = 0;
	    hgeti4 (header,"NAXIS",&naxis);
	    extend = 0;
	    hgetl (header,"EXTEND",&extend);
	    if (naxis == 0 && extend) {
		headnext = header;
		headend = NULL;
		irec = 0;
		}
	    else
		break;
	    }
	}

    return (fd);
}

static int offset1=0;
static int offset2=0;

int fitsrtline (fd, nbhead, lbuff, tbuff, irow, nbline, line)

int	fd;		/* File descriptor for FITS file */
int	nbhead;		/* Number of bytes in FITS header */
int	lbuff;		/* Number of bytes in table buffer */
char	*tbuff;		/* FITS table buffer */
int	irow;		/* Number of table row to read */
int	nbline;		/* Number of bytes to read for this line */
char	*line;		/* One line of FITS table (returned) */

{
    int nbuff,nlbuff,ipos,nbr,nbread;
    int offset, offend, ntry, ioff;
    char *tbuff1;

    offset = nbhead + (nbline * irow);
    offend = offset + nbline - 1;

/* Read a new buffer of the FITS table into memory if needed */
    if (offset < offset1 || offend > offset2) {
	nlbuff = lbuff / nbline;
	nbuff = nlbuff * nbline;
	ipos = 0;
	for (ntry = 0; ntry < 3; ntry++) {
	    ioff = lseek (fd, offset, L_SET);
	    if (ioff < offset) {
		if (ntry == 2)
		    return (0);
		else
		    continue;
		}
	    nbr = read (fd, tbuff, nbuff);
	    if (nbr < nbline) {
		if (verbose)
		    printf ("FITSRHEAD: %d / %d bytes read %d\n",
				nbr,nbuff,ntry);
		if (ntry == 2)
		    return (nbr);
		}
	    else
		break;
	    }
	offset1 = offset;
	offset2 = offset + nbr - 1;
	strncpy (line, tbuff, nbline);
	nbread = nbline;
	return (nbline);
	}
    else {
	tbuff1 = tbuff + (offset - offset1);
	strncpy (line, tbuff1, nbline);
	return (nbline);
	}
}

fitsrtlset ()
{
	offset1 = 0;
	offset2 = 0;
}

int
fitsrhead (filename, lhead, header)

char	*filename;	/* Name of FITS image file */
int	lhead;		/* Maximum length of image header in bytes */
char	*header;	/* FITS image header (filled) */

{
    float scale, bias;
    int fd;
    int nbhead, nbimage, naxis1, naxis2, bytepix, nbr;
    int bitpix, naxis, nblocks, nbytes;
    char *image;

    /* Open the image file and read the header */
    if ((fd = fitsropen (filename, lhead, header, &nbhead)) <= 0 )
	return (NULL);

    close (fd);
    return (nbhead);
}


char *
fitsrimage (filename, lhead, header)

char	*filename;	/* Name of IFTS image file */
int	lhead;		/* Maximum length of image header in bytes */
char	*header;	/* FITS image header (filled) */

{
    float scale, bias;
    int fd;
    int nbhead, nbimage, naxis1, naxis2, bytepix, nbr;
    int bitpix, naxis, nblocks, nbytes;
    char *image;

    /* Open the image file and read the header */
    if ((fd = fitsropen (filename, lhead, header, &nbhead)) <= 0 )
	return (NULL);

    /* Compute size of image in bytes using relevant header parameters */
    hgeti4 (header,"NAXIS",&naxis);
    hgeti4 (header,"NAXIS1",&naxis1);
    hgeti4 (header,"NAXIS2",&naxis2);
    hgeti4 (header,"BITPIX",&bitpix);
    bytepix = bitpix / 8;
    if (bytepix < 0) bytepix = -bytepix;
    nbimage = naxis1 * naxis2 * bytepix;

    /* Set number of bytes to integral number of 2880-byte blocks */
    nblocks = nbimage / FITSBLOCK;
    if (nblocks * FITSBLOCK < nbimage)
	nblocks = nblocks + 1;
    nbytes = nblocks * FITSBLOCK;

    /* Allocate and read image */
    image = malloc (nbytes);
    nbr = read (fd, image, nbytes);
    close (fd);
    if (nbr < nbimage)
	return (NULL);
    return (image);
}

int
fitswimage (filename, header, image)

char	*filename;	/* Name of IFTS image file */
char	*header;	/* FITS image header */
char	*image;		/* FITS image pixels */

{
    float scale, bias;
    int fd;
    int nbhead, nbimage, nblocks, bytepix;
    int bitpix, naxis, naxis1, naxis2, nbytes, nbw;

    /* Open the output file */
    fd = open (filename, O_RDWR, 0);
    if (fd < 3) {
	fprintf (stderr, "FITSWIMAGE:  cannot write file %s\n", filename);
	return (-1);
	}

    /* Write header to file */
    nbhead = ksearch (header,"END") - header + 80;
    nblocks = nbhead / FITSBLOCK;
    if (nblocks * FITSBLOCK < nbhead)
	nblocks = nblocks + 1;
    nbytes = nblocks * FITSBLOCK;
    nbw = write (fd, header, nbytes);
    if (nbw < nbhead) {
	fprintf (stderr, "FITSWIMAGE:  wrote %d / %d bytes of header to file %s\n",
		 nbw, nbytes, filename);
	close (fd);
	return (-1);
	}

    /* Compute size of image in bytes using relevant header parameters */
    hgeti4 (header,"NAXIS",&naxis);
    hgeti4 (header,"NAXIS1",&naxis1);
    hgeti4 (header,"NAXIS2",&naxis2);
    hgeti4 (header,"BITPIX",&bitpix);
    bytepix = bitpix / 8;
    if (bytepix < 0) bytepix = -bytepix;
    nbimage = naxis1 * naxis2 * bytepix;
    nblocks = nbimage / FITSBLOCK;
    if (nblocks * FITSBLOCK < nbimage)
	nblocks = nblocks + 1;
    nbytes = nblocks * FITSBLOCK;

    /* Write image to file */
    nbw = write (fd, image, nbytes);
    close (fd);
    if (nbw < nbimage) {
	fprintf (stderr, "FITSWIMAGE:  wrote %d / %d bytes of image to file %s\n",
		 nbw, nbytes, filename);
	return (-1);
	}
    return (nbw);
}
