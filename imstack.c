/* File imstack.c
 * August 6, 1998
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include "fitsfile.h"
#include "wcs.h"

static void usage();
static int verbose = 0;		/* verbose flag */
static int debug = 0;		/* debugging flag */
static int wfits = 0;		/* if 1, write FITS header before data */
static char *newname = "imstack.out";

static int nfiles = 0;
static int nbstack = 0;
static FILE *fstack = NULL;
static int StackImage();

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char filename[128];
    char *filelist[100];
    char *listfile;
    char *str;
    int readlist = 0;
    FILE *flist;
    int ifile, nblocks, nbytes, i, nbw;
    char *blanks;

    /* Crack arguments */
    for (av++; --ac > 0 && (*(str = *av)=='-' || *str == '@'); av++) {
	char c;
	if (*str == '@')
	    str = str - 1;
	while (c = *++str)
	switch (c) {
	case 'v':	/* more verbosity */
	    verbose++;
	    break;

	case 'f':	/* Write FITS header to output file */
	    wfits++;
	    strcpy (newname,"imstack.fit");
	    break;

	case 'k':	/* Print extra debugging information */
	    debug++;
	    break;

	case '@':	/* List of files to be read */
	    readlist++;
	    listfile = ++str;
	    str = str + strlen (str) - 1;
	    av++;
	    ac--;
	    break;

	default:
	    usage (progname);
	    break;
	}
    }

    /* Find number of images to stack  and leave listfile open for reading */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMSTACK: List file %s cannot be read\n",
		     listfile);
	    usage (progname);
	    }
	while (fgets (filename, 128, flist) != NULL)
	    nfiles++;
	fclose (flist);
	flist = NULL;
	if (nfiles > 0)
	    flist = fopen (listfile, "r");
	}

    /* If no arguments left, print usage */
    else if (ac == 0)
	usage (progname);

    /* Read ac remaining file names starting at av[0] */
    else {
	while (ac-- > 0) {
	    filelist[nfiles] = *av++;
	    nfiles++;
	    if (verbose)
		printf ("Reading %s\n",filelist[nfiles--]);
	    }
	}

    /* Stack images */
    for (ifile = 0;  ifile < nfiles; ifile++) {
	if (readlist) {
	    if (fgets (filename, 128, flist) != NULL) {
		filename[strlen (filename) - 1] = 0;
		if (StackImage (ifile, filename))
		    break;
		}
	    }
	else {
	    if (StackImage (ifile, filelist[ifile]))
		break;
	    }
	}

    if (readlist)
	fclose (flist);

    /* Pad out FITS file to 2880 blocks */
    if (wfits && fstack != NULL) {
	nblocks = nbstack / FITSBLOCK;
	if (nblocks * FITSBLOCK < nbstack)
	    nblocks = nblocks + 1;
	nbytes = (nblocks * FITSBLOCK) - nbstack;
	if (nbytes > 0) {
	    blanks = (char *) malloc ((size_t) nbytes);
	    for (i = 0;  i < nbytes; i++)
		blanks[i] = 0;
	    (void) fwrite (blanks, (size_t) 1, (size_t)nbytes, fstack);
	    free (blanks);
	    }
	}

    if (fstack != NULL)
	fclose (fstack);
    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"Stack FITS or IRAF images into single FITS image\n");
    fprintf(stderr,"%s: usage: [-vf] file1.fit file2.fit ... filen.fit\n", progname);
    fprintf(stderr,"%s: usage: [-vf] @filelist\n", progname);
    fprintf(stderr,"  -f: Print output FITS header, else do not \n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}


static int
StackImage (ifile, filename)

int	ifile;		/* Sequence number of input file */
char	*filename;	/* FITS or IRAF file filename */

{
    char *image;		/* FITS image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    char *irafheader;		/* IRAF image header */
    int nbimage, naxis, naxis1, naxis2, naxis3, naxis4, bytepix;
    int bitpix, nblocks, nbytes;
    int iraffile;
    int i;
    char pixname[128];

    /* Open IRAF header if .imh extension is present */
    if (strsrch (filename,".imh") != NULL) {
	iraffile = 1;
	if ((irafheader = irafrhead (filename, &lhead)) != NULL) {
	    nbhead = 0;
	    if ((header = iraf2fits (filename,irafheader,lhead, &nbhead)) == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",filename);
		free (irafheader);
		return (1);
		}
	    if ((image = irafrimage (header)) == NULL) {
		hgets (header,"PIXFILE", 64, pixname);
		fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
		free (irafheader);
		free (header);
		return (1);
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", filename);
	    free (header);
	    return (1);
	    }
	}

    /* Read FITS image header if .imh extension is not present */
    else {
	iraffile = 0;
	if ((header = fitsrhead (filename, &lhead, &nbhead)) != NULL) {
	    if ((image = fitsrimage (filename, nbhead, header)) == NULL) {
		fprintf (stderr, "Cannot read FITS image %s\n", filename);
		free (header);
		return (1);
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    return (1);
	    }
	}
    if (ifile < 1 && verbose)

    /* Compute size of input image in bytes from header parameters */
    hgeti4 (header,"NAXIS",&naxis);
    hgeti4 (header,"NAXIS1",&naxis1);
    naxis2 = 1;
    hgeti4 (header,"NAXIS2",&naxis2);
    naxis3 = 1;
    hgeti4 (header,"NAXIS3",&naxis3);
    hgeti4 (header,"BITPIX",&bitpix);
    bytepix = bitpix / 8;
    if (bytepix < 0) bytepix = -bytepix;
    nbimage = naxis1 * naxis2 * naxis3 * bytepix;
    nbstack = nbstack + nbimage;

    /* Set NAXIS2 to # of images stacked; pad out FITS header to 2880 blocks */
    if (ifile < 1 && wfits) {
	if (naxis == 1) {
	    hputi4 (header,"NAXIS", 2);
	    hputi4 (header,"NAXIS2", nfiles);
	    }
	else if (naxis == 2) {
	    hputi4 (header,"NAXIS", 3);
	    hputi4 (header,"NAXIS3", nfiles);
	    }
	else {
	    hputi4 (header,"NAXIS", 4);
	    hputi4 (header,"NAXIS4", nfiles);
	    }
	nbhead = strlen (header);
	nblocks = nbhead / FITSBLOCK;
	if (nblocks * FITSBLOCK < nbhead)
	    nblocks = nblocks + 1;
	nbytes = nblocks * FITSBLOCK;
	for (i = nbhead+1; i < nbytes; i++)
	    header[i] = ' ';
	nbhead = nbytes;
	}

    /* If first file open to write and, optionally, write FITS header */
    if (ifile == 0) {
	fstack = fopen (newname, "w");
	if (fstack == NULL) {
	    fprintf (stderr, "Cannot write image %s\n", newname);
	    return (1);
	    }
	if (wfits) {
	    if (fwrite (header, (size_t) 1, (size_t) nbhead, fstack)) {
		if (verbose) {
		    printf ("%d-byte FITS header from %s written to %s\n",
			    nbhead, filename, newname);
		    }
		}
	    else
		printf ("FITS file %s cannot be written.\n", newname);
	    }
	}
    else if (fstack == NULL) {
	fstack = fopen (newname, "a");
	if (fstack == NULL) {
	    fprintf (stderr, "Cannot write image %s\n", newname);
	    return (1);
	    }
	}
    if (fwrite (image, (size_t) 1, (size_t) nbimage, fstack)) {
	if (verbose) {
	    if (iraffile)
		printf ("IRAF file %s %d bytes added to %s[%d]\n",
			filename, nbimage, newname, ifile+1);
	    else
		printf ("FITS file %s %d bytes added to %s[%d]\n",
			filename, nbimage, newname, ifile+1);
	    }
	}
    else {
	if (iraffile)
	    printf ("IRAF file %s NOT added to %s[%d]\n",
		    filename, newname, ifile+1);
	else
	    printf ("FITS file %s NOT added to %s[%d]\n",
		    filename, newname, ifile+1);
	}

    if (ifile < 1) {
	fclose (fstack);
	fstack = NULL;
	}

    free (header);
    free (image);
    return (0);
}

/* May 15 1997	New program
 * May 30 1997	Fix FITS data padding to integral multiple of 2880 bytes
 *
 * May 28 1998	Include fitsio.h instead of fitshead.h
 * Jul 24 1998	Make irafheader char instead of int
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 */
