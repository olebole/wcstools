/* File imstack.c
 * April 9, 2002
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
#include "libwcs/fitsfile.h"
#include "libwcs/wcs.h"

static void usage();
static int StackImage();

static int verbose = 0;		/* verbose flag */
static int wfits = 1;		/* if 1, write FITS header before data */
static char *newname = NULL;
static int nfiles = 0;
static int nbstack = 0;
static FILE *fstack = NULL;
static int version = 0;		/* If 1, print only program name and version */
static char *outfile = NULL;	/* If not null, output filename */

main (ac, av)
int ac;
char **av;
{
    char filename[128];
    char *filelist[100];
    char *listfile;
    char *str;
    int readlist = 0;
    int ntimes = 1;
    FILE *flist;
    int ifile, nblocks, nbytes, i;
    char *blanks;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

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

	case 'i':	/* Write only image data to output file */
	    wfits = 0;
	    if (newname == NULL) {
		newname = calloc (16, 1);
		strcpy (newname,"imstack.out");
		}
	    break;

	case 'n':	/* Use each input file this many times */
	    if (ac < 2)
                usage ();
	    ntimes = (int) atof (*++av);
	    ac--;
	    break;

	case 'o':	/* Set output file name */
	    if (ac < 2)
                usage ();
	    if (newname != NULL)
		free (newname);
	    newname = *++av;
	    ac--;
	    break;

	case '@':	/* List of files to be read */
	    readlist++;
	    listfile = ++str;
	    str = str + strlen (str) - 1;
	    av++;
	    ac--;
	    break;

	default:
	    usage();
	    break;
	}
    }

    /* If output file name has not yet been set, set it */
    if (newname == NULL) {
	newname = calloc (16, 1);
	strcpy (newname,"imstack.fits");
	}

    /* Find number of images to stack  and leave listfile open for reading */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMSTACK: List file %s cannot be read\n",
		     listfile);
	    usage();
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
	usage();

    /* Read ac remaining file names starting at av[0] */
    else {
	while (ac-- > 0) {
	    filelist[nfiles] = *av++;
	    nfiles++;
	    if (verbose)
		printf ("Reading %s\n",filelist[nfiles-1]);
	    }
	}

    /* Stack images from list in file */
    if (readlist) {
	for (ifile = 0;  ifile < nfiles; ifile++) {
	    if (fgets (filename, 128, flist) != NULL) {
		filename[strlen (filename) - 1] = 0;
		if (StackImage (ifile, ntimes, filename))
		    break;
		}
	    }
	fclose (flist);
	}

    /* Stack images from list on command line */
    else {
	for (ifile = 0;  ifile < nfiles; ifile++) {
	    if (StackImage (ifile, ntimes, filelist[ifile]))
		break;
	    }
	}


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
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Stack FITS or IRAF images into single FITS image\n");
    fprintf(stderr,"Usage: imstack [-vi][-o filename][-n num] file1.fits file2.fits ... filen.fits\n");
    fprintf(stderr,"  or : imstack [-vi][-n num] @filelist\n");
    fprintf(stderr,"  -i: Do not put FITS header in output file\n");
    fprintf(stderr,"  -n: Use each file this many times\n");
    fprintf(stderr,"  -o: Output filename\n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}


static int
StackImage (ifile, ntimes, filename)

int	ifile;		/* Sequence number of input file */
int	ntimes;		/* Stack each image this many times */
char	*filename;	/* FITS or IRAF file filename */

{
    char *image;		/* FITS image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    char *irafheader;		/* IRAF image header */
    int nbimage, naxis, naxis1, naxis2, naxis3, bytepix;
    int bitpix, nblocks, nbytes;
    int iraffile;
    int i, itime, nout;
    char pixname[256];

    /* Open IRAF header */
    if (isiraf (filename)) {
	iraffile = 1;
	if ((irafheader = irafrhead (filename, &lhead)) != NULL) {
	    nbhead = 0;
	    if ((header = iraf2fits (filename,irafheader,lhead, &nbhead)) == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",filename);
		free (irafheader);
		return (1);
		}
	    if ((image = irafrimage (header)) == NULL) {
		hgetm (header,"PIXFIL", 255, pixname);
		fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
		free (irafheader);
		free (header);
		return (1);
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", filename);
	    return (1);
	    }
	}

    /* Read FITS image header */
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
    naxis = 0;
    naxis1 = 1;
    hgeti4 (header,"NAXIS1",&naxis1);
    if (naxis1 > 1)
	naxis = naxis + 1;
    naxis2 = 1;
    hgeti4 (header,"NAXIS2",&naxis2);
    if (naxis2 > 1)
	naxis = naxis + 1;
    naxis3 = 1;
    hgeti4 (header,"NAXIS3",&naxis3);
    if (naxis3 > 1)
	naxis = naxis + 1;
    hgeti4 (header,"BITPIX",&bitpix);
    bytepix = bitpix / 8;
    if (bytepix < 0) bytepix = -bytepix;
    nbimage = naxis1 * naxis2 * naxis3 * bytepix;
    nbstack = nbstack + nbimage;

    /* Set NAXIS2 to # of images stacked; pad out FITS header to 2880 blocks */
    if (ifile < 1 && wfits) {
	if (naxis == 1) {
	    hputi4 (header,"NAXIS", 2);
	    hputi4 (header,"NAXIS2", nfiles*ntimes);
	    }
	else if (naxis == 2) {
	    hputi4 (header,"NAXIS", 3);
	    hputi4 (header,"NAXIS3", nfiles*ntimes);
	    }
	else {
	    hputi4 (header,"NAXIS", 4);
	    hputi4 (header,"NAXIS4", nfiles*ntimes);
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
    for (itime = 0; itime < ntimes; itime++) {
	nout = (ifile * ntimes) + itime + 1;
	if (fwrite (image, (size_t) 1, (size_t) nbimage, fstack)) {
	    if (verbose) {
		if (iraffile)
		    printf ("IRAF %d bytes of file %s added to %s[%d]",
			    nbimage, filename, newname, nout);
		else
		    printf ("FITS %d bytes of file %s added to %s[%d]",
			    nbimage, filename, newname, nout);
		}
	    if (itime == 0 || itime == ntimes-1)
		printf ("\n");
	    else
		printf ("\r");
	    }
	else {
	    if (iraffile)
		printf ("IRAF file %s NOT added to %s[%d]\n",
		        filename, newname, nout);
	    else
		printf ("FITS file %s NOT added to %s[%d]\n",
		        filename, newname, nout);
	    }
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
 * Oct 14 1998	Use isiraf() to determine file type
 * Nov 30 1998	Add version and help commands for consistency
 *
 * Oct 22 1999	Drop unused variables after lint
 *
 * Feb  7 2000	Add option to repeat files in stack
 * Mar 23 2000	Use hgetm() to get the IRAF pixel file name, not hgets()
 * Sep  6 2000	Add -o option to set output filename
 * Sep  8 2000	Default ntimes to 1 so program works
 *
 * Apr  9 2002	Do not free unallocated header
 */
