/* File remap.c
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
static char *newname = "remap.fit";

static int nfiles = 0;
static int nbstack = 0;
static FILE *fstack = NULL;
static int RemapImage();
extern struct WorldCoor *GetFITSWCS();
static struct WorldCoor *wcsout = NULL;
static char outheader[14400];

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    char rastr[16];
    char decstr[16];
    char filename[128];
    char *filelist[100];
    char *listfile;
    int readlist = 0;
    FILE *flist;
    int ifile, nblocks, nbytes, i, nbw, nx, ny;
    char *blanks;

    setrot (0.0);

    /* Crack arguments */
    for (av++; --ac > 0 && (*(str = *av)=='-' || *str == '@'); av++) {
	char c;
	if (*str == '@')
	    str = str - 1;
	while (c = *++str)
	switch (c) {
    	case 'a':	/* Output rotation angle in degrees */
    	    if (ac < 2)
    		usage();
    	    setrot (atof (*++av));
    	    ac--;
    	    break;

    	case 'b':	/* Output image center on command line in B1950 */
    	    if (ac < 3)
    		usage();
	    setsys (WCS_B1950);
	    strcpy (rastr, *++av);
	    ac--;
	    strcpy (decstr, *++av);
	    ac--;
	    setcenter (rastr, decstr);
    	    break;

	case 'e':	/* Set WCS projection
	    if (ac < 2)
		usage();
	    setwcsproj (*++av);
	    ac--;
	    break; */
	    
    	case 'j':	/* Output image center on command line in J2000 */
    	    if (ac < 3)
    		usage();
	    setsys (WCS_J2000);
	    strcpy (rastr, *++av);
	    ac--;
	    strcpy (decstr, *++av);
	    ac--;
	    setcenter (rastr, decstr);
    	    break;

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

    	case 'o':	/* Specifiy output image filename */
    	    if (ac < 2)
    		usage();
	    strcpy (outname, *++av);
	    ac--;
	    if (outname[0] == '-')
		overwrite++;
	    else
		overwrite = 0;
    	    writeheader++;
    	    break;

    	case 'p':	/* Output image plate scale in arcseconds per pixel */
    	    if (ac < 2)
    		usage();
    	    setsecpix (atof (*++av));
    	    ac--;
    	    break;

	case 'x':	/* X and Y coordinates of output image reference pixel */
	    if (ac < 3)
		usage();
	    x = atof (*++av);
	    ac--;
	    y = atof (*++av);
	    ac--;
    	    setrefpix (x, y);
    	    break;

	case 'y':	/* Dimensions of output image in pixels */
	    if (ac < 3)
		usage();
	    nx = atoi (*++av);
	    ac--;
	    ny = atoi (*++av);
	    ac--;
    	    setnpix (nx, ny);
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

    /* Set up FITS header and WCS for output file */

    /* Remap images */
    for (ifile = 0;  ifile < nfiles; ifile++) {
	if (readlist) {
	    if (fgets (filename, 128, flist) != NULL) {
		filename[strlen (filename) - 1] = 0;
		if (RemapImage (wcs, ifile, filename))
		    break;
		}
	    }
	else {
	    if (RemapImage (ifile, filelist[ifile]))
		break;
	    }
	}

    if (readlist)
	fclose (flist);

    /* Pad out FITS file to 2880 blocks */
    if (fstack != NULL) {
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
    fprintf (stderr,"Remap FITS or IRAF images into single FITS image using WCS\n");
    fprintf(stderr,"%s: usage: [-vf] file1.fit file2.fit ... filen.fit\n", progname);
    fprintf(stderr,"%s: usage: [-vf] @filelist\n", progname);
    fprintf(stderr,"Usage: [-vodfl] [-m mag] [-n frac] [-s mode] [-g class] [-h maxref] [-i peak]\n");
    fprintf(stderr,"       [-c catalog] [-p scale] [-b ra dec] [-j ra dec] [-r deg] [-t tol] [-x x y] [-y frac]\n");
    fprintf(stderr,"       FITS or IRAF file(s)\n");
    fprintf(stderr,"  -a: output rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b: output center in B1950 (FK4) RA and Dec\n");
    fprintf(stderr,"  -e: WCS type (TAN default)\n");
    fprintf(stderr,"  -j: output center in J2000 (FK5) RA and Dec\n");
    fprintf(stderr,"  -o: name for output image, - to overwrite\n");
    fprintf(stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -x: X and Y coordinates of reference pixel (default is center)\n");
    fprintf(stderr,"  -y: output image dimensions (default is first input image)\n");
    exit (1);
}


static int
RemapImage (ifile, filename)

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
    int i, hpout, wpout;
    char pixname[128];
    struct WorldCoor *wcsin;

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

    /* Set output world coordinate system from command line and first image header */
    if (wcsout == NULL) {
	wcsout = GetFITSWCS (header, verbose, &cra, &cdec, &dra, &ddec, &secpix,
			     &wpout, &hpout, equinox);

    /* Set output header from command line and first image header */
	strcpy (headout, header);
	hputi4 (headout,"NAXIS", 2);
	hputi4 (headout,"NAXIS1", hpout);
	hputi4 (headout,"NAXIS2", wpout);
	SetFITSWCS (headout, wcsout);

    /* Set input world coordinate system for this image */
    wcsin = wcsinit (header);

    hgeti4 (header,"BITPIX",&bitpix);
    bytepix = bitpix / 8;
    if (bytepix < 0) bytepix = -bytepix;
    nbimage = naxis1 * naxis2 * naxis3 * bytepix;
    nbstack = nbstack + nbimage;

    /* Set NAXIS2 to # of images stacked; pad out FITS header to 2880 blocks */
    if (ifile < 1) {
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
	if (fwrite (header, (size_t) 1, (size_t) nbhead, fstack)) {
	    if (verbose)
		printf ("%d-byte FITS header from %s written to %s\n",
			nbhead, filename, newname);
	    }
	else
	    printf ("FITS file %s cannot be written.\n", newname);
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
 * Apr 24 1998	change coordinate setting to setsys() from setfk4()
 * Apr 28 1998	Change coordinate system flags to WCS_*
 * May 28 1998	Include fitsio.h instead of fitshead.h
 * Jul 24 1998	Make irafheader char instead of int
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 */
