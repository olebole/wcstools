/* File remap.c
 * April 29, 1999
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
static char *outname0 = "remap.fits";
static char *outname;

static int nfiles = 0;
static int nbstack = 0;
static int fitsout = 0;
static FILE *fstack = NULL;
static int RemapImage();
extern struct WorldCoor *GetFITSWCS();
static struct WorldCoor *wcsout = NULL;
static char outheader[14400];
static int eqsys = -1;
static double equinox = 0.0;
static int version = 0;		/* If 1, print only program name and version */

main (ac, av)
int ac;
char **av;
{
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
    double x, y;

    outname = outname0;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

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

	case 'f':	/* Force FITS output */
	    fitsout++;
	    break;
	    
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

	case 'k':	/* Print extra debugging information */
	    debug++;
	    break;

    	case 'o':	/* Specifiy output image filename */
    	    if (ac < 2)
    		usage();
	    outname = *++av;
	    ac--;
    	    break;

    	case 'p':	/* Output image plate scale in arcseconds per pixel */
    	    if (ac < 2)
    		usage();
    	    setsecpix (atof (*++av));
    	    ac--;
    	    break;

	case 'w':	/* Set WCS projection */
	    if (ac < 2)
		usage();
	    setwcsproj (*++av);
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
	    usage();
	    break;
	}
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
		printf ("Reading %s\n",filelist[nfiles--]);
	    }
	}

    /* Set up FITS header and WCS for output file */

    /* Remap images */
    for (ifile = 0;  ifile < nfiles; ifile++) {
	if (readlist) {
	    if (fgets (filename, 128, flist) != NULL) {
		filename[strlen (filename) - 1] = 0;
		if (RemapImage (ifile, filename))
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
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Remap FITS or IRAF images into single FITS image using WCS\n");
    fprintf(stderr,"usage: remap [-vf] file1.fit file2.fit ... filen.fit\n");
    fprintf(stderr,"       remap [-vf] @filelist\n");
    fprintf(stderr,"  -a: output rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b: output center in B1950 (FK4) RA and Dec\n");
    fprintf(stderr,"  -j: output center in J2000 (FK5) RA and Dec\n");
    fprintf(stderr,"  -o: name for output image\n");
    fprintf(stderr,"  -p: initial plate scale in arcsec per pixel (default 0)\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -w: WCS type (TAN default)\n");
    fprintf(stderr,"  -x: X and Y coordinates of reference pixel (default is center)\n");
    fprintf(stderr,"  -y: output image dimensions (default is first input image)\n");
    exit (1);
}


static int
RemapImage (ifile, filename)

int	ifile;		/* Sequence number of input file */
char	*filename;	/* FITS or IRAF file filename */

{
    char *image;		/* FITS input image */
    char *imout;		/* FITS output image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    char *irafheader;		/* IRAF image header */
    int nbimage, naxis, naxis1, naxis2, naxis3, naxis4, bytepix;
    int bitpix, nblocks, nbytes;
    double cra, cdec, dra, ddec, secpix;
    char *headout;
    int iraffile;
    int i, j, hpin, wpin, hpout, wpout, nbout;
    int offscl;
    char pixname[128];
    struct WorldCoor *wcsin;
    double bzin, bsin, bzout, bsout;
    double xout, yout, xin, yin, xpos, ypos, dpix;

    /* Read IRAF header and image */
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

    /* Read FITS image */
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

    /* Set input world coordinate system from first image header */
    wcsin = wcsinit (header);

    /* Set output WCS from command line and first image header */
    wcsout = GetFITSWCS (filename, header, verbose, &cra, &cdec, &dra,
			 &ddec, &secpix, &wpout, &hpout, &eqsys, &equinox);

    /* Set output header from command line and first image header */
    strcpy (headout, header);
    hputi4 (headout,"NAXIS", 2);
    hputi4 (headout,"NAXIS1", wpout);
    hputi4 (headout,"NAXIS2", hpout);
    SetFITSWCS (headout, wcsout);

    /* Allocate space for output image */
    hgeti4 (headout, "BITPIX", &bitpix);
    nbout = hpout * wpout * bitpix / 8;
    if (nbout < 0) nbout = -nbout;
    imout = (char * ) calloc (nbout, 1);

    /* Fill output image */
    bsin = 1.0;
    hgetr8 (header, "BSCALE", &bsin);
    bzin = 0.0;
    hgetr8 (header, "BZERO", &bzin);
    bsout = 1.0;
    hgetr8 (header, "BSCALE", &bsout);
    bzout = 0.0;
    hgetr8 (header, "BZERO", &bzout);
    for (i = 1; i <= wpout; i++) {
	xout = (double) i;
	for (j = 1; j <= hpout; j++) {
	    yout = (double) j;
	    pix2wcs (wcsout, xout, yout, &xpos, &ypos);
	    if (!wcsout->offscl) {
		wcs2pix (wcsin, xpos, ypos, &xin, &yin, &offscl);
		if (!offscl) {
		    dpix = getpix (image,bitpix,wpin,hpin,bzin,bsin,xin,yin);
		    putpix (imout,bitpix,wpout,hpout,bzout,bsout,xout,yout);
		    }
		}
	    }
	}

    /* Write output image */
    if (iraffile && !fitsout) {
        if (irafwimage (outname, lhead, irafheader, header, image) > 0 && verbose)
            printf ("%s: written successfully.\n", outname);
        }
    else {
        if (fitswimage (outname, header, image) > 0 && verbose)
            printf ("%s: written successfully.\n", outname);
        }

    free (header);
    free (image);
    wcsfree (wcsin);
    free (headout);
    wcsfree (wcsout);
    free (imout);
    return (0);
}

/* Apr 29 1999	New program
 */
