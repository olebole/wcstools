/* File remap.c
 * January 28, 2000
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

extern void setcenter();
extern void setsys();
extern void setrot();
extern void setsecpix();
extern void setsecpix2();
extern void setrefpix();
extern void setnpix();
extern void setwcsproj();

static void usage();
static int verbose = 0;		/* verbose flag */
static char *outname0 = "remap.fits";
static char *outname;

static int nfiles = 0;
static int nbstack = 0;
static int fitsout = 0;
static double secpix = 0;
static int bitpix0 = 0; /* Output BITPIX, =input if 0 */
static FILE *fstack = NULL;
static int RemapImage();
extern struct WorldCoor *GetFITSWCS();
static struct WorldCoor *wcsout = NULL;
static int eqsys = 0;
static double equinox = 0.0;
static int version = 0;		/* If 1, print only program name and version */
static int remappix=0;		/* Number of samples of input pixels */
static int nx = 0;
static int ny = 0;
static int nlog = 0;

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
    int ifile, nblocks, nbytes, i;
    char *blanks;
    double x, y;

    for (i = 0; i < 100; i++)
	filelist[i] = NULL;

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
	    
    	case 'g':	/* Output image center on command line in galactic */
    	    if (ac < 3)
    		usage();
	    setsys (WCS_GALACTIC);
	    strcpy (rastr, *++av);
	    ac--;
	    strcpy (decstr, *++av);
	    ac--;
	    setcenter (rastr, decstr);
    	    break;

	case 'i':	/* Bits per output pixel in FITS code */
    	    if (ac < 2)
    		usage();
    	    bitpix0 = atoi (*++av);
    	    ac--;
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

    	case 'l':	/* Logging interval for processing */
    	    if (ac < 2)
    		usage();
    	    nlog = atoi (*++av);
    	    ac--;
    	    break;

	case 'n':	/* Number of samples per linear input pixel */
    	    if (ac < 2)
    		usage();
    	    remappix = atoi (*++av);
    	    ac--;
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
	    secpix = atof (*++av);
    	    setsecpix (secpix);
    	    ac--;
	    if (ac > 1 && isnum (*(av+1))) {
		setsecpix2 (atof (*++av));
		ac--;
		}
    	    break;

	case 'v':	/* more verbosity */
	    verbose++;
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

	case 'z':       /* Use AIPS classic WCS */
	    setdefwcs (WCS_ALT);
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
	    if (verbose)
		printf ("Reading %s\n",filelist[nfiles]);
	    nfiles++;
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
    fprintf(stderr,"usage: remap [-vf][-a rot][[-b][-j] ra dec][-i bits][-l num] file1.fit file2.fit ... filen.fit\n");
    fprintf(stderr,"       remap [-vf][-a rot][[-b][-j] ra dec][-i bits][-l num] @filelist\n");
    fprintf(stderr,"  -a: Output rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b: Output center in B1950 (FK4) RA and Dec\n");
    fprintf(stderr,"  -f: Force FITS output\n");
    fprintf(stderr,"  -g: Output center in galactic longitude and latitude\n");
    fprintf(stderr,"  -i: Number of bits per output pixel (default is input)\n");
    fprintf(stderr,"  -j: Output center in J2000 (FK5) RA and Dec\n");
    fprintf(stderr,"  -l: Log every num rows of input image\n");
    fprintf(stderr,"  -n: Number of samples per linear input pixel\n");
    fprintf(stderr,"  -o: Name for output image\n");
    fprintf(stderr,"  -p: Output plate scale in arcsec/pixel (default =input)\n");
    fprintf(stderr,"  -v: Verbose\n");
    fprintf(stderr,"  -w: Output WCS type (TAN default)\n");
    fprintf(stderr,"  -x: Output image reference X and Y coordinates (default is center)\n");
    fprintf(stderr,"  -y: Output image dimensions (default is first input image)\n");
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
    int bitpix, bitpixout;
    double cra, cdec, dra, ddec;
    char *headout;
    int iraffile;
    int i, j, ii, jj, hpin, wpin, hpout, wpout, nbout, ix, iy, npout;
    int offscl, lblock;
    char pixname[128];
    struct WorldCoor *wcsin;
    double bzin, bsin, bzout, bsout;
    double dx, dy, dx0, dy0, pfrac, secpix0, secpix1;
    double xout, yout, xin, yin, xpos, ypos, dpixi;
    double pixratio, xrpix, yrpix;

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

    /* Change output plate scale to match output dimensions */
    secpix0 = wcsin->cdelt[1] * 3600.0;
    if (secpix0 < 0) secpix0 = -secpix0;
    if (secpix == 0)
	if (nx > 0 && ny > 0) {
	    secpix = secpix0 * (wcsin->nxpix / nx);
    	    setsecpix (secpix);
	    }
	else
	    secpix = secpix0;

    /* Change output dimensions to match output plate scale */
    pixratio = secpix0 / secpix;
    if (nx == 0 && ny == 0) {
	if (secpix > 0) {
	    nx = wcsin->nxpix * pixratio;
	    ny = wcsin->nypix * pixratio;
    	    setnpix (nx, ny);
	    }
	else {
	    nx = wcsin->nxpix;
	    ny = wcsin->nypix;
	    }
	}
    xrpix = wcsin->xrefpix * pixratio;
    yrpix = wcsin->yrefpix * pixratio;
    setrefpix (xrpix, yrpix);

    /* Set output header from command line and first image header */
    lhead = strlen (header);
    lblock = lhead / 2880;
    if (lblock * 2880  < lhead)
	lhead = (lblock+2) * 2880;
    else
	lhead = (lblock+1) * 2880;
    headout = (char *) calloc (lhead, sizeof (char));
    strcpy (headout, header);

    hputi4 (headout, "NAXIS1", nx);
    hputi4 (headout, "NAXIS2", ny);
    hputr8 (headout, "CRPIX1", xrpix);
    hputr8 (headout, "CRPIX2", yrpix);
    hputr8 (headout, "CDELT1", secpix/3600.0);
    hputr8 (headout, "CDELT2", secpix/3600.0);
    if (hgetr8 (headout, "SECPIX1", &secpix1)) {
	hputr8 (headout, "SECPIX1", secpix);
	hputr8 (headout, "SECPIX2", secpix);
	}
    else if (hgetr8 (headout, "SECPIX", &secpix1))
	hputr8 (headout, "SECPIX", secpix);

    hgeti4 (header, "BITPIX", &bitpix);
    if (bitpix0 != 0) {
	hputi4 (headout, "BITPIX", bitpix0);
	bitpixout = bitpix0;
	}
    else
	bitpixout = bitpix;

    /* Set output WCS from command line and first image header */
    wcsout = GetFITSWCS (filename, headout, verbose, &cra, &cdec, &dra,
			 &ddec, &secpix, &wpout, &hpout, &eqsys, &equinox);

    /* Warn if remapping is not acceptable */
    pixratio = wcsin->xinc / wcsout->xinc;
    if (remappix == 0)
	remappix = pixratio + 2;
    else if (pixratio > remappix) {
	fprintf (stderr, "REMAP: remapping %.1f pixels from 1; %d too small\n",
		 pixratio, remappix);
	}

    /* Set output coordinate system from input WCS to output coordinate system*/
    wcsin->sysout = wcsout->syswcs;
    strcpy (wcsin->radecout, wcsout->radecsys);

    /* Allocate space for output image */
    npout = hpout * wpout;
    nbout = npout * bitpixout / 8;
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
    wpin = wcsin->nxpix;
    hpin = wcsin->nypix;

    /* Redistribute flux proportionally from input image */
    pfrac = 1.0 / (float)(remappix * remappix);
    dx = 1.0 / remappix;
    dy = 1.0 / remappix;
    dx0 = 0.5 + (0.5 / remappix);
    dy0 = 0.5 + (0.5 / remappix);

    /* Loop through vertical pixels (image lines) */
    for (i = 1; i <= hpin; i++) {

	/* Loop through horizontal pixels */
	for (j = 1; j <= wpin; j++) {
	    dpixi = pfrac * getpix1 (image,bitpix,wpin,hpin,bzin,bsin,j,i);
	    yin = (double) i - dy0;

	    /* Break up each input pixel vertically */
	    for (ii = 0; ii < remappix; ii++) {
		yin = yin + dy;
		xin = (double) j - dx0;

		/* Break up each input pixel horizontally */
		for (jj = 0; jj < remappix; jj++) {
		    xin = xin + dx;

		    /* Get WCS coordinates of this subpixel in input image */
		    pix2wcs (wcsin, xin, yin, &xpos, &ypos);
		    if (!wcsin->offscl) {

			/* Get pixel coordinates in output image */
			xout = 0.0;
			yout = 0.0;
			wcs2pix (wcsout, xpos, ypos, &xout, &yout, &offscl);
			if (!offscl) {

			    /* Add fraction of input image pixel to output image */
			    ix = (int)(xout + 0.5);
			    iy = (int)(yout + 0.5);
			    addpix1 (imout,bitpixout,wpout,hpout,bzout,bsout,
				    ix, iy, dpixi);
			    }
			}
		    }
		}
	    }
	if (nlog > 0 && i%nlog == 0)
	    fprintf (stderr,"REMAP: Input image line %04d / %04d rebinned to %d / %d.\r",
		     i, hpin, iy, hpout);
	}
    if (nlog > 0)
	printf ("\n");

    /* Write output image */
    if (iraffile && !fitsout) {
        if (irafwimage (outname, lhead, irafheader, headout, imout) > 0 && verbose)
            printf ("%s: written successfully.\n", outname);
        }
    else {
        if (fitswimage (outname, headout, imout) > 0 && verbose)
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

/* Sep 28 1999	New program
 * Oct 22 1999	Drop unused variables after lint
 * Oct 22 1999	Add optional second plate scale argument
 * Nov 19 1999	Add galactic coordinate output option
 * Dec  3 1999	Add option to set output BITPIX
 *
 * Jan 28 2000	Call setdefwcs() with WCS_ALT instead of 1
 */
