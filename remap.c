/* File remap.c
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

#define REMAP_CLOSEST	1
#define REMAP_FLUX	2

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
static char *wcsproj;           /* WCS projection name */

static int nfiles = 0;
static int nbstack = 0;
static int fitsout = 0;
static double secpix = 0;
static double secpix2 = 0;
static int bitpix0 = 0; /* Output BITPIX, =input if 0 */
static FILE *fstack = NULL;
static int RemapImage();
extern struct WorldCoor *GetFITSWCS();
static struct WorldCoor *wcsout = NULL;
static char *headout;		/* FITS output header */
static char *imout;		/* FITS output image */
static int eqsys = 0;
static double equinox = 0.0;
static outsys = 0;
static int version = 0;		/* If 1, print only program name and version */
static int mode = REMAP_CLOSEST;
static int remappix=0;		/* Number of samples of input pixels */
static double xrpix = 0.0;
static double yrpix = 0.0;
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
    double x = 0.0;
    double y = 0.0;
    int lwcs;

    for (i = 0; i < 100; i++)
	filelist[i] = NULL;
    wcsproj = NULL;

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
	    outsys = WCS_B1950;
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
	    outsys = WCS_GALACTIC;
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
	    outsys = WCS_J2000;
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

	/* case 'n':	 Number of samples per linear input pixel
    	    if (ac < 2)
    		usage();
    	    remappix = atoi (*++av);
    	    ac--;
    	    break; */

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
		secpix2 = atof (*++av);
		setsecpix2 (secpix2);
		ac--;
		}
    	    break;

	case 'v':	/* more verbosity */
	    verbose++;
	    break;

	case 'w':	/* Set WCS projection */
	    if (ac < 2)
		usage();
	    wcsproj = *++av;
	    lwcs = strlen (wcsproj);
	    for (i =0; i < lwcs; i++) {
		if (wcsproj[i] > 96)
		    wcsproj[i] = wcsproj[i] - 64;
		}
	    ac--;
	    break;

	case 'x':	/* X and Y coordinates of output image reference pixel */
	    if (ac < 3)
		usage();
	    xrpix = atof (*++av);
	    ac--;
	    yrpix = atof (*++av);
	    ac--;
    	    setrefpix (xrpix, yrpix);
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
    free (headout);
    wcsfree (wcsout);
    free (imout);
    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Remap FITS or IRAF images into single FITS image using WCS\n");
    fprintf(stderr,"Usage: remap [-vf][-a rot][[-b][-j] ra dec][-i bits][-l num] file1.fit file2.fit ... filen.fit\n");
    fprintf(stderr,"  or : remap [-vf][-a rot][[-b][-j] ra dec][-i bits][-l num] @filelist\n");
    fprintf(stderr,"  -a: Output rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -b ra dec: Output center in B1950 (FK4) RA and Dec\n");
    fprintf(stderr,"  -f: Force FITS output\n");
    fprintf(stderr,"  -g long lat: Output center in galactic longitude and latitude\n");
    fprintf(stderr,"  -i num: Number of bits per output pixel (default is input)\n");
    fprintf(stderr,"  -j ra dec: center in J2000 (FK5) RA and Dec\n");
    fprintf(stderr,"  -l num: Log every num rows of output image\n");
    fprintf(stderr,"  -m mode: c closest pixel (more to come)\n");
    /* fprintf(stderr,"  -n: Number of samples per linear input pixel\n"); */
    fprintf(stderr,"  -o name: Name for output image\n");
    fprintf(stderr,"  -p secpix: Output plate scale in arcsec/pixel (default =input)\n");
    fprintf(stderr,"  -v: Verbose\n");
    fprintf(stderr,"  -w type: Output WCS type (input is default)\n");
    fprintf(stderr,"  -x x y: Output image reference X and Y coordinates (default is center)\n");
    fprintf(stderr,"  -y nx ny: Output image dimensions (default is first input image)\n");
    exit (1);
}


static int
RemapImage (ifile, filename)

int	ifile;		/* Sequence number of input file */
char	*filename;	/* FITS or IRAF file filename */

{
    char *image;		/* FITS input image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    char *irafheader;		/* IRAF image header */
    int bitpix, bitpixout;
    double cra, cdec, dra, ddec;
    int iraffile;
    int i, j, ii, jj, hpin, wpin, hpout, wpout, nbout, ix, iy, npout;
    int iin, iout, jin, jout;
    int iout1, iout2, jout1, jout2;
    int offscl, lblock;
    char pixname[256];
    struct WorldCoor *wcsin;
    double bzin, bsin, bzout, bsout;
    double dx, dy, dx0, dy0, secpixin1, secpixin2, secpix1;
    double xout, yout, xin, yin, xpos, ypos, dpixi;
    double xmin, xmax, ymin, ymax;
    double pixratio;
    char history[80];
    char wcstemp[16];

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

    if (!outsys)
	outsys = wcsin->syswcs;

    if (ifile < 1) {

	/* Set output plate scale */
	secpixin1 = wcsin->cdelt[1] * 3600.0;
	secpixin2 = wcsin->cdelt[2] * 3600.0;
	if (secpixin1 < 0)
	    secpixin1 = -secpixin1;
	if (secpix == 0)
	    secpix = secpixin1;
	if (secpix2 == 0)
	    secpix2 = secpix;

	/* Change output dimensions to match output plate scale */
	pixratio = secpixin1 / secpix;
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

	/* Set reference pixel to default of new image center */
	if (nx != 0 && xrpix == 0.0 && yrpix == 0.0) {
	    xrpix = 0.5 * (double) nx;
	    yrpix = 0.5 * (double) ny;
	    setrefpix (xrpix, yrpix);
	    }

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

	if (wcsproj != NULL || outsys != wcsin->syswcs) {
	    if (wcsproj == NULL)
		wcsproj = wcsin->ctype[0]+4;
	    if (outsys == WCS_GALACTIC)
		strcpy (wcstemp, "GLON-");
	    else if (outsys == WCS_ECLIPTIC)
		strcpy (wcstemp, "ELON-");
	    else
		strcpy (wcstemp, "RA---");
	    strcat (wcstemp, wcsproj);
	    hputs  (headout, "CTYPE1", wcstemp);

	    if (outsys == WCS_GALACTIC)
		strcpy (wcstemp, "GLAT-");
	    else if (outsys == WCS_ECLIPTIC)
		strcpy (wcstemp, "ELAT-");
	    else
		strcpy (wcstemp, "DEC--");
	    strcat (wcstemp, wcsproj);
	    hputs  (headout, "CTYPE2", wcstemp);
	    }
	hputr8 (headout, "CRPIX1", xrpix);
	hputr8 (headout, "CRPIX2", yrpix);
	hputr8 (headout, "CDELT1", -secpix/3600.0);
	hputr8 (headout, "CDELT2", secpix2/3600.0);
	if (hgetr8 (headout, "SECPIX1", &secpix1)) {
	    hputr8 (headout, "SECPIX1", secpix);
	    hputr8 (headout, "SECPIX2", secpix2);
	    }
	else if (hgetr8 (headout, "SECPIX", &secpix1)) {
	    if (secpix == secpix2)
		hputr8 (headout, "SECPIX", secpix);
	    else {
		hputr8 (headout, "SECPIX1", secpix);
		hputr8 (headout, "SECPIX2", secpix2);
		}
	    }

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
	    remappix = 1;
	else if (pixratio > remappix) {
	    fprintf (stderr, "REMAP: remapping %.1f pixels from 1; %d too small\n",
		 pixratio, remappix);
	    }

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
	}

    /* Set input WCS output coordinate system to output coordinate system*/
    wcsin->sysout = wcsout->syswcs;
    strcpy (wcsin->radecout, wcsout->radecsys);
    wpin = wcsin->nxpix;
    hpin = wcsin->nypix;

    /* Set output WCS output coordinate system to input coordinate system*/
    wcsout->sysout = wcsin->syswcs;

    sprintf (history, "REMAP input file %s", filename);
    hputc (headout, "HISTORY", history);

    /* Find limiting edges of input image in output image */
    pix2wcs (wcsin, 1.0, 1.0, &xpos, &ypos);
    wcs2pix (wcsout, xpos, ypos, &xout, &yout, &offscl);
    xmin = xout;
    xmax = xout;
    ymin = yout;
    ymax = yout;
    pix2wcs (wcsin, 1.0, (double)hpin, &xpos, &ypos);
    wcs2pix (wcsout, xpos, ypos, &xout, &yout, &offscl);
    if (xout < xmin) xmin = xout;
    if (xout > xmax) xmax = xout;
    if (yout < ymin) ymin = yout;
    if (yout > ymax) ymax = yout;
    pix2wcs (wcsin, (double)wpin, 1.0, &xpos, &ypos);
    wcs2pix (wcsout, xpos, ypos, &xout, &yout, &offscl);
    if (xout < xmin) xmin = xout;
    if (xout > xmax) xmax = xout;
    if (yout < ymin) ymin = yout;
    if (yout > ymax) ymax = yout;
    pix2wcs (wcsin, (double)wpin, (double)hpin, &xpos, &ypos);
    wcs2pix (wcsout, xpos, ypos, &xout, &yout, &offscl);
    if (xout < xmin) xmin = xout;
    if (xout > xmax) xmax = xout;
    if (yout < ymin) ymin = yout;
    if (yout > ymax) ymax = yout;
    iout1 = (int) (ymin + 0.5);
    if (iout1 < 1) iout1 = 1;
    iout2 = (int) (ymax + 0.5);
    if (iout2 > hpout) iout2 = hpout;
    jout1 = (int) (xmin + 0.5);
    if (jout1 < 1) jout1 = 1;
    jout2 = (int) (xmax + 0.5);
    if (jout2 > wpout) jout2 = wpout;

    if (verbose)
	printf ("REMAP: Output x: %d-%d, y: %d-%d\n",
		jout1,jout2,iout1,iout2);

    /* Loop through vertical pixels (output image lines) */
    for (iout = iout1; iout <= iout2; iout++) {
	yout = (double) iout;

	/* Loop through horizontal pixels (output image columns) */
	for (jout = jout1; jout <= jout2; jout++) {
	    xout = (double) jout;

	    /* Get WCS coordinates of this pixel in output image */
	    pix2wcs (wcsout, xout, yout, &xpos, &ypos);

	    /* Get image coordinates of this subpixel in input image */
	    wcs2pix (wcsin, xpos, ypos, &xin, &yin, &offscl);
	    if (!offscl) {
		iin = (int) (yin + 0.5);
		jin = (int) (xin + 0.5);

		/* Read pixel from input file */
		dpixi = getpix1 (image,bitpix,wpin,hpin,bzin,bsin,jin,iin);

		/* Write pixel to output file */
		addpix1 (imout,bitpixout,wpout,hpout,bzout,bsout,jout,iout,dpixi);
		}
	    }

	if (nlog > 0 && iout%nlog == 0)
	    fprintf (stderr,"REMAP: Output image line %04d / %04d filled from %d / %d.\r",
		     iout, iout2, iin, hpin);
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
    return (0);
}

/* Sep 28 1999	New program
 * Oct 22 1999	Drop unused variables after lint
 * Oct 22 1999	Add optional second plate scale argument
 * Nov 19 1999	Add galactic coordinate output option
 * Dec  3 1999	Add option to set output BITPIX
 *
 * Jan 28 2000	Call setdefwcs() with WCS_ALT instead of 1
 * Mar 23 2000	Use hgetm() to get the IRAF pixel file name, not hgets()
 *
 * Apr 13 2001	Set default center to center of output image, not input
 * Jul  2 2001	Fix recentering bug
 * Jul  5 2001	Simplify output by just grabbing the closest input pixel
 * Jul  9 2001	Write projection type into CRVALn of output file
 * Jul 10 2001	Fix -w help message so it correctly notes input proj is default
 * Jul 10 2001	Set CTYPEn correctly for ecliptic or galactic output
 *
 * Apr  9 2002	Do not free unallocated header
 */
