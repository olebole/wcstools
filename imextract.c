/* File imextract.c
 * May 18, 1999
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
#include "wcscat.h"

#define MAXKWD 50

static void usage();
static int ExtractImage();

static int verbose = 0;		/* verbose flag */
static int wfits = 0;		/* if 1, write FITS header before data */
static char *newname = "imstack.out";
static int krename = 0;
static char prefix[2];
static int nfiles = 0;
static int nbstack = 0;
static int version = 0;		/* If 1, print only program name and version */
static int fitsout = 0;		/* If 1, write FITS output for IRAF input */
static int nameout = 0;		/* If 1, write output file name */
static char *suffix = NULL;	/* Suffix if set on command line */
static char *outfile = NULL;	/* Output file name if set on command line */
static char *outdir = NULL;	/* Output directory if set on command line */

main (ac, av)
int ac;
char **av;
{
    char filename[128];
    char *filelist[100];
    char *listfile;
    char *str;
    char *kwd[MAXKWD];
    int nkwd = 0;
    int readlist = 0;
    char *temp;
    FILE *flist;
    int ifile, nblocks, nbytes, i, nbw;
    char *blanks;
    char *ranges = NULL;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    /* Crack arguments */
    for (av++; --ac > 0; av++) {
	char c;

	str = *av;

	/* Print help message */
	if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	    usage();

	/* Print version message */
	else if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	    version = 1;
	    usage();
	    }

	/* If filename preceded by @, it is a list of filenames */
	else if (*str == '@') {
	    readlist++;
	    listfile = ++str;
	    continue;
	    }

	/* Set range and make a list of extraction numbers from it */
	else if (strchr (*av + 1, '-') || strchr (*av + 1, ',')) {
	    if (ranges) {
		temp = ranges;
		ranges = (char *) calloc (strlen(ranges) + strlen(*av) + 2, 1);
		strcpy (ranges, temp);
		strcat (ranges, ",");
		strcat (ranges, *av);
		free (temp);
		}
	    else {
		ranges = (char *) calloc (strlen(*av) + 1, 1);
		strcpy (ranges, *av);
		}
	    continue;
	    }

	/* If numeric argument, set image to be extracted */
	else if (isnum (str)) {
	    if (ranges) {
		temp = ranges;
		ranges = (char *)calloc (strlen(ranges)+strlen(*av)+2, 1);
		strcpy (ranges, temp);
		strcat (ranges, ",");
		strcat (ranges, *av);
		free (temp);
		}
	    else {
		ranges = (char *) calloc (strlen(*av) + 1, 1);
		strcpy (ranges, *av);
		}
	    continue;
	    }

	/* If equal sign in argument, it is a header keyword assignment */
	else if (strsrch (*av,"=") != NULL) {
	    kwd[nkwd] = *av;
	    if (nkwd < MAXKWD)
		nkwd++;
	    continue;
	    }

	/* Otherwise, if there is no preceding -, it is a file name */
	else if (*str != '-') {
	    filelist[nfiles] = *av;
	    nfiles++;
	    if (verbose)
		printf ("Reading %s\n",filelist[nfiles-1]);
	    continue;
	    }
	while (c = *++str) {
	    switch (c) {
		case 'd':	/* Write to specific output directory */
		    if (ac < 2)
			usage();
		    outdir = *++av;
		    ac--;
		    break;
		case 'f':	/* FITS output for IRAF input */
		    fitsout++;
		    break;
		case 'n':	/* Echo output file name */
		    nameout++;
		    break;
		case 'o':	/* Write specific output file name */
		    if (ac < 2)
			usage();
		    outfile = *++av;
		    ac--;
		    break;
		case 'v':	/* more verbosity */
		    verbose++;
		    break;
		case 'x':	/* add specific suffix before extension */
		    if (ac < 2)
			usage();
		    suffix = *++av;
		    ac--;
		    break;
		default:
		    usage();
		    break;
		}
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

    if (nfiles < 1)
	usage();

    /* Extract images */
    for (ifile = 0;  ifile < nfiles; ifile++) {
	if (readlist) {
	    if (fgets (filename, 128, flist) != NULL) {
		filename[strlen (filename) - 1] = 0;
		ExtractImage (filename, ranges, ifile, nkwd, kwd);
		}
	    }
	else
	    ExtractImage (filelist[ifile], ranges, ifile, nkwd, kwd);
	}

    if (readlist)
	fclose (flist);

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Extract one FITS or IRAF image from compound image\n");
    fprintf(stderr,"usage: imextract [-vf] [-x ext] [kwn=valn] range file1.fit ... filen.fit\n");
    fprintf(stderr,"       imextract [-vf] [-o file] [kwn=valn] n filename\n");
    fprintf(stderr,"       imextract [-vf] [kwn=valn] n @filelist\n");
    fprintf(stderr,"  range: images to extract (by sequence number)\n");
    fprintf(stderr,"  -f: Write FITS out for IRAF input\n");
    fprintf(stderr,"  -n: Write out name of output file\n");
    fprintf(stderr,"  -o: Specify output file name (without extension) \n");
    fprintf(stderr,"  -v: Verbose\n");
    fprintf(stderr,"  -x: Add this extension instead of _n\n");
    exit (1);
}


static int
ExtractImage (filename, ranges, ifile, nkwd, kwd)

char	*filename;	/* FITS or IRAF file filename */
char	*ranges;	/* String with range of sequence numbers to extract */
int	ifile;		/* Number in list of files on which to operate */
int	nkwd;		/* Number of keywords for which to set values */
char	*kwd[];		/* Names and values of those keywords */

{
    char *image;		/* FITS image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    char *irafheader=NULL;	/* IRAF image header */
    char *outimage;
    int nbimage, naxis, naxis1, naxis2, naxis3, naxis4, bytepix;
    int bitpix, nblocks, nbytes;
    int iraffile;
    int i, lext, lroot, nbskip;
    char *kw, *kwv, *kwl, *kwv0;
    char *v, *vq0, *vq1;
    int ikwd, lkwd, lkwv;
    char pixname[128];
    char *fname, *ext, *imext, *imext1;
    char echar;
    char temp[64];
    struct Range *range; /* Range of sequence numbers to list */
    char newname[128];
    char newkey[10];
    char cval[24];
    int squote = 39;
    int dquote = 34;
    int nimages, nimage;
    int nidef = 1;

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
    if (verbose && !ifile)

    /* Compute size of input image in bytes from header parameters */
    hgeti4 (header,"NAXIS",&naxis);
    if (naxis == 1) {
	printf ("IMEXTRACT: Image %s has only one dimension\n", filename);
	return;
	}
    naxis1 = 1;
    hgeti4 (header,"NAXIS1",&naxis1);
    naxis2 = 1;
    hgeti4 (header,"NAXIS2",&naxis2);
    naxis3 = 1;
    hgeti4 (header,"NAXIS3",&naxis3);
    if (naxis1 * naxis2 == 1 || naxis1 * naxis3 == 1 || naxis2 * naxis3 == 1) {
        printf ("IMEXTRACT: Image %s has only one real dimension\n", filename);
        return;
        }
    hgeti4 (header,"BITPIX",&bitpix);
    bytepix = bitpix / 8;
    if (bytepix < 0) bytepix = -bytepix;

    /* Remove directory path and extension from file name */
    fname = strrchr (filename, '/');
    if (fname)
	fname = fname + 1;
    else
	fname = filename;
    ext = strrchr (fname, '.');

    /* Set output filename */
    if (ext != NULL) {
	lext = (fname + strlen (fname)) - ext;
	lroot = ext - fname;
	}
    else {
	lext = 0;
	lroot = strlen (fname);
	}
    imext = strchr (fname, ',');
    imext1 = NULL;
    if (imext == NULL) {
	imext = strchr (fname, '[');
	if (imext != NULL) {
	    imext1 = strchr (fname, ']');
	    *imext1 = (char) 0;
	    }
	}
    if (imext != NULL) {
	ranges = (char *) calloc (strlen(imext+1) + 1, 1);
	strcpy (ranges, imext);
	}

    /* Figure out how much to write out and where to start */
    range = RangeInit (ranges, nidef);
    nimages = rgetn (range);

    for (i = 0; i < nimages; i++) {
	nimage = rgeti4 (range);

    /* Set output directory path if given on command line */
    if (outdir != NULL) {
	strcpy (newname, outdir);
	strcat (newname, "/");
	}
    else
	newname[0] = (char) 0;

    /* Make up name for new FITS or IRAF output file */
    if (ext != NULL) {
	strncat (newname, fname, lroot);
	*(newname + lroot) = 0;
	}
    else
	strcat (newname, fname);

    /* If no section selected, copy entire image with header unchanged */
    if (nimage < 1)
	outimage = image;

    /* If only 2 axes, drop 2nd from output file header */
    else if (naxis < 3) {
	nbskip = (nimage - 1) * (naxis1 * bytepix);
	outimage = image + nbskip;
	hputi4 (header, "NAXIS", 1);
	hdel (header,"NAXIS2");
	}

    /* If only 2 populated axes, drop 2nd and 3rd from output file header */
    else if (naxis2 ==1 || naxis3 == 1) {
	nbskip = (nimage - 1) * (naxis1 * bytepix);
	outimage = image + nbskip;
	hputi4 (header, "NAXIS", 1);
	hdel (header,"NAXIS2");
	hdel (header,"NAXIS3");
	}

    /* If 3 populated axes, drop 3rd from output file header */
    else {
	nbskip = (nimage - 1) * (naxis1 * naxis2 * bytepix);
	outimage = image + nbskip;
	hputi4 (header, "NAXIS", 2);
	hdel (header,"NAXIS3");
	}

    /* Set keywords one at a time */
    for (ikwd = 0; ikwd < nkwd; ikwd++) {
        strcpy (cval,"                    ");
	kwv0 = strchr (kwd[ikwd],'=');
	*kwv0 = 0;
	lkwd = kwv0 - kwd[ikwd];
	kwv = kwv0 + 1;
	lkwv = strlen (kwv);

	/* Make keyword all upper case */
	kwl = kwd[ikwd] + lkwd;
	for (kw = kwd[ikwd]; kw < kwl; kw++) {
	    if (*kw > 96 && *kw < 123)
		*kw = *kw - 32;
	    }

	/* If keyword is already in header, krename it if requested */
	if (krename && ksearch (header, kwd[ikwd])) {
	    strcpy (newkey, prefix);
	    strcat (newkey, kwd[ikwd]);
	    if (strlen (newkey) > 8)
		newkey[8] = (char) 0;
	    hchange (header, kwd[ikwd], newkey);
	    }

	/* Write value to keyword */
	if ((vq0 = strchr (kwv,dquote))) {
	    vq0 = vq0 + 1;
	    vq1 = strchr (vq0,dquote);
	    if (vq0 && vq1) {
		kwv = vq0;
		*vq1 = 0;
		hputs (header, kwd[ikwd], kwv);
		}
	    else
		hputs (header, kwd[ikwd], kwv);
	    }
	else if ((vq0 = strchr (kwv,squote))) {
	    vq0 = vq0 + 1;
	    vq1 = strchr (vq0,squote);
	    if (vq0 && vq1) {
		kwv = vq0;
		*vq1 = 0;
		hputs (header, kwd[ikwd], kwv);
		}
	    else
		hputs (header, kwd[ikwd], kwv);
	    }
	else if (isnum (kwv)) {
	    i = 21 - lkwv;
	    for (v = kwv; v < kwv+lkwv; v++)
		cval[i++] = *v;
	    cval[21] = 0;
	    hputc (header, kwd[ikwd], cval);
	    }
	else if (!strcmp (kwv,"T") || !strcmp (kwv,"t"))
	    hputl (header, kwd[ikwd], 1);
	else if (!strcmp (kwv,"YES") || !strcmp (kwv,"yes"))
	    hputl (header, kwd[ikwd], 1);
	else if (!strcmp (kwv,"F") || !strcmp (kwv,"f"))
	    hputl (header, kwd[ikwd], 0);
	else if (!strcmp (kwv,"NO") || !strcmp (kwv,"no"))
	    hputl (header, kwd[ikwd], 0);
	else
	    hputs (header, kwd[ikwd], kwv);
	if (verbose)
	    printf ("%s = %s\n", kwd[ikwd], kwv);
	*kwv0 = '=';
	}

    /* Drop multispec suffix if output is 1-D file */
    hgeti4 (header, "NAXIS", &naxis);
    if (naxis == 1) {
	if ((imext1 = strstr (newname, ".ms")) != NULL) {
	    *imext1 = (char) 0;
	    *(imext1+1) = (char) 0;
	    *(imext1+2) = (char) 0;
	    }
	}

    /* Add suffix from command line, if one is present */
    if (nimage > 0 && suffix != NULL) {
	if (strlen (suffix) > 0) {
	    strcat (newname, "_");
	    strcat (newname, suffix);
	    }
	}

    /* Add suffix from file name image section, if one is present */
    else if (nimage > 0 && imext != NULL) {
	strcat (newname, "_");
	strcat (newname, imext+1);
	}

    /* Create suffix from number of extracted sub-image */
    else if (nimage > 0) {
	strcat (newname, "_");
	sprintf (temp,"%03d",nimage);
	strcat (newname, temp);
	}

    if (outfile != NULL)
	strcpy (newname, outfile);

    /* Create output IRAF file if input was an IRAF file */
    if (iraffile && !fitsout) {
	strcpy (pixname, newname);
	strcat (pixname, ".pix");
	hputs (header, "PIXFILE", pixname);
	strcat (newname, ".imh");
	}
    else if (lext > 0) {
	if (imext != NULL) {
	    echar = *imext;
	    *imext = (char) 0;
	    if (fitsout)
		strcat (newname, ".fits");
	    else
		strcat (newname, ext);
	    *imext = echar;
	    if (imext1 != NULL)
		*imext1 = ']';
	    }
	else if (fitsout)
	    strcat (newname, ".fits");
	else
	    strcat (newname, ext);
	}

    /* Write new IRAF or FITS file */
    if (iraffile && !fitsout) {
	if (!irafwimage (newname, lhead, irafheader, header, outimage ))
	    printf ("IRAF file %s not written\n", newname);
	else if (verbose)
	    printf ("IRAF file %s written from %s\n", newname, filename);
	}
    else {
	if (!fitswimage (newname, header, outimage))
	    printf ("FITS file %s not written\n", newname);
	else if (verbose)
	    printf ("FITS file %s written from %s\n", newname, filename);
	}
    if (nameout)
	printf ("%s\n", newname);

    }

    if (irafheader)
	free (irafheader);
    free (header);
    free (image);
    return (0);
}

/* Feb 11 1999	New program
 * Mar  2 1999	Add option to echo output file name
 * Mar  4 1999	Allow keyword assignments from command line
 * Mar  5 1999	Add command line output file name specification
 * Mar  8 1999	Extract range of images in one pass
 * Mar  9 1999	Add option to write to a specific directory
 * May 18 1999	If suffix is null string, do not add _
 */