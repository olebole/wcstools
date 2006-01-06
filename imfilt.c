/* File imfilt.c
 * October 25, 2005
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

#define MEDIAN 1
#define MEAN 2
#define GAUSSIAN 3

static void usage();
static void imFilt();
extern char *FiltFITS();

#define MAXFILES 1000
static int maxnfile = MAXFILES;

static int verbose = 0;		/* verbose/debugging flag */
char outname[128];		/* Name for output image */
static int fitsout = 0;		/* Output FITS file from IRAF input if 1 */
static int overwrite = 0;	/* allow overwriting of input image file */
static int version = 0;		/* If 1, print only program name and version */
static int xsize = 3;
static int ysize = 3;
static int filter = 0;		/* Filter code */
static int nlog = 100;		/* Number of lines between log messages */


int
main (ac, av)
int ac;
char **av;
{
    char *str;
    char c;
    int readlist = 0;
    char *lastchar;
    char filename[128];
    FILE *flist;
    char **fn, *fname;
    char *listfile;
    char *fni;
    int lfn, i, ifile, nfile;
    double dx, dy;

    nfile = 0;
    fn = (char **)calloc (maxnfile, sizeof(char *));
    outname[0] = 0;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    /* Loop through the arguments */
    for (av++; --ac > 0; av++) {
	str = *av;

	/* List of files to be read */
	if (*str == '@') {
	    readlist++;
	    listfile = ++str;
	    str = str + strlen (str) - 1;
	    av++;
	    ac--;
	    str = str - 1;
	    }

	/* Parameters */
	else if (str[0] == '-') {
	    while ((c = *++str) != 0) {
		switch (c) {

		case 'a':	/* Mean filter (average) */
		    if (ac < 3)
			usage (str);
		    filter = MEAN;
		    xsize = (int) atof (*++av);
		    ac--;
		    ysize = (int) atof (*++av);
		    ac--;
		    break;

		case 'g':	/* Gaussian filter */
		    if (ac < 3)
			usage (str);
		    filter = GAUSSIAN;
		    xsize = (int) atof (*++av);
		    ac--;
		    ysize = (int) atof (*++av);
		    ac--;
		    break;

		case 'l': /* Number of lines to log */
		    if (ac < 2)
			usage (str);
		    nlog = (int) atof (*++av);
		    ac--;
		    break;

		case 'm': /* Median filter */
		    if (ac < 3)
			usage (str);
		    filter = MEDIAN;
		    xsize = (int) atof (*++av);
		    ac--;
		    ysize = (int) atof (*++av);
		    ac--;
		    break;

		case 'v':	/* more verbosity */
		    verbose++;
		    break;

		default:
		    usage();
		    break;
		}
		}
	    }

        /* Image file */
        else if (isfits (str) || isiraf (str)) {
            if (nfile >= maxnfile) {
                maxnfile = maxnfile * 2;
                fn = (char **) realloc ((void *)fn, maxnfile);
                }
            fn[nfile] = str;
            nfile++;
            }

        else {
	    fprintf (stderr,"IMFILT: %s is not a FITS or IRAF file \n",str);
            usage();
            }

    }

    /* Process files in file of filenames */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMFILT: List file %s cannot be read\n",
		     listfile);
	    usage ();
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
	    imFilt (filename);
	    if (verbose)
		printf ("\n");
	    }
	fclose (flist);
	}

    /* Process files on command line */
    else if (nfile > 0) {
	for (ifile = 0; ifile < nfile; ifile++) {
	    fname = fn[ifile];
	    lfn = strlen (fname);
	    if (lfn < 8)
		lfn = 8;
	    else {
		if (verbose)
		    printf ("%s:\n", fname);
		imFilt (fname);
		if (verbose)
  		    printf ("\n");
		}
	    }
	}

    /* If no files processed, print usage */
    else
	usage ();

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Filter FITS and IRAF image files\n");
    fprintf(stderr,"Usage: [-v][-a dx[,dy]][-g dx[,dy]][-m dx[,dy]] file.fits ...\n");
    fprintf(stderr,"  -a dx dy: Mean filter dx x dy pixels\n");
    fprintf(stderr,"  -g dx dy: Gaussian filter dx x dy pixels\n");
    fprintf(stderr,"  -l num: Logging interval in lines\n");
    fprintf(stderr,"  -m dx dy: Median filter dx x dy pixels\n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}

static void
imFilt (name)
char *name;
{
    char *image;		/* FITS image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    int iraffile;		/* 1 if IRAF image */
    char *irafheader;		/* IRAF image header */
    char newname[256];		/* Name for revised image */
    char *ext;
    char *newimage;
    char *imext, *imext1;
    char *fname;
    char extname[16];
    int lext, lroot;
    int nx, ny, npix;
    char echar;
    char temp[8];
    char history[64];
    char pixname[256];
    double ctemp;
    double dmax, dpix;
    double bs, bz;

    /* If not overwriting input file, make up a name for the output file */
    if (!overwrite) {
	fname = strrchr (name, '/');
	if (fname)
	    fname = fname + 1;
	else
	    fname = name;
	ext = strrchr (fname, '.');
	if (ext != NULL) {
	    lext = (fname + strlen (fname)) - ext;
	    lroot = ext - fname;
	    strncpy (newname, fname, lroot);
	    *(newname + lroot) = 0;
	    }
	else {
	    lext = 0;
	    lroot = strlen (fname);
	    strcpy (newname, fname);
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
	}
    else
	strcpy (newname, name);

    /* Open IRAF image */
    if (isiraf (name)) {
	iraffile = 1;
	if ((irafheader = irafrhead (name, &lhead)) != NULL) {
	    if ((header = iraf2fits (name, irafheader, lhead, &nbhead)) == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",name);
		free (irafheader);
		return;
		}
	    if ((image = irafrimage (header)) == NULL) {
		hgetm (header,"PIXFIL", 255, pixname);
		fprintf (stderr, "Cannot read IRAF pixel file %s\n", pixname);
		free (irafheader);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", name);
	    return;
	    }
	}

    /* Open FITS file */
    else {
	iraffile = 0;
	if ((header = fitsrhead (name, &lhead, &nbhead)) != NULL) {
	    if ((image = fitsrimage (name, nbhead, header)) == NULL) {
		fprintf (stderr, "Cannot read FITS image %s\n", name);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", name);
	    return;
	    }
	}

    /* Use output filename if it is set on the command line */
    if (outname[0] > 0)
	strcpy (newname, outname);

    /* Create new file name */
    else if (!overwrite) {
	if (imext != NULL) {
	    if (hgets (header, "EXTNAME",8,extname)) {
		strcat (newname, ".");
		strcat (newname, extname);
		}
	    else {
		strcat (newname, "_");
		strcat (newname, imext+1);
		}
	    }
	if (filter == MEDIAN)
	    strcat (newname, "m");
	else if (filter == GAUSSIAN)
	    strcat (newname, "g");
	else
	    strcat (newname, "a");

	if (xsize < 10 && xsize > -1)
	    sprintf (temp,"%1d",xsize);
	else if (xsize < 100 && xsize > -10)
	    sprintf (temp,"%2d",xsize);
	else if (xsize < 1000 && xsize > -100)
	    sprintf (temp,"%3d",xsize);
	else
	    sprintf (temp,"%4d",xsize);
	strcat (newname, temp);
	strcat (newname, "x");
	if (ysize < 10 && ysize > -1)
	    sprintf (temp,"%1d",ysize);
	else if (ysize < 100 && ysize > -10)
	    sprintf (temp,"%2d",ysize);
	else if (ysize < 1000 && ysize > -100)
	    sprintf (temp,"%3d",ysize);
	else
	    sprintf (temp,"%4d",ysize);
	strcat (newname, temp);
	if (fitsout)
	    strcat (newname, ".fits");
	else if (lext > 0) {
	    if (imext != NULL) {
		echar = *imext;
		*imext = (char) 0;
		strcat (newname, ext);
		*imext = echar;
		if (imext1 != NULL)
		    *imext1 = ']';
		}
	    else
		strcat (newname, ext);
	    }
	}
    else
	strcpy (newname, name);

    if (verbose) {
	if (filter == GAUSSIAN)
	    fprintf (stderr,"%d x %d Gaussian Filter ", xsize, ysize);
	else if (filter == MEDIAN)
	    fprintf (stderr,"%d x %d Median Filter ", xsize, ysize);
	else
	    fprintf (stderr,"%d x %d Mean Filter ", xsize, ysize);
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s", name);
	else
	    fprintf (stderr,"FITS image file %s", name);
	fprintf (stderr, " -> %s\n", newname);
	}

    if ((newimage = FiltFITS (header,image,filter,xsize,ysize,nlog)) == NULL) {
	fprintf (stderr,"Cannot filter image %s; file is unchanged.\n",name);
	}
    else {
	if (iraffile && !fitsout) {
	    if (irafwimage (newname,lhead,irafheader,header,newimage) > 0) {
		if (verbose)
		    printf ("%s: written successfully.\n", newname);
		else
		    printf ("%s\n", newname);
		}
	    else if (verbose)
		printf ("IMFILT: File %s not written.\n", newname);
	    }
	else {
	    if (fitswimage (newname, header, newimage) > 0) {
		if (verbose)
		    printf ("%s: written successfully.\n", newname);
		else
		    printf ("%s\n", newname);
		}
	    else if (verbose)
		printf ("IMFILT: File %s not written.\n", newname);
	    }
	free (newimage);
	}

    free (header);
    if (iraffile)
	free (irafheader);
    free (image);
    return;
}
/* Oct 25 2005	New program
 */
