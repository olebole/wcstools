/* File imrot.c
 * July 20, 2000
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

static void usage();
static void imRot ();
extern char *RotFITS();

static int verbose = 0;	/* verbose/debugging flag */
static int mirror = 0;	/* reflect image across vertical axis */
static int rotate = 0;	/* rotation in degrees, degrees counter-clockwise */
static int bitpix = 0;	/* number of bits per pixel (FITS code) */
static int fitsout = 0;	/* Output FITS file from IRAF input if 1 */
static int nsplit = 0;	/* Output multiple FITS files from n-extension file */
static int overwrite = 0;	/* allow overwriting of input image file */
static int version = 0;		/* If 1, print only program name and version */

main (ac, av)
int ac;
char **av;
{
    char *str;
    int readlist = 0;
    char *lastchar;
    char filename[128];
    FILE *flist;
    char *listfile;
    char *fni;
    int lfn, i;

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    /* crack arguments */
    for (av++; --ac > 0 && (*(str = *av)=='-' || *str == '@'); av++) {
	char c;
	if (*str == '@')
	    str = str - 1;
	while (c = *++str)
    	switch (c) {
    	case 'f':	/* FITS file output */
    	    fitsout++;
    	    break;

    	case 'l':	/* image flipped around N-S axis */
	    mirror = 1;
    	    break;

	case 'o':	/* allow overwriting of existing image file */
	    overwrite++;
	    break;

    	case 'r':	/* Rotation angle in degrees */
    	    if (ac < 2)
    		usage ();
    	    rotate = (int) atof (*++av);
    	    ac--;
    	    break;

	case 's':	/* split input FITS multiextension image file */
    	    if (ac < 2)
    		usage ();
	    nsplit = atoi (*++av);
	    ac--;
	    break;

    	case 'x':	/* Number of bits per pixel */
    	    if (ac < 2)
    		usage ();
    	    bitpix = (int) atof (*++av);
    	    ac--;
    	    break;

    	case 'v':	/* more verbosity */
    	    verbose++;
    	    break;

	case '@':	/* List of files to be read */
	    readlist++;
	    listfile = ++str;
	    str = str + strlen (str) - 1;
	    av++;
	    ac--;

    	default:
    	    usage();
    	    break;
    	}
    }

    /* Process files in file of filenames */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMROT: List file %s cannot be read\n",
		     listfile);
	    usage ();
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
	    imRot (filename);
	    if (verbose)
		printf ("\n");
	    }
	fclose (flist);
	}

    /* If no arguments left, print usage */
    if (ac == 0)
	usage ();

    /* Process files on command line */
    else {
	while (ac-- > 0) {
    	    char *fn = *av++;
	    lfn = strlen (fn);
	    if (lfn < 8)
		lfn = 8;
	    if (nsplit > 0) {
		lfn = strlen (fn);
		if (lfn < 8)
		    lfn = 8;
		fni = (char *)calloc (lfn+4, 1);
		for (i = 1; i <= nsplit; i++) {
		    sprintf (fni, "%s,%d", fn, i);
    		    if (verbose)
    			printf ("%s:\n", fni);
		    imRot (fni);
		    if (verbose)
  			printf ("\n");
		    }
		free (fni);
		}
	    else {
    		if (verbose)
    		    printf ("%s:\n", fn);
    		imRot (fn);
		if (verbose)
  		    printf ("\n");
		}
	    }
	}

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Rotate and/or Reflect FITS and IRAF image files\n");
    fprintf(stderr,"Usage: [-vm][-r rot][-s num] file.fits ...\n");
    fprintf(stderr,"  -f: write FITS image from IRAF input\n");
    fprintf(stderr,"  -l: reflect image across vertical axis\n");
    fprintf(stderr,"  -o: allow overwriting of input image, else write new one\n");
    fprintf(stderr,"  -r: image rotation angle in degrees (default 0)\n");
    fprintf(stderr,"  -s: split n-extension FITS file\n");
    fprintf(stderr,"  -v: verbose\n");
    fprintf(stderr,"  -x: output pixel size in bits (FITS code, default=input)\n");
    exit (1);
}

static void
imRot (name)
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
    int lext, lroot;
    int bitpix0;
    char echar;
    char temp[8];
    char history[64];
    char pixname[256];

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
	if (imext != NULL) {
	    strcat (newname, "_");
	    strcat (newname, imext+1);
	    }
	if (mirror)
	    strcat (newname, "m");
	if (rotate != 0) {
	    strcat (newname, "r");
	    if (rotate < 10 && rotate > -1)
		sprintf (temp,"%1d",rotate);
	    else if (rotate < 100 && rotate > -10)
		sprintf (temp,"%2d",rotate);
	    else if (rotate < 1000 && rotate > -100)
		sprintf (temp,"%3d",rotate);
	    else
		sprintf (temp,"%4d",rotate);
	    strcat (newname, temp);
	    }
	if (bitpix == -64)
	    strcat (newname, "bn64");
	else if (bitpix == -32)
	    strcat (newname, "bn32");
	else if (bitpix == -16)
	    strcat (newname, "bn16");
	else if (bitpix == 32)
	    strcat (newname, "b32");
	else if (bitpix == 16)
	    strcat (newname, "b16");
	else if (bitpix == 8)
	    strcat (newname, "b8");
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
    if (verbose) {
	fprintf (stderr,"Rotate and/or reflect ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s", name);
	else
	    fprintf (stderr,"FITS image file %s", name);
	fprintf (stderr, " -> %s\n", newname);
	}

    if (bitpix != 0) {
	hgeti4 (header, "BITPIX", &bitpix0);
	if (verbose)
	    fprintf (stderr, "IMROT: %d bits/pixel -> %d bits/pixel\n",
		     bitpix0, bitpix);
	sprintf (history, "New copy of %s BITPIX %d -> %d",
		 name, bitpix0, bitpix);
	}
    else
	sprintf (history,"New copy of image %s", name);
    hputc (header,"HISTORY",history);

    if ((newimage = RotFITS (name,header,image,rotate,mirror,bitpix,verbose))
	== NULL) {
	fprintf (stderr,"Cannot rotate image %s; file is unchanged.\n",name);
	}
    else {
	if (bitpix != 0)
	    hputi4 (header, "BITPIX", bitpix);
	if (iraffile && !fitsout) {
	    if (irafwimage (newname,lhead,irafheader,header,newimage) > 0 &&
		verbose)
		printf ("%s: written successfully.\n", newname);
	    }
	else {
	    if (fitswimage (newname, header, newimage) > 0 && verbose)
		printf ("%s: written successfully.\n", newname);
	    }
	free (newimage);
	}

    free (header);
    if (iraffile)
	free (irafheader);
    free (image);
    return;
}
/* Apr 15 1996	New program
 * Apr 18 1996	Add option to write to current working directory
 * May  2 1996	Pass filename to rotFITS
 * May 28 1996	Change rotFITS to RotFITS
 * Jun  6 1996	Always write to current working directory
 * Jun 14 1996	Use single image buffer
 * Jul  3 1996	Allow optional overwriting of input image
 * Jul 16 1996	Update header reading and allocation
 * Aug 26 1996	Change HGETC call to HGETS; pass LHEAD in IRAFWIMAGE
 * Aug 27 1996	Remove unused variables after lint
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 * Aug 13 1997	Fix bug when overwriting an image
 * Dec 15 1997	Add capability of reading and writing IRAF 2.11 images
 *
 * Feb 24 1998	Add ext. to filename if writing part of multi-ext. file
 * May 26 1998	Fix bug when writing .imh images
 * May 27 1998	Include fitsio.h instead of fitshead.h
 * Jul 24 1998	Make irafheader char instead of int
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Oct 14 1998	Use isiraf() to determine file type
 * Nov 30 1998	Add version and help commands for consistency
 * Jun  8 1999  Return image pointer from RotFITS, not flag
 * Oct 22 1999	Drop unused variables after lint
 *
 * Jan 24 2000	Add to name if BITPIX is changed
 * Mar 23 2000	Use hgetm() to get the IRAF pixel file name, not hgets()
 * Jul 20 2000	Use .fits, not .fit, extension on output FITS file names
 * Jul 20 2000	Add -s option to split multi-extension FITS files
 */
