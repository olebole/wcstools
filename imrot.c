/* File imrot.c
 * October 17, 1996
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

#include "libwcs/fitshead.h"

static void usage();
static void imRot ();
extern int RotFITS();

static int verbose = 0;	/* verbose/debugging flag */
static int mirror = 0;	/* reflect image across vertical axis */
static int rotate = 0;	/* rotation in degrees, degrees counter-clockwise */
static int bitpix = 0;	/* number of bits per pixel (FITS code) */
static int fitsout = 0;	/* Output FITS file from IRAF input if 1 */
static int overwrite = 0;	/* allow overwriting of input image file */

main (ac, av)
int ac;
char **av;
{
    char *str;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
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

    	case 'x':	/* Number of bits per pixel */
    	    if (ac < 2)
    		usage ();
    	    bitpix = (int) atof (*++av);
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

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage ();

    else {
	while (ac-- > 0) {
    	    char *fn = *av++;
    	    if (verbose)
    		printf ("%s:\n", fn);
    	    imRot (fn);
    	    if (verbose)
    		printf ("\n");
	    }
	}

    return (0);
}

static void
usage ()
{
    fprintf (stderr,"Rotate and/or Reflect FITS and IRAF image files\n");
    fprintf(stderr,"Usage: [-vm [-r rot] file.fts ...\n");
    fprintf(stderr,"  -f: write FITS image from IRAF input\n");
    fprintf(stderr,"  -l: reflect image across vertical axis\n");
    fprintf(stderr,"  -o: allow overwriting of input image, else write new one\n");
    fprintf(stderr,"  -r: image rotation angle in degrees (default 0)\n");
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
    int *irafheader;		/* IRAF image header */
    char newname[64];		/* Name for revised image */
    char *ext;
    char *fname;
    int lext, lroot;
    int bitpix0;
    char temp[8];
    char history[64];
    char pixname[128];

    /* If not overwriting input file, make up a name for the output file */
    if (!overwrite) {
	ext = strrchr (name, '.');
	fname = strrchr (name, '/');
	if (fname)
	    fname = fname + 1;
	else
	    fname = name;
	if (ext) {
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
	if (rotate < 10 && rotate > -1)
	    sprintf (temp,"%1d",rotate);
	else if (rotate < 100 && rotate > -10)
	    sprintf (temp,"%2d",rotate);
	else if (rotate < 1000 && rotate > -100)
	    sprintf (temp,"%3d",rotate);
	else
	    sprintf (temp,"%4d",rotate);
	if (rotate != 0)
	    strcat (newname, temp);
	if (mirror)
	    strcat (newname, "m");
	if (fitsout)
	    strcat (newname, ".fit");
	else if (lext > 0)
	    strcat (newname, ext);
	}
    else
	strcpy (newname, name);

    /* Open IRAF image if .imh extension is present */
    if (strsrch (name,".imh") != NULL) {
	iraffile = 1;
	irafheader = irafrhead (name, &lhead);
	if (irafheader) {
	    header = iraf2fits (name, irafheader, lhead, &nbhead);
	    if (header == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",name);
		free (irafheader);
		return;
		}
	    image = irafrimage (header);
	    if (image == NULL) {
		hgets (header,"PIXFILE", 64, pixname);
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

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	header = fitsrhead (name, &lhead, &nbhead);
	if (header) {
	    image = fitsrimage (name, nbhead, header);
	    if (image == NULL) {
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
		 fname, bitpix0, bitpix);
	}
    else
	sprintf (history,"New copy of image %s", fname);
    hputc (header,"HISTORY",history);

    if (RotFITS (name, header, &image, rotate, mirror, bitpix, verbose)) {
	fprintf (stderr,"Cannot rotate image %s; file is unchanged.\n",name);
	return;
	}
    if (bitpix != 0)
	hputi4 (header, "BITPIX", bitpix);
    if (iraffile && !fitsout) {
	if (irafwimage (newname, lhead, irafheader, header, image) > 0 && verbose)
	    printf ("%s: written successfully.\n", newname);
	}
    else {
	if (fitswimage (newname, header, image) > 0 && verbose)
	    printf ("%s: written successfully.\n", newname);
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
 */
