/* File edhead.c
 * May 20, 1998
 * By Doug Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to dmink@cfa.harvard.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>
#include "libwcs/fitshead.h"
#include "libwcs/wcs.h"

static void usage();
static int newimage = 0;
static int verbose = 0;		/* verbose flag */
static void EditHead();

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;

    /* crack arguments */
    for (av++; --ac > 0 && *(str = *av) == '-'; av++) {
	char c;
	while (c = *++str)
	switch (c) {
	case 'v':	/* more verbosity */
	    verbose++;
	    break;
	case 'n':	/* ouput new file */
	    newimage++;
	    break;
	default:
	    usage (progname);
	    break;
	}
    }

    /* now there are ac remaining file names starting at av[0] */
    if (ac == 0)
	usage (progname);

    while (ac-- > 0) {
	char *fn = *av++;
	EditHead (fn);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"Edit header of FITS or IRAF image file\n");
    fprintf(stderr,"%s: usage: [-vn] file.fts file.imh...\n",
	    progname);
    fprintf(stderr,"  -n: write new file, else overwrite \n");
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
EditHead (filename)

char	*filename;	/* FITS or IRAF file filename */

{
    char *image;		/* FITS image */
    char *header;		/* FITS header */
    int lhead;			/* Maximum number of bytes in FITS header */
    int nbhead;			/* Actual number of bytes in FITS header */
    int iraffile;		/* 1 if IRAF image */
    int *irafheader;		/* IRAF image header */
    int i, nbytes, nhb, nhblk, lname, lext;
    char *head, *headend, *hlast;
    char headline[160];
    char newname[128];
    char tempname[128];
    FILE *fd;
    char *ext, *fname;
    char *editcom;
    char newline[1];

    newline[0] = 10;
    strcpy (tempname, "fitshead.temp");

    /* Open IRAF image and header if .imh extension is present */
    if (strsrch (filename,".imh") != NULL) {
	iraffile = 1;
	if ((irafheader = irafrhead (filename, &lhead)) != NULL) {
	    if ((header = iraf2fits (filename, irafheader, lhead, &nbhead)) == NULL) {;
		free (irafheader);
                fprintf (stderr, "Cannot translate IRAF header %s/n", filename);
                return;
                }
	    }

	else {
	    fprintf (stderr, "Cannot read IRAF header file %s\n", filename);
	    free (header);
	    return;
	    }
	}

    /* Read FITS image and header if .imh extension is not present */
    else {
	iraffile = 0;
	if ((header = fitsrhead (filename, &lhead, &nbhead)) != NULL) {
	    if ((image = fitsrimage (filename, nbhead, header)) == NULL) {
		fprintf (stderr, "Cannot read FITS image %s\n", filename);
		free (header);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    return;
	    }
	}
    if (verbose)

    /* Write current header to temporary file */
    if ((fd = fopen (tempname, "w"))) {
	headend = ksearch (header, "END") + 80;
	for (head = header; head < headend; head = head + 80) {
	    for (i = 0; i< 80; i++)
		headline[i] = 0;
	    strncpy (headline,head,80);
	    for (i = 79; i > 0; i--) {
		if (headline[i] == ' ')
		    headline[i] = 0;
		else
		    break;
		}
	    nbytes = i + 1;
	    (void) fwrite (headline, 1, nbytes, fd);
	    (void) fwrite (newline, 1, 1, fd);
	    }
	fclose (fd);
	free (header);
	}
    else {
	fprintf (stderr, "Cannot write temporary header file\n");
	free (header);
	if (iraffile)
	    free (irafheader);
	free (image);
	return;
	}

    /* Run an editor on the temporary header file */
    if (!(editcom = getenv ("EDITOR"))) {
	editcom = (char *)malloc (32);
	strcpy (editcom,"vi");
	}
    strcat (editcom," ");
    strcat (editcom,tempname);
    printf ("Edit command is '%s'\n",editcom);
    if (system (editcom)) {
	free (header);
	if (iraffile)
	    free (irafheader);
	free (image);
	return;
	}

    /* Read the new header from the temporary file */
    if ((fd = fopen (tempname, "r")) != NULL) {
	struct stat buff;
	if (stat (tempname, &buff))
            nbytes = -errno;
	else
            nbytes = ((((int) buff.st_size * 4) / 2880) + 1) * 2880;
	header = (char *) calloc (nbytes, 1);
	head = header;
	hlast = header + nbytes - 1;
	for (i = 0; i< 81; i++)
	    headline[i] = 0;
	while (fgets (headline,82,fd)) {
	    int i = 79;
	    while (headline[i] == 0 || headline[i] == 10)
		headline[i--] = ' ';
	    strncpy (head,headline,80);
	    head = head + 80;
	    if (head > hlast) {
		nhblk = (head - header) / 2880;
		nhb = (nhblk + 10) * 2880;
		header = (char *) realloc (header,nhb);
		head = header + nhb;
		hlast = hlast + 28800;
		}
	    for (i = 0; i< 80; i++)
		headline[i] = 0;
	    }
	fclose (fd);
	}
    else {
	fprintf (stderr, "Cannot read temporary header file %s\n", tempname);
	free (header);
	if (iraffile)
	    free (irafheader);
	free (image);
	return;
	}

    /* Make up name for new FITS or IRAF output file */
    if (newimage) {

    /* Remove directory path and extension from file name */
	ext = strrchr (filename, '.');
	fname = strrchr (filename, '/');
	if (fname)
	    fname = fname + 1;
	else
	    fname = filename;
	lname = strlen (fname);
	if (ext) {
	    lext = strlen (ext);
	    strncpy (newname, fname, lname - lext);
	    *(newname + lname - lext) = 0;
	    }
	else
	    strcpy (newname, fname);

    /* Add file extension preceded by a e */
	if (iraffile)
	    strcat (newname, "e.imh");
	else
	    strcat (newname, "e.fit");
	}
    else
	strcpy (newname, filename);

    /* Write fixed header to output file */
    if (iraffile) {
	if (irafwhead (newname, lhead, irafheader, header) > 0 && verbose)
	    printf ("%s rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
	free (irafheader);
	}
    else {
	if (fitswimage (newname, header, image) > 0 && verbose)
	    printf ("%s: rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
	free (image);
	}

    free (header);
    return;
}

/* Aug 15 1996	New program
 * Aug 26 1996	Change HGETC call to HGETS
 * Aug 27 1996	Read up to 82 characters per line to get newline
 * Aug 29 1996	Allow new file to be written
 * Oct 17 1996	Drop unused variables
 * Dec 11 1996	Allocate editcom if environment variable is not set
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 *
 * May 20 1998	Set reread buffer size based on temporary file size
 */
