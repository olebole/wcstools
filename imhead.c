/* File imhead.c
 * March 16, 1998
 * By Doug Mink Harvard-Smithsonian Center for Astrophysics)
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
static int PrintFITSHead();
static void PrintHead();
extern char *GetFITShead();

static int nskip = 0;		/* Number of bytes to skip */
static int nfiles = 0;		/* Nuber of files for headers */
static int verbose = 0;		/* verbose/debugging flag */

main (ac, av)
int ac;
char **av;
{
    char *progname = av[0];
    char *str;
    int readlist = 0;
    char *lastchar;
    char filename[128];
    FILE *flist;
    char *listfile;

    /* crack arguments */
    for (av++; --ac > 0 && (*(str = *av)=='-' || *str == '@'); av++) {
	char c;
	if (*str == '@')
	    str = str - 1;
	while (c = *++str)
	switch (c) {

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
	    usage(progname);
	    break;
	}
    }

    /* Find number of images to search and leave listfile open for reading */
    if (readlist) {
	if ((flist = fopen (listfile, "r")) == NULL) {
	    fprintf (stderr,"IMHEAD: List file %s cannot be read\n",
		     listfile);
	    usage (progname);
	    }
	while (fgets (filename, 128, flist) != NULL) {
	    lastchar = filename + strlen (filename) - 1;
	    if (*lastchar < 32) *lastchar = 0;
	    PrintHead (filename);
	    if (verbose)
		printf ("\n");
	    }
	fclose (flist);
	}

    /* If no arguments left, print usage */
    if (ac == 0)
	usage (progname);

    nfiles = ac;
    while (ac-- > 0) {
	char *fn = *av++;
	PrintHead (fn);
	if (verbose)
	    printf ("\n");
	}

    return (0);
}

static void
usage (progname)
char *progname;
{
    fprintf (stderr,"Print FITS or IRAF image header\n");
    fprintf(stderr,"%s: usage: [-v][-b num] file.fit ...\n", progname);
    fprintf(stderr,"  -v: verbose\n");
    exit (1);
}


static void
PrintHead (name)

char *name;

{
    char *header;	/* FITS image header */

    if ((header = GetFITShead (name)) == NULL)
	return;

    if (verbose)

    if (verbose || nfiles > 1) {
	if (strsrch (name,".imh") != NULL)
	    printf ("%s IRAF file header:\n", name);
	else
	    printf ("%s FITS file header:\n", name);
	}

    if (PrintFITSHead (header) && verbose)
	printf ("%s: no END of header found\n", name);

    free (header);
    return;
}


static int
PrintFITSHead (header)

char	*header;	/* Image FITS header */
{
    char line[80], *iline, *endhead;
    int i, nblank;

    endhead = ksearch (header, "END") + 80;
    if (endhead == NULL)
	return (1);

    nblank = 0;
    for (iline = header; iline < endhead; iline = iline + 80) {
	strncpy (line, iline, 80);
	i = 79;
	while (line[i] <= 32 && i > 0)
	    line[i--] = 0;
	if (i > 0) {
	    if (nblank > 1) {
		printf ("COMMENT   %d blank lines\n",nblank);
		nblank = 0;
		}
	    else if (nblank > 0) {
		printf ("COMMENT   %d blank line\n",nblank);
		nblank = 0;
		}
	    printf ("%s\n",line);
	    }
	else
	    nblank++;
	}

    return (0);
}
/* Jul 10 1996	New program
 * Jul 16 1996	Update header I/O
 * Aug 15 1996	Drop unnecessary reading of FITS image; clean up code
 * Aug 27 1996	Drop unused variables after lint
 * Nov 19 1996	Add linefeeds after filename in verbose mode
 * Dec  4 1996	Print "header" instead of "WCS" in verbose mode
 * Dec 17 1996	Add byte skipping before header
 *
 * Feb 21 1997  Get header from subroutine
 * May 28 1997  Add option to read a list of filenames from a file
 * Dec 12 1997	Read IRAF version 2 .imh files
 * Dec 15 1997	Note number of blank lines in header as comment
 *
 * Jan  5 1998	Print file name if multiple headers printed
 * Jan  5 1998	Print error message if no END is found in header
 * Jan 14 1998	Really get IRAF 2.11 files right on any architecture
 * Mar 16 1998	Print line instead of lines if there is only one blank line
 */
