/* File sethead.c
 * October 22, 1999
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
#include "libwcs/fitsfile.h"

#define MAXKWD 50
#define MAXFILES 1000
static int maxnkwd = MAXKWD;
static int maxnfile = MAXFILES;

static void usage();
static void SetValues ();

static int verbose = 0;		/* verbose/debugging flag */
static int newimage0 = 0;
static int keyset = 0;
static int histset = 0;
static int krename = 0;
static char prefix[2];
static int version = 0;		/* If 1, print only program name and version */


main (ac, av)
int ac;
char **av;
{
    char *str;
    char **kwd;
    int nkwd = 0;
    char **fn;
    int nfile = 0;
    int readlist = 0;
    int ifile;
    char filename[128];
    char *keybuff, *kw1, *kw2;
    FILE *flist;
    char *listfile;
    char *ilistfile;
    char *klistfile;
    int ikwd;
    char newline = 10;

    ilistfile = NULL;
    klistfile = NULL;
    keybuff = NULL;
    nkwd = 0;
    nfile = 0;
    fn = (char **)calloc (maxnfile, sizeof(char *));
    kwd = (char **)calloc (maxnkwd, sizeof(char *));

    /* Check for help or version command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help"))
	usage();
    if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	version = 1;
	usage();
	}

    /* crack arguments */
    for (av++; --ac > 0; av++) {
	if (*(str = *av) == '-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

		case 'h':	/* Set HISTORY */
		    histset++;
		    break;
	
		case 'k':	/* Set SETHEAD keyword */
		    keyset++;
		    break;

		case 'n':	/* Write new file */
		    newimage0++;
		    break;

		case 'r':	/* Rename keywords with replaced values */
		    krename++;
		    if (ac > 1) {
			strncpy (prefix, *++av, 1);
			ac--;
			}
		    else
			prefix[0] = 'X';
		    prefix[1] = (char) 0;
		    break;

		case 'v':	/* more verbosity */
		    verbose++;
		    break;

		default:
		    usage();
		    break;
		}
	    }

	/* File containing a list of keywords or files */
	else if (*av[0] == '@') {
	    readlist++;
	    listfile = *av + 1;
	    if (isimlist (listfile)) {
		ilistfile = listfile;
		nfile = getfilelines (ilistfile);
		}
	    else {
		klistfile = listfile;
		nkwd = getfilelines (klistfile);
		if (nkwd > 0) {
		    if (nkwd > maxnkwd) {
			kwd = realloc ((void *)kwd, nkwd);
			maxnkwd = nkwd;
			}
		    keybuff = getfilebuff (klistfile);
		    if (keybuff != NULL) {
			kw1 = keybuff;
			for (ikwd = 0; ikwd < nkwd; ikwd++) {
			    kwd[ikwd] = kw1;
			    kw2 = strchr (kw1, newline);
			    *kw2 = (char) 0;
			    if (ikwd < nkwd - 1)
				kw1 = kw2 + 1;
			    }
			}
		    else
			nkwd = 0;
		    }
		}
	    }

	/* Image file */
	else if (isfits (*av) || isiraf (*av)) {
	    if (nfile >= maxnfile) {
		maxnfile = maxnfile * 2;
		fn = realloc ((void *)fn, maxnfile);
		}
	    fn[nfile] = *av;
	    nfile++;
	    }

	/* Keyword */
	else {
	    if (nkwd >= maxnkwd) {
		maxnkwd = maxnkwd * 2;
		kwd = realloc ((void *)kwd, maxnkwd);
		}
	    kwd[nkwd] = *av;
	    nkwd++;
	    }
	}

    if (nkwd <= 0 || nfile <= 0 )
	usage ();

    /* Open file containing a list of images, if there is one */
    if (ilistfile != NULL) {
	if ((flist = fopen (ilistfile, "r")) == NULL) {
	    fprintf (stderr,"SETHEAD: Image list file %s cannot be read\n",
		     ilistfile);
	    usage ();
	    }
	}

    /* Read through headers of images */
    for (ifile = 0; ifile < nfile; ifile++) {
	if (ilistfile != NULL) {
	    first_token (flist, 254, filename);
	    SetValues (filename, nkwd, kwd);
	    }
	else
	    SetValues (fn[ifile], nkwd, kwd);

	if (verbose)
	    printf ("\n");
	}
    if (ilistfile != NULL)
	fclose (flist);

    if (keybuff != NULL)
	free (keybuff);

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Set FITS or IRAF header keyword values\n");
    fprintf(stderr,"Usage: [-nhkv][-f num][-m num][-r char] file1.fits [... filen.fits] kw1=val1 [ ... kwn=valuen]\n");
    fprintf(stderr,"  or : [-nhkv][-f num][-m num][-r char] file1.fits [... filen.fits] @keywordfile]\n");
    fprintf(stderr,"  or : [-nhkv][-f num][-m num][-r char] @listfile kw1=val1 [ ... kwn=valuen]\n");
    fprintf(stderr,"  or : [-nhkv][-f num][-m num][-r char] @listfile @keywordfile\n");
    fprintf(stderr,"  -h: Write HISTORY line\n");
    fprintf(stderr,"  -k: Write SETHEAD keyword\n");
    fprintf(stderr,"  -n: Write a new file (add e before the extension)\n");
    fprintf(stderr,"  -r: Rename reset keywords with prefix following\n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}


static void
SetValues (filename, nkwd, kwd)

char	*filename;	/* Name of FITS or IRAF image file */
int	nkwd;		/* Number of keywords for which to set values */
char	*kwd[];		/* Names and values of those keywords */

{
    char *header;	/* FITS image header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    char *irafheader;	/* IRAF image header */
    int iraffile;	/* 1 if IRAF image, 0 if FITS image */
    int newimage;	/* 1 to awrite new image file, else 0 */
    int i, lext, lroot;
    char *image;
    char newname[128];
    char *ext, *fname, *imext, *imext1;
    char *kw, *kwv, *kwl, *kwv0;
    char *v, *vq0, *vq1;
    char echar;
    int ikwd, lkwd, lkwv, lhist;
    int fdr, fdw, ipos, nbr, nbw;
    int squote = 39;
    int dquote = 34;
    int naxis = 1;
    int imageread = 0;
    char cval[24];
    char history[72];
    char *endchar;
    char *ltime;
    char newkey[10];
    char value[80];
    char keyroot[8];
    char ctemp;
    int lval, ii, lv;

    newimage = newimage0;

    /* Open IRAF image if .imh extension is present */
    if (isiraf (filename)) {
	iraffile = 1;
	if ((irafheader = irafrhead (filename, &lhead)) != NULL) {
	    if ((header = iraf2fits (filename, irafheader, lhead, &nbhead)) == NULL) {
		fprintf (stderr, "Cannot translate IRAF header %s/n",filename);
		free (irafheader);
		return;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read IRAF file %s\n", filename);
	    return;
	    }
	}

    /* Open FITS file if .imh extension is not present */
    else {
	iraffile = 0;
	if ((header = fitsrhead (filename, &lhead, &nbhead)) != NULL) {
	    hgeti4 (header,"NAXIS",&naxis);
	    if (naxis > 0) {
		if ((image = fitsrimage (filename, nbhead, header)) == NULL) {
		    if (verbose)
			fprintf (stderr, "No FITS image in %s\n", filename);
		    imageread = 0;
		    }
		else
		    imageread = 1;
		}
	    else {
		if (verbose)
		    fprintf (stderr, "Rewriting primary header, copying rest\n");
		newimage = 1;
		}
	    }
	else {
	    fprintf (stderr, "Cannot read FITS file %s\n", filename);
	    return;
	    }
	}
    if (verbose) {
	fprintf (stderr,"Set Header Parameter Values in ");
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", filename);
	else
	    fprintf (stderr,"FITS image file %s\n", filename);
	}

    if (nkwd < 1)
	return;

    /* Set keywords one at a time */
    for (ikwd = 0; ikwd < nkwd; ikwd++) {
        strcpy (cval,"                    ");
	kwv0 = strchr (kwd[ikwd],'=');
	if (kwv0 != NULL) {
	*kwv0 = 0;
	lkwd = kwv0 - kwd[ikwd];
	kwv = kwv0 + 1;
	lkwv = strlen (kwv);

	/* Get current length of header buffer */
	lhead = gethlength (header);

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

	/* Write numeric value to keyword */
	if (isnum (kwv)) {
	    i = 21 - lkwv;
	    for (v = kwv; v < kwv+lkwv; v++)
		cval[i++] = *v;
	    cval[21] = 0;
	    if (hputc (header, kwd[ikwd], cval)) {
		lhead = lhead + 28800;
		if ((header =
		    (char *)realloc(header,(unsigned int)lhead)) != NULL) {
		    hlength (header, lhead);
		    hputc (header,kwd[ikwd], cval);
		    }
		}
	    }

	/* Write boolean value to keyword */
	else if (!strcmp (kwv,"T") || !strcmp (kwv,"t") ||
		 !strcmp (kwv,"YES") || !strcmp (kwv,"yes")) {
	    if (hputl (header, kwd[ikwd], 1)) {
		lhead = lhead + 28800;
		if ((header =
		    (char *)realloc(header,(unsigned int)lhead)) != NULL) {
		    hlength (header, lhead);
		    hputl (header,kwd[ikwd], 1);
		    }
		}
	    }
	else if (!strcmp (kwv,"F") || !strcmp (kwv,"f") ||
		 !strcmp (kwv,"NO") || !strcmp (kwv,"no")) {
	    if (hputl (header, kwd[ikwd], 0)) {
		lhead = lhead + 28800;
		if ((header =
		    (char *)realloc(header,(unsigned int)lhead)) != NULL) {
		    hlength (header, lhead);
		    hputl (header,kwd[ikwd], 0);
		    }
		}
	    }

	/* Write character string to keyword */
	else {
	    if ((vq0 = strchr (kwv,dquote))) {
		vq0 = vq0 + 1;
		vq1 = strchr (vq0,dquote);
		if (vq0 && vq1) {
		    kwv = vq0;
		    *vq1 = 0;
		    }
		}
	    else if ((vq0 = strchr (kwv,squote))) {
		vq0 = vq0 + 1;
		vq1 = strchr (vq0,squote);
		if (vq0 && vq1) {
		    kwv = vq0;
		    *vq1 = 0;
		    }
		}
	    lval = strlen (kwv);
	    if (lval < 69)
		if (hputs (header, kwd[ikwd], kwv)) {
		    lhead = lhead + 14400;
		    if ((header =
			(char *)realloc(header,(unsigned int)lhead)) != NULL) {
			hlength (header, lhead);
			hputs (header, kwd[ikwd], kwv);
			}
		    }

	    /* If character string is longer than 68 characters, split it */
	    else {
		strcpy (keyroot, kwd[ikwd]);
		lroot = strlen (keyroot);
		if (lroot > 6) {
		    *(keyroot+6) = (char) 0;
		    lroot = 6;
		    }
		ii = '1';
		lkwv = strlen (kwv);
		while (lkwv > 0) {
		    if (lkwv > 67)
			lv = 67;
		    else
			lv = lkwv;
		    strncpy (value, kwv, lv);
		    ctemp = value[lv];
		    value[lv] = (char) 0;
		    strcpy (newkey, keyroot);
		    strcat (newkey, "_");
		    newkey[lroot+1] = ii;
		    newkey[lroot+2] = (char) 0;
		    ii++;
		    if (hputs (header, newkey, value)) {
			lhead = lhead + 28800;
			if ((header =
			  (char *)realloc(header,(unsigned int)lhead))!=NULL) {
			    hlength (header, lhead);
			    hputs (header, newkey, value);
			    }
			}
		    value[lv] = ctemp;
		    kwv = kwv + lv;
		    lkwv = lkwv - lv;
		    }
		kwv = kwv0 + 1;
		}
	    }
	if (verbose)
	    printf ("%s = %s\n", kwd[ikwd], kwv);
	*kwv0 = '=';
	}
	}

    /* Make up name for new FITS or IRAF output file */
    if (newimage) {

    /* Remove directory path and extension from file name */
	fname = strrchr (filename, '/');
	if (fname)
	    fname = fname + 1;
	else
	    fname = filename;
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
	if (imext != NULL && *(imext+1) != '0') {
	    strcat (newname, "_");
	    strcat (newname, imext+1);
	    }
	if (fname)
	    fname = fname + 1;
	else
	    fname = filename;
	strcat (newname, "e");
	if (lext > 0) {
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
	strcpy (newname, filename);

    /* Add history to header */
    if (keyset || histset) {
	if (hgets (header, "SETHEAD", 72, history))
	    hputc (header, "HISTORY", history);
	endchar = strchr (history, ',');
	*endchar = (char) 0;
	strcat (history, " ");
	ltime = getltime ();
	strcat (history, ltime);
	endchar = strrchr (history,':');
	*endchar = (char) 0;
	strcat (history, " ");
	for (ikwd = 0; ikwd < nkwd; ikwd++) {
	    kwv0 = strchr (kwd[ikwd],'=');
	    if (kwv0)
		*kwv0 = (char) 0;
	    lhist = strlen (history);
	    lkwd = strlen (kwd[ikwd]);

	    /* If too may keywords, start a second history line */
	    if (lhist + lkwd + 10 > 71) {
		if (histset) {
		    strcat (history, " updated");
		    hputc (header, "HISTORY", history);
		    endchar = strchr (history, ',');
		    *endchar = (char) 0;
		    strcat (history, " ");
		    ltime = getltime ();
		    strcat (history, ltime);
		    endchar = strrchr (history,':');
		    *endchar = (char) 0;
		    strcat (history, " ");
		    }
		else
		    break;
		}
	    strcat (history, kwd[ikwd]);
	    if (kwv0)
		*kwv0 = '=';
	    if (nkwd == 2 && ikwd < nkwd-1)
		strcat (history, " and ");
	    else if (ikwd < nkwd-1)
		strcat (history, ", ");
	    }
	strcat (history, " updated");
	if (keyset)
	    hputs (header, "SETHEAD", history);
	if (histset)
	    hputc (header, "HISTORY", history);
	}

    /* Write fixed header to output file */
    if (iraffile) {
	if (irafwhead (newname, lhead, irafheader, header) > 0 && verbose)
	    printf ("%s rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
	free (irafheader);
	}
    else if (naxis > 0 && imageread) {
	if (fitswimage (newname, header, image) > 0 && verbose)
	    printf ("%s: rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
	free (image);
	}
    else {
	if ((fdw = fitswhead (newname, header)) > 0) {
	    fdr = fitsropen (filename);
	    ipos = lseek (fdr, nbhead, SEEK_SET);
	    image = (char *) calloc (2880, 1);
	    while ((nbr = read (fdr, image, 2880)) > 0) {
		nbw = write (fdw, image, nbr);
		if (nbw < nbr)
		    fprintf (stderr,"SETHEAD: %d / %d bytes written\n",nbw,nbr);
		}
	    close (fdr);
	    close (fdw);
	    if (verbose)
		printf ("%s: rewritten successfully.\n", newname);
	    free (image);
	    }
	}

    free (header);
    return;
}

/* Oct 11 1996	New program
 * Dec 12 1996	Move ISNUM subroutine to hget.c
 *
 * Feb 21 1997  Check pointers against NULL explicitly for Linux
 *
 * May 28 1998	Include fitsio.h instead of fitshead.h
 * Jun  2 1998	Fix bug in hput()
 * Jul 24 1998	Deal coorectly with logical T or F
 * Jul 24 1998	Make irafheader char instead of int
 * Jul 30 1998	Allow use of list of files and multiple files
 * Aug  6 1998	Change fitsio.h to fitsfile.h
 * Aug 14 1998	Preserve extension when creating new file name
 * Aug 14 1998	If changing primary header, write out entire input file
 * Aug 31 1998	Add options to add HISTORY and/or SETHEAD keyword
 * Sep  1 1998	Add option to keep changed keywords with new names
 * Oct  5 1998	Allow header changes even if no data is present
 * Oct  5 1998	Determine assignment arguments by presence of equal sign
 * Oct  5 1998	Use isiraf() to determine if file is IRAF or FITS
 * Oct 29 1998	Fix history setting
 * Nov 30 1998	Add version and help commands for consistency
 * Dec 30 1998	Write header without image if no image is present
 *
 * Mar  4 1999	Reset = for each keyword after setting value in header
 * Apr  2 1999	Add warning if too many files or keywords on command line
 * Apr  2 1999	Add -f and -m to change maximum number of files or keywords
 * Jul 14 1999  Read lists of BOTH keywords and files simultaneously
 * Jul 14 1999  Reallocate keyword array if too many in file
 * Jul 15 1999	Add capability of writing multi-line keywords a la IRAF
 * Jul 15 1999	Reallocate keyword and file lists if default limits exceeded
 * Oct 14 1999	Reallocate header if length is exceeded
 * Oct 22 1999	Drop unused variables after lint
 */
