/* File cphead.c
 * June 8, 2000
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

#define MAXKWD 100
#define MAXFILES 2000
static int maxnkwd = MAXKWD;
static int maxnfile = MAXFILES;

static void usage();
static void CopyValues();
static void strclean();
extern char *GetFITShead();

static int verbose = 0;		/* verbose/debugging flag */
static int nfile = 0;
static int ndec = -9;
static int maxlfn = 0;
static int listall = 0;
static int listpath = 0;
static int newimage0 = 0;
static int keyset = 0;
static int histset = 0;
static int tabout = 0;
static int printhead = 0;
static int version = 0;		/* If 1, print only program name and version */
static int printfill=0;		/* If 1, print ___ for unfound keyword values */
static int printfile=1;		/* If 1, print filename first if >1 files */
static int keyeqval=0;		/* If 1, print keyword=value, not just value */
static char *rootdir=NULL;	/* Root directory for input files */

main (ac, av)
int ac;
char **av;
{
    char *str;
    char **kwd;
    int nkwd = 0;
    char **fn;
    int ifile;
    int lfn;
    char filename[256];
    char *infile;
    char *name;
    FILE *flist, *fdk;
    char *listfile;
    char *ilistfile;
    char *klistfile;
    int ikwd, lkwd, i;
    char *kw, *kwe;
    char string[80];

    ilistfile = NULL;
    klistfile = NULL;
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
	if (*(str = *av)=='-') {
	    char c;
	    while (c = *++str)
	    switch (c) {

		case 'd': /* Root directory for input */
		    if (ac < 2)
			usage();
		    rootdir = *++av;
		    ac--;
		    break;

		case 'h':	/* Set HISTORY */
		    histset++;
		    break;
	
		case 'k':	/* Set CPHEAD keyword */
		    keyset++;
		    break;

		case 'n':	/* Write new file */
		    newimage0++;
		    break;

		case 'p': /* Number of decimal places in output */
		    if (ac < 2)
			usage();
		    ndec = (int) (atof (*++av));
		    ac--;
		    break;
	
		case 'v': /* More verbosity */
		    verbose++;
		    break;

		default:
		    usage();
		    break;
		}
	    }

	/* File containing a list of keywords or files */
	else if (*av[0] == '@') {
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
			kwd = (char **) realloc ((void *)kwd, nkwd);
			maxnkwd = nkwd;
			}
		    if ((fdk = fopen (klistfile, "r")) == NULL) {
			fprintf (stderr,"GETHEAD: File %s cannot be read\n",
				 klistfile);
			nkwd = 0;
			}
		    else {
			for (ikwd = 0; ikwd < nkwd; ikwd++) {
			    kwd[ikwd] = (char *) calloc (32, 1);
			    first_token (fdk, 31, kwd[ikwd]);
			    }
			fclose (fdk);
			}
		    }
		}
	    }

	/* Image file */
	else if (isfits (*av) || isiraf (*av)) {
	    if (infile == NULL)
		infile = *av;
	    else {
		if (nfile >= maxnfile) {
		    maxnfile = maxnfile * 2;
		    fn = (char **) realloc ((void *)fn, maxnfile);
		    }
		fn[nfile] = *av;

		if (listpath || (name = strrchr (fn[nfile],'/')) == NULL)
		    name = fn[nfile];
		else
		    name = name + 1;
		lfn = strlen (name);
		if (lfn > maxlfn)
		    maxlfn = lfn;
		nfile++;
		}
	    }

	/* Keyword */
	else {
	    if (nkwd >= maxnkwd) {
		maxnkwd = maxnkwd * 2;
		kwd = (char **) realloc ((void *)kwd, maxnkwd);
		}
	    kwd[nkwd] = *av;
	    nkwd++;
	    }
	}

    if (nkwd <= 0 && nfile <= 0 )
	usage ();
    else if (nkwd <= 0) {
	fprintf (stderr, "CPHEAD: no keywords specified\n");
	exit (1);
	}
    else if (nfile <= 0 ) {
	fprintf (stderr, "CPHEAD: no files specified\n");
	exit (1);
	}

    if (nkwd > 1)
	printfill = 1;
    if (nfile < 2 && !listall)
	printfile = 0;


    /* Open file containing a list of images, if there is one */
    if (ilistfile != NULL) {
	if ((flist = fopen (ilistfile, "r")) == NULL) {
	    fprintf (stderr,"GETHEAD: Image list file %s cannot be read\n",
		     ilistfile);
	    usage ();
	    }
	}
    if (nfile < 1)
	usage();

    /* Read through headers of images */
    for (ifile = 0; ifile < nfile; ifile++) {
	if (ilistfile != NULL) {
	    first_token (flist, 254, filename);
	    CopyValues (infile, filename, nkwd, kwd);
	    }
	else
	    CopyValues (infile, fn[ifile], nkwd, kwd);

	if (verbose)
	    printf ("\n");
	}
    if (ilistfile != NULL)
	fclose (flist);

    return (0);
}

static void
usage ()
{
    if (version)
	exit (-1);
    fprintf (stderr,"Copy FITS or IRAF header keyword values\n");
    fprintf(stderr,"usage: cphead [-v][-d dir][-p num] file1.fit ... filen.fits kw1 kw2 ... kwn\n");
    fprintf(stderr,"       cphead [-v][-d dir][-p num] file1.fit @filelist kw1 kw2 ... kwn\n");
    fprintf(stderr,"       cphead [-v][-d dir][-p num] file1.fit @filelist @kwlist\n");
    fprintf(stderr,"  -d: Root directory for input files (default is cwd)\n");
    fprintf(stderr,"  -h: Write HISTORY line\n");
    fprintf(stderr,"  -k: Write CPHEAD keyword\n");
    fprintf(stderr,"  -n: Write a new file (add e before the extension)\n");
    fprintf(stderr,"  -p: Number of decimal places in numeric output\n");
    fprintf(stderr,"  -v: Verbose\n");
    exit (1);
}


static void
CopyValues (infile, filename, nkwd, kwd)

char	*infile;	/* FITS or IRAF image file from which to read */
char	*filename;	/* FITS or IRAF image file to which to write */
int	nkwd;		/* Number of keywords for which to print values */
char	*kwd[];		/* Names of keywords for which to print values */

{
    char *header;	/* FITS image header */
    char *headout;	/* FITS image header */
    int newimage;
    char *irafheader;	/* IRAF image header */
    char *image;	/* Input and output image buffer */
    double dval;
    int ival, nch;
    int iraffile;
    char fnform[8];
    char newname[128];
    char string[80];
    char temp[1028];
    char outline[1000];
    char mstring[800];
    char *kw, *kwe, *filepath, *fname, *ext, *imext, *imext1;
    char *kwv, *kwl, *kwv0, *knl;
    int ikwd, lkwd, nfound, notfound, lext, lroot, lhist, lhead;
    char *ltime;
    int naxis, ipos, nbhead, nbr, nbw;
    int fdr, fdw;
    char history[72];
    char echar;
    char *endchar;
    int imageread = 0;

    newimage = newimage0;

    if (rootdir) {
	nch = strlen (rootdir) + strlen (infile) + 2;
	filepath = (char *) calloc (1, nch);
	strcat (filepath, rootdir);
	strcat (filepath, "/");
	strcat (filepath, infile);
	}
    else
	filepath = infile;

    /* Retrieve FITS header from FITS or IRAF .imh file */
    if ((header = GetFITShead (filepath)) == NULL)
	return;

    /* Open IRAF image if .imh extension is present */
    if (isiraf (filename)) {
	iraffile = 1;
	if ((irafheader = irafrhead (filename, &lhead)) != NULL) {
	    if ((headout = iraf2fits (filename, irafheader, lhead, &nbhead)) == NULL) {
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
	if ((headout = fitsrhead (filename, &lhead, &nbhead)) != NULL) {
	    hgeti4 (headout,"NAXIS",&naxis);
	    if (naxis > 0) {
		if ((image = fitsrimage (filename, nbhead, headout)) == NULL) {
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
	fprintf (stderr,"Copy Header Parameter Values from ");
	hgeti4 (header, "IMHVER", &iraffile );
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s to ", infile);
	else
	    fprintf (stderr,"FITS image file %s to ", infile);
	hgeti4 (headout, "IMHVER", &iraffile );
	if (iraffile)
	    fprintf (stderr,"IRAF image file %s\n", filename);
	else
	    fprintf (stderr,"FITS image file %s\n", filename);
	}

    nfound = 0;

    notfound = 0;
    for (ikwd = 0; ikwd < nkwd; ikwd++) {
	lkwd = strlen (kwd[ikwd]);
	kwe = kwd[ikwd] + lkwd;
	for (kw = kwd[ikwd]; kw < kwe; kw++) {
	    if (*kw > 96 && *kw < 123)
		*kw = *kw - 32;
	    }
	if (hgets (header, kwd[ikwd], 80, string)) {
	    strclean (string);
	    if (isnum (string)) {
		if (strchr (string,'.')) {
		    hgetr8 (header, kwd[ikwd], &dval);
		    hputr8 (headout, kwd[ikwd], dval);
		    }
		else {
		    hgeti4 (header, kwd[ikwd], &ival);
		    hputi4 (headout, kwd[ikwd], ival);
		    }
		}
	    else
		hputs (headout, kwd[ikwd], string);

	    if (verbose)
		printf ("%s = %s", kwd[ikwd], string);
	    nfound++;
	    }
	else if (verbose)
	    printf ("%s not found\n", kwd[ikwd]);
	else
	    notfound = 1;
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
	if (hgets (headout, "CPHEAD", 72, history))
	    hputc (headout, "HISTORY", history);
	endchar = strchr (history, ',');
	*endchar = (char) 0;
	strcat (history, " ");
	ltime = lt2fd ();
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
		    hputc (headout, "HISTORY", history);
		    endchar = strchr (history, ',');
		    *endchar = (char) 0;
		    strcat (history, " ");
		    ltime = lt2fd ();
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
	    hputs (headout, "CPHEAD", history);
	if (histset)
	    hputc (headout, "HISTORY", history);
	}

    /* Write fixed header to output file */
    if (iraffile) {
	if (irafwhead (newname, lhead, irafheader, headout) > 0 && verbose)
	    printf ("%s rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
	free (irafheader);
	}
    else if (naxis > 0 && imageread) {
	if (fitswimage (newname, headout, image) > 0 && verbose)
	    printf ("%s: rewritten successfully.\n", newname);
	else if (verbose)
	    printf ("%s could not be written.\n", newname);
	free (image);
	}
    else {
	if ((fdw = fitswhead (newname, headout)) > 0) {
	    fdr = fitsropen (filename);
	    ipos = lseek (fdr, nbhead, SEEK_SET);
	    image = (char *) calloc (2880, 1);
	    while ((nbr = read (fdr, image, 2880)) > 0) {
		nbw = write (fdw, image, nbr);
		if (nbw < nbr)
		    fprintf (stderr,"CPHEAD: %d / %d bytes written\n",nbw,nbr);
		}
	    close (fdr);
	    close (fdw);
	    if (verbose)
		printf ("%s: rewritten successfully.\n", newname);
	    free (image);
	    }
	}

    free (header);
    free (headout);
    return;
}


/* Remove exponent and trailing zeroes, if reasonable */
static void
strclean (string)

char *string;

{
    char *sdot, *s;
    int ndek, lstr, i;

    /* Remove positive exponent if there are enough digits given */
    if (strsrch (string, "E+") != NULL) {
	lstr = strlen (string);
	ndek = (int) (string[lstr-1] - 48);
	ndek = ndek + (10 * ((int) (string[lstr-2] - 48)));
	if (ndek < lstr - 7) {
	    lstr = lstr - 4;
	    string[lstr] = (char) 0;
	    string[lstr+1] = (char) 0;
	    string[lstr+2] = (char) 0;
	    string[lstr+3] = (char) 0;
	    sdot = strchr (string, '.');
	    if (ndek > 0 && sdot != NULL) {
		for (i = 1; i <= ndek; i++) {
		    *sdot = *(sdot+1);
		    sdot++;
		    *sdot = '.';
		    }
		}
	    }
	}

    /* Remove trailing zeroes */
    if (strchr (string, '.') != NULL) {
	lstr = strlen (string);
	s = string + lstr - 1;
	while (*s == '0' && lstr > 1) {
	    if (*(s - 1) != '.') {
		*s = (char) 0;
		lstr --;
		}
	    s--;
	    }
	}

    /* Remove trailing decimal point */
    lstr = strlen (string);
    s = string + lstr - 1;
    if (*s == '.')
	*s = (char) 0;

    return;
	
}

/* Feb 24 2000	New program based on sethead and gethead
 * Mar 22 2000	Use lt2fd() instead of getltime()
 * Jun  8 2000	If no files or keywords specified, say so
 */
