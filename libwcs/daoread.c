/*** File libwcs/daoread.c
 *** March 20, 1997
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 */

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>

static int nlines;	/* Number of lines in catalog */
#define ABS(a) ((a) < 0 ? (-(a)) : (a))

static int daoopen();
static char *daoline();

static char newline = 10;
static char *daobuff;

/* DAOREAD -- Read DAOFIND file of star positions in an image */

int
daoread (daocat, xa, ya, ba, pa, nlog)

char	*daocat;	/* Name of DAOFIND catalog file */
double	**xa, **ya;	/* X and Y coordinates of stars, array returned */
double	**ba;		/* Fluxes of stars in counts, array returned */
int	**pa;		/* Peak counts of stars in counts, array returned */
int	nlog;		/* 1 to print each star's position */
{
    int nstars;
    double xi, yi, magi;
    double flux;
    int iline;
    char *line;

    line = 0;
    nstars = 0;

    if (daoopen (daocat) > 0) {
	line = daobuff;

    /* Loop through catalog */
	for (iline = 1; iline <= nlines; iline++) {
	    line = daoline (iline, line);
	    if (line == NULL) {
		fprintf (stderr,"DAOREAD: Cannot read line %d\n", iline);
		break;
		}
	    else if (line[0] != '#') {

		/* Extract X, Y, magnitude  */
		sscanf (line,"%lg %lg %lg", &xi, &yi, &magi);

		/* Save star position, scaled flux, and magnitude in table */
		nstars++;
		*xa= (double *) realloc(*xa, nstars*sizeof(double));
		*ya= (double *) realloc(*ya, nstars*sizeof(double));
		*ba= (double *) realloc(*ba, nstars*sizeof(double));
		*pa= (int *) realloc(*pa, nstars*sizeof(int));
		(*xa)[nstars-1] = xi;
		(*ya)[nstars-1] = yi;
		flux = 10000.0 * pow (10.0, (-magi / 2.5));
		(*ba)[nstars-1] = flux;
		(*pa)[nstars-1] = (int)(magi * 100.0);

		if (nlog == 1)
		    fprintf (stderr,"DAOREAD: %6d: %9.5f %9.5f %15.2f %6.2f\n",
			   nstars,xi,yi,flux,magi);
		}

	    /* Log operation */
	    if (nlog > 0 && iline%nlog == 0)
		fprintf (stderr,"DAOREAD: %5d / %5d / %5d stars from catalog %s\r",
			nstars, iline, nlines, daocat);

	    /* End of star loop */
	    }

	/* End of open catalog file */
	}

/* Summarize search */
    if (nlog > 0)
	fprintf (stderr,"DAOREAD: Catalog %s : %d / %d / %d found\n",
		 daocat, nstars, iline, nlines);

    free (daobuff);

    return (nstars);
}


/* DAOOPEN -- Open DAOFIND catalog, returning number of entries */

static int
daoopen (daofile)

char *daofile;	/* DAOFIND catalog file name */
{
    struct stat statbuff;
    FILE *fcat;
    int nr, lfile;
    char *daonew;
    
/* Find length of DAOFIND catalog */
    if (stat (daofile, &statbuff)) {
	fprintf (stderr,"DAOOPEN: DAOFIND catalog %s has no entries\n",daofile);
	return (0);
	}
    else
	lfile = (int) statbuff.st_size;

/* Open DAOFIND catalog */
    if (!(fcat = fopen (daofile, "r"))) {
	fprintf (stderr,"DAOOPEN: DAOFIND catalog %s cannot be read\n",daofile);
	return (0);
	}

/* Allocate buffer to hold entire catalog and read it */
    if ((daobuff = malloc (lfile)) != NULL) {
	nr = fread (daobuff, 1, lfile, fcat);
	if (nr < lfile) {
	    fprintf (stderr,"DAOOPEN: read only %d / %d bytes of file %s\n",
		     nr, lfile, daofile);
	    (void) fclose (fcat);
	    return (0);
	    }

    /* Enumerate entries in DAOFIND catalog by counting newlines */
	daonew = daobuff;
	nlines = 0;
	while ((daonew = strchr (daonew, newline)) != NULL) {
	    daonew = daonew + 1;
	    nlines = nlines + 1;
	    }
	}

    (void) fclose (fcat);
    return (nlines);
}


/* DAOLINE -- Get DAOFIND catalog entry for one star; return 0 if successful */

static char *
daoline (iline, line)

int iline;	/* Star sequence number in DAOFIND catalog */
char *line;	/* Pointer to iline'th entry (returned updated) */
{
    char *nextline;
    int i;

    if (iline > nlines) {
	fprintf (stderr, "DAOSTAR:  %d is not in catalog\n",iline);
	return (NULL);
	}
    else if (iline < 1 && line) {
	nextline = strchr (line, newline) + 1;
	}
    else {
	nextline = daobuff;
	for (i = 1; i < iline; i++) {
	    nextline = strchr (nextline, newline) + 1;
	    }
	}

    return (nextline);
}

/* Dec 11 1996	New subroutines
 *
 * Mar 20 1997	Removed unused variables, fixed logging after lint
 */
