/*** File libwcs/fileutil.c
 *** October 21, 1999
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics

 * Module:      fileutil.c (ASCII file utilities)
 * Purpose:     Find out things about ASCII files
 * Subroutine:	getfilelines (filename)
 *		Return number of lines in an ASCII file
 * Subroutine:	getfilebuff (filename)
 *		Return entire file contents in a character string
 * Subroutine:	getfilesize (filename)
 *		Return size of a binary or ASCII file
 * Subroutine:	isimlist (filename)
 *		Return 1 if file is list of FITS or IRAF image files, else 0

 * Copyright:   1999 Smithsonian Astrophysical Observatory
 *              You may do anything you like with this file except remove
 *              this copyright.  The Smithsonian Astrophysical Observatory
 *              makes no representations about the suitability of this
 *              software for any purpose.  It is provided "as is" without
 *              express or implied warranty.
 */

#include <stdlib.h>
#ifndef VMS
#include <unistd.h>
#endif
#include <stdio.h>
#include <fcntl.h>
#include <sys/file.h>
#include <errno.h>
#include <string.h>
#include "fitsfile.h"


/* GETFILELINES -- return number of lines in one file */

int
getfilelines (filename)

char    *filename;      /* Name of file for which to find number of lines */
{

    char *buffer, *bufline;
    int nlines = 0;
    char newline = 10;

    /* Read file */
    buffer = getfilebuff (filename);

    /* Count lines in file */
    if (buffer != NULL) {
	bufline = buffer;
	nlines = 0;
	while ((bufline = strchr (bufline, newline)) != NULL) {
            bufline = bufline + 1;
            nlines++;
	    }
	free (buffer);
	return (nlines);
	}
    else {
	return (0);
	}
}


/* GETFILEBUFF -- return entire file contents in one character string */

char *
getfilebuff (filename)

char    *filename;      /* Name of file for which to find number of lines */
{

    FILE *diskfile;
    int lfile, nr;
    char *buffer;

    /* Open file */
    if ((diskfile = fopen (filename, "r")) == NULL)
        return (NULL);

   /* Find length of file */
    if (fseek (diskfile, 0, 2) == 0)
        lfile = ftell (diskfile);
    else
        lfile = 0;
    if (lfile < 1) {
	fprintf (stderr,"GETFILEBUFF: File %s is empty\n", filename);
	fclose (diskfile);
	return (NULL);
	}

    /* Allocate buffer to hold entire file and read it */
    if ((buffer = calloc (1, lfile+1)) != NULL) {
 	fseek (diskfile, 0, 0);
        nr = fread (buffer, 1, lfile, diskfile);
	if (nr < lfile) {
	    fprintf (stderr,"GETFILEBUFF: File %s: read %d / %d bytes\n",
		     filename, nr, lfile);
	    free (buffer);
	    fclose (diskfile);
	    return (NULL);
	    }
	buffer[lfile] = (char) 0;
	fclose (diskfile);
	return (buffer);
	}
    else {
	fprintf (stderr,"GETFILEBUFF: File %s: no room for %d-byte buffer\n",
		 filename, lfile);
	fclose (diskfile);
	return (NULL);
	}
}


/* GETFILESIZE -- return size of one file in bytes */

int
getfilesize (filename)

char    *filename;      /* Name of file for which to find size */
{
    FILE *diskfile;
    long filesize;

    /* Open file */
    if ((diskfile = fopen (filename, "r")) == NULL)
        return (-1);

    /* Move to end of the file */
    if (fseek (diskfile, 0, 2) == 0)

        /* Position is the size of the file */
        filesize = ftell (diskfile);

    else
        filesize = -1;

    fclose (diskfile);

    return ((int) filesize);
}


/* ISIMLIST -- Return 1 if list of FITS or IRAF files, else 0 */
int
isimlist (filename)

char    *filename;      /* Name of file for which to find size */
{
    FILE *diskfile;
    char token[256];
    int ncmax = 254;

    if ((diskfile = fopen (filename, "r")) == NULL)
	return (0);
    else {
	first_token (diskfile, ncmax, token);
	fclose (diskfile);
	if (isfits (token) | isiraf (token))
	    return (1);
	else
	    return (0);
	}
}

static char *token1;


/* FIRST_TOKEN -- Return first token from the next line of an ASCII file */

int
first_token (diskfile, ncmax, token)

FILE	*diskfile;		/* File descriptor for ASCII file */
int	ncmax;			/* Maximum number of characters returned */
char	*token;			/* First token on next line (returned) */
{
    char *lastchar, *lspace;

    /* If line can be read, add null at the end of the first token */
    if (fgets (token, ncmax, diskfile) != NULL) {
	lastchar = token + strlen (token) - 1;

	/* Remove trailing spaces or control characters */
	while (*lastchar <= 32)
	    *lastchar-- = 0;

	if ((lspace = strchr (token, ' ')) != NULL) {
	    *lspace = (char) 0;
	    token1 = lspace + 1;
	    }
	else
	    token1 = NULL;
	return (1);
	}
    else
	return (0);
}

char *
next_token ()

{
    return (token1);
}


/*
 * Jul 14 1999	New subroutines
 * Jul 15 1999	Add getfilebuff()
 * Oct 15 1999	Fix format eror in error message
 * Oct 21 1999	Fix declarations after lint
 * Dec  9 1999	Add next_token(); set pointer to next token in first_token
 */
