/*** File iolib/parstring.c
 *** September 23, 1993
 *** By Doug Mink

 Subroutines to read, write and access ASCII packed FITS strings:

	pcread (pathname,parms)		Read parameter string from file
	pcwrite (pathname,parms)	Write parameter string to file
	pcgstr (parms,keyword,string)	Extract character string
	pcpstr (parms, keyword, string)	Implant character string
	pcgeti (parms, keyword, value)	Extract integer
	pcputi (parms, keyword, value)	Implant integer
	pcgetd (parms, keyword, value)	Extract double float
	pcputd (parms, keyword, value)	Implant double float
	pcgetb (parms, keyword, value)	Extract logical (y or n or t or f)
	pcputb (parms, keyword, value)	Implant logical (y or n or t or f)
	pcadd (parms, newvalue)		Add string to parameter string
	pcfind (parms, keyword)		Find keyword in parameter string
	char *pcsrch (a, b)		Find substring in string
*/

#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/file.h>
#include <fcntl.h>

int maxlparms = 2000;
char newline = 10;
char temp[80];		/* value of found parameter */
char *vstart;		/* pointer to start of found parameter */
char *vend;		/* pointer to end+1 of found parameter */


/* Read FITS-like parameter file, returning pointer to parameter
 * string or NULL if unsuccessful
 */

pcread (pathname,parms)

    char *pathname;
    char *parms;
{
    int fd, lparms;
    char *pcsrch();
    char *tparms, *last;
    int newfile = 0;
    char *malloc();

    /* Open file to read parameter string */
    if ((fd = open (pathname,O_RDONLY)) < 3) {
	fprintf (stderr,"Cannot open parameter file %s\n",pathname);
	newfile = 1;
	}

    /* Allocate maximum memory to read parameter string */
    else if ((tparms = malloc (maxlparms)) == NULL) {
	fprintf (stderr,"Can't malloc %d bytes for param string\n",maxlparms);
	exit (1);
	}

    /* Read parameter string from file */
    if ((lparms = read (fd,tparms,maxlparms)) < 4) {
	fprintf (stderr,"Error reading IRAF parameter file %s\n",pathname);
	newfile = 1;
	}

    /* Search for END<nl> in string */
    if ((last = pcsrch (tparms,"END\n")) == NULL) {
	fprintf (stderr,"No END found in parameter string\n");
	exit (1);
	}

    /* If no parameter string has been read, initialize a new one */
    if (newfile) {
	last = tparms;
	strcpy (tparms,"END\n");
	}

    /* Allocate new parameter string */
    lparms = last + 4 - tparms;
    if ((parms = malloc (lparms+1)) == NULL) {
	fprintf (stderr,"Can't malloc %d bytes for param string\n",lparms+1);
	exit (1);
	}

    strncpy (parms,tparms,lparms);
    free (tparms);

    close (fd);
    return;
}


/* Write parameter file */

pcwrite (pathname,parms)

    char *pathname;
    char *parms;

{
    int fd, nw, lparms;

    lparms = strlen (parms);

    if ((fd = open (pathname,O_WRONLY)) < 3) {
	fprintf (stderr,"Error opening IRAF parameter file %s\n",pathname);
	exit(1);
	}

    if ((nw = write (fd,parms,lparms)) < lparms) {
	fprintf (stderr,"Error writing IRAF parameter file %s\n",pathname);
	exit(1);
	}

    close (fd);
    return(0);
}


/* Extract character string from parameter string */

pcgstr (parms,keyword,string)

    char *parms;
    char *keyword;	/* Name of parameter to return */
    char *string;	/* Address of string returned */
{
    int slen;
    char *malloc();
    int pcfind ();

    /* Search for parameter name in parameter string */
    if (pcfind (parms, keyword) == 1) {
	slen = strlen (temp);

	if ((string = malloc (slen+1)) == NULL) {
	    fprintf (stderr,"Can't malloc %d bytes for string\n",slen+1);
	    exit(1);
	    }
	strncpy (string,temp,slen);
	string[slen] = 0;
	return (0);
	}

    return (1);
}


/* Implant character string into IRAF parameter file */

pcpstr (parms, keyword, string)

    char *parms;
    char *keyword;	/* Name of parameter to set */
    char *string;	/* Value of parameter to set */
{
    int pcfind (), atoi();

    /* Search for parameter name in parameter string */
    if (pcfind (parms, keyword) == 1) {
	pcadd (parms, string);
	return (0);
	}

    return (1);
}


/* Extract integer from parameter string */

pcgeti (parms, keyword, value)

    char *parms;
    char *keyword;	/* Name of parameter to return */
    int *value;		/* Value of parameter (returned) */
{
    int pcfind (), atoi();

    /* Search for parameter name in parameter string */
    if (pcfind (parms, keyword) == 1) {
	*value = atoi (temp);
	return (0);
	}

    return (1);
}


/* Implant integer into parameter string */

pcputi (parms, keyword, value)

    char *parms;
    char *keyword;	/* Name of parameter to set */
    int value;		/* Value of parameter to set */
{
    int pcfind ();
    char cvalue[20];

    /* Search for parameter name in parameter string */
    if (pcfind (parms, keyword) == 1) {
	sprintf (cvalue,"%d",value);
	pcadd (parms, cvalue);
	return (0);
	}

    return (1);
}


/* Extract double precision floating point value from parameter string */

pcgetd (parms, keyword, value)

    char *parms;
    char *keyword;	/* Name of parameter to return */
    double *value;
{
    int pcfind ();

    /* Search for parameter name in parameter string */
    if (pcfind (parms, keyword) == 1) {
	*value = atod (temp);
	return (0);
	}

    return (1);
}


/* Implant double precision floating point value into IRAF parameter file */

pcputd (parms, keyword, value)

    char *parms;
    char *keyword;	/* Name of parameter to set */
    double value;	/* Value of parameter to set */
{
    int pcfind ();
    char cvalue[20];

    /* Search for parameter name in parameter string */
    if (pcfind (parms, keyword) == 1) {
	bzero (cvalue,20);
	sprintf (cvalue,"%f",value);
	pcadd (parms, cvalue);
	return (0);
	}

    return (1);
}


/* Extract boolean value from IRAF parameter file */

pcgetb (parms, keyword, value)

    char *parms;
    char *keyword;	/* Name of parameter to return */
    int *value;
{
    int pcfind ();

    /* Search for parameter name in parameter string */
    if (pcfind (parms, keyword) == 1) {
	if (temp[0] == 'y')
	    *value = 1;
	else if (temp[0] == 't')
	    *value = 1;
	else if (temp[0] == 'f')
	    *value = 0;
	else
	    *value = 0;
	return (0);
	}

    return (1);
}


/* Implant boolean value into IRAF parameter file */

pcputb (parms, keyword, value)

    char *parms;
    char *keyword;	/* Name of parameter to set */
    int value;		/* Value of parameter to set */
{
    int pcfind (), atoi();

    /* Search for parameter name in parameter string */
    if (pcfind (parms, keyword) == 1) {
	if (value == 0)
	    pcadd (parms, "n");
	else
	    pcadd (parms, "y");
	return (0);
	}

    return (1);
}


pcadd (parms, newvalue)

    char *parms;
    char *newvalue;

{
    int lnew, lold;
    char *pold, *pnew;
    char *realloc();
    int lparms;

    lnew = strlen (newvalue);
    lold = vend - vstart;
    lparms = strlen (parms);

    /* If replacement shorter than old, add it and contract string */
    if (lnew < lold) {

	strncpy (vstart,newvalue,lnew);
	pold = vend;
	pnew = vstart + lnew;
	while (pold < parms + lparms) *pnew++ = *pold++;

	/* reallocate parameter string */
	lparms = lparms + lnew - lold;
	if ((parms = realloc (parms, lparms+1)) == NULL) {
	    fprintf (stderr,"Can't realloc %d bytes for param string\n",lparms+1);
	    exit (1);
	    }
	parms[lparms] = 0;
	}

    /* If replacement longer than old, lengthen string and add it */
    else {

	/* reallocate parameter string */
	lparms = lparms + lnew - lold;
	if ((parms = realloc (parms, lparms+1)) == NULL) {
	    fprintf (stderr,"Can't realloc %d bytes for param string\n",lparms+1);
	    exit (1);
	    }

	pold = parms + lparms;
	pnew = parms + lparms + lnew - lold;
	while (pold > vend) *pnew-- = *pold--;
	strncpy (vstart,newvalue,lnew);
	parms[lparms] = 0;
	}

    return;
}


/* Search for specific parameter in parameter string */

pcfind (parms, keyword)

    char *parms;
    char *keyword;	/* Name of parameter to return */
{
    int n, lpar;
    int strncmp();
    char *endline, *line, *lastchar;
    char *psearch, *nextchar;
    char *pcsrch();
    char quote1 = 39;
    char quote2 = 34;

    line = temp;
    psearch = parms;
    bzero (line,120);

    /* Search for first character of parameter name */
    while ((line = pcsrch (psearch,keyword)) != NULL) {

	/* Check to see if it's at the beginning of a line */
	if (line > parms) lastchar = line - 1;
	else lastchar = parms;
	if ((line == parms) || (*lastchar == newline)) {

	    /* If parameter in file matches parameter name, find value */
	    lpar = strlen (keyword);
	    nextchar = line + lpar;
	    if (*nextchar == ' ' || *nextchar == '=') {

		/* Find end of line */
		endline = strchr (line,newline);

		/* Look for single quotes around value first */
		vstart = strchr (line,quote1);
		if (vstart != NULL && vstart < endline) {
		    vstart = vstart + 1;
		    vend = strchr (vstart,quote1);
		    if (vend == NULL) {
			fprintf (stderr,"unmatched quotes in parameter string\n");
			vend = endline - 1;
			}
		    else
			vend = vend - 1;
		    }

		/* Then look for double quotes around value */
		if (vstart == NULL) {
		    vstart = strchr (line,quote2);
		    if (vstart != NULL && vstart < endline) {
			vstart = strchr (line,quote2) + 1;
			vend = strchr (vstart,quote2);
			if (vend == NULL) {
			    fprintf (stderr,"unmatched quotes in parameter string\n");
			    vend = endline - 1;
			    }
			else
			    vend = vend - 1;
			}
		    }

		/* If no quotes, value starts at = sign */
		if (vstart == NULL) {
		    vstart = strchr (line,'=') + 1;
		    while (*vstart == ' ' || *vstart == 9) vstart++;
		    }

		/* If no quotes, value ends at "/", or end of line */
		if (vend == NULL) {
		    vend = strchr (vstart,'/');
		    if (vend == NULL)
			vend = endline;
		    else
			vend = vend - 1;
		    while (*vend == ' ') vend--;
		    }

		/* Copy from vstart to vend into temp */
		n = vend - vstart + 1;
		if (n > 120) n = 120;
		strncpy (temp,vstart,n);
		return (1);
		}
	    }

	/* Otherwise start searching from this point */
	psearch = line + 1;
	}

    /* Return 0 if parameter name not found */
    return (0);
}


/* Search for string within string (borrowed from Unix index(3F)
 * Return pointer to start of substring; if not found, return NULL */

char
*pcsrch (a, b)

    char *a, *b;
{
    int la, lb;
    char *a1, *aend, *s, *t, *bend;

    la = strlen (a);
    aend = a + la;
    lb = strlen (b);
    bend = b + lb;

    for (a1 = a; a1 < aend; ++a1)
	{
	s = a1;
	t = b;
	while (t < bend)
		if (*s++ != *t++)
			goto no;
	return (a1);
	no: ;
	}
return (NULL);
}
/*
   Jul 11 1991  Wrote program
   May  4 1993	Call strchr instead of index for System 5 compatibility
   Sep 23 1993	Declare undeclared subroutine
*/
