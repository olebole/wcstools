/*** File fitslib/hget.c
 *** February 20, 1996
 *** By Doug Mink, Harvard-Smithsonian Center for Astrophysics

 * Module:	hget.c (Get FITS Header parameter values)
 * Purpose:	Extract values for variables from FITS header string
 * Subroutine:	hgeti2 (hstring,keyword) returns short integer
 * Subroutine:	hgeti4 (hstring,keyword) returns long integer
 * Subroutine:	hgetr4 (hstring,keyword) returns real
 * Subroutine:	hgetra (hstring,keyword) returns double RA in degrees
 * Subroutine:	hgetdec (hstring,keyword) returns double Dec in degrees
 * Subroutine:	hgetr8 (hstring,keyword) returns double
 * Subroutine:	hgetl  (hstring,keyword) returns logical int (0=F, 1=T)
 * Subroutine:	hgetc  (hstring,keyword) returns character string (drops quotes)
 * Subroutine:	hgets  (hstring,keyword, str) returns character string (drops quotes)
 * Subroutine:	ksearch (hstring,keyword) returns pointer to header string entry

 * Copyright:   1995, 1996 Smithsonian Astrophysical Observatory
 *              You may do anything you like with this file except remove
 *              this copyright.  The Smithsonian Astrophysical Observatory
 *              makes no representations about the suitability of this
 *              software for any purpose.  It is provided "as is" without
 *              express or implied warranty.

*/

#include <string.h>		/* NULL, strlen, strstr, strcpy */
#include <stdio.h>
#include "fitshead.h"	/* FITS header extraction subroutines */
#include <stdlib.h>

char *hgetc ();

char val[30];

/* Extract long value for variable from FITS header string */

int
hgeti4 (hstring,keyword,ival)

char *hstring;	/* character string containing FITS header information
		   in the format <keyword>= <value> {/ <comment>} */
char *keyword;	/* character string containing the name of the variable
		   to be returned.  hget searches for a line beginning
		   with this string.  if "[n]" is present, the n'th
		   token in the value is returned.
		   (the first 8 characters must be unique) */
int *ival;
{
char *value;

/* Get value and comment from header string */
	value = hgetc (hstring,keyword);

/* Translate value from ASCII to binary */
	if (value != NULL) {
	    strcpy (val, value);
	    *ival = (int) atol (val);
	    return (1);
	    }
	else {
	    return (0);
	    }
}


/* Extract integer*2 value for variable from fits header string */

int
hgeti2 (hstring,keyword,ival)

char *hstring;	/* character string containing FITS header information
		   in the format <keyword>= <value> {/ <comment>} */
char *keyword;	/* character string containing the name of the variable
		   to be returned.  hget searches for a line beginning
		   with this string.  if "[n]" is present, the n'th
		   token in the value is returned.
		   (the first 8 characters must be unique) */
short *ival;
{
char *value;

/* Get value and comment from header string */
	value = hgetc (hstring,keyword);

/* Translate value from ASCII to binary */
	if (value != NULL) {
	    strcpy (val, value);
	    *ival = (short) atoi (val);
	    return (1);
	    }
	else {
	    return (0);
	    }
}

/* Extract real value for variable from FITS header string */

int
hgetr4 (hstring,keyword,rval)

char *hstring;	/* character string containing FITS header information
		   in the format <keyword>= <value> {/ <comment>} */
char *keyword;	/* character string containing the name of the variable
		   to be returned.  hget searches for a line beginning
		   with this string.  if "[n]" is present, the n'th
		   token in the value is returned.
		   (the first 8 characters must be unique) */
float *rval;
{
	char *value;

/* Get value and comment from header string */
	value = hgetc (hstring,keyword);

/* translate value from ASCII to binary */
	if (value != NULL) {
	    strcpy (val, value);
	    *rval = (float) atof (val);
	    return (1);
	    }
	else {
	    return (0);
	    }
}


/* Extract real*8 right ascension in degrees from FITS header string */

int
hgetra (hstring,keyword,dval)

char *hstring;	/* character string containing FITS header information
		   in the format <keyword>= <value> {/ <comment>} */
char *keyword;	/* character string containing the name of the variable
		   to be returned.  hget searches for a line beginning
		   with this string.  if "[n]" is present, the n'th
		   token in the value is returned.
		   (the first 8 characters must be unique) */
double *dval;	/* Right ascension in degrees (returned) */
{
	int hgetdec();

	if (hgetdec (hstring,keyword,dval)) {
	    *dval = *dval * 15.0;
	    return (1);
	    }
	else
	    return (0);
}


/* Extract real*8 declination in degrees from FITS header string */

int
hgetdec (hstring,keyword,dval)

char *hstring;	/* character string containing FITS header information
		   in the format <keyword>= <value> {/ <comment>} */
char *keyword;	/* character string containing the name of the variable
		   to be returned.  hget searches for a line beginning
		   with this string.  if "[n]" is present, the n'th
		   token in the value is returned.
		   (the first 8 characters must be unique) */
double *dval;	/* Right ascension in degrees (returned) */
{
	double hval, rval, deg, min, sec, sign;
	char *value,val[30], *c1;
	char *strsrch();

/* Get value from header string */
	value = hgetc (hstring,keyword);

/* Translate value from ASCII colon-delimited string to binary */
	if (value != NULL) {
	    str2dec (value, dval);
	    return (1);
	    }
	else
	    return (0);
}


/* Extract real*8 value for variable from FITS header string */

int
hgetr8 (hstring,keyword,dval)

char *hstring;	/* character string containing FITS header information
		   in the format <keyword>= <value> {/ <comment>} */
char *keyword;	/* character string containing the name of the variable
		   to be returned.  hget searches for a line beginning
		   with this string.  if "[n]" is present, the n'th
		   token in the value is returned.
		   (the first 8 characters must be unique) */
double *dval;
{
	char *value,val[30];

/* Get value and comment from header string */
	value = hgetc (hstring,keyword);

/* Translate value from ASCII to binary */
	if (value != NULL) {
	    strcpy (val, value);
	    *dval = atof (val);
	    return (1);
	    }
	else {
	    return (0);
	    }
}


/* Extract logical value for variable from FITS header string */

int
hgetl (hstring,keyword,ival)

char *hstring;	/* character string containing FITS header information
		   in the format <keyword>= <value> {/ <comment>} */
char *keyword;	/* character string containing the name of the variable
		   to be returned.  hget searches for a line beginning
		   with this string.  if "[n]" is present, the n'th
		   token in the value is returned.
		   (the first 8 characters must be unique) */
int *ival;
{
	char *value;

/* Get value and comment from header string */
	value = hgetc (hstring,keyword);

/* Translate value from ASCII to binary */
	if (value != NULL) {
	    strcpy (val, value);
            value = &val[0];
	    if (value[0] == 't' || value[0] == 'T')
		*ival = 1;
	    else
		*ival = 0;
	    return (1);
	    }
	else {
	    return (0);
	    }
}


/* Extract string value for variable from FITS header string */

int
hgets (hstring, keyword, lcval, cval)

char *hstring;	/* character string containing FITS header information
		   in the format <keyword>= <value> {/ <comment>} */
char *keyword;	/* character string containing the name of the variable
		   to be returned.  hget searches for a line beginning
		   with this string.  if "[n]" is present, the n'th
		   token in the value is returned.
		   (the first 8 characters must be unique) */
int lcval;	/* Maximum length of cval in characters */
char *cval;
{
	char *value;
	int lval;

/* Get value and comment from header string */
	value = hgetc (hstring,keyword);

	if (value != NULL) {
	    lval = strlen (value);
	    if (lval < lcval)
		strcpy (cval, value);
	    else
		strncpy (cval, value, lval-1);
	    return (1);
	    }
	else
	    return (0);
}


/* Extract character value for variable from FITS header string */

char *
hgetc (hstring,keyword0)

char *hstring;
char *keyword0;
{
	static char cval[80];
	char *value;
	char cwhite[2];
	char squot[2],dquot[2],lbracket[2],rbracket[2],slash[2];
	char keyword[16];
	char line[100];
	char *vpos,*cpar;
	char *q1, *q2, *v1, *v2, *c1, *brack1, *brack2;
	char *strsrch();
	int ipar, i;

	squot[0] = 39;
	squot[1] = 0;
	dquot[0] = 34;
	dquot[1] = 0;
	lbracket[0] = 91;
	lbracket[1] = 0;
	rbracket[0] = 93;
	rbracket[1] = 0;
	slash[0] = 47;
	slash[1] = 0;

/* Find length of variable name */
	strcpy (keyword,keyword0);
	brack1 = strsrch (keyword,lbracket);
	if (brack1 != NULL) *brack1 = NULL;

/* Search header string for variable name */
	vpos = ksearch (hstring,keyword);

/* Exit if not found */
	if (vpos == NULL) {
	    return (NULL);
	    }

/* Initialize line to nulls */
	 for (i = 0; i < 100; i++)
	    line[i] = 0;

/* In standard FITS, data lasts until 80th character */

/* Extract entry for this variable from the header */
	strncpy (line,vpos,80);

/* check for quoted value */
	q1 = strsrch (line,squot);
	c1 = strsrch (line,slash);
	if (q1 != NULL) {
	    if (c1 != NULL && q1 < c1)
		q2 = strsrch (q1+1,squot);
	    else if (c1 == NULL)
		q2 = strsrch (q1+1,squot);
	    else
		q1 = NULL;
	    }
	else {
	    q1 = strsrch (line,dquot);
	    if (q1 != NULL) {
		if (c1 != NULL && q1 < c1)
		    q2 = strsrch (q1+1,dquot);
		else if (c1 == NULL)
		    q2 = strsrch (q1+1,dquot);
		else
		    q1 = NULL;
		}
	    else {
		q1 = NULL;
		q2 = line + 10;
		}
	    }

/* Extract value and remove excess spaces */
	if (q1 != NULL) {
	    v1 = q1 + 1;;
	    v2 = q2;
	    c1 = strsrch (q2,"/");
	    }
	else {
	    v1 = strsrch (line,"=") + 1;
	    c1 = strsrch (line,"/");
	    if (c1 != NULL)
		v2 = c1;
	    else
		v2 = line + 79;
	    }

/* Ignore leading spaces */
	while (*v1 == ' ' && v1 < v2) {
	    v1++;
	    }

/* Drop trailing spaces */
	*v2 = 0;
	v2--;
	while (*v2 == ' ' && v2 > v1) {
	    *v2 = 0;
	    v2--;
	    }

	strcpy (cval,v1);
	value = cval;

/* If keyword has brackets, extract appropriate token from value */
	if (brack1 != NULL) {
	    brack2 = strsrch (keyword,rbracket);
	    if (brack2 != NULL) {
		*brack2 = NULL;
		ipar = atoi (brack1);
		if (ipar > 0) {
		    cwhite[0] = 32;
		    cwhite[1] = NULL;
		    for (i = 1; i <= ipar; i++) {
			cpar = strtok (v1,cwhite);
			}
		    if (cpar != NULL) {
			strcpy (cval,cpar);
			}
		    else
			value = NULL;
		    }
		}
	    }

	return (value);
}


/* Find FITS header line containing specified keyword */

char *ksearch (hstring,keyword)

/* Find entry for keyword keyword in FITS header string hstring.
   (the keyword may have a maximum of eight letters)
   NULL is returned if the keyword is not found */

char *hstring;	/* character string containing fits-style header
		information in the format <keyword>= <value> {/ <comment>}
		the default is that each entry is 80 characters long;
		however, lines may be of arbitrary length terminated by
		nulls, carriage returns or linefeeds, if packed is true.  */
char *keyword;	/* character string containing the name of the variable
		to be returned.  ksearch searches for a line beginning
		with this string.  The string may be a character
		literal or a character variable terminated by a null
		or '$'.  it is truncated to 8 characters. */
{
char *loc, *headnext, *headlast, *pval;
int icol, icol0, nline, prevchar, nextchar, lkey;
char *strsrch();

	pval = 0;

/* Search header string for variable name */
	headlast = hstring + strlen (hstring);
	headnext = hstring;
	pval = NULL;
	while (headnext < headlast) {
	    loc = strsrch (headnext,keyword);

	/* Exit if keyword is not found */
	    if (loc == NULL) {
		break;
		}

	/* If parameter name in string entry is longer, keep searching */
	    lkey = strlen (keyword);
	    nextchar = (int) *(loc + lkey);
	    if (nextchar != 61 && nextchar > 32 && nextchar < 127) {
		headnext = loc + 1;
		}

	/* See if this is a legal position for a keyword.  if not, assume
	   it is in the comment or value  field, and keep searching. */
	    else {
		icol = loc - hstring + 1;
		nline = icol / 80;
		icol0 = (80 * nline) + 1;
		icol = icol - icol0;
		if (icol > 8) {
		    headnext = loc + 1;
		    }
		else if (icol > 1) {
		    prevchar = (int) *(loc - 1);
		    if (prevchar > 32 && prevchar < 127) {
			headnext = loc + 1;
			}
		    }
		else {
		    pval = loc;
		    break;
		    }
		}
	    }

/* Return pointer to calling program */
	return (pval);
}


/* Find string s2 within string s1 */

char *strsrch (s1, s2)

char *s1;	/* String to search */
char *s2;	/* String to look for */

{
    char *s,*s1e;
    char cfirst,clast;
    int i,ls1,ls2;

    /* Return null string if either pointer is NULL */
    if (s1 == NULL || s2 == NULL)
	return (NULL);

    /* A zero-length pattern is found in any string */
    ls2 = strlen (s2);
    if (ls2 ==0)
	return (s1);

    /* Only a zero-length string can be found in a zero-length string */
    ls1 = strlen (s1);
    if (ls1 ==0)
	return (NULL);

    cfirst = s2[0];
    clast = s2[ls2-1];
    s1e = s1 + ls1 - ls2 + 1;
    s = s1;
    while (s < s1e) { 

	/* Search for first character in pattern string */
	if (*s == cfirst) {

	    /* If single character search, return */
	    if (ls2 == 1)
		return (s);

	    /* Search for last character in pattern string if first found */
	    if (s[ls2-1] == clast) {

		/* If two-character search, return */
		if (ls2 == 2)
		    return (s);

		/* If 3 or more characters, check for rest of search string */
		i = 1;
		while (i++ < ls2 && s[i] == s2[i])
		    ;

		/* If entire string matches, return */
		if (i >= ls2)
		    return (s);
		}
	    }
	s++;
	}
    return (NULL);
}


/* Read the right ascension, ra, in sexagesimal hours from in[] */

void
str2ra (in, ra)

char	*in;	/* Character string */
double	*ra;	/* Right ascension in degrees (returned) */

{
    str2dec (in, ra);
    *ra = *ra * 15.0;

    return;
}


/* Read the declination, dec, in sexagesimal degrees from in[] */

void
str2dec (in, dec)

char	*in;	/* Character string */
double	*dec;	/* Declination in degrees (returned) */

{
    double deg, min, sec, sign;
    char *value, temp[4], *c1;
    char *strsrch();

    *dec = 0.0;

    /* Translate value from ASCII colon-delimited string to binary */
    if (in[0]) {
	value = in;
	if (strsrch (value,"-") == NULL)
	    sign = 1.0;
	else
	    sign = -1.0;
	if ((c1 = strsrch (value,":")) != NULL) {
	    *c1 = 0;
	    deg = (double) atoi (value);
	    *c1 = ':';
	    value = c1 + 1;
	    if ((c1 = strsrch (value,":")) != NULL) {
		*c1 = 0;
		min = (double) atoi (value);
		*c1 = ':';
		value = c1 + 1;
		sec = atof (value);
		}
	    else if ((c1 = strsrch (value,".")) != NULL) {
		min = atof (value);
		sec = 0.0;
		}
	    *dec = sign * (deg + (min / 60.0) + (sec / 3600.0));
	    }
	else if ((c1 = strsrch (value,".")) != NULL) {
	    *dec = atof (value);
	    }

	}
    return;
}
/* Oct 28 1994	New program
 *
 * Mar  1 1995	Search for / after second quote, not first one
 * May  2 1995	Initialize line in HGETC; deal with logicals in HGETL better
 * May  4 1995	Declare STRSRCH in KSEARCH
 * Aug  7 1995  Fix line initialization in HGETC
 * Dec 22 1995	Add HGETRA and HGETDEC to get degrees from xx:xx:xx.xxx string

 * Jan 26 1996	Fix HGETL to not crash when parameter is not present
 * Feb  1 1996	Fix HGETC to deal with quotes correctly
 * Feb  1 1996	Fix HGETDEG to deal with sign correctly
 * Feb  6 1996	Add HGETS to update character strings
 * Feb  8 1996	Fix STRSRCH to find final characters in string
 * Feb 23 1996	Add string to degree conversions
 */
