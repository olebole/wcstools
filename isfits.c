/* File isfits.c
 * August 12, 2013
 * By Jessica Mink, Harvard-Smithsonian Center for Astrophysics
 * Send bug reports to jmink@cfa.harvard.edu

   Copyright (C) 2008-2013
   Smithsonian Astrophysical Observatory, Cambridge, MA USA

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 *
 * Return 1 if argument is an integer, 2 if it is floating point, else 0
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libwcs/fitsfile.h"


int
main (ac, av)
int ac;
char **av;
{
    char *str;

    /* Check for version or help command first */
    str = *(av+1);
    if (!str || !strcmp (str, "help") || !strcmp (str, "-help")) {
	fprintf (stderr,"Usage: Return 1 if argument is a FITS file, else 0\n");
	exit (1);
	}
    else if (!strcmp (str, "version") || !strcmp (str, "-version")) {
	exit (1);
	}

    /* check to see if this is a FITS file */
    else
	printf ("%d\n", isfits (str));

    exit (0);
}
/* Apr 25 2008	New program
 *
 * Oct 29 2010	Include fitsfile.h instead of fitshead.h
 *
 * Dec 14 2011	Fix comments
 *
 * Aug 12 2013	Add linefeed to usage
 */
