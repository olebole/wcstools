/*** File webread.c
 *** January 3, 2001
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 *** (http code from John Roll)
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "wcs.h"
#include "wcscat.h"

#define CHUNK   8192
#define LINE    1024
#define MAXHOSTNAMELENGTH	256

#ifdef HAVE_FCNTL_H
#include <fcntl.h>
#else
#include <sys/fcntl.h>
#endif

#include <sys/time.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

static int FileINetParse(char *file, int port, struct sockaddr_in *adrinet);

static FILE *SokOpen();
#define XFREAD  1
#define XFWRITE 2
#define XFCREAT 4

#define File    FILE *
#define FileFd(fd)              fileno(fd)
static char newline = 10;
static struct TabTable *webopen();


/* WEBREAD -- Read a catalog over the web and return results */

int
webread (caturl,refcatname,distsort,cra,cdec,dra,ddec,drad,sysout,
                 eqout,epout,mag1,mag2,nstarmax,
		 unum,ura,udec,upra,updec,umag,umagb,utype,nlog)

char	*caturl;	/* URL of search engine */
char	*refcatname;	/* Name of catalog (UAC, USAC, UAC2, USAC2) */
int	distsort;	/* 1 to sort stars by distance from center */
double	cra;		/* Search center J2000 right ascension in degrees */
double	cdec;		/* Search center J2000 declination in degrees */
double	dra;		/* Search half width in right ascension in degrees */
double	ddec;		/* Search half-width in declination in degrees */
double	drad;		/* Limiting separation in degrees (ignore if 0) */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
double	mag1,mag2;	/* Limiting magnitudes (none if equal) */
int	nstarmax;	/* Maximum number of stars to be returned */
double	*unum;		/* Array of UA numbers (returned) */
double	*ura;		/* Array of right ascensions (returned) */
double	*udec;		/* Array of declinations (returned) */
double	*upra;		/* Array of right ascension proper motions (returned) */
double	*updec;		/* Array of declination proper motions (returned) */
double	*umag;		/* Array of red magnitudes (returned) */
double	*umagb;		/* Array of blue magnitudes (returned) */
int	*utype;		/* Array of plate numbers (returned) */
int	nlog;		/* Logging interval */
{
    char srchurl[LINE];
    char temp[64];
    struct TabTable *tabtable;
    double dtemp;
    int lurl;
    struct StarCat *starcat;
    char cstr[32];
    double ra, dec;

    /* Convert coordinate system to string */
    wcscstr (cstr, sysout, eqout, epout);

    /* Set up search query from arguments */
    lurl = strlen (caturl);

    /* Set up query for scat used as server */
    if (!strncmp (caturl+lurl-4,"scat",4)) {
	sprintf (srchurl, "?catalog=%s&ra=%.7f&dec=%.7f&system=%s&",
		 refcatname, cra, cdec, cstr);
	if (drad != 0.0) {
	    dtemp = drad * 3600.0;
	    sprintf (temp, "radius=%.3f&",dtemp);
	    strcat (srchurl, temp);
	    }
	else {
	    dtemp = dra * 3600.0;
	    sprintf (temp, "dra=%.3f&",dtemp);
	    strcat (srchurl, temp);
	    dtemp = ddec * 3600.0;
	    sprintf (temp, "ddec=%.3f&",dtemp);
	    strcat (srchurl, temp);
	    }
	if (mag1 != mag2) {
	    sprintf (temp, "mag1=%.2f&mag=%.2f&",mag1,mag2);
	    strcat (srchurl, temp);
	    }
	if (epout != 0.0) {
	    sprintf (temp, "epoch=%.5f&", epout);
	    strcat (srchurl, temp);
	    }
	if (nlog > 0)
	    fprintf (stderr, "%s%s\n", caturl, srchurl);
	}

    /* Set up query for ESO GSC server */
    if (!strncmp (caturl+lurl-10,"gsc-server",10)) {
	ra = cra;
	dec = cdec;
	if (sysout != WCS_J2000)
	    wcscon (sysout, WCS_J2000, eqout, 2000.0, &ra, &dec, epout);
	sprintf (srchurl, "?ra=%.7f&dec=%.7f&", ra, dec);
	if (drad != 0.0) {
	    dtemp = drad * 60.0;
	    sprintf (temp, "radius=%.3f&",dtemp);
	    }
	else {
	    dtemp = dra / cos (degrad(dec));
	    if (dtemp > ddec)
		sprintf (temp, "radius=%.3f&",dtemp*60.0);
	    else
		sprintf (temp, "radius=%.3f&",ddec*60.0);
	    }
	strcat (srchurl, temp);
	sprintf (temp, "nout=%d&mime=skycat", nstarmax);
	strcat (srchurl, temp);
	if (nlog > 0)
	    fprintf (stderr, "%s%s\n", caturl, srchurl);
	}

    /* Set up query for ESO USNO A server */
    if (!strncmp (caturl+lurl-12,"usnoa-server",12)) {
	ra = cra;
	dec = cdec;
	if (sysout != WCS_J2000)
	    wcscon (sysout, WCS_J2000, eqout, 2000.0, &ra, &dec, epout);
	sprintf (srchurl, "?ra=%.7f&dec=%.7f&", ra, dec);
	if (drad != 0.0) {
	    dtemp = drad * 3600.0;
	    sprintf (temp, "radius=%.3f&",dtemp);
	    }
	else {
	    if (dra > ddec)
		sprintf (temp, "radius=%.3f&",dra*3600.0);
	    else
		sprintf (temp, "radius=%.3f&",ddec*3600.0);
	    }
	strcat (srchurl, temp);
	if (mag1 != mag2) {
	    sprintf (temp, "mag=%.2f,%.2f&", mag1, mag2);
	    strcat (srchurl, temp);
	    }
	sprintf (temp, "format=8&sort=mr");
	strcat (srchurl, temp);
	if (nlog > 0)
	    fprintf (stderr,"%s%s\n", caturl, srchurl);
	}

    /* Run search across the web */
    if ((tabtable = webopen (caturl, srchurl, nlog)) == NULL) {
	if (nlog > 0)
	    fprintf (stderr, "WEBREAD: %s failed\n", srchurl);
	return (0);
	}

    /* Return if no data */
    if (tabtable->tabdata == NULL || strlen (tabtable->tabdata) == 0) {
	if (nlog > 0)
	    fprintf (stderr, "WEBRNUM: No data returned\n");
	return (0);
	}

    /* Open returned Starbase table as a catalog */
    if ((starcat = tabcatopen (caturl, tabtable)) == NULL) {
	if (nlog > 0)
	    fprintf (stderr, "WEBREAD: Could not open Starbase table as catalog\n");
	return (0);
	}

    /* Extract desired sources from catalog  and return them */
    return (tabread (caturl,distsort,cra,cdec,dra,ddec,drad,
	     sysout,eqout,epout,mag1,mag2,nstarmax,&starcat,
	     unum,ura,udec,upra,updec,umag,umagb,utype,NULL,nlog));
}


int
webrnum (caturl,refcatname,nnum,sysout,eqout,epout,
	 unum,ura,udec,upra,updec,umag,umagb,utype,nlog)

char	*caturl;	/* URL of search engine */
char	*refcatname;	/* Name of catalog (UAC, USAC, UAC2, USAC2) */
int	nnum;		/* Number of stars to find */
int	sysout;		/* Search coordinate system */
double	eqout;		/* Search coordinate equinox */
double	epout;		/* Proper motion epoch (0.0 for no proper motion) */
double	*unum;		/* Array of UA numbers to find */
double	*ura;		/* Array of right ascensions (returned) */
double	*udec;		/* Array of declinations (returned) */
double	*upra;		/* Array of right ascensions proper motion (returned) */
double	*updec;		/* Array of declination proper motions (returned) */
double	*umag;		/* Array of red magnitudes (returned) */
double	*umagb;		/* Array of blue magnitudes (returned) */
int	*utype;		/* Array of spectral types (returned) */
int	nlog;		/* Logging interval */
{
    char srchurl[LINE];
    char numlist[LINE];
    char numstr[32];
    char csys[32];
    struct TabTable *tabtable;
    int i, refcat;
    char title[64];	/* Description of catalog (returned) */
    int syscat;		/* Catalog coordinate system (returned) */
    double eqcat;	/* Equinox of catalog (returned) */
    double epcat;	/* Epoch of catalog (returned) */
    char cstr[32];
    char temp[64];
    struct StarCat *starcat;

    /* Make list of catalog numbers */
    for (i = 0; i < nnum; i++) {
	refcat = RefCat (refcatname, title, &syscat, &eqcat, &epcat);
	CatNum (refcat, 0, 0, unum[i], numstr);
	if (i > 0) {
	    strcat (numlist, ",");
	    strcat (numlist, numstr);
	    }
	else
	    strcpy (numlist, numstr);
	}

    /* Set up search query */
    wcscstr (cstr, sysout, eqout, epout);
    sprintf (srchurl, "?catalog=%s&num=%soutsys=%s&",refcatname,numlist,csys);
    if (epout != 0.0) {
	sprintf (temp, "epoch=%.5f&", epout);
	strcat (srchurl, temp);
	}

    /* Run search across the web */
    if ((tabtable = webopen (caturl, srchurl, nlog)) == NULL) {
	if (nlog > 0)
	    fprintf (stderr, "WEBRNUM: %s failed\n", srchurl);
	return (0);
	}

    /* Return if no data */
    if (tabtable->tabdata == NULL || strlen (tabtable->tabdata) == 0) {
	if (nlog > 0)
	    fprintf (stderr, "WEBRNUM: No data returned\n");
	return (0);
	}

    /* Open returned Starbase table as a catalog */
    if ((starcat = tabcatopen (caturl, tabtable)) == NULL) {
	if (nlog > 0)
	    fprintf (stderr, "WEBRNUM: Could not open Starbase table as catalog\n");
	return (0);
	}

    /* Extract desired sources from catalog  and return them */
    return (tabrnum (srchurl, nnum, sysout, eqout, epout, &starcat,
         unum, ura, udec, upra, updec, umag, umagb, utype, NULL, nlog));
}


static struct TabTable *
webopen (caturl, srchpar, nlog)

char	*caturl;	/* URL of search engine */
char	*srchpar;	/* Search engine parameters to append */
{
    File sok;
    char *server;
    char linebuff[LINE];
    char *buff;
    char *cgipart;
    char *srchurl;
    char *servurl;
    int	status;
    int lserver, lsrch;
    char *tabbuff, *newbuff;
    int	lbuff = 0;
    int lfname, lfa;
    char *tabnew, *tabline, *lastline;
    int formfeed = (char) 12;
    struct TabTable *tabtable;
    int nc;
    int chunked = 0;
    int lchunk, ltab, lline, lname;
    int diag;

    if (nlog == 1)
	diag = 1;
    else
	diag = 0;

    /* Extract server name from search engine URL */
    servurl = caturl;
    if (!strncmp(caturl, "http://", 7))
	servurl = servurl + 7;
    cgipart = strchr (servurl, '/');
    lserver = cgipart - servurl;
    if ((server = (char *) malloc (lserver+2)) == NULL)
	return (NULL);
    strncpy (server, servurl, lserver);
    server[lserver] = (char) 0;

    /* Combine CGI command and arguments */
    lsrch = strlen (srchpar) + strlen (cgipart) + 2;
    if ((srchurl = (char *) malloc (lsrch)) == NULL)
	return (NULL);
    strcpy (srchurl, cgipart);
    strcat (srchurl, srchpar);

    /* Open port to HTTP server */
    if ( !(sok = SokOpen(server, 80, XFREAD | XFWRITE)) )
	return (NULL);

    /* Send HTTP command */
    fprintf(sok, "GET %s HTTP/1.1\nHost: %s\n\n", srchurl, server);
    fflush(sok);

    fscanf(sok, "%*s %d %*s\n", &status);

    /* Skip continue lines
    if (status == 100) {
	while (status == 100)
	    fscanf(sok, "%*s %d %*s\n", &status);
	} */

    /* If status is not 200 return without data */
    if ( status != 200 )
	return (NULL);

    /* Skip over http header of returned stuff */
    while (fgets (linebuff, LINE, sok) ) {
	if (diag)
	    fprintf (stderr, "%s", linebuff);
	if (strsrch (linebuff, "chunked") != NULL)
	    chunked = 1;
	if (*linebuff == '\n') break;
	if (*linebuff == '\r') break;
	}

    /* Read table into buffer in memory a chunk at a time */
    tabbuff = NULL;
    nc = 0;
    if (chunked) {
	lchunk = 1;
	lbuff = 0;
	while (lchunk > 0) {
	    fgets (linebuff, LINE, sok);
	    lline = strlen (linebuff);
	    if (lline < 1)
		break;
	    if (linebuff[lline-1] < 32)
		linebuff[lline-1] = (char) 0;
	    if (linebuff[lline-2] < 32)
		linebuff[lline-2] = (char) 0;
	    lline = strlen (linebuff);
	    if (lline > 0)
		lchunk = (int) strtol (linebuff, NULL, 16);
	    else
		lchunk = 0;
	    if (diag)
		fprintf (stderr, "%s (=%d)\n", linebuff, lchunk);
	    if (lchunk > 0) {
		lbuff = lbuff + lchunk;
		if (tabbuff == NULL) {
		    tabbuff = (char *) malloc (lbuff+8);
		    buff = tabbuff;
		    }
		else {
		    newbuff = (char *) malloc (lbuff+8);
		    strcpy (newbuff, tabbuff);
		    free (tabbuff);
		    tabbuff = newbuff;
		    buff = tabbuff + lbuff - lchunk;
		    }
        	fread (buff, 1, lchunk, sok);
		buff[lchunk] = (char) 0;
		if (diag)
		    fprintf (stderr, "%s\n", buff);
		}
	    }
	}

    /* Read table into buffer in memory a line at a time */
    else {
	while (fgets (linebuff, LINE, sok)) {
	    if (diag)
        	fprintf (stderr, "%s", linebuff);
	    nc = nc + strlen (linebuff);
	    if (tabbuff == NULL) {
		lbuff = 100 * strlen (linebuff);
		tabbuff = (char *) calloc (lbuff, 1);
		strcpy (tabbuff, linebuff);
		}
	    else if (nc > lbuff) {
		lbuff = lbuff * 2;
		newbuff = (char *) calloc (lbuff, 1);
		strcpy (newbuff, tabbuff);
		free (tabbuff);
		tabbuff = newbuff;
		newbuff = NULL;
		strcat (tabbuff, linebuff);
		}
	    else
		strcat (tabbuff, linebuff);
	    }
	}
    (void) fclose (sok);

    /* Allocate tab table structure */
    ltab = sizeof (struct TabTable);
    if ((tabtable = (struct TabTable *) calloc (1, ltab)) == NULL) {
	fprintf (stderr,"WEBOPEN: cannot allocate Tab Table structure for %s",
	         srchurl);
	return (NULL);
	}

    /* Save pointers to file contents */
    tabtable->tabbuff = tabbuff;
    tabtable->tabheader = tabtable->tabbuff;

    /* Allocate space for and save catalog URL as filename */
    lname = strlen (caturl) + 2;
    if ((tabtable->filename = (char *) calloc (1, lname)) == NULL) {
	fprintf (stderr,"WEBOPEN: cannot allocate filename %s in structure",
	         caturl);
	tabclose (tabtable);
	return (NULL);
	}
    strcpy (tabtable->filename, caturl);

    /* Allocate space for and save search string as tabname */
    lname = strlen (srchpar) + 2;
    if ((tabtable->tabname = (char *) calloc (1, lname)) == NULL) {
	fprintf (stderr,"WEBOPEN: cannot allocate tabname %s in structure",
	         srchurl);
	tabclose (tabtable);
	return (NULL);
	}
    strcpy (tabtable->tabname, srchpar);

    /* Find column headings and start of data */
    tabline = tabtable->tabheader;
    lastline = NULL;
    while (*tabline != '-' && tabline < tabtable->tabbuff+lbuff) {
	lastline = tabline;
	tabline = strchr (tabline,newline) + 1;
	}
    if (*tabline != '-') {
	fprintf (stderr,"WEBOPEN: No - line in tab table %s",srchurl);
	tabclose (tabtable);
	return (NULL);
	}
    tabtable->tabhead = lastline;
    tabtable->tabdata = strchr (tabline, newline) + 1;

    /* Extract positions of keywords we will want to use */
    if (!tabparse (tabtable)) {
	fprintf (stderr,"TABOPEN: No columns in tab table %s\n",srchurl);
	tabclose (tabtable);
	return (NULL);
	}

    /* Enumerate entries in tab table catalog by counting newlines */
    tabnew = tabtable->tabdata;
    tabtable->nlines = 0;
    while ((tabnew = strchr (tabnew, newline)) != NULL) {
	tabnew = tabnew + 1;
	tabtable->nlines = tabtable->nlines + 1;
	if (*tabnew == formfeed)
	    break;
	}

    tabtable->tabline = tabtable->tabdata;
    tabtable->iline = 1;
    return (tabtable);
}


/* sokFile.c
 * copyright 1991, 1993, 1995, 1999 John B. Roll jr.
 */

static FILE *
SokOpen(name, port, mode)
	char *name;             /* "host:port" socket to open */
	int   port;
	int   mode;             /* mode of socket to open */
{
    int             xfd;        /* socket descriptor */
    int             type;       /* type returned from FileINetParse */
    struct sockaddr_in adrinet; /* socket structure parsed from name */
    int             reuse_addr = 1;


    File            f;          /* returned file descriptor */

    if (!(type = FileINetParse(name, port, &adrinet)))
	return NULL;

    if ( type == 1 
     || (mode & XFCREAT && mode & XFREAD && !(mode & XFWRITE)) ) {
	if ( ((xfd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
	  ||  setsockopt(xfd, SOL_SOCKET, SO_REUSEADDR,
	             (char *) &reuse_addr, sizeof(reuse_addr)) < 0
	  || (bind(xfd, (struct sockaddr *) & adrinet
	                 ,sizeof(adrinet)) != 0)
	  ||  listen(xfd, 5) ) {
	    close(xfd);
	    return NULL;
	}
      } else {
	if (((xfd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
	           || (connect(xfd, (struct sockaddr *) & adrinet
	                       ,sizeof(adrinet)) != 0)) {
	    close(xfd);
	    return NULL;
	}
    }

    f = fdopen (xfd, "r+");

    return f;
}


static int
FileINetParse(file, port, adrinet)
	char *file;             /* host/socket pair to parse? */
	int   port;
	struct sockaddr_in *adrinet; /* socket info structure to fill? */
{
    struct hostent *hp;         /* -> hostent structure for host */
    char            hostname[MAXHOSTNAMELENGTH + 12]; /* name of host */
    char           *portstr;    /* internet port number (ascii) */
    int             type = 2;   /* return code */
    extern int gethostname();

    if ( !strncmp(file, "http://", 7) ) {
	file += 7;
	if ( port == -1 ) port  = 80;
    }

    strcpy(hostname, file);

#ifdef msdos
    /* This is a DOS disk discriptor, not a machine name */
    if ((!(file[0] == '.')) && file[1] == ':')
	return 0;
#endif

    if ( portstr = strchr(hostname, '/') ) {
	*portstr = '\0';
    }

    if ( portstr = strchr(hostname, ':') ) {
	*portstr++ = '\0';

	if ((port = strtol(portstr, NULL, 0)) == 0) {
	    struct servent *service;

	    if ((service = getservbyname(portstr, NULL)) == NULL)
	        return 0;
	    port = service->s_port;
	}
    }

    if ( port == -1 ) return 0;

    if (hostname[0] == '\0')
	type = 1;
    if (hostname[0] == '\0' || hostname[0] == '.')
	if (gethostname(hostname, MAXHOSTNAMELENGTH) == -1)
	    return 0;

    if ((hp = gethostbyname(hostname)) == NULL)
	return 0;

    memset(adrinet, 0, sizeof(struct sockaddr_in));
    adrinet->sin_family = AF_INET;
    adrinet->sin_port = htons(port);
    memcpy(&adrinet->sin_addr, hp->h_addr, hp->h_length);

    return type;
}

/* Nov 29 2000	New subroutines
 * Dec 11 2000	Do not print messages unless nlog > 0
 * Dec 12 2000	Fix problems with return if no stars
 * Dec 18 2000	Clean up code after lint
 *
 * Jan  2 2001	Set MAXHOSTNAMELENGTH to 256, bypassing system constant
 * Jan  3 2001	Include string.h, not strings.h
 */
