/* File webbuff.c
 * by Doug Mink, after John Roll
 * July 12, 2001
 *
 * Return character buffer from given URL
 */

char *
webbuff (url, diag)

char	*url;	/* URL to read */
int	diag;	/* 1 to print diagnostic messages */
{
    File sok;
    char *server;
    char linebuff[LINE];
    char *buff;
    char *urlpath;
    char *servurl;
    int	status;
    int lserver;
    char *tabbuff, *newbuff;
    int	lbuff = 0;
    int lfname, lfa;
    int nc;
    int chunked = 0;
    int lchunk, ltab, lline, lname;

    /* Extract server name and path from URL */
    servurl = url;
    if (!strncmp(url, "http://", 7))
	servurl = servurl + 7;
    urlpath = strchr (servurl, '/');
    lserver = urlpath - servurl;
    if ((server = (char *) malloc (lserver+2)) == NULL)
	return (NULL);
    strncpy (server, servurl, lserver);
    server[lserver] = (char) 0;

    /* Open port to HTTP server */
    if ( !(sok = SokOpen(server, 80, XFREAD | XFWRITE)) )
	return (NULL);

    /* Send HTTP command */
    fprintf(sok, "GET %s HTTP/1.1\nHost: %s\n\n", urlpath, server);
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
	lline = 1;
	lbuff = 0;
	while (lline > 0) {
	    fgets (linebuff, LINE, sok);
	    lline = strlen (linebuff);
	    if (lline < 1)
		break;
	    if (linebuff[lline-1] < 32)
		linebuff[lline-1] = (char) 0;
	    if (linebuff[lline-2] < 32)
		linebuff[lline-2] = (char) 0;
	    if (strlen (linebuff) > 0) {
		lchunk = (int) strtol (linebuff, NULL, 16);
		if (lchunk < 1)
		    break;
		}
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

    return (tabbuff);
}
