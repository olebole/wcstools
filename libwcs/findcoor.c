/* File findcoor.c
 * February 16, 1996
 * By Joe Mohr, modified by Doug Mink
 *
 * This subroutine takes a list of reference star positions and a list
 * of object coordinates from an image and finds the image pixels
 * which correspond to each of the reference stars.  The program produces
 * an exhaustive list of triangles using the points from each of the two
 * files and them by comparing the distribution of scale changes and 
 * orientation twists it determines which triangles are correct and which
 * are just random similar triangles.  As a final step to exclude the
 * change similar triangle with the correct scale and correct orientation
 * the program sorts the list of matches points for each gsc star and chooses
 * the median value.  Thus, this puppy is pretty darned robust.
 * It should work regardless of whether the image is transposed and rotated
 * with respect to the gsc.  A list of the gsc stars with matches image
 * points is output... these can be fit using <CCD2sky>.
 *
 * Hey Brian... thanks for the neat triangle idea... but this code is
 * a heap "crisper" than yours!  Yeee haw!!!
 * October 17, 1994
 */

#include <stdio.h>
#include "wcs.h"

#define DECL(a,b,c) (a/fabs(a)*(3600.0*fabs(a)+60.0*b+c))
#define RA(a,b,c) (54000.0*a+900.0*b+15.0*c)
#define Squ(x) ((x)*(x))
#define DIST(a,b,c,d) (sqrt(Squ((a)-(c))+Squ((b)-(d))))
#define ORIENT(a,b,c,d) (atan(fabs(((b)-(d))/((a)-(c)))))

#define DELTA 1.015	/* agreement required for side ratios of triangles*/

#define LOW_LIMIT 0.05	/* minimum width of scale and orient dist is 5% */

#define MAXGSC	100	/* Maximum number of Guide Stars to use in an image */

static void hunt ();

typedef struct {unsigned a,b,c;} triad;

int
findCoords (cra, cdec, ng, gra, gdec, gm, gx, gy, ns, sx, sy, secpix, rot)

double	cra;	/* Center right ascension of reference star search */
double	cdec;	/* Center declination of reference star search */
int	ng;	/* Number of reference stars found in image */
double	*gra;	/* reference star right ascensions */
double	*gdec;	/* reference star declinations */
double	*gm;	/* reference star magnitudes */
double	*gx;	/* reference star X coordinates (returned) */
double	*gy;	/* reference star Y coordinates (returned) */
int	ns;	/* Number of stars found in image */
double	*sx; 	/* Image star X-coordinates */
double	*sy; 	/* Image star Y-coordinates */
double	*secpix;	/* Constrain scale to be within 5% of this if not 0.0 */
double	*rot;	/* Constrain twist to be within 2 deg of this if not 0.0 */

{
    float	skdec[MAXGSC],skra[MAXGSC], *skr,*focr, r12,r23,r31, val,
		*scale,*orient,width,median, ormedian,orwidth, *skmatches[50];
    int		skrat,forat,stnum,sknum,i,j,k,n=0,getline(),loc,*foin,
		tnum,p1,p2,p3,p4,match=0,sklist[MAXGSC];
    unsigned	sz=0;
    triad	*skpt,*fopt,*matchlist;
    void	indexx(),sort2(),sort(),locate(),calc_scales(), get_mode();

    /* Change input reference star positions to relative positions */
    if (ng <= MAXGSC)
	sknum = ng;
    else
	sknum = MAXGSC;

    /* Convert to relative coordinates */
    for (i = 0; i < sknum; i++) {
	skdec[i] = gdec[i] - cdec;
	skra[i] = cos (degrad (gdec[i])) * (gra[i] - cra);
	sklist[i] = 0;	/*mark the list of stars for later matching*/
	}

    /* Prepare enough memory for the calculations */
    /* Empirically determine how many triangle there will be */
    for (sz = 0, i = 0; i < sknum; i++)
	for(j = i + 1; j < sknum; j++)
	    for(k = j + 1; k < sknum; k++)
		sz++;
    skpt = (triad *) calloc (sz,sizeof(triad));
    skr = (float *) calloc (3*sz,sizeof(float));
    matchlist = (triad *)calloc (2000*sz,sizeof(triad));
    printf ("%d triads in skpt, %d in matchlist\n", sz, 20*sz);

    /* Calculate the side ratios of every combination*/
    for (n = loc = i = 0; i < sknum; i++) {
	for (j = i + 1; j < sknum; j++) {
	    r12 = DIST(skdec[i],skra[i],skdec[j],skra[j]);
	    for(k = j + 1; k < sknum; k++) {

		/* Calculate side lengths */
		r23 = DIST(skdec[j],skra[j],skdec[k],skra[k]); 
		r31 = DIST(skdec[k],skra[k],skdec[i],skra[i]); 

		/* Record the i,j,k of the three points */
		skpt[loc].a = i;
		skpt[loc].b = j;
		skpt[loc++].c = k;

		/* Calculate the three side ratios */
		skr[n++] = r12 / r23;
		skr[n++] = r23 / r31;
		skr[n++] = r31 / r12;
		}
	    }
	}
    skrat = n;
    printf("  %d ratios for %d triangles for %d stars\n", n,loc,sknum);

    /* Process the image object list */
    stnum = ns;

    /* Prepare enough memory for the calculations */
    /* Empirically determine how many triangle there will be */
    for (sz = 0, i = 0; i < stnum; i++)
	for (j = i + 1; j < stnum; j++)
	    for (k = j + 1; k < stnum; k++)
		sz++;
    fopt = (triad *) calloc (2*sz,sizeof(triad));
    foin = (int *) calloc (6*sz,sizeof(int));
    focr = (float *) calloc (6*sz,sizeof(float));

    /* Calculate the side ratios of every combination */
    for (n = loc = i = 0; i < stnum; i++) {
	for (j = i + 1; j < stnum; j++) {
	    r12 = DIST (sx[i], sy[i], sx[j], sy[j]);
	    for (k = j + 1; k < stnum; k++) {

		/* Calculate side lengths */
		r23 = DIST (sx[j], sy[j], sx[k], sy[k]);
		r31 = DIST (sx[k], sy[k], sx[i], sy[i]);

		/* Record the i,j,k of the three points */
		fopt[loc].a = i;
		fopt[loc].b = j;
		fopt[loc++].c = k;

		/* Produce inverted version of the same triangle */
		fopt[loc].a = i;
		fopt[loc].b = k;
		fopt[loc++].c = j;

		/* Calculate the three side ratios for both triangles */
		focr[n++] = r12 / r23;
		focr[n++] = r23 / r31;
		focr[n++] = r31 / r12;
		focr[n++] = r31 / r23;
		focr[n++] = r23 / r12;
		focr[n++] = r12 / r31;
		}
	    }
	}
    forat = n;
    indexx (forat,focr-1,foin-1);
    printf ("  sorted %d ratios for %d triangles for %d objects\n  ",
	    n,loc,stnum);

    /* Go through reference star triangle list looking for matches */
    n = 0;
    match = 0;
    do {
	r12 = skr[n]*DELTA;
	r23=skr[n+1]*DELTA;
	r31=skr[n+2];

	/* Search for matching triangles in the primary ratio list */
	hunt (focr, foin, forat, r12*(2./DELTA-1), &p1);
	p1++;
	while (p1 >= 0 && p1 < forat && focr[foin[p1]-1] < r12) {

	    /* Find out which point of which triangle we are considering */
	    matchlist[match].b = (foin[p1]-1) / 3;	/* image triangle */
	    matchlist[match].c = (foin[p1]-1) % 3;	/* focas point */

	    /* Find second point of triangle */
	    p3 = matchlist[match].c + 1;
	    if (p3 > 2)
		p3 = 3 * matchlist[match].b;
	    else
		p3 = 3 * matchlist[match].b + p3;

	    /* Find third point of triangle */
	    p4 = matchlist[match].c + 2;
	    if (p4 > 2)
		p4 -= 3;
	    p4 = 3 * matchlist[match].b + p4;

	    /* Get list of ratios consistent with second ratio */
	    val = r23 * ((2.0 / DELTA) - 1.0);
	    hunt (focr, foin, forat, val, &p2);
	    p2++;

	    /* Record if this match is the next point in the triangle */
	    while (p2 >= 0 && p2 < forat && focr[foin[p2]-1] < r23) 
		if (foin[p2++]-1 == p3 && fabs(focr[p4]-r31)/r31 < DELTA-1.0)
		    matchlist[match++].a = n/3; /* reference star triangle */
	    p1++;
	    }
	/* Next triangle in the reference star list */
	if (n%60 == 0) {
	    /* printf("."); */
	    printf ("%d %d\n",n,match);
	    fflush (stdout);
	    }
	n += 3;
	} while (n < skrat);
    printf("\n  %d triangle matches in list\n", match);
    fflush(stdout);
    if (match < 1) return (0);

    /* Use match lists to calculate distribution of scales*/
    scale = (float *) calloc(3*match,sizeof(float));
    orient = (float *) calloc(3*match,sizeof(float));
    calc_scales (skpt,skra,skdec,fopt,sx,sy,matchlist,match,scale,orient);

    /* Search for most probably scale value */
    get_mode (scale, match, &median, &width);
    sort2 (match,scale-1,orient-1);

    /* If scale is specified, change the median scale to that value */
    if (*secpix != 0.0) {
	if (fabs (median - *secpix) / *secpix > LOW_LIMIT)
	    median = *secpix;
	if (width > 0.05 * *secpix)
	    width = LOW_LIMIT * *secpix;
	}

    /* Use the scale constraints to work on the orientation constraints */
    locate (scale-1, match, median-width, &p1);
    locate (scale-1, match, median+width, &p2);

    /* Find the mode of the orientation distribution */
    if (*rot != 0.0)
	ormedian = *rot;
    else
	get_mode (&orient[p1], p2-p1, &ormedian, &orwidth);
    orwidth = 0.034906; /*set 2 degree width*/

    /* Print results for correct scale */
    printf("  median scale:   %8.5f (%8.5f)\n",median,width);
    printf("  median twist:   %8.5f (%8.5f)\n",ormedian*180./PI,orwidth*180./PI);
    fflush(stdout);
    *secpix = median;
    *rot = ormedian;

    /* Pick out correct matches */
    calc_scales(skpt,skra,skdec,fopt,sx,sy,matchlist,match,scale,orient);

    /* Count the number of "correct" matches */
    for (tnum = i = 0; i < match; i++) 
	if (fabs(scale[i]-median)<width && 
	    fabs(orient[i]-ormedian)<orwidth) tnum++;
    printf("  %d of %d triangles are probably correct matches\n",
	   tnum,match);

    /* Set up array to track the "correct" matches for reference points */
    for (i = 0; i < sknum; i++) 
	skmatches[i] = (float *)calloc(tnum,sizeof(float));
    for (i = 0; i < match; i++) { 

	/* For each triangle with correct scale and orientation */
	if (fabs(scale[i] - median) < width && 
	    fabs(orient[i] - ormedian) < orwidth) 

	    /* Cycle through all three points */
	    for (k = 0; k < 3; k++) {
		p1 = *(&skpt[matchlist[i].a].a+k);
		p2 = matchlist[i].c+k;if (p2>2) p2-=3;
		p2 = *(&fopt[matchlist[i].b].a+p2);
		*(skmatches[p1]+sklist[p1]++) = p2;
		}
	}

    /* Now sorting the match lists for each reference star-> taking median */
    for (tnum = j = i = 0; i < sknum; i++)
	if (sklist[i] > 0) {
	    tnum += sklist[i];
	    k++;
	    }
    tnum /= k;
    for (i = 0; i < sknum; i++)
	if (sklist[i] > 0.5 * tnum) {
	    if (sklist[i] > 1)
		sort(sklist[i], skmatches[i]-1);

	    /* Save the median as the correct match */
	    sklist[i] = *(skmatches[i] + sklist[i] / 2);
	    }
	else
	    sklist[i] = -1;
	
    /* Write out image position for each reference star */
    p1 = i = 0;
    for (i = 0; i < sknum; i++) {
	if (sklist[i] >= 0) {
	    gx[i] = sx[sklist[i]];
	    gy[i] = sy[sklist[i]];
	    printf ("   %8.3f %8.3f\n",gra[i],gdec[i],sx[sklist[i]],sy[sklist[i]]);
	    p1++;
	    }
	else {
	    gx[i] = -1.0;
	    gy[i] = -1.0;
	    }
	i++;
	}
    printf("  %d stars out of %d have matches in the image\n",p1,sknum);
    return (p1);
}


/* Hunt through a list which has a sorted index list */

static void
hunt (list, ind, len, val, location)

float *list;
int *ind;
int len;
float val;
int *location;

{
    int	ascnd,ju,jm,jl;

    jl = -1;
    ju = len;
    ascnd = list[ind[len-1]-1] > list[ind[0]-1]; /*ascending order?*/
    while (ju - jl > 1) {
	jm = (ju + jl) >> 1;
	if (val > list[ind[jm]-1] == ascnd)
	    jl = jm;
	else
	    ju = jm;
	}
    *location = jl;
}


/* Go down match list and calculate scales */

void
calc_scales (sobj, sposx, sposy, fobj, fposx, fposy, pt, length, scale, orient)

triad sobj[];
float sposx[];
float sposy[];
triad fobj[];
double fposx[];
double fposy[];
triad pt[];
int length;
float scale[];
float orient[];

{
    int	i,p1,p2;

    for (i = 0; i < length; i++) {

	/* Find distance between first two points in GSC triangle */
	p1 = sobj[pt[i].a].a;
	p2 = sobj[pt[i].a].b;
	scale[i] = DIST (sposx[p1], sposy[p1], sposx[p2], sposy[p2]);
	orient[i] = ORIENT (sposx[p1], sposy[p1], sposx[p2], sposy[p2]);

	/* Find distance between correct points in image triangle */
	p1 = *(&fobj[pt[i].b].a + pt[i].c);
	p2 = pt[i].c + 1;
	if (p2 > 2) p2-=3;
	p2 = *(&fobj[pt[i].b].a + p2);
	scale[i] /= DIST (fposx[p1], fposy[p1], fposx[p2], fposy[p2]);
	orient[i] -= ORIENT (fposx[p1], fposy[p1], fposx[p2], fposy[p2]);
	}
}


/* Get the mode of a distribution */

void
get_mode (list, num, mode, width)

float	*list;
int	num;
float	*mode;
float	*width;
{
    float	min, max, *slist, peak = 0.0;
    int		j, i, bins, dist[5000];
    void	sort(), locate();

    if (num == 0) {
	*mode = *width=0.0;
	return;
	}
    if (num ==1 || num == 2) {
	*mode = list[0];
	if (num == 1)
	    *width = LOW_LIMIT * fabs (*mode);
	else
	     *width = (1. + LOW_LIMIT) * fabs (list[1] - list[0]);
	return;
	}

    /* Find maximum and minimum values for distribution */
    slist = (float *) calloc(num,sizeof(float));
    min = max = list[0];
    for (i = 0; i < num; i++)
	if ((slist[i] = list[i]) > max)
	    max=list[i]; 
	else if (list[i] < min)
	    min = list[i];

    if (num < 10) {
	sort (num, slist-1);
	*mode = slist[num/2];
	*width = 0.0;
	i = 0;
	while (i < num)
	    *width += fabs (slist[i++] - *mode);
	*width /= (float) num;
	if (*width < LOW_LIMIT * fabs(*mode))
	    *width = LOW_LIMIT * fabs(*mode);
	return;
	}

    /* Make distributions only if there are 10 or more objects */
    bins = num / 3;

    /* Find the appropriate number of bins to find the mode */
    /* That is the number of bins which yield a peak 50% above the rest */
    do {
	bins *= 1.5;

	/* clean distribution */
	i = 0;
	while (i<bins)
	    dist[i++] = 0.0;

	/* build distribution */
	for (i = 0; i < num; i++)
	    dist[(int)(bins*(list[i]-min)/(max-min))] += 1.0;

        /* Search for the most probably value */
        sort (bins, dist-1);
        } while (dist[bins-2]/dist[bins-1]>0.6 && bins<5000);

    /* Rebuild the distribution*/
    i = 0;
    while (i < bins)
	dist[i++]=0.0;
    for (i = 0; i < num; i++)
	dist[(int)(bins*(list[i]-min)/(max-min))]+=1.0;
    i = 0;
    while (i++<bins) {
	if (dist[i-1]>peak) {
	    peak=dist[i-1];
	    j=i-1;
	    }
	}
    i = j;
    while (i++ < bins) {
	if (dist[i]/dist[i-1] > 1.0 || dist[i] == 0)
	    break;
	}
    *mode = min + ((float) j + 0.5) * (max - min) / (float) bins;
    *width = (float) (i - j) * (max - min) / (float) bins;
    sort (num, slist-1);
    locate (slist-1, num, *mode-*width, &i);
    locate (slist-1, num, *mode+*width, &j);
    *mode = slist[(i+j)/2];
    if (*width < LOW_LIMIT*fabs(*mode))
	*width = LOW_LIMIT*fabs(*mode);
    cfree ((char *)slist);
}


void sort (n,ra)

int n;
float ra[];

{
    int l,j,ir,i;
    float rra;

    l = (n >> 1)+1;
    ir = n;
    for (;;) {
    	if (l > 1)
    	    rra = ra[--l];
    	else {
    	    rra = ra[ir];
    	    ra[ir] = ra[1];
    	    if (--ir  == 1) {
    		ra[1] = rra;
    		return;
    		}
    	    }
    	i = l;
    	j = l << 1;
    	while (j <=  ir) {
    	    if (j < ir && ra[j] < ra[j+1]) ++j;
    	    if (rra < ra[j]) {
    		ra[i] = ra[j];
    		j += (i=j);
    		}
    	    else j = ir+1;
    	    }
    	ra[i] = rra;
	}
}


void sort2 (n, ra, rb)

int n;
float ra[],rb[];

{
    int l,j,ir,i;
    float rrb,rra;

    l = (n >> 1)+1;
    ir = n;
    for (;;) {
    	if (l > 1) {
    	    rra = ra[--l];
    	    rrb = rb[l];
    	    }
	else {
    	    rra = ra[ir];
    	    rrb = rb[ir];
    	    ra[ir] = ra[1];
    	    rb[ir] = rb[1];
    	    if (--ir  == 1) {
    		ra[1] = rra;
    		rb[1] = rrb;
    		return;
    		}
    	    }
    	i = l;
    	j = l << 1;
    	while (j <=  ir)	{
    	    if (j < ir && ra[j] < ra[j+1]) ++j;
    	    if (rra < ra[j]) {
    		ra[i] = ra[j];
    		rb[i] = rb[j];
    		j += (i=j);
    		}
    	    else j = ir+1;
    	    }
    	ra[i] = rra;
    	rb[i] = rrb;
	}
}

void locate (xx,n,x,j)
float xx[],x;
int n,*j;

{
    int ascnd,ju,jm,jl;

    jl = 0;
    ju = n+1;
    ascnd = xx[n] > xx[1];
    while (ju-jl > 1) {
    	jm = (ju+jl) >> 1;
    	if (x > xx[jm]  == ascnd)
    	    jl = jm;
    	else
    	    ju = jm;
	}
    *j = jl;
}

void indexx (n,arrin,indx)

int n,indx[];
float arrin[];

{
    int l,j,ir,indxt,i;
    float q;

    for (j = 1; j <= n; j++)
	indx[j]=j;
    l = (n >> 1) + 1;
    ir = n;
    for (;;) {
    	if (l > 1)
    	    q = arrin[(indxt=indx[--l])];
    	else {
    	    q = arrin[(indxt=indx[ir])];
    	    indx[ir] = indx[1];
    	    if (--ir  == 1) {
    		indx[1] = indxt;
    		return;
    		}
    	    }
    	i = l;
    	j = l << 1;
    	while (j <=  ir) {
    	    if (j < ir && arrin[indx[j]] < arrin[indx[j+1]])
		j++;
    	    if (q < arrin[indx[j]]) {
    		indx[i] = indx[j];
    		j += (i=j);
    		}
    	    else j = ir+1;
    	    }
    	indx[i] = indxt;
	}
}

