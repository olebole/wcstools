
/* The following are used in star finding (findstar.c) */
#define	NSTATPIX	25	/* Stats are computed for +- this many pixels */
#define	ISTATPIX	10	/* Stats are computed every this many pixels */
#define	MAXWALK		20	/* Farthest distance to walk from seed */
#define	BURNEDOUT	65535	/* Clamp pixels brighter than this */
#define NITERATE	5	/* Number of iterations for sigma clipping */
#define STARSIGMA	5.0	/* Stars must be this many sigmas above mean */
#define BORDER		10	/* Ignore this much of the edge */
#define MAXRAD		20	/* Maximum radius for a star */
#define MINRAD		1	/* Minimum radius for a star */
#define MINPEAK		10	/* Minimum peak for a star */
#define MINSEP		10	/* Minimum separations for stars */

/* The following are used in star matching (matchstar.c) */
#define	FTOL	0.000001	/* Fractional change of chisqr() to be done */
#define NMAX		500	/* Maximum number of minimization iterations */
#define	NPEAKS		20	/* Binning peak history */

/* The following are used in world coordinate system fitting (imsetwcs.c) */
#define MINSTARS	4	/* Minimum stars from reference and image */
#define MAXSTARS	25	/* Max star pairs we need to try using */
#define MAXREF		100	/* Max reference stars to use in image region */
#define MAGLIM1		-2	/* Faintest reference catalog magnitude to use*/
#define MAGLIM		0	/* Faintest reference catalog magnitude to use*/
#define PIXDIFF		10	/* +- this many pixels is a match */
#define PSCALE		0	/* Plate scale in arcsec/pixel */
				/* (if nonzero, this overrides image header) */
