/*** File wcscon.c
 *** Doug Mink, Harvard-Smithsonian Center for Astrophysics
 *** Based on Starlink subroutines by Patrick Wallace
 *** April 1, 1996

 * Module:	wcscon.c (World Coordinate System conversion)
 * Purpose:	Convert between various sky coordinate systems
 * Subroutine:	fk524 (ra,dec) Convert J2000(FK5) to B1950(FK4) coordinates
 * Subroutine:	fk524m (ra,dec,rapm,decpm)
 * Subroutine:	fk425 (ra,dec) Convert B1950(FK4) to J2000(FK5) coordinates
 * Subroutine:	fk425m (ra, dec, rapm, decpm)
 * Subroutine:	fk42gal (dtheta,dphi) Convert B1950(FK4) to galactic coordinates
 * Subroutine:	fk52gal (dtheta,dphi) Convert J2000(FK5) to galactic coordinates
 * Subroutine:	gal2fk4 (dtheta,dphi) Convert galactic coordinates to B1950(FK4)
 * Subroutine:	gal2fk5 (dtheta,dphi) Convert galactic coordinates to J2000<FK5)
 * Subroutine:  fk5prec (ep0, ep1, ra, dec) Precession ep0 to ep1, FK5 system
 * Subroutine:  fk4prec (ep0, ep1, ra, dec) Precession ep0 to ep1, FK4 system
 */

#include <math.h>
#ifndef VMS
#include <stdlib.h>
#endif
#include <ctype.h>
#include <string.h>

/*  Constant vector and matrix (by columns)
    These values were obtained by inverting C.Hohenkerk's forward matrix
    (private communication), which agrees with the one given in reference
    2 but which has one additional decimal place.  */

static double a[3] = {-1.62557e-6, -0.31919e-6, -0.13843e-6};
static double ad[3] = {1.245e-3,  -1.580e-3,  -0.659e-3};

/* FK524  convert J2000 FK5 star data to B1950 FK4
   based on Starlink sla_fk524 by P.T.Wallace 27 October 1987 */

static double emi[6][6] = {
	 0.999925679499910, -0.011181482788805, -0.004859004008828,
	-0.000541640798032, -0.237963047085011,  0.436218238658637,

	 0.011181482840782,  0.999937484898031, -0.000027155744957,
	 0.237912530551179, -0.002660706488970, -0.008537588719453,

	 0.004859003889183, -0.000027177143501,  0.999988194601879,
	-0.436101961325347,  0.012258830424865,  0.002119065556992,

	-0.000002423898405,  0.000000027105439,  0.000000011777422,
	 0.999904322043106, -0.011181451601069, -0.004858519608686,

	-0.000000027105439, -0.000002423927017,  0.000000000065851,
	 0.011181451608968,  0.999916125340107, -0.000027162614355,

	-0.000000011777422,  0.000000000065846, -0.000002424049954,
	 0.004858519590501, -0.000027165866691,  0.999966838131419};

void
fk524 (ra,dec)

double	*ra;		/* Right ascension in degrees (J2000 in, B1950 out) */
double	*dec;		/* Declination in degrees (J2000 in, B1950 out) */

{
double	rra;		/* Right ascension in radians (J2000 in, B1950 out) */
double	rdec;		/* Declination in radians (J2000 in, B1950 out) */
double	rapm;		/* Proper motion in right ascension */
double	decpm;		/* Proper motion in declination  */
			/* In:  rad/jul.yr.  Out: rad/trop.yr.  */
double	rcon = 1.74532925199433e-2; /* radians per degree */
void
fk524m();

	rra = *ra * rcon;
	rdec = *dec * rcon;
	rapm = (double) 0.0;
	decpm = (double) 0.0;
	fk524m (&rra, &rdec, &rapm, &decpm);
	*ra = rra / rcon;
	*dec = rdec / rcon;
	return;
}

void fk524m (ra,dec,rapm,decpm)

double	*ra;		/* Right ascension in radians (J2000 in, B1950 out) */
double	*dec;		/* Declination in radians (J2000 in, B1950 out) */
double	*rapm;		/* Proper motion in right ascension */
double	*decpm;		/* Proper motion in declination  */
			/* In:  rad/jul.yr.  Out: rad/trop.yr.  */

/*  This routine converts stars from the new, IAU 1976, FK5, Fricke
    system, to the old, Bessel-Newcomb, FK4 system, using Yallop's
    implementation (see ref 2) of a matrix method due to Standish
    (see ref 3).  The numerical values of ref 2 are used canonically.

 *  Notes:

      1)  The proper motions in ra are dra / dt rather than
 	    cos(dec) * dra / dt, and are per year rather than per century.
 
      2)  Note that conversion from Julian epoch 2000.0 to Besselian
 	    epoch 1950.0 only is provided for.  Conversions involving
 	    other epochs will require use of the appropriate precession,
 	    proper motion, and e-terms routines before and/or after
 	    fk524 is called.
 
      3)  In the fk4 catalogue the proper motions of stars within
 	    10 degrees of the poles do not embody the differential
 	    e - term effect and should, strictly speaking, be handled
 	    in a different manner from stars outside these regions.
 	    however, given the general lack of homogeneity of the star
 	    data available for routine astrometry, the difficulties of
 	    handling positions that may have been determined from
 	    astrometric fields spanning the polar and non - polar regions,
 	    the likelihood that the differential e - terms effect was not
 	    taken into account when allowing for proper motion in past
 	    astrometry, and the undesirability of a discontinuity in
 	    the algorithm, the decision has been made in this routine to
 	    include the effect of differential e - terms on the proper
 	    motions for all stars, whether polar or not.  at epoch 2000,
 	    and measuring on the sky rather than in terms of dra, the
 	    errors resulting from this simplification are less than
 	    1 milliarcsecond in position and 1 milliarcsecond per
 	    century in proper motion.

   References:

      1  "Mean and apparent place computations in the new IAU System.
          I. The transformation of astrometric catalog systems to the
 	  equinox J2000.0." Smith, C.A.; Kaplan, G.H.; Hughes, J.A.;
	  Seidelmann, P.K.; Yallop, B.D.; Hohenkerk, C.Y.
 	  Astronomical Journal vol. 97, Jan. 1989, p. 265-273.

      2  "Mean and apparent place computations in the new IAU System.
	  II. Transformation of mean star places from FK4 B1950.0 to
 	  FK5 J2000.0 using matrices in 6-space."  Yallop, B.D.;
	  Hohenkerk, C.Y.; Smith, C.A.; Kaplan, G.H.; Hughes, J.A.;
	  Seidelmann, P.K.; Astronomical Journal vol. 97, Jan. 1989,
	  p. 274-279.

      3  "Conversion of positions and proper motions from B1950.0 to the
	  IAU system at J2000.0", Standish, E.M.  Astronomy and
	  Astrophysics, vol. 115, no. 1, Nov. 1982, p. 20-22.

   P.T.Wallace   Starlink   27 October 1987
   Doug Mink     Smithsonian Astrophysical Observatory  7 June 1995 */

{
double	r2000,d2000;	/* J2000.0 ra,dec (radians) */
double	dr2000,dd2000;	/* J2000.0 proper motions (rad/jul.yr)*/
double	r1950,d1950;	/* B1950.0 ra,dec (rad) */
double	dr1950,dd1950;	/* B1950.0 proper motions (rad/trop.yr) */

/* Miscellaneous */
double	ur,ud;
double	sr,cr,sd,cd,x,y,z,w;
double	v1[6],v2[6];
double	xd,yd,zd;
double	rxyz,rxysq,rxy;
double	dra,ddec,scon,tcon;
int	i,j;
int	diag = 0;

/* Constants */
double	d2pi = 6.283185307179586476925287;	/* two PI */
double	pmf;	/* radians per year to arcsec per century */
double	tiny = 1.e-30; /* small number to avoid arithmetic problems */
double	zero = (double) 0.0;
double	rcon = 1.74532925199433e-2; /* radians per degree */

	pmf = 100. * 60. * 60. * 360. / d2pi;

/* Pick up J2000 data (units radians and arcsec / jc) */
	r2000 = *ra;
	d2000 = *dec;
	dr2000 = *rapm;
	dd2000 = *decpm;
	ur = dr2000 * pmf;
	ud = dd2000 * pmf;

/* Spherical to Cartesian */
	sr = sin (r2000);
	cr = cos (r2000);
	sd = sin (d2000);
	cd = cos(d2000);

	x = cr * cd;
	y = sr * cd;
	z = sd;

	v1[0] = x;
	v1[1] = y;
	v1[2] = z;

	if (ur != zero || ud != zero) {
	    v1[3] = -(ur*y) - (cr*sd*ud);
	    v1[4] =  (ur*x) - (sr*sd*ud);
	    v1[5] =	         (cd*ud);
	    }
	else {
	    v1[3] = zero;
	    v1[4] = zero;
	    v1[5] = zero;
	    }

/* Convert position + velocity vector to bn system */
	for (i = 0; i < 6; i++) {
	    w = zero;
	    for (j = 0; j < 6; j++) {
		w = w + emi[j][i] * v1[j];
		}
	    v2[i] = w;
	    }

/* Vector components */
	x = v2[0];
	y = v2[1];
	z = v2[2];
	xd = v2[3];
	yd = v2[4];
	zd = v2[5];

/* Magnitude of position vector */
	rxyz = sqrt (x*x + y*y + z*z);

/* Include e-terms */
	x = x + a[0] * rxyz;
	y = y + a[1] * rxyz;
	z = z + a[2] * rxyz;
	xd = xd + ad[0] * rxyz;
	yd = yd + ad[1] * rxyz;
	zd = zd + ad[2] * rxyz;

/* Convert to spherical */
	rxysq = x*x + y*y;
	rxy = sqrt (rxysq);

	if (x == zero && y == zero)
	    r1950 = zero;
	else {
	    r1950 = atan2 (y,x);
	    if (r1950 < zero)
		r1950 = r1950 + d2pi;
	    }
	d1950 = atan2 (z,rxy);

	if (rxy > tiny) {
	    ur = (x*yd - y*xd) / rxysq;
	    ud = (zd*rxysq - z * (x*xd + y*yd)) / ((rxysq + z*z) * rxy);
	    }
	dr1950 = ur / pmf;
	dd1950 = ud / pmf;

/* Return results */
	*ra = r1950;
	*dec = d1950;
	*rapm = dr1950;
	*decpm = dd1950;

	if (diag) {
	    scon = 3.6e3 / rcon;
	    tcon = 2.4e2 / rcon;
	    dra = tcon * (r1950 - r2000);
	    ddec = scon * (d1950 - d2000);
	    printf("B1950-J2000: dra= %11.5f sec  ddec= %f11.5f arcsec\n",
		dra, ddec);
	    }

	return;
}


/* Convert B1950.0 fk4 star data to J2000.0 fk5 */
static double em[6][6] = {
	 0.999925678186902,  0.011182059571766,  0.004857946721186,
	-0.000541652366951,  0.237917612131583, -0.436111276039270,

	-0.011182059642247,  0.999937478448132, -0.000027147426498,
	-0.237968129744288, -0.002660763319071,  0.012259092261564,

	-0.004857946558960, -0.000027176441185,  0.999988199738770,
	 0.436227555856097, -0.008537771074048,  0.002119110818172,

	 0.000002423950176,  0.000000027106627,  0.000000011776559,
	 0.999947035154614,  0.011182506007242,  0.004857669948650,

	-0.000000027106627,  0.000002423978783, -0.000000000065816,
	-0.011182506121805,  0.999958833818833, -0.000027137309539,

	-0.000000011776558, -0.000000000065874,  0.000002424101735,
	-0.004857669684959, -0.000027184471371,  1.000009560363559};

void
fk425 (ra,dec)

double	*ra;		/* Right ascension in degrees (B1950 in, J2000 out) */
double	*dec;		/* Declination in degrees (B1950 in, J2000 out) */

{
double	rra;		/* Right ascension in radians (B1950 in, J2000 out) */
double	rdec;		/* Declination in radians (B1950 in, J2000 out) */
double	rapm;		/* Proper motion in right ascension */
double	decpm;		/* Proper motion in declination  */
			/* In: rad/trop.yr.  Out:  rad/jul.yr. */
double	rcon = 1.74532925199433e-2; /* radians per degree */
void	fk425m();

	rra = *ra * rcon;
	rdec = *dec * rcon;
	rapm = (double) 0.0;
	decpm = (double) 0.0;
	fk425m (&rra, &rdec, &rapm, &decpm);
	*ra = rra / rcon;
	*dec = rdec / rcon;
	return;
}

void
fk425m (ra, dec, rapm, decpm)

double	*ra, *dec;	/* Right ascension and declination in radians
			   input:  B1950.0,fk4	returned:  J2000.0,fk5 */
double	*rapm, *decpm;	/* Proper motion in right ascension and declination
			   input:  B1950.0,fk4	returned:  J2000.0,fk5
			           rad/trop.yr.            rad/jul.yr.  */

/* This routine converts stars from the old, Bessel-Newcomb, FK4
   system to the new, IAU 1976, FK5, Fricke system, using Yallop's
   implementation (see ref 2) of a matrix method due to Standish
   (see ref 3).  The numerical values of ref 2 are used canonically.

   Notes:

      1)  The proper motions in ra are dra/dt rather than
 	   cos(dec)*dra/dt, and are per year rather than per century.

      2)  Conversion from besselian epoch 1950.0 to Julian epoch
 	   2000.0 only is provided for.  Conversions involving other
 	   epochs will require use of the appropriate precession,
 	   proper motion, and e-terms routines before and/or
 	   after fk425 is called.

      3)  In the FK4 catalogue the proper motions of stars within
 	   10 degrees of the poles do not embody the differential
 	   e-term effect and should, strictly speaking, be handled
 	   in a different manner from stars outside these regions.
 	   However, given the general lack of homogeneity of the star
 	   data available for routine astrometry, the difficulties of
 	   handling positions that may have been determined from
 	   astrometric fields spanning the polar and non-polar regions,
 	   the likelihood that the differential e-terms effect was not
 	   taken into account when allowing for proper motion in past
 	   astrometry, and the undesirability of a discontinuity in
 	   the algorithm, the decision has been made in this routine to
 	   include the effect of differential e-terms on the proper
 	   motions for all stars, whether polar or not.  At epoch 2000,
 	   and measuring on the sky rather than in terms of dra, the
 	   errors resulting from this simplification are less than
 	   1 milliarcsecond in position and 1 milliarcsecond per
 	   century in proper motion.

   References:

      1  "Mean and apparent place computations in the new IAU System.
          I. The transformation of astrometric catalog systems to the
 	  equinox J2000.0." Smith, C.A.; Kaplan, G.H.; Hughes, J.A.;
	  Seidelmann, P.K.; Yallop, B.D.; Hohenkerk, C.Y.
 	  Astronomical Journal vol. 97, Jan. 1989, p. 265-273.

      2  "Mean and apparent place computations in the new IAU System.
	  II. Transformation of mean star places from FK4 B1950.0 to
 	  FK5 J2000.0 using matrices in 6-space."  Yallop, B.D.;
	  Hohenkerk, C.Y.; Smith, C.A.; Kaplan, G.H.; Hughes, J.A.;
	  Seidelmann, P.K.; Astronomical Journal vol. 97, Jan. 1989,
	  p. 274-279.

      3  "Conversion of positions and proper motions from B1950.0 to the
	  IAU system at J2000.0", Standish, E.M.  Astronomy and
	  Astrophysics, vol. 115, no. 1, Nov. 1982, p. 20-22.

   P.T.Wallace   Starlink   27 October 1987
   Doug Mink     Smithsonian Astrophysical Observatory  7 June 1995 */

{
double	r1950,d1950;	/* B1950.0 ra,dec (rad) */
double	dr1950,dd1950;	/* B1950.0 proper motions (rad/trop.yr) */
double	r2000,d2000;	/* J2000.0 ra,dec (rad) */
double	dr2000,dd2000;	/*J2000.0 proper motions (rad/jul.yr) */

/* Miscellaneous */
double	ur,ud,sr,cr,sd,cd,w,wd;
double	x,y,z,xd,yd,zd, dra,ddec,rcon,scon,tcon;
double	rxysq,rxyzsq,rxy,spxy;
int	i,j;
int	diag = 0;

double	r0[3],r1[3];	/* star position and velocity vectors */
double	v1[6],v2[6];	/* combined position and velocity vectors */

/* Constants */
double	d2pi = 6.283185307179586476925287;	/* two PI */
double	pmf;	/* radians per year to arcsec per century */
double	tiny = 1.e-30; /* small number to avoid arithmetic problems */
double	zero = (double) 0.0;

	pmf = 100 * 60 * 60 * 360 / d2pi;

/* Pick up B1950 data (units radians and arcsec / tc) */
	r1950 = *ra;
	d1950 = *dec;
	dr1950 = *rapm;
	dd1950 = *decpm;
	ur = dr1950 * pmf;
	ud = dd1950 * pmf;

/* Spherical to cartesian */
	sr = sin (r1950);
	cr = cos (r1950);
	sd = sin (d1950);
	cd = cos (d1950);

	r0[0] = cr * cd;
	r0[1] = sr * cd;
	r0[2] =      sd;

	r1[0] = -sr*cd*ur - cr*sd*ud;
	r1[1] =  cr*cd*ur - sr*sd*ud;
	r1[2] = 		cd*ud;

/* Allow for e-terms and express as position + velocity 6-vector */
	w = r0[0] * a[0] + r0[1] * a[1] + r0[2] * a[2];
	wd = r0[0] * ad[0] + r0[1] * ad[1] + r0[2] * ad[2];
	for (i = 0; i < 3; i++) {
	    v1[i] = r0[i] - a[i] + w*r0[i];
	    v1[i+3] = r1[i] - ad[i] + wd*r0[i];
	    }

/* Convert position + velocity vector to Fricke system */
	for (i = 0; i < 6; i++) {
	    w = zero;
	    for (j = 0; j < 6; j++) {
		w = w + em[j][i] * v1[j];
		}
	    v2[i] = w;
	    }

/* Revert to spherical coordinates */
	x = v2[0];
	y = v2[1];
	z = v2[2];
	xd = v2[3];
	yd = v2[4];
	zd = v2[5];

	rxysq = x*x + y*y;
	rxyzsq = rxysq + z*z;
	rxy = sqrt (rxysq);

	spxy = x*xd + y*yd;
	
	if (x == zero && y == zero)
	    r2000 = zero;
	else {
	    r2000 = atan2 (y,x);
	    if (r2000 < zero) 
		r2000 = r2000 + d2pi;
	    }
	d2000 = atan2 (z,rxy);

	if (rxy > tiny) {
	    ur = (x * yd - y * xd) / rxysq;
	    ud = (zd * rxysq - z * spxy) / (rxyzsq * rxy);
	    }
	dr2000 = ur / pmf;
	dd2000 = ud / pmf;

/* Return results */
	*ra = r2000;
	*dec = d2000;
	*rapm = dr2000;
	*decpm = dd2000;

	if (diag) {
	    rcon = 1.74532925199433e-2;
	    scon = 3.6e3 / rcon;
	    tcon = 2.4e2 / rcon;
	    dra = tcon * (r2000 - r1950);
	    ddec = scon * (d2000 - d1950);
	    printf("J2000-B1950: dra= %11.5f sec  ddec= %f11.5f arcsec\n",
		dra, ddec);
	    }
	return;
}

int	idg=0;

/*  l2,b2 system of galactic coordinates
/*  p = 192.25       ra of galactic north pole (mean b1950.0)
/*  q =  62.6        inclination of galactic to mean b1950.0 equator
/*  r =  33          longitude of ascending node
/*  p,q,r are degrees

/*  Equatorial to galactic rotation matrix
    (The Eulerian angles are p, q, 90-r)

	+cp.cq.sr-sp.cr     +sp.cq.sr+cp.cr     -sq.sr
	-cp.cq.cr-sp.sr     -sp.cq.cr+cp.sr     +sq.cr
	cp.sq              +sp.sq              +cq */

static double bgal[3][3] = {
	-0.066988739415,-0.872755765852,-0.483538914632,
	 0.492728466075,-0.450346958020, 0.744584633283,
	-0.867600811151,-0.188374601723, 0.460199784784};

/*---  Transform b1950.0 'fk4' equatorial coordinates to
/*     IAU 1958 galactic coordinates */

fk42gal (dtheta,dphi)

double *dtheta;	/* b1950.0 'fk4' ra in radians
		   Galactic longitude (l2) in radians (returned) */
double *dphi;	/* b1950.0 'fk4' dec in radians
		   Galactic latitude (b2) in radians (returned) */

/*  Note:   The equatorial coordinates are b1950.0 'fk4'.  use the
	    routine jpgalj if conversion from j2000.0 coordinates
	    is required.
	    Reference: blaauw et al, MNRAS,121,123 (1960) */
{
	double pos[3],pos1[3],r,dl,db,rl,rb,rra,rdec,dra,ddec;
	double rcon = 0.0174532925199433;
	void jpcon(),jpcop();
	int i;
	char *eqcoor, *eqstrn();

	dra = *dtheta;
	ddec = *dphi;
	rra = dra * rcon;
	rdec = ddec * rcon;

/*  remove e-terms */
/*	call jpabe (rra,rdec,-1,idg) */

/*  Spherical to Cartesian */
	r = 1.;
	jpcop (rra,rdec,r,pos);

/*  rotate to galactic */
	for (i = 0; i<3; i++) {
	    pos1[i] = pos[0]*bgal[i][0] + pos[1]*bgal[i][1] + pos[2]*bgal[i][2];
	    }

/*  Cartesian to spherical */
        jpcon (pos1,&rl,&rb,&r);

	dl = rl / rcon;
	db = rb / rcon;
	*dtheta = dl;
	*dphi = db;

/*  Print result if in diagnostic mode */
	if (idg) {
	    eqcoor = eqstrn (dra,ddec);
	    printf ("FK42GAL: B1950 RA,Dec= %s\n",eqcoor);
	    printf ("FK42GAL: long = %.5f lat = %.5f\n",dl,db);
	    free (eqcoor);
	    }

	return;
}


/*--- Transform IAU 1958 galactic coordinates to B1950.0 'fk4'
 *    equatorial coordinates */

void gal2fk4 (dtheta,dphi)

double *dtheta;	/* Galactic longitude (l2) in radians
		   B1950 FK4 RA in radians (returned) */
double *dphi;	/* Galactic latitude (b2) in radians
		   B1950 FK4 Dec in radians (returned) */

/*  Note:
       The equatorial coordinates are B1950.0 FK4.  Use the
       routine GAL2FK5 if conversion to J2000 coordinates
       is required.
    Reference:  Blaauw et al, MNRAS,121,123 (1960) */

{
	double pos[3],pos1[3],r,dl,db,rl,rb,rra,rdec,dra,ddec;
	double rcon = 0.0174532925199433;
	void jpcon(),jpcop();
	char *eqcoor, *eqstrn();
	int i;

/*  spherical to cartesian */
	dl = *dtheta;
	db = *dphi;
	rl = dl * rcon;
	rb = db * rcon;
	r = 1.0;
	jpcop (rl,rb,r,pos);

/*  rotate to equatorial coordinates */
	for (i = 0; i < 3; i++) {
	    pos1[i] = pos[0]*bgal[0][i] + pos[1]*bgal[1][i] + pos[2]*bgal[2][i];
	    }

/*  cartesian to spherical */
        jpcon (pos1,&rra,&rdec,&r);

/*  introduce e-terms */
/*	jpabe (rra,rdec,-1,idg); */

	dra = rra / rcon;
	ddec = rdec / rcon;
	*dtheta = dra;
	*dphi = ddec;

/*  print result if in diagnostic mode */
	if (idg) {
	    printf ("GAL2FK4: long = %.5f lat = %.5f\n",dl,db);
	    eqcoor = eqstrn (dra,ddec);
	    printf ("GAL2FK4: B1950 RA,Dec= %s\n",eqcoor);
	    free (eqcoor);
	    }

	return;
}


/*  l2,b2 system of galactic coordinates
    p = 192.25       ra of galactic north pole (mean b1950.0)
    q =  62.6        inclination of galactic to mean b1950.0 equator
    r =  33          longitude of ascending node
    p,q,r are degrees */

/*  Equatorial to galactic rotation matrix
    The eulerian angles are p, q, 90-r
	+cp.cq.sr-sp.cr     +sp.cq.sr+cp.cr     -sq.sr
	-cp.cq.cr-sp.sr     -sp.cq.cr+cp.sr     +sq.cr
	+cp.sq              +sp.sq              +cq		*/

static double jgal[3][3] = {
	-0.054875539726,-0.873437108010,-0.483834985808,
	 0.494109453312,-0.444829589425, 0.746982251810,
	-0.867666135858,-0.198076386122, 0.455983795705};

/* Transform J2000 equatorial coordinates to IAU 1958 galactic coordinates */

void fk52gal (dtheta,dphi)

double *dtheta;	/* J2000 right ascension in degrees
		   Galactic longitude (l2) in degrees (returned) */
double *dphi;	/* J2000 declination in degrees
		   Galactic latitude (b2) in degrees (returned) */

/* Rotation matrices by P.T.Wallace, Starlink eqgal and galeq, March 1986 */
/*  Note:
	The equatorial coordinates are J2000 FK5.  Use the routine
	GAL2FK4 if conversion from B1950 FK4 coordinates is required.
    Reference: Blaauw et al, MNRAS,121,123 (1960) */
{
	double pos[3],pos1[3],r,dl,db,rl,rb,rra,rdec,dra,ddec;
	double rcon = 0.0174532925199433;
	void jpcon(),jpcop();
	char *eqcoor, *eqstrn();
	int i;

/*  Spherical to cartesian */
	dra = *dtheta;
	ddec = *dphi;
	rra = dra * rcon;
	rdec = ddec * rcon;
	r = 1.0;
	(void)jpcop (rra,rdec,r,pos);

/*  Rotate to galactic */
	for (i = 0; i < 3; i++) {
	    pos1[i] = pos[0]*jgal[i][0] + pos[1]*jgal[i][1] + pos[2]*jgal[i][2];
	    }

/*  Cartesian to spherical */
        jpcon (pos1,&rl,&rb,&r);

	dl = rl / rcon;
	db = rb / rcon;
	*dtheta = dl;
	*dphi = db;

/*  Print result if in diagnostic mode */

	if (idg) {
	    eqcoor = eqstrn (dra,ddec);
	    printf ("FK52GAL: J2000 RA,Dec= %s\n",eqcoor);
	    printf ("FK52GAL: long = %.5f lat = %.5f\n",dl,db);
	    free (eqcoor);
	    }

	return;
}


/*--- Transform IAU 1958 galactic coordinates to J2000 equatorial coordinates */

void gal2fk5 (dtheta,dphi)

double *dtheta;	/* Galactic longitude (l2) in degrees
		   J2000.0 ra in degrees (returned) */
double *dphi;	/* Galactic latitude (b2) in degrees
		   J2000.0 dec in degrees (returned) */

/*  Note:
       The equatorial coordinates are J2000.  Use the routine FK42GAL
       if conversion to J2000 coordinates is required.
    Reference: Blaauw et al, MNRAS,121,123 (1960) */

{
	double pos[3],pos1[3],r,dl,db,rl,rb,rra,rdec,dra,ddec;
	double rcon = 0.0174532925199433;
	void jpcon(),jpcop();
	int i;
	char *eqcoor, *eqstrn();

/*  Spherical to Cartesian */
	dl = *dtheta;
	db = *dphi;
	rl = dl * rcon;
	rb = db * rcon;
	r = 1.0;
	jpcop (rl,rb,r,pos);

/*  Rotate to equatorial coordinates */
	for (i = 0; i < 3; i++) {
	    pos1[i] = pos[0]*jgal[0][i] + pos[1]*jgal[1][i] + pos[2]*jgal[2][i];
	    }

/*  Cartesian to Spherical */
        jpcon (pos1,&rra,&rdec,&r);
	dra = rra / rcon;
	ddec = rdec / rcon;
	*dtheta = dra;
	*dphi = ddec;

/*  Print result if in diagnostic mode */
	if (idg) {
	    printf ("GAL2FK5: long = %.5f lat = %.5f\n",dl,db);
	    eqcoor = eqstrn (dra,ddec);
	    printf ("GAL2FK5: J2000 RA,Dec= %s\n",eqcoor);
	    free (eqcoor);
	    }

	return;
}


/* Return string with right ascension in hours and declination in degrees */

char *eqstrn (dra, ddec)

double	dra;		/* Right ascension in degrees */
double	ddec;		/* Declination in degrees */

{
char	*eqcoor;	/* ASCII character string of position (returned) */
char	decp;
int	rah,irm,decd,decm;
double	xpos,ypos,xp,yp,ras,decs;

/*  Right ascension to hours, minutes, and seconds */
	xpos = dra / 15.0;
	rah = (int) xpos;
	xp = (double) 60.0 * (xpos - (double) rah);
	irm = (int) xp;
	ras = (double) 60.0 * (xp - (double) irm);

/* Declination to degrees, minutes, seconds */
	if (ddec < 0) {
	    ypos = -ddec;
	    decp = '-';
	    }
	else {
	    decp = '+';
	    ypos = ddec;
	    }
	decd = (int) ypos;
	yp = (double) 60.0 * (ypos - (double) decd);
	decm = (int) yp;
	decs = (double) 60.0 * (yp - (double) decm);

	eqcoor = malloc (32);
	(void)sprintf (eqcoor,"%02d:%02d:%06.3f %c%02d:%02d:%05.2f",
		       rah,irm,ras,decp,decd,decm,decs);
	if (eqcoor[6] == ' ')
	    eqcoor[6] = '0';
	if (eqcoor[20] == ' ')
	    eqcoor[20] = '0';

	return (eqcoor);
}


/* Convert geocentric equatorial rectangular coordinates to
   right ascension and declination, and distance */

void jpcon (pos,rra,rdec,r)

double pos[3];	/* x,y,z geocentric equatorial position of object */
double *rra;	/* Right ascension in radians */
double *rdec;	/* Declination in radians */
double *r;	/* Distance to object in same units as pos */

{
	double x,y,z,rxy,rxy2,z2;

	x = pos[0];
	y = pos[1];
	z = pos[2];

	*rra = atan2 (y, x);
	if (*rra < 0.) *rra = *rra + 6.283185307179586;

	rxy2 = x*x + y*y;
	rxy = sqrt (rxy2);
	*rdec = atan2 (z, rxy);

	z2 = z * z;
	*r = sqrt (rxy2 + z2);

	return;
}


/* Convert right ascension, declination, and distance to
   geocentric equatorial rectangular coordinates */

void jpcop (rra,rdec,r,pos)

double rra;	/* Right ascension in radians */
double rdec;	/* Declination in radians */
double r;	/* Distance to object in same units as pos */
double pos[3];	/* x,y,z geocentric equatorial position of object */
{
	pos[0] = r * cos (rra) * cos (rdec);
	pos[1] = r * sin (rra) * cos (rdec);
	pos[2] = r * sin (rdec);

	return;
}

/* The following routines are almost verbatim from Patrick Wallace's SLALIB */

void fk4prec (ep0, ep1, ra, dec)

double ep0;	/* Starting Besselian epoch */
double ep1;	/* Ending Besselian epoch */
double *ra;	/* RA in radians mean equator & equinox of epoch ep0
		      mean equator & equinox of epoch ep1 (returned) */
double *dec;	/* Dec in radians mean equator & equinox of epoch ep0
		       mean equator & equinox of epoch ep1 (returned) */
/*
**  slaPreces:
**  Precession - FK4 (Bessel-Newcomb, pre-IAU1976)
**
**  Note:
**      This routine will not correctly convert between the old and
**      the new systems - for example conversion from B1950 to J2000.
**      For these purposes see fk425, fk524, fk45m and fk54m.
**
**  P.T.Wallace   Starlink   22 December 1993
*/
{
    double pm[3][3], v1[3], v2[3];
    void mprecfk4(), slaDcs2c(), slaDmxv(), slaDcc2s();
    double slaDranrm();
 
    /* Generate appropriate precession matrix */
    mprecfk4 ( ep0, ep1, pm );
 
    /* Convert RA,Dec to x,y,z */
    slaDcs2c ( *ra, *dec, v1 );
 
    /* Precess */
    slaDmxv ( pm, v1, v2 );
 
    /* Back to RA,Dec */
    slaDcc2s ( v2, ra, dec );
    *ra = slaDranrm ( *ra );
}

void fk5prec (ep0, ep1, ra, dec)

double ep0;	/* Starting epoch */
double ep1;	/* Ending epoch */
double *ra;	/* RA in radians mean equator & equinox of epoch ep0
		      mean equator & equinox of epoch ep1 (returned) */
double *dec;	/* Dec in radians mean equator & equinox of epoch ep0
		       mean equator & equinox of epoch ep1 (returned) */
/*
**  slaPreces:
**  Precession -  FK5 (Fricke, post-IAU1976)
**
**  Note:
**      This routine will not correctly convert between the old and
**      the new systems - for example conversion from B1950 to J2000.
**      For these purposes see fk425, fk524, fk45m and fk54m.
**
**  P.T.Wallace   Starlink   22 December 1993
*/
{
    double pm[3][3], v1[3], v2[3];
    void mprecfk5(), slaDcs2c(), slaDmxv(), slaDcc2s();
    double slaDranrm();
 
    /* Generate appropriate precession matrix */
    mprecfk5 ( ep0, ep1, pm );
 
    /* Convert RA,Dec to x,y,z */
    slaDcs2c ( *ra, *dec, v1 );
 
    /* Precess */
    slaDmxv ( pm, v1, v2 );
 
    /* Back to RA,Dec */
    slaDcc2s ( v2, ra, dec );
    *ra = slaDranrm ( *ra );
}


void slaDcs2c (a, b, v)

double a;	/* Right ascension in radians */
double b;	/* Declination in radians */
double *v;	/* x,y,z unit vector (returned) */

/*
**  slaDcs2c: Spherical coordinates to direction cosines.
**
**  The spherical coordinates are longitude (+ve anticlockwise
**  looking from the +ve latitude pole) and latitude.  The
**  Cartesian coordinates are right handed, with the x axis
**  at zero longitude and latitude, and the z axis at the
**  +ve latitude pole.
**
**  P.T.Wallace   Starlink   31 October 1993
*/
{
    double cosb;
 
    cosb = cos ( b );
    v[0] = cos ( a ) * cosb;
    v[1] = sin ( a ) * cosb;
    v[2] = sin ( b );
}


void slaDmxv (dm, va, vb)

double (*dm)[3];	/* 3x3 Matrix */
double *va;		/* Vector */
double *vb;		/* Result vector (returned) */

/*
**  slaDmxv:
**  Performs the 3-d forward unitary transformation:
**     vector vb = matrix dm * vector va
**
**  P.T.Wallace   Starlink   31 October 1993
*/
{
    int i, j;
    double w, vw[3];
 
    /* Matrix dm * vector va -> vector vw */
    for ( j = 0; j < 3; j++ ) {
	w = 0.0;
	for ( i = 0; i < 3; i++ ) {
	    w += dm[j][i] * va[i];
	    }
	vw[j] = w;
	}
 
    /* Vector vw -> vector vb */
    for ( j = 0; j < 3; j++ ) {
	vb[j] = vw[j];
	}
}


void slaDcc2s (v, a, b)

double *v;	/* x,y,z vector */
double *a;	/* Right ascension in radians */
double *b;	/* Declination in radians */

/*
**  slaDcc2s:
**  Direction cosines to spherical coordinates.
**
**  Returned:
**     *a,*b  double      spherical coordinates in radians
**
**  The spherical coordinates are longitude (+ve anticlockwise
**  looking from the +ve latitude pole) and latitude.  The
**  Cartesian coordinates are right handed, with the x axis
**  at zero longitude and latitude, and the z axis at the
**  +ve latitude pole.
**
**  If v is null, zero a and b are returned.
**  At either pole, zero a is returned.
**
**  P.T.Wallace   Starlink   31 October 1993
*/
{
    double x, y, z, r;
 
    x = v[0];
    y = v[1];
    z = v[2];
    r = sqrt ( x * x + y * y );
 
    *a = ( r != 0.0 ) ? atan2 ( y, x ) : 0.0;
    *b = ( z != 0.0 ) ? atan2 ( z, r ) : 0.0;
}

 
/* 2pi */
#define D2PI 6.2831853071795864769252867665590057683943387987502

double slaDranrm (angle)

double angle;	/* angle in radians */

/*
**  slaDranrm:
**  Normalize angle into range 0-2 pi.
**  The result is angle expressed in the range 0-2 pi (double).
**  Defined in slamac.h:  D2PI
**
**  P.T.Wallace   Starlink   30 October 1993
*/
{
    double w;
 
    w = fmod ( angle, D2PI );
    return ( w >= 0.0 ) ? w : w + D2PI;
}


/* pi/(180*3600):  arcseconds to radians */
#define DAS2R 4.8481368110953599358991410235794797595635330237270e-6

void mprecfk4 (bep0, bep1, rmatp)

double bep0;		/* Beginning Besselian epoch */
double bep1;		/* Ending Besselian epoch */
double (*rmatp)[3];	/* 3x3 Precession matrix (returned) */

/*
**  slaPrebn:
**  Generate the matrix of precession between two epochs,
**  using the old, pre-IAU1976, Bessel-Newcomb model, using
**  Kinoshita's formulation (double precision)
**
**  The matrix is in the sense   v(bep1)  =  rmatp * v(bep0)
**
**  Reference:
**     Kinoshita, H. (1975) 'Formulas for precession', SAO Special
**     Report No. 364, Smithsonian Institution Astrophysical
**     Observatory, Cambridge, Massachusetts.
**
**  P.T.Wallace   Starlink   30 October 1993
*/
{
    double bigt, t, tas2r, w, zeta, z, theta;
    void slaDeuler();
 
    /* Interval between basic epoch B1850.0 and beginning epoch in TC */
    bigt  = ( bep0 - 1850.0 ) / 100.0;
 
    /* Interval over which precession required, in tropical centuries */
    t = ( bep1 - bep0 ) / 100.0;
 
    /* Euler angles */
    tas2r = t * DAS2R;
    w = 2303.5548 + ( 1.39720 + 0.000059 * bigt ) * bigt;
    zeta = (w + ( 0.30242 - 0.000269 * bigt + 0.017996 * t ) * t ) * tas2r;
    z = (w + ( 1.09478 + 0.000387 * bigt + 0.018324 * t ) * t ) * tas2r;
    theta = ( 2005.1125 + ( - 0.85294 - 0.000365* bigt ) * bigt +
	    ( - 0.42647 - 0.000365 * bigt - 0.041802 * t ) * t ) * tas2r;
 
    /* Rotation matrix */
    slaDeuler ( "ZYZ", -zeta, theta, -z, rmatp );
}


void mprecfk5 (ep0, ep1, rmatp)

double ep0;		/* Beginning epoch */
double ep1;		/* Ending epoch */
double (*rmatp)[3];	/* 3x3 Precession matrix (returned) */

/*
**  slaPrec:
**  Form the matrix of precession between two epochs (IAU 1976, FK5).
**  Notes:
**  1)  The epochs are TDB (loosely ET) Julian epochs.
**  2)  The matrix is in the sense   v(ep1)  =  rmatp * v(ep0) .
**
**  References:
**     Lieske,J.H., 1979. Astron. Astrophys.,73,282.
**          equations (6) & (7), p283.
**     Kaplan,G.H., 1981. USNO circular no. 163, pa2.
**
**  P.T.Wallace   Starlink   31 October 1993
*/
{
    double t0, t, tas2r, w, zeta, z, theta;
    void slaDeuler();
 
    /* Interval between basic epoch J2000.0 and beginning epoch (JC) */
    t0 = ( ep0 - 2000.0 ) / 100.0;
 
    /* Interval over which precession required (JC) */
    t =  ( ep1 - ep0 ) / 100.0;
 
    /* Euler angles */
    tas2r = t * DAS2R;
    w = 2306.2181 + ( ( 1.39656 - ( 0.000139 * t0 ) ) * t0 );
    zeta = (w + ( ( 0.30188 - 0.000344 * t0 ) + 0.017998 * t ) * t ) * tas2r;
    z = (w + ( ( 1.09468 + 0.000066 * t0 ) + 0.018203 * t ) * t ) * tas2r;
    theta = ( ( 2004.3109 + ( - 0.85330 - 0.000217 * t0 ) * t0 )
	  + ( ( -0.42665 - 0.000217 * t0 ) - 0.041833 * t ) * t ) * tas2r;
 
    /* Rotation matrix */
    slaDeuler ( "ZYZ", -zeta, theta, -z, rmatp );
}


void slaDeuler (order, phi, theta, psi, rmat)

char *order;		/* specifies about which axes the rotations occur */
double phi;		/* 1st rotation (radians) */
double theta;		/* 2nd rotation (radians) */
double psi;		/* 3rd rotation (radians) */
double (*rmat)[3];	/* 3x3 Rotation matrix (returned) */

/*
**  slaDeuler:
**  Form a rotation matrix from the Euler angles - three successive
**  rotations about specified Cartesian axes.
**
**  A rotation is positive when the reference frame rotates
**  anticlockwise as seen looking towards the origin from the
**  positive region of the specified axis.
**
**  The characters of order define which axes the three successive
**  rotations are about.  A typical value is 'zxz', indicating that
**  rmat is to become the direction cosine matrix corresponding to
**  rotations of the reference frame through phi radians about the
**  old z-axis, followed by theta radians about the resulting x-axis,
**  then psi radians about the resulting z-axis.
**
**  The axis names can be any of the following, in any order or
**  combination:  x, y, z, uppercase or lowercase, 1, 2, 3.  Normal
**  axis labelling/numbering conventions apply;  the xyz (=123)
**  triad is right-handed.  Thus, the 'zxz' example given above
**  could be written 'zxz' or '313' (or even 'zxz' or '3xz').  Order
**  is terminated by length or by the first unrecognised character.
**
**  Fewer than three rotations are acceptable, in which case the later
**  angle arguments are ignored.  Zero rotations produces a unit rmat.
**
**  P.T.Wallace   Starlink   17 November 1993
*/
{
   int j, i, l, n, k;
   double result[3][3], rotn[3][3], angle, s, c , w, wm[3][3];
   char axis;
 
/* Initialize result matrix */
   for ( j = 0; j < 3; j++ ) {
      for ( i = 0; i < 3; i++ ) {
         result[i][j] = ( i == j ) ? 1.0 : 0.0;
      }
   }
 
/* Establish length of axis string */
   l = strlen ( order );
 
/* Look at each character of axis string until finished */
   for ( n = 0; n < 3; n++ ) {
      if ( n <= l ) {
 
      /* Initialize rotation matrix for the current rotation */
         for ( j = 0; j < 3; j++ ) {
            for ( i = 0; i < 3; i++ ) {
               rotn[i][j] = ( i == j ) ? 1.0 : 0.0;
            }
         }
 
      /* Pick up the appropriate Euler angle and take sine & cosine */
         switch ( n ) {
         case 0 :
           angle = phi;
           break;
         case 1 :
           angle = theta;
           break;
         case 2 :
           angle = psi;
           break;
         }
         s = sin ( angle );
         c = cos ( angle );
 
      /* Identify the axis */
         axis =  order[n];
         if ( ( axis == 'X' ) || ( axis == 'x' ) || ( axis == '1' ) ) {
 
         /* Matrix for x-rotation */
            rotn[1][1] = c;
            rotn[1][2] = s;
            rotn[2][1] = -s;
            rotn[2][2] = c;
         }
         else if ( ( axis == 'Y' ) || ( axis == 'y' ) || ( axis == '2' ) ) {
 
         /* Matrix for y-rotation */
            rotn[0][0] = c;
            rotn[0][2] = -s;
            rotn[2][0] = s;
            rotn[2][2] = c;
         }
         else if ( ( axis == 'Z' ) || ( axis == 'z' ) || ( axis == '3' ) ) {
 
         /* Matrix for z-rotation */
            rotn[0][0] = c;
            rotn[0][1] = s;
            rotn[1][0] = -s;
            rotn[1][1] = c;
         } else {
 
         /* Unrecognized character - fake end of string */
            l = 0;
         }
 
      /* Apply the current rotation (matrix rotn x matrix result) */
         for ( i = 0; i < 3; i++ ) {
            for ( j = 0; j < 3; j++ ) {
               w = 0.0;
               for ( k = 0; k < 3; k++ ) {
                  w += rotn[i][k] * result[k][j];
               }
               wm[i][j] = w;
            }
         }
         for ( j = 0; j < 3; j++ ) {
            for ( i= 0; i < 3; i++ ) {
               result[i][j] = wm[i][j];
            }
         }
      }
   }
 
/* Copy the result */
   for ( j = 0; j < 3; j++ ) {
      for ( i = 0; i < 3; i++ ) {
         rmat[i][j] = result[i][j];
      }
   }
}
/*
 * Nov  6 1995	Include stdlib.h instead of malloc.h
 * Apr  1 1996	Add arbitrary epoch precession
 */
