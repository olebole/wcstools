#include "slalib.h"
#include "slamac.h"
void slaHfk5z ( double rh, double dh, double epoch,
                double *r5, double *d5, double *dr5, double *dd5 )
/*
**  - - - - - - - - -
**   t p t H f k 5 z
**  - - - - - - - - -
**
**  Transform a Hipparcos star position into FK5 J2000, assuming
**  zero Hipparcos proper motion.
**
**  (double precision)
**
**  Given:
**     rh      double    Hipparcos RA (radians)
**     dh      double    Hipparcos Dec (radians)
**     epoch   double    Julian epoch (TDB)
**
**  Returned (all FK5, equinox J2000, epoch EPOCH):
**     r5      double    RA (radians)
**     d5      double    Dec (radians)
**
**  Called:  slaDcs2c, slaDav2m, slaDmxv, slaDav2m, slaDmxm,
**           slaDimxv, slaDvxv, slaDc62s
**
**  Notes:
**
**  1)  The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.
**
**  2)  The FK5 to Hipparcos transformation consists of a pure
**      rotation and spin;  zonal errors in the FK5 catalogue are
**      not taken into account.
**
**  3)  The published orientation and spin components are interpreted
**      as "axial vectors".  An axial vector points at the pole of the
**      rotation and its length is the amount of rotation in radians.
**
**  4)  It was the intention that Hipparcos should be a close
**      approximation to an inertial frame, so that distant objects
**      have zero proper motion;  such objects have (in general)
**      non-zero proper motion in FK5, and this routine returns those
**      fictitious proper motions.
**
**  5)  The position returned by this routine is in the FK5 J2000
**      reference frame but at the specified epoch.
**
**  6)  See also slaFk52h, slaH2fk5, slaFk5zhz.
**
**  Reference:
**
**     M.Feissel & F.Mignard, Astron. Astrophys. 331, L33-L36 (1998).
**
**  Last revision:   8 October 1998
**
**  Copyright P.T.Wallace.  All rights reserved.
**
**  Copyright (C) 1998 Rutherford Appleton Laboratory
*/

#define AS2R 0.484813681109535994e-5    /* arcseconds to radians */

{
/* FK5 to Hipparcos orientation and spin (radians, radians/year) */
   static double ortn[3] = { -19.9e-3 * AS2R,
                              -9.1e-3 * AS2R,
                              22.9e-3 * AS2R },
                   s5[3] = { -0.30e-3 * AS2R,
                              0.60e-3 * AS2R,
                              0.70e-3 * AS2R };

   double ph[6], r5h[3][3], sh[3], t, vst[3], rst[3][3], r5ht[3][3],
          pv5e[6], vv[3], r, v;
   int i;


/* Hipparcos barycentric position vector (normalized). */
   slaDcs2c ( rh, dh, ph );

/* FK5 to Hipparcos orientation matrix. */
   slaDav2m ( ortn, r5h );

/* Rotate the spin vector into the Hipparcos frame. */
   slaDmxv ( r5h, s5, sh );

/* Time interval from J2000 to epoch. */
   t = epoch - 2000.0;

/* Axial vector:  accumulated Hipparcos wrt FK5 spin over that interval. */
   for ( i = 0; i < 3; i++ ) {
      vst [ i ] = s5 [ i ] * t;
   }

/* Express the accumulated spin as a rotation matrix. */
   slaDav2m ( vst, rst );

/* Rotation matrix:  accumulated spin, then FK5 to Hipparcos. */
   slaDmxm ( r5h, rst, r5ht );

/* De-orient & de-spin the vector into FK5 J2000 at epoch. */
   slaDimxv ( r5ht, ph, pv5e );
   slaDvxv ( sh, ph, vv );
   slaDimxv ( r5ht, vv, pv5e + 3 );

/* FK5 position/velocity 6-vector to spherical. */
   slaDc62s( pv5e, r5, d5, &r, dr5, dd5, &v );
}
