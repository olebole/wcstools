#include "slalib.h"
#include "slamac.h"
void slaFk52h ( double r5, double d5, double dr5, double dd5,
                double *rh, double *dh, double *drh, double *ddh )
/*
**  - - - - - - - - -
**   t p t F k 5 2 h
**  - - - - - - - - -
**
**  Transform FK5 (J2000) star data into the Hipparcos frame.
**
**  (double precision)
**
**  This routine transforms FK5 star positions and proper motions
**  into the frame of the Hipparcos catalogue.
**
**  Given (all FK5, equinox J2000, epoch J2000):
**     r5      double    RA (radians)
**     d5      double    Dec (radians)
**     dr5     double    proper motion in RA (dRA/dt, rad/Jyear)
**     dd5     double    proper motion in Dec (dDec/dt, rad/Jyear)
**
**  Returned (all Hipparcos, epoch J2000):
**     rh      double    RA (radians)
**     dh      double    Dec (radians)
**     drh     double    proper motion in RA (dRA/dt, rad/Jyear)
**     ddh     double    proper motion in Dec (dDec/dt, rad/Jyear)
**
**  Called:  slaDs2c6, slaDav2m, slaDmxv, slaDvxv, slaDc62s
**
**  Notes:
**
**  1)  The proper motions in RA are dRA/dt rather than
**      cos(Dec)*dRA/dt, and are per year rather than per century.
**
**  2)  The FK5 to Hipparcos transformation consists of a pure
**      rotation and spin;  zonal errors in the FK5 catalogue are
**      not taken into account.
**
**  3)  The published orientation and spin components are interpreted
**      as "axial vectors".  An axial vector points at the pole of the
**      rotation and its length is the amount of rotation in radians.
**
**  4)  See also slaH2fk5, slaFk5hz, slaHfk5z.
**
**  Reference:
**
**     M.Feissel & F.Mignard, Astron. Astrophys. 331, L33-L36 (1998).
**
**  Last revision:   8 October 1998
**
**  Copyright P.T.Wallace.  All rights reserved.
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

   double pv5[6], r5h[3][3], vv[3], pvh[6], r, v;
   int i;


/* FK5 barycentric position/velocity 6-vector (normalized). */
   slaDs2c6 ( r5, d5, 1.0, dr5, dd5, 0.0, pv5 );

/* FK5 to Hipparcos orientation matrix. */
   slaDav2m ( ortn, r5h );

/* Rotate & spin the 6-vector into the Hipparcos frame. */
   slaDmxv ( r5h, pv5, pvh );
   slaDvxv ( pv5, s5, vv );
   for ( i = 0; i < 3; i++ ) {
      vv [ i ] = pv5 [ i + 3 ] + vv [ i ];
   }
   slaDmxv ( r5h, vv, pvh + 3 );

/* Hipparcos 6-vector to spherical. */
   slaDc62s ( pvh, rh, dh, &r, drh, ddh, &v );
}
