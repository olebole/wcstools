#include "slalib.h"
#include "slamac.h"
void slaH2fk5 ( double rh, double dh, double drh, double ddh,
                double *r5, double *d5, double *dr5, double *dd5 )
/*
**  - - - - - - - - -
**   t p t H 2 f k 5
**  - - - - - - - - -
**
**  Transform Hipparcos star data into the FK5 (J2000) system.
**
**  (double precision)
**
**  This routine transforms Hipparcos star positions and proper
**  motions into FK5 J2000.
**
**  Given (all Hipparcos, epoch J2000):
**     rh      double    RA (radians)
**     dh      double    Dec (radians)
**     drh     double    proper motion in RA (dRA/dt, rad/Jyear)
**     ddh     double    proper motion in Dec (dDec/dt, rad/Jyear)
**
**  Returned (all FK5, equinox J2000, Epoch J2000):
**     r5      double    RA (radians)
**     d5      double    Dec (radians)
**     dr5     double    proper motion in RA (dRA/dt, rad/Jyear)
**     dd5     double    proper motion in Dec (dDec/dt, rad/Jyear)
**
**  Called:  slaDs2c6, slaDav2m, slaDmxv, slaDimxv, slaDvxv, slaDc62s
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
**  4)  See also sla_FK52H, sla_FK5HZ, sla_HFK5Z.
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

   double pvh[6], r5h[3][3], sh[3], vv[3], pv5[6], r, v;
   int i;


/* Hipparcos barycentric position/velocity 6-vector (normalized). */
   slaDs2c6 ( rh, dh, 1.0, drh, ddh, 0.0, pvh );

/* FK5 to Hipparcos orientation matrix. */
   slaDav2m ( ortn, r5h );

/* Rotate the spin vector into the Hipparcos frame. */
   slaDmxv ( r5h, s5, sh );

/* De-orient & de-spin the 6-vector into FK5 J2000. */
   slaDimxv ( r5h, pvh, pv5 );
   slaDvxv ( pvh, sh, vv );
   for ( i = 0; i < 3; i++ ) {
      vv [ i ] = pvh [ i + 3 ] - vv [ i ];
   }
   slaDimxv ( r5h, vv, pv5 + 3 );

/* FK5 6-vector to spherical. */
   slaDc62s ( pv5, r5, d5, &r, dr5, dd5, &v );
}
