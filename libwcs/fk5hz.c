#include "slalib.h"
#include "slamac.h"
void slaFk5hz ( double r5, double d5, double epoch, double *rh, double *dh )
/*
**  - - - - - - - - -
**   t p t F k 5 h z
**  - - - - - - - - -
**
**  Transform an FK5 (J2000) star position into the frame of the
**  Hipparcos catalogue, assuming zero Hipparcos proper motion.
**
**  (double precision)
**
**  This routine converts a star position from the FK5 system to
**  the Hipparcos system, in such a way that the Hipparcos proper
**  motion is zero.  Because such a star has, in general, a non-zero
**  proper motion in the FK5 system, the routine requires the epoch
**  at which the position in the FK5 system was determined.
**
**  Given:
**     r5      double    FK5 RA (radians), equinox J2000, epoch EPOCH
**     d5      double    FK5 Dec (radians), equinox J2000, epoch EPOCH
**     epoch   double    Julian epoch (TDB)
**
**  Returned (all Hipparcos):
**     rh      double    RA (radians)
**     dh      double    Dec (radians)
**
**  Called:  slaDcs2c, slaDav2m, slaDimxv, slaDmxv, slaDcc2s
**
**  Notes:
**
**  1)  The FK5 to Hipparcos transformation consists of a pure
**      rotation and spin;  zonal errors in the FK5 catalogue are
**      not taken into account.
**
**  2)  The published orientation and spin components are interpreted
**      as "axial vectors".  An axial vector points at the pole of the
**      rotation and its length is the amount of rotation in radians.
**
**  3)  See also sla_FK52H, sla_H2FK5, sla_HFK5Z.
**
**  Reference:
**
**     M.Feissel & F.Mignard, Astron. Astrophys. 331, L33-L36 (1998).
**
**  P.T.Wallace   Starlink   8 October 1998
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

   double p5e[3], r5h[3][3], t, vst[3], rst[3][3], p5[3], ph[3];
   int i;


/* FK5 barycentric position vector. */
   slaDcs2c ( r5, d5, p5e );

/* FK5 to Hipparcos orientation matrix. */
   slaDav2m ( ortn, r5h );

/* Time interval from epoch to J2000. */
   t = 2000.0 - epoch;

/* Axial vector:  accumulated Hipparcos wrt FK5 spin over that interval. */
   for ( i = 0; i < 3; i++ ) {
      vst [ i ] = s5 [ i ] * t;
   }

/* Express the accumulated spin as a rotation matrix. */
   slaDav2m ( vst, rst );

/* Derotate the vector's FK5 axes back to epoch. */
   slaDimxv ( rst, p5e, p5 );

/* Rotate the vector into the Hipparcos frame. */
   slaDmxv ( r5h, p5, ph );

/* Hipparcos vector to spherical. */
   slaDcc2s ( ph, rh, dh );
}
