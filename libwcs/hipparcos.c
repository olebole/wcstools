/* Convert from FK5 to Hipparcos from if no proper motion */

void
fk52h (ra,dec)

double	*ra;		/* Right ascension in degrees (J2000 in, B1950 out) */
double	*dec;		/* Declination in degrees (J2000 in, B1950 out) */
{
    double *rapm;	/* Proper motion in right ascension */
    double *decpm;	/* Proper motion in declination */

    *rapm = 0.0;
    *decpm = 0.0;

    fk52hm (ra, dec, rapm, decpm);
}


void
fk52hm (ra,dec,rapm,decpm)

double	*ra;		/* Right ascension in degrees (J2000 in, B1950 out) */
double	*dec;		/* Declination in degrees (J2000 in, B1950 out) */
double	*rapm;		/* Proper motion in right ascension */
double	*decpm;		/* Proper motion in declination */
{
    double theta, phi, ct, st, dt, cp, sp, dp;
    double v[6];
    int diag = 0;

    /* convert to vector */
    s2v6 (*ra, *dec, *rapm, *decpm, v);

    /* Compute FK5 to Hipparcos orientation matrix */

    /* Rotate and spin vector into Hipparcos frame */

    /* Convert back to spherical coordinates */
    v2s6 (v, ra, dec, rapm, decpm);

    if (diag) {
	double scon, tcon, dra, ddec;
	scon = raddeg (3.6e3);
	tcon = raddeg (2.4e2);
	dra = tcon * (r1950 - r2000);
	ddec = scon * (d1950 - d2000);
	printf("B1950-J2000: dra= %11.5f sec  ddec= %f11.5f arcsec\n",
		dra, ddec);
	}

}


/* Convert spherical coordinates to vector */

void
s2v6 (ra, dec, rapm, decpm, v)

double ra;	/* Right ascension in degrees */
double dec;	/* Declination in degrees */
double rapm;	/* Right ascension proper motion in degrees/year */
double decpm;	/* Declination in proper motion in degrees/year */
double *v;	/* 6-vector of position and velocity (returned) */
{
    double theta, phi, ct, st, dt, cp, sp, dp;

    /* convert to vector */
    theta = degrad(ra);
    ct = cos (theta);
    st = sin (theta);
    dt = degrad (rapm);
    phi = degrad(dec);
    cp = cos (phi);
    sp = sin (phi);
    dp = degrad (decpm);
    v[0] = ct * cp;
    v[1] = st * cp;
    v[2] = sp;
    v[3] = (-v[1] * dt) - (dp * sp * ct);
    v[4] = (v[0] * dt) - (dp * sp * st);
    v[5] = dp * cp;

    return;
}


/* Convert vector to spherical coordinates */

void
v2s6 (v, ra, dec, rapm, decpm)

double *v;	/* 6-vector of position and velocity */
double *ra;	/* Right ascension in degrees (returned) */
double *dec;	/* Declination in degrees (returned) */
double *rapm;	/* Right ascension proper motion in degrees/year (returned) */
double *decpm;	/* Declination in proper motion in degrees/year (returned) */
{
    double rxy2, rxy, r2, theta, phi, dt, dp;

    rxy2 = v[0]*v[0] + [v[1]*v[1];
    rxy = sqrt (rxy2);
    r2 = rxy2 + v[2] * v[2]

    theta = atan2 (v[1], v[0]);
    phi = atan2 (v[2], rxy);
    dt = (v[0] * v[3] - v[1] * v[4]) / rxy2;
    dp = (v[5] * RXY2 - v[2] * (v[0] * v[3] + v[1] * v[4]) / (r2 * rxy);

    /* Convert results to degrees */
    *ra = raddeg (theta);
    *dec = raddeg (phi);
    *rapm = raddeg (dt);
    *decpm = raddeg (dp);

    return;
}
