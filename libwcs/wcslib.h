#ifndef wcslib_h_
#define wcslib_h_

#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(__STDC__) && !defined(__cplusplus)
#ifndef const
#define const
#endif
#endif

#ifdef TRIGD
#include <TRIGD>
#else /* not TRIGD */
#if __STDC__ || defined(__cplusplus)
   double cosdeg(const double);
   double sindeg(const double);
   double tandeg(const double);
   double acosdeg(const double);
   double asindeg(const double);
   double atandeg(const double);
   double atan2deg(const double, const double);
#else
   double cosdeg();
   double sindeg();
   double tandeg();
   double acosdeg();
   double asindeg();
   double atandeg();
   double atan2deg();
#endif

/* Domain tolerance for asin and acos functions. */
#define WCSTRIG_TOL 1e-10
#endif /* TRIGD */

#ifdef __cplusplus
};
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct prjprm {
   int flag;
   int n;
   double r0;
   double p[10];
   double w[10];
};

#if __STDC__ || defined(__cplusplus)
   int azpset(struct prjprm *);
   int azpfwd(const double, const double, struct prjprm *, double *, double *);
   int azprev(const double, const double, struct prjprm *, double *, double *);
   int tanset(struct prjprm *);
   int tanfwd(const double, const double, struct prjprm *, double *, double *);
   int tanrev(const double, const double, struct prjprm *, double *, double *);
   int sinset(struct prjprm *);
   int sinfwd(const double, const double, struct prjprm *, double *, double *);
   int sinrev(const double, const double, struct prjprm *, double *, double *);
   int stgset(struct prjprm *);
   int stgfwd(const double, const double, struct prjprm *, double *, double *);
   int stgrev(const double, const double, struct prjprm *, double *, double *);
   int arcset(struct prjprm *);
   int arcfwd(const double, const double, struct prjprm *, double *, double *);
   int arcrev(const double, const double, struct prjprm *, double *, double *);
   int zpnset(struct prjprm *);
   int zpnfwd(const double, const double, struct prjprm *, double *, double *);
   int zpnrev(const double, const double, struct prjprm *, double *, double *);
   int zeaset(struct prjprm *);
   int zeafwd(const double, const double, struct prjprm *, double *, double *);
   int zearev(const double, const double, struct prjprm *, double *, double *);
   int airset(struct prjprm *);
   int airfwd(const double, const double, struct prjprm *, double *, double *);
   int airrev(const double, const double, struct prjprm *, double *, double *);
   int cypset(struct prjprm *);
   int cypfwd(const double, const double, struct prjprm *, double *, double *);
   int cyprev(const double, const double, struct prjprm *, double *, double *);
   int carset(struct prjprm *);
   int carfwd(const double, const double, struct prjprm *, double *, double *);
   int carrev(const double, const double, struct prjprm *, double *, double *);
   int merset(struct prjprm *);
   int merfwd(const double, const double, struct prjprm *, double *, double *);
   int merrev(const double, const double, struct prjprm *, double *, double *);
   int ceaset(struct prjprm *);
   int ceafwd(const double, const double, struct prjprm *, double *, double *);
   int cearev(const double, const double, struct prjprm *, double *, double *);
   int copset(struct prjprm *);
   int copfwd(const double, const double, struct prjprm *, double *, double *);
   int coprev(const double, const double, struct prjprm *, double *, double *);
   int codset(struct prjprm *);
   int codfwd(const double, const double, struct prjprm *, double *, double *);
   int codrev(const double, const double, struct prjprm *, double *, double *);
   int coeset(struct prjprm *);
   int coefwd(const double, const double, struct prjprm *, double *, double *);
   int coerev(const double, const double, struct prjprm *, double *, double *);
   int cooset(struct prjprm *);
   int coofwd(const double, const double, struct prjprm *, double *, double *);
   int coorev(const double, const double, struct prjprm *, double *, double *);
   int bonset(struct prjprm *);
   int bonfwd(const double, const double, struct prjprm *, double *, double *);
   int bonrev(const double, const double, struct prjprm *, double *, double *);
   int pcoset(struct prjprm *);
   int pcofwd(const double, const double, struct prjprm *, double *, double *);
   int pcorev(const double, const double, struct prjprm *, double *, double *);
   int glsset(struct prjprm *);
   int glsfwd(const double, const double, struct prjprm *, double *, double *);
   int glsrev(const double, const double, struct prjprm *, double *, double *);
   int parset(struct prjprm *);
   int parfwd(const double, const double, struct prjprm *, double *, double *);
   int parrev(const double, const double, struct prjprm *, double *, double *);
   int aitset(struct prjprm *);
   int aitfwd(const double, const double, struct prjprm *, double *, double *);
   int aitrev(const double, const double, struct prjprm *, double *, double *);
   int molset(struct prjprm *);
   int molfwd(const double, const double, struct prjprm *, double *, double *);
   int molrev(const double, const double, struct prjprm *, double *, double *);
   int cscset(struct prjprm *);
   int cscfwd(const double, const double, struct prjprm *, double *, double *);
   int cscrev(const double, const double, struct prjprm *, double *, double *);
   int qscset(struct prjprm *);
   int qscfwd(const double, const double, struct prjprm *, double *, double *);
   int qscrev(const double, const double, struct prjprm *, double *, double *);
   int tscset(struct prjprm *);
   int tscfwd(const double, const double, struct prjprm *, double *, double *);
   int tscrev(const double, const double, struct prjprm *, double *, double *);
#else
   int azpset(), azpfwd(), azprev();
   int tanset(), tanfwd(), tanrev();
   int sinset(), sinfwd(), sinrev();
   int stgset(), stgfwd(), stgrev();
   int arcset(), arcfwd(), arcrev();
   int zpnset(), zpnfwd(), zpnrev();
   int zeaset(), zeafwd(), zearev();
   int airset(), airfwd(), airrev();
   int cypset(), cypfwd(), cyprev();
   int carset(), carfwd(), carrev();
   int merset(), merfwd(), merrev();
   int ceaset(), ceafwd(), cearev();
   int copset(), copfwd(), coprev();
   int codset(), codfwd(), codrev();
   int coeset(), coefwd(), coerev();
   int cooset(), coofwd(), coorev();
   int bonset(), bonfwd(), bonrev();
   int pcoset(), pcofwd(), pcorev();
   int glsset(), glsfwd(), glsrev();
   int parset(), parfwd(), parrev();
   int aitset(), aitfwd(), aitrev();
   int molset(), molfwd(), molrev();
   int cscset(), cscfwd(), cscrev();
   int qscset(), qscfwd(), qscrev();
   int tscset(), tscfwd(), tscrev();
#endif

extern const char *prjset_errmsg[];
extern const char *prjfwd_errmsg[];
extern const char *prjrev_errmsg[];

#ifndef PI
#define PI      3.141592653589793238462643
#endif
#define D2R	PI/180.0
#define R2D	180.0/PI
#define SQRT2	1.4142135623730950488
#define SQRT2INV 1.0/SQRT2

#define PRJSET 137

#ifdef __cplusplus
};
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern int npcode;
extern char pcodes[25][4];

struct celprm {
   int flag;
   double ref[4];
   double euler[5];

#if __STDC__  || defined(__cplusplus)
   int (*prjfwd)(const double, const double,
                 struct prjprm *,
                 double *, double *);
   int (*prjrev)(const double, const double,
                 struct prjprm *,
                 double *, double *);
#else
   int (*prjfwd)();
   int (*prjrev)();
#endif
};

#if __STDC__  || defined(__cplusplus)
   int celset(const char *, struct celprm *, struct prjprm *);
   int celfwd(const char *,
              const double, const double,
              struct celprm *,
              double *, double *,
              struct prjprm *,
              double *, double *);
   int celrev(const char *,
              const double, const double,
              struct prjprm *,
              double *, double *,
              struct celprm *,
              double *, double *);
#else
   int celset(), celfwd(), celrev();
#endif

extern const char *celset_errmsg[];
extern const char *celfwd_errmsg[];
extern const char *celrev_errmsg[];

#define CELSET 137

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(__STDC__) && !defined(__cplusplus)
#ifndef const
#define const
#endif
#endif

struct linprm {
   int flag;
   int naxis;
   double *crpix;
   double *pc;
   double *cdelt;

   /* Intermediates. */
   double *piximg;
   double *imgpix;
};

#if __STDC__  || defined(__cplusplus)
   int linset(struct linprm *);
   int linfwd(const double[], struct linprm *, double[]);
   int linrev(const double[], struct linprm *, double[]);
   int matinv(const int, const double [], double []);
#else
   int linset(), linfwd(), linrev(), matinv();
#endif

extern const char *linset_errmsg[];
extern const char *linfwd_errmsg[];
extern const char *linrev_errmsg[];

#define LINSET 137

#ifdef __cplusplus
};
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct wcsprm {
   int flag;
   char pcode[4];
   char lngtyp[5], lattyp[5];
   int lng, lat;
   int cubeface;
};

#if __STDC__ || defined(__cplusplus)
   int wcsset(const int,
              const char[][9],
              struct wcsprm *);

   int wcsfwd(const char[][9],
              struct wcsprm *,
              const double[],
              const double[],
              struct celprm *, 
              double *,
              double *, 
              struct prjprm *, 
              double[], 
              struct linprm *,
              double[]);

   int wcsrev(const char[][9],
              struct wcsprm *,
              const double[], 
              struct linprm *,
              double[], 
              struct prjprm *, 
              double *,
              double *, 
              const double[], 
              struct celprm *, 
              double[]);

   int wcsmix(const char[][9],
              struct wcsprm *,
              const int,
              const int,
              const double[],
              const double,
              int,
              double[],
              const double[],
              struct celprm *,
              double *,
              double *,
              struct prjprm *,
              double[], 
              struct linprm *,
              double[]);

#else
   int wcsset(), wcsfwd(), wcsrev(), wcsmix();
#endif

extern const char *wcsset_errmsg[];
extern const char *wcsfwd_errmsg[];
extern const char *wcsrev_errmsg[];
extern const char *wcsmix_errmsg[];

#define WCSSET 137

#ifdef __cplusplus
};
#endif

#endif /* wcslib_h_ */

/* Apr 15 1998	"deg" added to function names by Doug Mink, SAO
 * May 27 1998	ifndef WCSLIB changed to ifndef wcslib_h_
 * May 27 1998	ifndef WCSTRIG changed to ifndef wcstrig_h_
 * May 27 1998	ifndef CEL changed to ifndef cel_h_
 * May 27 1998	ifndef LIN changed to ifndef lin_h_
 */