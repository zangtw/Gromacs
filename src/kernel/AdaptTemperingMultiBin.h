/* multiple-bin estimators of thermodynamical quantities */
#ifndef _ADAPTIVE_TEMPERING_H__
#define _ADAPTIVE_TEMPERING_H__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "typedefs.h"
#include "physics.h"

#define ZCOM_PICK
#define ZCOM_CFG
#include "zcom1.h"

#ifndef HAVE_BOOL
#define bool int
#endif

/* raw data */
typedef struct {
  double    s;
  double    se;
  double    se2;
  double    se3;
} sm_t;

/* multiple-bin estimator parameters */
typedef struct {
  double    bmin;     /* minimal beta (highest temperature) */
  double    bmax;     /* maximal beta (lowest temperature) */
  double    bdel;     /* bin size of beta */
  int       n;        /* number of temperature bins */
  double    *barr;    /* temperature array */
  double    beta;     /* current value of beta */
  int       m;        /* maximal number of bins in a window */
  int       order;    /* order, should be 1 */
  unsigned  flags;    /* combination of flags */
  int       bwmod;    /* 0: d(beta) 1: dT/T  2: d(kT) */
  double    bwdel;    /* delta lnT */
  int       *js;          /* lower  boundary of asym. beta windows for ehat (asym.) */
  int       *jt;          /* higher boundary of asym. beta windows for ehat (asym.) */
  int       *jset;        /* lower  boundary of beta windows for et (usu. sym.), i - jset[i] == jtet[i] - (i+1) */
  int       *jtet;        /* higher boundary of beta windows for et (usu. sym.), i - jset[i] == jtet[i] - (i+1) */
  int       nstrefresh;   /* interval of recalculating et for all temperature */
  int       av_nstsave;   /* interval of writing mbav and ze files */
  int       av_binary;    /* use binary format in mbav file */
  char      *av_file;     /* name of mbav file */
  char      *ze_file;     /* name of ze file */
  int       wze_reps;     /* number of iterations before writing ze file */
  double    *vis;         /* number of visits */
  double    totvis;       /* total number of visits, number of tempering */
  double    *winstot;     /* total of sum.s over a multiple-bin temperature window */
  /* langevin equation */
  double    lgv_dt;           /* time step for the temperature Langevin eq */
  double    lgv_dTmax;        /* maximal amount of temperature change in a step */
  double    lgv_rej;          /* number of attempts of langevin equation trying to change beta too drastically */
  double    lgv_tot;          /* total number of attempts of using langevin equation */
  int       regl;             /* average within a bin first */
  double    fracmin;          /* minimal allowable coefficient during left/right combination */
  double    cvshiftmax;       /* maximal fraction for shift energy fluct. if cv is monotonic, it should be 0.0, for ising model, it can restrain the magnitude */
  /* parameter for temperature distribution perference */
  double    ens_exp;  /* ensemble exponent of beta */
	int       mode;         /* mode 0: disable 1: Gaussian 2:exponential */
	double    beta0;        /* beta0 in Gaussian distribution */
	double    invsigma2;    /* 1/sigma^2 in Gaussian distribution */
  double    c;            /* c in exponential distribution */
  double    *ens_w;   /* array of ensemble weights at bin boundaries */
  /* shrink parameters */
  double    shk_base;     /* current generic shrink amplitude */
  int       shk_winadj;   /* adjust shrink according to temperature window width */
  double    shk_max;      /* initial and maximal shrink (adjusted) */
  double    *shk_gauge;   /* array used of modulation shrinking factors */
  int       shk_mode;     /* 0: const, 1: amp/t, 2: amp/t^exp */
  double    shk_min;      /* minimal value for enforcing acc. sampling */
  int       shk_stop;     /* stop shrinking after this number of steps */
  double    shk_amp;      /* amp t^(-exp) */
  double    shk_exp;      /* amp t^(-exp) */
  /* reconstructed averages */
  double    *lnz;     /* logarithm of the partition function */
  double    *ehat;    /* internal energy */
  double    *cvhat;   /* heat capacity */
  double    *et;      /* bin-averaged internal energy */
  double    *imbal;   /* |a+ - a-| / (a+ + a-) for left-right combination */
  int       *haset;   /* current et[i] is reasonably good */
  unsigned  *qua;     /* bits represent whether estimated values are unbiased */
  double    *ampf;    /* currently amplification factor for adaptive averaging */
  sm_t  *sums;    /* normal data */
  int   has_xsums;
  int   cnt_int;      /* number of additional integer variables to be written to binary file */
  int   cnt_dbl;      /* number of additional double variables to be written to binary file */
  int   *idxcnt;      /* index count */
  int   *midx;        /* index look-up table */
  sm_t  *xsums;       /* multiple-bin damping data */
  /* energy histogram stuff, a new fold */
  int   eh_mode;  /* 0: disable; 1: simple histogram */
  int       eh_skip;      /* interval of reconstructing energy histograms */
  int       eh_bwmod;     /* 0: d(beta) 1: dT/T  2: d(kT) */
  double    eh_bwdel;     /* delta lnT */
  double    eh_min;       /* minimal energy */
  double    eh_max;       /* maximal energy */
  double    eh_del;       /* energy bin size */
  int       eh_cnt;       /* number of energy bins */
  int       eh_binary;    /* binary format for ehist file */
  int       eh_nstsave;   /* interval of writing histogrm files */
  char      *eh_file;     /* name of ehist file */
  char      *eh_rfile;    /* name of reconstructed energy histogram */
  double    *eh_his;      /* energy histogram data */
  double    *eh_recon;    /* temporary space for reconstructing histogram */
  int       *eh_is;       /* indices for temperature windows (lower) */
  int       *eh_it;       /* indices for temperature windows (higher) */
} mb_t;

#define   MB_DAMP       0x00000001    /* use adaptive averaging */
#define   MB_CV         0x00000002    /* compute heat capacity */
#define   MB_SYMWIN     0x00000004    /* use symmetrical window */
#define   MB_ONEBIN     0x00000020    /* use single bin estimator */
#define   MB_VERBOSE    0x00001000    /* being verbose */
#define   MB_SBCORR     0x00002000    /* include energy fluctuation correction due to a small bin width for internal energy, etc */
#define   MB_EH_ADDAHALF    0x00010000    /* add a half energy bin width in output */
#define   MB_EH_KEEPEDGE    0x00020000    /* keep zero edge at sides */
#define   MB_EH_NOZEROES    0x00040000    /* do not output zeroes */

#define MBQ_ET    0x00000001  /*  et quality bit   */
#define MBQ_EHAT  0x00000002  /*  ehat quality bit */
#define MBQ_CV    0x00000004  /*  cv quality bit   */
#define MBQ_LNZ   0x00000008  /*  lnz quality bit  */
  
#define MB_LOOSE  0x00000010  /*  temporarily allow empty temperature windows */

#define mb_close(mb) { mb_close_low(mb); free(mb); mb = NULL; }

/* estimate `ergt' at the current temperature `bet'
 * calculate the new bet after integrating Langevin equation,
 * dt is the time step for kT or 1/beta; ndlnwfdbeta = -d(lnwf)/d(beta) */
double mb_move(mb_t *mb, double erg, double bet, int ib,
		      double ndlnwfdbeta, double noise, double *ergt);

/* write various averages to ze_file */
int mb_wze(mb_t *mb, const char *fname);

/* use the integral identity to reconstruct an unbiased histogram (robust method) */
int mb_eh_recon(mb_t *mb, const char *fname);

/* return the reciprocal weight of multicanonical ensemble */
double mb_ensinvwf(mb_t *mb, double beta, double *f, double *ndfdbeta);

/* add energy and bet */
void mb_add(mb_t *mb, double e, const double v[], double bet,
	    int *pib, double *pinvwf, double *ndlnwfdbeta);

/* compute the average energy Et at current bin ib */
double mb_calc_et(mb_t *mb, int ib, int flags);

/* recompute all average energy */
void mb_refresh_et(mb_t *mb, int reps);

/* prepare and write mb data  */
int mb_write(mb_t *mb);

/* mb_cfgopen_low: initialize members of mb_t from configuration
 * file `cfg', or if unavailable, from default values */
int mb_cfgopen_low(mb_t *mb, cfgdata_t *cfg, double bmin, double bmax, char suffix);

/* mb_cfgopen: return a pointer of an initialized mb_t
 * if possible, initial values are taken from configuration
 * file `cfg', otherwise default values are assumed */
mb_t *mb_cfgopen(cfgdata_t *cfg, double bmin, double bmax, char suffix);

/* mb_close_low: close a pointer to mb_t */
void mb_close_low(mb_t *mb);

/* mb_clear: clear mb_t data */
void mb_clear(mb_t *mb);

/* mb_readbin_low: read mb_t data as binary */
int mb_readbin_low(mb_t *mb, FILE *fp, int ver, int endn);

/* mb_readbin: read mb_t data as binary */
int mb_readbin(mb_t *mb, const char *fname, int *pver);

/* mb_writebin_low: write mb_t data as binary */
int mb_writebin_low(mb_t *mb, FILE *fp, int ver);

/* mb_writebin: write mb_t data as binary */
int mb_writebin(mb_t *mb, const char *fname, int ver);

/* mb_eh_clear: clear mb_t data/eh */
void mb_eh_clear(mb_t *mb);

/* mb_eh_readbin_low: read mb_t/eh data as binary */
int mb_eh_readbin_low(mb_t *mb, FILE *fp, int ver, int endn);

/* mb_eh_readbin: read mb_t/eh data as binary */
int mb_eh_readbin(mb_t *mb, const char *fname, int *pver);

/* mb_eh_writebin_low: write mb_t/eh data as binary */
int mb_eh_writebin_low(mb_t *mb, FILE *fp, int ver);

/* mb_eh_writebin: write mb_t/eh data as binary */
int mb_eh_writebin(mb_t *mb, const char *fname, int ver);

/* mb_manifest: clear mb_t data */
void mb_manifest(mb_t *mb, FILE *fp, int arrmax);

#undef ZCOM_PICK
#undef ZCOM_CFG
#endif /* FILE */
