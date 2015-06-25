#include "AdaptTemperingMultiBin.h"

#define ZCOM_PICK
#define ZCOM_ENDN
#include "zcom1.h"

/*################################*/
/* raw data operations -- all static */

/* multiply everything in sm by a factor `f' */
static void sm_mul(sm_t *sm, double f)
{
  sm->s *= f;
  sm->se *= f;
  sm->se2 *= f;
  sm->se3 *= f;
}

/* add a new energy to sm_t with a weight `w' */
static void sm_addE(sm_t *sm, double w, double e, bool bCv)
{
  double s, sp, se, sep, see, seep, de, dep, de2;

  sm->s = sp = (s = sm->s) + w;
  sm->se = sep = (se = sm->se) + e*w;
  if (s <= 0.0) return;
  de = e - se/s;     /* old sum */
  dep = e - sep/sp;  /* new sum */
  sm->se2 = seep = (see = sm->se2) + (de2 = de * dep) * w; 
  if (bCv) sm->se3 += ((de2 -seep/s) - 2.0*see/s)* dep *w;
}

static int sm_cfgopen_low(sm_t *sm, cfgdata_t *cfg)
{
  if (sm == NULL) {
    fprintf(stderr, "null pointer to sm_t\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  sm->s = 0.0;
  sm->se = 0.0;
  sm->se2 = 0.0;
  sm->se3 = 0.0;
  (void) cfg;
  return 0;
ERR:
  return -1;
}

static void sm_close_low(sm_t *sm)
{
  memset(sm, 0, sizeof(*sm));
}

static void sm_clear(sm_t *sm)
{
  sm->s = 0.0;
  sm->se = 0.0;
  sm->se2 = 0.0;
  sm->se3 = 0.0;
}

static int sm_readbin_low(sm_t *sm, FILE *fp, int ver, int endn)
{
  if (sm == NULL) {
    fprintf(stderr, "passing null pointer to sm_readbin_low\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    return -1;
  }
  /* clear data before reading */
  sm_clear(sm);

  /* s */
  if (endn_fread(&sm->s, sizeof(sm->s), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading sm->s\n");
    goto ERR;
  }
  /* se */
  if (endn_fread(&sm->se, sizeof(sm->se), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading sm->se\n");
    goto ERR;
  }
  /* se2 */
  if (endn_fread(&sm->se2, sizeof(sm->se2), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading sm->se2\n");
    goto ERR;
  }
  /* se3 */
  if (endn_fread(&sm->se3, sizeof(sm->se3), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading sm->se3\n");
    goto ERR;
  }
  (void) ver;
  return 0;
ERR:
  sm_clear(sm);
  return -1;
}

static int sm_writebin_low(sm_t *sm, FILE *fp, int ver)
{
  if (sm == NULL) {
    fprintf(stderr, "passing null pointer to sm_writebin_low\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    return -1;
  }
  /* s */
  if (endn_fwrite(&sm->s, sizeof(sm->s), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing sm->s\n");
    goto ERR;
  }
  /* se */
  if (endn_fwrite(&sm->se, sizeof(sm->se), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing sm->se\n");
    goto ERR;
  }
  /* se2 */
  if (endn_fwrite(&sm->se2, sizeof(sm->se2), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing sm->se2\n");
    goto ERR;
  }
  /* se3 */
  if (endn_fwrite(&sm->se3, sizeof(sm->se3), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing sm->se3\n");
    goto ERR;
  }
  (void) ver;
  return 0;
ERR:
  return -1;
}

static void sm_manifest(sm_t *sm, FILE *fp, int arrmax)
{
  /* s */
  fprintf(fp, "sm->s: double, %g\n", sm->s);

  /* se */
  fprintf(fp, "sm->se: double, %g\n", sm->se);

  /* se2 */
  fprintf(fp, "sm->se2: double, %g\n", sm->se2);

  /* se3 */
  fprintf(fp, "sm->se3: double, %g\n", sm->se3);

  (void) arrmax;
}

/*################################*/
/* mb operations, static part*/

/* check if mb->barr is arranged in an ascending order */
static int mb_check_barr(mb_t *mb)
{
  int i;

  for (i = 0; i <= mb->n; i++)
    if (i > 0 && mb->barr[i] <= mb->barr[i-1]) {
      fprintf(stderr, "barr should ascend: barr[%d] = %g, barr[%d] = %g\n",
          i, mb->barr[i], i-1, mb->barr[i-1]);
      return 1;
    }
  return 0;
}

/* setup temperature windows, [ajs[i], ajt[i]) for each i in [0..n],
 * and symmetrical ones [ajset[i], ajtet[i]) for each i in [0..n) */
static int mb_mkwin(mb_t *mb, int bwmod, double bwdel,
    int ajs[], int ajt[], int ajset[], int ajtet[])
{
  int i, n, js, jt, idel, di1, di2, di;
  double bet, dbet = 0.0;

  die_if (mb == NULL || ajs == NULL || ajt == NULL,
      "null pointer mb: %p, ajs: %p, ajt: %p", mb, ajs, ajt);
  n = mb->n;
  for (i = 0; i <= n; i++) {
    bet = mb->barr[i];
    switch (bwmod) {
      case 0: dbet = bwdel; break;
      case 1: dbet = bwdel * bet; break;
      case 2: dbet = bwdel * (bet * bet); break;
      default: die_if(1, "bad bwmod=%d\n", bwmod); break;
    }
    idel = (int)(dbet/mb->bdel + 0.50000001);
    if (idel < 1) idel = 1; /* although idel = 0 should be fine */
    if ((js = i - idel) < 0) js = 0;
    if ((jt = i + idel) > n) jt = n;
    die_if (i < js || i > jt, "bad window %d (%d,%d)\n", i, js, jt);
    ajs[i] = js;
    ajt[i] = jt;
  }

  if (ajset == NULL || ajtet == NULL) return 0;
  /* calculate jset and jtet */
  for (i = 0; i < n; i++) {
    if (mb->flags & MB_SYMWIN) {
      di1 = i - ajs[i];
      di2 = ajt[i+1] - (i + 1);
      /* choose the smaller one from di1 and di2 */
      di  = (di1 > di2) ? di2 : di1;
      ajset[i] = i - di;
      ajtet[i] = i + 1 + di;
    } else {
      ajset[i] = ajs[ i ];
      ajtet[i] = ajt[ i + 1 ];
    }
  }
  return 0;
}

/* construct a more symmetrical bin for temperature mb->barr[i] than
 * (mb->js[i], mb->jt[i]) for quantities at a bin boundary, e.g. ehat */
static void mb_mksymwin(mb_t *mb, int i, int *js, int *jt)
{
  int j;

  j = (i < mb->n) ? i : (i - 1);
  (*js) = mb->jset[j];
  (*jt) = mb->jtet[j];
  if (i > 0 && i < mb->n) {
    if (i < mb->n/2) (*jt)--;
    else (*js)--;
  }
  die_if (*jt - *js <= 0, "empty window (%d,%d) for %d\n", *js, *jt, i);
}

/* compute the largest span of temperature windows in number of bins */
static int mb_maxwinspan(mb_t *mb)
{
  int i, cnt, max;

  /* first calculate the maximal # of bins, (attempted value) */
  for (max = 0, i = 0; i < mb->n; i++) {
    cnt = mb->jt[i + 1] - mb->js[i];
    if (cnt > max) max = cnt;
  }
  return mb->m = max;
}

/* adjust mb->m as the max. number of estimators affected by a single bin,
 * idxcnt needs to be allocated before calling */
static int mb_maxidxcnt(mb_t *mb)
{
  int i, j, js, jt, max;

  if (!(mb->flags & MB_DAMP)) return mb->m;
  /* mb->idxcnt[j]: number of estimators affected by stat. in bin j */
  for (i = 0; i < mb->n; i++) {
    js = mb->js[i];
    jt = mb->jt[i + 1];
    for (j = js; j < jt; j++)
      mb->idxcnt[j]++;
  }
  for (max = 0, j = mb->n; j; j--) /* compute max { mb->idxcnt[j] } */
    if (mb->idxcnt[j] > max) max = mb->idxcnt[j];
  if (max > mb->m) {
    if (mb->flags & MB_VERBOSE) printf("mb->m: %d => %d\n", mb->m, max);
    mb->m = max;
  }
  return mb->m;
}

/* set up matrix indices mb->midx, need mb->js, mb->jt */
static int mb_mkidx(mb_t *mb)
{
  int i, j, cnt, js, jt;

  if (!(mb->flags & MB_DAMP)) return 0;
  /* clear mb->idxcnt, whose space is to be used
   * to dynamically construct the matrix indices
   * mb->idxcnt will be restored at the end */
  for (j = 0; j < mb->n; j++) mb->idxcnt[j] = 0;
  for (i = 0; i < mb->n; i++) {
    js = mb->js[ i ];
    jt = mb->jt[ i + 1 ];
    for (j = js; j < jt; j++) {
      cnt = mb->idxcnt[j];
      die_if (cnt >= mb->m, 
        "cnt: %d, m: %d, i: %d, j: %d, (js, jt) = (%d, %d).\n",
          cnt, mb->m, i, j, js, jt);
      die_if (j - js >= mb->m,
        "j: %d, js: %d, m: %d\n", j, js, mb->m);
      mb->midx[ j * mb->m + cnt ] = i * mb->m + j - js;
      mb->idxcnt[ j ] = cnt + 1;
    }
  }
  return 0;
}

/* normalize damping weight back to 1.0, averages are not affected 
 * if i < 0, do for all estimators */
static void mb_normalize(mb_t *mb, int i)
{
  int i0, i1, j, js, jt;
  double fac;

  if (!(mb->flags & MB_DAMP)) return;
  if (i >= 0) { i0 = i1 = i; }
  else { i0 = 0; i1 = mb->n; } /* all estimators, if i < 0 */
  for (i = i0; i <= i1; i++) { /* loop over estimators */
    if (fabs(mb->ampf[i] - 1.0) < 1e-12) continue;
    fac = 1.0 / mb->ampf[i];
    mb->ampf[i] = 1.0;
    js = mb->js[i];
    jt = mb->jt[i+1];
    for (j = js; j < jt; j++)
      sm_mul(mb->xsums + i * mb->m + j - js, fac);
  }
}

/* compute the temperature-independent shrinking factor */
static double mb_shkbase(mb_t *mb)
{
  double x, shk;

  if (mb->shk_stop >= 0 && mb->totvis > mb->shk_stop)
    return 0.0;
  else if (mb->shk_mode == 0)
    return mb->shk_min;
  shk = mb->shk_max;
  if (mb->totvis > 10 * mb->n) {
    x = mb->totvis / mb->n;
    if (mb->shk_mode == 1) {
      x = mb->shk_amp / x;
    } else if (mb->shk_mode == 2) {
      x = mb->shk_amp / pow(x, mb->shk_exp);
    } else die_if(1, "invalid shk_mode: %d\n", mb->shk_mode);
    if (x < shk) shk = x;
    if (shk < mb->shk_min) shk = mb->shk_min;
  }
  return mb->shk_base = shk;
}

/* compute the adjusted shk and invgam */
static double mb_invgam(mb_t *mb, int ib)
{
  double shk;

  shk = mb_shkbase(mb); /* compute the unadjusted */
  if (mb->shk_winadj) { /* multiply the gauge */
    die_if (mb->shk_gauge == NULL, "gauge is null\n");
    die_if (ib < 0 || ib >= mb->n, "index %d out of range\n", ib);
    shk *= mb->shk_gauge[ib];
    if (shk > mb->shk_max) shk = mb->shk_max;
  }
  return 1.0 / (1.0 - shk);
}

/* compute total weighted number of visits to T. window of bin i */
static void mb_winstot(mb_t *mb)
{
  int i, j, js, jt, damp;
  double tot;
  sm_t *sm0;

  damp = (mb->flags & MB_DAMP);
  for (i = 0; i < mb->n; i++) {
    js = mb->jset[i];
    jt = mb->jtet[i];
    die_if (js < mb->js[i] || i < js || i >= jt,
      "bad window (%d, %d), js=%d, i=%d, n=%d\n", js, jt, mb->js[i], i, mb->n);
    sm0 = damp ? (mb->xsums + i*mb->m - js) : (mb->sums);
    for (tot = 0.0, j = js; j < jt; j++) tot += sm0[j].s;
    mb->winstot[i] = damp ? (tot/mb->ampf[i]) : tot;
  }
}

/* set quality bit */
static void mb_setqbit(unsigned *ptr, unsigned mask, int on)
{
  if (mask) { if (on) *ptr |= mask; else *ptr &= ~mask; }
}

/* translate quality bits into a 0-1 string */
static char *mb_qua2s(unsigned i)
{
  static char buf[64]; /* has to be static to be the return value */
  int  cnt = 0;

  buf[cnt++] = (char)((i & MBQ_ET  ) ? '1' : '0');
  buf[cnt++] = (char)((i & MBQ_EHAT) ? '1' : '0');
  buf[cnt++] = (char)((i & MBQ_CV  ) ? '1' : '0');
  buf[cnt++] = (char)((i & MBQ_LNZ ) ? '1' : '0');
  buf[cnt] = '\0';
  return buf;
}

/* collect moments from i's left and i's right */
static void mb_lrcol(mb_t *mb, int i, int js, int jt,
    double t1[], double *tb, double s0[], double s1[],
    sm_t *sm0, int bcorr)
{
  double del, et0, el, var, x;
  int j, jx, lr, regl = mb->regl;

  die_if (0 > js || js >= jt || jt > mb->n, "bad window (%d, %d)", js, jt);
  s0[0] = s0[1] = s1[0] = s1[1] = t1[0] = t1[1] = 0.0;
  if (tb != NULL) *tb = 0.0;
  for (j = js; j < jt; j++) { /* loop over bins */
    sm_t *sm = sm0 + j;
    if (fabs(sm->s) < 1e-6) continue; /* empty bin */
    if (j < i) { jx = js; lr = 0; } 
    else       { jx = jt; lr = 1; }
    /* correct energy fluctuation for small bin size: 
     *   < (E - el)^2 > + (el - et0)^2,
     * el is the average energy of the bin */
    el = sm->se/sm->s;
    et0 = mb->haset[j] ? mb->et[j] : el;
    del = bcorr ? (el - et0) : 0.0;
    var = sm->se2/sm->s + del * del;
    if (regl) {
      s0[lr] += 1;
      s1[lr] += sm->se / sm->s;
      x = var;
    } else {
      s0[lr] += sm->s;
      s1[lr] += sm->se;
      x = var * sm->s;
    }
    t1[lr] += x * (j - jx + 0.5);
    if (j == i && tb != NULL) *tb = 0.5 * x;
  }
}

/* collect second order moments from i's left and i's right */
static void mb_lrcol2(mb_t *mb, int i, int js, int jt,
    double t2[], double s0[], double s2[], sm_t *sm0)
{
  double x;
  int j, jx, lr, regl = mb->regl;

  die_if (0 > js || js >= jt || jt > mb->n, "bad window (%d, %d)", js, jt);
  s0[0] = s0[1] = s2[0] = s2[1] = t2[0] = t2[1] = 0.0;
  for (j = js; j < jt; j++) { /* loop over bins */
    sm_t *sm = sm0+j;
    if (fabs(sm->s) < 1e-6) continue;
    if (j < i) { jx = js; lr = 0; }
    else       { jx = jt; lr = 1; }
    if (regl) {
      s0[lr] += 1;
      s2[lr] += sm->se2/sm->s;
      x = sm->se3/sm->s;
    } else {
      s0[lr] += sm->s;
      s2[lr] += sm->se2;
      x = sm->se3;
    }
    t2[lr] += x * (j - jx + 0.5);
  }
  t2[0] *= -mb->bdel;
  t2[1] *= -mb->bdel;
}

/* return an estimate from combining data from left and right
 * *imbal returns the difference between two coefficients
 * *good returns if an estimate is successful */
static double mb_lrbal(mb_t *mb, int ib, int js, int jt, int loose,
    const double t1[], double tb, const double s0[], const double s1[],
    double *imbal, int *good, unsigned qbit)
{
  double a, b;  /* combination coefficients, a for left, b for right */
  double num, den, del, tmp;
  int ip, im;

  a = b = 1.0;
  mb_setqbit(&mb->qua[ib], qbit, 0);

  if (fabs(s0[0] + s0[1]) < 1e-3) { /* no data in the entire window */
    if (good != NULL) *good = 0;
    die_if (!loose, "an empty window ib = %d\n", ib);
    return 0.0;
  }
  ip = jt - ib;
  im = ib - js;
  num = t1[0] + t1[1] + tb * (jt - js);
  den = t1[0] * ip - t1[1] * im; /* t1[0] >= 0 and t1[1] <= 0, so den >= 0 */
  if (fabs(den) < 1e-8) goto FALLBACK; /* visited at most once */
  del = num / den;
  if ((a = 1.0 - del*ip) < 0.0) a = 0.0;
  if ((b = 1.0 + del*im) < 0.0) b = 0.0;
  if (fabs(a + b) < 1e-8) goto FALLBACK;
  tmp = 1.0/(a + b);
  a *= tmp;
  b *= tmp;
  if (a < mb->fracmin) {
    b = 1.0 - (a = mb->fracmin);
  } else if (b < mb->fracmin) {
    a = 1.0 - (b = mb->fracmin);
  } else {
    mb_setqbit(&mb->qua[ib], qbit, 1);
  }
FALLBACK:
  /* in case denominator is 0.0, e.g., a*s0[0] + b*s0[1] = 1*0 + 0*3  */
  if (fabs(a * s0[0] + b * s0[1]) < 1e-6) {
    a = b = 1.0;
    mb_setqbit(&mb->qua[ib], qbit, 0);
  }
  if (imbal != NULL) *imbal = (a - b) / (a + b);
  if (good != NULL) *good = 1.0;
  return (a * s1[0] + b * s1[1]) / (a * s0[0] + b * s0[1]);
}

/* update the estimated energy et of the given ib
 * using a combination of data from left and right
 * note `flags' is different from mb->flags
 * if MB_LOOSE is in `flags', failure due to empty data is allowed */
static double mb_etlr(mb_t *mb, int ib, int flags)
{
  sm_t *sm0;
  int js, jt;
  double s0[2], s1[2], t1[2], tb;  /* first order */
  static int once;

  if (once++ == 0) fprintf(stderr, "etlr: fracmin = %g\n", mb->fracmin);

  js = mb->jset[ib];
  jt = mb->jtet[ib];
  die_if (js < 0 || jt <= js || jt > mb->n, 
      "invalid indices %d, %d, ib = %d/%d", js, jt, ib, mb->n);

  /* choose xsums or sums */
  sm0 = (mb->flags & MB_DAMP) ? (mb->xsums + ib*mb->m-mb->js[ib]) : mb->sums;

  if ((mb->flags & MB_ONEBIN) || jt == js + 1) /* window reduced to a bin */
    return mb->et[ib] = (fabs(sm0[ib].s)>1e-6) ? (sm0[ib].se/sm0[ib].s) : 0.;

  /* using stat. from bins (js, jt) to form two estimates, left & right */
  mb_lrcol(mb, ib, js, jt, t1, &tb, s0, s1, sm0, 1);
  /* compute linear combination coefficients and et */
  return mb->et[ib] = mb_lrbal(mb, ib, js, jt, (flags & MB_LOOSE),
      t1, tb, s0, s1, &mb->imbal[ib], &mb->haset[ib], MBQ_ET);
}

static void mb_calc_ehat(mb_t *mb)
{
  int i, js, jt, needcv;
  sm_t *sm0;
  double bet, del, s0[2], s1[2], s2[2], t1[2], t2[2];
  static int once;

  if (once++ == 0) fprintf(stderr, "calc_ehat: fracmin = %g, cvshiftmax = %g\n", 
        mb->fracmin, mb->cvshiftmax);

  needcv = (mb->flags & MB_CV);
  for (i = 0; i <= mb->n; i++) {
    mb_mksymwin(mb, i, &js, &jt);
    bet = mb->barr[i];
    if (mb->flags & MB_DAMP) { /* try to use damp sum if possible */
      int ip = (i == mb->n) ? (i-1) : i;
      sm0 = mb->xsums + ip * mb->m - mb->js[ip];
    } else sm0 = mb->sums;
    mb_lrcol(mb, i, js, jt, t1, NULL, s0, s1, sm0, mb->flags & MB_SBCORR);
    if (fabs(s0[0] + s0[1]) < 1e-3) { /* no data in the entire window */
      mb_setqbit(&mb->qua[i], MBQ_EHAT, 0);
      mb_setqbit(&mb->qua[i], MBQ_CV, 0);
      continue;
    }
    mb->ehat[i] = mb_lrbal(mb, i, js, jt, 1,
        t1, 0.0, s0, s1, NULL, NULL, MBQ_EHAT);

    if (needcv) { /* calculate Cv */
      mb_lrcol2(mb, i, js, jt, t2, s0, s2, sm0);
      if (t2[0] * t2[1] > 0) { /* t2[0] and t2[1] share the same sign */
        del = (t2[0] + t2[1]) / (s2[0] + s2[1]);
        if (del < -1) del = -1;
        if (fabs(del) > mb->cvshiftmax)
          del = mb->cvshiftmax * ((del > 0) ? 1.0 : -1.0);
        mb_setqbit(&mb->qua[i], MBQ_CV, 1);
        mb->cvhat[i] = BOLTZ * (bet * bet) * 
                       (s2[0]+s2[1]) * (1+del) / (s0[0]+s0[1]);
        /* = (s2[0]+s2[1]+t2[0]+t2[1])/(s0[0]+s0[1])*bet*bet; */
      } else /* normal case */
        mb->cvhat[i] = BOLTZ * (bet * bet) * 
          mb_lrbal(mb, i, js, jt, 1, t2, 0.0, s0, s2, NULL, NULL, MBQ_CV);
    }
  }
  /* estimate the partition function */
  for (mb->lnz[0] = 0.0, i = 0; i < mb->n; i++) 
    mb->lnz[i+1] = mb->lnz[i] + mb->et[i] * (mb->barr[i] - mb->barr[i+1]);
}

/*################################*/
/* mb operations, non-static part*/

double mb_move(mb_t *mb, double erg, double bet, int ib,
          double ndlnwfdbeta, double noise, double *ergt)
{
  double dkt, dktmax, kt1, kt2, dt, bet2, rndmag;

  *ergt = mb_calc_et(mb, ib, 0);
  dt = mb->lgv_dt * BOLTZ;
  dktmax = mb->lgv_dTmax * BOLTZ;
  kt1 = 1.0/bet;
  rndmag = kt1*sqrt(2.0*dt);
  
  /* Langevin integration */
  dkt  = (erg - *ergt + ndlnwfdbeta)*dt + rndmag * noise;
  
  if (dkt > dktmax) {
    dkt = dktmax;
    mb->lgv_rej += 1.0;
  } else if (dkt < -dktmax) {
    dkt = -dktmax;
    mb->lgv_rej += 1.0;
  }
  mb->lgv_tot += 1.0;
  kt2 = kt1 + dkt;
  bet2 = 1.0 / kt2;
  mb->beta = (bet2 < mb->bmax && bet2 > mb->bmin) ? bet2 : bet;
  mb->beta = (mb->beta >= mb->bmax - 1e-5) ? (mb->bmax - 1e-5) : mb->beta;
  mb->beta = (mb->beta <= mb->bmin + 1e-5) ? (mb->bmin + 1e-5) : mb->beta;
  return mb->beta;
}

int mb_wze(mb_t *mb, const char *fname)
{
  int i, ip;
  FILE *fp;

  if (fname == NULL) fname = mb->ze_file;
  die_if (fname == NULL, "file name is NULL");

  mb_winstot(mb);
  mb_calc_ehat(mb);
  if ((fp = fopen(fname, "w")) == NULL) {
    fprintf(stderr, "cannot open %s.\n", fname);
    return 1;
  }
  for (i = 0; i <= mb->n; i++) {
    fprintf(fp, "%16.10f %20.10f %22.10f %22.10f ",
      mb->barr[i], mb->lnz[i], mb->ehat[i], mb->cvhat[i]);
    ip = (i < mb->n) ? i : (i-1); /* for quantities with no [mb->n] */
    fprintf(fp, " %22.10f %s %+10.6f %22.10e %22.10e %22.10e %22.10e",
      mb->et[ip], mb_qua2s(mb->qua[i]), mb->imbal[ip], mb->sums[ip].s, 
      mb->vis[ip], mb->shk_gauge[ip], mb->winstot[ip]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

int mb_eh_recon(mb_t *mb, const char *fname)
{
  FILE *fp;
  int ib, j, js, jt, ie, imin, imax, cols, full, keep0;
  double eav, db, x;
  double num, den;
  double del, base, inc;

  if (mb->eh_mode == 0) return 0;
  die_if (mb->eh_mode != 1, "invalid eh_mode %d\n", mb->eh_mode);
  if ((fp = fopen((fname != NULL) ? fname : mb->eh_rfile, "w")) == NULL) {
    fprintf(stderr, "cannot write reconstructed histogram [%s].\n",
        mb->eh_rfile);
    return 1;
  }
  full = mb->flags & MB_EH_KEEPEDGE;
  keep0 = !(mb->flags & MB_EH_NOZEROES);
  del = (mb->flags & MB_EH_ADDAHALF) ? 0.5 : 0; /* for continuous system */
  cols = mb->eh_cnt;
  base = mb->eh_min;
  inc = mb->eh_del;
  db = mb->bdel;
  for (mb->lnz[0] = 0.0, ib = 0; ib < mb->n; ib++) /* build lnZ */
    mb->lnz[ib+1] = mb->lnz[ib] + mb->et[ib]*(mb->barr[ib] - mb->barr[ib+1]);
  
  /* loop over temperatures, and skip a few intermediate temperatures */
  for (ib = 0; ib <= mb->n; ib += mb->eh_skip) {
    /* reconstruct energy histogram at beta = mb->barr[ib] */
    js = mb->eh_is[ib];
    jt = mb->eh_it[ib];
    die_if(js < 0 || jt > mb->n || js >= jt, "bad window (%d, %d)\n", js, jt);
    /* loop over energy bins */
    for (ie = 0; ie < cols; ie++) {
      eav = base + (ie+del)*inc;
      for (den = 0, j = js; j <= jt; j++) { /* denominator */
        x = mb->ens_w[j] * exp(-eav*db*(j - ib) - mb->lnz[j] + mb->lnz[ib]);
        if (j == js || j == jt) x *= 0.5;
        den += x;
      }
      for (num = 0.0, j = js; j < jt; j++) { /* numerator */
        x = mb->eh_his[j*cols + ie];
        if (fabs(x) < 0.5) continue;
        num += x;
      }
      mb->eh_recon[ie] = num/den;
    }
    /* determine the output range */
    if (full) {
      imin = 0;
      imax = cols;
    } else { /* only print the smallest non-zero region */
      for (ie = cols-1; ie >= 0; ie--) if (mb->eh_recon[ie] > 0) break;
      if ((imax=ie+1) == 0) continue;
      for (imin = 0; imin < imax; imin++) if (mb->eh_recon[imin] > 0) break;
    }
    /* normalize the histogram, and output */
    for (x = 0, ie = imin; ie < imax; ie++) x += mb->eh_recon[ie];
    x = (fabs(x) > 1e-6) ? 1.0/(x*inc) : 1.0;
    for (ie = imin; ie < imax; ie++)
      if (keep0 || mb->eh_recon[ie] > 1e-6)
        fprintf(fp, "%g %.14E %g\n", 
          base + (ie+del)*inc, mb->eh_recon[ie]*x, mb->barr[ib]);
    fprintf(fp, "\n");
  }
  fclose(fp);
  return 0;
}

double mb_ensinvwf(mb_t *mb, double beta, double *pf, double *pndfdbeta)
{
  const double eps = 1e-5;
  double dif, invwf, f, ndfdbeta;
  int ifac;

  /* invwf: 1/w(beta)/f(beta);
   * f: f(beta);
   * ndfdbeta: -df(beta)/dbeta;
   */

  if(mb->mode == 1)
  {
    dif    = beta - mb->beta0;
    f      = exp(-0.5 * (dif * dif) * mb->invsigma2);
    ndfdbeta = f * dif * mb->invsigma2;
  }
  else if(mb->mode ==2)
  {
    f      = exp(-beta * mb->c);
    ndfdbeta = f * mb->c;
  }
  else 
  {
    f = 1.0;
    ndfdbeta = 0.0;
  }
  if (pf   != NULL) *pf   = f;
  if (pndfdbeta != NULL) *pndfdbeta = ndfdbeta;

  beta /= mb->bmax; /* to relative beta */
  ifac = (int)(mb->ens_exp + 0.5); /* round to nearest int */
  if (fabs(mb->ens_exp - ifac) < 1e-5 && ifac >= 0) {
    for (invwf = 1.0/f; ifac > 0; ifac--)
      invwf *= beta;
  } else  /* unable to avoid exponential */
    invwf = exp(mb->ens_exp * log(beta)) / f;

  die_if (invwf > 1e6 || invwf < 1e-6, "bad invwf=%g, beta=%g\n", invwf, beta);
  
  return invwf;
}

void mb_add(mb_t *mb, double e, const double v[], double bet,
    int *pib, double *pinvwf, double *ndlnwfdbeta)
{
  double ginvwf, invwf, f = 1.0, ndfdbeta = 0.0;
  int j, bCv = mb->flags & MB_CV;

  *pib = j = (int)((bet - mb->bmin)/mb->bdel);
  die_if (j < 0 || j >= mb->n, "beta = %d, %g, range = (%g, %g, %g)\n",
      j, bet, mb->bmin, mb->bdel, mb->bmax);
  mb->vis[j] += 1.0;
  mb->totvis += 1.0;

  /* f: f(beta);
   * ndfdbeta: -df/dbeta;
   * invwf: 1/w(beta)/f(beta);
   * ndlnwfdbeta: -dln(w(beta)f(beta))/dbeta;
   * ginvwf: adaptive weight = ampf * invwf;
   */
  *pinvwf = invwf = mb_ensinvwf(mb, bet, &f, &ndfdbeta); /* get weight */
  *ndlnwfdbeta = mb->ens_exp/bet + ndfdbeta/f;
  sm_addE(mb->sums + j, invwf, e, bCv);

  die_if(v != NULL, "v[] must be null\n");

  if (mb->eh_mode > 0) { /* add to energy histogram */
    int ie = (int)((e - mb->eh_min)/mb->eh_del);
    if (ie >= 0 && ie < mb->eh_cnt) /* no invw for histogram */
      mb->eh_his[j*mb->eh_cnt+ie] += 1.0;
  }

  if (mb->flags & MB_DAMP) { /* add to damping data */
    int i, l, mid, cnt, upd;

    cnt = mb->idxcnt[j];
    for (l = 0; l < cnt; l++) { /* loop over affected estimators */
      mid = mb->midx[j * mb->m + l];
      i = mid / mb->m;
      die_if (i * mb->m + j - mb->js[i] != mid, /* check index */
        "index corruption, i=%d, m=%d, j=%d, js=%d, mid=%d(%d)\n",
            i, mb->m, j, mb->js[i], mid, i * mb->m + j - mb->js[i]);
      /* et is computed from (jset, jtet), which is a subset of (js, jt), 
       * avoid weight updating when j lies outside of the former */
      if (mb->flags & MB_ONEBIN)
        upd = (j == i);
      else if (mb->flags & MB_SYMWIN)
        upd = (j >= mb->jset[i] && j < mb->jtet[i]);
      else upd = 1;
      /* apply adaptive averaging */
      if (upd) mb->ampf[i] *= mb_invgam(mb, i);
      ginvwf = mb->ampf[i] * invwf; /* multiply accumulated 1/gamma */
      sm_addE(mb->xsums+mid, ginvwf, e, bCv);
      /* we call normalization when the weight starts to blow up */
      if (ginvwf > 2.0) mb_normalize(mb, i);
    }
  }
}

double mb_calc_et(mb_t *mb, int ib, int flags)
{
  die_if (ib < 0 || ib >= mb->n, "bad ib %d [0, %d).\n", ib, mb->n);
  return mb_etlr(mb, ib, flags);
}

void mb_refresh_et(mb_t *mb, int reps)
{
  int i, rep;

  for (rep = 0; rep < reps; rep++) {
    for (i = 0; i < mb->n; i++)
      mb_calc_et(mb, i, MB_LOOSE);
  }
}

int mb_write(mb_t *mb)
{
  return mb_writebin(mb, mb->av_file, 1);
}

int mb_cfgopen_low(mb_t *mb, cfgdata_t *cfg, double bmin, double bmax, char suffix)
{
  int i;
  char buf[100];
  char suf[2];
  suf[0] = suffix;
  suf[1] = '\0';

  if (mb == NULL) {
    fprintf(stderr, "null pointer to mb_t\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  /* bmin: minimal beta (highest temperature) */
  mb->bmin = bmin;
  /* bmax: maximal beta (lowest temperature) */
  mb->bmax = bmax;
  /* bdel: bin size of beta */
  mb->bdel = 0.0001;
  if (cfg != NULL && 0 != cfgget(cfg, &mb->bdel, "beta_del", "%lf")) {
    fprintf(stderr, "missing var: mb->bdel, key: beta_del, fmt: %%lf\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  if ( !(mb->bdel > 1e-6) ) {
    fprintf(stderr, "mb->bdel: failed validation: mb->bdel > 1e-6\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }

  /* n: number of temperature bins */
  mb->n = (int)((mb->bmax - mb->bmin)/mb->bdel - 1e-5) + 1;
  /* barr: temperature array */
  if ((mb->barr = calloc((mb->n + 2), sizeof(double))) == NULL) {
    fprintf(stderr, "no memory! var: mb->barr, type: double\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n+1; i++)
    mb->barr[i] = mb->bmin + i * mb->bdel;
  /* check beta array */
  if ( !(mb_check_barr(mb) == 0) ) {
    fprintf(stderr, "check beta array\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  /* fix bmax to a bin boundary */
  mb->bmax = mb->bmin + mb->bdel * mb->n;
  /* beta: current value of beta */

  //just for debug use
  mb->beta = mb->bmax;
  /*mb->beta = mb->bmin;*/
  /*mb->beta = 0.5 * mb->bmin + 0.5 * mb->bmax;*/
  
  mb->beta = (mb->beta >= mb->bmax - 1e-5) ? (mb->bmax - 1e-5) : mb->beta;
  mb->beta = (mb->beta <= mb->bmin + 1e-5) ? (mb->bmin + 1e-5) : mb->beta;
  
  /* m: maximal number of bins in a window */
  mb->m = 0;
  /* order: order, should be 1 */
  mb->order = 1;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->order, "mbest_order", "%d")) {
    fprintf(stderr, "assuming default: mb->order = 1, key: mbest_order\n");
  }
  if ( !(mb->order == 1) ) {
    fprintf(stderr, "mb->order: failed validation: mb->order == 1\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }

  /* flags: combination of flags */
  mb->flags = 0;
  /* bwmod: 0: d(beta) 1: dT/T  2: d(kT) */
  mb->bwmod = 1;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->bwmod, "mbest_mbin_mode", "%d")) {
    fprintf(stderr, "assuming default: mb->bwmod = 1, key: mbest_mbin_mode\n");
  }
  if ( !(mb->bwmod >= 0 && mb->bwmod <= 2) ) {
    fprintf(stderr, "mb->bwmod: failed validation: mb->bwmod >= 0 && mb->bwmod <= 2\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }

  /* bwdel: delta lnT */
  mb->bwdel = 0.05;
  if (mb->bwmod == 1) {
    if (cfg != NULL && 0 != cfgget(cfg, &mb->bwdel, "mbest_delta_lnT", "%lf")) {
      fprintf(stderr, "missing var: mb->bwdel, key: mbest_delta_lnT, fmt: %%lf\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    if ( !(mb->bwdel > mb->bdel/pow(mb->bmin, 1.0)) ) {
      fprintf(stderr, "mb->bwdel: failed validation: mb->bwdel > mb->bdel/pow(mb->bmin, 1.0)\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }

  /* bwdel: delta beta */
  if (mb->bwmod == 0) {
    mb->bwdel = 0.02;
    if (cfg != NULL && 0 != cfgget(cfg, &mb->bwdel, "mbest_delta_beta", "%lf")) {
      fprintf(stderr, "missing var: mb->bwdel, key: mbest_delta_beta, fmt: %%lf\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    if ( !(mb->bwdel > mb->bdel/pow(mb->bmin, 0.0)) ) {
      fprintf(stderr, "mb->bwdel: failed validation: mb->bwdel > mb->bdel/pow(mb->bmin, 0.0)\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }

  /* bwdel: delta kT */
  if (mb->bwmod == 2) {
    mb->bwdel = 0.1;
    if (cfg != NULL && 0 != cfgget(cfg, &mb->bwdel, "mbest_delta_kT", "%lf")) {
      fprintf(stderr, "missing var: mb->bwdel, key: mbest_delta_kT, fmt: %%lf\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    if ( !(mb->bwdel > mb->bdel/pow(mb->bmin, 2.0)) ) {
      fprintf(stderr, "mb->bwdel: failed validation: mb->bwdel > mb->bdel/pow(mb->bmin, 2.0)\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }

  /* MB_DAMP: use adaptive averaging */
  i = 1;
  if (cfg == NULL || 0 != cfgget(cfg, &i, "mbest_damp", "%u")) {
    fprintf(stderr, "assuming default: MB_DAMP = 1, key: mbest_damp\n");
  }
  if ( !(i == 0 || i == 1) ) {
    fprintf(stderr, "MB_DAMP: failed validation: i == 0 || i == 1\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  if (i) {
    mb->flags |= MB_DAMP;
  } else {
    mb->flags &= ~MB_DAMP;
  }

  /* MB_CV: compute heat capacity */
  i = 1;
  if (cfg == NULL || 0 != cfgget(cfg, &i, "mbest_needcv", "%u")) {
    fprintf(stderr, "assuming default: MB_CV = 1, key: mbest_needcv\n");
  }
  if ( !(i == 0 || i == 1) ) {
    fprintf(stderr, "MB_CV: failed validation: i == 0 || i == 1\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  if (i) {
    mb->flags |= MB_CV;
  } else {
    mb->flags &= ~MB_CV;
  }

  /* MB_SYMWIN: use symmetrical window */
  i = 1;
  if (cfg == NULL || 0 != cfgget(cfg, &i, "mbest_sym_mbin", "%u")) {
    fprintf(stderr, "assuming default: MB_SYMWIN = 1, key: mbest_sym_mbin\n");
  }
  if ( !(i == 0 || i == 1) ) {
    fprintf(stderr, "MB_SYMWIN: failed validation: i == 0 || i == 1\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  if (i) {
    mb->flags |= MB_SYMWIN;
  } else {
    mb->flags &= ~MB_SYMWIN;
  }

  /* MB_ONEBIN: use single bin estimator */
  i = 0;
  if (cfg == NULL || 0 != cfgget(cfg, &i, "mbest_single_bin", "%u")) {
    fprintf(stderr, "assuming default: MB_ONEBIN = 0, key: mbest_single_bin\n");
  }
  if ( !(i == 0 || i == 1) ) {
    fprintf(stderr, "MB_ONEBIN: failed validation: i == 0 || i == 1\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  if (i) {
    mb->flags |= MB_ONEBIN;
  } else {
    mb->flags &= ~MB_ONEBIN;
  }

  /* MB_VERBOSE: being verbose */
  i = 1;
  if (cfg == NULL || 0 != cfgget(cfg, &i, "mbest_verbose", "%u")) {
    fprintf(stderr, "assuming default: MB_VERBOSE = 1, key: mbest_verbose\n");
  }
  if ( !(i == 0 || i == 1) ) {
    fprintf(stderr, "MB_VERBOSE: failed validation: i == 0 || i == 1\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  if (i) {
    mb->flags |= MB_VERBOSE;
  } else {
    mb->flags &= ~MB_VERBOSE;
  }

  /* MB_SBCORR: include energy fluctuation correction due to a small bin width for internal energy, etc */
  i = 1;
  if (cfg == NULL || 0 != cfgget(cfg, &i, "mbest_sbcorr", "%u")) {
    fprintf(stderr, "assuming default: MB_SBCORR = 1, key: mbest_sbcorr\n");
  }
  if ( !(i == 0 || i == 1) ) {
    fprintf(stderr, "MB_SBCORR: failed validation: i == 0 || i == 1\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  if (i) {
    mb->flags |= MB_SBCORR;
  } else {
    mb->flags &= ~MB_SBCORR;
  }

  /* js: lower  boundary of asym. beta windows for ehat (asym.) */
  if ((mb->js = calloc((mb->n + 2), sizeof(int))) == NULL) {
    fprintf(stderr, "no memory! var: mb->js, type: int\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n+1; i++)
    mb->js[i] = 0;
  /* jt: higher boundary of asym. beta windows for ehat (asym.) */
  if ((mb->jt = calloc((mb->n + 2), sizeof(int))) == NULL) {
    fprintf(stderr, "no memory! var: mb->jt, type: int\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n+1; i++)
    mb->jt[i] = 0;
  /* jset: lower  boundary of beta windows for et (usu. sym.), i - jset[i] == jtet[i] - (i+1) */
  if ((mb->jset = calloc((mb->n + 1), sizeof(int))) == NULL) {
    fprintf(stderr, "no memory! var: mb->jset, type: int\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++)
    mb->jset[i] = 0;
  /* jtet: higher boundary of beta windows for et (usu. sym.), i - jset[i] == jtet[i] - (i+1) */
  if ((mb->jtet = calloc((mb->n + 1), sizeof(int))) == NULL) {
    fprintf(stderr, "no memory! var: mb->jtet, type: int\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++)
    mb->jtet[i] = 0;
  /* setup the temperature windows (js, jt, jset, jtet) need to setup damp and symwin before calling! */
  if ( !(mb_mkwin(mb, mb->bwmod, mb->bwdel, mb->js, mb->jt, mb->jset, mb->jtet) == 0) ) {
    fprintf(stderr, "setup the temperature windows (js, jt, jset, jtet) need to setup damp and symwin before calling!\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  /* compute m as the maximal window span */
  mb->m = mb_maxwinspan(mb);
  /* nstrefresh: interval of recalculating et for all temperature */
  mb->nstrefresh = 10000;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->nstrefresh, "nstrefresh", "%d")) {
    fprintf(stderr, "assuming default: mb->nstrefresh = 10000, key: nstrefresh\n");
  }

  /* av_nstsave: interval of writing mbav and ze files */
  mb->av_nstsave = 10000;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->av_nstsave, "nstav", "%d")) {
    fprintf(stderr, "assuming default: mb->av_nstsave = 10000, key: nstav\n");
  }

  /* av_binary: use binary format in mbav file */
  mb->av_binary = 1;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->av_binary, "mbav_binary", "%d")) {
    fprintf(stderr, "assuming default: mb->av_binary = 1, key: mbav_binary\n");
  }

  /* av_file: name of mbav file */
  strcpy(buf, "mb");
  strcat(buf, suf);
  strcat(buf,".av");
  mb->av_file = ssdup(buf);
  if (cfg == NULL || 0 != cfgget(cfg, &mb->av_file, "mbav_file", "%s")) {
    fprintf(stderr, "assuming default: mb->av_file = \"mb.av\", key: mbav_file\n");
  }

  /* ze_file: name of ze file */
  strcpy(buf, "ZE");
  strcat(buf, suf);
  mb->ze_file = ssdup(buf);
  if (cfg == NULL || 0 != cfgget(cfg, &mb->ze_file, "ze_file", "%s")) {
    fprintf(stderr, "assuming default: mb->ze_file = \"ZE\", key: ze_file\n");
  }

  /* wze_reps: number of iterations before writing ze file */
  mb->wze_reps = 5;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->wze_reps, "mbest_wze_reps", "%d")) {
    fprintf(stderr, "assuming default: mb->wze_reps = 5, key: mbest_wze_reps\n");
  }

  /* vis: number of visits */
  if ((mb->vis = calloc((mb->n + 1), sizeof(double))) == NULL) {
    fprintf(stderr, "no memory! var: mb->vis, type: double\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++)
    mb->vis[i] = 0.0;
  /* totvis: total number of visits, number of tempering */
  mb->totvis = 0.0;
  /* winstot: total of sum.s over a multiple-bin temperature window */
  if ((mb->winstot = calloc((mb->n + 1), sizeof(double))) == NULL) {
    fprintf(stderr, "no memory! var: mb->winstot, type: double\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++)
    mb->winstot[i] = 0.0;
  /* lgv_dt: time step for the temperature Langevin eq */
  mb->lgv_dt = 1e-5;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->lgv_dt, "Tdt", "%lf")) {
    fprintf(stderr, "assuming default: mb->lgv_dt = 1e-5, key: Tdt\n");
  }

  /* lgv_dTmax: maximal amount of temperature change in a step */
  mb->lgv_dTmax = 25.0;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->lgv_dTmax, "dTmax", "%lf")) {
    fprintf(stderr, "assuming default: mb->lgv_dTmax = 25.0, key: dTmax\n");
  }

  /* lgv_rej: number of attempts of langevin equation trying to change beta too drastically */
  mb->lgv_rej = 0.0;
  /* lgv_tot: total number of attempts of using langevin equation */
  mb->lgv_tot = 0.0;
  /* regl: average within a bin first */
  mb->regl = 2;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->regl, "mbest_regularize", "%d")) {
    fprintf(stderr, "assuming default: mb->regl = 2, key: mbest_regularize\n");
  }

  /* fracmin: minimal allowable coefficient during left/right combination */
  mb->fracmin = 0.0;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->fracmin, "mbest_fracmin", "%lf")) {
    fprintf(stderr, "assuming default: mb->fracmin = 0.0, key: mbest_fracmin\n");
  }

  /* cvshiftmax: maximal fraction for shift energy fluct. if cv is monotonic, it should be 0.0, for ising model, it can restrain the magnitude */
  mb->cvshiftmax = 1.0;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->cvshiftmax, "mbest_cvshiftmax", "%lf")) {
    fprintf(stderr, "assuming default: mb->cvshiftmax = 1.0, key: mbest_cvshiftmax\n");
  }

  /* ens_exp: ensemble exponent of beta */
  mb->ens_exp = 1.0;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->ens_exp, "ensemble_factor", "%lf")) {
    fprintf(stderr, "assuming default: mb->ens_exp = 1.0, key: ensemble_factor\n");
  }

  /* multicanonical ensemble mode */
  mb->mode = 0;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->mode, "ensemble_mode", "%d"))
    fprintf(stderr, "assuming default: mb->mode = 1.0, key: ensemble_mode\n");
  
  /* default values */
  mb->beta0 = 0.5 * (mb->bmax + mb->bmin);
  mb->invsigma2 = 1.0;
  mb->c = 0.0;

  if(mb->mode == 1)
  {
    /* beta0 */
    if (cfg == NULL || 0 != cfgget(cfg, &mb->beta0, "ensemble_beta0", "%lf"))
      fprintf(stderr, "assuming default: mb->beta0 = 0.5 * (mb->bmax + mb->bmin), key: ensemble_beta0\n");
    if (mb->beta0 >= mb->bmax || mb->beta0 <= mb->bmin)
      fprintf(stderr, "WARNING: beta0 is not in the temperature range!\n");

    /* sigma */
    double sigma = 1.0;
    if (cfg == NULL || 0 != cfgget(cfg, &sigma, "ensemble_sigma", "%lf"))
      fprintf(stderr, "assuming default: mb->sigma = 1.0, key: ensemble_sigma\n");
    if (sigma == 0)
    {
      fprintf(stderr, "ERROR: sigma cannot be zero!\n");
      goto ERR;
    }
    mb->invsigma2 = 1.0/sigma/sigma;
  }
  else if(mb->mode == 2)
  {
    /* c */
    mb->c = 0.0;
    if (cfg == NULL || 0 != cfgget(cfg, &mb->c, "ensemble_c", "%lf"))
      fprintf(stderr, "assuming default: mb->c = 0.0, key: ensemble_c\n");
  }
  else if(mb->mode != 0)
  {
    fprintf(stderr, "invalid multicanonical ensemble parameter\n");
    goto ERR;
  }

  /* ens_w: array of ensemble weights at bin boundaries */
  if ((mb->ens_w = calloc((mb->n + 2), sizeof(double))) == NULL) {
    fprintf(stderr, "no memory! var: mb->ens_w, type: double\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n+1; i++)
    mb->ens_w[i] = 1.0/(mb_ensinvwf(mb, mb->barr[i], NULL, NULL));
  /* shk_base: current generic shrink amplitude */
  mb->shk_base = 0.0;
  /* shk_winadj: adjust shrink according to temperature window width */
  mb->shk_winadj = 1;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->shk_winadj, "shrink_mbin_adjust", "%d")) {
    fprintf(stderr, "assuming default: mb->shk_winadj = 1, key: shrink_mbin_adjust\n");
  }

  /* shk_max: initial and maximal shrink (adjusted) */
  mb->shk_max = 0.01;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->shk_max, "shrink0", "%lf")) {
    fprintf(stderr, "assuming default: mb->shk_max = 0.01, key: shrink0\n");
  }
  if ( !(mb->shk_max < 0.9 && mb->shk_max >= 0.0) ) {
    fprintf(stderr, "mb->shk_max: failed validation: mb->shk_max < 0.9 && mb->shk_max >= 0.0\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }

  /* shk_gauge: array used of modulation shrinking factors */
  if ((mb->shk_gauge = calloc((mb->n + 1), sizeof(double))) == NULL) {
    fprintf(stderr, "no memory! var: mb->shk_gauge, type: double\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++)
    mb->shk_gauge[i] = (mb_ensinvwf(mb, 0.5*(mb->barr[i] + mb->barr[i+1]), NULL, NULL) * mb->m) / (mb->jtet[i] - mb->jset[i]);
  /* shk_mode: 0: const, 1: amp/t, 2: amp/t^exp */
  mb->shk_mode = 1;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->shk_mode, "shrink_mode", "%d")) {
    fprintf(stderr, "assuming default: mb->shk_mode = 1, key: shrink_mode\n");
  }
  if ( !(mb->shk_mode >= 0 && mb->shk_mode <= 2) ) {
    fprintf(stderr, "mb->shk_mode: failed validation: mb->shk_mode >= 0 && mb->shk_mode <= 2\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }

  /* shk_min: minimal value for enforcing acc. sampling */
  mb->shk_min = 0.0;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->shk_min, "shrinkmin", "%lf")) {
    fprintf(stderr, "assuming default: mb->shk_min = 0.0, key: shrinkmin\n");
  }

  /* shk_stop: stop shrinking after this number of steps */
  mb->shk_stop = -1;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->shk_stop, "shrinkstop", "%d")) {
    fprintf(stderr, "assuming default: mb->shk_stop = -1, key: shrinkstop\n");
  }

  /* shk_amp: amp t^(-exp) */
  mb->shk_amp = 0.1;
  if (mb->shk_mode >= 1) {
    if (cfg == NULL || 0 != cfgget(cfg, &mb->shk_amp, "shrinkamp", "%lf")) {
      fprintf(stderr, "assuming default: mb->shk_amp = 0.1, key: shrinkamp\n");
    }
  }

  /* shk_exp: amp t^(-exp) */
  mb->shk_exp = 1.0;
  if (mb->shk_mode >= 2) {
    if (cfg == NULL || 0 != cfgget(cfg, &mb->shk_exp, "shrinkexp", "%lf")) {
      fprintf(stderr, "assuming default: mb->shk_exp = 1.0, key: shrinkexp\n");
    }
  }

  /* lnz: logarithm of the partition function */
  if ((mb->lnz = calloc((mb->n + 2), sizeof(double))) == NULL) {
    fprintf(stderr, "no memory! var: mb->lnz, type: double\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n+1; i++)
    mb->lnz[i] = 0.0;
  /* ehat: internal energy */
  if ((mb->ehat = calloc((mb->n + 2), sizeof(double))) == NULL) {
    fprintf(stderr, "no memory! var: mb->ehat, type: double\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n+1; i++)
    mb->ehat[i] = 0.0;
  /* cvhat: heat capacity */
  if ((mb->cvhat = calloc((mb->n + 2), sizeof(double))) == NULL) {
    fprintf(stderr, "no memory! var: mb->cvhat, type: double\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n+1; i++)
    mb->cvhat[i] = 0.0;
  /* et: bin-averaged internal energy */
  if ((mb->et = calloc((mb->n + 1), sizeof(double))) == NULL) {
    fprintf(stderr, "no memory! var: mb->et, type: double\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++)
    mb->et[i] = 0.0;
  /* imbal: |a+ - a-| / (a+ + a-) for left-right combination */
  if ((mb->imbal = calloc((mb->n + 2), sizeof(double))) == NULL) {
    fprintf(stderr, "no memory! var: mb->imbal, type: double\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n+1; i++)
    mb->imbal[i] = 0.0;
  /* haset: current et[i] is reasonably good */
  if ((mb->haset = calloc((mb->n + 1), sizeof(int))) == NULL) {
    fprintf(stderr, "no memory! var: mb->haset, type: int\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++)
    mb->haset[i] = 0;
  /* qua: bits represent whether estimated values are unbiased */
  if ((mb->qua = calloc((mb->n + 2), sizeof(unsigned))) == NULL) {
    fprintf(stderr, "no memory! var: mb->qua, type: unsigned\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n+1; i++)
    mb->qua[i] = 0;
  /* ampf: currently amplification factor for adaptive averaging */
  if ((mb->ampf = calloc((mb->n + 1), sizeof(double))) == NULL) {
    fprintf(stderr, "no memory! var: mb->ampf, type: double\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++)
    mb->ampf[i] = 1.0;
  /* sums: normal data */
  if ((mb->sums = calloc((mb->n + 1), sizeof(sm_t))) == NULL) {
    fprintf(stderr, "no memory! var: mb->sums, type: sm_t\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++) {
    if (0 != sm_cfgopen_low(mb->sums+i, cfg)) {
      fprintf(stderr, "failed to initialize mb->sums[%d]\n\n", i);
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }
  /* has_xsums */
  mb->has_xsums = mb->flags & MB_DAMP;
  /* cnt_int: number of additional integer variables to be written to binary file */
  mb->cnt_int = 0;
  /* cnt_dbl: number of additional double variables to be written to binary file */
  mb->cnt_dbl = 5;
  /* idxcnt: index count */
  if ((mb->idxcnt = calloc((mb->n + 1), sizeof(int))) == NULL) {
    fprintf(stderr, "no memory! var: mb->idxcnt, type: int\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++)
    mb->idxcnt[i] = 0;
  /* compute idxcnt, and adjust m to the maximal value among idxcnt[j].  the call is needed before allocating midx */
  mb->m = mb_maxidxcnt(mb);
  /* midx: index look-up table */
  if ((mb->midx = calloc((mb->n*mb->m + 1), sizeof(int))) == NULL) {
    fprintf(stderr, "no memory! var: mb->midx, type: int\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n*mb->m; i++)
    mb->midx[i] = 0;
  /* setup indices */
  if ( !(mb_mkidx(mb) == 0) ) {
    fprintf(stderr, "setup indices\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  /* xsums: multiple-bin damping data */
  if ((mb->xsums = calloc((mb->n*mb->m + 1), sizeof(sm_t))) == NULL) {
    fprintf(stderr, "no memory! var: mb->xsums, type: sm_t\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n*mb->m; i++) {
    if (0 != sm_cfgopen_low(mb->xsums+i, cfg)) {
      fprintf(stderr, "failed to initialize mb->xsums[%d]\n\n", i);
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }
  /* eh_mode: 0: disable; 1: simple histogram */
  mb->eh_mode = 0;
  if (cfg == NULL || 0 != cfgget(cfg, &mb->eh_mode, "ehist_mode", "%d")) {
    fprintf(stderr, "assuming default: mb->eh_mode = 0, key: ehist_mode\n");
  }
  if ( !(mb->eh_mode == 0 || mb->eh_mode == 1) ) {
    fprintf(stderr, "mb->eh_mode: failed validation: mb->eh_mode == 0 || mb->eh_mode == 1\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }

  /* eh_skip: interval of reconstructing energy histograms */
  mb->eh_skip = 10;
  if (mb->eh_mode) {
    if (cfg == NULL || 0 != cfgget(cfg, &mb->eh_skip, "ehist_skip", "%d")) {
      fprintf(stderr, "assuming default: mb->eh_skip = 10, key: ehist_skip\n");
    }
    if ( !(mb->eh_skip > 0) ) {
      fprintf(stderr, "mb->eh_skip: failed validation: mb->eh_skip > 0\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }

  /* eh_bwmod: 0: d(beta) 1: dT/T  2: d(kT) */
  mb->eh_bwmod = 1;
  if (mb->eh_mode) {
    if (cfg == NULL || 0 != cfgget(cfg, &mb->eh_bwmod, "ehist_mbin_mode", "%d")) {
      fprintf(stderr, "assuming default: mb->eh_bwmod = 1, key: ehist_mbin_mode\n");
    }
    if ( !(mb->eh_bwmod >= 0 && mb->eh_bwmod <= 2) ) {
      fprintf(stderr, "mb->eh_bwmod: failed validation: mb->eh_bwmod >= 0 && mb->eh_bwmod <= 2\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }

  /* eh_bwdel: delta lnT */
  mb->eh_bwdel = 0.05;
  if (mb->eh_mode && mb->eh_bwmod == 1) {
    if (cfg != NULL && 0 != cfgget(cfg, &mb->eh_bwdel, "ehist_delta_lnT", "%lf")) {
      fprintf(stderr, "missing var: mb->eh_bwdel, key: ehist_delta_lnT, fmt: %%lf\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    if ( !(mb->eh_bwdel > mb->bdel/pow(mb->bmin, 1.0)) ) {
      fprintf(stderr, "mb->eh_bwdel: failed validation: mb->eh_bwdel > mb->bdel/pow(mb->bmin, 1.0)\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }

  /* eh_bwdel: delta beta */
  if (mb->eh_mode && mb->eh_bwmod == 0) {
    mb->eh_bwdel = 0.02;
    if (cfg != NULL && 0 != cfgget(cfg, &mb->eh_bwdel, "ehist_delta_beta", "%lf")) {
      fprintf(stderr, "missing var: mb->eh_bwdel, key: ehist_delta_beta, fmt: %%lf\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    if ( !(mb->eh_bwdel > mb->bdel/pow(mb->bmin, 0.0)) ) {
      fprintf(stderr, "mb->eh_bwdel: failed validation: mb->eh_bwdel > mb->bdel/pow(mb->bmin, 0.0)\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }

  /* eh_bwdel: delta kT */
  if (mb->eh_mode && mb->eh_bwmod == 2) {
    mb->eh_bwdel = 0.10;
    if (cfg != NULL && 0 != cfgget(cfg, &mb->eh_bwdel, "ehist_delta_kT", "%lf")) {
      fprintf(stderr, "missing var: mb->eh_bwdel, key: ehist_delta_kT, fmt: %%lf\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    if ( !(mb->eh_bwdel > mb->bdel/pow(mb->bmin, 2.0)) ) {
      fprintf(stderr, "mb->eh_bwdel: failed validation: mb->eh_bwdel > mb->bdel/pow(mb->bmin, 2.0)\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }

  /* eh_min: minimal energy */
  mb->eh_min = -12.6e4;
  if (mb->eh_mode) {
    if (cfg == NULL || 0 != cfgget(cfg, &mb->eh_min, "ehist_min", "%lf")) {
      fprintf(stderr, "assuming default: mb->eh_min = -12.6e4, key: ehist_min\n");
    }
  }

  /* eh_max: maximal energy */
  mb->eh_max = -9.0e4;
  if (mb->eh_mode) {
    if (cfg == NULL || 0 != cfgget(cfg, &mb->eh_max, "ehist_max", "%lf")) {
      fprintf(stderr, "assuming default: mb->eh_max = -9.0e4, key: ehist_max\n");
    }
    if ( !(mb->eh_max > mb->eh_min) ) {
      fprintf(stderr, "mb->eh_max: failed validation: mb->eh_max > mb->eh_min\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }

  /* eh_del: energy bin size */
  mb->eh_del = 20.0;
  if (mb->eh_mode) {
    if (cfg == NULL || 0 != cfgget(cfg, &mb->eh_del, "ehist_del", "%lf")) {
      fprintf(stderr, "assuming default: mb->eh_del = 20.0, key: ehist_del\n");
    }
    if ( !(mb->eh_del > 0) ) {
      fprintf(stderr, "mb->eh_del: failed validation: mb->eh_del > 0\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }

  /* eh_cnt: number of energy bins */
  mb->eh_cnt = (int)((mb->eh_max-mb->eh_min)/mb->eh_del - 1e-5 + 1);
  /* eh_binary: binary format for ehist file */
  mb->eh_binary = 1;
  if (mb->eh_mode) {
    if (cfg == NULL || 0 != cfgget(cfg, &mb->eh_binary, "ehist_binary", "%d")) {
      fprintf(stderr, "assuming default: mb->eh_binary = 1, key: ehist_binary\n");
    }
  }

  /* eh_nstsave: interval of writing histogrm files */
  mb->eh_nstsave = 100000;
  if (mb->eh_mode) {
    if (cfg == NULL || 0 != cfgget(cfg, &mb->eh_nstsave, "nsthist", "%d")) {
      fprintf(stderr, "assuming default: mb->eh_nstsave = 100000, key: nsthist\n");
    }
  }

  /* eh_file: name of ehist file */
  mb->eh_file = NULL;
  if (mb->eh_mode) {
    strcpy(buf, "hist");
    strcat(buf, suf);
    strcat(buf,".bin");
    mb->eh_file = ssdup(buf);
    if (cfg == NULL || 0 != cfgget(cfg, &mb->eh_file, "ehist_file", "%s")) {
      fprintf(stderr, "assuming default: mb->eh_file = \"hist.bin\", key: ehist_file\n");
    }
  }

  /* eh_rfile: name of reconstructed energy histogram */
  mb->eh_rfile = NULL;
  if (mb->eh_mode) {
    strcpy(buf, "HMB");
    strcat(buf, suf);
    mb->eh_rfile = ssdup(buf);
    if (cfg == NULL || 0 != cfgget(cfg, &mb->eh_rfile, "ehist_mbin_file", "%s")) {
      fprintf(stderr, "assuming default: mb->eh_rfile = \"HMB\", key: ehist_mbin_file\n");
    }
  }

  /* eh_his: energy histogram data */
  mb->eh_his = NULL;
  if (mb->eh_mode) {
    if ((mb->eh_his = calloc((mb->n*mb->eh_cnt + 1), sizeof(double))) == NULL) {
      fprintf(stderr, "no memory! var: mb->eh_his, type: double\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      exit(1);
    }
    for (i = 0; i < mb->n*mb->eh_cnt; i++)
      mb->eh_his[i] = 0.0;
  }
  /* eh_recon: temporary space for reconstructing histogram */
  mb->eh_recon = NULL;
  if (mb->eh_mode) {
    if ((mb->eh_recon = calloc((mb->eh_cnt + 1), sizeof(double))) == NULL) {
      fprintf(stderr, "no memory! var: mb->eh_recon, type: double\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      exit(1);
    }
    for (i = 0; i < mb->eh_cnt; i++)
      mb->eh_recon[i] = 0.0;
  }
  /* eh_is: indices for temperature windows (lower) */
  mb->eh_is = NULL;
  if (mb->eh_mode) {
    if ((mb->eh_is = calloc((mb->n + 2), sizeof(int))) == NULL) {
      fprintf(stderr, "no memory! var: mb->eh_is, type: int\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      exit(1);
    }
    for (i = 0; i < mb->n + 1; i++)
      mb->eh_is[i] = 0;
  }
  /* eh_it: indices for temperature windows (higher) */
  mb->eh_it = NULL;
  if (mb->eh_mode) {
    if ((mb->eh_it = calloc((mb->n + 2), sizeof(int))) == NULL) {
      fprintf(stderr, "no memory! var: mb->eh_it, type: int\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      exit(1);
    }
    for (i = 0; i < mb->n + 1; i++)
      mb->eh_it[i] = 0;
    if ( !(0 == mb_mkwin(mb, mb->eh_bwmod, mb->eh_bwdel, mb->eh_is, mb->eh_it, NULL, NULL)) ) {
      fprintf(stderr, "mb->eh_it: failed validation: 0 == mb_mkwin(mb, mb->eh_bwmod, mb->eh_bwdel, mb->eh_is, mb->eh_it, NULL, NULL)\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }
  /* MB_EH_ADDAHALF: add a half energy bin width in output */
  if (mb->eh_mode) {
    i = 1;
    if (cfg == NULL || 0 != cfgget(cfg, &i, "ehist_addahalf", "%u")) {
      fprintf(stderr, "assuming default: MB_EH_ADDAHALF = 1, key: ehist_addahalf\n");
    }
    if ( !(i == 0 || i == 1) ) {
      fprintf(stderr, "MB_EH_ADDAHALF: failed validation: i == 0 || i == 1\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    if (i) {
      mb->flags |= MB_EH_ADDAHALF;
    } else {
      mb->flags &= ~MB_EH_ADDAHALF;
    }
  }

  /* MB_EH_KEEPEDGE: keep zero edge at sides */
  if (mb->eh_mode) {
    i = 0;
    if (cfg == NULL || 0 != cfgget(cfg, &i, "ehist_keepedge", "%u")) {
      fprintf(stderr, "assuming default: MB_EH_KEEPEDGE = 0, key: ehist_keepedge\n");
    }
    if ( !(i == 0 || i == 1) ) {
      fprintf(stderr, "MB_EH_KEEPEDGE: failed validation: i == 0 || i == 1\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    if (i) {
      mb->flags |= MB_EH_KEEPEDGE;
    } else {
      mb->flags &= ~MB_EH_KEEPEDGE;
    }
  }

  /* MB_EH_NOZEROES: do not output zeroes */
  if (mb->eh_mode) {
    i = 0;
    if (cfg == NULL || 0 != cfgget(cfg, &i, "ehist_nozeroes", "%u")) {
      fprintf(stderr, "assuming default: MB_EH_NOZEROES = 0, key: ehist_nozeroes\n");
    }
    if ( !(i == 0 || i == 1) ) {
      fprintf(stderr, "MB_EH_NOZEROES: failed validation: i == 0 || i == 1\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    if (i) {
      mb->flags |= MB_EH_NOZEROES;
    } else {
      mb->flags &= ~MB_EH_NOZEROES;
    }
  }

  return 0;
ERR:
  return -1;
}

mb_t *mb_cfgopen(cfgdata_t *cfg, double bmin, double bmax, char suffix)
{
  mb_t *mb;

  /* allocate memory for mb_t */
  if ((mb = calloc(1, sizeof(*mb))) == NULL) {
    fprintf(stderr, "no memory for mb_t (size %u)\n", (unsigned) sizeof(*mb));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  
  /* call low level function */
  if (0 != mb_cfgopen_low(mb, cfg, bmin, bmax, suffix)) {
    fprintf(stderr, "mb_t: error while reading configuration file\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    free(mb);
    return NULL;
  }

  return mb;
}

void mb_close_low(mb_t *mb)
{
  int i;

  if (mb->barr      != NULL) free(mb->barr);
  if (mb->js        != NULL) free(mb->js);
  if (mb->jt        != NULL) free(mb->jt);
  if (mb->jset      != NULL) free(mb->jset);
  if (mb->jtet      != NULL) free(mb->jtet);
  if (mb->av_file   != NULL) ssdelete(mb->av_file);
  if (mb->ze_file   != NULL) ssdelete(mb->ze_file);
  if (mb->vis       != NULL) free(mb->vis);
  if (mb->winstot   != NULL) free(mb->winstot);
  if (mb->ens_w     != NULL) free(mb->ens_w);
  if (mb->shk_gauge != NULL) free(mb->shk_gauge);
  if (mb->lnz       != NULL) free(mb->lnz);
  if (mb->ehat      != NULL) free(mb->ehat);
  if (mb->cvhat     != NULL) free(mb->cvhat);
  if (mb->et        != NULL) free(mb->et);
  if (mb->imbal     != NULL) free(mb->imbal);
  if (mb->haset     != NULL) free(mb->haset);
  if (mb->qua       != NULL) free(mb->qua);
  if (mb->ampf      != NULL) free(mb->ampf);
  if (mb->sums != NULL) {
    for (i = 0; i < mb->n; i++) {
      sm_close_low(mb->sums + i);
    }
    free(mb->sums);
  }
  if (mb->idxcnt    != NULL) free(mb->idxcnt);
  if (mb->midx      != NULL) free(mb->midx);
  if (mb->xsums != NULL) {
    for (i = 0; i < mb->n*mb->m; i++) {
      sm_close_low(mb->xsums + i);
    }
    free(mb->xsums);
  }
  if (mb->eh_mode) {
    if (mb->eh_file   != NULL) ssdelete(mb->eh_file);
    if (mb->eh_rfile  != NULL) ssdelete(mb->eh_rfile);
    if (mb->eh_his    != NULL) free(mb->eh_his);
    if (mb->eh_recon  != NULL) free(mb->eh_recon);
    if (mb->eh_is     != NULL) free(mb->eh_is);
    if (mb->eh_it     != NULL) free(mb->eh_it);
  }
  memset(mb, 0, sizeof(*mb));
}

void mb_clear(mb_t *mb)
{
  int i;

  mb->lgv_rej = 0.0;
  mb->lgv_tot = 0.0;
  for (i = 0; i < mb->n+1; i++)
    mb->lnz[i] = 0.0;
  for (i = 0; i < mb->n+1; i++)
    mb->ehat[i] = 0.0;
  for (i = 0; i < mb->n+1; i++)
    mb->cvhat[i] = 0.0;
  for (i = 0; i < mb->n; i++)
    mb->et[i] = 0.0;
  for (i = 0; i < mb->n+1; i++)
    mb->imbal[i] = 0.0;
  for (i = 0; i < mb->n; i++)
    mb->haset[i] = 0;
  for (i = 0; i < mb->n+1; i++)
    mb->qua[i] = 0;
  for (i = 0; i < mb->n; i++)
    mb->ampf[i] = 1.0;
  if (mb->sums != NULL) {
    for (i = 0; i < mb->n; i++)
      sm_clear(mb->sums+i);
  }
  if (mb->xsums != NULL) {
    for (i = 0; i < mb->n*mb->m; i++)
      sm_clear(mb->xsums+i);
  }
}

int mb_readbin_low(mb_t *mb, FILE *fp, int ver, int endn)
{
  double lgv_rate;
  int i;
  int j;
  int itmp;
  int jmax;
  int jmin;

  if (mb == NULL) {
    fprintf(stderr, "passing null pointer to mb_readbin_low\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    return -1;
  }
  /* clear data before reading */
  mb_clear(mb);

  /* n: number of temperature bins */
  if (endn_fread(&itmp, sizeof(itmp), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading itmp\n");
    goto ERR;
  }
  if (itmp != mb->n) {
    fprintf(stderr, "mb->n mismatch, expect: %d, read: %d, pos: %#lx\n",
        mb->n, itmp, (unsigned long) ftell(fp));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  /* m: maximal number of bins in a window */
  if (endn_fread(&itmp, sizeof(itmp), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading itmp\n");
    goto ERR;
  }
  if (itmp != mb->m) {
    fprintf(stderr, "mb->m mismatch, expect: %d, read: %d, pos: %#lx\n",
        mb->m, itmp, (unsigned long) ftell(fp));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  /* order: order, should be 1 */
  if (endn_fread(&mb->order, sizeof(mb->order), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading mb->order\n");
    goto ERR;
  }
  if ( !(mb->order == 1) ) {
    fprintf(stderr, "mb->order: failed validation: mb->order == 1\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  /* flags: combination of flags */
  if (endn_fread(&mb->flags, sizeof(mb->flags), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading mb->flags\n");
    goto ERR;
  }
  /* has_xsums */
  if (endn_fread(&mb->has_xsums, sizeof(mb->has_xsums), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading mb->has_xsums\n");
    goto ERR;
  }
  /* cnt_int: number of additional integer variables to be written to binary file */
  if (endn_fread(&itmp, sizeof(itmp), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading itmp\n");
    goto ERR;
  }
  if (itmp != mb->cnt_int) {
    fprintf(stderr, "mb->cnt_int mismatch, expect: %d, read: %d, pos: %#lx\n",
        mb->cnt_int, itmp, (unsigned long) ftell(fp));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  /* cnt_dbl: number of additional double variables to be written to binary file */
  if (endn_fread(&itmp, sizeof(itmp), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading itmp\n");
    goto ERR;
  }
  if (itmp != mb->cnt_dbl) {
    fprintf(stderr, "mb->cnt_dbl mismatch, expect: %d, read: %d, pos: %#lx\n",
        mb->cnt_dbl, itmp, (unsigned long) ftell(fp));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  /* beta: current value of beta */
  if (endn_fread(&mb->beta, sizeof(mb->beta), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading mb->beta\n");
    goto ERR;
  }
  if ( !(mb->beta >= mb->bmin && mb->beta <= mb->bmax) ) {
    fprintf(stderr, "mb->beta: failed validation: mb->beta >= mb->bmin && mb->beta <= mb->bmax\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  mb->beta = (mb->beta >= mb->bmax - 1e-5) ? (mb->bmax - 1e-5) : mb->beta;
  mb->beta = (mb->beta <= mb->bmin + 1e-5) ? (mb->bmin + 1e-5) : mb->beta;
  /* totvis: total number of visits, number of tempering */
  if (endn_fread(&mb->totvis, sizeof(mb->totvis), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading mb->totvis\n");
    goto ERR;
  }
  if ( !(mb->totvis > 0) ) {
    fprintf(stderr, "mb->totvis: failed validation: mb->totvis > 0\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  /* shk_base: current generic shrink amplitude */
  if (endn_fread(&mb->shk_base, sizeof(mb->shk_base), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading mb->shk_base\n");
    goto ERR;
  }
  /* lgv_rate */
  lgv_rate = 0.0;
  if (endn_fread(&lgv_rate, sizeof(lgv_rate), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading lgv_rate\n");
    goto ERR;
  }
  /* lgv_tot: total number of attempts of using langevin equation */
  if (endn_fread(&mb->lgv_tot, sizeof(mb->lgv_tot), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading mb->lgv_tot\n");
    goto ERR;
  }
  if ( !((mb->lgv_rej=mb->lgv_tot*lgv_rate) >= 0.0) ) {
    fprintf(stderr, "mb->lgv_tot: failed validation: (mb->lgv_rej=mb->lgv_tot*lgv_rate) >= 0.0\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  /* barr: temperature array */
  if ((mb->n > 0)
    && (endn_fread(mb->barr, sizeof(*(mb->barr)), mb->n, fp, endn) != (size_t) (mb->n))) {
    fprintf(stderr, "error in reading mb->barr, n = mb->n(%d)\n", mb->n);
    goto ERR;
  }
  /* et: bin-averaged internal energy */
  if ((mb->n > 0)
    && (endn_fread(mb->et, sizeof(*(mb->et)), mb->n, fp, endn) != (size_t) (mb->n))) {
    fprintf(stderr, "error in reading mb->et, n = mb->n(%d)\n", mb->n);
    goto ERR;
  }
  /* sums: normal data */
  if (mb->sums == NULL) {
    fprintf(stderr, "mb->sums is null\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++) {
    if (endn_fread(&itmp, sizeof(itmp), 1, fp, endn) != 1) {
      fprintf(stderr, "error in reading itmp\n");
      goto ERR;
    }
    if (itmp != i) {
      fprintf(stderr, "i mismatch, expect: %d, read: %d, pos: %#lx\n",
          i, itmp, (unsigned long) ftell(fp));
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    if (0 != sm_readbin_low(mb->sums+i, fp, ver, endn)) {
      fprintf(stderr, "error reading object mb->sums+i\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }
  /* xsums: multiple-bin damping data */
  if (mb->xsums == NULL) {
    fprintf(stderr, "mb->xsums is null\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++) {
    if (endn_fread(&itmp, sizeof(itmp), 1, fp, endn) != 1) {
      fprintf(stderr, "error in reading itmp\n");
      goto ERR;
    }
    if (itmp != i) {
      fprintf(stderr, "i mismatch, expect: %d, read: %d, pos: %#lx\n",
          i, itmp, (unsigned long) ftell(fp));
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    jmin = mb->js[i];
    jmax = mb->jt[i+1];
    for (j = jmin; j < jmax; j++) {
      if (endn_fread(&itmp, sizeof(itmp), 1, fp, endn) != 1) {
        fprintf(stderr, "error in reading itmp\n");
        goto ERR;
      }
      if (itmp != j) {
        fprintf(stderr, "j mismatch, expect: %d, read: %d, pos: %#lx\n",
            j, itmp, (unsigned long) ftell(fp));
        fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
        goto ERR;
      }
      if (0 != sm_readbin_low(mb->xsums+i*mb->m-jmin+j, fp, ver, endn)) {
        fprintf(stderr, "error reading object mb->xsums+i*mb->m-jmin+j\n");
        fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
        goto ERR;
      }
    }
  }
  return 0;
ERR:
  mb_clear(mb);
  return -1;
}

int mb_readbin(mb_t *mb, const char *fname, int *pver)
{
  FILE *fp;
  int ver;
  int itmp;
  int i;
  int endn;

  if ((fp = fopen(fname, "rb")) == NULL) {
    fprintf(stderr, "cannot read binary file [%s].\n", fname);
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    return -1;
  }

  /* determine file endian */
  if ((endn = endn_rmatchi(&itmp, sizeof(int), fp)) < 0) {
    fprintf(stderr, "itmp 0x%X cannot match sizeof(int) 0x%X\n",
        (unsigned) itmp, (unsigned) sizeof(int));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  if (endn_fread(&itmp, sizeof(itmp), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading itmp\n");
    goto ERR;
  }
  if (itmp != (int) sizeof(double)) {
    fprintf(stderr, "(int) sizeof(double) mismatch, expect: %d, read: %d, pos: %#lx\n",
        (int) sizeof(double), itmp, (unsigned long) ftell(fp));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  if (endn_fread(&ver, sizeof(ver), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading ver\n");
    goto ERR;
  }
  if (pver != NULL) *pver = ver;

  /* call low level read function for members */
  i = mb_readbin_low(mb, fp, ver, endn);
  fclose(fp);
  return i;
ERR:
  fclose(fp);
  return -1;
}

int mb_writebin_low(mb_t *mb, FILE *fp, int ver)
{
  int i;
  int jmax;
  int jmin;
  int j;
  double lgv_rate;

  if (mb == NULL) {
    fprintf(stderr, "passing null pointer to mb_writebin_low\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    return -1;
  }
  /* n: number of temperature bins */
  if (endn_fwrite(&mb->n, sizeof(mb->n), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->n\n");
    goto ERR;
  }
  /* m: maximal number of bins in a window */
  if (endn_fwrite(&mb->m, sizeof(mb->m), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->m\n");
    goto ERR;
  }
  /* order: order, should be 1 */
  if (endn_fwrite(&mb->order, sizeof(mb->order), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->order\n");
    goto ERR;
  }
  /* flags: combination of flags */
  if (endn_fwrite(&mb->flags, sizeof(mb->flags), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->flags\n");
    goto ERR;
  }
  /* has_xsums */
  if (endn_fwrite(&mb->has_xsums, sizeof(mb->has_xsums), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->has_xsums\n");
    goto ERR;
  }
  /* cnt_int: number of additional integer variables to be written to binary file */
  if (endn_fwrite(&mb->cnt_int, sizeof(mb->cnt_int), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->cnt_int\n");
    goto ERR;
  }
  /* cnt_dbl: number of additional double variables to be written to binary file */
  if (endn_fwrite(&mb->cnt_dbl, sizeof(mb->cnt_dbl), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->cnt_dbl\n");
    goto ERR;
  }
  /* beta: current value of beta */
  if (endn_fwrite(&mb->beta, sizeof(mb->beta), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->beta\n");
    goto ERR;
  }
  /* totvis: total number of visits, number of tempering */
  if (endn_fwrite(&mb->totvis, sizeof(mb->totvis), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->totvis\n");
    goto ERR;
  }
  /* shk_base: current generic shrink amplitude */
  if (endn_fwrite(&mb->shk_base, sizeof(mb->shk_base), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->shk_base\n");
    goto ERR;
  }
  /* lgv_rate */
  lgv_rate = (mb->lgv_tot > 1.0) ? (mb->lgv_rej/mb->lgv_tot) : 0.0;
  if (endn_fwrite(&lgv_rate, sizeof(lgv_rate), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing lgv_rate\n");
    goto ERR;
  }
  /* lgv_tot: total number of attempts of using langevin equation */
  if (endn_fwrite(&mb->lgv_tot, sizeof(mb->lgv_tot), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->lgv_tot\n");
    goto ERR;
  }
  /* barr: temperature array */
  if ((mb->n > 0)
    && (endn_fwrite(mb->barr, sizeof(*(mb->barr)), mb->n, fp, 1) != (size_t) (mb->n))) {
    fprintf(stderr, "error in writing mb->barr, n = mb->n(%d)\n", mb->n);
    goto ERR;
  }
  /* et: bin-averaged internal energy */
  if ((mb->n > 0)
    && (endn_fwrite(mb->et, sizeof(*(mb->et)), mb->n, fp, 1) != (size_t) (mb->n))) {
    fprintf(stderr, "error in writing mb->et, n = mb->n(%d)\n", mb->n);
    goto ERR;
  }
  /* sums: normal data */
  if (mb->sums == NULL) {
    fprintf(stderr, "mb->sums is null\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++) {
    if (endn_fwrite(&i, sizeof(i), 1, fp, 1) != 1) {
      fprintf(stderr, "error in writing i\n");
      goto ERR;
    }
    if (0 != sm_writebin_low(mb->sums+i, fp, ver)) {
      fprintf(stderr, "error writing object mb->sums+i\n");
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
  }
  /* xsums: multiple-bin damping data */
  mb_normalize(mb, -1);
  if (mb->xsums == NULL) {
    fprintf(stderr, "mb->xsums is null\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  for (i = 0; i < mb->n; i++) {
    if (endn_fwrite(&i, sizeof(i), 1, fp, 1) != 1) {
      fprintf(stderr, "error in writing i\n");
      goto ERR;
    }
    jmin = mb->js[i];
    jmax = mb->jt[i+1];
    for (j = jmin; j < jmax; j++) {
      if (endn_fwrite(&j, sizeof(j), 1, fp, 1) != 1) {
        fprintf(stderr, "error in writing j\n");
        goto ERR;
      }
      if (0 != sm_writebin_low(mb->xsums+i*mb->m-jmin+j, fp, ver)) {
        fprintf(stderr, "error writing object mb->xsums+i*mb->m-jmin+j\n");
        fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
        goto ERR;
      }
    }
  }
  return 0;
ERR:
  return -1;
}

int mb_writebin(mb_t *mb, const char *fname, int ver)
{
  FILE *fp;
  int i;
  int size;

  if ((fp = fopen(fname, "wb")) == NULL) {
    fprintf(stderr, "cannot write binary file [%s].\n", fname);
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    return -1;
  }

  size = (int) sizeof(int);
  if (endn_fwrite(&size, sizeof(size), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing size\n");
    goto ERR;
  }
  size = (int) sizeof(double);
  if (endn_fwrite(&size, sizeof(size), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing size\n");
    goto ERR;
  }
  if (endn_fwrite(&ver, sizeof(ver), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing ver\n");
    goto ERR;
  }

  /* call low level write function for members */
  i = mb_writebin_low(mb, fp, ver);
  fclose(fp);
  return i;
ERR:
  fclose(fp);
  return -1;
}

void mb_eh_clear(mb_t *mb)
{
  int i;

  if ( !(mb->eh_mode > 0) ) return ;
  if ( !(mb->eh_mode == 1) ) {
    fprintf(stderr, "mb_eh_clear: failed validation: mb->eh_mode == 1\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  if (mb->eh_mode) {
    for (i = 0; i < mb->n*mb->eh_cnt; i++)
      mb->eh_his[i] = 0.0;
    for (i = 0; i < mb->eh_cnt; i++)
      mb->eh_recon[i] = 0.0;
  }
}

int mb_eh_readbin_low(mb_t *mb, FILE *fp, int ver, int endn)
{
  int i;
  int itmp;
  double *pd;
  int jmin;
  double dtmp;
  int size;

  if (mb == NULL) {
    fprintf(stderr, "passing null pointer to mb_eh_readbin_low\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    return -1;
  }
  /* clear data before reading */
  mb_eh_clear(mb);

  /* n: number of temperature bins */
  if (endn_fread(&itmp, sizeof(itmp), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading itmp\n");
    goto ERR;
  }
  if (itmp != mb->n) {
    fprintf(stderr, "mb->n mismatch, expect: %d, read: %d, pos: %#lx\n",
        mb->n, itmp, (unsigned long) ftell(fp));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  /* eh_cnt: number of energy bins */
  if (endn_fread(&itmp, sizeof(itmp), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading itmp\n");
    goto ERR;
  }
  if (itmp != mb->eh_cnt) {
    fprintf(stderr, "mb->eh_cnt mismatch, expect: %d, read: %d, pos: %#lx\n",
        mb->eh_cnt, itmp, (unsigned long) ftell(fp));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  /* eh_min: minimal energy */
  if (endn_fread(&dtmp, sizeof(dtmp), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading dtmp\n");
    goto ERR;
  }
  if (fabs(dtmp - mb->eh_min) > 1e-5) {
    fprintf(stderr, "mb->eh_min mismatch, expect: %g, read: %g, pos: %#lx\n",
        mb->eh_min, dtmp, (unsigned long) ftell(fp));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  /* eh_del: energy bin size */
  if (endn_fread(&dtmp, sizeof(dtmp), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading dtmp\n");
    goto ERR;
  }
  if (fabs(dtmp - mb->eh_del) > 1e-5) {
    fprintf(stderr, "mb->eh_del mismatch, expect: %g, read: %g, pos: %#lx\n",
        mb->eh_del, dtmp, (unsigned long) ftell(fp));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  if ( !(mb->eh_del > 0) ) {
    fprintf(stderr, "mb->eh_del: failed validation: mb->eh_del > 0\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  /* eh_his: energy histogram data */
  for (i = 0; i < mb->n; i++) {
    if (endn_fread(&itmp, sizeof(itmp), 1, fp, endn) != 1) {
      fprintf(stderr, "error in reading itmp\n");
      if (feof(fp)) break;
    }
    if (itmp > i && itmp < mb->n) {
      i = itmp;
    } else if (i != itmp) {
      fprintf(stderr, "mb->eh_his bad major index, i: %d, read %d\n", i, itmp);
      goto ERR;
    }
    pd = mb->eh_his + i * mb->eh_cnt;
    if (endn_fread(&jmin, sizeof(jmin), 1, fp, endn) != 1) {
      fprintf(stderr, "error in reading jmin\n");
      goto ERR;
    }
    if ( !(jmin >= 0 && jmin < mb->eh_cnt) ) {
      fprintf(stderr, "mb->eh_his: base index %d out of boudary [0, mb->eh_cnt=%d)\n",
          jmin, mb->eh_cnt);
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    if (endn_fread(&size, sizeof(size), 1, fp, endn) != 1) {
      fprintf(stderr, "error in reading size\n");
      goto ERR;
    }
    if ( !(size > 0 && jmin + size <= mb->eh_cnt) ) {
      fprintf(stderr, "mb->eh_his: invalid size %d, jmin=%d, [0, mb->eh_cnt=%d)\n",
          size, jmin, mb->eh_cnt);
      fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
      goto ERR;
    }
    if ((size > 0)
      && (endn_fread(pd+jmin, sizeof(*(pd+jmin)), size, fp, endn) != (size_t) (size))) {
      fprintf(stderr, "error in reading pd+jmin, n = size(%d)\n", size);
      goto ERR;
    }
  }
  (void) ver;
  return 0;
ERR:
  mb_eh_clear(mb);
  return -1;
}

int mb_eh_readbin(mb_t *mb, const char *fname, int *pver)
{
  FILE *fp;
  int ver;
  int itmp;
  int i;
  int endn;

  if ( !(mb->eh_mode > 0) ) return 0;
  if ( !(mb->eh_mode == 1) ) {
    fprintf(stderr, "mb_eh_readbin: failed validation: mb->eh_mode == 1\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  if ((fp = fopen(fname, "rb")) == NULL) {
    fprintf(stderr, "cannot read binary file [%s].\n", fname);
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    return -1;
  }

  /* determine file endian */
  if ((endn = endn_rmatchi(&itmp, sizeof(int), fp)) < 0) {
    fprintf(stderr, "itmp 0x%X cannot match sizeof(int) 0x%X\n",
        (unsigned) itmp, (unsigned) sizeof(int));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  if (endn_fread(&itmp, sizeof(itmp), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading itmp\n");
    goto ERR;
  }
  if (itmp != (int) sizeof(double)) {
    fprintf(stderr, "(int) sizeof(double) mismatch, expect: %d, read: %d, pos: %#lx\n",
        (int) sizeof(double), itmp, (unsigned long) ftell(fp));
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    goto ERR;
  }
  if (endn_fread(&ver, sizeof(ver), 1, fp, endn) != 1) {
    fprintf(stderr, "error in reading ver\n");
    goto ERR;
  }
  if (pver != NULL) *pver = ver;

  /* call low level read function for members */
  i = mb_eh_readbin_low(mb, fp, ver, endn);
  fclose(fp);
  return i;
ERR:
  fclose(fp);
  return -1;
}

int mb_eh_writebin_low(mb_t *mb, FILE *fp, int ver)
{
  int i;
  int jmin;
  int jmax;
  double *pd;
  int size;

  if (mb == NULL) {
    fprintf(stderr, "passing null pointer to mb_eh_writebin_low\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    return -1;
  }
  /* n: number of temperature bins */
  if (endn_fwrite(&mb->n, sizeof(mb->n), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->n\n");
    goto ERR;
  }
  /* eh_cnt: number of energy bins */
  if (endn_fwrite(&mb->eh_cnt, sizeof(mb->eh_cnt), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->eh_cnt\n");
    goto ERR;
  }
  /* eh_min: minimal energy */
  if (endn_fwrite(&mb->eh_min, sizeof(mb->eh_min), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->eh_min\n");
    goto ERR;
  }
  /* eh_del: energy bin size */
  if (endn_fwrite(&mb->eh_del, sizeof(mb->eh_del), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing mb->eh_del\n");
    goto ERR;
  }
  /* eh_his: energy histogram data */
  for (i = 0; i < mb->n; i++) {
    pd = mb->eh_his + i * mb->eh_cnt;
    for (jmax = mb->eh_cnt; jmax > 0 && pd[jmax-1] <= 0.0; jmax--) ;
    for (jmin = 0; jmin < jmax && pd[jmin] <= 0.0; jmin++) ;
    if ((size = jmax - jmin) <= 0) continue;
    if (endn_fwrite(&i, sizeof(i), 1, fp, 1) != 1) {
      fprintf(stderr, "error in writing i\n");
      goto ERR;
    }
    if (endn_fwrite(&jmin, sizeof(jmin), 1, fp, 1) != 1) {
      fprintf(stderr, "error in writing jmin\n");
      goto ERR;
    }
    if (endn_fwrite(&size, sizeof(size), 1, fp, 1) != 1) {
      fprintf(stderr, "error in writing size\n");
      goto ERR;
    }
    if ((size > 0)
      && (endn_fwrite(pd+jmin, sizeof(*(pd+jmin)), size, fp, 1) != (size_t) (size))) {
      fprintf(stderr, "error in writing pd+jmin, n = size(%d)\n", size);
      goto ERR;
    }
  }
  (void) ver;
  return 0;
ERR:
  return -1;
}

int mb_eh_writebin(mb_t *mb, const char *fname, int ver)
{
  FILE *fp;
  int i;
  int size;

  if ( !(mb->eh_mode > 0) ) return 0;
  if ( !(mb->eh_mode == 1) ) {
    fprintf(stderr, "mb_eh_writebin: failed validation: mb->eh_mode == 1\n");
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    exit(1);
  }
  if ((fp = fopen(fname, "wb")) == NULL) {
    fprintf(stderr, "cannot write binary file [%s].\n", fname);
    fprintf(stderr, "FILE: %s, LINE: %d\n", __FILE__, __LINE__);
    return -1;
  }

  size = (int) sizeof(int);
  if (endn_fwrite(&size, sizeof(size), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing size\n");
    goto ERR;
  }
  size = (int) sizeof(double);
  if (endn_fwrite(&size, sizeof(size), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing size\n");
    goto ERR;
  }
  if (endn_fwrite(&ver, sizeof(ver), 1, fp, 1) != 1) {
    fprintf(stderr, "error in writing ver\n");
    goto ERR;
  }

  /* call low level write function for members */
  i = mb_eh_writebin_low(mb, fp, ver);
  fclose(fp);
  return i;
ERR:
  fclose(fp);
  return -1;
}

void mb_manifest(mb_t *mb, FILE *fp, int arrmax)
{
  int i;
  int pacnt;

  /* bdel: bin size of beta */
  fprintf(fp, "mb->bdel: double, %g\n", mb->bdel);

  /* n: number of temperature bins */
  fprintf(fp, "mb->n: int, %4d\n", mb->n);

  /* barr: temperature array */
  fprintf(fp, "mb->barr: dynamic array of mb->n+1: ");
  for (i = mb->n+1-1; i >= 0; i--) if (fabs(mb->barr[i]) > 1e-30) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n+1 > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n+1; i++) {
      if (i == arrmax && i < mb->n+1-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n+1-arrmax)) continue;
      fprintf(fp, "%g, ", mb->barr[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* beta: current value of beta */
  fprintf(fp, "mb->beta: double, %g\n", mb->beta);

  /* m: maximal number of bins in a window */
  fprintf(fp, "mb->m: int, %4d\n", mb->m);

  /* order: order, should be 1 */
  fprintf(fp, "mb->order: int, %4d\n", mb->order);

  /* flags: combination of flags */
  fprintf(fp, "mb->flags: unsigned, 0x%X\n", mb->flags);

  /* bwmod: 0: d(beta) 1: dT/T  2: d(kT) */
  fprintf(fp, "mb->bwmod: int, %4d\n", mb->bwmod);

  /* bwdel: delta lnT */
  fprintf(fp, "mb->bwdel: double, %g\n", mb->bwdel);

  /* MB_DAMP: use adaptive averaging */
  fprintf(fp, "mb->flags & MB_DAMP (0x%X, mbest_damp): %s\n",
    MB_DAMP, (mb->flags & MB_DAMP) ? "on" : "off");

  /* MB_CV: compute heat capacity */
  fprintf(fp, "mb->flags & MB_CV (0x%X, mbest_needcv): %s\n",
    MB_CV, (mb->flags & MB_CV) ? "on" : "off");

  /* MB_SYMWIN: use symmetrical window */
  fprintf(fp, "mb->flags & MB_SYMWIN (0x%X, mbest_sym_mbin): %s\n",
    MB_SYMWIN, (mb->flags & MB_SYMWIN) ? "on" : "off");

  /* MB_ONEBIN: use single bin estimator */
  fprintf(fp, "mb->flags & MB_ONEBIN (0x%X, mbest_single_bin): %s\n",
    MB_ONEBIN, (mb->flags & MB_ONEBIN) ? "on" : "off");

  /* MB_VERBOSE: being verbose */
  fprintf(fp, "mb->flags & MB_VERBOSE (0x%X, mbest_verbose): %s\n",
    MB_VERBOSE, (mb->flags & MB_VERBOSE) ? "on" : "off");

  /* MB_SBCORR: include energy fluctuation correction due to a small bin width for internal energy, etc */
  fprintf(fp, "mb->flags & MB_SBCORR (0x%X, mbest_sbcorr): %s\n",
    MB_SBCORR, (mb->flags & MB_SBCORR) ? "on" : "off");

  /* js: lower  boundary of asym. beta windows for ehat (asym.) */
  fprintf(fp, "mb->js: dynamic array of mb->n+1: ");
  for (i = mb->n+1-1; i >= 0; i--) if (mb->js[i]) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n+1 > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n+1; i++) {
      if (i == arrmax && i < mb->n+1-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n+1-arrmax)) continue;
      fprintf(fp, "%4d, ", mb->js[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* jt: higher boundary of asym. beta windows for ehat (asym.) */
  fprintf(fp, "mb->jt: dynamic array of mb->n+1: ");
  for (i = mb->n+1-1; i >= 0; i--) if (mb->jt[i]) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n+1 > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n+1; i++) {
      if (i == arrmax && i < mb->n+1-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n+1-arrmax)) continue;
      fprintf(fp, "%4d, ", mb->jt[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* jset: lower  boundary of beta windows for et (usu. sym.), i - jset[i] == jtet[i] - (i+1) */
  fprintf(fp, "mb->jset: dynamic array of mb->n: ");
  for (i = mb->n-1; i >= 0; i--) if (mb->jset[i]) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n; i++) {
      if (i == arrmax && i < mb->n-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n-arrmax)) continue;
      fprintf(fp, "%4d, ", mb->jset[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* jtet: higher boundary of beta windows for et (usu. sym.), i - jset[i] == jtet[i] - (i+1) */
  fprintf(fp, "mb->jtet: dynamic array of mb->n: ");
  for (i = mb->n-1; i >= 0; i--) if (mb->jtet[i]) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n; i++) {
      if (i == arrmax && i < mb->n-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n-arrmax)) continue;
      fprintf(fp, "%4d, ", mb->jtet[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* nstrefresh: interval of recalculating et for all temperature */
  fprintf(fp, "mb->nstrefresh: int, %4d\n", mb->nstrefresh);

  /* av_nstsave: interval of writing mbav and ze files */
  fprintf(fp, "mb->av_nstsave: int, %4d\n", mb->av_nstsave);

  /* av_binary: use binary format in mbav file */
  fprintf(fp, "mb->av_binary: int, %4d\n", mb->av_binary);

  /* av_file: name of mbav file */
  fprintf(fp, "mb->av_file: char *, %s\n", mb->av_file);

  /* ze_file: name of ze file */
  fprintf(fp, "mb->ze_file: char *, %s\n", mb->ze_file);

  /* wze_reps: number of iterations before writing ze file */
  fprintf(fp, "mb->wze_reps: int, %4d\n", mb->wze_reps);

  /* vis: number of visits */
  fprintf(fp, "mb->vis: dynamic array of mb->n: ");
  for (i = mb->n-1; i >= 0; i--) if (fabs(mb->vis[i]) > 1e-30) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n; i++) {
      if (i == arrmax && i < mb->n-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n-arrmax)) continue;
      fprintf(fp, "%g, ", mb->vis[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* totvis: total number of visits, number of tempering */
  fprintf(fp, "mb->totvis: double, %g\n", mb->totvis);

  /* winstot: total of sum.s over a multiple-bine temperature window */
  fprintf(fp, "mb->winstot: dynamic array of mb->n: ");
  for (i = mb->n-1; i >= 0; i--) if (fabs(mb->winstot[i]) > 1e-30) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n; i++) {
      if (i == arrmax && i < mb->n-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n-arrmax)) continue;
      fprintf(fp, "%g, ", mb->winstot[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* lgv_dt: time step for the temperature Langevin eq */
  fprintf(fp, "mb->lgv_dt: double, %g\n", mb->lgv_dt);

  /* lgv_dTmax: maximal amount of temperature change in a step */
  fprintf(fp, "mb->lgv_dTmax: double, %g\n", mb->lgv_dTmax);

  /* lgv_rej: number of attempts of langevin equation trying to change beta too drastically */
  fprintf(fp, "mb->lgv_rej: double, %g\n", mb->lgv_rej);

  /* lgv_tot: total number of attempts of using langevin equation */
  fprintf(fp, "mb->lgv_tot: double, %g\n", mb->lgv_tot);

  /* regl: average within a bin first */
  fprintf(fp, "mb->regl: int, %4d\n", mb->regl);

  /* fracmin: minimal allowable coefficient during left/right combination */
  fprintf(fp, "mb->fracmin: double, %g\n", mb->fracmin);

  /* cvshiftmax: maximal fraction for shift energy fluct. if cv is monotonic, it should be 0.0, for ising model, it can restrain the magnitude */
  fprintf(fp, "mb->cvshiftmax: double, %g\n", mb->cvshiftmax);

  /* ens_exp: ensemble exponent of beta */
  fprintf(fp, "mb->ens_exp: double, %g\n", mb->ens_exp);
  
  /* mode */
  fprintf(fp, "mb->mode: int, %d\n", mb->mode);
  if(mb->mode == 1)
  {
    /* beta0 */
    fprintf(fp, "mb->beta0: double, %g\n", mb->beta0);
    
    /* invsigma2 */
    fprintf(fp, "mb->invsigma2: double, %g\n", mb->invsigma2);
  }
  else if(mb->mode == 2)
  {
    /* c */
    fprintf(fp, "mb->c: double, %g\n", mb->c);
  }

  /* ens_w: array of ensemble weights at bin boundaries */
  fprintf(fp, "mb->ens_w: dynamic array of mb->n+1: ");
  for (i = mb->n+1-1; i >= 0; i--) if (fabs(mb->ens_w[i]) > 1e-30) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n+1 > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n+1; i++) {
      if (i == arrmax && i < mb->n+1-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n+1-arrmax)) continue;
      fprintf(fp, "%g, ", mb->ens_w[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* shk_base: current generic shrink amplitude */
  fprintf(fp, "mb->shk_base: double, %g\n", mb->shk_base);

  /* shk_winadj: adjust shrink according to temperature window width */
  fprintf(fp, "mb->shk_winadj: int, %4d\n", mb->shk_winadj);

  /* shk_max: initial and maximal shrink (adjusted) */
  fprintf(fp, "mb->shk_max: double, %g\n", mb->shk_max);

  /* shk_gauge: array used of modulation shrinking factors */
  fprintf(fp, "mb->shk_gauge: dynamic array of mb->n: ");
  for (i = mb->n-1; i >= 0; i--) if (fabs(mb->shk_gauge[i]) > 1e-30) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n; i++) {
      if (i == arrmax && i < mb->n-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n-arrmax)) continue;
      fprintf(fp, "%g, ", mb->shk_gauge[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* shk_mode: 0: const, 1: amp/t, 2: amp/t^exp */
  fprintf(fp, "mb->shk_mode: int, %4d\n", mb->shk_mode);

  /* shk_min: minimal value for enforcing acc. sampling */
  fprintf(fp, "mb->shk_min: double, %g\n", mb->shk_min);

  /* shk_stop: stop shrinking after this number of steps */
  fprintf(fp, "mb->shk_stop: int, %4d\n", mb->shk_stop);

  if (mb->shk_mode >= 1) {
    /* shk_amp: amp t^(-exp) */
    fprintf(fp, "mb->shk_amp: double, %g\n", mb->shk_amp);
  }
  if (mb->shk_mode >= 2) {
    /* shk_exp: amp t^(-exp) */
    fprintf(fp, "mb->shk_exp: double, %g\n", mb->shk_exp);
  }
  /* lnz: logarithm of the partition function */
  fprintf(fp, "mb->lnz: dynamic array of mb->n+1: ");
  for (i = mb->n+1-1; i >= 0; i--) if (fabs(mb->lnz[i]) > 1e-30) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n+1 > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n+1; i++) {
      if (i == arrmax && i < mb->n+1-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n+1-arrmax)) continue;
      fprintf(fp, "%g, ", mb->lnz[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* ehat: internal energy */
  fprintf(fp, "mb->ehat: dynamic array of mb->n+1: ");
  for (i = mb->n+1-1; i >= 0; i--) if (fabs(mb->ehat[i]) > 1e-30) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n+1 > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n+1; i++) {
      if (i == arrmax && i < mb->n+1-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n+1-arrmax)) continue;
      fprintf(fp, "%g, ", mb->ehat[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* cvhat: heat capacity */
  fprintf(fp, "mb->cvhat: dynamic array of mb->n+1: ");
  for (i = mb->n+1-1; i >= 0; i--) if (fabs(mb->cvhat[i]) > 1e-30) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n+1 > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n+1; i++) {
      if (i == arrmax && i < mb->n+1-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n+1-arrmax)) continue;
      fprintf(fp, "%g, ", mb->cvhat[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* et: bin-averaged internal energy */
  fprintf(fp, "mb->et: dynamic array of mb->n: ");
  for (i = mb->n-1; i >= 0; i--) if (fabs(mb->et[i]) > 1e-30) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n; i++) {
      if (i == arrmax && i < mb->n-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n-arrmax)) continue;
      fprintf(fp, "%g, ", mb->et[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* imbal: |a+ - a-| / (a+ + a-) for left-right combination */
  fprintf(fp, "mb->imbal: dynamic array of mb->n+1: ");
  for (i = mb->n+1-1; i >= 0; i--) if (fabs(mb->imbal[i]) > 1e-30) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n+1 > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n+1; i++) {
      if (i == arrmax && i < mb->n+1-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n+1-arrmax)) continue;
      fprintf(fp, "%g, ", mb->imbal[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* haset: current et[i] is reasonably good */
  fprintf(fp, "mb->haset: dynamic array of mb->n: ");
  for (i = mb->n-1; i >= 0; i--) if (mb->haset[i]) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n; i++) {
      if (i == arrmax && i < mb->n-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n-arrmax)) continue;
      fprintf(fp, "%4d, ", mb->haset[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* qua: bits represent whether estimated values are unbiased */
  fprintf(fp, "mb->qua: dynamic array of mb->n+1: ");
  for (i = mb->n+1-1; i >= 0; i--) if (mb->qua[i]) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n+1 > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n+1; i++) {
      if (i == arrmax && i < mb->n+1-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n+1-arrmax)) continue;
      fprintf(fp, "0x%X, ", mb->qua[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* ampf: currently amplification factor for adaptive averaging */
  fprintf(fp, "mb->ampf: dynamic array of mb->n: ");
  for (i = mb->n-1; i >= 0; i--) if (fabs(mb->ampf[i]) > 1e-30) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n; i++) {
      if (i == arrmax && i < mb->n-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n-arrmax)) continue;
      fprintf(fp, "%g, ", mb->ampf[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* sums: normal data */
  if (mb->sums != NULL) {
    fprintf(fp, "mb->sums: sm_t array of mb->n:");
    for (i = mb->n*sizeof(sm_t)-1; i >= 0; i--) if (*((char *)mb->sums + i)) break;
    if (i >= 0) {
      fprintf(fp, "\n");
      for (i = 0; i < mb->n; i++) {
        if (i == arrmax && i < mb->n-arrmax) {
          fprintf(fp, "\n...\n");
        }
        if (arrmax >= 0 && i >= arrmax && i < (mb->n-arrmax)) continue;
        fprintf(fp, "mb->sums[%d]:\n", i);
        sm_manifest(mb->sums+i, fp, arrmax);
      }
    } else {
      fprintf(fp, " {0}\n");
    }
  }

  /* has_xsums */
  fprintf(fp, "mb->has_xsums: int, %4d\n", mb->has_xsums);

  /* cnt_int: number of additional integer variables to be written to binary file */
  fprintf(fp, "mb->cnt_int: int, %4d\n", mb->cnt_int);

  /* cnt_dbl: number of additional double variables to be written to binary file */
  fprintf(fp, "mb->cnt_dbl: int, %4d\n", mb->cnt_dbl);

  /* idxcnt: index count */
  fprintf(fp, "mb->idxcnt: dynamic array of mb->n: ");
  for (i = mb->n-1; i >= 0; i--) if (mb->idxcnt[i]) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n; i++) {
      if (i == arrmax && i < mb->n-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n-arrmax)) continue;
      fprintf(fp, "%4d, ", mb->idxcnt[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* midx: index look-up table */
  fprintf(fp, "mb->midx: dynamic array of mb->n*mb->m: ");
  for (i = mb->n*mb->m-1; i >= 0; i--) if (mb->midx[i]) break;
  if (i >= 0) {
    if ((arrmax < 0 || arrmax > 3) && mb->n*mb->m > 6)
      fprintf(fp, "\n");
    for (pacnt = 0, i = 0; i < mb->n*mb->m; i++) {
      if (i == arrmax && i < mb->n*mb->m-arrmax) {
        if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
        fprintf(fp, "..., ");
        if (arrmax > 3) fprintf(fp, "\n");
      }
      if (arrmax >= 0 && i >= arrmax && i < (mb->n*mb->m-arrmax)) continue;
      fprintf(fp, "%4d, ", mb->midx[i]);
      if (++pacnt % 10 == 0) fprintf(fp, "\n");
    }
    if (pacnt % 10 != 0) fprintf(fp, "\n");
  } else {
    fprintf(fp, " {0}\n");
  }

  /* xsums: multiple-bin damping data */
  if (mb->xsums != NULL) {
    fprintf(fp, "mb->xsums: sm_t array of mb->n*mb->m:");
    for (i = mb->n*mb->m*sizeof(sm_t)-1; i >= 0; i--) if (*((char *)mb->xsums + i)) break;
    if (i >= 0) {
      fprintf(fp, "\n");
      for (i = 0; i < mb->n*mb->m; i++) {
        if (i == arrmax && i < mb->n*mb->m-arrmax) {
          fprintf(fp, "\n...\n");
        }
        if (arrmax >= 0 && i >= arrmax && i < (mb->n*mb->m-arrmax)) continue;
        fprintf(fp, "mb->xsums[%d]:\n", i);
        sm_manifest(mb->xsums+i, fp, arrmax);
      }
    } else {
      fprintf(fp, " {0}\n");
    }
  }

  /* eh_mode: 0: disable; 1: simple histogram */
  fprintf(fp, "mb->eh_mode: int, %4d\n", mb->eh_mode);

  if (mb->eh_mode) {
    /* eh_skip: interval of reconstructing energy histograms */
    fprintf(fp, "mb->eh_skip: int, %4d\n", mb->eh_skip);

    /* eh_bwmod: 0: d(beta) 1: dT/T  2: d(kT) */
    fprintf(fp, "mb->eh_bwmod: int, %4d\n", mb->eh_bwmod);

    /* eh_bwdel: delta lnT */
    fprintf(fp, "mb->eh_bwdel: double, %g\n", mb->eh_bwdel);

    /* eh_min: minimal energy */
    fprintf(fp, "mb->eh_min: double, %g\n", mb->eh_min);

    /* eh_max: maximal energy */
    fprintf(fp, "mb->eh_max: double, %g\n", mb->eh_max);

    /* eh_del: energy bin size */
    fprintf(fp, "mb->eh_del: double, %g\n", mb->eh_del);

    /* eh_cnt: number of energy bins */
    fprintf(fp, "mb->eh_cnt: int, %4d\n", mb->eh_cnt);

    /* eh_binary: binary format for ehist file */
    fprintf(fp, "mb->eh_binary: int, %4d\n", mb->eh_binary);

    /* eh_nstsave: interval of writing histogrm files */
    fprintf(fp, "mb->eh_nstsave: int, %4d\n", mb->eh_nstsave);

    /* eh_file: name of ehist file */
    fprintf(fp, "mb->eh_file: char *, %s\n", mb->eh_file);

    /* eh_rfile: name of reconstructed energy histogram */
    fprintf(fp, "mb->eh_rfile: char *, %s\n", mb->eh_rfile);

    /* eh_his: energy histogram data */
    fprintf(fp, "mb->eh_his: dynamic array of mb->n*mb->eh_cnt: ");
    for (i = mb->n*mb->eh_cnt-1; i >= 0; i--) if (fabs(mb->eh_his[i]) > 1e-30) break;
    if (i >= 0) {
      if ((arrmax < 0 || arrmax > 3) && mb->n*mb->eh_cnt > 6)
        fprintf(fp, "\n");
      for (pacnt = 0, i = 0; i < mb->n*mb->eh_cnt; i++) {
        if (i == arrmax && i < mb->n*mb->eh_cnt-arrmax) {
          if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
          fprintf(fp, "..., ");
          if (arrmax > 3) fprintf(fp, "\n");
        }
        if (arrmax >= 0 && i >= arrmax && i < (mb->n*mb->eh_cnt-arrmax)) continue;
        fprintf(fp, "%g, ", mb->eh_his[i]);
        if (++pacnt % 10 == 0) fprintf(fp, "\n");
      }
      if (pacnt % 10 != 0) fprintf(fp, "\n");
    } else {
      fprintf(fp, " {0}\n");
    }

    /* eh_recon: temporary space for reconstructing histogram */
    fprintf(fp, "mb->eh_recon: dynamic array of mb->eh_cnt: ");
    for (i = mb->eh_cnt-1; i >= 0; i--) if (fabs(mb->eh_recon[i]) > 1e-30) break;
    if (i >= 0) {
      if ((arrmax < 0 || arrmax > 3) && mb->eh_cnt > 6)
        fprintf(fp, "\n");
      for (pacnt = 0, i = 0; i < mb->eh_cnt; i++) {
        if (i == arrmax && i < mb->eh_cnt-arrmax) {
          if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
          fprintf(fp, "..., ");
          if (arrmax > 3) fprintf(fp, "\n");
        }
        if (arrmax >= 0 && i >= arrmax && i < (mb->eh_cnt-arrmax)) continue;
        fprintf(fp, "%g, ", mb->eh_recon[i]);
        if (++pacnt % 10 == 0) fprintf(fp, "\n");
      }
      if (pacnt % 10 != 0) fprintf(fp, "\n");
    } else {
      fprintf(fp, " {0}\n");
    }

    /* eh_is: indices for temperature windows (lower) */
    fprintf(fp, "mb->eh_is: dynamic array of mb->n + 1: ");
    for (i = mb->n + 1-1; i >= 0; i--) if (mb->eh_is[i]) break;
    if (i >= 0) {
      if ((arrmax < 0 || arrmax > 3) && mb->n + 1 > 6)
        fprintf(fp, "\n");
      for (pacnt = 0, i = 0; i < mb->n + 1; i++) {
        if (i == arrmax && i < mb->n + 1-arrmax) {
          if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
          fprintf(fp, "..., ");
          if (arrmax > 3) fprintf(fp, "\n");
        }
        if (arrmax >= 0 && i >= arrmax && i < (mb->n + 1-arrmax)) continue;
        fprintf(fp, "%4d, ", mb->eh_is[i]);
        if (++pacnt % 10 == 0) fprintf(fp, "\n");
      }
      if (pacnt % 10 != 0) fprintf(fp, "\n");
    } else {
      fprintf(fp, " {0}\n");
    }

    /* eh_it: indices for temperature windows (higher) */
    fprintf(fp, "mb->eh_it: dynamic array of mb->n + 1: ");
    for (i = mb->n + 1-1; i >= 0; i--) if (mb->eh_it[i]) break;
    if (i >= 0) {
      if ((arrmax < 0 || arrmax > 3) && mb->n + 1 > 6)
        fprintf(fp, "\n");
      for (pacnt = 0, i = 0; i < mb->n + 1; i++) {
        if (i == arrmax && i < mb->n + 1-arrmax) {
          if (arrmax > 3 && pacnt % 10 != 0) fprintf(fp, "\n");
          fprintf(fp, "..., ");
          if (arrmax > 3) fprintf(fp, "\n");
        }
        if (arrmax >= 0 && i >= arrmax && i < (mb->n + 1-arrmax)) continue;
        fprintf(fp, "%4d, ", mb->eh_it[i]);
        if (++pacnt % 10 == 0) fprintf(fp, "\n");
      }
      if (pacnt % 10 != 0) fprintf(fp, "\n");
    } else {
      fprintf(fp, " {0}\n");
    }

    /* MB_EH_ADDAHALF: add a half energy bin width in output */
    fprintf(fp, "mb->flags & MB_EH_ADDAHALF (0x%X, ehist_addahalf): %s\n",
      MB_EH_ADDAHALF, (mb->flags & MB_EH_ADDAHALF) ? "on" : "off");

    /* MB_EH_KEEPEDGE: keep zero edge at sides */
    fprintf(fp, "mb->flags & MB_EH_KEEPEDGE (0x%X, ehist_keepedge): %s\n",
      MB_EH_KEEPEDGE, (mb->flags & MB_EH_KEEPEDGE) ? "on" : "off");

    /* MB_EH_NOZEROES: do not output zeroes */
    fprintf(fp, "mb->flags & MB_EH_NOZEROES (0x%X, ehist_nozeroes): %s\n",
      MB_EH_NOZEROES, (mb->flags & MB_EH_NOZEROES) ? "on" : "off");
  }
}

#undef ZCOM_PICK
#undef ZCOM_ENDN
