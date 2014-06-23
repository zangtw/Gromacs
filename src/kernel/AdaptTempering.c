#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <unistd.h>

#include "AdaptTemperingMultiBin.h"
#include "AdaptTempering.h"

#define ZCOM_PICK
#define ZCOM_RNG
#include "zcom1.h"

/*################# TYPE DEFINATIONS##################*/

/*################################*/
/* at operations, static part*/

#define Beta2T(x) (1.0/(BOLTZ * x))
#define T2Beta(x) (1.0/(BOLTZ * x))

/* load previous data */
static int at_loaddata(at_t *at, bool bCPT)
{
  mb_t *mb = at->mb;
  int ver = 0, readmb = 1;

  if (!bCPT) /* initial run */
    return 0;
	/*readmb = (at->premode > 0) ? 0 : bCPT;*/
	readmb = bCPT;
  if (readmb) {
    if (mb_readbin(mb, mb->av_file, &ver) != 0) { /* read stat. */
      fprintf(stderr, "cannot load mb data from %s\n", mb->av_file);
      return 1;
    }
    fprintf(stderr, "%s version: %d\n", mb->av_file, ver);
    mb_wze(mb, "ze_init");
		at->beta = at->mb->beta; /* update the current temperature */

    if (mb_eh_readbin(mb, mb->eh_file, &ver) != 0) { /* read E-histograms */
      fprintf(stderr, "cannot load energy histogram from %s\n", mb->eh_file);
      return 2;
    }
  }

  fprintf(stderr,"successfully load previous data\n");
  return 0;
}

/* don't write at the first step; write at the last step */
static int at_doevery(llong_t step, int nst, bool bfirst, bool blast)
{
  return !bfirst && (blast || (nst > 0 && step % nst == 0));
}

/* write various output files */
static void at_output(at_t *at, llong_t step,
		      int ib, double invw, double t1, double t2, double Eav,
		      bool bfirst, bool blast, bool btr, bool bflush)
{
  bool btrace;

  /* write the trace file */
  if (at->nsttrace > 0)
    btrace = (step % at->nsttrace == 0) || bfirst || blast;
  else /* tracing is disabled if at->nsttrace == 0 */
    btrace = (at->nsttrace == 0) ? btr : 0;

	if(bflush)
		at->log->flags |= LOG_FLUSHAFTER;
	else
		at->log->flags ^= LOG_FLUSHAFTER;

  if (btrace) {
    log_printf(at->log, "%10.3f %5d %10.6f %12.3f %12.3f %10.6f %8.3f %8.5f",
      step * at->tmstep, ib, t2 - t1, at->Ea, Eav, at->beta, t1, invw);    
    log_printf(at->log, "\n");
  }

  if (at_doevery(step, at->mb->av_nstsave, bfirst, blast)) { /* save averages */
    mb_write(at->mb);
    mb_wze(at->mb, NULL);
    mtsave(at->rng_file);
  }
  if (at_doevery(step, at->mb->eh_nstsave, bfirst, blast)) { /* save energy histograms */
    mb_eh_writebin(at->mb, at->mb->eh_file, 0);
    mb_eh_recon(at->mb, NULL);
  }
}

/* The only reason we copy the whole beta array from 'mb' to 'at' is to reduce the MPI broadcast time when doing at_move. */
static double *at_barr_init(mb_t *mb)
{
  double *barr;
  int n = mb->n;

  if ((barr = calloc(n, sizeof(double)) )== NULL) {
    fprintf(stderr, "no memory for barr (size %u)\n", (unsigned) n*sizeof(double));
    exit(1);
  }

  int i;

  for(i=0; i<n; i++)
		barr[i] = mb->barr[i];

  return barr;
}

static double *at_kappa_init(mb_t *mb, double Tref, double kappa0)
{
  double *kappa;
  int n = mb->n;
  double Tmin = Beta2T(mb->bmax);
  double Tmax = Beta2T(mb->bmin);
  double T;
  double kappa1 = 1 - kappa0;

  if ((kappa = calloc(n, sizeof(double)) )== NULL) {
    fprintf(stderr, "no memory for kappa (size %u)\n", (unsigned) n*sizeof(double));
    exit(1);
  }

  int i;

  for(i=0; i<n; i++){

    T = Beta2T(mb->barr[i]);
    
    if (T > Tref)
      kappa[i] = 1 - (T - Tref) * kappa1 / (Tmax - Tref);
    else
      kappa[i] = 1;
  }

  return kappa;
}

static double *at_epsilon_init(mb_t *mb, double Tref, double epsilon0)
{
  double *epsilon;
  int n = mb->n;
  double Tmin = Beta2T(mb->bmax);
  double Tmax = Beta2T(mb->bmin);
  double T;

  if ((epsilon = calloc(n, sizeof(double)) ) == NULL) {
    fprintf(stderr, "no memory for epsilon (size %u)\n", (unsigned) n*sizeof(double));
    exit(1);
  }

  int i;

  for(i=0; i<n; i++){

    T = Beta2T(mb->barr[i]);
    
    if (T > Tref)
      epsilon[i] = epsilon0 * (T - Tref) / (Tmax - Tref);
    else
      epsilon[i] = 0;
  }

  return epsilon;
}

static void at_manifest(at_t *at, FILE *fp, int arrmax)
{
  fprintf(fp, "at->bmin: double, %g\n", at->bmin);
  fprintf(fp, "at->bmax: double, %g\n", at->bmax);
  fprintf(fp, "at->T0: double, %g\n", at->T0);
  fprintf(fp, "at->beta: double, %g\n", at->beta);
  fprintf(fp, "at->nsttemp: int, %4d\n", at->nsttemp);
  fprintf(fp, "at->mvreps: int, %4d\n", at->mvreps);
  fprintf(fp, "at->tmstep: double, %g\n", at->tmstep);
  fprintf(fp, "at->rng_file: char *, %s\n", at->rng_file);
  fprintf(fp, "at->nsttrace: int, %4d\n", at->nsttrace);
  fprintf(fp, "at->grand: function pointer, %p\n", at->grand);
  fprintf(fp, "at->trace_file: char *, %s\n", at->trace_file);
  fprintf(fp, "at->log->fname: char *, %s\n", at->log->fname);
  fprintf(fp, "at->bTH: bool, %s",(at->bTH)?"true":"false");
  if (at->bTH) {
    fprintf(fp, "at->th_Tref: double, %g\n", at->TH_Tref);
    fprintf(fp, "at->kappa0: double, %g\n", at->kappa0);
    fprintf(fp, "at->epsilon0: double, %g\n", at->epsilon0);
  }
  if (at->mb != NULL) {
    fprintf(fp, "at->mb: object pointer to mb_t\n");
    mb_manifest(at->mb, fp, arrmax);
  }
  fprintf(fp, "at->Ea: double, %g\n", at->Ea);
}

/* at_cfgopen_low: initialize members of at_t from configuration
* file `cfg', or if unavailable, from default values */
static int at_cfgopen_low(at_t *at, cfgdata_t *cfg, double tmstep)
{
	char buf1[100], buf2[100];
	char suf[2];
	suf[0] = at->suffix;
	suf[1] = '\0';
	/* NOTE: only at->beta won't be initialized in this function */

  die_if((at == NULL),"null pointer to at_t\n");
  die_if((cfg == NULL),"cannot load *.cfg file\n");
  
	/* bmin: minimal beta (highest temperature) */
  die_if(cfgget(cfg, &at->bmin, "beta_min", "%lf"), "missing var: at->bmin, key: beta_min, fmt: %%lf\n");
  die_if((at->bmin <= 0.0), "at->bmin should be positive!\n");
	
	/* bmax: maximum beta (lowest temperature) */
  die_if(cfgget(cfg, &at->bmax, "beta_max", "%lf"), "missing var: at->bmax, key: beta_min, fmt: %%lf\n");
  die_if((at->bmax <= 0.0), "at->bmax: should be positive!\n");
  
	die_if(!(at->bmax > at->bmin), "at->bmax: failed validation: at->bmax > at->bmin\n");

  /* T0: thermostat temperature */
  at->T0 = 300.0;
  if (cfgget(cfg, &at->T0, "T0", "%lf"))
    fprintf(stderr, "assuming default: at->T0 = 300.0, key: T0\n");

  /* nsttemp: frequency of tempering, 0: disable, -1: only ns */
  at->nsttemp = -1;
  if (cfgget(cfg, &at->nsttemp, "nsttemp", "%d"))
		fprintf(stderr, "assuming default: at->nsttemp = -1, key: nsttemp\n");

  /* mvreps: number of repeating Langevin eq */
  at->mvreps = 1;
  if (cfgget(cfg, &at->mvreps, "move_repeats", "%d"))
		fprintf(stderr, "assuming default: at->mvreps = 1, key: move_repeats\n");

  /* tmstep: MD integration step, for convenience */
  at->tmstep = tmstep;

  /* nsttrace: interval of writing trace file; -1: only when doing neighbor search, 0: disable */
  at->nsttrace = -1;
  if (cfgget(cfg, &at->nsttrace, "nsttrace", "%d"))
		fprintf(stderr, "assuming default: at->nsttrace = -1, key: nsttrace\n");

  /* grand: function pointer to a gaussian random number generator */
  at->grand = &grand0;
  
	/* rng_file: file name of random number state */
  if (cfgget(cfg, &buf1, "rng_file", "%s"))
	{
		fprintf(stderr, "assuming default: at->rng_file = \"MTSEED\", key: rng_file\n");
		strcpy(buf1, "MTSEED");
	}
	strcat(buf1, suf);
  at->rng_file = ssdup(buf1);

  /* trace_file: name of trace file */
  if (cfgget(cfg, &buf2, "trace_file", "%s"))
	{
		fprintf(stderr, "assuming default: at->trace_file = \"TRACE\", key: trace_file\n");
		strcpy(buf2, "TRACE");
	}
	strcat(buf2, suf);
  at->trace_file = ssdup(buf2);

  /* log: logfile */
  at->log = NULL;

  /* bTH : 0: disable; 1:enable */
  at->bTH = 0;
  if (cfgget(cfg, &at->bTH, "boost_mode", "%d")) 
		fprintf(stderr, "assuming default: at->th_mode = 0, key: boost_mode\n");

  /* TH_Tref */
  at->TH_Tref = 300.0;
  if (at->bTH)
    if (cfgget(cfg, &at->TH_Tref, "boost_Tref", "%lf"))
			fprintf(stderr, "assuming default: at->th_Tref = 300.0, key: boost_Tref\n");

  /* kappa0 */
  at->kappa0 = 1.0;
  if (at->bTH) 
    if (cfgget(cfg, &at->kappa0, "kappa0", "%lf"))
      fprintf(stderr, "assuming default: at->kappa0 = 1.0, key: kappa0\n");

  /* epsilon0 */
  at->epsilon0 = 0.0;
  if (at->bTH)
    if (cfgget(cfg, &at->epsilon0, "epsilon0", "%lf"))
			fprintf(stderr, "assuming default: at->epsilon0 = 0.0, key: epsilon0\n");
  
  /* mb: handle for multiple-bin estimator */
	at->mb = mb_cfgopen(cfg, at->bmin, at->bmax, at->suffix);
  die_if((at->mb == NULL), "failed to initialize at->mb\n\n");
	
  /* Ea: total potential energy */
  at->Ea = 0.0;

	return 1;
}

/* return a pointer of an initialized at_t
 * if possible, initial values are taken from configuration
 * file `cfg', otherwise default values are assumed */
static at_t *AdaptTempering_OpenCfg(const char *cfgname, double tmstep, int suffix)
{
  cfgdata_t *cfg;
  at_t *at;
	bool bLoaded;
	char *p;
	int delay;

  /* open configuration file */
  die_if(!(cfg = cfgopen(cfgname)), "at_t: cannot open config. file %s.\n", cfgname);

  /* allocate memory for at_t */
	xnew(at, 1);

	/* Get the file suffix first */
	die_if(suffix >= 10, "do not support # of simulations > 10 currently\n");
	at->suffix = (char)(((int)'0') + suffix);
  
	/* call low level function */
  die_if (!(bLoaded = at_cfgopen_low(at,cfg,tmstep)), "at_t: error while reading configuration file %s\n", cfgname);
	
	printf("Successfully loaded cfg data!\n");
	
	/* generate different random seeds in multi-simulation */
	delay = suffix * 2;
	sleep(delay);

	/* load random number generator */
	mtload(at->rng_file, 0);
	
  /* close handle to configuration file */
  cfgclose(cfg);
  return at;
}

/* at_close_low: close a pointer to at_t */
static void at_close_low(at_t *at)
{
  if (at->mpi_rank == 0) mtsave(at->rng_file);
  if (at->rng_file   != NULL) ssdelete(at->rng_file);
  if (at->mpi_rank == 0) log_close(at->log);
  if (at->trace_file != NULL) ssdelete(at->trace_file);
  if (at->mb         != NULL) mb_close(at->mb);
  memset(at, 0, sizeof(*at));
}

/*########################END OF STATIC FUNCTIONS ########################*/

real AdaptTempering_CurrentBeta(at_t *at)
{
	return at->beta;
}

real AdaptTempering_CurrentT(at_t *at)
{
	return Beta2T(at->beta);
}

real AdaptTempering_CurrentPara(at_t *at)
{
	if (at->mb->mode == 0)
		return 1;
	else if (at->mb->mode == 1)
		return at->mb->beta0 * at->mb->invsigma2;
	else if (at->mb->mode == 2)
		return at->mb->c;
}

real AdaptTempering_CurrentSecondPara(at_t *at)
{
	if (at->mb->mode == 0)
		return 1;
	else if (at->mb->mode == 1)
		return -0.5 * at->mb->invsigma2;
	else if (at->mb->mode == 2)
		return 1;
}

real AdaptTempering_ReferenceTemperature(at_t *at)
{
  return at->T0;
}

void AdaptTempering_ForceChangeBeta(at_t *at, double newbeta)
{
	at->beta = at->mb->beta = newbeta;

	AdaptTempering_UpdateTemperature(at);
}

at_t *AdaptTempering_MasterCreate(const char *fname, bool bCPT, double tmstep, int suffix)
{
  at_t *at;

  /* this will also initialize settings for member objects such as at->mb */
  at = AdaptTempering_OpenCfg((fname != NULL) ? fname : "at.cfg", tmstep, suffix); 
	die_if(at == NULL, "failed to load configuration file.\n"); 
	
	/* Make the initial temperature = T0 */
	AdaptTempering_ForceChangeBeta(at, Beta2T(at->T0));
  
	/* we only load previous data if it's continuation */
  if (at_loaddata(at, bCPT) != 0) {
		fprintf(stderr, "Warning: This simulation is started from checkpoint, while some files are missing. Will assume no previous simulation data is available.\n");

		/* AdaptTempering_Close(at);
		return NULL; */
  }

	return at;
}

at_t *AdaptTempering_NonMasterCreate()
{
	at_t *at;

	xnew(at, 1);

	return at;
}

void AdaptTempering_OpenLog(at_t *at)
{
	die_if(at == NULL, "failed to load adaptive tempering data.\n");
	die_if(at->trace_file == NULL, "failed to get the name of trace file. Is this function called by a non-master node?\n");
	
	at->log = log_open(at->trace_file); /* set an attempt value for at->beta before reading the value from mb.av */
  
  AdaptTempering_UpdateTemperature(at); /* update temperature */
}

int AdaptTempering_DumpToFile(at_t *at, const char *fname, int arrmax)
{
  FILE *fp;
	char buf[100];
	char suf[2];
	suf[0] = at->suffix;
	suf[1] = '\0';

	strcpy(buf, fname);
	if (at->suffix)
		strcat(buf, suf);
	
	if ((fp = fopen(buf, "w")) == NULL)
	{
			fprintf(stderr, "cannot write %s\n", fname);
			return -1;
  }
  at_manifest(at, fp, arrmax);
  fclose(fp);
  return 0;
}

int AdaptTempering_Langevin(at_t *at, llong_t step, bool bfirst, bool blast, bool btr, bool bflush)
{
  double invwf = 1.0, T1 = 0., T2 = 0., Eav = 0., ndlnwfdbeta;
  int ib, rep;
  double *varr = NULL;
	double noise;
  
	noise = (*at->grand)();
	
	die_if (at->grand == NULL, "no gaussian RNG\n");

  /* update energy data, change at->beta */
  /* repeat several times to change the temperature */
  for (rep = 0; rep < at->mvreps; rep++) {
    /* 1. deposit the current energy and temperature */
    mb_add(at->mb, at->Ea, varr, at->beta, &ib, &invwf, &ndlnwfdbeta);

    /* 2. use Langevin equation to update the temperature */
    T1 = Beta2T(at->beta);
    at->beta = mb_move(at->mb, at->Ea, at->beta, ib, ndlnwfdbeta, noise, &Eav);
    T2 = Beta2T(at->beta);
  }

  if (at_doevery(step, at->mb->nstrefresh, 0, blast))
    mb_refresh_et(at->mb, 1);
  at_output(at, step, ib, invwf, T1, T2, Eav, bfirst, blast, btr, bflush);
  return 0;
}

void AdaptTempering_UpdateTemperature(at_t *at)
{
	at->beta = at->mb->beta;
}

real AdaptTempering_ForceScaleFactor(at_t *at)
{
	return (at->T0 / Beta2T(at->beta));
}

void AdaptTempering_Close(at_t *at)
{ at_close_low(at); free(at); at=NULL; }

#ifdef GMX_MPI
int AdaptTempering_SyncAllNodes(at_t *at, MPI_Comm comm)
{
  int mpisize = 1;
  int mpirank = 0;

  die_if (at == NULL, "null pointer at to at_t\n");
  
	if (comm != MPI_COMM_NULL) 
    die_if((MPI_SUCCESS != MPI_Comm_size(comm, &mpisize)), "cannot even get MPI size\n");
	
  if (mpisize > 1) {
    die_if((MPI_SUCCESS != MPI_Comm_rank(comm, &mpirank)), "cannot get MPI rank\n");

		/* Broadcast at_t */
		die_if((MPI_SUCCESS != MPI_Bcast(at, sizeof(*at), MPI_BYTE, 0, comm)), "%3d/%3d: failed to bcast at (%p), type = *at, size = 1 (%d), comm = 0x%lX\n", mpirank, mpisize, at, 1, (unsigned long) comm);
	}

  at->mpi_comm = comm;
  at->mpi_size = mpisize;
  at->mpi_rank = mpirank;

  if (at->mpi_rank != 0) {
    at->rng_file = NULL;
    at->trace_file = NULL;
    at->mb = NULL;
  }
  return 0;
}
#endif

#undef ZCOM_PICK
#undef ZCOM_RNG
