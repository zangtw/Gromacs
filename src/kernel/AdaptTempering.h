/* Accelerated tempering types */
#ifndef ADAPTIVE_TEMPERING_H__
#define ADAPTIVE_TEMPERING_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "AdaptTemperingMultiBin.h"

#define ZCOM_PICK
#define ZCOM_LOG
#include "zcom1.h"

#ifdef GMX_MPI
#include "mpi.h"
#endif

#ifndef BOLTZ
#define BOLTZ  8.314511212e-3 /* Boltzmann constant */
#endif

#ifndef HAVE_BOOL
#define bool int
#endif

/* define a long long int type */
#ifndef I32ONLY
typedef long long  llong_t;
#define llong_pfmt "%lld"
#else
typedef long llong_t;
#define llong_pfmt "%ld"
#endif

typedef struct AdaptTemperingType at_t;

struct AdaptTemperingType{
  /*------Temperature bins-------------------------*/
  double    bmin;     // minimal beta 
  double    bmax;     // maximal beta
	double    beta;     // current beta 
  double    T0;       // thermostat temperature 
	int			  index;    // index of current beta

  /*-----Integration of MD and Langevin related ---*/
  int       mvreps;    // number of repeating Langevin eq 
  double    tmstep;    // MD integration step
	int       nsttemp;   // interval of tempering, 0: disable, -1: only when doing neighbor searching
  int       nsttrace;  // interval of writing trace file, 0: disable, -1: only when doing neighbor searching
	double    (*grand)(void);  /* function pointer to a gaussian random number genrator */

  /*------Files------------------------------------*/
  char      *rng_file;    // random number file name
  char      *trace_file;  // trace file name
  logfile_t *log;         // log file
	char      suffix;       // file suffix

  /*------Combined Hamiltonian:
    H = kappa* H0 + epsilon * H1
    kappa = 1-(T-Tref)*(1-kappa0)/(Tmax-Tref) if T>Tref; kappa=1 if T<Tref
    epsilon= epsilon0*(T-Tref)/(Tmax-Tref) if T>Tref; epsilon=0 if T<Tref */
  bool      bTH;	
  double    TH_Tref;
  double    *kappa, *epsilon;
  double    kappa0, epsilon0;

  /*-----Multiple bin estimator of thermodynamical quantities------*/
  mb_t      *mb;      // multiple-bin estimator

	/*-----Thermodynamical staffs------------------------------------*/
  double    Ea;       // total potential energy

  /*-----Parallel staffs-------------------------------------------*/
  int mpi_rank, mpi_size;
#ifdef GMX_MPI
  MPI_Comm mpi_comm;
#endif
}; /* types for accelerated tempering */

/* Returns the current beta */
real AdaptTempering_CurrentBeta(at_t *at);

/* Returns the current temperature */
real AdaptTempering_CurrentT(at_t *at);

/* Returns the current paramters in macrocanonical ensemble. */
real AdaptTempering_CurrentPara(at_t *at);

/* Returns the current 2nd paramters in macrocanonical ensemble. */
real AdaptTempering_CurrentSecondPara(at_t *at);

/* Returns the reference temperature. Should be always equal to the thermostat temperature. */
real AdaptTempering_ReferenceTemperature(at_t *at);

/* Force update beta without doing Langevin equation. Might be called during initialization in the beginning of simulation or after every parameter exchange. */
void AdaptTempering_ForceChangeBeta(at_t *at, double newbeta);

/* Create at_t *at for master nodes. */
at_t *AdaptTempering_MasterCreate(const char *fname, bool bCPT, double tmstep, int suffix);

/* Create at_t *at for non-master nodes. */
at_t *AdaptTempering_NonMasterCreate();

/* Open log files. */
void AdaptTempering_OpenLog(at_t *at);

/* Dump all the information of at_t *at to file. */
int AdaptTempering_DumpToFile(at_t *at, const char *fname, int arrmax);

/* Integrate Langevin equation and update the temperature. */
int AdaptTempering_Langevin(at_t *at, llong_t step, bool bfirst, bool blast, bool btr, bool bflush);

/* Update the current temperature according to the multiple bin data. */
void AdaptTempering_UpdateTemperature(at_t *at);

/* Get the rescaling factor of forces */
real AdaptTempering_ForceScaleFactor(at_t *at);

/* Close at_t *at and clean data. */
void AdaptTempering_Close(at_t *at);

/* Sync the information of at_t *at to all nodes. */
#ifdef GMX_MPI
int AdaptTempering_SyncAllNodes(at_t *at, MPI_Comm comm);
#endif

#undef ZCOM_PICK
#undef ZCOM_LOG
#endif
