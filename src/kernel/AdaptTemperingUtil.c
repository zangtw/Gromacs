#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "AdaptTempering.h"
#include "AdaptTemperingUtil.h"
#include "typedefs.h"
#include "main.h"
#include "random.h"
#include "gmx_fatal.h"
#include "smalloc.h"
#include "vec.h"
#include "force.h"
#include "domdec.h"
#include "mvdata.h"

#define PROBABILITYCUTOFF 100

struct AdaptTemperingParaExchangeType
{
    int      repl;   
		int      nrepl; 
    real     *q0;
    real     *q1;     
    int      *ind;   
		int      *allswaps;
    int      nst;
    int      seed;
    int      nattempt[2];
    real     *prob_sum;
    int      **nmoves;
    int      *nexchange;

    /* these are helper arrays for replica exchange; allocated here so they
       don't have to be allocated each time */
    int      *destinations;
    int     **cyclic;
    int     **order;
    int      *tmpswap;
    gmx_bool *incycle;
    gmx_bool *bEx;

    /* helper arrays to hold the quantities that are exchanged */
    real  *prob;
    real  *beta;
};

/*################ STATIC FUNCTIONS: PARAMETER EXCHANGE RELATED #########################*/
static void AdaptTemperingCollectQuantities(at_t *at, at_repl_ex_t *re,
		const gmx_multisim_t *ms, real *para_array, real *para2_array,
		real *beta_array)
{
	int i;

	if(para_array)
	{
		for(i=0; i<re->nrepl; i++)
			para_array[i] = 0;
		para_array[re->repl] = AdaptTempering_CurrentPara(at);
		gmx_sum_sim(ms->nsim, para_array, ms);
	}

	if(para2_array)
	{
		for(i=0; i<re->nrepl; i++)
			para2_array[i] = 0;
		para2_array[re->repl] = AdaptTempering_CurrentSecondPara(at);
		gmx_sum_sim(ms->nsim, para2_array, ms);
	}

	if(beta_array)
	{
		for(i=0; i<re->nrepl; i++)
			beta_array[i] = 0;
		beta_array[re->repl] = AdaptTempering_CurrentBeta(at);
		gmx_sum_sim(ms->nsim, beta_array, ms);
	}
}

static void AdaptTemperingComputeExchangeOrder(FILE *fplog, int **cyclic,
		int **order, const int nrepl, const int maxswap)
{
    int i, j;

    for (j = 0; j < maxswap; j++)
    {
        for (i = 0; i < nrepl; i++)
        {
            if (cyclic[i][j+1] >= 0)
            {
                order[cyclic[i][j+1]][j] = cyclic[i][j];
                order[cyclic[i][j]][j]   = cyclic[i][j+1];
            }
        }
        for (i = 0; i < nrepl; i++)
        {
            if (order[i][j] < 0)
            {
                order[i][j] = i; /* if it's not exchanging, it should stay this round*/
            }
        }
    }

    if (debug)
    {
        fprintf(fplog, "Replica Exchange Order\n");
        for (i = 0; i < nrepl; i++)
        {
            fprintf(fplog, "Replica %d:", i);
            for (j = 0; j < maxswap; j++)
            {
                if (order[i][j] < 0)
                {
                    break;
                }
                fprintf(debug, "%2d", order[i][j]);
            }
            fprintf(fplog, "\n");
        }
        fflush(fplog);
    }
}

static void AdaptTemperingCyclicDecomposition(FILE *fplog, 
		const int *destinations, int **cyclic, gmx_bool *incycle, const int nrepl,
		int *nswap)
{

    int i, j, c, p;
    int maxlen = 1;
    for (i = 0; i < nrepl; i++)
    {
        incycle[i] = FALSE;
    }
    for (i = 0; i < nrepl; i++)  /* one cycle for each replica */
    {
        if (incycle[i])
        {
            cyclic[i][0] = -1;
            continue;
        }
        cyclic[i][0] = i;
        incycle[i]   = TRUE;
        c            = 1;
        p            = i;
        for (j = 0; j < nrepl; j++) /* potentially all cycles are part, but we will break first */
        {
            p = destinations[p];    /* start permuting */
            if (p == i)
            {
                cyclic[i][c] = -1;
                if (c > maxlen)
                {
                    maxlen = c;
                }
                break; /* we've reached the original element, the cycle is complete, and we marked the end. */
            }
            else
            {
                cyclic[i][c] = p;  /* each permutation gives a new member of the cycle */
                incycle[p]   = TRUE;
                c++;
            }
        }
    }
    *nswap = maxlen - 1;

    if (debug)
    {
        for (i = 0; i < nrepl; i++)
        {
            fprintf(debug, "Cycle %d:", i);
            for (j = 0; j < nrepl; j++)
            {
                if (cyclic[i][j] < 0)
                {
                    break;
                }
                fprintf(debug, "%2d", cyclic[i][j]);
            }
            fprintf(debug, "\n");
        }
        fflush(debug);
    }
}

static void AdaptTemperingExchangePrepare(FILE *fplog,
		const int *destinations, const int replica_id, const int nrepl, int *maxswap,
		int **order, int **cyclic, int *incycle, gmx_bool *bThisReplicaExchanged)
{
  int i, j;
  /* Hold the cyclic decomposition of the (multiple) replica
   * exchange. */
  gmx_bool bAnyReplicaExchanged = FALSE;
  *bThisReplicaExchanged = FALSE;

  for (i=0; i<nrepl; i++)
  {
		if (destinations[i] != i)
		{
			/* only mark as exchanged if the index has been shuffled */
			bAnyReplicaExchanged = TRUE;
			break;
		}
  }
  if (bAnyReplicaExchanged)
  {
		/* reinitialize the placeholder arrays */
		for (i=0; i<nrepl; i++)
		{
			for (j = 0; j < nrepl; j++)
			{
				cyclic[i][j] = -1;
				order[i][j]  = -1;
			}
		}

    /* Identify the cyclic decomposition of the permutation (very
     * fast if neighbor replica exchange). */
    AdaptTemperingCyclicDecomposition(fplog, destinations, cyclic, incycle, nrepl, maxswap);

    /* Now translate the decomposition into a replica exchange
     * order at each step. */
    AdaptTemperingComputeExchangeOrder(fplog, cyclic, order, nrepl, *maxswap);

    /* Did this replica do any exchange at any point? */
    for (j = 0; j < *maxswap; j++)
    {
			if (replica_id != order[replica_id][j])
			{
				*bThisReplicaExchanged = TRUE;
				break;
			}
		}
  }
}

static void AdaptTemperingExchangeReals(const gmx_multisim_t *ms, int b, 
		real *v, int n)
{
  real *buf;
  int   i;

  if (v)
  {
		snew(buf, n);
#ifdef GMX_MPI
    /*
       MPI_Sendrecv(v,  n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
       buf,n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
       ms->mpi_comm_masters,MPI_STATUS_IGNORE);
     */
		{
			MPI_Request mpi_req;
			
			MPI_Isend(v, n*sizeof(real), MPI_BYTE, MSRANK(ms, b), 0,
                    ms->mpi_comm_masters, &mpi_req);
			MPI_Recv(buf, n*sizeof(real), MPI_BYTE, MSRANK(ms, b), 0,
                   ms->mpi_comm_masters, MPI_STATUS_IGNORE);
			MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
    }
#endif
		for (i=0; i<n; i++)
		{
			v[i] = buf[i];
		}
		sfree(buf);
  }
}

static void AdaptTemperingExchangeInts(const gmx_multisim_t *ms, int b, 
		int *v, int n)
{
  int *buf;
  int  i;

  if (v)
  {
		snew(buf, n);
#ifdef GMX_MPI
    /*
       MPI_Sendrecv(v,  n*sizeof(int),MPI_BYTE,MSRANK(ms,b),0,
       buf,n*sizeof(int),MPI_BYTE,MSRANK(ms,b),0,
       ms->mpi_comm_masters,MPI_STATUS_IGNORE);
     */
		{
			MPI_Request mpi_req;
			
			MPI_Isend(v, n*sizeof(int), MPI_BYTE, MSRANK(ms, b), 0,
                    ms->mpi_comm_masters, &mpi_req);
			MPI_Recv(buf, n*sizeof(int), MPI_BYTE, MSRANK(ms, b), 0,
                   ms->mpi_comm_masters, MPI_STATUS_IGNORE);
			MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
		}
#endif
		for (i=0; i<n; i++)
		{
			v[i] = buf[i];
		}
		sfree(buf);
  }
}

static void AdaptTemperingExchangeDoubles(const gmx_multisim_t *ms, int b, 
		double *v, int n)
{
  double *buf;
  int     i;

  if (v)
  {
		snew(buf, n);
#ifdef GMX_MPI
		/*
       MPI_Sendrecv(v,  n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
       buf,n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
       ms->mpi_comm_masters,MPI_STATUS_IGNORE);
     */
		{
			MPI_Request mpi_req;
			
			MPI_Isend(v, n*sizeof(double), MPI_BYTE, MSRANK(ms, b), 0,
                    ms->mpi_comm_masters, &mpi_req);
			MPI_Recv(buf, n*sizeof(double), MPI_BYTE, MSRANK(ms, b), 0,
                   ms->mpi_comm_masters, MPI_STATUS_IGNORE);
			MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
		}
#endif
		for (i=0; i<n; i++)
		{
			v[i] = buf[i];
		}
		sfree(buf);
  }
}

static void AdaptTemperingExchangeRvecs(const gmx_multisim_t *ms, int b, 
		rvec *v, int n)
{
  rvec *buf;
  int   i;

  if (v)
  {
		snew(buf, n);
#ifdef GMX_MPI
		/*
       MPI_Sendrecv(v[0],  n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
       buf[0],n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
       ms->mpi_comm_masters,MPI_STATUS_IGNORE);
		 */
		{
			MPI_Request mpi_req;
			
			MPI_Isend(v[0], n*sizeof(rvec), MPI_BYTE, MSRANK(ms, b), 0,
                    ms->mpi_comm_masters, &mpi_req);
			MPI_Recv(buf[0], n*sizeof(rvec), MPI_BYTE, MSRANK(ms, b), 0,
                   ms->mpi_comm_masters, MPI_STATUS_IGNORE);
			MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
		}
#endif
		for (i=0; i<n; i++)
		{
			copy_rvec(buf[i], v[i]);
		}
		sfree(buf);
  }
}

static void exchange_reals(const gmx_multisim_t *ms, int b, real *v, int n)
{
    real *buf;
    int   i;

    if (v)
    {
        snew(buf, n);
#ifdef GMX_MPI
        /*
           MPI_Sendrecv(v,  n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
           buf,n*sizeof(real),MPI_BYTE,MSRANK(ms,b),0,
           ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v, n*sizeof(real), MPI_BYTE, MSRANK(ms, b), 0,
                      ms->mpi_comm_masters, &mpi_req);
            MPI_Recv(buf, n*sizeof(real), MPI_BYTE, MSRANK(ms, b), 0,
                     ms->mpi_comm_masters, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
#endif
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}

static void exchange_ints(const gmx_multisim_t *ms, int b, int *v, int n)
{
    int *buf;
    int  i;

    if (v)
    {
        snew(buf, n);
#ifdef GMX_MPI
        /*
           MPI_Sendrecv(v,  n*sizeof(int),MPI_BYTE,MSRANK(ms,b),0,
             buf,n*sizeof(int),MPI_BYTE,MSRANK(ms,b),0,
             ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v, n*sizeof(int), MPI_BYTE, MSRANK(ms, b), 0,
                      ms->mpi_comm_masters, &mpi_req);
            MPI_Recv(buf, n*sizeof(int), MPI_BYTE, MSRANK(ms, b), 0,
                     ms->mpi_comm_masters, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
#endif
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}

static void exchange_doubles(const gmx_multisim_t *ms, int b, double *v, int n)
{
    double *buf;
    int     i;

    if (v)
    {
        snew(buf, n);
#ifdef GMX_MPI
        /*
           MPI_Sendrecv(v,  n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
           buf,n*sizeof(double),MPI_BYTE,MSRANK(ms,b),0,
           ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v, n*sizeof(double), MPI_BYTE, MSRANK(ms, b), 0,
                      ms->mpi_comm_masters, &mpi_req);
            MPI_Recv(buf, n*sizeof(double), MPI_BYTE, MSRANK(ms, b), 0,
                     ms->mpi_comm_masters, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
#endif
        for (i = 0; i < n; i++)
        {
            v[i] = buf[i];
        }
        sfree(buf);
    }
}

static void exchange_rvecs(const gmx_multisim_t *ms, int b, rvec *v, int n)
{
    rvec *buf;
    int   i;

    if (v)
    {
        snew(buf, n);
#ifdef GMX_MPI
        /*
           MPI_Sendrecv(v[0],  n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
           buf[0],n*sizeof(rvec),MPI_BYTE,MSRANK(ms,b),0,
           ms->mpi_comm_masters,MPI_STATUS_IGNORE);
         */
        {
            MPI_Request mpi_req;

            MPI_Isend(v[0], n*sizeof(rvec), MPI_BYTE, MSRANK(ms, b), 0,
                      ms->mpi_comm_masters, &mpi_req);
            MPI_Recv(buf[0], n*sizeof(rvec), MPI_BYTE, MSRANK(ms, b), 0,
                     ms->mpi_comm_masters, MPI_STATUS_IGNORE);
            MPI_Wait(&mpi_req, MPI_STATUS_IGNORE);
        }
#endif
        for (i = 0; i < n; i++)
        {
            copy_rvec(buf[i], v[i]);
        }
        sfree(buf);
    }
}

static void AdaptTemperingExchangeState(at_t *at, 
		at_repl_ex_t *re, const gmx_multisim_t *ms, 
		int b, t_state *state)
{
    /* Exchange the state */
    int ngtc, nnhpres;
    ngtc    = state->ngtc * state->nhchainlength;
    nnhpres = state->nnhpres* state->nhchainlength;
    exchange_rvecs(ms, b, state->box, DIM);
    exchange_rvecs(ms, b, state->box_rel, DIM);
    exchange_rvecs(ms, b, state->boxv, DIM);
    exchange_reals(ms, b, &(state->veta), 1);
    exchange_reals(ms, b, &(state->vol0), 1);
    exchange_rvecs(ms, b, state->svir_prev, DIM);
    exchange_rvecs(ms, b, state->fvir_prev, DIM);
    exchange_rvecs(ms, b, state->pres_prev, DIM);
    exchange_doubles(ms, b, state->nosehoover_xi, ngtc);
    exchange_doubles(ms, b, state->nosehoover_vxi, ngtc);
    exchange_doubles(ms, b, state->nhpres_xi, nnhpres);
    exchange_doubles(ms, b, state->nhpres_vxi, nnhpres);
    exchange_doubles(ms, b, state->therm_integral, state->ngtc);
    exchange_rvecs(ms, b, state->x, state->natoms);
    exchange_rvecs(ms, b, state->v, state->natoms);
    exchange_rvecs(ms, b, state->sd_X, state->natoms);
		
		/* Exchange the beta */
		double mybeta=(double)(re->beta[re->repl]);
		exchange_doubles(ms, b, &mybeta, 1);
		AdaptTempering_ForceChangeBeta(at, mybeta);
}

static void copy_rvecs(rvec *s, rvec *d, int n)
{
    int i;

    if (d != NULL)
    {
        for (i = 0; i < n; i++)
        {
            copy_rvec(s[i], d[i]);
        }
    }
}

static void copy_doubles(const double *s, double *d, int n)
{
    int i;

    if (d != NULL)
    {
        for (i = 0; i < n; i++)
        {
            d[i] = s[i];
        }
    }
}

static void copy_reals(const real *s, real *d, int n)
{
    int i;

    if (d != NULL)
    {
        for (i = 0; i < n; i++)
        {
            d[i] = s[i];
        }
    }
}

static void copy_ints(const int *s, int *d, int n)
{
    int i;

    if (d != NULL)
    {
        for (i = 0; i < n; i++)
        {
            d[i] = s[i];
        }
    }
}

#define scopy_rvecs(v, n)   copy_rvecs(state->v, state_local->v, n);
#define scopy_doubles(v, n) copy_doubles(state->v, state_local->v, n);
#define scopy_reals(v, n) copy_reals(state->v, state_local->v, n);
#define scopy_ints(v, n)   copy_ints(state->v, state_local->v, n);
static void AdaptTemperingUpdateLocalState(t_state *state, t_state *state_local)
{
    /* When t_state changes, this code should be updated. */
    int ngtc, nnhpres;
    ngtc    = state->ngtc * state->nhchainlength;
    nnhpres = state->nnhpres* state->nhchainlength;
    scopy_rvecs(box, DIM);
    scopy_rvecs(box_rel, DIM);
    scopy_rvecs(boxv, DIM);
    state_local->veta = state->veta;
    state_local->vol0 = state->vol0;
    scopy_rvecs(svir_prev, DIM);
    scopy_rvecs(fvir_prev, DIM);
    scopy_rvecs(pres_prev, DIM);
    scopy_doubles(nosehoover_xi, ngtc);
    scopy_doubles(nosehoover_vxi, ngtc);
    scopy_doubles(nhpres_xi, nnhpres);
    scopy_doubles(nhpres_vxi, nnhpres);
    scopy_doubles(therm_integral, state->ngtc);
    scopy_rvecs(x, state->natoms);
    scopy_rvecs(v, state->natoms);
    scopy_rvecs(sd_X, state->natoms);
    copy_ints(&(state->fep_state), &(state_local->fep_state), 1);
    scopy_reals(lambda, efptNR);
}

static void print_count(FILE *fplog, const char *leg, int n, int *count)
{
    int i;

    fprintf(fplog, "Repl %2s ", leg);
    for (i = 1; i < n; i++)
    {
        fprintf(fplog, " %4d", count[i]);
    }
    fprintf(fplog, "\n");
}

static void print_ind(FILE *fplog, const char *leg, int n, int *ind, gmx_bool *bEx)
{
    int i;

    fprintf(fplog, "Repl %2s %2d", leg, ind[0]);
    for (i = 1; i < n; i++)
    {
        fprintf(fplog, " %c %2d", (bEx != 0 && bEx[i]) ? 'x' : ' ', ind[i]);
    }
    fprintf(fplog, "\n");
}

static void print_prob(FILE *fplog, const char *leg, int n, real *prob)
{
    int  i;
    char buf[8];

    fprintf(fplog, "Repl %2s ", leg);
    for (i = 1; i < n; i++)
    {
        if (prob[i] >= 0)
        {
            sprintf(buf, "%4.2f", prob[i]);
            fprintf(fplog, "  %3s", buf[0] == '1' ? "1.0" : buf+1);
        }
        else
        {
            fprintf(fplog, "     ");
        }
    }
    fprintf(fplog, "\n");
}

static void print_transition_matrix(FILE *fplog, const char *leg, int n, int **nmoves, int *nattempt)
{
    int   i, j, ntot;
    float Tprint;

    ntot = nattempt[0] + nattempt[1];
    fprintf(fplog, "\n");
    fprintf(fplog, "Repl");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "    ");  /* put the title closer to the center */
    }
    fprintf(fplog, "Empirical Transition Matrix\n");

    fprintf(fplog, "Repl");
    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "%8d", (i+1));
    }
    fprintf(fplog, "\n");

    for (i = 0; i < n; i++)
    {
        fprintf(fplog, "Repl");
        for (j = 0; j < n; j++)
        {
            Tprint = 0.0;
            if (nmoves[i][j] > 0)
            {
                Tprint = nmoves[i][j]/(2.0*ntot);
            }
            fprintf(fplog, "%8.4f", Tprint);
        }
        fprintf(fplog, "%3d\n", i);
    }
}

static void AdaptTemperingExchangeTest(at_t *at, FILE *fplog, 
		const gmx_multisim_t *ms, at_repl_ex_t *re, gmx_large_int_t step, 
		real time)
{
  int       m,i,j,a,b,tmp;
  real      delta = 0;
  gmx_bool  bPrint;
  gmx_bool *bEx      = re->bEx;
	real     *beta     = re->beta;
  real     *prob     = re->prob;
  int      *pind     = re->destinations; /* permuted index */

  fprintf(fplog, "Parameter exchange at step " gmx_large_int_pfmt " time %g\n", step, time);

	AdaptTemperingCollectQuantities(at, re, ms, NULL, NULL, beta);

  /* make a duplicate set of indices for shuffling */
  for (i=0; i<re->nrepl; i++)
  {
		pind[i] = re->ind[i];
  }

  /* standard nearest neighbor replica exchange */
  m = (step / re->nst) % 2;
  for (i=1; i<re->nrepl; i++)
  {
		a = re->ind[i-1];
		b = re->ind[i];

		bPrint = (re->repl == a || re->repl == b);
		if (i % 2 == m)
		{
			delta = (beta[b] - beta[a]) * (re->q0[b] - re->q0[a]);
			delta += (beta[b] * beta[b] - beta[a] * beta[a]) * (re->q1[b] - re->q1[a]);
			
			if (bPrint)
			{
				fprintf(fplog, "Repl %d <-> %d  Delta = %10.3e \n", a, b, delta);
				fprintf(fplog, "T[b] = %f, T[a] = %f, T0[b] = %f, T0[a] = %f\n", 1.0/beta[b]/BOLTZ, 1.0/beta[a]/BOLTZ, re->q1[b]/re->q0[b]*(-2.0)/BOLTZ, re->q1[a]/re->q0[a]*(-2.0)/BOLTZ);
				//fprintf(fplog, "beta[b] = %f, beta[a] = %f, q0[b] = %f, q0[a] = %f, q1[b] = %f, q1[a] = %f\n", beta[b], beta[a], re->q0[b], re->q0[a], re->q1[b], re->q1[a]);
			}

			if (delta <= 0)
			{
				/* accepted */
				prob[i] = 1;
				bEx[i]  = TRUE;
			}
			else
			{
				if (delta > PROBABILITYCUTOFF)
				{
					prob[i] = 0;
				}
				else
				{
					prob[i] = exp(-delta);
				}
				
				/* roll a number to determine if accepted */
        bEx[i] = (rando(&(re->seed)) < prob[i]);
      }
			re->prob_sum[i] += prob[i];

			if (bEx[i])
			{
				/* swap these two */
				tmp       = pind[i-1];
				pind[i-1] = pind[i];
				pind[i]   = tmp;
				re->nexchange[i]++;  /* statistics for back compatibility */
			}
		}
		else
		{
			prob[i] = -1;
			bEx[i]  = FALSE;
		}
  }
  /* print some statistics */
  print_ind(fplog, "ex", re->nrepl, re->ind, bEx);
  print_prob(fplog, "pr", re->nrepl, prob);
  fprintf(fplog, "\n");
  re->nattempt[m]++;
  
  /* record which moves were made and accepted */
  for (i = 0; i < re->nrepl; i++)
  {
      re->nmoves[re->ind[i]][pind[i]] += 1;
      re->nmoves[pind[i]][re->ind[i]] += 1;
  }
  fflush(fplog); /* make sure we can see what the last exchange was */
}

/*################ STATIC FUNCTIONS: ADAPTIVE TEMPERING RELATED #########################*/
static void AdaptTemperingCheckTemperatureCoupling(at_t *at, t_inputrec *ir)
{
	int i;

	for(i=0; i<ir->opts.ngtc; i++)
	{
		if(ir->opts.ref_t[i] != AdaptTempering_ReferenceTemperature(at))
			gmx_fatal(FARGS,"Coupling temperature of different groups are different from the one in *.cfg file. \
				         This is not allowed for adaptive tempering. \
								 Modify your *.mdp file and try the simulation again.");
	}
}

/*############################## END OF STATIC FUNCTIONS ###############################*/

char *AdaptTemperingGetCfgFileName(char *arg, int nfile, const t_filenm fnm[], t_commrec *cr)
{
  int i;
	char mysim[2];
  char *fname, *p;

  for (i = 0; i < nfile; i++) 
	{
    if (strcmp(arg, fnm[i].opt) == 0) 
		{
			fname = fnm[i].fns[0];

      if (fnm[i].ftp == efMDP) 
			{
        /* modify the extension from .mdp to .cfg */
        if (strcmp(fname, "grompp.mdp") == 0) 
					return NULL;
				else
				{
					/* founded the *.cfg file */
					if(cr->ms != NULL)
					{
						if (cr->ms->nsim >= 10)
							gmx_fatal(FARGS,"Dont support number of simulations larger than 9.\n");
						sprintf(mysim,"%d",cr->ms->sim);
					}
					else mysim[0]='\0';

					p = strstr(fname, ".cfg.mdp");
				  if(mysim[0] != '\0')
				  {
				    strncpy(p, mysim, 1);
				    p++;
				  }
				  strncpy(p, ".cfg",4);
				  p[4] = '\0';
				}
      }
      return fname;
    }
  }
  return NULL;
}

at_t *AdaptTemperingInit(char *cfg_fnm, gmx_bool bCPT, t_inputrec *ir, 
		t_commrec *cr) 
{
	at_t *at = NULL;

	if (SIMMASTER(cr))
	{
		at = AdaptTempering_MasterCreate(cfg_fnm, bCPT, ir->delta_t);
		if (at == NULL)
			gmx_fatal(FARGS,"Error from node %d during initializing adaptive tempering. This maybe caused by a failure to allocate memory.\n", cr->sim_nodeid);
	}
	else
	{
		at = AdaptTempering_NonMasterCreate();
		if (at == NULL)
			gmx_fatal(FARGS,"Error from node %d during initializing adaptive tempering. This maybe caused by a failure to allocate memory.\n", cr->sim_nodeid);
	}

#ifdef GMX_MPI
	if (PAR(cr))
		AdaptTempering_SyncAllNodes(at, cr->mpi_comm_mygroup);
#endif
	
	/* Check the temperature coupling for all nodes */
	AdaptTemperingCheckTemperatureCoupling(at, ir);

	if (SIMMASTER(cr))
	{
		AdaptTempering_OpenLog(at);
		
		/*if (!AdaptTempering_CreateMultipleBin(at))*/
		/*gmx_fatal(FARGS,"Error during creating multiple bins. This maybe caused by a failure to allocate memory.\n");*/
	}

	/* Update the current temperature; update the force rescaling factor */
	if (SIMMASTER(cr))
		AdaptTempering_UpdateTemperature(at);
	if (PAR(cr))
		gmx_bcast(sizeof(at->beta), &at->beta, cr);
	
	if (SIMMASTER(cr))
		AdaptTempering_DumpToFile(at, "at.manifest", 3);

	return at;
}

gmx_bool AdaptTemperingDoTemperingThisStep(at_t *at, gmx_large_int_t step,
		gmx_bool bFirstStep, gmx_bool bLastStep, gmx_bool bNS, t_commrec *cr)
{
	gmx_bool bTempering;

	if(SIMMASTER(cr))
	{
		/* nsttemp < 0 means do tempering at an NS step */
		bTempering = (at->nsttemp > 0) || bNS || bLastStep;

		/* no tempering during prerun, temperature is fixed */
		if (at->nsttemp > 0 && ((int)step % at->nsttemp) != 0 && !bLastStep)
			bTempering = 0; /* if nsttemp is set, do tempering at a regular interval */
	}
  
	if (PAR(cr))
		gmx_bcast(sizeof(gmx_bool), &bTempering, cr);

	return bTempering;
}

gmx_bool AdaptTemperingUpdate(at_t *at, gmx_large_int_t step, 
		gmx_bool bTempering, gmx_bool bFirstStep, gmx_bool bLastStep, 
		gmx_bool bTotE, gmx_bool bXTC, t_commrec *cr, gmx_enerdata_t *enerd)
{
	if(!bTempering)
		return 0;
	
	if(!bTotE)
		gmx_fatal(FARGS,"Total potential energy is absent when doing tempering.");
	if(SIMMASTER(cr))
		at->Ea = enerd->term[F_EPOT];

  /* change temperature, and regularly write output files */
  if (SIMMASTER(cr)) {
    if (AdaptTempering_Langevin(at, (llong_t)step, bFirstStep, bLastStep, bXTC))
			gmx_fatal(FARGS,"node %d, step: " llong_pfmt ", error during moving master\n", cr->nodeid, step);
  }
	if (SIMMASTER(cr))
		AdaptTempering_UpdateTemperature(at);
	if (PAR(cr))
		gmx_bcast(sizeof(at->beta), &at->beta, cr);

  return 1;
}
			
at_repl_ex_t *AdaptTemperingInitParaExchange(at_t *at, FILE *fplog, 
		gmx_multisim_t *ms, t_state *state, t_inputrec *ir, int nst, int seed)
{
  int  i,j,k;
  at_repl_ex_t *re;
  gmx_bool    bDiff;

  fprintf(fplog, "\nInitializing Parameter Exchange for Adapt Tempering\n");

  if (ms == NULL || ms->nsim == 1)
  {
		gmx_fatal(FARGS, "Nothing to exchange with only one replica, maybe you forgot to set the -multi option of mdrun?");
  }

  snew(re, 1);

  re->repl     = ms->sim;
  re->nrepl    = ms->nsim;
  snew(re->q0, re->nrepl); 
  snew(re->q1, re->nrepl); 
	
	fprintf(fplog, "Repl  There are %d replicas:\n", re->nrepl);

  check_multi_int(fplog, ms, state->natoms, "the number of atoms", FALSE);
  check_multi_int(fplog, ms, ir->eI, "the integrator", FALSE);
  check_multi_large_int(fplog, ms, ir->init_step+ir->nsteps, "init_step+nsteps", FALSE);
  check_multi_large_int(fplog, ms, (ir->init_step+nst-1)/nst,
                        "first exchange step: init_step/-replex", FALSE);
  check_multi_int(fplog, ms, ir->etc, "the temperature coupling", FALSE);
  check_multi_int(fplog, ms, ir->opts.ngtc,
                  "the number of temperature coupling groups", FALSE);
  check_multi_int(fplog, ms, ir->epc, "the pressure coupling", FALSE);
  check_multi_int(fplog, ms, ir->efep, "free energy", FALSE);
  check_multi_int(fplog, ms, ir->fepvals->n_lambda, "number of lambda states", FALSE);

	AdaptTemperingCollectQuantities(at, re, ms, re->q0, NULL, NULL);
	AdaptTemperingCollectQuantities(at, re, ms, NULL, re->q1, NULL);

  bDiff = FALSE;
  for (i=1; i<ms->nsim; i++)
  {
		if (re->q0[i] != re->q0[0] || re->q1[i] != re->q1[0])
    {
			bDiff = TRUE;
    }
  }

	if (!bDiff)
		gmx_fatal(FARGS,"The parameters of the %d systems are all the same, there is nothing to exchange",re->nrepl);
      
  /* Make an index for increasing replica order */
  /* only makes sense if one or the other is varying, not both!
     if both are varying, we trust the order the person gave. */
  snew(re->ind, re->nrepl);
  for(i=0; i<re->nrepl; i++)
  {
    re->ind[i] = i;
  }

	/* Sort: currently use q0/q1(propotional to minus beta) as the index sorting criteria */
  for(i=0; i<re->nrepl; i++)
  {
		for(j=i+1; j<re->nrepl; j++)
    {
			if (re->q0[re->ind[j]]/re->q1[re->ind[j]] < re->q0[re->ind[i]]/re->q1[re->ind[i]])
			{
				k = re->ind[i];
				re->ind[i] = re->ind[j];
				re->ind[j] = k;
			}
		}
  }

  /* keep track of all the swaps, starting with the initial placement. */
  snew(re->allswaps, re->nrepl);
  for(i=0; i<re->nrepl; i++)
  {
      re->allswaps[i] = re->ind[i];
  }

	fprintf(fplog, "\nParameter exchange in adaptive tempering\n");
  for(i=0; i<re->nrepl; i++)
	{
		fprintf(fplog, " %5.1f", re->q0[re->ind[i]]);
	}
	fprintf(fplog, "\n");
  
	re->nst = nst;
  if (seed == -1)
  {
		if (MASTERSIM(ms))
		{
			re->seed = make_seed();
    }
		else
		{
			re->seed = 0;
		}
		gmx_sumi_sim(1, &(re->seed), ms);
  }
  else
  {
      re->seed = seed;
  }
  fprintf(fplog, "\nParameter exchange interval: %d\n", re->nst);
  fprintf(fplog, "\nParameter random seed: %d\n", re->seed);

  re->nattempt[0] = 0;
  re->nattempt[1] = 0;

  snew(re->prob_sum, re->nrepl);
  snew(re->nexchange, re->nrepl);
  snew(re->nmoves, re->nrepl);
  for (i=0; i<re->nrepl; i++)
  {
		snew(re->nmoves[i], re->nrepl);
  }
  fprintf(fplog, "Replica exchange information below: x=exchange, pr=probability\n");

  /* generate space for the helper functions so we don't have to snew each time */

  snew(re->destinations, re->nrepl);
  snew(re->incycle, re->nrepl);
  snew(re->tmpswap, re->nrepl);
  snew(re->cyclic, re->nrepl);
  snew(re->order, re->nrepl);
  for (i=0; i<re->nrepl; i++)
  {
		snew(re->cyclic[i], re->nrepl);
    snew(re->order[i], re->nrepl);
  }
  /* allocate space for the functions storing the data for the replicas */
  /* not all of these arrays needed in all cases, but they don't take
     up much space, since the max size is nrepl**2 */
  snew(re->prob, re->nrepl);
  snew(re->bEx, re->nrepl);
  snew(re->beta, re->nrepl);
  return re;
}
							
gmx_bool AdaptTemperingDoParaExchange(at_t *at, FILE *fplog, t_commrec *cr,
		at_repl_ex_t *re, t_state *state, gmx_enerdata_t *enerd, t_state *state_local,
		gmx_large_int_t step, real time)
{
  int i, j;
  int replica_id = 0;
  int exchange_partner;
  int maxswap = 0;
  /* Number of rounds of exchanges needed to deal with any multiple
   * exchanges. */
  /* Where each replica ends up after the exchange attempt(s). */
  /* The order in which multiple exchanges will occur. */
  gmx_bool bThisReplicaExchanged = FALSE;

  if (MASTER(cr))
  {
      replica_id  = re->repl;
      AdaptTemperingExchangeTest(at, fplog, cr->ms, re, step, time);
      AdaptTemperingExchangePrepare(fplog, re->destinations, 
					replica_id, re->nrepl, &maxswap, re->order, re->cyclic, re->incycle, 
					&bThisReplicaExchanged);
  }
  /* Do intra-simulation broadcast so all processors belonging to
   * each simulation know whether they need to participate in
   * collecting the state. Otherwise, they might as well get on with
   * the next thing to do. */
  if (PAR(cr))
  {
#ifdef GMX_MPI
    MPI_Bcast(&bThisReplicaExchanged, sizeof(gmx_bool), MPI_BYTE, MASTERRANK(cr),
                cr->mpi_comm_mygroup);
#endif
  }

  if (bThisReplicaExchanged)
  {
		/* Exchange the states */
		if (PAR(cr))
		{
			/* Collect the global state on the master node */
			if (DOMAINDECOMP(cr))
			{
				dd_collect_state(cr->dd, state_local, state);
			}
    }

		if (MASTER(cr))
    {
       /* There will be only one swap cycle with standard replica
        * exchange, but there may be multiple swap cycles if we
        * allow multiple swaps. */

			for (j=0; j<maxswap; j++)
			{
				exchange_partner = re->order[replica_id][j];

				if (exchange_partner != replica_id)
				{
					/* Exchange the global states between the master nodes */
					fprintf(fplog, "Exchanging %d with %d\n", replica_id, exchange_partner);
					AdaptTemperingExchangeState(at, re, cr->ms, exchange_partner, state);
        }
      }
    }
      
		/* With domain decomposition the global state is distributed later */
    if (!DOMAINDECOMP(cr))
		{
			/* Copy the global state to the local state data structure */
			AdaptTemperingUpdateLocalState(state, state_local);

			if (PAR(cr))
			{
				bcast_state(cr, state, FALSE);
			}
    }
  }

  return bThisReplicaExchanged;
}

void AdaptTemperingPrintExchangeStatistics(at_t * at, FILE *fplog, at_repl_ex_t *re)
{
    int  i;

    fprintf(fplog, "\nParameter exchange statistics\n");

    fprintf(fplog, "Repl  %d attempts, %d odd, %d even\n",
            re->nattempt[0]+re->nattempt[1], re->nattempt[1], re->nattempt[0]);

    fprintf(fplog, "Repl  average probabilities:\n");
    for (i = 1; i < re->nrepl; i++)
    {
        if (re->nattempt[i%2] == 0)
        {
            re->prob[i] = 0;
        }
        else
        {
            re->prob[i] =  re->prob_sum[i]/re->nattempt[i%2];
        }
    }
    print_ind(fplog, "", re->nrepl, re->ind, NULL);
    print_prob(fplog, "", re->nrepl, re->prob);

    fprintf(fplog, "Repl  number of exchanges:\n");
    print_ind(fplog, "", re->nrepl, re->ind, NULL);
    print_count(fplog, "", re->nrepl, re->nexchange);

    fprintf(fplog, "Repl  average number of exchanges:\n");
    for (i = 1; i < re->nrepl; i++)
    {
        if (re->nattempt[i%2] == 0)
        {
            re->prob[i] = 0;
        }
        else
        {
            re->prob[i] =  ((real)re->nexchange[i])/re->nattempt[i%2];
        }
    }
    print_ind(fplog, "", re->nrepl, re->ind, NULL);
    print_prob(fplog, "", re->nrepl, re->prob);

    fprintf(fplog, "\n");
    
    /* print the transition matrix */
    print_transition_matrix(fplog, "", re->nrepl, re->nmoves, re->nattempt);
}

gmx_bool AdaptTemperingCalcForce(FILE *fplog, t_commrec *cr, 
		t_inputrec *inputrec, gmx_large_int_t step, t_nrnb *nrnb, 
		gmx_wallcycle_t wcycle, gmx_localtop_t *top, gmx_mtop_t *mtop,
		gmx_groups_t *groups, matrix box, rvec x[], history_t *hist, rvec f[], 
		tensor vir_force, t_mdatoms *mdatoms, gmx_enerdata_t *enerd, t_fcdata *fcd, 
		real *lambda, t_graph *graph, t_forcerec *fr, gmx_vsite_t *vsite, 
		rvec mu_tot, double t, FILE *field, gmx_edsam_t ed, gmx_bool bBornRadii,
		int flags, at_t *at, gmx_bool bFirstStep)
{
	gmx_bool bUpdated;

	/* For future use */
	/* if(bUpdated || bFirstStep)
	 * {
	 *
	 * }
	 */

	do_force(fplog, cr, inputrec, step, nrnb, wcycle, top, mtop, groups, box, x, hist,
           f, vir_force, mdatoms, enerd, fcd, lambda, graph, fr, vsite, mu_tot, t, 
	  			 field, ed, bBornRadii, flags);

	/* Rescale the forces */
  if (at != NULL) {
		int i;
		real scale;
		
		scale = AdaptTempering_ForceScaleFactor(at);
    
    /* scale the force */
    for (i = mdatoms->start; i < mdatoms->start + mdatoms->homenr; i++) {
      f[i][0] *= scale; 
      f[i][1] *= scale; 
      f[i][2] *= scale; 
    }
  }

	return bUpdated;
}
