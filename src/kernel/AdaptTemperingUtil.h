#ifndef ADAPTIVE_TEMPERING_UTILTIES___
#define ADAPTIVE_TEMPERING_UTILTIES___

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "vsite.h"

#include "AdaptTempering.h"

typedef struct AdaptTemperingParaExchangeType at_repl_ex_t;

char *AdaptTemperingGetCfgFileName(char *arg, int nfile, const t_filenm fnm[], t_commrec *cr);

at_t *AdaptTemperingInit(char *cfg_fnm, gmx_bool bCPT, t_inputrec *ir, t_commrec *cr);

gmx_bool AdaptTemperingDoTemperingThisStep(at_t *at, gmx_large_int_t step,
		gmx_bool bFirstStep, gmx_bool bLastStep, gmx_bool bNS, t_commrec *cr);

gmx_bool AdaptTemperingUpdate(at_t *at, gmx_large_int_t step, 
		gmx_bool bTempering, gmx_bool bFirstStep, gmx_bool bLastStep, 
		gmx_bool bTotE, gmx_bool bXTC, t_commrec *cr, gmx_enerdata_t *enerd);
			
at_repl_ex_t *AdaptTemperingInitParaExchange(at_t *at, FILE *fplog, 
		gmx_multisim_t *ms, t_state *state, t_inputrec *ir, int nst, int seed);
							
gmx_bool AdaptTemperingDoParaExchange(at_t *at, FILE *fplog, t_commrec *cr,
		at_repl_ex_t *re, t_state *state, gmx_enerdata_t *enerd, t_state *state_local,
		gmx_large_int_t step, real time);

void AdaptTemperingPrintExchangeStatistics(at_t * at, FILE *fplog, at_repl_ex_t *re);

gmx_bool AdaptTemperingCalcForce(FILE *fplog, t_commrec *cr, 
		t_inputrec *inputrec, gmx_large_int_t step, t_nrnb *nrnb, 
		gmx_wallcycle_t wcycle, gmx_localtop_t *top, gmx_mtop_t *mtop,
		gmx_groups_t *groups, matrix box, rvec x[], history_t *hist, rvec f[], 
		tensor vir_force, t_mdatoms *mdatoms, gmx_enerdata_t *enerd, t_fcdata *fcd, 
		real *lambda, t_graph *graph, t_forcerec *fr, gmx_vsite_t *vsite, 
		rvec mu_tot, double t, FILE *field, gmx_edsam_t ed, gmx_bool bBornRadii,
		int flags, at_t *at, gmx_bool bFirstStep);

#endif
