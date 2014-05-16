#ifndef MULTOP_H__
#define MULTOP_H__

#include "typedefs.h"
#include "types/commrec.h"

#define MT_MAX_BLOCK 10 /* max block per topology */
#define MT_MAX_ATOM_PER_MOL 100000 /* max atom # per molecule */
#define MT_MAX_PDIHS 10000 /* max proper dihedrals */

typedef struct MulTop_GlobalTops_t mt_gtops_t;
typedef struct MulTop_LocalTops_t  mt_ltops_t;

/* Get the names of input topologies. Return # of total topology. */
int MulTop_Global_GetInputFileName(char **fns[], const char *opt, int nfile, const t_filenm fnm[], t_commrec *cr);

/* Initializaion. Skip gtops->gr for non-master nodes. */ 
mt_gtops_t *MulTop_Global_Init(int ntop, real Tref, real Tmax, real Wmax, gmx_bool bMASTER);

/* Set reference topology */
void MulTop_Global_SetReferenceTopology(mt_gtops_t *gtops, gmx_mtop_t *mtop);

/* Get non-reference topologies from states */
void MulTop_Global_GetOtherTopologies(mt_gtops_t *gtops, char **files, t_commrec *cr, t_state **states);

/* Load gtops->gr from file. */
void MulTop_Global_LoadData(mt_gtops_t *gtops);

/* Calculate gtops->gr based on the current states. */
void MulTop_Global_CalcData(mt_gtops_t *gtops, t_state **states_tot, t_state *state_ref);

/* Reconstruct all topologies with a common format. After that, save and clean gtops->gr data. */ 
void MulTop_Global_RefreshForceFieldParameters(mt_gtops_t *gtops);

/* Broadcast gtops information to non-master node. */
void MulTop_Global_Bcast(mt_gtops_t *gtops, t_commrec *cr);

/* Initialize ltops from gtops */
mt_ltops_t *MulTop_Local_Init(mt_gtops_t *gtops);

/* Set local reference topology */
void MulTop_Local_SetReferenceTopology(mt_ltops_t *ltops, gmx_localtop_t *top);

/* Get non-reference local topologies */
void MulTop_Local_GetOtherTopologies(mt_ltops_t *ltops, mt_gtops_t *gtops, t_commrec *cr, t_inputrec *ir);

/* Update the final local topology record (after re-partition the system). */
void MulTop_Local_UpdateRecords(mt_ltops_t *ltops, gmx_domdec_t *dd);

/* Init the final topology */
void MulTop_Local_InitFinalTopology(mt_ltops_t *ltops);

/* Update the final topology (every step) */
void MulTop_Local_UpdateFinalTopologyBasic(mt_ltops_t *ltops);

/* Update the final topology force parameters (firststep or after a termperature change). */
void MulTop_Local_UpdateFinalTopologyParameters(mt_ltops_t *ltops, real *pot, real T);

/* Return the final topology (for force calculation). */
gmx_localtop_t *MulTop_Local_GetFinalTopology(mt_ltops_t *ltops);

#endif
