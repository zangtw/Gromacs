#include "typedefs.h"
#include "gmx_fatal.h"
#include "physics.h"
#include "macros.h"
#include "smalloc.h"
#include "string.h"
#include "vec.h"
#include "tpxio.h"
#include "mvdata.h"
#include "network.h"
#include "domdec.h"
#include "filenm.h"
#include "mtop_util.h"
#include "pbc.h"
#include "bondf.h"
#include "MultiTopo.h"

/* Records of the relationship between various global topologies.
 * There is only ONE reference topology, which is the topology 0.
 * The molecule/atom sets in any other topologies MUST be subsets of the ones in the reference topology.
 * This type is only available on the master node. 
 * Information will be saved as file and discarded after finishing initialization and loaded
 * from file at continuation. */
typedef struct MulTop_GlobalRecords_t{

	/* Objective: Find the corresponding atom index in the reference topology for an atom in other topologies.
	 * Usage: find_atom[i][j], where i=topology index, j=atom index.
	 * Properties: find_atom[0][j] always gives j. */
	int **find_atom; 

	/* Objective: Find the corresponding molecule block index between topologies.
	 * Usage: find_block[i][j][k], where i=current topology index, j=current molblock index, k=target topology index.
	 * Properties:
	 * 1. find_block[i][j][i] always equals to j.
	 * 2. find_block[i][j][k] always equals to find_block[m][find_block[i][j][m]][k] for any topology index m.
	 * 3. find_block[i][j][k]=-1 tells that no corresponding molblock found.
	 * 4. find_block[i][j][0] should never be -1.
	 * 5. use function "t1b1t2" to return find_block[i][j][k] conveniently.(t1=topology 1, b1=block1, t2=topology 2) */
	int ***find_block;

	/* name of the file recording the information above */
	char file_nm[10];
}MulTop_GlobalRecords;

/* Records of the relationship between local and global topology.
 * This type is available on every node. 
 * Information will be saved as file after finishing initialization and loaded
 * from file at continuation.
 * LocalRecords are not used in this version, where non-bonded interaction is not considered in additional topologies. It will be used in the next version. */
typedef struct MulTop_LocalRecords_t{
	
	/* Objective: Find the corresponding atom index in the local topology for a global topology.
	 * Usage: index_g2l[i], i=global atom index.
	 * Properties: Result=-1 means that no such atom in the local topology. */
	int *index_g2l;
		
	/* Objective: Find the corresponding atom index in the global topology for a local topology.
	 * Usage: index_l2g[i], i=local atom index.
	 * Properties: Result cannot be -1. */
	int *index_l2g;
}MulTop_LocalRecords;

/* Records of some parameters (related to  dihedral angels) to avoid repeat calculation. */
typedef struct MulTop_DihedralRecords_t{
	real **PDIHcos, **PDIHsin;
	real **PIDIHcos, **PIDIHsin;
	int PDIH_start_index, PIDIH_start_index; /* start index of periodic proper/improper dihedral interaction in ffparams.functype and ffparams.iparams */
}MulTop_DihedralRecords;

struct MulTop_GlobalTops_t{
	int ntop; /* # of topologies */
	real Tref; /* when T>Tref, the additionaled topologies will be added. */
	real Tmax; /* largest reachable temperature */
	real Wmax; /* largest weight of additional topologies */
	gmx_mtop_t **mtops; /* all the global topologies */
	MulTop_DihedralRecords *dr;
	MulTop_GlobalRecords   *gr;
}; /* global topologies */

struct MulTop_LocalTops_t{
	int ntop; /* # of topologies */
	real Tref; /* when T>Tref, the additionaled topologies will be added. */
	real Tmax; /* largest reachable temperature */
	real Wmax; /* largest weight of additional topologies */
	gmx_localtop_t **tops; /* all the local topologies */
	gmx_localtop_t *top_final; /* final topology for force calculation */
	MulTop_DihedralRecords *dr;
	MulTop_LocalRecords    *lr;
}; /* local topologies */
	
/* For convenience */
#define GREC_FA     gr->find_atom
#define GREC_FB     gr->find_block
#define LREC_L2G    lr->index_l2g
#define LREC_G2L    lr->index_g2l
#define DREC_Pcos   dr->PDIHcos
#define DREC_Psin   dr->PDIHsin
#define DREC_PIcos  dr->PIDIHcos
#define DREC_PIsin  dr->PIDIHsin
#define DREC_Pind   dr->PDIH_start_index
#define DREC_PIind  dr->PIDIH_start_index

/*################################################
 * for debug use */

static void PDIH_manifest(mt_gtops_t *gtops)
{
	FILE *PDIHcos_manifest0 = fopen("./debug/Pcos_0.txt","w");
	FILE *PDIHsin_manifest0 = fopen("./debug/Psin_0.txt","w");
	FILE *PDIHcos_manifest1 = fopen("./debug/Pcos_1.txt","w");
	FILE *PDIHsin_manifest1 = fopen("./debug/Psin_1.txt","w");
	int i,j;
	MulTop_DihedralRecords *dr = gtops->dr;

	for(i=DREC_Pind; gtops->mtops[0]->ffparams.functype[i] == F_PDIHS ; i++)
	{
		fprintf(PDIHcos_manifest0,"PDIH[%d]: phiA=%f\tcpA=%f\tmult=%d\n",i,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.phiA,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.cpA,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.mult);
		fprintf(PDIHcos_manifest0,"Res=%lf\n",DREC_Pcos[i-DREC_Pind][0]);
		fprintf(PDIHsin_manifest0,"PDIH[%d]: phiA=%f\tcpA=%f\tmult=%d\n",i,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.phiA,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.cpA,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.mult);
		fprintf(PDIHsin_manifest0,"Res=%lf\n",DREC_Psin[i-DREC_Pind][0]);
		fprintf(PDIHcos_manifest1,"PDIH[%d]: phiA=%f\tcpA=%f\tmult=%d\n",i,
				   gtops->mtops[1]->ffparams.iparams[i].pdihs.phiA,
				   gtops->mtops[1]->ffparams.iparams[i].pdihs.cpA,
				   gtops->mtops[1]->ffparams.iparams[i].pdihs.mult);
		fprintf(PDIHcos_manifest1,"Res=%lf\n",DREC_Pcos[i-DREC_Pind][1]);
		fprintf(PDIHsin_manifest1,"PDIH[%d]: phiA=%f\tcpA=%f\tmult=%d\n",i,
				   gtops->mtops[1]->ffparams.iparams[i].pdihs.phiA,
				   gtops->mtops[1]->ffparams.iparams[i].pdihs.cpA,
				   gtops->mtops[1]->ffparams.iparams[i].pdihs.mult);
		fprintf(PDIHsin_manifest1,"Res=%lf\n",DREC_Psin[i-DREC_Pind][1]);
	}

	fclose(PDIHcos_manifest0);
	fclose(PDIHsin_manifest0);
	fclose(PDIHcos_manifest1);
	fclose(PDIHsin_manifest1);
}

static void PIDIH_manifest(mt_gtops_t *gtops)
{
	FILE *PIDIHcos_manifest0 = fopen("./debug/PIcos_0.txt","w");
	FILE *PIDIHsin_manifest0 = fopen("./debug/PIsin_0.txt","w");
	FILE *PIDIHcos_manifest1 = fopen("./debug/PIcos_1.txt","w");
	FILE *PIDIHsin_manifest1 = fopen("./debug/PIsin_1.txt","w");
	int i;
	MulTop_DihedralRecords *dr = gtops->dr;

	for(i=DREC_PIind; gtops->mtops[0]->ffparams.functype[i] == F_PIDIHS ; i++)
	{
		fprintf(PIDIHcos_manifest0,"PDIH[%d]: phiA=%f\tcpA=%f\tmult=%d\n",i,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.phiA,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.cpA,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.mult);
		fprintf(PIDIHcos_manifest0,"Res=%lf\n",DREC_PIcos[i-DREC_PIind][0]);
		fprintf(PIDIHsin_manifest0,"PDIH[%d]: phiA=%f\tcpA=%f\tmult=%d\n",i,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.phiA,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.cpA,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.mult);
		fprintf(PIDIHsin_manifest0,"Res=%lf\n",DREC_PIsin[i-DREC_PIind][0]);
		fprintf(PIDIHcos_manifest1,"PDIH[%d]: phiA=%f\tcpA=%f\tmult=%d\n",i,
				   gtops->mtops[1]->ffparams.iparams[i].pdihs.phiA,
				   gtops->mtops[1]->ffparams.iparams[i].pdihs.cpA,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.mult);
		fprintf(PIDIHcos_manifest1,"Res=%lf\n",DREC_PIcos[i-DREC_PIind][1]);
		fprintf(PIDIHsin_manifest1,"PDIH[%d]: phiA=%f\tcpA=%f\tmult=%d\n",i,
				   gtops->mtops[1]->ffparams.iparams[i].pdihs.phiA,
				   gtops->mtops[1]->ffparams.iparams[i].pdihs.cpA,
				   gtops->mtops[0]->ffparams.iparams[i].pdihs.mult);
		fprintf(PIDIHsin_manifest1,"Res=%lf\n",DREC_PIsin[i-DREC_PIind][1]);
	}

	fclose(PIDIHcos_manifest0);
	fclose(PIDIHsin_manifest0);
	fclose(PIDIHcos_manifest1);
	fclose(PIDIHsin_manifest1);
}

/*################################################
 * static array */

static int *booltable; /* just for speed up the judge_samebonds */

/*################################################
 * static functions */

/* Return find_block[t1][b1][t2].
 * The level of recursion is always < 2. */
static int t1b1t2(int t1, int b1, int t2, MulTop_GlobalRecords *gr)
{
	if (GREC_FB[t1][b1][t2] != -1)
		return GREC_FB[t1][b1][t2];
	else 
	{
		/* find_block[i][j][0] should not be -1. See the "Properties" section above. */
		if (t2 == 0)
		{
			gmx_fatal(FARGS,"error! block %d in your topology %d is not found in reference topology\n", b1, t1);
			return -1;
		}
		else if (t1 == 0)  /* find_block[0][b1][t2] is allowed to be -1. */
			return -1; 
		else /* try to initialize find_block[t1][b1][t2] if it is -1. See the "Properties" section above. */
			return GREC_FB[t1][b1][t2] = t1b1t2(0, t1b1t2(t1, b1, 0, gr), t2, gr);
	}
}

/* fill start_tar with atom indexs corresponding to the ones in start_src in reference topology. The topology index of start_src is n. */
static void fillin_bonds(int *start_tar, int *start_src, int nr, 
		MulTop_GlobalRecords *gr, int n)
{
	int i, j;

	for(i = 0; i < nr; i++)
	{
		j = GREC_FA[n][*(start_src+i)];

		*(start_tar+i) = j;
	}
}

/* judge whether start_tar and start_src represents the same bond */
static gmx_bool judge_samebonds(int *start_tar, int *start_src, 
		MulTop_GlobalRecords *gr, int n, int nr)
{
	int i, j, k;
	gmx_bool same = 1;

	for(i=0; i<nr; i++)
	{
	  j = *(start_tar+i);

	  if (j > MT_MAX_ATOM_PER_MOL)
	    gmx_fatal(FARGS, "# of atoms larger than a preset value! You can change this value in MulTop.h\n");

	  booltable[j] = 1;
	}

	for(j=0; j<nr; j++)
	{
		k = GREC_FA[n][*(start_src+j)];

		if (k > MT_MAX_ATOM_PER_MOL)
			gmx_fatal(FARGS, "# of atoms larger than a preset value! You can change this value in MulTop.h\n");

		if(!booltable[k])
		{
			same = 0;
			break;
		}
	}

	for(i=0; i<nr; i++)
		booltable[*(start_tar+i)] = 0;

	return same;
}

static void iparams_cal(t_iparams *tar, t_iparams *src0, t_iparams *src1,
		int ftype, real *w, int index, real *pot, MulTop_DihedralRecords *dr)
{
	static const real DEG2RAD2 = DEG2RAD * DEG2RAD;
	gmx_bool Zpara0 = !(w[0]);
	gmx_bool Zpara1 = !(w[1]);
	
	switch(ftype)
	{
	case F_LJ:
		if (!Zpara0)
		{
			tar->lj.c6 = src0->lj.c6 * w[0];
			tar->lj.c12 = src0->lj.c12 * w[0];
		}
		if(!Zpara1)
		{
			tar->lj.c6 += src1->lj.c6 * w[1];	
  	  tar->lj.c12 += src1->lj.c12 * w[1];
		}
  	break;
  case F_BONDS: case F_ANGLES: case F_IDIHS:
  	{
			if(Zpara0)
			{
				tar->harmonic.krA = src1->harmonic.krA * w[1];
				tar->harmonic.rA = src1->harmonic.rA;
				tar->harmonic.krB = tar->harmonic.krA;
				tar->harmonic.rB = tar->harmonic.rA;
			}
			else if(Zpara1)
			{
				tar->harmonic.krA = src0->harmonic.krA * w[0];
				tar->harmonic.rA = src0->harmonic.rA;
				tar->harmonic.krB = tar->harmonic.krA;
				tar->harmonic.rB = tar->harmonic.rA;

			}
			else
			{
				tar->harmonic.krA = src0->harmonic.krA * w[0] + src1->harmonic.krA * w[1];
				tar->harmonic.rA = src0->harmonic.krA * src0->harmonic.rA * w[0];
				tar->harmonic.rA += src1->harmonic.krA * src1->harmonic.rA * w[1];

				if(tar->harmonic.krA == 0)
					gmx_fatal(FARGS,"Error! Some bonded parameters are zero both in topology 0 and topology 1. Check your topology settings.");
			
				tar->harmonic.rA /= tar->harmonic.krA;
			
				if(ftype != F_BONDS)
				{
					if(tar->harmonic.rA >= 1800 || tar->harmonic.rA < -1800)
						fprintf(stderr, "Large calculated angle: %lg\n", tar->harmonic.rA);
					
					while (tar->harmonic.rA >= 180)
						tar->harmonic.rA -= 360;
					while (tar->harmonic.rA < -180)
						tar->harmonic.rA += 360;
				}
				
				pot[ftype] += src0->harmonic.krA * src1->harmonic.krA / tar->harmonic.krA 
											* (src0->harmonic.rA - src1->harmonic.rA)
											* (src0->harmonic.rA - src1->harmonic.rA)
					            * w[0] * w[1] * 0.5 * ((ftype == F_BONDS) ? 1 : DEG2RAD2);
				tar->harmonic.rB = tar->harmonic.rA;
				tar->harmonic.krB = tar->harmonic.krA;
			}
  	}
  	break;
  case F_PDIHS:
  {
  	int pdih_ind = index - DREC_Pind;

		if (Zpara0 || !(DREC_Pcos[pdih_ind][0] || DREC_Psin[pdih_ind][0]))
		{
			tar->pdihs.phiA = src1->pdihs.phiA;
			tar->pdihs.phiB = tar->pdihs.phiA;
			tar->pdihs.cpA = src1->pdihs.cpA * w[1];
			tar->pdihs.cpB = tar->pdihs.cpA;
			tar->pdihs.mult = src1->pdihs.mult;
		}
		else if (Zpara1 || !(DREC_Pcos[pdih_ind][1] || DREC_Psin[pdih_ind][1]))
		{
			tar->pdihs.phiA = src0->pdihs.phiA;
			tar->pdihs.phiB = tar->pdihs.phiA;
			tar->pdihs.cpA = src0->pdihs.cpA * w[0];
			tar->pdihs.cpB = tar->pdihs.cpA;
			tar->pdihs.mult = src0->pdihs.mult;
		}
		else
  	{
  		real A = DREC_Pcos[pdih_ind][0] * w[0] + DREC_Pcos[pdih_ind][1] * w[1];
  		real B = DREC_Psin[pdih_ind][0] * w[0] + DREC_Psin[pdih_ind][1] * w[1];
			if(!A)
			{
				tar->pdihs.cpA = B;
				if(B > 0)
					tar->pdihs.phiA = 90.0;
				else if (B < 0)
					tar->pdihs.phiA = -90.0;
				else
					tar->pdihs.phiA = 0;
			}
			else if (!B)
			{
				tar->pdihs.cpA = A;
				tar->pdihs.phiA = 0;
			}
			else	
			{
				tar->pdihs.cpA = sqrt(A * A + B * B);
				tar->pdihs.phiA = RAD2DEG * atan2(B, A); // atan2(y,x)=atan(y/x)
			}

  		tar->pdihs.phiB = tar->pdihs.phiA;
  		tar->pdihs.cpB = tar->pdihs.cpA;
  		tar->pdihs.mult = max(src0->pdihs.mult, src1->pdihs.mult); //Here either they are equal to each other, or one of them is zero (according to the logics in refresh_mtops)
			pot[ftype] += src0->pdihs.cpA * w[0] + src1->pdihs.cpA * w[1] - tar->pdihs.cpA;
  	}
  	break;
  }
  case F_PIDIHS:
  {
		int pdih_ind = index - DREC_PIind;
		if (Zpara0 || !(DREC_PIcos[pdih_ind][0] || DREC_PIsin[pdih_ind][0]))
		{
			tar->pdihs.phiA = src1->pdihs.phiA;
			tar->pdihs.phiB = tar->pdihs.phiA;
			tar->pdihs.cpA = src1->pdihs.cpA;
			tar->pdihs.cpB = tar->pdihs.cpB * w[1];
			tar->pdihs.mult = src1->pdihs.mult;
		}
		if (Zpara1 || !(DREC_PIcos[pdih_ind][1] || DREC_PIsin[pdih_ind][1]))
		{
			tar->pdihs.phiA = src0->pdihs.phiA;
			tar->pdihs.phiB = tar->pdihs.phiA;
			tar->pdihs.cpA = src0->pdihs.cpA * w[0];
			tar->pdihs.cpB = tar->pdihs.cpA;
			tar->pdihs.mult = src0->pdihs.mult;
		}
		else
  	{
  		real A = DREC_PIcos[pdih_ind][0] * w[0] + DREC_PIcos[pdih_ind][1] * w[1];
  		real B = DREC_PIsin[pdih_ind][0] * w[0] + DREC_PIsin[pdih_ind][1] * w[1];
			if(!A)
			{
				tar->pdihs.cpA = B;
				if(B > 0)
					tar->pdihs.phiA = 90.0;
				else if (B < 0)
					tar->pdihs.phiA = -90.0;
				else
					tar->pdihs.phiA = 0;
			}
			else if (!B)
			{
				tar->pdihs.cpA = A;
				tar->pdihs.phiA = 0;
			}
			else	
			{
				tar->pdihs.cpA = sqrt(A * A + B * B);
				tar->pdihs.phiA = RAD2DEG * atan2(B, A); // atan2(y,x)=tan(y/x)
			}

  		tar->pdihs.phiB = tar->pdihs.phiA;
  		tar->pdihs.cpB = tar->pdihs.cpA;
  		tar->pdihs.mult = max(src0->pdihs.mult, src1->pdihs.mult); //Here either they are equal to each other, or one of them is zero (according to the logics in redo_mtops)
			pot[ftype] += src0->pdihs.cpA * w[0] + src1->pdihs.cpA * w[1] - tar->pdihs.cpA;
  	}
  	break;
  }
  case F_RBDIHS: case F_FOURDIHS:
  {
  	int i;
  	for(i=0; i<NR_RBDIHS; i++)
  	{
			if(!Zpara0)
			{
				tar->rbdihs.rbcA[i] = src0->rbdihs.rbcA[i] * w[0];
				tar->rbdihs.rbcB[i] = tar->rbdihs.rbcA[i];
			}
			if(!Zpara1)
			{
				tar->rbdihs.rbcA[i] += src1->rbdihs.rbcA[i] * w[1];
				tar->rbdihs.rbcB[i] = tar->rbdihs.rbcA[i];
			}
  	}
		break;
  }
	case F_LJ14: case F_VDW14:
		if(!Zpara0)
		{
			tar->lj14.c6A = src0->lj14.c6A * w[0];
			tar->lj14.c12A = src0->lj14.c12A * w[0];
			tar->lj14.c6B = tar->lj14.c6A;
			tar->lj14.c12B = tar->lj14.c12A;
		}
		if(!Zpara1)
		{
			tar->lj14.c6A += src1->lj14.c6A * w[1];
			tar->lj14.c12A += src1->lj14.c12A * w[1];
			tar->lj14.c6B = tar->lj14.c6A;
			tar->lj14.c12B = tar->lj14.c12A;
		}
  	break;
  case F_CONSTR:
  	tar->constr.dA = src0->constr.dA;
  	tar->constr.dB = src0->constr.dB;
  	break;
  case F_SETTLE:
  	tar->settle.doh = src0->settle.doh;
  	tar->settle.dhh = src0->settle.dhh;
  	break;
  default:
  	gmx_fatal(FARGS, "currently don't support forces for ftype %d\n",ftype);
  	break;
	}
}

static void iparams_copy(t_iparams *tar, t_iparams *src, int ftype)
{
	switch(ftype)
	{
	case F_LJ:
		tar->lj.c6 = src->lj.c6;
		tar->lj.c12 = src->lj.c12;
		break;
	case F_BONDS: case F_ANGLES: case F_IDIHS:
		tar->harmonic.rA = src->harmonic.rA;
		tar->harmonic.krA = src->harmonic.krA;
		tar->harmonic.rB = src->harmonic.rB;
		tar->harmonic.krB = src->harmonic.krB;
		break;
	case F_PDIHS: case F_PIDIHS:
		tar->pdihs.phiA = src->pdihs.phiA;
		tar->pdihs.cpA = src->pdihs.cpA;
		tar->pdihs.phiB = src->pdihs.phiB;
		tar->pdihs.cpB = src->pdihs.cpB;
		tar->pdihs.mult = src->pdihs.mult;
		break;
	case F_RBDIHS: case F_FOURDIHS:
		{
			int i;
			for(i=0; i<NR_RBDIHS; i++)
			{
				tar->rbdihs.rbcA[i] = src->rbdihs.rbcA[i];
				tar->rbdihs.rbcB[i] = src->rbdihs.rbcB[i];
			}
			break;
		}
	case F_LJ14: case F_VDW14:
		{
			tar->lj14.c6A = src->lj14.c6A;
			tar->lj14.c12A = src->lj14.c12A;
			tar->lj14.c6B = src->lj14.c6B;
			tar->lj14.c12B = src->lj14.c12B;
			break;
		}
	case F_CONSTR:
		tar->constr.dA = src->constr.dA;
		tar->constr.dB = src->constr.dB;
		break;
	case F_SETTLE:
		tar->settle.doh = src->settle.doh;
		tar->settle.dhh = src->settle.dhh;
		break;
	default:
		gmx_fatal(FARGS, "currently don't support forces for ftype %d\n",ftype);
		break;
	}
}

/* Note: gr->find_atom[i] will be re-allocated in MulTop_Global_LoadData or MulTop_Global_CalcData). */
static void MulTop_Global_InitRecords(mt_gtops_t *gtops)
{
	int i, j, k;
	int ntop = gtops->ntop;
	MulTop_GlobalRecords *gr = gtops->gr;

	snew(GREC_FA, ntop);
	for(i=0; i<ntop; i++)
		snew(GREC_FA[i], 1);

	snew(GREC_FB, ntop);
	for(i=0; i<ntop; i++)
	{
		snew(GREC_FB[i], MT_MAX_BLOCK);

		for(j=0; j < MT_MAX_BLOCK; j++)
		{
			snew(GREC_FB[i][j], ntop);

			for(k=0; k<ntop; k++)
				GREC_FB[i][j][k] = -1;
		}
	}

	sprintf(gtops->gr->file_nm, "GREC");
}

static void MulTop_Global_InitDihedral(mt_gtops_t *gtops)
{
	int i;
	int ntop = gtops->ntop;
	MulTop_DihedralRecords *dr = gtops->dr;

	snew(DREC_Pcos, MT_MAX_PDIHS);
	snew(DREC_Psin, MT_MAX_PDIHS);
	snew(DREC_PIcos, MT_MAX_PDIHS);
	snew(DREC_PIsin, MT_MAX_PDIHS);
		
	for(i=0; i < MT_MAX_PDIHS; i++)
	{
		snew(DREC_Pcos[i], ntop);
		snew(DREC_Psin[i], ntop);
		snew(DREC_PIcos[i], ntop);
		snew(DREC_PIsin[i], ntop);
	}
}

static gmx_bool MulTop_Global_OpenFile(mt_gtops_t *gtops, FILE **fp)
{
	FILE *file;
	MulTop_GlobalRecords *gr = gtops->gr;

	if((file = fopen(gr->file_nm, "r")) == NULL)
	{
		if((file = fopen(gr->file_nm, "w")) == NULL)
			gmx_fatal(FARGS,"Error opening file: %s.", gr->file_nm);
		*fp = file;
		return FALSE;
	}

	*fp = file;
	return TRUE; 
}

/* Save the global record to file. */
static void MulTop_Global_SaveData(mt_gtops_t *gtops, int t, FILE *fp)
{
	int i,j,k;
	int natoms = gtops->mtops[t]->natoms;
	MulTop_GlobalRecords *gr = gtops->gr;
	
	/* write header */
	if(t!=0)
		fprintf(fp, " ");
	fprintf(fp, "#");
	fprintf(fp, " %d", t);

	/* write GREC_FA */
	for(i=0; i<natoms; i++)
		fprintf(fp, " %d", GREC_FA[t][i] + 1); /* '+1' here is to make index corresponding to *.gro */

	/* write GREC_FB after last topology is processed */
	if(t == gtops->ntop - 1)
		for(i=0; i< gtops->ntop; i++)
			for(j=0; j< MT_MAX_BLOCK; j++)
				for(k=0; k< gtops->ntop; k++)
					fprintf(fp, " %d", GREC_FB[i][j][k]);
}

/* Save the global record to file. */
static void MulTop_Global_CleanRecords(MulTop_GlobalRecords *gr, int ntop)
{
	int i,j;
	
	/* free GREC_FA */
	for(i=0; i<ntop; i++)
		sfree(GREC_FA[i]);
	sfree(GREC_FA);

	/* clean GREC_FB */
	for(i=0; i< ntop; i++)
		for(j=0; j< MT_MAX_BLOCK; j++)
			sfree(GREC_FB[i][j]);
	for(i=0; i< ntop; i++)
		sfree(GREC_FB[i]);
	sfree(GREC_FB);

	/* clean gr */
	sfree(gr);
}

/* Calculate gr data for the (t)th topology. Note the 0th topology is the reference one. */
static void MulTop_Global_CalcDataForEachTopology(mt_gtops_t *gtops, int t, t_state *state_tar, t_state *state_ref, FILE *fp)
{
  int b0, bi; /* block index for topology ref and tar */
  int m0, mi; /* molecule index within a block */    
  int a0, ai; /* atom index within a molecule  */  
  int bloc0, bloci;  /* location of the current block in mtop->mols.index */
  int nmol0, nmoli;  /* # of molecule in a block  */
  int res0, resi;       /* residue index in the current molecule  */
  char **mnm0, **mnmi;  /* molecule name */
  char *anm0, *anmi;    /* atom name     */
  gmx_bool *blockdone0, *moldone0;
	gmx_mtop_t *mtop0 = gtops->mtops[0];  
	gmx_mtop_t *mtopi = gtops->mtops[t];
	MulTop_GlobalRecords *gr = gtops->gr;

	srenew(GREC_FA[t], mtopi->natoms);

  snew(blockdone0, mtop0->nmolblock);

  /* Begin to find the cooresponding atoms in block bi of the current topology in the reference topology */
  for(bi = 0, bloci = 0, nmoli = 0; bi < mtopi->nmolblock; bi++, bloci+=nmoli)
  {
  	mnmi = mtopi->moltype[bi].name;
 		nmoli = mtopi->molblock[bi].nmol;
  	gmx_bool BlockbFinded = 0;

		/* We need to find the corresponding b0 first.. */
		for(b0 = 0, bloc0 = 0, nmol0 = 0; b0 < mtop0->nmolblock; b0++, bloc0+=nmol0)
		{
			mnm0 = mtop0->moltype[b0].name;
			nmol0 = mtop0->molblock[b0].nmol;		
			
			if(blockdone0[b0]) /* check whether b0 has already been processed. */
				continue;
			
			if(strcmp(mnm0[0],mnmi[0])) /* molecule name must be the same */
				continue;

			BlockbFinded = 1;
			break;
		}

		if(!BlockbFinded)
  		gmx_fatal(FARGS,"Cannot find the corresponding molecule block in the reference topology for the molecule '%s' in topology %d.\n",mnmi[0], t);
  	
		/* number of molecules/atoms CANNOT be larger than the reference topology */
    if(nmol0 < nmoli)
    	gmx_fatal(FARGS,"Number of molecules in molecule block: '%s' in topology %d is larger than the one in the reference topology, which is banned.\n",mnmi[0],t);

    if(mtop0->moltype[b0].atoms.nr < mtopi->moltype[bi].atoms.nr)
     	gmx_fatal(FARGS,"Number of atoms in molecule '%s' in topology %d is larger than the one in the reference topology, which is banned.\n",mnmi[0],t);

		if(b0 > MT_MAX_BLOCK || bi > MT_MAX_BLOCK)
			gmx_fatal(FARGS,"Number of blocks in your topologies is larger than a preset value %d. You can modify this value (named MT_MAX_BLOCK) in MultiTopo.h", MT_MAX_BLOCK);
		GREC_FB[0][b0][t] = bi;
		GREC_FB[t][bi][0] = b0;
		GREC_FB[t][bi][t] = bi;

    snew(moldone0, nmol0);
		
		/* Next, find the corresponding molecules in the current block */
    for(mi = 0; mi < nmoli; mi++)
    {
     	gmx_bool MolbFinded = 0;
    	ai = 0;
    	resi = mtopi->moltype[bi].atoms.atom[ai].resind;
    	anmi = *(mtopi->moltype[bi].atoms.atomname[ai]);

      for(m0 = 0; m0 < nmol0; m0++)
      {
      	if(moldone0[m0])
      		continue;

      	gmx_bool AtombFinded = 0;
      	for(a0 = 0; a0 < mtop0->moltype[b0].atoms.nr; a0++)
				{
       		res0 = mtop0->moltype[b0].atoms.atom[a0].resind;
       		anm0 = *(mtop0->moltype[b0].atoms.atomname[a0]);

					/* Followings are the criteria */

       		if (res0!=resi)
       			continue;

       		/* use fuzzy search, e.g. CG1==CG2 but CG1!=CB1 */
       		if (*anm0 != *anmi)
       			continue;
       		if(*(anm0+1) != *(anmi+1))
       			continue;

       		if(distance2(state_ref->x[mtop0->mols.index[bloc0+m0]+a0],
       								 state_tar->x[mtopi->mols.index[bloci+mi]]) > 1e-6)
       			continue;

    	    AtombFinded = 1;
    	    break;
       	}

       	if(!AtombFinded)
       		continue;

       	MolbFinded = 1;
       	break;
      }

      if(!MolbFinded)
      	gmx_fatal(FARGS,"cannot find the corresponding %dth molecule in block '%s' in topology %d",mi,mnmi[0],t);

      /* Link the first atom in this molecule */
			GREC_FA[t][mtopi->mols.index[bloci+mi]] = mtop0->mols.index[bloc0+m0]+a0;

      /* Now link other atoms in current molecule */
      for(ai = 1; ai < mtopi->moltype[bi].atoms.nr; ai++)
      {
      	resi = mtopi->moltype[bi].atoms.atom[ai].resind;
       	anmi = *(mtopi->moltype[bi].atoms.atomname[ai]);

       	gmx_bool AtombFinded = 0;
       	for(a0 = 0; a0 < mtop0->moltype[b0].atoms.nr; a0++)
				{
					res0 = mtop0->moltype[b0].atoms.atom[a0].resind;
					anm0 = *(mtop0->moltype[b0].atoms.atomname[a0]);

					if (res0!=resi)
						continue;
					if (*anm0 != *anmi)
						continue;
					if(*(anm0+1) != *(anmi+1))
						continue;

					if(distance2(state_ref->x[mtop0->mols.index[bloc0+m0]+a0],
						  				 state_tar->x[mtopi->mols.index[bloci+mi]+ai]) > 1e-6)
						continue;

					AtombFinded = 1;
					break;
       	}

    	  if(!AtombFinded)
    	  	gmx_fatal(FARGS,"cannot find the corresponding atom '%s' of the %dth molecule in block '%s' in topology %d",anmi[0], mi, mnmi[0],t);

    	  GREC_FA[t][mtopi->mols.index[bloci+mi]+ai] = mtop0->mols.index[bloc0+m0]+a0;
      }

      moldone0[m0] = 1;
    }

		sfree(moldone0);

    blockdone0[b0] = 1;
	}

	sfree(blockdone0);

	MulTop_Global_SaveData(gtops, t, fp);
}

/* Currently not useful */
/*#define MA_INDEX (*(int**)((char*)(dd->ma)+16))*/
/*#define MA_CG    (*(int**)((char*)(dd->ma)+24))*/
/*#define MA_NAT   (*(int**)((char*)(dd->ma)+32))*/
/*#define MA_IBUF  (*(int**)((char*)(dd->ma)+40))*/
/*#define COMM_CGS ((t_block*)((char*)(dd->comm)+120))*/

/* Init local records. */ 
static void MulTop_Local_InitRecords(mt_ltops_t *ltops, mt_gtops_t *gtops)
{
	int i;
	int t_atom = gtops->mtops[0]->natoms;
	MulTop_LocalRecords *lr = ltops->lr;

	snew(LREC_L2G, t_atom);
	snew(LREC_G2L, t_atom);

	for(i=0; i<t_atom; i++)
	{
		LREC_L2G[i] = -1;
		LREC_G2L[i] = -1;
	}
}

/* Update the local records for single processor. */ 
static void MulTop_Local_UpdateRecordsSingleProcessor(MulTop_LocalRecords *lr, int natom)
{
	int i;

	for(i=0; i<natom; i++)
	{
		LREC_L2G[i] = i;
		LREC_G2L[i] = i;
	}
}

static real MulTop_Local_CalcBonds(int nbonds, const t_iatom atoms[], const t_iparams params[], const rvec coord[], const t_pbc *pbc)
{
	real vtot, k, x;
	int i, type, ai, aj;
	rvec dd;
	real dr, dr2, dx, dx2;

	for(i=0, vtot=0; i<nbonds;)
	{
		type = atoms[i++];
		ai   = atoms[i++];
		aj   = atoms[i++];
		k    = params[type].harmonic.krA;
		x    = params[type].harmonic.rA;
		
		if(k == 0)
			continue;

		if(pbc)
      pbc_dx_aiuc(pbc, coord[ai], coord[aj], dd);
		else 
      rvec_sub(coord[ai], coord[aj], dd);

    dr2  = iprod(dd, dd);       
    dr   = dr2*gmx_invsqrt(dr2);
		dx   = dr - x;
		dx2  = dx * dx;

		vtot += 0.5 * k * dx2;
	}

	return vtot;
}

static real MulTop_Local_CalcAngles(int nbonds, const t_iatom atoms[], const t_iparams params[], const rvec coord[], const t_pbc *pbc)
{
	real vtot, k, x;
	int i, type, ai, aj, ak;
	real dx, dx2;
	rvec r_ij, r_kj;
	real cos_theta, theta;
	int t1, t2;

	for(i=0, vtot=0; i<nbonds;)
	{
		type = atoms[i++];
		ai   = atoms[i++];
		aj   = atoms[i++];
		ak   = atoms[i++];
		k    = params[type].harmonic.krA;
		x    = params[type].harmonic.rA * DEG2RAD;
		
		if(k == 0)
			continue;

    theta  = bond_angle(coord[ai], coord[aj], coord[ak], pbc,
                            r_ij, r_kj, &cos_theta, &t1, &t2);  /*  41		*/

		dx   = theta - x;
		dx2  = dx * dx;

		vtot += 0.5 * k * dx2;
	}

	return vtot;
}

static real MulTop_Local_CalcIdihs(int nbonds, const t_iatom atoms[], const t_iparams params[], const rvec coord[], const t_pbc *pbc)
{
	real vtot, k, x;
	int i, type, ai, aj, ak, al;
	real dx, dx2;
	rvec r_ij, r_kj, r_kl, m, n;
	real phi, sign;
	int t1, t2, t3;
	static real ANGLE_WARNING = M_PI * 10;

	for(i=0, vtot=0; i<nbonds;)
	{
		type = atoms[i++];
		ai   = atoms[i++];
		aj   = atoms[i++];
		ak   = atoms[i++];
		al   = atoms[i++];
		k    = params[type].harmonic.krA;
		x    = params[type].harmonic.rA * DEG2RAD;
		
		if(k == 0)
			continue;
    
		phi = dih_angle(coord[ai], coord[aj], coord[ak], coord[al], pbc, 
				r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);  /*  84		*/

		dx   = phi - x;
					
		if(dx >= ANGLE_WARNING || dx < -1 * ANGLE_WARNING)
			fprintf(stderr, "Large calculated angle diff: %lg\n", dx * RAD2DEG);
		while (dx >= M_PI)
      dx -= 2*M_PI;
    while (dx < -M_PI)
      dx += 2*M_PI;

		dx2  = dx * dx;

		vtot += 0.5 * k * dx2;
	}

	return vtot;
}

static real MulTop_Local_CalcPdihs(int nbonds, const t_iatom atoms[], const t_iparams params[], const rvec coord[], const t_pbc *pbc)
{
	real vtot, c0, phi0;
	int i, type, ai, aj, ak, al;
	rvec r_ij, r_kj, r_kl, m, n;
	real phi, sign;
	int t1, t2, t3;
	real multphi;

	for(i=0, vtot=0; i<nbonds;)
	{
		type = atoms[i++];
		ai   = atoms[i++];
		aj   = atoms[i++];
		ak   = atoms[i++];
		al   = atoms[i++];
		c0   = params[type].pdihs.cpA;
		
		if(c0 == 0)
			continue;
		
		phi0 = params[type].pdihs.phiA * DEG2RAD;
    
		phi = dih_angle(coord[ai], coord[aj], coord[ak], coord[al], pbc, 
				r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);  /*  84		*/
		multphi = params[type].pdihs.mult * phi;

		vtot += c0 + c0 * cos(multphi - phi0);
	}

	return vtot;
}

static real MulTop_Local_CalcPIdihs(int nbonds, const t_iatom atoms[], const t_iparams params[], const rvec coord[], const t_pbc *pbc)
{
	return MulTop_Local_CalcPdihs(nbonds, atoms, params, coord, pbc);
}

static real MulTop_Local_CalcRBdihs(int nbonds, const t_iatom atoms[], const t_iparams params[], const rvec coord[], const t_pbc *pbc)
{
	real vtot, c0, phi0;
	int i, j, type, ai, aj, ak, al;
	rvec r_ij, r_kj, r_kl, m, n;
	real phi, sign;
	int t1, t2, t3;
  real param[NR_RBDIHS];
	gmx_bool skip;
	real v, cosfac, cos_phi;

	for(i=0, vtot=0; i<nbonds;)
	{
		type = atoms[i++];
		ai   = atoms[i++];
		aj   = atoms[i++];
		ak   = atoms[i++];
		al   = atoms[i++];

		for(j=0, skip = TRUE; j<NR_RBDIHS; j++)
		{
			param[j] = params[type].rbdihs.rbcA[j];
			
			if(param[j] != 0)
				skip = FALSE;
		}
		if(skip)
			continue;

		phi = dih_angle(coord[ai], coord[aj], coord[ak], coord[al], pbc, 
				r_ij, r_kj, r_kl, m, n, &sign, &t1, &t2, &t3);  /*  84		*/
		cos_phi = cos(phi);

		for(j = 1, v = param[0], cosfac = cos_phi; j<NR_RBDIHS; j++)
		{
			cosfac *= cos_phi;
			v += param[j] * cosfac;
		}

		vtot += v;
	}

	return vtot;
}

static real MulTop_Local_CalcFdihs(int nbonds, const t_iatom atoms[], const t_iparams params[], const rvec coord[], const t_pbc *pbc)
{
	return MulTop_Local_CalcRBdihs(nbonds, atoms, params, coord, pbc);
}

/*static void update_lr_send(gmx_domdec_t *dd, MulTop_LocalRecords *lr)*/
/*{*/
/*int i,j,n;*/
/*int curr, nalloc=0;*/
/*int *buf=NULL;*/
/*int ms_rank = dd->masterrank;*/
/*int my_rank = dd->rank;*/
/*int my_nr = dd->nat_tot;*/
/*int tot_nodes = dd->nnodes;*/

/*if(my_rank == ms_rank)*/
/*{*/
/*for(n=0; n< tot_nodes; n++)*/
/*{*/
/*if (n == my_rank)*/
/*continue;*/

/*if (MA_NAT[n] > nalloc)*/
/*{*/
/*nalloc = over_alloc_dd(MA_NAT[n]);*/
/*srenew(buf,nalloc);*/
/*}*/

/*curr = 0;*/
/*for(i=MA_INDEX[n]; i<MA_INDEX[n+1]; i++)*/
/*{*/
/*for(j=COMM_CGS->index[MA_CG[i]]; j<COMM_CGS->index[MA_CG[i]+1]; j++)*/
/*buf[curr++] = j;*/
/*}*/

/*#ifdef GMX_MPI*/
/*MPI_Send(buf,curr,MPI_INT,n,n,dd->mpi_comm_all);*/
/*#endif*/
/*}*/
/*sfree(buf);*/

/*Last do master node */
/*n = my_rank;*/
/*curr = 0;*/
/*for(i=MA_INDEX[n]; i<MA_INDEX[n+1]; i++)*/
/*{*/
/*for(j=COMM_CGS->index[MA_CG[i]]; j<COMM_CGS->index[MA_CG[i]+1]; j++)*/
/*{*/
/*LREC_G2L[j] = curr;*/
/*LREC_L2G[curr++] = j;*/
/*}*/
/*}*/
/*}*/
/*else*/
/*{*/
/*#ifdef GMX_MPI*/
/*MPI_Recv(LREC_L2G,my_nr,MPI_INT,ms_rank,MPI_ANY_TAG,dd->mpi_comm_all,MPI_STATUS_IGNORE);*/
/*#endif*/
/*for(i=0; i< my_nr; i++)*/
/*LREC_G2L[LREC_L2G[i]] = i;*/
/*}*/
/*}*/

/*static void update_lr_scatter(gmx_domdec_t *dd, MulTop_LocalRecords *lr)*/
/*{*/
/*int i,j,n;*/
/*int curr, nalloc=0;*/
/*int *scounts = NULL, *disps = NULL, *buf=NULL;*/
/*int ms_rank = dd->masterrank;*/
/*int my_rank = dd->rank;*/
/*int my_nr = dd->nat_tot;*/
/*int tot_nodes = dd->nnodes;*/

/*if(my_rank == ms_rank)*/
/*{*/
/*scounts = MA_IBUF;*/
/*disps = scounts + tot_nodes;*/

/*for(n=0; n < tot_nodes; n++)*/
/*{*/
/*scounts[n] = MA_NAT[n];*/
/*disps[n] = n == 0 ? 0: disps[n-1] + scounts[n-1];*/
/*}*/

/*snew(buf,COMM_CGS->index[COMM_CGS->nr]);*/
/*curr = 0;*/
/*for(n=0; n < tot_nodes; n++)*/
/*{*/
/*for(i=MA_INDEX[n]; i<MA_INDEX[n+1]; i++)*/
/*{*/
/*for(j=COMM_CGS->index[MA_CG[i]]; j<COMM_CGS->index[MA_CG[i]+1]; j++)*/
/*buf[curr++] = j;*/
/*}*/
/*}*/
/*}*/

/*#ifdef GMX_MPI*/
/*MPI_Scatterv(buf,scounts,disps,MPI_INT,LREC_L2G,my_nr,MPI_INT,ms_rank,dd->mpi_comm_all);*/
/*#endif*/

/*for(i=0; i < my_nr; i++)*/
/*LREC_G2L[LREC_L2G[i]] = i;*/
/*}*/

/*################################################
 * non-static functions */

int MulTop_Global_GetInputFileName(char **fns[], const char *opt, int nfile, const t_filenm fnm[], t_commrec *cr)
{
  int ntop, i;
	char mysim[2];
  char *fname, *p;

	ntop = opt2fns(fns, opt, nfile, fnm);

	/* replace the  */
	if(cr->ms != NULL)
	{
		if (cr->ms->nsim >= 10)
			gmx_fatal(FARGS,"Dont support number of simulations larger than 9.\n");
		sprintf(mysim,"%d",cr->ms->sim);
	}
	else mysim[0]='\0';

	for(i=0; i<ntop; i++)
	{
		fname = *(fns[i]);
	
		p = strstr(fname, ".tpr.mdp");
	
		if(mysim[0] != '\0')
		{
			strncpy(p, mysim, 1);
			p++;
		}
		strncpy(p, ".tpr",4);
		p[4] = '\0';
	}

	return ntop+1;
}

mt_gtops_t *MulTop_Global_Init(int ntop, real Tref, real Tmax, real Wmax, gmx_bool bMASTER)
{
	int i;
	mt_gtops_t *gtops;

	snew(gtops, 1);

	gtops->ntop = ntop;
	gtops->Tref = Tref;
	gtops->Tmax = Tmax;
	gtops->Wmax = Wmax;

	snew(gtops->mtops, ntop);
	for(i=1; i<ntop; i++)
		snew(gtops->mtops[i], 1);

	snew(gtops->dr, 1);
	MulTop_Global_InitDihedral(gtops);

	gtops->gr=NULL;
	if(bMASTER)
	{
		snew(gtops->gr, 1);
		MulTop_Global_InitRecords(gtops);
	}

	/* Init the static array. The memory will be released after writing the GREC file. */
	snew(booltable, MT_MAX_ATOM_PER_MOL);

	return gtops;
}

void MulTop_Global_SetReferenceTopology(mt_gtops_t *gtops, gmx_mtop_t *mtop)
{
	gtops->mtops[0] = mtop;
}

void MulTop_Global_GetOtherTopologies(mt_gtops_t *gtops, char **files, t_commrec *cr, t_state **states)
{
	int i;
	gmx_mtop_t *mtop;
	t_inputrec *ir;

	for(i=1; i<gtops->ntop; i++)
	{
		mtop = gtops->mtops[i];

		snew(ir, 1);

		read_tpx_state(files[i-1], ir, states[i], NULL, mtop);

		sfree(ir);
	}
}

void MulTop_Global_LoadData(mt_gtops_t *gtops)
{
	int natoms;
	int ntop = gtops->ntop;
	int i,j,k,t,index;
	char header;
	int header_index, header_natoms;
	MulTop_GlobalRecords *gr = gtops->gr;
	FILE *fp;

	if(!MulTop_Global_OpenFile(gtops, &fp))
		gmx_fatal(FARGS,"File type is 'w'. No data to load!");

	for(t=0; t<ntop; t++)
	{
		natoms = gtops->mtops[t]->natoms;
		srenew(GREC_FA[t], natoms);	

		/* Read header information */
		if(t!=0)
			fscanf(fp, "%c", &header);
		fscanf(fp, "%c", &header);
		if(header != '#')
			gmx_fatal(FARGS,"Error reading record files: cannot find locations for topology %d in file.", t);
	
		fscanf(fp, "%d", &header_index);
		if(header_index != t)
			gmx_fatal(FARGS,"Error reading record files: topology index in file is not consistent with the current topology (index = %d).",t);

		fscanf(fp, "%d", &header_natoms);
		if(header_natoms != natoms)
			gmx_fatal(FARGS,"Error reading record files: number of atoms in file is not consistent with one in topology %d which is %d.",t, natoms);

		/* Read GREC_FA information */
		for(i=0; i<natoms; i++)
		{
			if(feof(fp))
				gmx_fatal(FARGS,"Error reading record files: number of atoms in file is not consistent with one in topology %d which is %d.",t, natoms);
		
			fscanf(fp, "%d", &index); 
			*(GREC_FA[t]+i) = index - 1; /* index is the one in *.gro so we need to reduce it by one. */
		}
		
		/* Read GREC_FB information after the last topology is processed */
		if(t == ntop - 1)
		{
			for(i=0; i< ntop; i++)
			for(j=0; j< MT_MAX_BLOCK; j++)
				for(k=0; k< ntop; k++)
				{
					if(feof(fp))
						gmx_fatal(FARGS,"Error reading record files: file content is not consistent with molecule block information in topology %d.",i);

					fscanf(fp, "%d", GREC_FB[i][j]+k);
				}
		}
	}

	fclose(fp);
}

void MulTop_Global_CalcData(mt_gtops_t *gtops, t_state **states_tot, t_state *state_ref)
{
	int ntop = gtops->ntop;
	FILE *fp;
	t_state *state_tar;
	int i;

	if(MulTop_Global_OpenFile(gtops, &fp))
		gmx_fatal(FARGS,"File type is 'r'. No data to be calculated!");

	for(i=0; i<ntop; i++)
	{
		state_tar = states_tot[i]; /* state of the (i)th topology */
	
		MulTop_Global_CalcDataForEachTopology(gtops, i, state_tar, state_ref, fp);
	}

	fclose(fp);
}

void MulTop_Global_RefreshForceFieldParameters(mt_gtops_t *gtops)
{
	int b0, bi;
	int a0, ai;
	int n, ftype;
	t_iatom ***iatoms_buf;
	t_iparams **iparams_buf;
	t_functype *functype_buf;
	int ff_bufsize, ff_actualsize;
	int iatoms_bufsize, iatoms_actualsize;
	int ntop=gtops->ntop;
	gmx_mtop_t *mtop0 = gtops->mtops[0];
	MulTop_GlobalRecords   *gr = gtops->gr;
	MulTop_DihedralRecords *dr = gtops->dr;

/*####################################*/
/* First, make all mtops share the same set of atomtypes.*/

	int atnr = mtop0->ffparams.atnr;
	int atnr2 = atnr * atnr;
	int atnr_n, atnr_n2;

	for(n=1; n<ntop; n++)
	{
		atnr_n = gtops->mtops[n]->ffparams.atnr;
		atnr_n2 = atnr_n * atnr_n; 
			
		if(atnr_n2 == atnr2)
			fprintf(stdout,"# of atomtypes in topology %d is consistent with the reference topology. Will neglect the order of atomtype.\n", n);
		else if(atnr_n2 == 1)
			fprintf(stdout,"# of atomtypes in topology %d is 1. Will add the non-bonded term to every LJ pairs in the reference topology\n", n);
		else
			gmx_fatal(FARGS,"# of atomtypes in topology %d is %d, which is neither 1 nor %d(# of atomtypes in reference topology). Don't how to deal with the non-bonded interation.", n, atnr_n, atnr);
	}

/*####################################*/
/* Estimate how much memory is needed (at maximum), and allocate memory. */

	snew(iatoms_buf, mtop0->nmolblock);
	
	for(b0=0, ff_bufsize=atnr2; b0<mtop0->nmolblock; b0++)
	{
		snew(iatoms_buf[b0], F_NRE);
		
		for(ftype=0, iatoms_bufsize=0; ftype<F_NRE; ftype++)
		{
			if(ftype == F_LJ) /* # of the non-bonded interaction has alrady been counted. (atnr2) */
				continue;

			for(n=0 ; n<ntop; n++)
			{
				bi = t1b1t2(0, b0, n, gtops->gr);
	
				if(bi == -1)
					continue;
			
				ff_bufsize += (gtops->mtops[n]->moltype[bi].ilist[ftype].nr)/
											(interaction_function[ftype].nratoms + 1);
				iatoms_bufsize += gtops->mtops[n]->moltype[bi].ilist[ftype].nr;
			}
	
			if(iatoms_bufsize)
				snew(iatoms_buf[b0][ftype], iatoms_bufsize);
		}
	}
	snew(functype_buf, ff_bufsize);
	
	snew(iparams_buf, ntop);
	for(n=0; n<ntop; n++)
		snew(iparams_buf[n], ff_bufsize);

/*####################################*/
/* Fill in the non-bonded terms (for both functype and iparams). */
	
	int ff_curr; 
	
	for(ff_curr=0; ff_curr<atnr2; ff_curr++)
	{
		functype_buf[ff_curr] = mtop0->ffparams.functype[ff_curr];
			
		iparams_copy(&(iparams_buf[0][ff_curr]),&(mtop0->ffparams.iparams[ff_curr]),F_LJ);

		for(n=1; n<ntop; n++)
		{
			atnr_n = gtops->mtops[n]->ffparams.atnr;
			atnr_n2 = atnr_n * atnr_n; 
		
			if(atnr_n2 == 1)
				iparams_copy(&(iparams_buf[n][ff_curr]), &(gtops->mtops[n]->ffparams.iparams[0]), F_LJ);
			else
				iparams_copy(&(iparams_buf[n][ff_curr]), &(gtops->mtops[n]->ffparams.iparams[ff_curr]), F_LJ);
		}
	}

	for(n=1; n<ntop; n++)
	{
		gtops->mtops[n]->atomtypes.nr = atnr;
		gtops->mtops[n]->ffparams.atnr = atnr;
	}

/*####################################*/
/* Do the bonded terms. Fill in functype, iparams and ilists. */
/* Here F_VDW14 is a new interaction type introduced in MulTop. Unlike F_LJ14 in which both Vdw and Coulomb force will be calculated, F_VDW14 only calculate the Vdw force. */

	t_iatom *iatoms;
	t_ilist *reflist;

	for(ftype=0; ftype<F_NRE; ftype++)
	{
		if(ftype == F_VDW14) 
			continue;

		gmx_bool bCONSTRAINT = interaction_function[ftype].flags & IF_CONSTRAINT;
		int nratoms = interaction_function[ftype].nratoms;
		
		for(b0=0; b0<mtop0->nmolblock; b0++)
		{
			if(iatoms_buf[b0][ftype] == NULL)
				continue;

			int curr = 0; /* current length of the new list(iatoms) */
			int VDW14_curr = 0; /* current length of VDW14 list(iatoms) */
			
			iatoms = iatoms_buf[b0][ftype];
			reflist=&(mtop0->moltype[b0].ilist[ftype]);
				
			/* First copy the reflist to the new list */
			int constr_start_index;
			int constr_end_index;
			int settle_start_index;
			int settle_end_index;
			int ff_curr_record = ff_curr;
			if(reflist->nr > 0)
			{
				constr_start_index = reflist->iatoms[0];
				constr_end_index = reflist->iatoms[0] - 1;
				settle_start_index = reflist->iatoms[0];
				settle_end_index = reflist->iatoms[0] - 1;
			}
				
			for(;	curr<reflist->nr; curr++)
			{
				gmx_bool bAddNewFunc;

				if(!(curr % (nratoms+1)))
				{
					if(ftype == F_CONSTR)
					{
						if(reflist->iatoms[curr] > constr_end_index)
							{ bAddNewFunc = 1; constr_end_index ++;}
						else bAddNewFunc = 0; /* As constraint and settle are not allowed in the additional topologies, we don't need to add repeating interactions here. */
					}
					else if(ftype == F_SETTLE)
					{
						if(reflist->iatoms[curr] > settle_end_index)
							{ bAddNewFunc = 1; settle_end_index ++;}
						else bAddNewFunc = 0;
					}
					else bAddNewFunc = 1;

					if(bAddNewFunc)
					{
						functype_buf[ff_curr] = mtop0->ffparams.functype[reflist->iatoms[curr]];

						iparams_copy(&(iparams_buf[0][ff_curr]), &(mtop0->ffparams.iparams[reflist->iatoms[curr]]), ftype);

						iatoms[curr] = ff_curr;

						ff_curr++;
					}
					else if(ftype == F_CONSTR)
						iatoms[curr] = ff_curr_record + reflist->iatoms[curr] - constr_start_index;
					else  /* if(ftype == F_SETTLE) */
						iatoms[curr] = ff_curr_record + reflist->iatoms[curr] - settle_start_index;
				}
				else
					iatoms[curr] = reflist->iatoms[curr];
			}

			/* Next fill in the new ilist with datas in mtops[i] */
			t_ilist *mylist;
			int mylist_cur;  /* cursor in mylist(mtops[n]) */
			int newlist_cur; /* cursor in new lists */
			gmx_bool bFound;
		
			for(n=1; n<ntop; n++)
			{
				gmx_mtop_t *mytop = gtops->mtops[n];
				bi = t1b1t2(0, b0, n, gtops->gr);

				if(bi == -1)
					continue;

				mylist=&(mytop->moltype[bi].ilist[ftype]);

				for(mylist_cur=0; mylist_cur < mylist->nr; mylist_cur+=nratoms+1)
				{
					bFound = 0;

					for(newlist_cur=0; newlist_cur < curr && !bFound ; newlist_cur+=nratoms+1)
					{
						if(judge_samebonds(&(iatoms[newlist_cur+1]), 
									             &(mylist->iatoms[mylist_cur+1]), 
															 gr, n, nratoms))
						{
							if(bCONSTRAINT)
								gmx_fatal(FARGS,"More than one constraints are included to the same bond, which is forbidden. Check your topologies.\n");

							/* Special case: even if they represents the same bonds, we still need to separate them because of different multiplicities.*/
							if(ftype == F_PDIHS || ftype == F_PIDIHS)
							{
								if (mytop->ffparams.iparams[mylist->iatoms[mylist_cur]].pdihs.mult != iparams_buf[0][iatoms[newlist_cur]].pdihs.mult )
										continue;
							}

							iparams_copy(&(iparams_buf[n][iatoms[newlist_cur]]), 
													 &(mytop->ffparams.iparams[mylist->iatoms[mylist_cur]]), ftype);
							bFound = 1;
						}
					}

					if(!bFound) /* need to increase the new list length: new bond founded. also increase ff_curr */
					{
						
						/* Special for MulTop: add VDW14 interaction (see nonbonded.c). 
						 * For non-reference topologies, LJ14 pairs CANNOT have Coulomb interaction, which is the reason we define VDW14 here.*/
						if(ftype != F_LJ14)
						{
							iatoms[curr] = ff_curr;
							fillin_bonds(&(iatoms[curr+1]), 
												 &(mylist->iatoms[mylist_cur+1]),nratoms,gr,n);
							functype_buf[ff_curr] = mytop->ffparams.functype[mylist->iatoms[mylist_cur]];
							iparams_copy(&(iparams_buf[n][ff_curr]),
													 &(mytop->ffparams.iparams[mylist->iatoms[mylist_cur]]), ftype);
							curr += nratoms+1;
						}
						else
						{
							srenew(iatoms_buf[b0][F_VDW14],VDW14_curr+nratoms+1);
							iatoms_buf[b0][F_VDW14][VDW14_curr] = ff_curr;
							fillin_bonds(&(iatoms_buf[b0][F_VDW14][VDW14_curr+1]), 
												 &(mylist->iatoms[mylist_cur+1]),nratoms,gr,n);
							functype_buf[ff_curr] = F_VDW14;
							iparams_copy(&(iparams_buf[n][ff_curr]),
													 &(mytop->ffparams.iparams[mylist->iatoms[mylist_cur]]), F_VDW14);
							VDW14_curr += nratoms+1;
						}

						ff_curr++;	
					}
				}
			}

			iatoms_actualsize = curr;

			for(n=0; n<ntop; n++)
			{
				bi = t1b1t2(0, b0, n, gtops->gr);
				if(bi == -1)
					continue;
				
				gtops->mtops[n]->moltype[bi].ilist[ftype].nr = iatoms_actualsize;
				sfree(gtops->mtops[n]->moltype[bi].ilist[ftype].iatoms);
				gtops->mtops[n]->moltype[bi].ilist[ftype].iatoms = iatoms_buf[b0][ftype];
				
				if(ftype == F_LJ14)
				{
					gtops->mtops[n]->moltype[bi].ilist[F_VDW14].nr = VDW14_curr;
					sfree(gtops->mtops[n]->moltype[bi].ilist[F_VDW14].iatoms);
					gtops->mtops[n]->moltype[bi].ilist[F_VDW14].iatoms = iatoms_buf[b0][F_VDW14];
				}
			}
		}
	}

	ff_actualsize = ff_curr;
	
	for(n=0; n<ntop; n++)
	{
		sfree(gtops->mtops[n]->ffparams.functype);
		sfree(gtops->mtops[n]->ffparams.iparams);
		gtops->mtops[n]->ffparams.ntypes = ff_actualsize;
		gtops->mtops[n]->ffparams.functype  = functype_buf;
		gtops->mtops[n]->ffparams.iparams = iparams_buf[n];
	}

/*####################################*/
/* The last thing to do is filling in the pdih/pidih cos and sin datas(static) */

	if(DREC_Pind !=0 )
		DREC_Pind = 0;

	while(mtop0->ffparams.functype[DREC_Pind] != F_PDIHS && DREC_Pind<mtop0->ffparams.ntypes)
		DREC_Pind += 1;
	
	int i = DREC_Pind;
	while(mtop0->ffparams.functype[i] == F_PDIHS)
	{
		for(n=0; n<ntop; n++)
		{
			DREC_Pcos[i-DREC_Pind][n] = gtops->mtops[n]->ffparams.iparams[i].pdihs.cpA
					                             * cos(DEG2RAD * gtops->mtops[n]->ffparams.iparams[i].pdihs.phiA);
			DREC_Psin[i-DREC_Pind][n] = gtops->mtops[n]->ffparams.iparams[i].pdihs.cpA
					                             * sin(DEG2RAD * gtops->mtops[n]->ffparams.iparams[i].pdihs.phiA);
		}

		i++;
	}

	/* just for debug use */
	/*PDIH_manifest(gtops);*/

	/* for future use. Currently we don't have pidih interactions. */
	if(DREC_PIind !=0 )
		DREC_PIind = 0;

	while(mtop0->ffparams.functype[DREC_PIind] != F_PIDIHS && DREC_PIind<mtop0->ffparams.ntypes)
		DREC_PIind += 1;
	
	i = DREC_PIind;
	while(mtop0->ffparams.functype[i] == F_PIDIHS)
	{
		for(n=0; n<ntop; n++)
		{
			DREC_PIcos[i-DREC_PIind][n] = gtops->mtops[n]->ffparams.iparams[i].pdihs.cpA
					                             * cos(DEG2RAD * gtops->mtops[n]->ffparams.iparams[i].pdihs.phiA);
			DREC_PIsin[i-DREC_PIind][n] = gtops->mtops[n]->ffparams.iparams[i].pdihs.cpA
					                             * sin(DEG2RAD * gtops->mtops[n]->ffparams.iparams[i].pdihs.phiA);
		}

		i++;
	}

	/* just for debug use */
	/*PIDIH_manifest(gtops);*/
	
	/* Because gr has already been saved and will NOT be used in the rest of the current simulation, we will clean the records (free and memory). */
	MulTop_Global_CleanRecords(gtops->gr, ntop);
}

void MulTop_Global_Bcast(mt_gtops_t *gtops, t_commrec *cr)
{
	int i;
	int ntop = gtops->ntop;
	t_inputrec *ir_dummy;

	/* broadcast topologies */
	for(i=0; i<ntop; i++)
	{
		snew(ir_dummy, 1);

		bcast_ir_mtop(cr, ir_dummy, gtops->mtops[i]);
		
		sfree(ir_dummy);
	}
	
	/* broadcast dihedral records */
	MulTop_DihedralRecords *dr = gtops->dr;

	gmx_bcast(sizeof(int), &DREC_Pind, cr);
	gmx_bcast(sizeof(int), &DREC_PIind, cr);
	for(i=0; i<MT_MAX_PDIHS; i++)
	{
		gmx_bcast(ntop*sizeof(real), DREC_Pcos[i], cr);
		gmx_bcast(ntop*sizeof(real), DREC_Psin[i], cr);
		gmx_bcast(ntop*sizeof(real), DREC_PIcos[i], cr);
		gmx_bcast(ntop*sizeof(real), DREC_PIsin[i], cr);
	}
}

real MulTop_Global_GetEnergy(real v, t_commrec *cr)
{
	gmx_sum(1, &v, cr);

	return v;
}

mt_ltops_t *MulTop_Local_Init(mt_gtops_t *gtops)
{
	int ntop = gtops->ntop;
	int i;
	mt_ltops_t *ltops;

	snew(ltops, 1);

	ltops->ntop = ntop;
	snew(ltops->top_final, 1);
	snew(ltops->tops, ntop);

	ltops->Tref = gtops->Tref;
	ltops->Tmax = gtops->Tmax;
	ltops->Wmax = gtops->Wmax;
	
	snew(ltops->lr, 1);
	MulTop_Local_InitRecords(ltops, gtops);
	
	ltops->dr = gtops->dr;

	return ltops;
}

void MulTop_Local_SetReferenceTopology(mt_ltops_t *ltops, gmx_localtop_t *top)
{
	ltops->tops[0] = top;
}
 
void MulTop_Local_GetOtherTopologies(mt_ltops_t *ltops, mt_gtops_t *gtops, t_commrec *cr, t_inputrec *ir)
{
	int i;
	int ntop = ltops->ntop;

	for(i=1; i<ntop; i++)
	{
		if(DOMAINDECOMP(cr))
			ltops->tops[i] = dd_init_local_top(gtops->mtops[i]);
		else if(PAR(cr))
			gmx_fatal(FARGS, "MulTop does not support multithread scheme currently!\n");
		else
			ltops->tops[i] = gmx_mtop_generate_local_top(gtops->mtops[i], ir);
	}
}

void MulTop_Local_UpdateRecords(mt_ltops_t *ltops, gmx_domdec_t *dd)
{
	MulTop_LocalRecords *lr = ltops->lr;
	int i, i_gl, i_glmax=0;
	
	if(dd == NULL)
	{
		memset(LREC_L2G, -1, ltops->tops[0]->excls.nr);
		memset(LREC_G2L, -1, ltops->tops[0]->excls.nr);
		
		MulTop_Local_UpdateRecordsSingleProcessor(lr, ltops->tops[0]->excls.nr);
	}
	else
	{
		memset(LREC_L2G, -1, dd->nat_tot);

		for(i=0; i<dd->nat_tot; i++)
		{
			i_gl = dd->gatindex[i];
			LREC_L2G[i] = i_gl;
			i_glmax = i_gl > i_glmax ? i_gl : i_glmax;
		}
		
		memset(LREC_G2L, -1, i_glmax);

		for(i=0; i<dd->nat_tot; i++)
			LREC_G2L[LREC_L2G[i]] = i;

		/* old method */
		/*if(dd->nnodes <= 4)*/
		/*update_lr_send(dd, lr);*/
		/*else*/
		/*update_lr_scatter(dd, lr);*/
	}
}

void MulTop_Local_InitFinalTopology(mt_ltops_t *ltops)
{
	int i;

	gmx_localtop_t *top_final = ltops->top_final;
	gmx_localtop_t **tops = ltops->tops;

	top_final->idef.ntypes = tops[0]->idef.ntypes;
	top_final->idef.atnr = tops[0]->idef.atnr;
	top_final->idef.fudgeQQ = tops[0]->idef.fudgeQQ;
	top_final->idef.iparams_posres_nalloc = tops[0]->idef.iparams_posres_nalloc;
	top_final->idef.ilsort = tops[0]->idef.ilsort;
	top_final->idef.cmap_grid.ngrid = tops[0]->idef.cmap_grid.ngrid;
	top_final->idef.cmap_grid.grid_spacing = tops[0]->idef.cmap_grid.ngrid;
	top_final->idef.cmap_grid.cmapdata = tops[0]->idef.cmap_grid.cmapdata;

	snew(top_final->idef.functype, top_final->idef.ntypes);
	for(i=0; i<top_final->idef.ntypes; i++)
		top_final->idef.functype[i] = tops[0]->idef.functype[i];
	snew(top_final->idef.iparams, top_final->idef.ntypes);
	snew(top_final->idef.iparams_posres, top_final->idef.ntypes);
	top_final->atomtypes.nr = tops[0]->atomtypes.nr;
	top_final->atomtypes.radius = tops[0]->atomtypes.radius;
	top_final->atomtypes.vol = tops[0]->atomtypes.vol;
	top_final->atomtypes.surftens = tops[0]->atomtypes.surftens;
	top_final->atomtypes.gb_radius = tops[0]->atomtypes.gb_radius;
	top_final->atomtypes.S_hct = tops[0]->atomtypes.S_hct;
	top_final->atomtypes.atomnumber = tops[0]->atomtypes.atomnumber;
}

void MulTop_Local_UpdateFinalTopologyBasic(mt_ltops_t *ltops)
{	
	int i;
	gmx_localtop_t *top_final = ltops->top_final;
	gmx_localtop_t **tops = ltops->tops;

	for(i=0; i<F_NRE; i++)
	{
		top_final->idef.il[i].nr = tops[0]->idef.il[i].nr;
		top_final->idef.il[i].nr_nonperturbed= tops[0]->idef.il[i].nr_nonperturbed;
		top_final->idef.il[i].nalloc = tops[0]->idef.il[i].nalloc;
		top_final->idef.il[i].iatoms = tops[0]->idef.il[i].iatoms;
	}
	top_final->cgs.nr = tops[0]->cgs.nr;
	top_final->cgs.nalloc_index = tops[0]->cgs.nalloc_index;
	top_final->cgs.index = tops[0]->cgs.index;
	top_final->excls.nr = tops[0]->excls.nr;
	top_final->excls.nra = tops[0]->excls.nra;
	top_final->excls.nalloc_a = tops[0]->excls.nalloc_a;
	top_final->excls.nalloc_index = tops[0]->excls.nalloc_index;
	top_final->excls.index = tops[0]->excls.index;
	top_final->excls.a = tops[0]->excls.a;
}
	
/* Currently this function only support two topologies... */
void MulTop_Local_UpdateFinalTopologyParameters(mt_ltops_t *ltops, real *pot, real T)
{
	gmx_localtop_t *top_final = ltops->top_final;
	gmx_localtop_t **tops = ltops->tops;
	int i, j, k, ftype;
	real weight[2];
	static int count = 0;
	static real parameter;
	
	if(!count)
	{
		if(ltops->Tmax != ltops->Tref)
			parameter = ltops->Wmax / (ltops->Tmax - ltops->Tref);
		else parameter = ltops->Wmax;
	}

	for(i=0; i<F_EPOT; i++)
		pot[i] = 0;

	weight[0] = 1;
	if(T < ltops->Tref)
		weight[1] = 0;
	else if(ltops->Tmax != ltops->Tref)
		weight[1] = parameter * (T - ltops->Tref);
	else weight[1] = parameter;
	if(weight[1] < 1e-6)
		weight[1] = 0;

	for(j=0; j<tops[0]->idef.ntypes; j++)
	{
		ftype = top_final->idef.functype[j];

		iparams_cal(&(top_final->idef.iparams[j]), &(tops[0]->idef.iparams[j]), &(tops[1]->idef.iparams[j]), ftype, weight, j, pot, ltops->dr);
	}

	count ++;
}

gmx_localtop_t *MulTop_Local_GetFinalTopology(mt_ltops_t *ltops)
{
	return ltops->top_final;
}

real MulTop_Local_OnlyCalcAdditionalEnergy(mt_ltops_t *ltops, const t_forcerec *fr, const t_state *state, const gmx_enerdata_t *enerd, real T)
{
	real v, vi;
	real *pot;
	int thread;
	int ntop = ltops->ntop;
	t_ilist *ilist;
	t_iatom *iatoms;
	t_iparams *iparams;
	int i;
	t_pbc *pbc;
	rvec *coord = state->x;
	real weight;
	static int count = 0;
	static real parameter;
	
	if(!count)
	{
		if(ltops->Tmax != ltops->Tref)
			parameter = ltops->Wmax / (ltops->Tmax - ltops->Tref);
		else parameter = ltops->Wmax;
	}
	if(T < ltops->Tref)
		weight = 0;
	else if(ltops->Tmax != ltops->Tref)
		weight = parameter * (T - ltops->Tref);
	else weight = parameter;
	if(weight < 1e-6)
		weight = 0;
	
	if(fr->bMolPBC)
	{
		snew(pbc, 1);
		set_pbc(pbc, fr->ePBC, state->box);
	}
	else pbc = NULL;
	  
	ilist = &(ltops->top_final->idef.il);

	for(i=1, v=0; i<ntop; i++)
	{
		iparams = ltops->tops[i]->idef.iparams;
		
		snew(pot, fr->nthreads);
		
		vi = 0;

#pragma omp parallel for num_threads(fr->nthreads) schedule(static)
    for (thread = 0; thread < fr->nthreads; thread++)
    {
      int  ftype, nbonds, nat1;
      int  nb0, nbn;
	  	real *ener;

	  	ener = pot + thread;

      /* Loop over all bonded force types to calculate the bonded forces */
      for (ftype = 0; (ftype < F_NRE); ftype++)
      {
				if(ilist[ftype].nr == 0)
					continue;

        nat1   = interaction_function[ftype].nratoms + 1;
        nbonds = ilist[ftype].nr/nat1;
				iatoms = ilist[ftype].iatoms;
        nb0 = ((nbonds* thread   )/(fr->nthreads))*nat1;
        nbn = ((nbonds*(thread+1))/(fr->nthreads))*nat1 - nb0;
          
				if (ftype == F_BONDS)
          *ener += MulTop_Local_CalcBonds(nbn, iatoms+nb0, iparams, coord, pbc);
				else if(ftype == F_ANGLES)
          *ener += MulTop_Local_CalcAngles(nbn, iatoms+nb0, iparams, coord, pbc);
				else if(ftype == F_IDIHS)
          *ener += MulTop_Local_CalcIdihs(nbn, iatoms+nb0, iparams, coord, pbc);
				else if(ftype == F_PDIHS)
          *ener += MulTop_Local_CalcPdihs(nbn, iatoms+nb0, iparams, coord, pbc);
				else if(ftype == F_PIDIHS)
          *ener += MulTop_Local_CalcPIdihs(nbn, iatoms+nb0, iparams, coord, pbc);
				else if(ftype == F_RBDIHS)
          *ener += MulTop_Local_CalcRBdihs(nbn, iatoms+nb0, iparams, coord, pbc);
				else if(ftype == F_FOURDIHS)
          *ener += MulTop_Local_CalcFdihs(nbn, iatoms+nb0, iparams, coord, pbc);
				else continue;
      }
    }

    for (thread = 0; thread < fr->nthreads; thread++)
			vi += pot[thread];

		v += vi * weight;
	}

	v += enerd->term[F_VDW14];

	if(pbc != NULL)
		sfree(pbc);

	count++;

	return v;
}
