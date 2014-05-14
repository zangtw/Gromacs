#ifndef MT_DEBUG_H__
#define MT_DEBUG_H__

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "MultiTopo.h"

void mt_debug_p_vec_before_update(rvec *x, rvec *f, t_commrec *cr, t_mdatoms *md, 
		int step, int totstep, gmx_bool bVirtual);

void mt_debug_p_vec_after_update(rvec *x, rvec *f, t_commrec *cr, t_mdatoms *md, 
		int step, int totstep, gmx_bool bVirtual, gmx_bool TERMINATE);

void mt_debug_p_top(gmx_localtop_t *top, t_blocka *excls, t_commrec *cr, int ltopid);

void mt_debug_p_mtop(gmx_mtop_t *mtop, t_commrec *cr, int mtopid);

void mt_debug_p_force(rvec *f, rvec *f_pme, t_commrec *cr, char *label, t_mdatoms *md, gmx_bool bVirtual);

void mt_debug_p_x(t_commrec *cr, rvec *x, int nr, gmx_bool bVirtual);

void mt_debug_p_mdatoms(t_commrec *cr, t_mdatoms *md, gmx_bool bVirtual);

#endif
