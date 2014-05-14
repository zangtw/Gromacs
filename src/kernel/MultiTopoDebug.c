#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "gmx_fatal.h"
#include "txtdump.h"

#include "MultiTopoDebug.h"

void mt_debug_p_vec_before_update(rvec *x, rvec *f, t_commrec *cr, t_mdatoms *md, 
		int step, int totstep, gmx_bool bMT)
{
	static int count = 0;
	static int warn_once = 0;
	static char my_name[]="mt_debug_p_vec_before_update";

	if(count>=10)
	{
		if(warn_once == 0)
		{
			fprintf(stdout, "function %s called more than 10 times!\n", my_name);
			warn_once = 1;
		}
		return;
	}

	if(step<totstep)
		return;

	FILE *fp_x, *fp_f;
	char fx[30] = "./debug/x_B_";
	char ff[30] = "./debug/f_B_";
	char nodeid[3];
	char countid[2];
	int i;

	sprintf(nodeid,"%d",cr->nodeid);
	strcat(fx,nodeid);
	strcat(ff,nodeid);
	
	sprintf(countid,"%d",count);
	strcat(fx,countid);
	strcat(ff,countid);

	if((fp_x=fopen(fx,"w"))==NULL)
		gmx_fatal(FARGS,"cannot open file %s\n",fx);
	if((fp_f=fopen(ff,"w"))==NULL)
		gmx_fatal(FARGS,"cannot open file %s\n",ff);

	if(!DOMAINDECOMP(cr))
	{
		for(i=0; i< md->nr; i++)
		{
			fprintf(fp_x, "x[%d]:\t%f\t%f\t%f\tgindex%d\n",i,x[i][0],x[i][1],x[i][2], i);
			fprintf(fp_f, "f[%d]:\t%f\t%f\t%f\tgindex%d\n",i,f[i][0],f[i][1],f[i][2], i);
		}
	}
	else
	{
		fprintf(fp_x,"nat_home=%d,nat_tot=%d\n",cr->dd->nat_home,cr->dd->nat_tot);
		fprintf(fp_f,"nat_home=%d,nat_tot=%d\n",cr->dd->nat_home,cr->dd->nat_tot);
		for(i=0; i< cr->dd->nat_tot; i++)
		{
			fprintf(fp_x, "x[%d]:\t%f\t%f\t%f\tgindex%d\n",i,x[i][0],x[i][1],x[i][2], cr->dd->gatindex[i]);
			fprintf(fp_f, "f[%d]:\t%f\t%f\t%f\tgindex%d\n",i,f[i][0],f[i][1],f[i][2], cr->dd->gatindex[i]);
		}
	}

	count++;

	fclose(fp_x);
	fclose(fp_f);
}

void mt_debug_p_vec_after_update(rvec *x, rvec *f, t_commrec *cr, t_mdatoms *md, 
		int step, int totstep, gmx_bool bMT, gmx_bool TERMINATE)
{
	static int count = 0;
	static int warn_once = 0;
	static char my_name[]="mt_debug_p_vec_after_update";

	if(TERMINATE)
		gmx_fatal(FARGS,"debug finished.\n");
	
	if(count>=10)
	{
		if(warn_once == 0)
		{
			fprintf(stdout, "function %s called more than 10 times!\n", my_name);
			warn_once = 1;
		}
		return;
	}

	if(step<totstep)
		return;

	FILE *fp_x, *fp_f;
	char fx[30] = "./debug/x_A_";
	char ff[30] = "./debug/f_A_";
	char nodeid[3];
	char countid[2];
	int i;

	sprintf(nodeid,"%d",cr->nodeid);
	strcat(fx,nodeid);
	strcat(ff,nodeid);
	
	sprintf(countid,"%d",count);
	strcat(fx,countid);
	strcat(ff,countid);

	if((fp_x=fopen(fx,"w"))==NULL)
		gmx_fatal(FARGS,"cannot open file %s\n",fx);
	if((fp_f=fopen(ff,"w"))==NULL)
		gmx_fatal(FARGS,"cannot open file %s\n",ff);

	if(!DOMAINDECOMP(cr))
	{
		for(i=0; i< md->nr; i++)
		{
			fprintf(fp_x, "x[%d]:\t%f\t%f\t%f\tgindex%d\n",i,x[i][0],x[i][1],x[i][2], i);
			fprintf(fp_f, "f[%d]:\t%f\t%f\t%f\tgindex%d\n",i,f[i][0],f[i][1],f[i][2], i);
		}
	}
	else
	{
		fprintf(fp_x,"nat_home=%d,nat_tot=%d\n",cr->dd->nat_home,cr->dd->nat_tot);
		fprintf(fp_f,"nat_home=%d,nat_tot=%d\n",cr->dd->nat_home,cr->dd->nat_tot);
		for(i=0; i< cr->dd->nat_tot; i++)
		{
			fprintf(fp_x, "x[%d]:\t%f\t%f\t%f\tgindex%d\n",i,x[i][0],x[i][1],x[i][2], cr->dd->gatindex[i]);
			fprintf(fp_f, "f[%d]:\t%f\t%f\t%f\tgindex%d\n",i,f[i][0],f[i][1],f[i][2], cr->dd->gatindex[i]);
		}
	}

	count++;

	fclose(fp_x);
	fclose(fp_f);
}

void mt_debug_p_top(gmx_localtop_t *top, t_blocka *excls, t_commrec *cr, int ltopid)
{
	static int count = 0;
	static int warn_once = 0;
	static char my_name[]="mt_debug_p_top";

	if(count>=10)
	{
		if(warn_once == 0)
		{
			fprintf(stdout, "function %s called more than 10 times!\n", my_name);
			warn_once = 1;
		}
		return;
	}

	FILE *fp;
	char filenm[30] = "./debug/top_";
	char nodeid[3];
	char countid[2];
	char topid[2];

	if(cr!=NULL)
	{
		sprintf(nodeid,"%d",cr->nodeid);
		strcat(filenm,nodeid);
	}
	
	sprintf(countid,"%d",count);
	strcat(filenm,countid);
	
	sprintf(topid,"%d",ltopid);
	strcat(filenm,topid);

	if((fp=fopen(filenm,"w"))==NULL)
		gmx_fatal(FARGS,"cannot open file %s\n",filenm);
	
  pr_title(fp,0,"topology");
	if(excls != NULL)
		pr_blocka(fp,3,"excls",excls,1);
  pr_idef(fp,6,"idef",&top->idef,1);
	
	count++;

	fclose(fp);
}

void mt_debug_p_mtop(gmx_mtop_t *mtop, t_commrec *cr, int mtopid)
{
	static int count = 0;
	static int warn_once = 0;
	static char my_name[]="mt_debug_p_mtop";

	if(count>=10)
	{
		if(warn_once == 0)
		{
			fprintf(stdout, "function %s called more than 10 times!\n", my_name);
			warn_once = 1;
		}
		return;
	}

	FILE *fp;
	char filenm[30] = "./debug/mtop_";
	char nodeid[3];
	char countid[2];
	char topid[2];

	if(cr!=NULL)
	{
		sprintf(nodeid,"%d",cr->nodeid);
		strcat(filenm,nodeid);
	}
	
	sprintf(countid,"%d",count);
	strcat(filenm,countid);
	
	sprintf(topid,"%d",mtopid);
	strcat(filenm,topid);

	if((fp=fopen(filenm,"w"))==NULL)
		gmx_fatal(FARGS,"cannot open file %s\n",filenm);
	
	pr_mtop(fp,0,"topology",mtop,1);
	
	count++;

	fclose(fp);
}

void mt_debug_p_force(rvec *f, rvec *f_pme, t_commrec *cr, char *label, t_mdatoms *md, gmx_bool bMT)
{
	static int count = 0;
	static int warn_once = 0;
	static char my_name[]="mt_debug_p_force";

	if(count>=10)
	{
		if(warn_once == 0)
		{
			fprintf(stdout, "function %s called more than 10 times!\n", my_name);
			warn_once = 1;
		}
		return;
	}
	
	FILE *fp;
	char filenm[30] = "./debug/force_A_";
	char nodeid[3];
	char countid[2];
	int i;

	strcat(filenm,label);
	strcat(filenm,"_");
	
	sprintf(nodeid,"%d",cr->nodeid);
	strcat(filenm,nodeid);
	
	sprintf(countid,"%d",count);
	strcat(filenm,countid);
	
	if((fp=fopen(filenm,"w"))==NULL)
		gmx_fatal(FARGS,"cannot open file %s\n",filenm);
	
	if(f_pme == NULL)
	{
		for(i=0; i<md->nalloc; i++)
		{
			fprintf(fp,"f[%d]:\t",i);
			if(!bMT)
				fprintf(fp,"%f\t%f\t%f\n",f[i][0],f[i][1],f[i][2]);
			else
				fprintf(fp,"%f\t%f\t%f\n",f[i][0]+f[i+md->nalloc][0]+f[i+2*md->nalloc][0],
						f[i][1]+f[i+md->nalloc][1]+f[i+2*md->nalloc][1],
						f[i][2]+f[i+md->nalloc][2]+f[i+2*md->nalloc][2]);
		}
	}
	else
	{
		for(i=0; i<md->nalloc; i++)
		{
			fprintf(fp,"f[%d]:\t",i);
			if(!bMT)
				fprintf(fp,"%f\t%f\t%f\n",f[i][0]+f_pme[i][0],f[i][1]+f_pme[i][1],f[i][2]+f_pme[i][2]);
			else
				fprintf(fp,"%f\t%f\t%f\n",f[i][0]+f[i+md->nalloc][0]+f[i+2*md->nalloc][0]+f_pme[i][0],
						f[i][1]+f[i+md->nalloc][1]+f[i+2*md->nalloc][1]+f_pme[i][1],
						f[i][2]+f[i+md->nalloc][2]+f[i+2*md->nalloc][2]+f_pme[i][2]);
		}
	}

	count++;

	fclose(fp);
}

void mt_debug_p_x(t_commrec *cr, rvec *x, int nr, gmx_bool bMT)
{
	static int count = 0;
	static int warn_once = 0;
	static char my_name[]="mt_debug_p_state";

	if(count>=10)
	{
		if(warn_once == 0)
		{
			fprintf(stdout, "function %s called more than 10 times!\n", my_name);
			warn_once = 1;
		}
		return;
	}

	FILE *fp;
	char filenm[30] = "./debug/state_";
	char countid[2];
	char nodeid[3];
	int i;

	sprintf(nodeid,"%d",cr->nodeid);
	strcat(filenm,nodeid);

	sprintf(countid,"%d",count);
	strcat(filenm,countid);

	if((fp=fopen(filenm,"w"))==NULL)
		gmx_fatal(FARGS,"cannot open file %s\n",filenm);

	for(i=0; i<nr+(bMT)*2*nr; i++)
		fprintf(fp,"x[%d]={%f,\t%f,\t%f}\n",i, x[i][0], x[i][1], x[i][2]);

	count++;

	fclose(fp);
}

void mt_debug_p_mdatoms(t_commrec *cr, t_mdatoms *md, gmx_bool bMT)
{
	static int count = 0;
	static int warn_once = 0;
	static char my_name[]="mt_debug_p_mdatoms";

	if(count>=10)
	{
		if(warn_once == 0)
		{
			fprintf(stdout, "function %s called more than 10 times!\n", my_name);
			warn_once = 1;
		}
		return;
	}

	FILE *fp;
	char filenm[30] = "./debug/mdatoms_";
	char countid[2];
	char nodeid[3];
	int i;

	sprintf(nodeid,"%d",cr->nodeid);
	strcat(filenm,nodeid);

	sprintf(countid,"%d",count);
	strcat(filenm,countid);

	if((fp=fopen(filenm,"w"))==NULL)
		gmx_fatal(FARGS,"cannot open file %s\n",filenm);

	fprintf(fp,"nr=%d, nalloc=%d, start=%d, homenr=%d\n",md->nr,md->nalloc,md->start,md->homenr);
	if(DOMAINDECOMP(cr))
		fprintf(fp,"nat_home=%d,nat_tot=%d\n",cr->dd->nat_home,cr->dd->nat_tot);
	
	for(i=0; i<md->nalloc+(bMT)*2*md->nalloc; i++)
	{
		if(cr->dd != NULL)
		{
			fprintf(fp,"charge[%d]=%f\t\t\ttype[%d]=%d\t\t\t",i,(md->chargeA)[i],i,(md->typeA)[i]);
			if(i<md->nalloc)
				fprintf(fp,"global_ind[%d]=%d\n",i,cr->dd->gatindex[i]);
			else fprintf(fp,"\n");
		}
		else
			fprintf(fp,"charge[%d]=%f\t\t\ttype[%d]=%d\t\t\tglobal_ind[%d]=%d\n",i,(md->chargeA)[i],i,(md->typeA)[i],i,i);
	}

	count++;

	fclose(fp);
}
