/*!----------------------------------------------------------------------
\file
\brief setting free surface conditions for fluid and ale

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID
*//*! @{ (documentation module open)*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "fluid_prototypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid2/fluid2_prototypes.h"
#include "../fluid3/fluid3.h"
#include "../fluid3/fluid3_prototypes.h"
INT cmp_int(const void *a, const void *b );
void fluid_calelm_hf(	      
                      PARTITION           *actpart,
                      INTRA               *actintra,  
                      CALC_ACTION         *action,
                      CONTAINER           *container,
                      OLL                 *hfmat_oll,
                      DIST_VECTOR         *rhs,
                      INT                  init
                     );
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 7/01    |
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;    /* element stiffness matrix       */
extern struct _ARRAY etforce_global;  /* element Time RHS               */
extern struct _ARRAY eiforce_global;  /* element Iter RHS               */

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;  

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h                                                  
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
enum _CALC_ACTION calc_action[MAXFIELD];

/*!--------------------------------------------------------------------- 
\brief seperat treatment of heightfunction

<pre>                                                         genk 11/03

In this function the equation for the heightfunction = position of the
free surface is solved seperately from the fluid problem. This is done
during the nonlinear iteration.

</pre>
\param    hctrl       INT                control flag
\param   *grat        DOUBLE             convergence ration
\param   *fdyn        FLUID_DYNAMIC 
\param   *actfield    FIELD               actual field
\param   *actpart     PARTITION           actual partition
\param   *actintra    INTRA               actual intra-com.
\param   *action      CALC_ACTION
\param   *container   CONTAINER
\param    converged   INT                 convergence flag

\return void                                            

------------------------------------------------------------------------*/
void fluid_heightfunc(INT                  hctrl,
                      DOUBLE              *grat, 
		      FIELD               *actfield,
                      PARTITION           *actpart,
                      INTRA               *actintra,
                      CALC_ACTION         *action,
                      CONTAINER           *container,
                      INT                  converged)
{
#ifdef D_FSI
static INT             numnp_total;
static INT             numele_total;
static INT             numeq_total;
static INT             numeq;
static INT             numdf;
static INT            *update;
INT                    counter;
INT                    i;
INT                    init;
INT                    dof;
static INT             phipos;
static INT             myrank;
static SOLVAR         *hfsolv; 
static OLL            *hfmat_oll;  
static DIST_VECTOR    *rhs;
static DIST_VECTOR    *sol;
static DISCRET        *actdis;
NODE                  *actnode;
static ARRAY           result_a;
static DOUBLE         *result;
static ARRAY           result_old_a;
static DOUBLE         *result_old;
DOUBLE                 gnorm=ZERO;
DOUBLE                 dgnorm=ZERO;
static FLUID_DYNAMIC  *fdyn;

#ifdef DEBUG 
dstrc_enter("fluid_heightfunc");
#endif


switch (hctrl)
{
/*======================================================================*
 |                      I N I T I A L I S A T I O N                     |
 *======================================================================*/
case 1: 
/*----------------------------------------------------------------------*/
fdyn   = alldyn[genprob.numff].fdyn;
fdyn->hf_stab=2;
numnp_total  = actfield->dis[0].numnp;
numele_total = actfield->dis[0].numele;
numeq_total  = fdyn->hf_numeq_total;
numeq        = fdyn->hf_numeq;
numdf        = fdyn->numdf;
phipos       = numdf-2;
myrank       = actintra->intra_rank;
actdis       = &(actfield->dis[0]);
result       = amdef("result",&result_a,numeq_total,1,"DV");
result_old   = amdef("result_old",&result_old_a,numeq_total,1,"DV");

printf("          | HEIGHTFUNCTION  | total number of equations: %10d \n",numeq_total);
printf("\n");

/*-------------------------------- allocate and initialise solver (OLL) */
hfsolv = (SOLVAR*)CCACALLOC(1,sizeof(SOLVAR));
hfsolv->fieldtyp = fluid;
#ifdef PARALLEL
#ifndef SPOOLES_PACKAGE
dserror("SPOOLES package is not compiled in");
#endif
hfsolv->solvertyp = SPOOLES_nonsym;
#else
hfsolv->solvertyp = umfpack;
#endif
hfsolv->parttyp = cut_elements;
hfsolv->matrixtyp = oll_matrix;

/* ------------------------------------------------------ create matrix */
hfsolv->nsysarray = 1;
hfsolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(hfsolv->nsysarray,sizeof(SPARSE_TYP));
hfsolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(hfsolv->nsysarray,sizeof(SPARSE_ARRAY));

hfsolv->sysarray_typ[0] = oll;
hfsolv->sysarray[0].oll = (OLL*)CCACALLOC(1,sizeof(OLL));
hfmat_oll = hfsolv->sysarray[0].oll;

/*--------------------------------------------------- initialise matrix */
hfmat_oll->numeq       = numeq; 
hfmat_oll->numeq_total = numeq_total;
hfmat_oll->nnz         =0;
hfmat_oll->is_copied   =0;
hfmat_oll->is_masked   =0;
hfmat_oll->rdim = hfmat_oll->numeq;
hfmat_oll->cdim = hfmat_oll->numeq_total;
hfmat_oll->row  = (MATENTRY**)CCACALLOC(hfmat_oll->rdim,sizeof(MATENTRY*));
hfmat_oll->col  = (MATENTRY**)CCACALLOC(hfmat_oll->cdim,sizeof(MATENTRY*));
update = amdef("update", &(hfmat_oll->update), numeq, 1, "IV");
counter=0;
/*-------------------------------------------------- fill update vector */
for (i=0;i<numnp_total;i++)
{
   actnode=&(actdis->node[i]);
   if (actnode->hfdof==NULL) continue;
   if (actnode->proc!=myrank) continue;
   update[counter]=actnode->hfdof[0];
   counter++;   
}
dsassert(counter==numeq,"number of local dofs wrong for height function!\n");
/*---------------------------------------------- sort the vector update */
qsort(update, counter, sizeof(INT), cmp_int);
hfmat_oll->total = numeq;
hfmat_oll->used  = 0;
hfmat_oll->spare = (MATENTRY*)CCACALLOC(hfmat_oll->total,sizeof(MATENTRY));

/*--------------------------------------- allocate 1 dist. vector 'rhs' */
solserv_create_vec(&rhs,1,numeq_total,numeq,"DV");
solserv_zero_vec(rhs);
/*--------------------------------------- allocate 1 dist. vector 'sol' */
solserv_create_vec(&sol,1,numeq_total,numeq,"DV");
solserv_zero_vec(sol);

/*----------------------------------------- initialise element function */
fluid_calelm_hf(actpart,actintra,action,container,hfmat_oll,rhs,0);

/*------------------------------------------ initialise solver function */
init=1;
solver_control(hfsolv, actintra,
               &(hfsolv->sysarray_typ[0]),
               &(hfsolv->sysarray[0]),
               sol,
               &(rhs[0]),
               init);

/*---------------------------------------------- fill old result vector */
for (i=0;i<numnp_total;i++)
{
   actnode = &(actdis->node[i]);
   if (actnode->hfdof==NULL) continue;
   dof = actnode->hfdof[0];
   result[dof]= actnode->sol_increment.a.da[3][numdf];
}
break; /* break for initialisation phase */

/*======================================================================*
 |                     S O L U T I O N    P H A S E                     |
 *======================================================================*/
case 2:

/*---------------------------- intitialise global matrix and global rhs */
solserv_zero_vec(&(rhs[0]));
solserv_zero_mat(actintra,&(hfsolv->sysarray[0]),
                 &(hfsolv->sysarray_typ[0]));

/*--------------------------------------------------- call the elements */
*action = calc_fluid_heightfunc;
container->dvec         = NULL;
container->ftimerhs     = NULL;
container->fiterhs      = NULL;
container->global_numeq = 0;
container->nii          = 0;
container->nif          = 0;
container->nif          = 0;
container->nim          = 0;
container->kstep        = 0;
container->fieldtyp     = fluid;
fluid_calelm_hf(actpart,actintra,action,container,hfmat_oll,rhs,1);

/*------------------------------------ solve for height function values */
init=0;
solver_control(hfsolv, actintra,
               &(hfsolv->sysarray_typ[0]),
               &(hfsolv->sysarray[0]),
               sol,
               &(rhs[0]),
               init);

/*----------------------------------------- copy solution to the nodes */
solserv_reddistvec(
                      sol,
                      &(hfsolv->sysarray[0]),
                      &(hfsolv->sysarray_typ[0]),
                      result,
                      numeq_total,
                      actintra
                     );
for (i=0;i<numnp_total;i++)
{
   actnode = &(actdis->node[i]);
   if (actnode->hfdof==NULL) continue;
   dof = actnode->hfdof[0];
   actnode->sol_increment.a.da[3][numdf] = result[dof];
}

/*-------------------------------------------------- check convergence */
switch (fdyn->itnorm)
{
case fncc_L2:
   for (i=0;i<numeq_total;i++)
   {
      dgnorm += DSQR(result[i]-result_old[i]);
       gnorm += DSQR(result[i]);
   }
   dgnorm = sqrt(dgnorm);
    gnorm = sqrt(gnorm);
break;
default:
   dserror("parameter itchk out of range!\n"); 
}
if (gnorm<EPS5) gnorm=ONE;
*grat = dgnorm/gnorm;

break; /* break for solution phase */
default:
   dserror("Parameter hctrl out of range!\n");
} /* end of switch (hctrl) */

/*----------------------------------- actual solution gets old solution */
amcopy(&result_a, &result_old_a); 

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

#else
dserror("FSI-functions not compiled in!\n");
#endif
return;
} /* end of fluid_heighfunc*/

/*!--------------------------------------------------------------------- 
\brief seperat treatment of heightfunction

<pre>                                                         genk 11/03

In this function the equation for the heightfunction = position of the
free surface is solved seperately from the fluid problem. This is done
during the nonlinear iteration.

</pre>
\param   *actpart     PARTITION           actual partition
\param   *actintra    INTRA               actual intra-com.
\param   *action      CALC_ACTION
\param   *container   CONTAINER
\param    converged   INT                 convergence flag

\return void                                            

------------------------------------------------------------------------*/
void fluid_calelm_hf(	      
                      PARTITION           *actpart,
                      INTRA               *actintra,  
                      CALC_ACTION         *action,
                      CONTAINER           *container,
                      OLL                 *hfmat_oll,
                      DIST_VECTOR         *rhs,
                      INT                  init
                     )
{
#ifdef D_FSI
INT                  i,j,jj,nd,counter,node;
INT                  rindex;
static INT           init_nnz;
static INT           myrank;
static INT           numele;
static INT           numeq;
static INT          *update;
INT                  lm[MAXDOFPERELE];         /* location vector for this element */
#ifdef PARALLEL
INT                  owner[MAXDOFPERELE];      /* the owner of every dof */
#endif
static DOUBLE        **estif;
static DOUBLE         *etforce;
static DOUBLE         *eiforce;
ELEMENT               *actele;
MATENTRY             **row;          /* matrix column                               */
MATENTRY              *actentry;     /* actual matrix entry                         */


#ifdef DEBUG 
dstrc_enter("fluid_calelm_hf");
#endif

if (init==0)
{
   init_nnz=0;
   numele       = actpart->pdis[0].numele;
   update       = hfmat_oll->update.a.iv;
   numeq        = hfmat_oll->numeq;
   myrank       = actintra->intra_rank;
   estif        = estif_global.a.da;
   etforce      = etforce_global.a.dv;
   eiforce      = eiforce_global.a.dv;
   goto end;
}

/*--------------------------------------- loop elements at free surface */
for (i=0;i<numele;i++)
{
   actele = actpart->pdis[0].element[i];
   switch(actele->eltyp)
   {
#ifdef D_FLUID2   
   case el_fluid2: 
      if (actele->e.f2->fs_on!=3) continue;
      fluid2(actpart,actintra,actele,NULL,
             &estif_global,NULL,
             &etforce_global,&eiforce_global,NULL,
	     action,NULL,NULL,container);
   break;
#endif
#ifdef D_FLUID3   
   case el_fluid3:
      if (actele->e.f3->fs_on!=3) continue;
      fluid3(actpart,actintra,actele,
             &estif_global,NULL,
             &etforce_global,&eiforce_global,NULL,
	     action,NULL,NULL,container);   
   break;
#endif
   default:
      dserror("eltyp unknown!\n");
   }
   /*--------------------------------------------------------- assemble */
   counter=0;
   for (j=0; j<container->ngnode; j++)
   {
      node = container->iedgnod[j];
      dsassert(actele->node[node]->hfdof!=NULL,"cannot read hfdof\n");
      lm[counter] = actele->node[node]->hfdof[0];
#ifdef PARALLEL 
      owner[counter] = actele->node[node]->proc;
#endif
      counter++;
   }/* end of loop over element nodes */
   /*--------------------------------------- now start looping the dofs */
   nd = counter;
   /*------------------------------------ loop over i (the element row) */
   for (j=0; j<nd; j++)
   {
      jj = lm[j];
      /*----------------------------------------- loop only my own rows */
#ifdef PARALLEL 
      if (owner[j]!=myrank) continue;
#endif
      /* --------------------------------------------- add complete row */
      oll_addrow(hfmat_oll, jj, lm, estif[j], nd);
      /*-------------------------------------------------- assemble rhs */
      rindex = find_index(jj,update,numeq);
      dsassert(rindex>=0,"dof does not exist in update\n");
      rhs[0].vec.a.dv[rindex]+=eiforce[j];
   }/* end loop over i */   
} 

/*----------------------------------------------------------- count nnz */
if (init_nnz==0)
{
   counter=0;
   row   = hfmat_oll->row;
   for (i=0;i<hfmat_oll->rdim;i++)
   {
      actentry = row[i];
      while (actentry!=NULL)
      {
         counter++;
	 actentry = actentry->rnext;
      }
   }
   hfmat_oll->nnz=counter;
   init_nnz++;
}

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#else
dserror("FSI-functions not compiled in!\n");
#endif

return;
} /* end of fluid_heighfunc*/

#endif
/*! @} (documentation module close)*/
