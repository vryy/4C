/*!---------------------------------------------------------------------
\file
\brief contains assembly routines for multilevel preconditioned cg

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 9/01    |
 | (defined in global_calelm.c, so they are extern here)                |                
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;
extern struct _ARRAY emass_global;
/*! 
\addtogroup MLPCG 
*//*! @{ (documentation module open)*/
/*!---------------------------------------------------------------------
\brief assemble element matrices to compressed sparse row matrix                                              

<pre>                                                        m.gee 9/02 
assemble element matrices to compressed sparse row matrix
</pre>
\param actsolv    SOLVAR*      (i)   general structure of solver informations                   
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\param bdcsr      DBCSR_ROOT*  (i)   the dbcsr matrix                 
\param sol        DIST_VECTOR* (o)   the distributed solution vector
\param rhs        DIST_VECTOR* (i)   the distributed right hand side vector
\return void                                               

------------------------------------------------------------------------*/
void  add_bdcsr(struct _PARTITION     *actpart,
                struct _SOLVAR        *actsolv,
                struct _INTRA         *actintra,
                struct _ELEMENT       *actele,
                struct _DBCSR         *bdcsr1,
                struct _DBCSR         *bdcsr2)
{
int         i,j,k,l,counter;          /* some counter variables */
int         start,index,lenght;       /* some more special-purpose counters */
int         ii,jj;                    /* counter variables for system matrix */
int         ii_owner;                 /* who is owner of dof ii -> procnumber */
int         ii_index;                 /* place of ii in dmsr format */
int         jj_index;                 /* place of jj in dmsr format */
int         nd,ndnd;                  /* size of estif */
int         nnz;                      /* number of nonzeros in sparse system matrix */
int         numeq_total;              /* total number of equations */
int         numeq;                    /* number of equations on this proc */
int         lm[MAXDOFPERELE];         /* location vector for this element */
int         owner[MAXDOFPERELE];      /* the owner of every dof */
int         myrank;                   /* my intra-proc number */
int         nprocs;                   /* my intra- number of processes */
double    **estif;                    /* element matrix to be added to system matrix */
double    **emass;                    /* element matrix to be added to system matrix */
int        *update;                   /* csr-vector update see AZTEC manual */
int        *ja;                       /*    "       ja           "         */
int        *ia;                       /*    "       ia           "         */
double     *a1,*a2;                   /*    "       a            "         */
#ifdef DEBUG 
dstrc_enter("add_bdcsr");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------- set some pointers and variables */
myrank     = actintra->intra_rank;
nprocs     = actintra->intra_nprocs;
estif      = estif_global.a.da;
if (bdcsr2) emass = emass_global.a.da;
else        emass = NULL;
nd         = actele->numnp * actele->node[0]->numdf;
ndnd       = nd*nd;
nnz        = bdcsr1->nnz;
numeq_total= bdcsr1->numeq_total;
numeq      = bdcsr1->numeq;
update     = bdcsr1->update.a.iv;
ja         = bdcsr1->ja.a.iv;
ia         = bdcsr1->ia.a.iv;
a1         = bdcsr1->a.a.dv;
if (bdcsr2) a2 = bdcsr2->a.a.dv;
else        a2 = NULL;
/*---------------------------------------------- make location vector lm*/
counter=0;
for (i=0; i<actele->numnp; i++)
{
   for (j=0; j<actele->node[i]->numdf; j++)
   {
      lm[counter]    = actele->node[i]->dof[j];
#ifdef PARALLEL 
      owner[counter] = actele->node[i]->proc;
#endif
      counter++;
   }
}/* end of loop over element nodes */
if (counter != nd) dserror("assemblage failed due to wrong dof numbering");
/*========================================== now start looping the dofs */
/*======================================= loop over i (the element row) */
for (i=0; i<nd; i++)
{
   ii = lm[i];
   /*-------------------------------------------- loop only my own rows */
#ifdef PARALLEL 
   if (owner[i]!=myrank) continue;
#endif
   /*------------------------------------- check for boundary condition */
   if (ii>=numeq_total) continue;
   ii_index    = find_index(ii,update,numeq);
   if (ii_index==-1) dserror("dof ii not found on this proc");
   /*================================= loop over j (the element column) */
   /*                            This is the full unsymmetric version ! */
   for (j=0; j<nd; j++)
   {
      jj = lm[j];
      /*---------------------------------- check for boundary condition */
      if (jj>=numeq_total) continue;
      start       = ia[ii_index];
      lenght      = ia[ii_index+1]-ia[ii_index];
      index       = find_index(jj,&(ja[start]),lenght);
      if (index==-1) dserror("dof jj not found in this row ii");
      index      += start;
      dsassert(ja[index]==jj,"Column indize wrong in assembly");
      a1[index]  += estif[i][j];
      if (bdcsr2)
      a2[index]  += emass[i][j];

   } /* end loop over j */
}/* end loop over i */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of add_bdcsr */

/*! @} (documentation module close)*/





