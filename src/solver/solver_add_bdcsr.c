/*!----------------------------------------------------------------------
\file
\brief contains assembly routines for multilevel preconditioned cg

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
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
                struct _DBCSR         *bdcsr2,
                struct _ARRAY         *elearray1,
                struct _ARRAY         *elearray2)
{
INT         i,j,counter;          /* some counter variables */
INT         start,index,lenght;       /* some more special-purpose counters */
INT         ii,jj;                    /* counter variables for system matrix */
INT         ii_index;                 /* place of ii in dmsr format */
INT         nd,ndnd;                  /* size of estif */
INT         nnz;                      /* number of nonzeros in sparse system matrix */
INT         numeq_total;              /* total number of equations */
INT         numeq;                    /* number of equations on this proc */
INT         lm[MAXDOFPERELE];         /* location vector for this element */
#ifdef PARALLEL
INT         owner[MAXDOFPERELE];      /* the owner of every dof */
#endif
INT         myrank;                   /* my intra-proc number */
INT         nprocs;                   /* my intra- number of processes */
DOUBLE    **estif;                    /* element matrix to be added to system matrix */
DOUBLE    **emass;                    /* element matrix to be added to system matrix */
INT        *update;                   /* csr-vector update see AZTEC manual */
INT        *ja;                       /*    "       ja           "         */
INT        *ia;                       /*    "       ia           "         */
DOUBLE     *a1,*a2;                   /*    "       a            "         */
#ifdef DEBUG
dstrc_enter("add_bdcsr");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------- set some pointers and variables */
myrank     = actintra->intra_rank;
nprocs     = actintra->intra_nprocs;
estif      = elearray1->a.da;
if (bdcsr2) emass = elearray2->a.da;
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





