/*!----------------------------------------------------------------------
\file
\brief 

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

*----------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 9/01    |
 | (defined in global_calelm.c, so they are extern here)                |                
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;
extern struct _ARRAY emass_global;
/*----------------------------------------------------------------------*
 |  routine to assemble element array to global DENSE-matrix            |
 |  in parallel and sequentiell,taking care of coupling conditions      |
 |                                                                      |
 |                                                                      |
 |                                                         m.gee 9/01   |
 *----------------------------------------------------------------------*/
void  add_dense(struct _PARTITION     *actpart,
                  struct _SOLVAR        *actsolv,
                  struct _INTRA         *actintra,
                  struct _ELEMENT       *actele,
                  struct _DENSE         *dense1,
                  struct _DENSE         *dense2)
{
INT               i,j,counter;           /* some counter variables */
INT               ii,jj;                 /* counter variables for system matrix */
INT               nd;                    /* size of estif */
INT               nnz;                   /* number of nonzeros in sparse system matrix */
INT               numeq_total;           /* total number of equations */
INT               numeq;                 /* number of equations on this proc */
INT               lm[MAXDOFPERELE];      /* location vector for this element */
#ifdef PARALLEL 
INT               owner[MAXDOFPERELE];   /* the owner of every dof */
#endif
INT               myrank;                /* my intra-proc number */
INT               nprocs;                /* my intra- number of processes */
enum _SOLVER_TYP *solvertyp;
DOUBLE           **estif;                 /* element matrix to be added to system matrix */
DOUBLE           **emass;                 /* second element matrix to be assembled */
DOUBLE           **A;                     /* the dense matrix */
DOUBLE           **B;                     /* the second dense matrix */
/*
INT       **cdofs;
INT         ncdofs;
*/
#ifdef DEBUG 
dstrc_enter("add_dense");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------- set some pointers and variables */
myrank     = actintra->intra_rank;
nprocs     = actintra->intra_nprocs;
solvertyp  = &(actsolv->solvertyp);
estif      = estif_global.a.da;
if (dense2)
emass      = emass_global.a.da;
else
emass      = NULL;
nd         = actele->numnp * actele->node[0]->numdf;
numeq_total= dense1->numeq_total;
numeq      = dense1->numeq;
A          = dense1->A.a.da;
/* check for assembly of 2 matrices */
if (dense2)
B          = dense2->A.a.da;
else
B          = NULL;
/*
cdofs      = actpart->coupledofs.a.ia;
ncdofs     = actpart->coupledofs.fdim;
*/
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
/*
   NOTE:
   I don't have to care for coupling at all in this case, because
   system matrix is redundant on all procs, every proc adds his part
   (also slave and master owners of a coupled dof) and the system matrix
   is then allreduced. This makes things very comfortable for the moment.
*/
/*======================================= loop over i (the element row) */
if (*solvertyp==lapack_nonsym)/*-------------- unsymmetric lapack solve */
{
   for (i=0; i<nd; i++)
   {
      ii = lm[i];
      /*---------------------------------- check for boundary condition */
      if (ii>=numeq_total) continue;
      /*--------------------------------- check for ownership of row ii */
#ifdef PARALLEL 
      if (owner[i]!=myrank) continue;
#endif
      /*============================== loop over j (the element column) */
      /*                         This is the full unsymmetric version ! */
      for (j=0; j<nd; j++)
      {
         jj = lm[j];
         /*------------------------------- check for boundary condition */
         if (jj>=numeq_total) continue;
         /*--------------------------------------- add to system matrix */
         A[jj][ii] += estif[i][j];
         if (B)
         B[jj][ii] += emass[i][j];
      } /* end loop over j */
   }/* end loop over i */
}
else/*------------------------------------------ symmetric lapack solve */
{
   for (i=0; i<nd; i++)
   {
      ii = lm[i];
      /*---------------------------------- check for boundary condition */
      if (ii>=numeq_total) continue;
      /*--------------------------------- check for ownership of row ii */
#ifdef PARALLEL 
      if (owner[i]!=myrank) continue;
#endif
      /*============================== loop over j (the element column) */
      /*                         This is the full unsymmetric version ! */
      for (j=i; j<nd; j++)
      {
         jj = lm[j];
         /*------------------------------- check for boundary condition */
         if (jj>=numeq_total) continue;
         if (i==j)/*------------------------------- main diagonal entry */
         {
            A[jj][ii] += estif[i][j];
            if (B)
            B[jj][ii] += emass[i][j];
         }
         else/*------------------------------------- off diagonal entry */
         {
            A[ii][jj] += estif[i][j];
            A[jj][ii] += estif[i][j];
            if (B) {
            B[ii][jj] += emass[i][j];
            B[jj][ii] += emass[i][j];
            }
         }
      } /* end loop over j */
   }/* end loop over i */
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of add_dense */



/*----------------------------------------------------------------------*
 |  make redundant dense matrix on all procs                 m.gee 11/01|
 *----------------------------------------------------------------------*/
void redundant_dense(
                        PARTITION     *actpart,
                        SOLVAR        *actsolv,
                        INTRA         *actintra,
                        DENSE         *dense1,
                        DENSE         *dense2
                        )
{
INT      i;
INT      imyrank;
INT      inprocs;

ARRAY    recv_a;
DOUBLE **recv;

#ifdef DEBUG 
dstrc_enter("redundant_dense");
#endif
/*----------------------------------------------------------------------*/
/*  NOTE:
          In this routine, for a relatively short time the system matrix
          exists 2 times. This takes a lot of memory and may be a 
          bottle neck!
          In MPI2 there exists a flag for in-place-Allreduce:
          
          MPI_Allreduce(MPI_IN_PLACE,
                        ucchb->a.a.dv,
                        (ucchb->a.fdim)*(ucchb->a.sdim),
                        MPI_DOUBLE,
                        MPI_SUM,
                        actintra->MPI_INTRA_COMM);
          
          But there is no MPI2 in for HP, yet.
*/
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*--- very comfortable: the only thing to do is to alreduce the ARRAY a */
/*                      (all coupling conditions are done then as well) */
/*--------------------------------------------------- allocate recvbuff */
recv = amdef("recv_a",&recv_a,dense1->A.fdim,dense1->A.sdim,"DA");
/*----------------------------------------------------------- Allreduce */  
MPI_Allreduce(dense1->A.a.da[0],
              recv[0],
              (dense1->A.fdim)*(dense1->A.sdim),
              MPI_DOUBLE,
              MPI_SUM,
              actintra->MPI_INTRA_COMM);
/*----------------------------------------- copy reduced data back to a */
amcopy(&recv_a,&(dense1->A));
/* check for presence of second system matrix and reduce this one as well */
if (dense2)
{
   /*-------------------------------------------------------- Allreduce */  
   MPI_Allreduce(dense2->A.a.da[0],
                 recv[0],
                 (dense2->A.fdim)*(dense2->A.sdim),
                 MPI_DOUBLE,
                 MPI_SUM,
                 actintra->MPI_INTRA_COMM);

   /*-------------------------------------- copy reduced data back to a */
   amcopy(&recv_a,&(dense2->A));
}
/*----------------------------------------------------- delete recvbuff */
amdel(&recv_a);
#endif /*---------------------------------------------- end of PARALLEL */ 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of redundant_dense */
