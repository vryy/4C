#include "../headers/standardtypes.h"
#include "../headers/solution.h"
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 9/01    |
 | (defined in global_calelm.c, so they are extern here)                |                
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;
extern struct _ARRAY emass_global;
/*----------------------------------------------------------------------*
 |  routine to assemble element array to global UCCHB-matrix            |
 |  in parallel,taking care of coupling conditions                      |
 |                                                                      |
 |                                                                      |
 |                                                        m.gee 10/01   |
 *----------------------------------------------------------------------*/
void  add_ucchb(struct _PARTITION     *actpart,
                  struct _SOLVAR        *actsolv,
                  struct _INTRA         *actintra,
                  struct _ELEMENT       *actele,
                  struct _UCCHB         *ucchb)
{
#ifdef PARSUPERLU_PACKAGE
int         i,j,counter;              /* some counter variables */
int         ii,jj;                    /* counter variables for system matrix */
int         ii_index;
int         jj_height;
int         jj_start;
int         nd;                       /* size of estif */
int         nnz;                      /* number of nonzeros in sparse system matrix */
int         numeq_total;              /* total number of equations */
int         numeq;                    /* number of equations on this proc */
int         lm[MAXDOFPERELE];         /* location vector for this element */
int         owner[MAXDOFPERELE];      /* the owner of every dof */
int         myrank;                   /* my intra-proc number */
int         nprocs;                   /* my intra- number of processes */
double    **estif;                    /* element matrix to be added to system matrix */
double     *a;                        /* the ucchb matrix */
int        *asub;                     /* the ucchb matrix */
int        *xa;                       /* the ucchb matrix */
/*
int       **cdofs;
int         ncdofs;
*/
#ifdef DEBUG 
dstrc_enter("add_ucchb");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------- set some pointers and variables */
myrank     = actintra->intra_rank;
nprocs     = actintra->intra_nprocs;
estif      = estif_global.a.da;
nd         = actele->numnp * actele->node[0]->numdf;
nnz        = ucchb->nnz;
numeq_total= ucchb->numeq_total;
numeq      = ucchb->numeq;
a          = ucchb->a.a.dv;
asub       = ucchb->asub.a.iv;
xa         = ucchb->xa.a.iv;
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
for (i=0; i<nd; i++)
{
   ii = lm[i];
   /*------------------------------------- check for boundary condition */
   if (ii>=numeq_total) continue;
   /*------------------------------------ check for ownership of row ii */
   if (owner[i]!=myrank) continue;
   /*================================= loop over j (the element column) */
   /*                            This is the full unsymmetric version ! */
   for (j=0; j<nd; j++)
   {
      jj = lm[j];
      /*---------------------------------- check for boundary condition */
      if (jj>=numeq_total) continue;
      /*------------------ entry [ii][jj] is added to the system matrix */
      /*                       (ii here denotes the row, jj the column) */
      /*--------------------------------------------height of column jj */
      jj_height = xa[jj+1]-xa[jj];
      /*--------------------------------- start of indizes of column jj */
      jj_start  = xa[jj];
      /*--------------------------------------------------- find row ii */
      ii_index = find_index(ii,&(asub[jj_start]),jj_height);
      if (ii_index==-1) dserror("dofs severely mixed up");
      /*--------------------------------------- add value to the matrix */
      a[jj_start+ii_index] += estif[i][j];
   } /* end loop over j */
}/* end loop over i */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef PARSUPERLU_PACKAGE */
return;
} /* end of add_ucchb */



/*----------------------------------------------------------------------*
 |  make redundant ucchb matrix on all procs                 m.gee 11/01|
 *----------------------------------------------------------------------*/
void redundant_ucchb(
                        PARTITION     *actpart,
                        SOLVAR        *actsolv,
                        INTRA         *actintra,
                        UCCHB         *ucchb
                        )
{
int     i;
int     imyrank;
int     inprocs;

ARRAY   recv_a;
double *recv;

#ifdef DEBUG 
dstrc_enter("redundant_ucchb");
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
          
          But there is no MPI2 in here, yet.
*/
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
/*--- very comfortable: the only thing to do is to alreduce the ARRAY a */
/*                      (all coupling conditions are done then as well) */
/*--------------------------------------------------- allocate recvbuff */
recv = amdef("recv_a",&recv_a,ucchb->a.fdim,ucchb->a.sdim,"DV");
/*----------------------------------------------------------- Allreduce */  
MPI_Allreduce(ucchb->a.a.dv,
              recv,
              (ucchb->a.fdim)*(ucchb->a.sdim),
              MPI_DOUBLE,
              MPI_SUM,
              actintra->MPI_INTRA_COMM);
/*----------------------------------------- copy reduced data back to a */
amcopy(&recv_a,&(ucchb->a));
/*----------------------------------------------------- delete recvbuff */
amdel(&recv_a);
#endif /*---------------------------------------------- end of PARALLEL */ 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of redundant_ucchb */
