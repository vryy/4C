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
#include "../solver/solver.h"
/*----------------------------------------------------------------------*
 | global dense matrices for element routines             m.gee 9/01    |
 | (defined in global_calelm.c, so they are extern here)                |                
 *----------------------------------------------------------------------*/
extern struct _ARRAY estif_global;
extern struct _ARRAY emass_global;
/*----------------------------------------------------------------------*
 |  routine to assemble element array to global rcptr-matrix            |
 |  in parallel,taking care of coupling conditions                      |
 |                                                                      |
 |                                                                      |
 |                                                         m.gee 4/02   |
 *----------------------------------------------------------------------*/
void  add_spo(struct _PARTITION     *actpart,
              struct _SOLVAR        *actsolv,
              struct _INTRA         *actintra,
              struct _ELEMENT       *actele,
              struct _SPOOLMAT      *spo1,
              struct _SPOOLMAT      *spo2)
{
#ifdef SPOOLES_PACKAGE
INT         i,j,k,l,counter;          /* some counter variables */
INT         istwo=0;
INT         start,index,lenght;       /* some more special-purpose counters */
INT         ii,jj;                    /* counter variables for system matrix */
INT         ii_iscouple;              /* flag whether ii is a coupled dof */
INT         ii_owner;                 /* who is owner of dof ii -> procnumber */
INT         ii_index;                 /* place of ii in dmsr format */
INT         jj_index;                 /* place of jj in dmsr format */
INT         nd,ndnd;                  /* size of estif */
INT         nnz;                      /* number of nonzeros in sparse system matrix */
INT         numeq_total;              /* total number of equations */
INT         numeq;                    /* number of equations on this proc */
INT         lm[MAXDOFPERELE];         /* location vector for this element */
INT         owner[MAXDOFPERELE];      /* the owner of every dof */
INT         myrank;                   /* my intra-proc number */
INT         nprocs;                   /* my intra- number of processes */
DOUBLE    **estif;                    /* element matrix to be added to system matrix */
DOUBLE    **emass;                    /* element matrix to be added to system matrix */
INT        *update;                   /* vector update see AZTEC manual */
DOUBLE     *A_loc;                    /*    "       A_loc see MUMPS manual */
DOUBLE     *B_loc;                    /*    "       A_loc see MUMPS manual */
INT        *irn;                      /*    "       irn see MUMPS manual */
INT        *jcn;                      /*    "       jcn see MUMPS manual */
INT        *rowptr;                   /*    "       rowptr see rc_ptr structure */
INT       **cdofs;                    /* list of coupled dofs and there owners, see init_assembly */
INT         ncdofs;                   /* total number of coupled dofs */
INT       **isend1;                   /* pointer to sendbuffer to communicate coupling conditions */
DOUBLE    **dsend1;                   /* pointer to sendbuffer to communicate coupling conditions */
INT       **isend2;                   /* pointer to sendbuffer to communicate coupling conditions */
DOUBLE    **dsend2;                   /* pointer to sendbuffer to communicate coupling conditions */
INT         nsend;
#ifdef DEBUG 
dstrc_enter("add_spo");
#endif
/*----------------------------------------------------------------------*/
/*----------------------- check whether to assemble one or two matrices */
if (spo2) istwo=1;
/*------------------------------------- set some pointers and variables */
myrank     = actintra->intra_rank;
nprocs     = actintra->intra_nprocs;
estif      = estif_global.a.da;
emass      = emass_global.a.da;
nd         = actele->numnp * actele->node[0]->numdf;
ndnd       = nd*nd;
nnz        = spo1->nnz;
numeq_total= spo1->numeq_total;
numeq      = spo1->numeq;
update     = spo1->update.a.iv;
A_loc      = spo1->A_loc.a.dv;
if (istwo)
B_loc      = spo2->A_loc.a.dv;
irn        = spo1->irn_loc.a.iv;
jcn        = spo1->jcn_loc.a.iv;
rowptr     = spo1->rowptr.a.iv;
cdofs      = actpart->pdis[0].coupledofs.a.ia;
ncdofs     = actpart->pdis[0].coupledofs.fdim;
/*---------------------------------- put pointers to sendbuffers if any */
#ifdef PARALLEL 
if (spo1->couple_i_send) 
{
   isend1 = spo1->couple_i_send->a.ia;
   dsend1 = spo1->couple_d_send->a.da;
   nsend  = spo1->couple_i_send->fdim;
   if (istwo)
   {
      isend2 = spo2->couple_i_send->a.ia;
      dsend2 = spo2->couple_d_send->a.da;
   }
}
#endif
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
/* end of loop over element nodes *//* this check is not possible any more for fluid element with implicit 
free surface condition: nd not eqaual numnp*numdf!!!                    */
#if 0
if (counter != nd) dserror("assemblage failed due to wrong dof numbering");
#endif
nd = counter;
/*========================================== now start looping the dofs */
/*======================================= loop over i (the element row) */
ii_iscouple = 0;
ii_owner    = myrank;
for (i=0; i<nd; i++)
{
   ii = lm[i];
   /*-------------------------------------------- loop only my own rows */
#ifdef PARALLEL 
   if (owner[i]!=myrank) continue;
#endif
   /*------------------------------------- check for boundary condition */
   if (ii>=numeq_total) continue;
   /*------------------------------------- check for coupling condition */
#ifdef PARALLEL 
   if (ncdofs)
   {
      ii_iscouple = 0;
      ii_owner    = -1;
      add_msr_checkcouple(ii,cdofs,ncdofs,&ii_iscouple,&ii_owner,nprocs);
   }
#endif
   /*-------------------- ii is not a coupled dofs or I am master owner */
   ii_index      = find_index(ii,update,numeq);
#ifndef D_CONTACT
   if (!ii_iscouple || ii_owner==myrank)
   {

      if (ii_index==-1) dserror("dof ii not found on this proc");
      start         = rowptr[ii_index];
      lenght        = rowptr[ii_index+1]-rowptr[ii_index];

   }
#endif   
   /*================================= loop over j (the element column) */
   /*                            This is the full unsymmetric version ! */
   for (j=0; j<nd; j++)
   {
      jj = lm[j];
      /*---------------------------------- check for boundary condition */
      if (jj>=numeq_total) continue;
      /*---------------------------------- check for coupling condition */
      /* 
        coupling condition for jj is not checked, because we add to 
        row ii here, which must also hold the coupled columns jj
      */ 
      /*======================================== do main-diagonal entry */
      /*                (either not a coupled dof or I am master owner) */
      if (!ii_iscouple || ii_owner==myrank)
      {
#ifdef D_CONTACT
         add_val_spo(ii,ii_index,jj,spo1,estif[i][j],actintra);
         if (istwo)
         add_val_spo(ii,ii_index,jj,spo2,emass[i][j],actintra);
#else
         index         = find_index(jj,&(jcn[start]),lenght);
         if (index==-1) dserror("dof jj not found in this row ii");
         index        += start;
         A_loc[index] += estif[i][j];
         if (istwo)
         B_loc[index] += emass[i][j];
#endif
      }
      /*======================================== do main-diagonal entry */
      /*                           (a coupled dof and I am slave owner) */
      else
      {
         add_spo_sendbuff(ii,jj,i,j,ii_owner,isend1,dsend1,estif,nsend);
         if (istwo)
         add_spo_sendbuff(ii,jj,i,j,ii_owner,isend2,dsend2,emass,nsend);
      }


   } /* end loop over j */
}/* end loop over i */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef SPOOLES_PACKAGE */
return;
} /* end of add_spo */


/*----------------------------------------------------------------------*
 |  add value to spooles matrix (with enlargment)            m.gee 11/02|
 *----------------------------------------------------------------------*/
void add_val_spo(INT ii,INT index, INT jj, struct _SPOOLMAT *spo, DOUBLE val, INTRA *actintra)
{
#ifdef SPOOLES_PACKAGE
INT     i,j,k,l,colstart,colend,foundit;
INT     counter,hasmoved;
INT    *irn,*jcn,*update,*rptr,numeq;
DOUBLE *A;
INT     nnz,nnz_new;
INT     rsize;
INT     move_index;
INT     sf_col,ef_col,st_col,et_col;
#ifdef DEBUG 
dstrc_enter("add_val_spo");
#endif
/*----------------------------------------------------------------------*/
irn    = spo->irn_loc.a.iv;
jcn    = spo->jcn_loc.a.iv;
update = spo->update.a.iv;
numeq  = spo->numeq;
rptr   = spo->rowptr.a.iv;
A      = spo->A_loc.a.dv;
/*----------------------------------------------------------------------*/
/*index  = find_index(ii,update,numeq);*/
if (index==-1) dserror("Cannot find local dof");
colstart = rptr[index];
colend   = rptr[index+1];
/* check whether entry jj already exists */
foundit=find_index(jj,&jcn[colstart],colend-colstart);
/* found the entry jj */
if (foundit != -1)
{
   A[colstart+foundit] += val;
}
/* entry does not exists in this row */
else
{
   if (jcn[colstart]==-1) /* there is still room in this row */
   {
      irn[colstart] = ii;
      jcn[colstart] = jj;
      A[colstart] = val;
      mg_sort(&(jcn[colstart]),colend-colstart,&(irn[colstart]),&(A[colstart]));
   }
   else /* there is no room in matrix */
   {
      startenlarge:
      nnz = spo->A_loc.fdim;
      nnz_new = (INT)(1.5*nnz);
      rsize   = (INT)(nnz_new/numeq);
      nnz_new = rsize*numeq;
      irn = amredef(&(spo->irn_loc),nnz_new,1,"IV");
      jcn = amredef(&(spo->jcn_loc),nnz_new,1,"IV");
      A   = amredef(&(spo->A_loc)  ,nnz_new,1,"DV");
      /* make sure, the new rowsize is larger then the old one */
      for (l=0; l<numeq; l++)
         if (rptr[l+1]-rptr[l] > rsize-1)
            goto startenlarge;
      /* init the new part of irn and jcn */
      for (i=nnz; i<nnz_new; i++)
      {
         irn[i]=-1;
         jcn[i]=-1;
      }
      /* loop row from back to front and move them */
      for (k=numeq-1; k>=0; k--)
      {
         move_index = k;
         /* get old column range */
         sf_col = rptr[move_index];
         ef_col = rptr[move_index+1];
         /* get new column range */
         st_col = move_index * rsize;
         et_col = (move_index+1) *rsize;
         /* loop the old row and copy to new location */
         counter=0;
         hasmoved=0;
         for (l=sf_col; l<ef_col; l++)
         {
            if (jcn[l]==-1) continue;
            hasmoved++;
            jcn[st_col+counter] = jcn[l];
            irn[st_col+counter] = irn[l];
            A[st_col+counter]   = A[l];
            counter++;
         }
         /* make sure, there is no old stuff in the new row */
         for (l=st_col+counter; l<et_col; l++)
         {
            jcn[l] = -1;
            irn[l] = -1;
            A[l]   = 0.0;
         }
         /* transfer complete, now sort the values to the back of new row */
         mg_sort(&(jcn[st_col]),et_col-st_col,&(irn[st_col]),&(A[st_col]));
         /* set new ptr in rptr */
         rptr[move_index+1] = et_col;
      }
      /* now there is room in the row, add value */
      /* the row has moved, so get the range again */
      colstart = rptr[index];
      colend   = rptr[index+1];
      if (jcn[colstart]==-1)
      {
         jcn[colstart] = jj;
         irn[colstart] = ii;
         A[colstart]   = val;
         mg_sort(&(jcn[colstart]),colend-colstart,&(irn[colstart]),&(A[colstart]));
      }
      else
         dserror("Fatal error in redimensioning of system matrix");
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef SPOOLES_PACKAGE */
return;
} /* end of add_val_spo */



/*----------------------------------------------------------------------*
 |  set value to spooles matrix (with enlargment)            m.gee 11/02|
 *----------------------------------------------------------------------*/
void set_val_spo(INT ii,INT index, INT jj, struct _SPOOLMAT *spo, DOUBLE val, INTRA *actintra)
{
#ifdef SPOOLES_PACKAGE
INT     i,j,k,l,colstart,colend,foundit;
INT     counter,hasmoved;
INT    *irn,*jcn,*update,*rptr,numeq;
DOUBLE *A;
INT     nnz,nnz_new;
INT     rsize;
INT     move_index;
INT     sf_col,ef_col,st_col,et_col;
#ifdef DEBUG 
dstrc_enter("set_val_spo");
#endif
/*----------------------------------------------------------------------*/
irn    = spo->irn_loc.a.iv;
jcn    = spo->jcn_loc.a.iv;
update = spo->update.a.iv;
numeq  = spo->numeq;
rptr   = spo->rowptr.a.iv;
A      = spo->A_loc.a.dv;
/*----------------------------------------------------------------------*/
/*index  = find_index(ii,update,numeq);*/
if (index==-1) dserror("Cannot find local dof");
colstart = rptr[index];
colend   = rptr[index+1];
/* check whether entry jj already exists */
foundit=find_index(jj,&jcn[colstart],colend-colstart);
/* found the entry jj */
if (foundit != -1)
{
   A[colstart+foundit] = val;
}
/* entry does not exists in this row */
else
{
   if (jcn[colstart]==-1) /* there is still room in this row */
   {
      irn[colstart] = ii;
      jcn[colstart] = jj;
      A[colstart] = val;
      mg_sort(&(jcn[colstart]),colend-colstart,&(irn[colstart]),&(A[colstart]));
   }
   else /* there is no room in matrix */
   {
      startenlarge:
      printf("MESSAGE: Enlargment of system matrix\n");
      nnz = spo->A_loc.fdim;
      nnz_new = (INT)(1.5*nnz);
      rsize   = (INT)(nnz_new/numeq);
      nnz_new = rsize*numeq;
      irn = amredef(&(spo->irn_loc),nnz_new,1,"IV");
      jcn = amredef(&(spo->jcn_loc),nnz_new,1,"IV");
      A   = amredef(&(spo->A_loc)  ,nnz_new,1,"DV");
      /* make sure, the new rowsize is larger then the old one */
      for (l=0; l<numeq; l++)
         if (rptr[l+1]-rptr[l] > rsize-1)
            goto startenlarge;
      /* init the new part of irn and jcn */
      for (i=nnz; i<nnz_new; i++)
      {
         irn[i]=-1;
         jcn[i]=-1;
      }
      /* loop row from back to front and move them */
      for (k=numeq-1; k>=0; k--)
      {
         move_index = k;
         /* get old column range */
         sf_col = rptr[move_index];
         ef_col = rptr[move_index+1];
         /* get new column range */
         st_col = move_index * rsize;
         et_col = (move_index+1) *rsize;
         /* loop the old row and copy to new location */
         counter=0;
         hasmoved=0;
         for (l=sf_col; l<ef_col; l++)
         {
            if (jcn[l]==-1) continue;
            hasmoved++;
            jcn[st_col+counter] = jcn[l];
            irn[st_col+counter] = irn[l];
            A[st_col+counter]   = A[l];
            counter++;
         }
         /* make sure, there is no old stuff in the new row */
         for (l=st_col+counter; l<et_col; l++)
         {
            jcn[l] = -1;
            irn[l] = -1;
            A[l]   = 0.0;
         }
         /* transfer complete, now sort the values to the back of new row */
         mg_sort(&(jcn[st_col]),et_col-st_col,&(irn[st_col]),&(A[st_col]));
         /* set new ptr in rptr */
         rptr[move_index+1] = et_col;
      }
      /* now there is room in the row, add value */
      /* the row has moved, so get the range again */
      colstart = rptr[index];
      colend   = rptr[index+1];
      if (jcn[colstart]==-1)
      {
         jcn[colstart] = jj;
         irn[colstart] = ii;
         A[colstart]   = val;
         mg_sort(&(jcn[colstart]),colend-colstart,&(irn[colstart]),&(A[colstart]));
      }
      else
         dserror("Fatal error in redimensioning of system matrix");
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef SPOOLES_PACKAGE */
return;
} /* end of set_val_spo */



 
/*----------------------------------------------------------------------*
 |  finalize the assembly to the spooles matrix              m.gee 11/02|
 *----------------------------------------------------------------------*/
void close_spooles_matrix(struct _SPOOLMAT *spo, INTRA *actintra)
{
#ifdef SPOOLES_PACKAGE
INT     i,j,k,index,colstart,colend,foundit,actrow,offset;
INT    *irn,*jcn,*update,*rptr,numeq;
DOUBLE *A;

#ifdef DEBUG 
dstrc_enter("close_spooles_matrix");
#endif
/*----------------------------------------------------------------------*/
irn    = spo->irn_loc.a.iv;
jcn    = spo->jcn_loc.a.iv;
update = spo->update.a.iv;
numeq  = spo->numeq;
rptr   = spo->rowptr.a.iv;
A      = spo->A_loc.a.dv;
/*----------------------------------------------------------------------*/
/*-------------------------- loop the rows and move values to the front */
for (i=0; i<numeq; i++)
{
   actrow = i;
   colstart = rptr[actrow];
   colend   = rptr[actrow+1];
   offset   = 0;
   /* count number of unused entries in this row */
   for (j=colstart; j<colend; j++)
   {
      if (jcn[j]==-1) offset++;
      else            break;
   }
   /* do nothing, if there is no offset */
   if (offset==0) continue;
   /* move all values offset to the front */
   for (k=j; k<colend; k++)
   {
      jcn[k-offset] = jcn[k];
      irn[k-offset] = irn[k];
      A[k-offset]   = A[k];
      jcn[k]        = -1;
      irn[k]        = -1;
      A[k]          = 0.0;
   }
   /* resize the row */
   rptr[actrow+1] = colend-offset;
}
/* redefine the array */
if (spo->jcn_loc.fdim > (INT)(1.2 * rptr[numeq]))
{
   amredef(&(spo->jcn_loc),rptr[numeq],1,"IV");
   amredef(&(spo->irn_loc),rptr[numeq],1,"IV");
   amredef(&(spo->A_loc)  ,rptr[numeq],1,"DV");
}
/* set correct number of nonzeros */
spo->nnz = rptr[numeq];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef SPOOLES_PACKAGE */
return;
} /* end of close_spooles_matrix */



/*----------------------------------------------------------------------*
 |  add one spooles matrix to another                        m.gee 11/02|
 |  init=0 : to = to + from * factor                                    |
 |  init=1:  to =      from * factor                                    |
 *----------------------------------------------------------------------*/
void add_spooles_matrix(struct _SPOOLMAT *to, struct _SPOOLMAT *from,
                       DOUBLE factor, INT init, INTRA *actintra)
{
#ifdef SPOOLES_PACKAGE
INT     i,j,k,index,foundit,actrow,offset;
INT     colstart_to,colend_to;
INT     colstart_from,colend_from;
INT    *irn_to,*jcn_to,*update_to,*rptr_to,numeq_to;
DOUBLE *A_to;
INT    *irn_from,*jcn_from,*update_from,*rptr_from,numeq_from;
DOUBLE *A_from;

#ifdef DEBUG 
dstrc_enter("add_spooles_matrix");
#endif
/*----------------------------------------------------------------------*/
irn_to    = to->irn_loc.a.iv;
jcn_to    = to->jcn_loc.a.iv;
update_to = to->update.a.iv;
numeq_to  = to->numeq;
rptr_to   = to->rowptr.a.iv;
A_to      = to->A_loc.a.dv;
irn_from    = from->irn_loc.a.iv;
jcn_from    = from->jcn_loc.a.iv;
update_from = from->update.a.iv;
numeq_from  = from->numeq;
rptr_from   = from->rowptr.a.iv;
A_from      = from->A_loc.a.dv;
/*----------------------------------------------------------------------*/
if (numeq_to != numeq_from) dserror("Number of rows not equal");
/*-------------------------- loop the rows and move values to the front */
for (i=0; i<numeq_to; i++)
{
   actrow        = i;
   colstart_to   = rptr_to[actrow];
   colend_to     = rptr_to[actrow+1];
   colstart_from = rptr_from[actrow];
   colend_from   = rptr_from[actrow+1];
   for (j=colstart_from; j<colend_from; j++)
   {
       index = find_index(irn_from[j],update_to,numeq_to);
       if (!init)
       add_val_spo(irn_from[j],index,jcn_from[j],to,A_from[j]*factor,actintra);
       else
       set_val_spo(irn_from[j],index,jcn_from[j],to,A_from[j]*factor,actintra);
   }
}
close_spooles_matrix(to,actintra);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef SPOOLES_PACKAGE */
return;
} /* end of add_spooles_matrix */






















/*----------------------------------------------------------------------*
 |  fill sendbuffer isend and dsend                           m.gee 1/02|
 *----------------------------------------------------------------------*/
void add_spo_sendbuff(INT ii,INT jj,INT i,INT j,INT ii_owner,INT **isend,
                      DOUBLE **dsend,DOUBLE **estif, INT numsend)
{
#ifdef SPOOLES_PACKAGE
INT         k,l;
#ifdef DEBUG 
dstrc_enter("add_spo_sendbuff");
#endif
/*----------------------------------------------------------------------*/
for (k=0; k<numsend; k++)
{
   if (isend[k][0]==ii) break;
}
isend[k][1]  = ii_owner;
dsend[k][jj]+= estif[i][j];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef SPOOLES_PACKAGE */
return;
} /* end of add_spo_sendbuff */



/*----------------------------------------------------------------------*
 |  exchange coupled dofs and add to row/column ptr matrix    m.gee 1/02|
 *----------------------------------------------------------------------*/
void exchange_coup_spo(
                         PARTITION     *actpart,
                         SOLVAR        *actsolv,
                         INTRA         *actintra,
                         SPOOLMAT      *spo
                        )
{
#ifdef SPOOLES_PACKAGE
INT            i,j;
INT            ii,jj,ii_index;
INT            start;
INT            lenght;
INT            index;
INT            tag;
INT            source;
INT            numeq,numeq_total;
INT            numsend;
INT            numrecv;
INT           *update;
INT          **isend;
DOUBLE       **dsend;
INT          **irecv;
DOUBLE       **drecv;
INT            imyrank;
INT            inprocs;
DOUBLE        *A_loc;                    /*    "       A_loc see MUMPS manual */
INT           *irn;                      /*    "       irn see MUMPS manual */
INT           *jcn;                      /*    "       jcn see MUMPS manual */
INT           *rowptr;                   /*    "       rowptr see rc_ptr structure */

#ifdef PARALLEL 
MPI_Status    *irecv_status;
MPI_Status    *drecv_status;

MPI_Request   *isendrequest;
MPI_Request   *dsendrequest;

MPI_Comm      *ACTCOMM;
#endif

#ifdef DEBUG 
dstrc_enter("exchange_coup_spo");
#endif
/*----------------------------------------------------------------------*/
#ifdef PARALLEL 
/*----------------------------------------------------------------------*/
imyrank = actintra->intra_rank;
inprocs = actintra->intra_nprocs;
ACTCOMM = &(actintra->MPI_INTRA_COMM);
/*---------------------------------------- set some pointers and values */
numsend     = spo->numcoupsend;
numrecv     = spo->numcouprecv;
A_loc      = spo->A_loc.a.dv;
irn        = spo->irn_loc.a.iv;
jcn        = spo->jcn_loc.a.iv;
rowptr     = spo->rowptr.a.iv;
update      = spo->update.a.iv;
numeq_total = spo->numeq_total;
numeq       = spo->numeq;
if (spo->couple_i_send) isend = spo->couple_i_send->a.ia;
if (spo->couple_d_send) dsend = spo->couple_d_send->a.da;
if (spo->couple_i_recv) irecv = spo->couple_i_recv->a.ia;
if (spo->couple_d_recv) drecv = spo->couple_d_recv->a.da;
/*--------------------------------------------- allocate some envelopes */
if (numrecv)
{
   irecv_status = (MPI_Status*)CCACALLOC(numrecv,sizeof(MPI_Status));
   drecv_status = (MPI_Status*)CCACALLOC(numrecv,sizeof(MPI_Status));
   if (!irecv_status || !drecv_status) dserror("Allocation of memory failed");
}
if (numsend)
{
   isendrequest = (MPI_Request*)CCACALLOC(numsend,sizeof(MPI_Request));
   dsendrequest = (MPI_Request*)CCACALLOC(numsend,sizeof(MPI_Request));
   if ( !isendrequest || !dsendrequest) dserror("Allocation of memory failed");
}
/*-------------------------------------------- loop the dofs to be send */
/* do all non-blocking sends and don't care about arrival (wird scho' klappe)*/
/*     the only thing to care for is the order in which things are send */
for (i=0; i<numsend; i++)
{
/*            sendbuffer       lenght    typ        dest.        tag          comm      request-handle */
   MPI_Isend(&(isend[i][0]),          2,MPI_INT   ,isend[i][1],isend[i][0],(*ACTCOMM),&(isendrequest[i]));
   MPI_Isend(&(dsend[i][0]),numeq_total,MPI_DOUBLE,isend[i][1],isend[i][0],(*ACTCOMM),&(dsendrequest[i]));
}/*------------------------------------------------ end of sending loop */
/*------------------------------- now loop over the dofs to be received */
/* 
   do blocking receives, 'cause one can't add something to the system
   matrix, which has not yet arrived, easy, isn't it?
*/
for (i=0; i<numrecv; i++)
{
   /*--------------------------- use wildcards to receive first to come */
   /*          recv-buf  lenght typ     source           tag       comm              status */
   MPI_Recv(&(irecv[i][0]),2,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,(*ACTCOMM),&(irecv_status[i]));
   if (irecv_status[i].MPI_ERROR) dserror("An error in MPI - communication occured !");

   /*---------------------- the dof number was sent as tag and as entry */
   tag    = irecv_status[i].MPI_TAG;
   if (tag != irecv[i][0]) dserror("MPI messages somehow got mixed up");
   source = irecv_status[i].MPI_SOURCE;
   
   /* do not use wildcards for second recv, we know now where it should come from */
   MPI_Recv(&(drecv[i][0]),numeq_total,MPI_DOUBLE,source,tag,(*ACTCOMM),&(drecv_status[i]));   
   if (drecv_status[i].MPI_ERROR) dserror("An error in MPI - communication occured !");

   /* now add the received data properly to my own piece of sparse matrix */
   ii = tag;
   ii_index = find_index(ii,update,numeq);
   if (ii_index==-1) dserror("dof ii not found on this proc");
   start         = rowptr[ii_index];
   lenght        = rowptr[ii_index+1]-rowptr[ii_index];
   for (j=0; j<lenght; j++)
   {
      index         = start+j;
      jj            = jcn[index];
      if (jj==-1) continue;
#ifdef D_CONTACT
      add_val_spo(ii,ii_index,jj,spo,drecv[i][jj],actintra);
#else
      A_loc[index] += drecv[i][jj];
#endif
   }
}/*---------------------------------------------- end of receiving loop */
/*-------------------------------------------- free allocated MPI-stuff */
if (numrecv){CCAFREE(irecv_status);CCAFREE(drecv_status);}
if (numsend){CCAFREE(isendrequest);CCAFREE(dsendrequest);}
/*----------------------------------------------------------------------
  do a barrier, because this is the end of the assembly, the msr matrix
  is now ready for solve
*/ 
MPI_Barrier(*ACTCOMM);
#endif /*---------------------------------------------- end of PARALLEL */ 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef SPOOLES_PACKAGE */
return;
} /* end of exchange_coup_spo */
