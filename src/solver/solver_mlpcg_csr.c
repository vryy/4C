/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/prototypes_mlpcg.h"
int cmp_int(const void *a, const void *b );
double cmp_double(const void *a, const void *b );
/*! 
\addtogroup MLPCG 
*//*! @{ (documentation module open)*/
/*!----------------------------------------------------------------------
\brief the multilevel preconditioner main structure

<pre>                                                         m.gee 09/02    
defined in solver_mlpcg.c
</pre>

*----------------------------------------------------------------------*/
extern struct _MLPRECOND mlprecond;


/*!---------------------------------------------------------------------
\brief extract a local column from a   DBCSR  matrix                                         

<pre>                                                        m.gee 10/02 

</pre>
\param P          DBCSR*       (i/o) the Prolongator
\param actcol     int          (i)   the column to be extracted
\param block      double[][500](i)   working matrix
\param rindex     int*         (i)   row indize of block
\param nrow       int*         (i)   row dimension of block
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_extractcollocal(DBCSR *P, int actcol, double *col, 
                           int *rcol, int *nrow)
{
int           i,j,counter;
int           *ia,*ja,*update,numeq;
double        *a;
int            actrow,colstart,colend;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_extractcollocal");
#endif
/*----------------------------------------------------------------------*/
*nrow   = 0;
counter = 0;
numeq   = P->numeq;
ia      = P->ia.a.iv;
ja      = P->ja.a.iv;
update  = P->update.a.iv;
a       = P->a.a.dv;
/*----------------------------------------------------------------------*/
for (i=0; i<numeq; i++)
{
   actrow   = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   for (j=colstart; j<colend; j++)
   {
      if (ja[j]==actcol)
      {
         col[counter]  = a[j];
         rcol[counter] = actrow;
         dsassert(counter<1000,"aggblock row dimension too small");
         counter++;
         break;
      }
   }
}
*nrow = counter;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_extractcollocal */

/*!---------------------------------------------------------------------
\brief opens a matrix for adding                                          

<pre>                                                        m.gee 9/02 

</pre>
\param matrix      DBCSR*    (i/o) the active level in the ml-precond.                   
\param firstdof    int       (i)   first dof updated on this proc
\param lastdof     int       (i)   last dof updated on this proc
\param numeq_total int       (i)   total number of rows
\param nnz_guess   int       (i)   guess for number of nonzero entries 
\param actintra    INTRA*    (i)   the intra-communicator of this field  
\return void                                               
\sa mlpcg_csr_setblock
------------------------------------------------------------------------*/
void mlpcg_csr_open(DBCSR*  matrix,
                    int     firstdof,
                    int     lastdof,
                    int     numeq_total,
                    int     nnz_guess,
                    INTRA  *actintra)
{
int        i,j;
int        dofrangesend[MAXPROC*2],dofrangerecv[MAXPROC*2];
int        myrank,nproc;
int       *update;
int        counter;
int        ione=-1;
int        isnew;
int        rowguess;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_csr_open");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*------------------------------------------------- check the dof range */
for (i=0; i<nproc*2; i++) dofrangesend[i]=0;
dofrangesend[myrank*2]   = firstdof;
dofrangesend[myrank*2+1] = lastdof;
#ifdef PARALLEL
MPI_Allreduce(dofrangesend,dofrangerecv,nproc*2,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
for (i=0; i<nproc*2; i++) dofrangerecv[i] = dofrangesend[i];
#endif
/*------------------------------------------------------- test last dof */
if (dofrangerecv[(nproc-1)*2+1] != numeq_total-1)
dserror("inconsistent total number of rows in csr creation");
/*------------------------------------------------------- test firstdof */
if (dofrangerecv[0] != 0) dserror("Numbering of dofs doesn't start with zero");
/*---------------------------------------- test monotony over the procs */
   if (dofrangerecv[1] <= dofrangerecv[0])
   dserror("Monotony inside proc wrong");
for (i=1; i<nproc; i++)
{
   if (dofrangerecv[i*2] != dofrangerecv[(i-1)*2+1]+1)
   dserror("Monotony over procs wrong");
   if (dofrangerecv[i*2+1] <= dofrangerecv[i*2])
   dserror("Monotony inside proc wrong");
}
/*------------------------------------------ check the sums of the dofs */
counter=0;
for (i=0; i<nproc; i++) 
   counter += dofrangerecv[i*2+1]-dofrangerecv[i*2]+1;
if (counter != numeq_total)
dserror("total sum of dofs wrong");
/*-------------------------------- check, whether matrix already exists */
if (matrix->update.Typ != cca_XX) isnew=0;
else                              isnew=1;
/*------------------------------------------------------------ put size */
matrix->numeq       = dofrangerecv[myrank*2+1]-dofrangerecv[myrank*2]+1;
matrix->numeq_total = counter; 
/*----------------------------------------------------- allocate update */
if (isnew==1)
{
   update = amdef("update",&(matrix->update),matrix->numeq,1,"IV");
   for (i=0; i<matrix->numeq; i++) update[i] = firstdof+i;
   /*------------------------------------------------------ allocate ia */
   amdef("ia",&(matrix->ia),matrix->numeq+1,1,"IV");
   /*------------------------------------------------------ allocate ja */
   amdef("ja",&(matrix->ja),nnz_guess,1,"IV");
   aminit(&(matrix->ja),(void*)(&ione));
   /*------------------------------------------------------- allocate a */
   amdef("a",&(matrix->a),nnz_guess,1,"DV");
   /* split the nnz_guess nonzeros into matrix->numeq pieces and make row ptrs*/
   rowguess = (int)(nnz_guess/(matrix->numeq));
   /*----------------------- make the approximately equal size row ptrs */
   counter=0;
   for (i=0; i<=matrix->numeq; i++)
   {
      matrix->ia.a.iv[i] = counter;
      counter += rowguess;
   }
   /*------ there is a rest in ja and a which is not assigned to a line */
} /* end of if (isnew==1) */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_csr_open */


/*!---------------------------------------------------------------------
\brief close a DBCSR matrix                                         

<pre>                                                        m.gee 10/02 
remove all the remaining unneeded 'ellbow-space' from the csr format.
If the matrix is then 20% too large, redefine the size, if not, leave the
allocated memory unchanged (not worth the work).
NOTE:
After mlpcg_csr_close the matrix is an exact matching csr matrix. But
it does not have to be opened again before adding to it again. The routine
mlpcg_csr_setblock is clever enough to add new 'ellbow-space' to the matrix
if needed.
</pre>
\param matrix      DBCSR*    (i/o) the active level in the ml-precond.                   
\return void                                               
\sa mlpcg_csr_open mlpcg_csr_setblock
------------------------------------------------------------------------*/
void mlpcg_csr_close(DBCSR*   matrix)
{
int        i,j,k,offset;
int        numeq;
int       *update;
int       *ja,*ia;
double    *a;
int        size,nnz_guess;
int        colstart,colend,actrow;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_csr_close");
#endif
/*----------------------------------------------------------------------*/
update    = matrix->update.a.iv;
ja        = matrix->ja.a.iv;
ia        = matrix->ia.a.iv;
a         = matrix->a.a.dv;
numeq     = matrix->numeq;
size      = matrix->ja.fdim;
nnz_guess = ia[numeq];
/*----------------------------- loop rows and move entries to the front */
for (i=0; i<numeq; i++)
{
   actrow = i;
   colstart = ia[actrow];
   colend   = ia[actrow+1];
   offset   = 0;
   /* count number of unused entries in row */
   for (j=colstart; j<colend; j++)
   {
      if (ja[j]==-1)
         offset++;
      else
         break;
   }
   /* move all values offset to the front */
   if (offset==0) continue;
   for (k=j; k<colend; k++)
   {
      ja[k-offset] = ja[k];
       a[k-offset] = a[k];
      ja[k]        = -1;
#ifdef DEBUG /* this is probably unnecessary, test this later..... */
       a[k]        = 0.0;
#endif
   }
   /* resize the row */
   ia[actrow+1] = colend-offset;
}
/*------------------------------ redefine, if too much memory is wasted */
/*------------------------- if the true size is more then 20% too large */
if ((int)(1.2*ia[numeq]) < size)
{
   amredef(&(matrix->ja),ia[numeq],1,"IV");
   amredef(&(matrix->a) ,ia[numeq],1,"DV");
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_csr_close */



/*!---------------------------------------------------------------------
\brief destroy a DBCSR matrix                                         

<pre>                                                        m.gee 10/02 
check for alocated memory and free everything
</pre>
\param matrix      DBCSR*    (i/o) the active level in the ml-precond.                   
\return void                                               
\sa mlpcg_csr_open mlpcg_csr_setblock mlpcg_csr_close
------------------------------------------------------------------------*/
void mlpcg_csr_destroy(DBCSR*   matrix)
{
int        i,j,k;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_csr_destroy");
#endif
/*----------------------------------------------------------------------*/
if (matrix->update.Typ != cca_XX) 
   amdel(&(matrix->update));
if (matrix->a.Typ != cca_XX) 
   amdel(&(matrix->a));
if (matrix->ja.Typ != cca_XX) 
   amdel(&(matrix->ja));
if (matrix->ia.Typ != cca_XX) 
   amdel(&(matrix->ia));
if (matrix->blocks.Typ != cca_XX) 
   amdel(&(matrix->blocks));
if (matrix->gdofrecv.Typ != cca_XX) 
   amdel(&(matrix->gdofrecv));
if (matrix->recvbuff.Typ != cca_XX) 
   amdel(&(matrix->recvbuff));
if (matrix->computebuff.Typ != cca_XX) 
   amdel(&(matrix->computebuff));
if (matrix->gdofsend.Typ != cca_XX) 
   amdel(&(matrix->gdofsend));
if (matrix->sendbuff.Typ != cca_XX) 
   amdel(&(matrix->sendbuff));
if (matrix->csc)
{
   mlpcg_csr_destroy(matrix->csc);
   matrix->csc = FREE(matrix->csc);
}
#ifdef PARALLEL
if (matrix->status != NULL) 
   FREE(matrix->status);
if (matrix->request != NULL) 
   FREE(matrix->request);
#endif   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_csr_destroy */




/*!---------------------------------------------------------------------
\brief sets a dense block into a DBCSR matrix                                         

<pre>                                                        m.gee 10/02 
puts a dense block into a BDCSR matrix in a processor local manner.
This means, that the given row indizes in rindex belong to this processor.
The DBCSR matrix has to be opened by mlpcg_csr_open before, the vectors
DBCSR.update, DBCSR.a, DBCSR.ja, DBCSR.ia have to be alloacted.
</pre>
\param matrix      DBCSR*        (i/o) the active level in the ml-precond.                   
\param block       double[][500] (i)   block of values
\param rindex      int*          (i)   row indize of block in matrix
\param cindex      int*          (i)   col indize of block in matrix
\param nrow        int           (i)   number of rows in block
\param ncol        int           (i)   number of oclumns in block
\param actintra    INTRA*        (i)   the intra-communicator of this field  
\return void                                               
\sa mlpcg_csr_open
------------------------------------------------------------------------*/
void mlpcg_csr_setblock(DBCSR*   matrix,
                        double   block[][500],
                        int     *rindex,
                        int     *cindex,
                        int      nrow, 
                        int      ncol,
                        INTRA   *actintra)
{
int        i,j,k,l,m,n,counter;
int        numeq;
int        nnz_guess;
int        new_nnz_guess;
int        new_row_guess;
int        move_row,move_index;
int        sf_col,ef_col,st_col,et_col;
int       *update;
int       *ja,*ia;
double    *a;
int        actrow,row_index;
int        actcol;
int        colstart,colend;
int        foundit;
int        rowstart,rowend;
int        hasmoved;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_csr_setblock");
#endif
/*----------------------------------------------------------------------*/
update    = matrix->update.a.iv;
ja        = matrix->ja.a.iv;
ia        = matrix->ia.a.iv;
a         = matrix->a.a.dv;
numeq     = matrix->numeq;
nnz_guess = ia[numeq];
/*------------------------------------ loop the rows in the dense block */
for (i=0; i<nrow; i++)
{
   actrow    = rindex[i];
   row_index = mlpcg_getindex(actrow,update,numeq);
   if (row_index==-1) dserror("Cannot find row on local piece of BDCSR matrix");
   /*------------------------------ get the col range of this row in ia */
   colstart = ia[row_index];
   colend   = ia[row_index+1];
   /*------------------------------ loop the columns and put all values */
   for (j=0; j<ncol; j++)
   {
      actcol = cindex[j];
      /*------- check whether actcol exists between colstart and colend */
      foundit=0;
      for (k=colstart; k<colend; k++)
         if (ja[k]==actcol)
         {
            foundit++;
            break;
         }
      /*--------------------------- if foundit, set value at this point */
      if (foundit==1)
         a[k] = block[i][j];
      /*------ if not found, create a new column entry if there is room */
      /*-------------------------------- and then resort the entire row */
      else
      {
         if (ja[colstart]==-1) /*------------------- yes, there is room */
         {
            /*------------ set the new column entry to the start of row */
            ja[colstart] = actcol;
            /*------------------------------- set the approbiate values */
            a[colstart]  = block[i][j];
            /*--------------------------- resort the row in ja AND in a */
            mg_sort(&(ja[colstart]),colend-colstart,NULL,&(a[colstart]));
         }
         else /*--- no, there is no room, row is full, needs enlargment */
         {
            printf("RANK %d: Enlargment of csr matrix\n",actintra->intra_rank);
            new_nnz_guess = (int)(2.0*nnz_guess);
            new_row_guess = (int)(new_nnz_guess/numeq);
            new_nnz_guess = new_row_guess*numeq;
            ja = amredef(&(matrix->ja),new_nnz_guess,1,"IV");
            a  = amredef(&(matrix->a) ,new_nnz_guess,1,"DV");
            /*------------------------- have to init the new part of ja */
            for (k=nnz_guess; k<new_nnz_guess; k++) 
            ja[k]=-1;
            /*- loop the rows from the back and move values to the back */
            for (k=numeq-1; k>=0; k--)
            {
               move_index = k;
               /*------------------------------------- old column range */
               sf_col = ia[move_index];
               ef_col = ia[move_index+1];
               /*------------------------------------- new column range */
               st_col =   k    * new_row_guess;
               et_col =  (k+1) * new_row_guess;
               /* make sure, the new rowsize is larger then the old one */
               dsassert(et_col-st_col>ef_col-sf_col,"row enlargment failed");
               /*--------- assure, that from and to area do not overlap */
/*
               if (ef_col >= st_col) 
               printf("Warning: move area overlaps!!!!!!!!!!!!!!!!!!!!\n");
*/
               /*------------ loop the old row and copy to new location */
               counter  = 0;
               hasmoved = 0;
               for (l=sf_col; l<ef_col; l++)
               {
                  /* nothing to move */
                  if (ja[l]==-1) continue; 
                  /* flag, that a move happens */
                  hasmoved++;
                  /* assert target is free */
/*
                  if (ja[st_col+counter]!=-1)
                  printf("Warning: override happens!!!!!!!!!!!!!!!!!!!!\n");
*/
                  /* copy column index to target */
                  ja[st_col+counter] = ja[l];
                  /* copy value to target */
                   a[st_col+counter] =  a[l];
                  /* increase counter */
                  counter++;
               }
               /* everything has moved now, now we have to assure, that
                  there is no old stuff in the new row */
               for (l=st_col+counter; l<et_col; l++)
               {
                  ja[l]= -1;
                  a[l] =  0.0;
               }
               /* transfer of row complete, now sort it again */
               if (hasmoved != 0)
               mg_sort(&(ja[st_col]),et_col-st_col,NULL,&(a[st_col]));
               /* movement complete, set new ptr in ia */
               ia[move_index+1] = et_col;
            }
            /*----------------------------- now ther is room in the row */
            /* the row has moved, so get range again */
            colstart = ia[row_index];
            colend   = ia[row_index+1];
            if (ja[colstart]==-1) /*------------ yes, there is room now */
            {
               /*--------- set the new column entry to the start of row */
               ja[colstart] = actcol;
               /*---------------------------- set the approbiate values */
               a[colstart]  = block[i][j];
               /*------------------------ resort the row in ja AND in a */
               mg_sort(&(ja[colstart]),colend-colstart,NULL,&(a[colstart]));
            }
            else 
               dserror("Sever fatal error in redimension of BDCSR matrix");
            /*--------------------- movement finished, adjust nnz_guess */
            nnz_guess = new_nnz_guess;
         } /* end of else */
      } /* end of else */
   } /* and of for (j=0; j<ncol; j++) */
} /* end of for (i=0; i<nrow; i++) */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_csr_setblock */


/*!---------------------------------------------------------------------
\brief adds a dense block into a DBCSR matrix                                         

<pre>                                                        m.gee 10/02 
puts a dense block into a BDCSR matrix in a processor local manner.
This means, that the given row indizes in rindex belong to this processor.
The DBCSR matrix has to be opened by mlpcg_csr_open before, the vectors
DBCSR.update, DBCSR.a, DBCSR.ja, DBCSR.ia have to be alloacted.
</pre>
\param matrix      DBCSR*        (i/o) the active level in the ml-precond.                   
\param block       double[][500] (i)   block of values
\param rindex      int*          (i)   row indize of block in matrix
\param cindex      int*          (i)   col indize of block in matrix
\param nrow        int           (i)   number of rows in block
\param ncol        int           (i)   number of oclumns in block
\param actintra    INTRA*        (i)   the intra-communicator of this field  
\return void                                               
\sa mlpcg_csr_open
------------------------------------------------------------------------*/
void mlpcg_csr_addblock(DBCSR*   matrix,
                        double   block[][500],
                        int     *rindex,
                        int     *cindex,
                        int      nrow, 
                        int      ncol,
                        INTRA   *actintra)
{
int        i,j,k,l,m,n,counter;
int        numeq;
int        nnz_guess;
int        new_nnz_guess;
int        new_row_guess;
int        move_row,move_index;
int        sf_col,ef_col,st_col,et_col;
int       *update;
int       *ja,*ia;
double    *a;
int        actrow,row_index;
int        actcol;
int        colstart,colend;
int        foundit;
int        rowstart,rowend;
int        hasmoved;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_csr_addblock");
#endif
/*----------------------------------------------------------------------*/
update    = matrix->update.a.iv;
ja        = matrix->ja.a.iv;
ia        = matrix->ia.a.iv;
a         = matrix->a.a.dv;
numeq     = matrix->numeq;
nnz_guess = ia[numeq];
/*------------------------------------ loop the rows in the dense block */
for (i=0; i<nrow; i++)
{
   actrow    = rindex[i];
   row_index = mlpcg_getindex(actrow,update,numeq);
   if (row_index==-1) dserror("Cannot find row on local piece of BDCSR matrix");
   /*------------------------------ get the col range of this row in ia */
   colstart = ia[row_index];
   colend   = ia[row_index+1];
   /*------------------------------ loop the columns and put all values */
   for (j=0; j<ncol; j++)
   {
      actcol = cindex[j];
      /*------- check whether actcol exists between colstart and colend */
      foundit=0;
      for (k=colstart; k<colend; k++)
         if (ja[k]==actcol)
         {
            foundit++;
            break;
         }
      /*--------------------------- if foundit, set value at this point */
      if (foundit==1)
         a[k] += block[i][j];
      /*------ if not found, create a new column entry if there is room */
      /*-------------------------------- and then resort the entire row */
      else
      {
         if (ja[colstart]==-1) /*------------------- yes, there is room */
         {
            /*------------ set the new column entry to the start of row */
            ja[colstart] = actcol;
            /*------------------------------- set the approbiate values */
            a[colstart]  = block[i][j];
            /*--------------------------- resort the row in ja AND in a */
            mg_sort(&(ja[colstart]),colend-colstart,NULL,&(a[colstart]));
         }
         else /*--- no, there is no room, row is full, needs enlargment */
         {
            printf("RANK %d: Enlargment of csr matrix\n",actintra->intra_rank);
            new_nnz_guess = (int)(2.0*nnz_guess);
            new_row_guess = (int)(new_nnz_guess/numeq);
            new_nnz_guess = new_row_guess*numeq;
            ja = amredef(&(matrix->ja),new_nnz_guess,1,"IV");
            a  = amredef(&(matrix->a) ,new_nnz_guess,1,"DV");
            /*------------------------- have to init the new part of ja */
            for (k=nnz_guess; k<new_nnz_guess; k++) 
            ja[k]=-1;
            /*- loop the rows from the back and move values to the back */
            for (k=numeq-1; k>=0; k--)
            {
               move_index = k;
               /*------------------------------------- old column range */
               sf_col = ia[move_index];
               ef_col = ia[move_index+1];
               /*------------------------------------- new column range */
               st_col =   k    * new_row_guess;
               et_col =  (k+1) * new_row_guess;
               /* make sure, the new rowsize is larger then the old one */
               dsassert(et_col-st_col>ef_col-sf_col,"row enlargment failed");
               /*--------- assure, that from and to area do not overlap */
/*
               if (ef_col >= st_col) 
               printf("Warning: move area overlaps!!!!!!!!!!!!!!!!!!!!\n");
*/
               /*------------ loop the old row and copy to new location */
               counter  = 0;
               hasmoved = 0;
               for (l=sf_col; l<ef_col; l++)
               {
                  /* nothing to move */
                  if (ja[l]==-1) continue; 
                  /* flag, that a move happens */
                  hasmoved++;
                  /* assert target is free */
/*
                  if (ja[st_col+counter]!=-1)
                  printf("Warning: override happens!!!!!!!!!!!!!!!!!!!!\n");
*/
                  /* copy column index to target */
                  ja[st_col+counter] = ja[l];
                  /* copy value to target */
                   a[st_col+counter] =  a[l];
                  /* increase counter */
                  counter++;
               }
               /* everything has moved now, now we have to assure, that
                  there is no old stuff in the new row */
               for (l=st_col+counter; l<et_col; l++)
               {
                  ja[l]= -1;
                  a[l] =  0.0;
               }
               /* transfer of row complete, now sort it again */
               if (hasmoved != 0)
               mg_sort(&(ja[st_col]),et_col-st_col,NULL,&(a[st_col]));
               /* movement complete, set new ptr in ia */
               ia[move_index+1] = et_col;
            }
            /*----------------------------- now ther is room in the row */
            /* the row has moved, so get range again */
            colstart = ia[row_index];
            colend   = ia[row_index+1];
            if (ja[colstart]==-1) /*------------ yes, there is room now */
            {
               /*--------- set the new column entry to the start of row */
               ja[colstart] = actcol;
               /*---------------------------- set the approbiate values */
               a[colstart]  = block[i][j];
               /*------------------------ resort the row in ja AND in a */
               mg_sort(&(ja[colstart]),colend-colstart,NULL,&(a[colstart]));
            }
            else 
               dserror("Sever fatal error in redimension of BDCSR matrix");
            /*--------------------- movement finished, adjust nnz_guess */
            nnz_guess = new_nnz_guess;
         } /* end of else */
      } /* end of else */
   } /* and of for (j=0; j<ncol; j++) */
} /* end of for (i=0; i<nrow; i++) */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_csr_addblock */



/*!---------------------------------------------------------------------
\brief adds an entry into a DBCSR matrix                                         

<pre>                                                        m.gee 10/02 
puts a value into a BDCSR matrix in a processor local manner.
This means, that the given row indizes in rindex belong to this processor.
The DBCSR matrix has to be opened by mlpcg_csr_open before, the vectors
DBCSR.update, DBCSR.a, DBCSR.ja, DBCSR.ia have to be alloacted.
</pre>
\param matrix      DBCSR*    (i/o) the active level in the ml-precond.                   
\param val         double    (i)   value to be added to matrix
\param rindex      int       (i)   row index of value
\param cindex      int       (i)   column index of value
\param actintra    INTRA*    (i)   the intra-communicator of this field  
\return void                                               
------------------------------------------------------------------------*/
void mlpcg_csr_addentry(DBCSR*   matrix,
                        double  val,
                        int     rindex,
                        int     cindex,
                        INTRA   *actintra)
{
int        i,j,k,l,m,n,counter;
int        numeq;
int        nnz_guess;
int        new_nnz_guess;
int        new_row_guess;
int        move_row,move_index;
int        sf_col,ef_col,st_col,et_col;
int       *update;
int       *ja,*ia;
double    *a;
int        actrow,row_index;
int        actcol;
int        colstart,colend;
int        foundit;
int        rowstart,rowend;
int        hasmoved;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_csr_addentry");
#endif
/*----------------------------------------------------------------------*/
update    = matrix->update.a.iv;
ja        = matrix->ja.a.iv;
ia        = matrix->ia.a.iv;
a         = matrix->a.a.dv;
numeq     = matrix->numeq;
nnz_guess = ia[numeq];
/*------------------------------------ loop the rows in the dense block */
   actrow    = rindex;
   row_index = mlpcg_getindex(actrow,update,numeq);
   if (row_index==-1) dserror("Cannot find row on local piece of BDCSR matrix");
   /*------------------------------ get the col range of this row in ia */
   colstart = ia[row_index];
   colend   = ia[row_index+1];
   /*------------------------------ loop the columns and put all values */
      actcol = cindex;
      /*------- check whether actcol exists between colstart and colend */
      foundit=0;
      for (k=colstart; k<colend; k++)
         if (ja[k]==actcol)
         {
            foundit++;
            break;
         }
      /*--------------------------- if foundit, set value at this point */
      if (foundit==1)
         a[k] += val;
      /*------ if not found, create a new column entry if there is room */
      /*-------------------------------- and then resort the entire row */
      else
      {
         if (ja[colstart]==-1) /*------------------- yes, there is room */
         {
            /*------------ set the new column entry to the start of row */
            ja[colstart] = actcol;
            /*------------------------------- set the approbiate values */
            a[colstart]  = val;
            /*--------------------------- resort the row in ja AND in a */
            mg_sort(&(ja[colstart]),colend-colstart,NULL,&(a[colstart]));
         }
         else /*--- no, there is no room, row is full, needs enlargment */
         {
            printf("RANK %d: Enlargment of csr matrix\n",actintra->intra_rank);
            new_nnz_guess = (int)(2.0*nnz_guess);
            new_row_guess = (int)(new_nnz_guess/numeq);
            new_nnz_guess = new_row_guess*numeq;
            ja = amredef(&(matrix->ja),new_nnz_guess,1,"IV");
            a  = amredef(&(matrix->a) ,new_nnz_guess,1,"DV");
            /*------------------------- have to init the new part of ja */
            for (k=nnz_guess; k<new_nnz_guess; k++) 
            ja[k]=-1;
            /*- loop the rows from the back and move values to the back */
            for (k=numeq-1; k>=0; k--)
            {
               move_index = k;
               /*------------------------------------- old column range */
               sf_col = ia[move_index];
               ef_col = ia[move_index+1];
               /*------------------------------------- new column range */
               st_col =   k    * new_row_guess;
               et_col =  (k+1) * new_row_guess;
               /* make sure, the new rowsize is larger then the old one */
               dsassert(et_col-st_col>ef_col-sf_col,"row enlargment failed");
               /*--------- assure, that from and to area do not overlap */
/*
               if (ef_col >= st_col) 
               printf("Warning: move area overlaps!!!!!!!!!!!!!!!!!!!!\n");
*/
               /*------------ loop the old row and copy to new location */
               counter  = 0;
               hasmoved = 0;
               for (l=sf_col; l<ef_col; l++)
               {
                  /* nothing to move */
                  if (ja[l]==-1) continue; 
                  /* flag, that a move happens */
                  hasmoved++;
                  /* assert target is free */
/*
                  if (ja[st_col+counter]!=-1)
                  printf("Warning: override happens!!!!!!!!!!!!!!!!!!!!\n");
*/
                  /* copy column index to target */
                  ja[st_col+counter] = ja[l];
                  /* copy value to target */
                   a[st_col+counter] =  a[l];
                  /* increase counter */
                  counter++;
               }
               /* everything has moved now, now we have to assure, that
                  there is no old stuff in the new row */
               for (l=st_col+counter; l<et_col; l++)
               {
                  ja[l]= -1;
                  a[l] =  0.0;
               }
               /* transfer of row complete, now sort it again */
               if (hasmoved != 0)
               mg_sort(&(ja[st_col]),et_col-st_col,NULL,&(a[st_col]));
               /* movement complete, set new ptr in ia */
               ia[move_index+1] = et_col;
            }
            /*----------------------------- now ther is room in the row */
            /* the row has moved, so get range again */
            colstart = ia[row_index];
            colend   = ia[row_index+1];
            if (ja[colstart]==-1) /*------------ yes, there is room now */
            {
               /*--------- set the new column entry to the start of row */
               ja[colstart] = actcol;
               /*---------------------------- set the approbiate values */
               a[colstart]  = val;
               /*------------------------ resort the row in ja AND in a */
               mg_sort(&(ja[colstart]),colend-colstart,NULL,&(a[colstart]));
            }
            else 
               dserror("Sever fatal error in redimension of BDCSR matrix");
            /*--------------------- movement finished, adjust nnz_guess */
            nnz_guess = new_nnz_guess;
         } /* end of else */
      } /* end of else */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_csr_addentry */


/*!---------------------------------------------------------------------
\brief sets an entry into a DBCSR matrix                                         

<pre>                                                        m.gee 10/02 
puts a value into a BDCSR matrix in a processor local manner.
This means, that the given row indizes in rindex belong to this processor.
The DBCSR matrix has to be opened by mlpcg_csr_open before, the vectors
DBCSR.update, DBCSR.a, DBCSR.ja, DBCSR.ia have to be alloacted.
</pre>
\param matrix      DBCSR*    (i/o) the active level in the ml-precond.                   
\param val         double    (i)   value to be added to matrix
\param rindex      int       (i)   row index of value
\param cindex      int       (i)   column index of value
\param actintra    INTRA*    (i)   the intra-communicator of this field  
\return void                                               
------------------------------------------------------------------------*/
void mlpcg_csr_setentry(DBCSR*   matrix,
                        double  val,
                        int     rindex,
                        int     cindex,
                        INTRA   *actintra)
{
int        i,j,k,l,m,n,counter;
int        numeq;
int        nnz_guess;
int        new_nnz_guess;
int        new_row_guess;
int        move_row,move_index;
int        sf_col,ef_col,st_col,et_col;
int       *update;
int       *ja,*ia;
double    *a;
int        actrow,row_index;
int        actcol;
int        colstart,colend;
int        foundit;
int        rowstart,rowend;
int        hasmoved;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_csr_setentry");
#endif
/*----------------------------------------------------------------------*/
update    = matrix->update.a.iv;
ja        = matrix->ja.a.iv;
ia        = matrix->ia.a.iv;
a         = matrix->a.a.dv;
numeq     = matrix->numeq;
nnz_guess = ia[numeq];
/*------------------------------------ loop the rows in the dense block */
   actrow    = rindex;
   row_index = mlpcg_getindex(actrow,update,numeq);
   if (row_index==-1) dserror("Cannot find row on local piece of BDCSR matrix");
   /*------------------------------ get the col range of this row in ia */
   colstart = ia[row_index];
   colend   = ia[row_index+1];
   /*------------------------------ loop the columns and put all values */
      actcol = cindex;
      /*------- check whether actcol exists between colstart and colend */
      foundit=0;
      for (k=colstart; k<colend; k++)
         if (ja[k]==actcol)
         {
            foundit++;
            break;
         }
      /*--------------------------- if foundit, set value at this point */
      if (foundit==1)
         a[k] = val;
      /*------ if not found, create a new column entry if there is room */
      /*-------------------------------- and then resort the entire row */
      else
      {
         if (ja[colstart]==-1) /*------------------- yes, there is room */
         {
            /*------------ set the new column entry to the start of row */
            ja[colstart] = actcol;
            /*------------------------------- set the approbiate values */
            a[colstart]  = val;
            /*--------------------------- resort the row in ja AND in a */
            mg_sort(&(ja[colstart]),colend-colstart,NULL,&(a[colstart]));
         }
         else /*--- no, there is no room, row is full, needs enlargment */
         {
            printf("RANK %d: Enlargment of csr matrix\n",actintra->intra_rank);
            new_nnz_guess = (int)(2.0*nnz_guess);
            new_row_guess = (int)(new_nnz_guess/numeq);
            new_nnz_guess = new_row_guess*numeq;
            ja = amredef(&(matrix->ja),new_nnz_guess,1,"IV");
            a  = amredef(&(matrix->a) ,new_nnz_guess,1,"DV");
            /*------------------------- have to init the new part of ja */
            for (k=nnz_guess; k<new_nnz_guess; k++) 
            ja[k]=-1;
            /*- loop the rows from the back and move values to the back */
            for (k=numeq-1; k>=0; k--)
            {
               move_index = k;
               /*------------------------------------- old column range */
               sf_col = ia[move_index];
               ef_col = ia[move_index+1];
               /*------------------------------------- new column range */
               st_col =   k    * new_row_guess;
               et_col =  (k+1) * new_row_guess;
               /* make sure, the new rowsize is larger then the old one */
               dsassert(et_col-st_col>ef_col-sf_col,"row enlargment failed");
               /*--------- assure, that from and to area do not overlap */
/*
               if (ef_col >= st_col) 
               printf("Warning: move area overlaps!!!!!!!!!!!!!!!!!!!!\n");
*/
               /*------------ loop the old row and copy to new location */
               counter  = 0;
               hasmoved = 0;
               for (l=sf_col; l<ef_col; l++)
               {
                  /* nothing to move */
                  if (ja[l]==-1) continue; 
                  /* flag, that a move happens */
                  hasmoved++;
                  /* assert target is free */
/*
                  if (ja[st_col+counter]!=-1)
                  printf("Warning: override happens!!!!!!!!!!!!!!!!!!!!\n");
*/
                  /* copy column index to target */
                  ja[st_col+counter] = ja[l];
                  /* copy value to target */
                   a[st_col+counter] =  a[l];
                  /* increase counter */
                  counter++;
               }
               /* everything has moved now, now we have to assure, that
                  there is no old stuff in the new row */
               for (l=st_col+counter; l<et_col; l++)
               {
                  ja[l]= -1;
                  a[l] =  0.0;
               }
               /* transfer of row complete, now sort it again */
               if (hasmoved != 0)
               mg_sort(&(ja[st_col]),et_col-st_col,NULL,&(a[st_col]));
               /* movement complete, set new ptr in ia */
               ia[move_index+1] = et_col;
            }
            /*----------------------------- now ther is room in the row */
            /* the row has moved, so get range again */
            colstart = ia[row_index];
            colend   = ia[row_index+1];
            if (ja[colstart]==-1) /*------------ yes, there is room now */
            {
               /*--------- set the new column entry to the start of row */
               ja[colstart] = actcol;
               /*---------------------------- set the approbiate values */
               a[colstart]  = val;
               /*------------------------ resort the row in ja AND in a */
               mg_sort(&(ja[colstart]),colend-colstart,NULL,&(a[colstart]));
            }
            else 
               dserror("Sever fatal error in redimension of BDCSR matrix");
            /*--------------------- movement finished, adjust nnz_guess */
            nnz_guess = new_nnz_guess;
         } /* end of else */
      } /* end of else */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_csr_setentry */




/*!---------------------------------------------------------------------
\brief adds a row into a DBCSR matrix                                         

<pre>                                                        m.gee 10/02 
puts a dense block into a BDCSR matrix in a processor local manner.
This means, that the given row indizes in rindex belong to this processor.
The DBCSR matrix has to be opened by mlpcg_csr_open before, the vectors
DBCSR.update, DBCSR.a, DBCSR.ja, DBCSR.ia have to be alloacted.
</pre>
\param matrix      DBCSR*        (i/o) the active level in the ml-precond.                   
\param block       double[][500] (i)   block of values
\param rindex      int*          (i)   row indize of block in matrix
\param cindex      int*          (i)   col indize of block in matrix
\param nrow        int           (i)   number of rows in block
\param ncol        int           (i)   number of oclumns in block
\param actintra    INTRA*        (i)   the intra-communicator of this field  
\return void                                               
\sa mlpcg_csr_open
------------------------------------------------------------------------*/
void mlpcg_csr_addrow(DBCSR*   matrix,
                      int       rownum,
                      double   *row,
                      int     *cindex,
                      int      ncol,
                      INTRA   *actintra)
{
int        i,j,k,l,m,n,counter;
int        numeq;
int        nnz_guess;
int        new_nnz_guess;
int        new_row_guess;
int        move_row,move_index;
int        sf_col,ef_col,st_col,et_col;
int       *update;
int       *ja,*ia;
double    *a;
int        actrow,row_index;
int        actcol;
int        colstart,colend;
int        foundit;
int        rowstart,rowend;
int        hasmoved;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_csr_addrow");
#endif
/*----------------------------------------------------------------------*/
update    = matrix->update.a.iv;
ja        = matrix->ja.a.iv;
ia        = matrix->ia.a.iv;
a         = matrix->a.a.dv;
numeq     = matrix->numeq;
nnz_guess = ia[numeq];
/*------------------------------------ loop the rows in the dense block */
actrow    = rownum;
row_index = mlpcg_getindex(actrow,update,numeq);
if (row_index==-1) dserror("Cannot find row on local piece of BDCSR matrix");
/*------------------------------ get the col range of this row in ia */
colstart = ia[row_index];
colend   = ia[row_index+1];
/*------------------------------ loop the columns and put all values */
for (j=0; j<ncol; j++)
{
   actcol = cindex[j];
   /*------- check whether actcol exists between colstart and colend */
   foundit=0;
   for (k=colstart; k<colend; k++)
      if (ja[k]==actcol)
      {
         foundit++;
         break;
      }
   /*--------------------------- if foundit, set value at this point */
   if (foundit==1)
      a[k] += row[j];
   /*------ if not found, create a new column entry if there is room */
   /*-------------------------------- and then resort the entire row */
   else
   {
      if (ja[colstart]==-1) /*------------------- yes, there is room */
      {
         /*------------ set the new column entry to the start of row */
         ja[colstart] = actcol;
         /*------------------------------- set the approbiate values */
         a[colstart]  = row[j];
         /*--------------------------- resort the row in ja AND in a */
         mg_sort(&(ja[colstart]),colend-colstart,NULL,&(a[colstart]));
      }
      else /*--- no, there is no room, row is full, needs enlargment */
      {
         printf("RANK %d: Enlargment of csr matrix\n",actintra->intra_rank);
         new_nnz_guess = (int)(2.0*nnz_guess);
         new_row_guess = (int)(new_nnz_guess/numeq);
         new_nnz_guess = new_row_guess*numeq;
         ja = amredef(&(matrix->ja),new_nnz_guess,1,"IV");
         a  = amredef(&(matrix->a) ,new_nnz_guess,1,"DV");
         /*------------------------- have to init the new part of ja */
         for (k=nnz_guess; k<new_nnz_guess; k++) 
         ja[k]=-1;
         /*- loop the rows from the back and move values to the back */
         for (k=numeq-1; k>=0; k--)
         {
            move_index = k;
            /*------------------------------------- old column range */
            sf_col = ia[move_index];
            ef_col = ia[move_index+1];
            /*------------------------------------- new column range */
            st_col =   k    * new_row_guess;
            et_col =  (k+1) * new_row_guess;
            /* make sure, the new rowsize is larger then the old one */
            dsassert(et_col-st_col>ef_col-sf_col,"row enlargment failed");
            /*--------- assure, that from and to area do not overlap */
/*
               if (ef_col >= st_col) 
               printf("Warning: move area overlaps!!!!!!!!!!!!!!!!!!!!\n");
*/
            /*------------ loop the old row and copy to new location */
            counter  = 0;
            hasmoved = 0;
            for (l=sf_col; l<ef_col; l++)
            {
               /* nothing to move */
               if (ja[l]==-1) continue; 
               /* flag, that a move happens */
               hasmoved++;
               /* assert target is free */
/*
               if (ja[st_col+counter]!=-1)
               printf("Warning: override happens!!!!!!!!!!!!!!!!!!!!\n");
*/
               /* copy column index to target */
               ja[st_col+counter] = ja[l];
               /* copy value to target */
                a[st_col+counter] =  a[l];
               /* increase counter */
               counter++;
            }
            /* everything has moved now, now we have to assure, that
               there is no old stuff in the new row */
            for (l=st_col+counter; l<et_col; l++)
            {
               ja[l]= -1;
               a[l] =  0.0;
            }
            /* transfer of row complete, now sort it again */
            if (hasmoved != 0)
            mg_sort(&(ja[st_col]),et_col-st_col,NULL,&(a[st_col]));
            /* movement complete, set new ptr in ia */
            ia[move_index+1] = et_col;
         }
         /*----------------------------- now ther is room in the row */
         /* the row has moved, so get range again */
         colstart = ia[row_index];
         colend   = ia[row_index+1];
         if (ja[colstart]==-1) /*------------ yes, there is room now */
         {
            /*--------- set the new column entry to the start of row */
            ja[colstart] = actcol;
            /*---------------------------- set the approbiate values */
            a[colstart]  = row[j];
            /*------------------------ resort the row in ja AND in a */
            mg_sort(&(ja[colstart]),colend-colstart,NULL,&(a[colstart]));
         }
         else 
            dserror("Sever fatal error in redimension of BDCSR matrix");
         /*--------------------- movement finished, adjust nnz_guess */
         nnz_guess = new_nnz_guess;
      } /* end of else */
   } /* end of else */
} /* and of for (j=0; j<ncol; j++) */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_csr_addrow */
/*!---------------------------------------------------------------------
\brief makes compressed sparse column matrix from csr                                         

<pre>                                                        m.gee 10/02 
allocates csc inside the csr matrix and creates a compressed sparse column
matrix from the compressed sparse row matrix. the csr matrix has to be
closed before (no 'ellbow-room' inside)
</pre>
\param matrix          DBCSR*       (i/o) the DBCSR to be transformed
\param actintra    INTRA*    (i) the communicator                 
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_csr_csrtocsc(DBCSR *matrix, INTRA *actintra)
{
int           i,j,k,counter,num,foundit;
int           *ia,*ja,*update,numeq;
double        *a;
int           *csc_ia,*csc_ja,*csc_update,csc_numeq;
double        *csc_a;
int            actcol;
DBCSR         *csc;
int            actrow,colstart,colend;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_csr_csrtocsc");
#endif
/*----------------------------------------------------------------------*/
numeq   = matrix->numeq;
ia      = matrix->ia.a.iv;
ja      = matrix->ja.a.iv;
update  = matrix->update.a.iv;
a       = matrix->a.a.dv;
/*----------------------------------------------------------------------*/
/*------------------------------------------------- allocate csc matrix */
if (!matrix->csc)
{
   matrix->csc = (DBCSR*)CALLOC(1,sizeof(DBCSR));
   if (!(matrix->csc)) dserror("Allocation of memory failed");
}
else dserror("CSC matrix already exists in DBCSR matrix");
csc = matrix->csc;
/*---------------------------- find all columns, which exists in matrix */
/*------------------------- make an initial guess for number of columns */
csc_update = amdef("update",&(csc->update),numeq,1,"IV");
counter = 0;
num = ia[numeq];
for (i=0; i<num; i++)
{
   actcol = ja[i];
   /* look, whether actcol was found before */
   foundit = 0;
   foundit = find_index(actcol,csc_update,counter);
   /* this column already exists, continue */
   if (foundit!=-1) continue;
   /* this column is new */
   /* check whether update is large enough */
   if (counter >= csc->update.fdim)
      csc_update = amredef(&(csc->update),csc->update.fdim+500,1,"IV");
   /* add as new and resort */
   csc_update[counter] = actcol;
   counter++;
   mg_sort(csc_update,counter,NULL,NULL);
}
if (counter != csc->update.fdim)
   csc_update = amredef(&(csc->update),counter,1,"IV");
/*------------------------------- set number of columns in this matrix */
csc->numeq = csc_numeq = counter;
/*----------- set numeq_total to -1, 'cause one cannot add the columns */
csc->numeq_total = -1;
/*------------------------------------------------- allocate vector ia */
csc_ia = amdef("ia",&(csc->ia),csc_numeq+1,1,"IV");
/*------------------------------------------ allocate vectors ja and a */
csc_ja = amdef("ja",&(csc->ja),matrix->ja.fdim,1,"IV");
csc_a  = amdef("a" ,&(csc->a) ,matrix->a.fdim ,1,"DV");
/*------------------------------------------- loop and fill csc matrix */
counter   = 0;
for (i=0; i<csc_numeq; i++)
{
   csc_ia[i] = counter;
   actcol    = csc_update[i];
   /* loop csr matrix and fill values */
   for (j=0; j<numeq; j++)/* loop the row */
   {
        actrow   = update[j];
        colstart = ia[j];
        colend   = ia[j+1];
        for (k=colstart; k<colend; k++) /* loop the columns in the row */
        {
           if (ja[k]  > actcol) 
              break; /* won't find actcol anymore in this row */
           if (ja[k] != actcol) 
              continue; /* not the right column */
           csc_ja[counter] = actrow;
           csc_a[counter]  = a[k];
           counter++;
           break; /* actcol cannot be twice in the same row */
        }
   }
}
/* csr and csc have to have the same size */
dsassert(counter==matrix->ia.a.iv[matrix->numeq],"csr to csc went wrong in size check");
csc_ia[i] = counter;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_csr_csrtocsc */
/*!---------------------------------------------------------------------
\brief extract a local column from a CSC matrix                                          

<pre>                                                        m.gee 10/02 
there must not be any ellbow room in the matrix!
</pre>
\param col    int       (i) column to extract
\param update int*      (i) update vector of csc matrix 
\param ia     int*      (i) ia vector of csc matrix    
\param ja     int*      (i) ja vector of csc matrix    
\param a      double*   (i) a vector of csc matrix    
\param col    double[]  (o) column vector, maxlength 1000 !             
\param rcol   int[]     (o) row indizes vector, maxlength 1000 !             
\param nrow   int*      (o) number of rows in column
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_extractcolcsc(int     col,
                         int     numeq,
                         int    *update,
                         int    *ia,
                         int    *ja,
                         double *a,
                         double  col_out[],
                         int     rcol_out[],
                         int     *nrow)
{
int           i,j,k;
int           index,rowstart,rowend;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_extractcolcsc");
#endif
/*----------------------------------------------------------------------*/
index = find_index(col,update,numeq);
if (index==-1)
{
   *nrow = 0;
   goto exit;
}
/*-------------------------------------------------- extract the column */
rowstart = ia[index];
rowend   = ia[index+1];
*nrow    = rowend-rowstart;
if (*nrow > 1000) dserror("Column too large to extract (>1000)");
k=0;
for (i=rowstart; i<rowend; i++)
{
   col_out[k] = a[i];
   rcol_out[k] = ja[i];
   k++;
}
/*----------------------------------------------------------------------*/
exit:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_extractcolcsc */
/*!---------------------------------------------------------------------
\brief extract a local row from a CSR matrix                                          

<pre>                                                        m.gee 10/02 
there must not be any ellbow room in the matrix!
</pre>
\param col    int       (i) column to extract
\param update int*      (i) update vector of csc matrix 
\param ia     int*      (i) ia vector of csc matrix    
\param ja     int*      (i) ja vector of csc matrix    
\param a      double*   (i) a vector of csc matrix    
\param col    double[]  (o) column vector, maxlength 1000 !             
\param rcol   int[]     (o) row indizes vector, maxlength 1000 !             
\param nrow   int*      (o) number of rows in column
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_extractrowcsr(int     row,
                         int     numeq,
                         int    *update,
                         int    *ia,
                         int    *ja,
                         double *a,
                         double  row_out[],
                         int     rrow_out[],
                         int     *ncol)
{
int           i,j,k;
int           index,rowstart,rowend;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_extractrowcsr");
#endif
/*----------------------------------------------------------------------*/
index = find_index(row,update,numeq);
if (index==-1)
{
   *ncol = 0;
   goto exit;
}
/*-------------------------------------------------- extract the column */
rowstart = ia[index];
rowend   = ia[index+1];
*ncol    = rowend-rowstart;
if (*ncol > 1000) dserror("Column too large to extract (>1000)");
k=0;
for (i=rowstart; i<rowend; i++)
{
   row_out[k] = a[i];
   rrow_out[k] = ja[i];
   k++;
}
/*----------------------------------------------------------------------*/
exit:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_extractrowcsr */




/*!---------------------------------------------------------------------
\brief extract the inverse of the diagonal from csr matrix                                         

<pre>                                                        m.gee 10/02 

</pre>
\param Dinv       double*         (o)   vector to hold inverse of D=diag(csr)
\param csr        DBCSR*          (i)   the csr matrix
\param numeq      int             (i)   local number of rows in csr matrix and length of Dinv
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_csr_getdinv(double *Dinv, DBCSR *csr, int numeq)
{
int           i,j,foundit;
int           *ia,*ja,*update;
double        *a;
int            actrow,colstart,colend;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_extractcollocal");
#endif
/*----------------------------------------------------------------------*/
ia      = csr->ia.a.iv;
ja      = csr->ja.a.iv;
update  = csr->update.a.iv;
a       = csr->a.a.dv;
/*----------------------------------------------------------------------*/
foundit = 1;
for (i=0; i<numeq; i++)
{
   if (!foundit) dserror("Cannot find diagonal matrix element");
   actrow   = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   foundit  = 0;
   for (j=colstart; j<colend; j++)
   {
      if (ja[j]==actrow)
      {
         if (FABS(a[j])<EPS14) dserror("Zero diagonal detected");
         foundit = 1;
         Dinv[i] = 1.0/a[j];
         break;
      }
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_extractcollocal */







/*! @} (documentation module close)*/
