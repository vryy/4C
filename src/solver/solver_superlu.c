#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
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
 |  control solver lib SuperLU_MPI                       m.gee 11/01    |
 *----------------------------------------------------------------------*/
void solver_psuperlu_ucchb( 
                      struct _SOLVAR         *actsolv,
                      struct _INTRA          *actintra,
                      struct _UCCHB          *ucchb,
                      struct _DIST_VECTOR    *sol,
                      struct _DIST_VECTOR    *rhs,
                      int                     option
                     )
{
#ifdef PARSUPERLU_PACKAGE
int            i;
int            dof;
int            info;
int            nprow;
int            npcol;
int            imyrank;
int            inprocs;
double         berr[1]; 

PSUPERLUVARS  *psuperluvar;
ARRAY          b_a;
double        *b;
ARRAY          tmp_a;
double        *tmp;
 
#ifdef DEBUG 
dstrc_enter("solver_psuperlu_ucchb");
#endif
/*----------------------------------------------------------------------*/
imyrank     = actintra->intra_rank;
inprocs     = actintra->intra_nprocs;
psuperluvar = actsolv->psuperluvars;
/*----------------------------------------------------------------------*/
switch(option)
{
/*----------------------------------------------------------------------*/
/*                                                           init phase */
/*----------------------------------------------------------------------*/
case 1:
   /*------------------------ construct a 2D processor grid for SuperLU */
   switch(inprocs)
   {
   case 1:  nprow=1;   npcol=1; break;
   case 2:  nprow=1;   npcol=2; break;
   case 3:  nprow=1;   npcol=3; break;
   case 4:  nprow=2;   npcol=2; break;
   case 5:  dserror("Number of procs in INTRA_COMM not supported, choose other"); break;
   case 6:  nprow=1;   npcol=6; break;
   case 7:  dserror("Number of procs in INTRA_COMM not supported, choose other"); break;
   case 8:  nprow=2;   npcol=4; break;
   case 9:  nprow=3;   npcol=3; break;
   case 10: nprow=5;   npcol=2; break;
   case 11: dserror("Number of procs in INTRA_COMM not supported, choose other"); break;
   case 12: nprow=4;   npcol=3; break;
   case 13: dserror("Number of procs in INTRA_COMM not supported, choose other"); break;
   case 14: nprow=7;   npcol=2; break;
   case 15: nprow=5;   npcol=3; break;
   case 16: nprow=4;   npcol=4; break;
   case 17: dserror("Number of procs in INTRA_COMM not supported, choose other"); break;
   case 18: nprow=6;   npcol=3; break;
   case 19: dserror("Number of procs in INTRA_COMM not supported, choose other"); break;
   case 20: nprow=5;   npcol=4; break;
   case 21: nprow=7;   npcol=3; break;
   case 22: nprow=11;  npcol=2; break;
   case 23: dserror("Number of procs in INTRA_COMM not supported, choose other"); break;
   case 24: nprow=6;   npcol=4; break;
   case 25: nprow=5;   npcol=5; break;
   case 26: nprow=13;  npcol=2; break;
   case 27: nprow=9;   npcol=3; break;
   case 28: nprow=7;   npcol=4; break;
   case 29: dserror("Number of procs in INTRA_COMM not supported, choose other"); break;
   case 30: nprow=6;   npcol=5; break;
   case 31: dserror("Number of procs in INTRA_COMM not supported, choose other"); break;
   case 32: nprow=16;  npcol=2; break;
   case 33: nprow=11;  npcol=3; break;
   case 34: dserror("Number of procs in INTRA_COMM not supported, choose other"); break;
   case 35: nprow=7;   npcol=5; break;
   case 36: nprow=6;   npcol=6; break;
   default:     
        dserror("Number of procs in INTRA_COMM not supported, choose other"); break;
   } 
   superlu_gridinit(actintra->MPI_INTRA_COMM,nprow,npcol,&(ucchb->grid));
   /*---- make a backup copy of xa and asub, because they could be destroyed */
   /*                                                        during solution */
   am_alloc_copy(&(ucchb->asub),&(ucchb->asub_backup));
   am_alloc_copy(&(ucchb->xa),&(ucchb->xa_backup));
   am_alloc_copy(&(ucchb->asub),&(ucchb->asub_perm_backup));
   am_alloc_copy(&(ucchb->xa),&(ucchb->xa_perm_backup));
   /*--------------------------------- allocate and set the matrix structure */
   dCreate_CompCol_Matrix(&(ucchb->A), 
                          ucchb->numeq_total,
                          ucchb->numeq_total,
                          ucchb->nnz_total,
                          ucchb->a.a.dv,
                          ucchb->asub.a.iv,
                          ucchb->xa.a.iv,
                          NC,_D_super,GE);
   /*--------------------------------------------------- set default options */
   set_default_options(&(ucchb->options));
   /*------------------------------------- Initialize ScalePerm and LUstruct */
   ScalePermstructInit(ucchb->numeq_total,ucchb->numeq_total,&(ucchb->ScalePerm));
   LUstructInit(ucchb->numeq_total,ucchb->numeq_total,&(ucchb->LUstruct));
   /*------------------------------------------------------- Init statistics */
   PStatInit(&(ucchb->stat));
   /* set flag, that this matrix has been initialized and is ready for solve */   
   ucchb->is_init    = 1;
   ucchb->ncall      = 0;
   ucchb->is_factored = 0;
break;
/*----------------------------------------------------------------------*/
/*                                                    end of init phase */
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                    calculation phase */
/*----------------------------------------------------------------------*/
case 0:
   /*--------------- check whether matrix has been initialized properly */
   if (!(ucchb->is_init)) dserror("Try to solve uninitialized matrix");
   /*--------------------------------- allocate rhs and solution vector */
   b   = amdef("b",&b_a,ucchb->numeq_total,1,"DV");
         amzero(&b_a);
   tmp = amdef("tmp",&tmp_a,ucchb->numeq_total,1,"DV");
         amzero(&tmp_a);
   /*--------------------------------------------------------- fill rhs */
   for (i=0; i<rhs->numeq; i++)
   {
      dof      = ucchb->update.a.iv[i];
      tmp[dof] = rhs->vec.a.dv[i]; 
   }
   /*--------------------------------------------- allreduce rhs vector */
   #ifdef PARALLEL 
   MPI_Allreduce(tmp,b,ucchb->numeq_total,MPI_DOUBLE,MPI_SUM,actintra->MPI_INTRA_COMM);
   #endif
   amdel(&tmp_a);
   /*------------------------------------------- set additional options */
   if (ucchb->is_factored==1)/*-------- matrix has been factored before */
   {
      ucchb->options.Fact = FACTORED;
      amcopy(&(ucchb->asub_perm_backup),&(ucchb->asub));
      amcopy(&(ucchb->xa_perm_backup),&(ucchb->xa));           
   }
   else/*-------------------------- matrix has not been factored before */
   {
      /*----- destroy and create the structures ScalePerm LUstruct stat */
      if (ucchb->ncall!=0)
      {
         Destroy_LU(ucchb->numeq_total,&(ucchb->grid),&(ucchb->LUstruct));
         LUstructInit(ucchb->numeq_total,ucchb->numeq_total,&(ucchb->LUstruct));
      }
      /*------------------------------------------- factor from scratch */
      ucchb->options.Fact = DOFACT;
   }
   /*-------------------------------------------------- call the solver */
   pdgssvx_ABglobal(
                      &(ucchb->options),
                      &(ucchb->A),
                      &(ucchb->ScalePerm),
                      b,
                      ucchb->numeq_total,
                      1,
                      &(ucchb->grid),
                      &(ucchb->LUstruct),
                       berr,
                      &(ucchb->stat),
                      &info
                   );
   /*------------------------------------------------- print statistics */
   /*PStatPrint(&(ucchb->stat),&(ucchb->grid));*/
   /*
   fprintf(allfiles.out_err,"call %d buffersizes %d %d %d %d %d\n",ucchb->ncall,
                                              ucchb->LUstruct.Llu->bufmax[0],
                                              ucchb->LUstruct.Llu->bufmax[1],
                                              ucchb->LUstruct.Llu->bufmax[2],
                                              ucchb->LUstruct.Llu->bufmax[3],
                                              ucchb->LUstruct.Llu->bufmax[4]);
   fflush(allfiles.out_err);*/                                          
   /*---------------------------------------------- check for breakdown */
   if (info) 
   {
      printf("WARNING: Numerical Breakdown occured in call to SuperLU_MPI\n");
      printf("         results may be nonsense, continue....\n");
      fprintf(allfiles.out_err,"WARNING: Numerical Breakdown occured in call to SuperLU_MPI\n");
      fprintf(allfiles.out_err,"         results may be nonsense, continue....\n");
      fflush(stdout);
   }       
   /* this solution could have permuted the row and column vectors asub */
   /* and xa. For the assembly routines it is necessary to have the     */
   /* original unpermuted vectors. So the vectors are overwritten here  */
   /* the backup-copy of them, which was built in the init phase of     */
   /* this routine                                                      */ 
   amcopy(&(ucchb->asub),&(ucchb->asub_perm_backup));
   amcopy(&(ucchb->xa),&(ucchb->xa_perm_backup));           
   amcopy(&(ucchb->asub_backup),&(ucchb->asub));
   amcopy(&(ucchb->xa_backup),&(ucchb->xa));           
   /*---------- set flag, that this matrix is factored and count solves */
   ucchb->ncall++;
   ucchb->is_factored=1;
   /*------------------------------------------------------- get result */
   for (i=0; i<sol->numeq; i++)
   {
      dof      = ucchb->update.a.iv[i];
      sol->vec.a.dv[i] = b[dof];
   }
   /*----------------------- delete temporary rhs and solution vector b */
   amdel(&b_a);
break;
/*----------------------------------------------------------------------*/
/*                                             end of calculation phase */
/*----------------------------------------------------------------------*/
default:
   dserror("Unknown option for solver call to SuperLU");
break;   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef PARSUPERLU_PACKAGE */
return;
} /* end of solver_psuperlu_ucchb */




