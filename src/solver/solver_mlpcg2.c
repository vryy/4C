/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0771 - 685-6572
</pre>

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/prototypes_mlpcg.h"
#include "../shell8/shell8.h"
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
\brief smoothed aggregation algebraic multigrid preconditioner                                              

<pre>                                                        m.gee 10/02 
on input r is the residuum, it is not altered on output !
on input z is zero!         it is     altered on output
This is the V-cycle preconditioner or *gamma=1
This is the W-cycle preconditioner or *gamma=2
</pre>
\param actlevel     INT          (i)   the active level in the grid hierachy
\param z_a          ARRAY*       (o)   the correction on this active level
\param r_a          ARRAY*       (i)   the residuum on this active level
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_amgVW(INT    level,
                        ARRAY *z_a,
                        ARRAY *r_a,
                        INTRA *actintra,
                        INT   *gamma)
{
INT      i,j;
DBCSR   *stiff;
MLLEVEL *actlevel;
MLLEVEL *nextlevel;
INT      nlevel;
INT      numeq;

DOUBLE  *z,*r;
ARRAY    rwork_a;
DOUBLE  *rwork;
ARRAY    zwork_a;
DOUBLE  *zwork;

ARRAY    rcwork_a;
DOUBLE  *rcwork;
ARRAY    zcwork_a;
DOUBLE  *zcwork;

DOUBLE   done=1.0;
INT      ione=1;
INT      izero=0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_amgVW");
#endif
/*----------------------------------------------------------------------*/
nlevel   = mlprecond.numlev;
actlevel = &(mlprecond.level[level]);
stiff    = actlevel->csr;
numeq    = stiff->numeq;
z        = z_a->a.dv;
r        = r_a->a.dv;
/*----------------------------------------------------------------------*/
if (level==0) 
   j = 1;
else          
   j = *gamma;
/*----------------------------------------------------------------------*/
if (level != nlevel-1)
{
   /*----------------------------------allocate working copy of r and z */
   rwork = amdef("rwork",&rwork_a,numeq,1,"DV");
   zwork = amdef("zwork",&zwork_a,numeq,1,"DV");
   /*-------------------------------------------------- copy r to rwork */
   mlpcgupdvec(rwork,r,&done,&ione,&numeq);
}
/*=================================== on coarsest level do coarse solve */
if (level == nlevel-1)
{
   mlpcg_precond_coarsesolv(z,r,actlevel,actintra);
   goto exit;
}
/*======================================================================*/


/*================================ do presmoothing if not coarsest grid */
#if 1
if (level < nlevel-1)
{
   amzero(&zwork_a);
   /* make zwork from r */
   mlpcg_precond_presmo(zwork,rwork,stiff,actlevel,actintra,level);
   /* make z += zwork */
   mlpcgupdvec(z,zwork,&done,&izero,&numeq);
   /* make rwork += -A*zwork */
   mlpcg_matvec(rwork,stiff,zwork,-1.0,0,actintra);
}
#endif
/*======================================================================*/


/*============================================ go to next coarser level */
#if 1
if (level < nlevel-1)
{
   /* set next level */
   nextlevel = actlevel+1;
   /* define r and z for use on next level */
   rcwork = amdef("rcwork",&rcwork_a,nextlevel->csr->numeq,1,"DV");
   zcwork = amdef("zcwork",&zcwork_a,nextlevel->csr->numeq,1,"DV");
   /* zero zc */
   amzero(&zcwork_a);
   /* restrict the residuum rwork to rc */
   mlpcg_precond_restrictr(rcwork,rwork,actlevel->P,nextlevel->csr,actintra);
   /* call next level */
   for (i=0; i<j; i++)
   {
      mlpcg_precond_amgVW(level+1,&zcwork_a,&rcwork_a,actintra,gamma);
   }
   /* prolongue zc to this level */
   mlpcg_precond_prolongz(zcwork,zwork,actlevel->P,nextlevel->csr,actintra);
   /* update z += zwork */
   mlpcgupdvec(z,zwork,&done,&izero,&numeq);
   /* update rwork += -A*z */
   mlpcg_matvec(rwork,stiff,zwork,-1.0,0,actintra);
   /* tidy up */
   amdel(&rcwork_a);
   amdel(&zcwork_a);
}
#endif
/*======================================================================*/

/*=============================== do postsmoothing if not coarsest grid */
#if 1
if (level < nlevel-1)
{
   amzero(&zwork_a);
   /* make zwork from rwork */
   mlpcg_precond_postsmo(zwork,rwork,stiff,actlevel,actintra,level);
   /* make z += zwork */
   mlpcgupdvec(z,zwork,&done,&izero,&numeq);
}
#endif
/*======================================================================*/


/*================================================== tidy up this level */
amdel(&rwork_a);
amdel(&zwork_a);
/*======================================================================*/

exit:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_amgVW */


/*!---------------------------------------------------------------------
\brief smoothed aggregation algebraic multigrid preconditioner                                              

<pre>                                                        m.gee 10/02 
on input r is the residuum, it is not altered on output !
on input z is zero!         it is     altered on output
This is the F-cycle preconditioner or *gamma=2
</pre>
\param actlevel     INT          (i)   the active level in the grid hierachy
\param z_a          ARRAY*       (o)   the correction on this active level
\param r_a          ARRAY*       (i)   the residuum on this active level
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_amgF(INT    level,
                        ARRAY *z_a,
                        ARRAY *r_a,
                        INTRA *actintra,
                        INT   *gamma)
{
INT      i,j;
DBCSR   *stiff;
MLLEVEL *actlevel,*nextlevel;
INT      nlevel;
INT      numeq;

DOUBLE  *z,*r;
ARRAY    rwork_a;
DOUBLE  *rwork;
ARRAY    zwork_a;
DOUBLE  *zwork;

ARRAY    rcwork_a;
DOUBLE  *rcwork;
ARRAY    zcwork_a;
DOUBLE  *zcwork;

DOUBLE   done=1.0;
INT      ione=1;
INT      izero=0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_amgF");
#endif
/*----------------------------------------------------------------------*/
nlevel   = mlprecond.numlev;
actlevel = &(mlprecond.level[level]);
stiff    = actlevel->csr;
numeq    = stiff->numeq;
z        = z_a->a.dv;
r        = r_a->a.dv;
/*----------------------------------------------------------------------*/
if (level==0) 
   j = 1;
else          
   j = *gamma;
/*-------------------------------------allocate working copy of r and z */
rwork = amdef("rwork",&rwork_a,numeq,1,"DV");
zwork = amdef("zwork",&zwork_a,numeq,1,"DV");
/*----------------------------------------------------- copy r to rwork */
mlpcgupdvec(rwork,r,&done,&ione,&numeq);

/*=================================== on coarsest level do coarse solve */
if (level == nlevel-1)
{
   mlpcg_precond_coarsesolv(z,r,actlevel,actintra);
   *gamma = 1;
   goto exit;
}
/*======================================================================*/


/*================================ do presmoothing if not coarsest grid */
#if 1
if (level < nlevel-1)
{
   /* make zwork from r */
   mlpcg_precond_presmo(zwork,rwork,stiff,actlevel,actintra,level);
   /* make rwork += -A*zwork */
   mlpcg_matvec(rwork,stiff,zwork,-1.0,0,actintra);
   /* make z += zwork */
   mlpcgupdvec(z,zwork,&done,&izero,&numeq);
}
#endif
/*======================================================================*/


/*============================================ go to next coarser level */
#if 1
if (level < nlevel-1)
{
   /* set next level */
   nextlevel = actlevel+1;
   /* define r and z for use on next level */
   rcwork = amdef("rcwork",&rcwork_a,nextlevel->csr->numeq,1,"DV");
   zcwork = amdef("zcwork",&zcwork_a,nextlevel->csr->numeq,1,"DV");
   /* zero zc */
   amzero(&zcwork_a);
   /* restrict the residuum rwork to rc */
   mlpcg_precond_restrictr(rcwork,rwork,actlevel->P,nextlevel->csr,actintra);
   for (i=0; i<j; i++)
   {
      /* call next level */
      mlpcg_precond_amgF(level+1,&zcwork_a,&rcwork_a,actintra,gamma);
   }
   /* prolongue zc to this level */
   mlpcg_precond_prolongz(zcwork,zwork,actlevel->P,nextlevel->csr,actintra);
   /* update rwork += -A*z */
   mlpcg_matvec(rwork,stiff,zwork,-1.0,0,actintra);
   /* update z += zwork */
   mlpcgupdvec(z,zwork,&done,&izero,&numeq);
   /* tidy up */
   amdel(&rcwork_a);
   amdel(&zcwork_a);
}
#endif
/*======================================================================*/

/*=============================== do postsmoothing if not coarsest grid */
#if 1
if (level < nlevel-1)
{
   /* make zwork from rwork */
   mlpcg_precond_postsmo(zwork,rwork,stiff,actlevel,actintra,level);
   /* make z += zwork */
   mlpcgupdvec(z,zwork,&done,&izero,&numeq);
   /* make rwork += -A*zwork */
   if (*gamma==2 && i==0)
   mlpcg_matvec(rwork,stiff,zwork,-1.0,0,actintra);
}
#endif
/*======================================================================*/


/*---------------------------------------------------tidy up this level */
exit:
amdel(&rwork_a);
amdel(&zwork_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_amgF */




/*!---------------------------------------------------------------------
\brief make presmoothing                                              

<pre>                                                        m.gee 10/02 

</pre>
\param z            DOUBLE*      (o)   the smoothed residuum, on input zero!
\param r            DOUBLE*      (i)   the residuum
\param csr          DBCSR*       (i)   the matrix on the active level 
\param lev          MLLEVEL*     (i)   the structure of the active level
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_presmo(DOUBLE *z, DOUBLE *r, DBCSR *csr, MLLEVEL *lev, INTRA *actintra, INT level)
{
INT j;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_presmo");
#endif
/*----------------------------------------------------------------------*/
switch(lev->presmoother)
{
case pre_none:
break;
case pre_fwdGS:
break;
case pre_Jacobi:
   if (lev->presweep < 0) dserror("number of sweeps negative for jacobi smoother");
   mlpcg_precond_smoJacobi(z,r,csr,lev->presweep,actintra);
break;
case pre_ilu:
   if (lev->presweep < 0)
   {
      j = -(lev->presweep)-level;
   if (j<0) j=0;
   }
   else
      j = lev->presweep;
   if (mlprecond.overlap && actintra->intra_nprocs>1 && level==0)
      mlpcg_precond_smo_ILUn_overlap(z,r,csr,j,actintra);
   else
      mlpcg_precond_smo_ILUn(z,r,csr,j,actintra);
break;
default:
   dserror("Unknown type of presmoother");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_presmo */
/*!---------------------------------------------------------------------
\brief make postsmoothing                                              

<pre>                                                        m.gee 10/02 

</pre>
\param z            DOUBLE*      (o)   the smoothed residuum, on input zero!
\param r            DOUBLE*      (i)   the residuum
\param csr          DBCSR*       (i)   the matrix on the active level 
\param lev          MLLEVEL*     (i)   the structure of the active level
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_postsmo(DOUBLE *z, DOUBLE *r, DBCSR *csr, MLLEVEL *lev, INTRA *actintra, INT level)
{
INT j;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_postsmo");
#endif
/*----------------------------------------------------------------------*/
switch(lev->postsmoother)
{
case post_none:
break;
case post_bckGS:
break;
case post_Jacobi:
   if (lev->postsweep < 0) dserror("number of sweeps negative for jacobi smoother");
   mlpcg_precond_smoJacobi(z,r,csr,lev->postsweep,actintra);
break;
case post_ilu:
   if (lev->presweep < 0)
   {
      j = -(lev->presweep)-level;
      if (j<0) j=0;
   }
   else
      j = lev->presweep;
   if (mlprecond.overlap && actintra->intra_nprocs>1 && level==0)
      mlpcg_precond_smo_ILUn_overlap(z,r,csr,j,actintra);
   else
      mlpcg_precond_smo_ILUn(z,r,csr,j,actintra);
break;
default:
   dserror("Unknown type of postsmoother");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_postsmo */
/*!---------------------------------------------------------------------
\brief make coarsest solution                                              

<pre>                                                        m.gee 11/02 

</pre>
\param z            DOUBLE*      (o)   the smoothed residuum, on input zero!
\param r            DOUBLE*      (i)   the residuum
\param lev          MLLEVEL*     (i)   the structure of the active level
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_coarsesolv(DOUBLE *z, DOUBLE *r, MLLEVEL *lev, INTRA *actintra)
{
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_coarsesolv");
#endif
/*----------------------------------------------------------------------*/
switch(lev->coarsesolv)
{
case co_none:
break;
case co_ilu:
   mlpcg_precond_smo_ILUn(z,r,lev->csr,lev->co_ilu_n,actintra);
break;
case co_lapack:
   mlpcg_precond_lapacksolve(z,r,lev->csr,actintra);
break;
case co_spooles:
   mlpcg_precond_spoolessolve(z,r,lev->csr,actintra);
break;
default:
   dserror("Unknown type of coarse solver");
break;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_coarsesolv */

 

/*!---------------------------------------------------------------------
\brief create multilevel preconditioner                                              

<pre>                                                        m.gee 9/02 

</pre>
\param mlpcgvars    MLPCGVARS*   (i)   variables needed for mlpcg                   
\param bdcsr        DBCSR*       (i)   the distributed  csr matrix                   
\param actintra     INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_create(DBCSR     *bdcsr, 
                          MLPCGVARS *mlpcgvars,
                          INTRA     *actintra)
{
INT        i;
MLLEVEL   *actlev;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_create");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------------- this has been inited before */
if (mlprecond.fielddis==mlpcgvars->fielddis) goto exit;
/*--------------------------------------------- set field and partition */
mlprecond.reuse        = mlpcgvars->reuse;
mlprecond.fielddis     = mlpcgvars->fielddis;
mlprecond.partdis      = mlpcgvars->partdis;
mlprecond.omega        = mlpcgvars->p_omega;
mlprecond.overlap      = mlpcgvars->overlap;
mlprecond.typ          = mlpcgvars->typ;
mlprecond.gamma        = mlpcgvars->gamma;
mlprecond.numdf        = mlpcgvars->numdf;
/*----------------------------------------------------- allocate levels */
mlprecond.numlev = mlpcgvars->numlev;
mlprecond.level = (MLLEVEL*)CCACALLOC(mlprecond.numlev,sizeof(MLLEVEL));
if (!mlprecond.level) dserror("Allocation of MLPRECOND failed");
/*-------------- loop levels and set smoothers and prolongation damping */
for (i=0; i<mlprecond.numlev-1; i++)
{
   actlev               = &(mlprecond.level[i]);
   actlev->coarsesolv   = co_none;
   actlev->co_ilu_n     = 0;
   actlev->presmoother  = mlpcgvars->presmoother;
   actlev->presweep     = mlpcgvars->presweep;
   actlev->postsmoother = mlpcgvars->postsmoother;
   actlev->postsweep    = mlpcgvars->postsweep;
}
   actlev               = &(mlprecond.level[mlprecond.numlev-1]);
   actlev->coarsesolv   = mlpcgvars->coarsesolv;
   actlev->co_ilu_n     = mlpcgvars->co_ilu_n;
   actlev->presmoother  = pre_none;
   actlev->presweep     = 0;
   actlev->postsmoother = post_none;
   actlev->postsweep    = 0;
/*-- loop levels again and allocate a csr matrix on all levels except 0 */
for (i=1; i<mlprecond.numlev; i++)
{
   actlev               = &(mlprecond.level[i]);
   actlev->csr = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
   if (!actlev->csr) dserror("Allocation of DBCSR failed");
}
/*--------------------- set pointer in lowest level to the bdcsr matrix */
mlprecond.level[0].csr = bdcsr;
/*----------------------------------------------------------------------*/
exit:;
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_create */



/*!---------------------------------------------------------------------
\brief init multilevel preconditioner                                              

<pre>                                                        m.gee 9/02 

</pre>
\param bdcsr      DBCSR*       (i)   the distributed  csr matrix                   
\param actintra   INTRA*       (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_init(DBCSR  *bdcsr,MLPCGVARS *mlpcgvars, INTRA *actintra)
{
INT        i;
MLLEVEL   *actlev;
DOUBLE     t1,t2;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_init");
#endif
/*-------------------------------- prepare the levels and the smoothers */
actlev = &(mlprecond.level[0]);
/*============================================== make level0 separately */
/*------------------------------ init the smoothers on the active level */

/*------------------------------------------------------------take time */
t1 = ds_cputime();
/* make aggregates on the lowest level, using the grid and the partition */
if (mlprecond.ncall==0)
   mlpcg_precond_agg(actlev,actintra);
/*------- set dof numbers in the domain decomposition of the aggregates */
if (mlprecond.ncall==0 && mlprecond.typ==2)
   mlpcg_precond_aggsetdofs(actlev,mlpcgvars->numdf,actintra);
/*------------------------------------------------------------take time */
t2 = ds_cputime();
/*------------------------------------------ print time for aggregation */
if (actintra->intra_rank==0)
{
   printf("Time aggregation : %20.10f\n",t2-t1);
   fflush(stdout);
}
/*------------------------------------------------------------take time */
t1 = ds_cputime();
/*---------------------------------- create prolongator for finest grid */
if (mlprecond.mod==0 && mlprecond.typ==2) /* prolongator vanek-style */
   mlpcg_precond_P0(actlev,actintra);
if (mlprecond.mod==0 && mlprecond.typ==1) /* prolongator fish-style */
   mlpcg_precond_P_fish(actlev,actintra);
/*------------------------------------------------------------take time */
t2 = ds_cputime();
/*------------------------------------------ print time for prolongator */
if (actintra->intra_rank==0)
{
   printf("Time P lev 0     : %20.10f\n",t2-t1);
   fflush(stdout);
}
/*------------------------------------------------------------take time */
t1 = ds_cputime();
/*------------------------------- restrict the matrix to the next level */
if (mlprecond.mod==0)
   mlpcg_precond_restrictK(actlev,actlev+1,actintra);
/*------------------------------------------------------------take time */
t2 = ds_cputime();
/*------------------------------------------------- print time for PTKP */
if (actintra->intra_rank==0)
{
   printf("Time PTKP        : %20.10f\n",t2-t1);
   fflush(stdout);
}
/*------------------------------------------------------- make printout */
if (actintra->intra_rank==0)
{
   printf("level 0: size %d nnz %d\n",mlprecond.level[0].csr->numeq_total,
                                      mlprecond.level[0].csr->ja.fdim);
   fflush(stdout);
}				   
#if 0
mlpcg_printfmatrix(mlprecond.level[0].csr,actintra);
#endif                                   
/*-------------------------------------------- loop level 1 to numlev-2 */
for (i=1; i<mlprecond.numlev-1; i++)
{
   actlev = &(mlprecond.level[i]);
   /*---------------------------------------------------------take time */
   t1 = ds_cputime();
   /*---------------------------------- make aggregates on active level */
   if (mlprecond.ncall==0)
      mlpcg_precond_agg(actlev,actintra);
   /*-------------------------------- set dof numbers of the aggregates */
   if (mlprecond.ncall==0 && mlprecond.typ==2)
      mlpcg_precond_aggsetdofs(actlev,mlpcgvars->numdf,actintra);
   /*---------------------------------------------------------take time */
   t2 = ds_cputime();
   /*--------------------------------------- print time for aggregation */
   if (actintra->intra_rank==0)
   {
      printf("Time aggregation : %20.10f\n",t2-t1);
      fflush(stdout);
   }
   /*---------------------------------------------------------take time */
   t1 = ds_cputime();
   /*----------------------------------------------- create prolongator */
   if (mlprecond.mod==0 && mlprecond.typ==2) /* prolongator vanek-style */
      mlpcg_precond_P(actlev,actintra);
   if (mlprecond.mod==0 && mlprecond.typ==1) /* prolongator fish-style */
      mlpcg_precond_P_fish(actlev,actintra);
   /*---------------------------------------------------------take time */
   t2 = ds_cputime();
   /*--------------------------------------- print time for prolongator */
   if (actintra->intra_rank==0)
   {
      printf("Time P lev %d     : %20.10f\n",i,t2-t1);
      fflush(stdout);
   }
   /*---------------------------------------------------------take time */
   t1 = ds_cputime();
   /*-------------------------------- restrict matrix to the next level */
   if (mlprecond.mod==0)
      mlpcg_precond_restrictK(actlev,actlev+1,actintra);
   /*---------------------------------------------------------take time */
   t2 = ds_cputime();
   /*---------------------------------------------- print time for PTKP */
   if (actintra->intra_rank==0)
   {
      printf("Time PTKP        : %20.10f\n",t2-t1);
      fflush(stdout);
   }
   /*---------------------------------------------------- make printout */
   if (actintra->intra_rank==0)
   {
      printf("level %d: size %d nnz %d\n",i,mlprecond.level[i].csr->numeq_total,
                                            mlprecond.level[i].csr->ja.fdim);
      fflush(stdout);
   }					    
#if 0
   mlpcg_printfmatrix(mlprecond.level[i].csr,actintra);
#endif                                   
}
/*------------------------ make level numlev-1 (coarse grid) separately */
/*------------------------------------------------------- make printout */
if (actintra->intra_rank==0)
{
   printf("level %d: size %d nnz %d\n",i,mlprecond.level[i].csr->numeq_total,
                                         mlprecond.level[i].csr->ja.fdim);
   fflush(stdout);
}					 
#if 0
mlpcg_printfmatrix(mlprecond.level[i].csr,actintra);
#endif                                   
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_init */


/*!---------------------------------------------------------------------
\brief prolonge z down the finer level                                              

<pre>                                                        m.gee 11/02 
z = P * zc
</pre>
\param zc         DOUBLE*     (o)   residuum in the dimension of the next level                   
\param z          DOUBLE*     (i)   residuum on the active level
\param P          DBCSR*      (i)   the prolongator matrix
\param coarsecsr  DBCSR*      (i)   the fine level stiffness matrix
\param actintra   INTRA*      (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_prolongz(DOUBLE *zc, DOUBLE *z, DBCSR *P, DBCSR *coarsecsr, 
                            INTRA *actintra)
{
INT        i,j,n,m,counter;
INT        actrow,actcol,colstart,colend,index,owner,tag,length;
INT        myrank,nproc;
INT       *update,*ia,*ja,numeq;
DOUBLE    *a;
INT        ineeds[MAXPROC][MAXPROC],ineedr[MAXPROC][MAXPROC];
INT        nsend=0,nrecv=0;
ARRAY      irecv_a,drecv_a;
INT       *irecv;
DOUBLE    *drecv;
#ifdef PARALLEL
MPI_Request *request;
MPI_Status   status;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_prolongz");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
numeq  = P->numeq;
update = P->update.a.iv;
ia     = P->ia.a.iv;
ja     = P->ja.a.iv;
a      = P->a.a.dv;
/*----------------------------------------------------------------------*/
for (i=0; i<numeq; i++) z[i] = 0.0;

/*================================ make interproc multiplication part I */
#ifdef PARALLEL
if (nproc>1)
{
/*------------------------------ make list of owners I need from remote */
   for (n=0; n<nproc; n++)
   for (m=0; m<nproc; m++) ineeds[n][m] = 0;
   index = mlpcg_getindex(P->firstcoupledof,update,numeq);
   for (i=ia[index]; i<ia[numeq]; i++)
/*
   for (i=0; i<ia[numeq]; i++)
*/
   {
      actcol = ja[i];
      owner  = mlpcg_getowner(actcol,coarsecsr->owner,nproc);
      if (owner==myrank) 
         continue;
      ineeds[myrank][owner] = 1;
   }
   MPI_Allreduce(&(ineeds[0][0]),&(ineedr[0][0]),MAXPROC*MAXPROC,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
/*------------------------------------------- count the nrecv and nsend */
   for (n=0; n<nproc; n++)
   {
      nrecv += ineedr[myrank][n];
      nsend += ineedr[n][myrank];
   }
/*------------------------------------------- allocate nsend*2 requests */
   request = (MPI_Request*)CCAMALLOC(nsend*2*sizeof(MPI_Request));
/*----------------------- send my zc and update to the requesting procs */
   counter=0;
   for (n=0; n<nproc; n++)
   {
      if (ineedr[n][myrank]==0) 
         continue;
      MPI_Isend(coarsecsr->update.a.iv,
                coarsecsr->numeq,
                MPI_INT,
                n,
                counter,
                actintra->MPI_INTRA_COMM,
                &(request[counter]));
      counter++;
      MPI_Isend(zc,
                coarsecsr->numeq,
                MPI_DOUBLE,
                n,
                counter,
                actintra->MPI_INTRA_COMM,
                &(request[counter]));
      counter++;
   }
   dsassert(counter==2*nsend,"Number of sends wrong");
}
#endif
/*======================================================================*/



/*============================== make the local piece of multiplication */
for (i=0; i<numeq; i++)
{
   actrow   = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol = ja[j];
      owner = mlpcg_getowner(actcol,coarsecsr->owner,nproc);
      if (owner != myrank)
         continue;
      index = mlpcg_getindex(actcol,coarsecsr->update.a.iv,coarsecsr->numeq);
      if (index==-1) dserror("Cannot find local dof");
      z[i] += a[j] * zc[index];
   }
}
/*======================================================================*/

/*=============================== make interproc multiplication part II */
/*----------------------------- allocate a guess of the receive buffers */
#ifdef PARALLEL
if (nproc>1)
{
   irecv = amdef("tmp",&irecv_a,(INT)(1.2*coarsecsr->numeq),1,"IV");
   drecv = amdef("tmp",&drecv_a,(INT)(1.2*coarsecsr->numeq),1,"DV");
   while(nrecv != 0)
   {
      /* probe for incoming message */
/*
      flag=0;
      while(!flag)
         MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,actintra->MPI_INTRA_COMM,&flag,&status);
*/
      MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,actintra->MPI_INTRA_COMM,&status);
      /* get sender and tag and length */
      n = status.MPI_SOURCE;
      tag = status.MPI_TAG;
      MPI_Get_count(&status,MPI_INT,&length);
      /* resize receivebuffers, if necessary */
      if (length>irecv_a.fdim)
      {
         amdel(&irecv_a);
         amdel(&drecv_a);
         irecv = amdef("tmp",&irecv_a,length,1,"IV");
         drecv = amdef("tmp",&drecv_a,length,1,"DV");
      }
      /* receive the messages */
      MPI_Recv(irecv,length,MPI_INT   ,n,tag  ,actintra->MPI_INTRA_COMM,&status);
      MPI_Recv(drecv,length,MPI_DOUBLE,n,tag+1,actintra->MPI_INTRA_COMM,&status);
      nrecv--;
      /* process multiplication with this message */
      for (i=0; i<numeq; i++)
      {
         actrow   = update[i];
         colstart = ia[i];
         colend   = ia[i+1];
         for (j=colstart; j<colend; j++)
         {
            actcol = ja[j];
            owner = mlpcg_getowner(actcol,coarsecsr->owner,nproc);
            if (owner != n)
               continue;
            index = mlpcg_getindex(actcol,irecv,length);
            if (index==-1) dserror("Cannot find dof in received message");
            z[i] += a[j] * drecv[index];
         }
      }
   }
/*------------------------------------ wait for sent messages to finish */
for (i=0; i<nsend*2; i++)
   MPI_Wait(&(request[i]),&status);
}
#endif
/*======================================================================*/




/*============================================================= tidy up */
#ifdef PARALLEL
if (nproc>1)
{
   CCAFREE(request);
   amdel(&irecv_a);
   amdel(&drecv_a);
}
#endif
/*======================================================================*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_prolongz */



/*!---------------------------------------------------------------------
\brief restrict the residuum to the next level                                              

<pre>                                                        m.gee 11/02 
rc = Ptransposed * r
</pre>
\param rc         DOUBLE*     (o)   residuum in the dimension of the next level                   
\param r          DOUBLE*     (i)   residuum on the active level
\param P          DBCSR*      (i)   the prolongator matrix
\param coarsecsr  DBCSR*      (i)   the coarse level stiffness matrix
\param actintra   INTRA*      (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_restrictr(DOUBLE *rc, DOUBLE *r, DBCSR *P, DBCSR *coarsecsr, 
                            INTRA *actintra)
{
INT        i,j,n,m,counter;
INT        myrank,nproc;
INT        sender,tag,length;
INT        numeq,*update,*ia,*ja;
DOUBLE    *a;
INT        actrow,actcol,colstart,colend,owner;
INT        cindex;
INT        sendtos[MAXPROC][MAXPROC],sendtor[MAXPROC][MAXPROC];
INT        nsend=0,nrecv=0,sendsize;
ARRAY      send_rindex_a,send_val_a,recv_rindex_a,recv_val_a;
INT       *send_rindex,*recv_rindex;
DOUBLE    *send_val,*recv_val;
#ifdef PARALLEL
MPI_Request *request=NULL;
MPI_Status   status;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_restrictr");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------- if the column version of P doesn't exist, create it */
numeq       = P->numeq;
update      = P->update.a.iv;
ia          = P->ia.a.iv;
ja          = P->ja.a.iv;
a           = P->a.a.dv;
/*--------------------------- make the matvec init of the coarse matrix */
mlpcg_matvec_init(coarsecsr,actintra);
/*----------------------------------------------------------------------*/
for (i=0; i<coarsecsr->numeq; i++) rc[i] = 0.0;
/*---------------------- loop my interproc columns, and check the owner */
/*--------------------------------------- allocate guess of sendbuffers */
#ifdef PARALLEL
if (nproc>1)
{
   send_rindex = amdef("tmp",&send_rindex_a,2000,1,"IV");
   send_val    = amdef("tmp",&send_val_a   ,2000,1,"DV");
                 amzero(&send_val_a);
   recv_rindex = amdef("tmp",&recv_rindex_a,1,1,"IV");
   recv_val    = amdef("tmp",&recv_val_a   ,1,1,"DV");
}
for (n=0; n<nproc; n++)
for (m=0; m<nproc; m++) sendtos[n][m] = 0;
counter=0;
if (nproc>1)
for (i=0; i<numeq; i++)
{
   actrow   = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol = ja[j];
      owner  = mlpcg_getowner(actcol,coarsecsr->owner,nproc);
      if (owner==myrank) continue;
      sendtos[myrank][owner]=1;
      /* enlarge sendbuffer, if necessary */
      if (counter>=send_rindex_a.fdim)
      {
         send_rindex = amredef(&send_rindex_a,send_rindex_a.fdim+2000,1,"IV");
         send_val    = amredef(&send_val_a   ,send_val_a.fdim+2000   ,1,"DV");
      }
      send_rindex[counter] = actcol;
      send_val[counter]   += a[j] * r[i];
      counter++;
   }
}
sendsize = counter;
MPI_Allreduce(sendtos,sendtor,MAXPROC*MAXPROC,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
/*-------------------------------------- count my sends and my receives */
for (n=0; n<nproc; n++)
{
   nsend += sendtor[myrank][n];
   nrecv += sendtor[n][myrank];
}
/*--------------------------------------------------- allocate requests */
if (nproc>1)
   request = (MPI_Request*)CCAMALLOC(2*nsend*sizeof(MPI_Request));
/*------------------------------------------------- send my sendbuffers */
counter=0;
for (n=0; n<nproc; n++)
{
   if (n==myrank) continue;
   if (sendtor[myrank][n]==0) continue;
   MPI_Isend(send_rindex,sendsize,MPI_INT,n,counter,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
   MPI_Isend(send_val,sendsize,MPI_DOUBLE,n,counter,actintra->MPI_INTRA_COMM,&(request[counter]));
   counter++;
}
dsassert(counter==nsend*2,"Number of sends wrong");
#endif
/*=========================make the local multiplication in daxpy style */
for (i=0; i<numeq; i++) /* loop rows of P */
{
   actrow   = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol = ja[j];
      owner  = mlpcg_getowner(actcol,coarsecsr->owner,nproc);
      if (owner != myrank)
         continue;
      cindex = mlpcg_getindex(actcol,coarsecsr->update.a.iv,coarsecsr->numeq);
      if (cindex==-1) dserror("Cannot find local dof");
      rc[cindex] += a[j] * r[i];
   }
}
/*================= make the interprocessor part of the multiplication */
#ifdef PARALLEL
while (nrecv!=0)
{
   /* probe for any incoming message */
/*
   flag=0;
   while(!flag)
      MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,actintra->MPI_INTRA_COMM,&flag,&status);
*/
   MPI_Probe(MPI_ANY_SOURCE,MPI_ANY_TAG,actintra->MPI_INTRA_COMM,&status);
   /* get size of message */
   MPI_Get_count(&status,MPI_INT,&length);
   /* get tag and sender */
   sender = status.MPI_SOURCE;
   tag    = status.MPI_TAG;
   /* resize the recvbuffers if necessary */
   if (length > recv_rindex_a.fdim)
   {
      recv_rindex = amredef(&recv_rindex_a,length,1,"IV");
      recv_val    = amredef(&recv_val_a   ,length,1,"DV");
   }
   /* receive this message */
   MPI_Recv(recv_rindex,length,MPI_INT   ,sender,tag  ,actintra->MPI_INTRA_COMM,&status);
   MPI_Recv(recv_val   ,length,MPI_DOUBLE,sender,tag+1,actintra->MPI_INTRA_COMM,&status);
   nrecv--;
   /* process this message */
   for (i=0; i<length; i++)
   {
      actcol = recv_rindex[i];
      owner  = mlpcg_getowner(actcol,coarsecsr->owner,nproc);
      if (owner != myrank) 
         continue;
      cindex = mlpcg_getindex(actcol,coarsecsr->update.a.iv,coarsecsr->numeq);
      if (cindex==-1) dserror("Cannot find local dof");
      rc[cindex] += recv_val[i];
   }
}
for (n=0; n<nsend*2; n++) 
   MPI_Wait(&(request[n]),&status);
#endif
/*==============================================================tidy up */
#ifdef PARALLEL
if (nproc>1)
{
   amdel(&send_rindex_a);
   amdel(&send_val_a);
   amdel(&recv_rindex_a);
   amdel(&recv_val_a);
   if (request)
      CCAFREE(request);
}
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_restrictr */






/*!---------------------------------------------------------------------
\brief restrict the system matrix to the next level                                              

<pre>                                                        m.gee 9/02 

</pre>
\param actlev      MLLEVEL*     (i)   the active level in the ml-precond.                   
\param nextlev     MLLEVEL*     (i/o) the next coarser level.
\param actintra   INTRA*        (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_restrictK(MLLEVEL  *actlev, MLLEVEL *nextlev, INTRA *actintra)
{
INT        i,j,counter;
INT        myrank,nproc;
INT        startdof,enddof;
INT        numeq_total,nnz_guess;
DBCSR     *work=NULL;
INT      **blocks;
INT        max=0;
AGG       *actagg;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_restrictK");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*------------ get the row numbers of the nextlev->csr stiffness matrix */
/*-------- these correspond to the columns of the prolongator actlev->P */
/* the blocks of the prolongator actlev->P are rowblocks (correspond to */
/* the rows and columns of actlev->csr) so we have to compute the column */
/* blocks from the aggregates themselfs.As the BDCSR is continous on proc*/
/* it is o.k. to take first and last dof from aggregates                */
startdof = actlev->agg[0].dof[0];
enddof   = actlev->agg[actlev->nagg-1].dof[actlev->agg[actlev->nagg-1].numdf-1];
/*------------- now we have to collect the unknown total number of dofs */
counter=0;
for (i=0; i<actlev->nagg; i++) 
{
   counter += actlev->agg[i].numdf;
   if (max < actlev->agg[i].numdf) 
      max = actlev->agg[i].numdf;
}
#ifdef PARALLEL
MPI_Allreduce(&counter,&numeq_total,1,MPI_INT,MPI_SUM,actintra->MPI_INTRA_COMM);
#else
numeq_total=counter;
#endif
/*------- as guess for the nnz serves the bigger finer grid matrix size */
nnz_guess = actlev->csr->nnz;
/*---------------- set numeq and numeq_total in the nextlev->csr matrix */
/*------------------------------------------ now we can open the matrix */
/*--- note that the variable nextlev->csr->firscoupledof is not yet set */
if (mlprecond.ncall==0)
   mlpcg_csr_open(nextlev->csr,startdof,enddof,numeq_total,nnz_guess,actintra);
else
   mlpcg_csr_zero(nextlev->csr,actintra);
/*----------------------------------------------------------------------*/
/* 
   for the restriction, we need a working csr matrix of dimensions      
   nrow        = nrow of nextlev->csr
   ncol        = nrow of actlev->csr
   numeq_total = numeq_total of nextlev->csr
   The size can be estimated to nnz_guess from actlev->csr
   The splitting of the dofs on the processores correponds to nextlev->csr
*/
work = (DBCSR*)CCACALLOC(1,sizeof(DBCSR));
mlpcg_csr_open(work,startdof,enddof,numeq_total,nnz_guess,actintra);
/*------------------------------------- we can now make the restriction */
mlpcg_precond_PtKP(actlev->P,actlev->csr,nextlev->csr,work,actlev->agg,
                   actlev->nagg,actintra);
/*------------------------------------ destroy the temporary csr matrix */
mlpcg_csr_destroy(work);
work = CCAFREE(work);
/*----------------------------------------------------------------------*/
/* 
  due to boundary conditions there may be empty columns in the prolongator.
  This leads to zero rows and columns and zero diagonal element in the
  coarse grid stiffness matrices. Without influence on the physical behaviour,
  a ONE is set on the diagonal of such rows 
*/
mlpcg_precond_checkdirich(nextlev->csr,actintra);
/*--------------------------- close the stiffness matrix on coarse grid */
mlpcg_csr_close(nextlev->csr);
/*----------------------- check for the first coupled row in the matrix */
mlpcg_precond_check_fcd(nextlev->csr,actintra);
/*--------------------------------------------------- set nnz in matrix */
nextlev->csr->nnz = nextlev->csr->ia.a.iv[nextlev->csr->numeq];
/*-------------------------- make the aggregates blocks in nextlev->csr */
if (mlprecond.ncall==0)
{
   blocks = amdef("blocks",&(nextlev->csr->blocks),actlev->nagg,max+1,"IA");
   for (i=0; i<actlev->nagg; i++)
   {
      actagg  = &(actlev->agg[i]);
      counter = 0;
      for (j=0; j<actagg->numdf; j++)
         blocks[i][1+counter++] = actagg->dof[j];
      blocks[i][0] = counter;
   }
}
/*----------------------------------------------------------------------*/
fflush(stdout);
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_restrictK */




/*!---------------------------------------------------------------------
\brief check for the first coupled row in the matrix                                             

<pre>                                                        m.gee 10/02 
build the owner array inside csr and check for the first coupled 
row in the matrix
</pre>
\param matrix     DBCSR* (i/o) the csr matrix to be ckecked
\param actintra   INTRA* (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_check_fcd(DBCSR *matrix, INTRA *actintra)
{
INT           i,j;
INT           myrank,nproc;
INT           owner;
INT           numeq,*update,*ia,*ja;
DOUBLE       *a;
DOUBLE        one=1.0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_check_fcd");
#endif
/*----------------------------------------------------------------------*/
myrank = actintra->intra_rank;
nproc  = actintra->intra_nprocs;
/*----------------------------------------------------------------------*/
numeq  = matrix->numeq;
update = matrix->update.a.iv;
ia     = matrix->ia.a.iv;
ja     = matrix->ja.a.iv;
a      = matrix->a.a.dv;
/*----------------------------------------------------------------------*/
matrix->owner[myrank][0] = update[0];
matrix->owner[myrank][1] = update[numeq-1];
/*
#ifdef PARALLEL
for (n=0; n<nproc; n++)
   MPI_Bcast(matrix->owner[n],2,MPI_INT,n,actintra->MPI_INTRA_COMM);
#endif
*/
/*----------------------------------------------------------------------*/
owner = myrank;
for (i=0; i<numeq; i++)
{
   for (j=ia[i]; j<ia[i+1]; j++)
   {
      /* most times i am the owner myself */
      if (matrix->owner[myrank][0]<=ja[j] && matrix->owner[myrank][1]>= ja[j])
         continue;
      else
      {
         matrix->firstcoupledof = update[i];
         goto exit;
      }
   }
}
/*----------------------------------------------------------------------*/
exit:
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_check_fcd */

/*!---------------------------------------------------------------------
\brief check for zero diagonal entries                                             

<pre>                                                        m.gee 10/02 
  due to boundary conditions there may be empty columns in the prolongator.
  This leads to zero rows and columns and zero diagonal element in the
  coarse grid stiffness matrices. Without influence on the physical behaviour,
  a ONE is set on the diagonal of such rows 
</pre>
\param matrix     DBCSR* (i/o) the csr matrix to be ckecked
\param actintra   INTRA* (i)   the intra-communicator of this field                  
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_checkdirich(DBCSR *matrix, INTRA *actintra)
{
INT           i,j,actrow,actcol,colstart,colend;
INT           foundit;
INT           numeq,*update,*ia,*ja;
DOUBLE       *a;
DOUBLE        one=1.0;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_checkdirich");
#endif
/*----------------------------------------------------------------------*/
numeq  = matrix->numeq;
update = matrix->update.a.iv;
ia     = matrix->ia.a.iv;
ja     = matrix->ja.a.iv;
a      = matrix->a.a.dv;
/*----------------------------------------------------------------------*/
for (i=0; i<numeq; i++)
{
   foundit  = 0;
   actrow   = update[i];
   colstart = ia[i];
   colend   = ia[i+1];
   if (colstart==colend)
   {
      printf("No entry in row with index %d\n",i);
      mlpcg_csr_addentry(matrix,one,actrow,actrow,actintra);
      continue;
   }
   for (j=colstart; j<colend; j++)
   {
      actcol = ja[j];
      if (actcol == actrow)
      {
         foundit=1;
         if (FABS(a[j])<EPS14) 
         a[j] = 1.0;
         break;
      }
   }
   if (!foundit)
   {
      mlpcg_csr_addentry(matrix,one,actrow,actrow,actintra);
      numeq  = matrix->numeq;
      update = matrix->update.a.iv;
      ia     = matrix->ia.a.iv;
      ja     = matrix->ja.a.iv;
      a      = matrix->a.a.dv;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_checkdirich */


/*!---------------------------------------------------------------------
\brief get directors from the shell8 element for rbm's                                              

<pre>                                                        m.gee 9/02 

</pre>
\return void                                               

------------------------------------------------------------------------*/
void mlpcg_precond_getdirs(void)
{
INT           i,j,counter;
PARTDISCRET  *actpdis;
NODE         *actnode;
ELEMENT      *actele;
DOUBLE      **director;
DOUBLE        fac;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_getdirs");
#endif
/*----------------------------------------------------------------------*/
/*----------------------------------- check presence of shell8 elements */
#ifndef D_SHELL8
dserror("SHELL8 not compiled in! MLPCG is a SHELL8 solver!");
#endif
/*-------------------------------------------- get pointer to partition */   
actpdis = mlprecond.partdis;
/*------------------------- allocate pointers to directors in mlprecond */
if (mlprecond.director.Typ==cca_XX)
{
   director = amdef("direct",&(mlprecond.director),actpdis->numnp,3,"DA");
   mlprecond.node = (NODE**)CCAMALLOC((actpdis->numnp)*sizeof(NODE*));
}   
/*------------------- loop all nodes in actpdis and get there directors */
/* 
IMPORTANT NOTE:
in the nonlinear dynamic or static case, it could be advisable to use the
deformed configuration (total lagrange!) to create the rigid body modes !
*/
counter=0;
for (i=0; i<actpdis->numnp; i++)
{
   actnode = actpdis->node[i];
   mlprecond.node[counter] = actnode;
   /* get element connected to this node */
   actele = actnode->element[0];
   /* loop nodes on element to find right one */
   for (j=0; j<actele->numnp; j++)
      if (actele->node[j] == actnode) break;
   dsassert(j!=actele->numnp,"Couldn't find node in element");
#ifdef D_SHELL8
   /* h/2 */
   fac = actele->e.s8->thick_node.a.dv[j] * actele->e.s8->sdc / 2.0 ;
   /* averaged director of unit length bischoff style */
   director[counter][0] = actele->e.s8->a3ref.a.da[0][j]*fac;
   director[counter][1] = actele->e.s8->a3ref.a.da[1][j]*fac;
   director[counter][2] = actele->e.s8->a3ref.a.da[2][j]*fac;
#endif
   counter++;
}
dsassert(counter==actpdis->numnp,"Number of nodes wrong");
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_getdirs */

/*! @} (documentation module close)*/
