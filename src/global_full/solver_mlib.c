#include "../headers/standardtypes.h"
#include "../headers/solution.h"
#ifdef MLIB_PACKAGE
#include "/opt/mlib/include/veclib.h" 
#endif
/* prototypes */
/*----------------------------------------------------------------------*
 |  control solver lib MLIB                               al   10/01    |
 *----------------------------------------------------------------------*/
void solver_mlib( 
                      struct _SOLVAR         *actsolv,
                      struct _INTRA          *actintra,
                      struct _ML_ARRAY_MDS   *mds,
                      struct _DIST_VECTOR    *sol,
                      struct _DIST_VECTOR    *rhs,
                      int                     option
                     )
{
/*-----------------------------------------------------------------------|
| option = 1: initialize sparse matrix package                           |
| option = 2: initialize sparse matrix with zero, keep matrix structure  |
| option = 0:  calculation phase                                         |
|-----------------------------------------------------------------------*/
#ifdef MLIB_PACKAGE
  int i,j,k;
  int symm;
  char order[4];
  char usymm[3];
  double      *vz;
  double      *vzh;
  MLVAR        *mlvar;
  #ifdef DEBUG 
  dstrc_enter("solver_mlib");
  #endif
/*----------------------------------------------------------------------*/
  mlvar = actsolv->mlvar;
  symm = mlvar->symm;
/*----------------------------------------------------------------------*/
switch(option)
{ 
case 1:/*========================= INITIALIZE THE SPARSE MATRIX PACKAGE */
  
  mds->output = 6;
  mds->ierr   = 0;
  
  dslein (&mds->numeq,&mlvar->msglvl,&mds->output,mds->global,&mds->ierr);

  usymm[0]='S';
  usymm[1]='U';
  
  if(!symm) dslema (usymm,mds->global,&mds->ierr,2);
/*--------------------------------------- INPUT THE MATRIX STRUCTURE ---*/
  dsleim(&mds->colstr.a.iv[0],&mds->rowind.a.iv[0],mds->global,&mds->ierr);
/*----------------------------------------------- REORDER THE MATRIX ---*/
/*-------------------------------- optional, kann man auch weglassen ---*/
switch(mlvar->order)
{
/* default by DSLEIN: use multiple minimum degree ordering */
/*     ORDER='MMD'
/* Use natural ordering */
/*     ORDER='NAT'
/* use constrained minimum degree ordering */
/*      ORDER='CMD'
/* geht nicht:
/*      ORDER='RCM'
/*      ORDER='1WD'
/*      ORDER='GND'
/*      ORDER='MET'
/* */
case 0:
break;
case 1:
   order[0]='M';
   order[1]='M';
   order[2]='D';
break;
case 2:
   order[0]='N';
   order[1]='A';
   order[2]='T';
break;
case 3:
   order[0]='C';
   order[1]='M';
   order[2]='D';
break;
default:
   dserror("Unknown typ of ordering");
break;   
}
if(mlvar->order>0)   dsleop (order, mds->global, &mds->ierr,mds->numeq);
/*      IF ( IER .NE. 0 ) GO TO 8000*/

  dsleor (&mlvar->maxzer,mds->global,&mds->ierr);

break;
case 2:/*======================= INITIALIZE THE SPARSE MATRIX WITH ZERO */
  vz = (double*)calloc(mds->nnz ,sizeof(double));
  if (!vz)  dserror("Allocation of memory int 'add_mds' failed");
  vzh = vz;
  for (i=0; i<mds->nnz; i++) *(vzh++) = 0.0;
  
  dslevm (&mds->colstr.a.iv[0],&mds->rowind.a.iv[0]  ,
          vz                  ,mds->global,&mds->ierr);
          
  free(vz);
break;
case 0:/*============================================ calculation phase */
/**/
/*-------------- FACTOR THE MATRIX AND ESTIMATE ITS CONDITION NUMBER ---*/
  if(!symm) dslefa (&mlvar->pvttol, mds->inrtia, mds->global, &mds->ierr);
  else
  dsleco (&mlvar->pvttol, &mds->cond, mds->inrtia,mds->global,&mds->ierr);
/*-- print additional information, depends on the stage of execution ---*/
  if(mlvar->msglvl==4) dsleps (mds->global);
/*-------------------------------- SOLVE FOR A GIVEN RIGHT HAND SIDE ---*/
  dslesl(&actsolv->nrhs, &(rhs->vec.a.dv[0]), &mds->numeq,
                                                 mds->global, &mds->ierr);
  if(mds->ierr!=0) exit; 
/*-- print additional information, depends on the stage of execution ---*/
   if(mlvar->msglvl==4) dsleps (mds->global);
/*----------------------------- USE THE SOLUTION STORED IN ARRAY RHS ---*/
   for (i=0; i<rhs->numeq; i++)
   {
      sol->vec.a.dv[i] = rhs->vec.a.dv[i];
   }
break;/*=============================================================== */
default:
   dserror("Unknown option for solver call to hp's mlib");
break;   
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#endif /* end of ifdef MLIB_PACKAGE */
return;
} /* end of solver_mlib */




