/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/prototypes_mlpcg.h"
#include "../shell8/shell8.h"
/*!----------------------------------------------------------------------
\brief the multilevel preconditioner main structure

<pre>                                                         m.gee 09/02    
defined in solver_mlpcg.c
</pre>

*----------------------------------------------------------------------*/
extern struct _MLPRECOND mlprecond;



#if 0
int dense_factor(DBCSR *asm, DBCSR *ilu)
{
int i,j,k,l,n,ii,jj;
int numeq;
int *update,*ia,*ja;
int actrow,actcol,colstart,colend;
double *a;
int    info;

double **dense;
int     *ipiv;

double   diff;

numeq  = asm->numeq;
update = asm->update.a.iv;
ia     = asm->ia.a.iv;
ja     = asm->ja.a.iv;
a      = asm->a.a.dv;

/* allocate dense matrix in ilu matrix */
ilu->dense = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
dense      = amdef("dense",ilu->dense,numeq,numeq,"DA");
amzero(ilu->dense);
/* allocate ipiv */
ilu->ipiv  = (ARRAY*)CCACALLOC(1,sizeof(ARRAY));
ipiv       = amdef("ipiv",ilu->ipiv,numeq,1,"IV");
amzero(ilu->ipiv);


/* copy the matrix to the dense matrix */
for (i=0; i<numeq; i++)
{
   actrow   = i;
   colstart = ia[i];
   colend   = ia[i+1];
   for (j=colstart; j<colend; j++)
   {
      actcol                = ja[j];
      dense[actrow][actcol] = a[j];
   }
}
/* check for syymetry of the matrix */
for (i=0; i<numeq; i++)
{
   for (j=0; j<numeq; j++)
   {
      diff = FABS(dense[i][j]-dense[j][i]);
      if (diff>EPS14)
      {
         printf("i %d j %d diff = %30.20f\n",i,j,diff);
      }
   }
}
/* factor the matrix */
info = 1;
dgetrf(&numeq,&numeq,dense[0],&numeq,ipiv,&info);
if (info != 0) dserror("Lapack dgetrf returned info nonzero");
return;
} 


int dense_solve(DBCSR *ilu, double *r, double *z)
{
int i,j;
int numeq;
int    info=1;

double **dense;
int     *ipiv;
char    trans[1];
double *b;
int     ione=1;
trans[0] = 'N';
numeq  = ilu->numeq;
dense  = ilu->dense->a.da;
ipiv   = ilu->ipiv->a.iv;

b = (double*)CCAMALLOC(numeq*sizeof(double));
for (i=0; i<numeq; i++) b[i] = r[i];

dgetrs(trans,&numeq,&ione,dense[0],&numeq,ipiv,b,&numeq,&info);
if (info != 0) dserror("Lapack dgetrf returned info nonzero");

for (i=0; i<numeq; i++) z[i] = b[i];
CCAFREE(b);

return;
}
#endif

#if 0
int mlpcg_precond_oneP_vanek_test(AGG     *actagg,
                             double   aggblock[][500],
                             int      rindex[],
                             int      cindex[],
                             int     *nrow,
                             int     *ncol,
                             DBCSR   *actstiff,
                             MLLEVEL *prevlevel)
{
int           i,j,k,l,counter;
AGG          *prevagg[200];
AGG          *actprevagg;
NODE         *node[200];
NODE         *actnode;
PARTDISCRET  *actpdis;
int          *actblock;
int           dof;
double        x0,y0,z0,x,y,z,a1,a2,a3;
int           index;
int           foundit;
int           shift,bins[5000];
double      **R;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("mlpcg_precond_oneP_vanek");
#endif
/*----------------------------------------------------------------------*/
/*------------------------ get the size of the tentative prolong. block */
*nrow=0;
for (i=0; i<actagg->nblock; i++)
   (*nrow) += actagg->block[i][0];
   
*ncol = actagg->numdf;
if (*nrow >= 1000 || *ncol >= 500)
dserror("Local variable aggblock[1000][500] too small");
/*--------------------------------------- get the global column indizes */
for (i=0; i<actagg->numdf; i++)
   cindex[i] = actagg->dof[i];
/*------------------------------------------ get the global row indizes */   
counter=0;
for (i=0; i<actagg->nblock; i++)
   for (j=0; j<actagg->block[i][0]; j++)
   {
      rindex[counter] = actagg->block[i][j+1];
      counter++;
   }
dsassert(counter==*nrow,"Number of dofs in prolongator wrong");
qsort((int*)rindex,*nrow,sizeof(int),cmp_int);
/*------------------------------ get the nodal patch from the partition */
/*======================================================================*/
/*======================================================================*/
#if 1 /*------------------ this is with R factors from rigid body modes */
/*======================================================================*/
/*======================================================================*/
dsassert(actagg->nblock<=200,"Local variable prevagg[200] too small");
dsassert(*ncol==6,"number of rbm's has to be 6 in this case");
/* 
the index of the agg in prevlev corresponds to the index of the block in
blocks
*/
counter=0;
for (i=0; i<actagg->nblock; i++)
{
   actblock = actagg->block[i];
   dof      = actblock[1];
   for (j=0; j<prevlevel->nagg; j++)
   {
       actprevagg = &(prevlevel->agg[j]);
       for (k=0; k<actprevagg->numdf; k++)
       {
          if (dof == actprevagg->dof[k])
          {
             prevagg[counter] = actprevagg;
             counter++;
             goto nextblock;
          }
       }
   }
   nextblock:;
}
dsassert(counter==actagg->nblock,"Cannot find aggregates on lower level");
/*----------------------- loop the nodes and put values to the aggblock */
if (*nrow > 18000) dserror("local bins too small");
init_quick_find(rindex,*nrow,&shift,bins);
for (i=0; i<actagg->nblock; i++)
{
   /* these are the aggregates on prevlev which are part of actagg */
   actprevagg = prevagg[i];
   R          = actprevagg->R->a.da;
   /* loop the dofs of the previous aggregate */
   for (j=0; j<actprevagg->numdf; j++)
   {
      dof = actprevagg->dof[j];
      /* find the dof in the row index */
      index = quick_find(dof,rindex,*nrow,shift,bins);
      if (index==-1) dserror("Cannot find aggregate-local dof");
      /* put a row of R to the aggblock */
      for (k=0; k<actprevagg->numdf; k++)
         aggblock[index][k] = R[j][k];
   }/* end of for (j=0; j<actprevagg->numdf; j++) */
} /* end of for (i=0; i<actagg->nblock; i++) */
/*======================================================================*/
/*======================================================================*/
#endif /*- this is the calculation of the rigid body mode for SHELL8 !!!*/
/*======================================================================*/
/*======================================================================*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of mlpcg_precond_oneP_vanek */
#endif
