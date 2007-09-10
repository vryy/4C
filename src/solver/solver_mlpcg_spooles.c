#ifndef CCADISCRET
#ifdef MLPCG

/*!---------------------------------------------------------------------
\file
\brief contains the multilevel preconditioner for shells

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
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

/*!------------------------------------------------------------------------
\brief stuff needed by the spooles solver

m.gee 6/01

matrix needed by the MLPCG solver on the finest grid only

-------------------------------------------------------------------------*/
typedef struct _SPO_DATA
{
#ifdef SPOOLES_PACKAGE
FrontMtx               *frontmtx;
InpMtx                 *mtxA;
DenseMtx               *mtxY;
DenseMtx               *mtxX;
Graph                  *graph;
IVL                    *adjIVL;
ETree                  *frontETree;
IVL                    *symbfacIVL;
SubMtxManager          *mtxmanager;
SubMtxManager          *solvemanager;
ChvManager             *chvmanager;
Chv                    *rootchv;
IV                     *oldToNewIV;
IV                     *newToOldIV;
IV                     *ownersIV;
IV                     *vtxmapIV;
DV                     *cumopsDV;
IV                     *ownedColumnsIV;
InpMtx                 *newA;
DenseMtx               *newY;
SolveMap               *solvemap;
#else
INT i; /* just a dummy to prevent empty structure */
#endif
} SPO_DATA;
static struct _SPO_DATA sdata;


/*!---------------------------------------------------------------------
\brief lapack solver

<pre>                                                        m.gee 11/02

</pre>
\param z            DOUBLE*      (o)   the solution of the solve
\param r            DOUBLE*      (i)   the rhs
\param csr          DBCSR*       (i)   the matrix to be solved with
\param actintra     INTRA*       (i)   the intra-communicator of this field
\return void

------------------------------------------------------------------------*/
void mlpcg_precond_spoolessolve(DOUBLE *z, DOUBLE *r, DBCSR *csr, INTRA *actintra)
{
#ifdef SPOOLES_PACKAGE
INT              i,j;
INT              myrank,nproc;
INT              pivotingflag;
INT             *update,*ia,*ja;
DOUBLE          *a;
INT              numeq,numeq_total,nnz;
DOUBLE           cpus[20];
INT              stats[20];
INT             *rowind1;
INT              nrow;
INT              msglvl = 0;
FILE            *msgFile = NULL;
INT              nedges;
INT              seed = 10101;
DOUBLE          *opcounts,minops,cutoff,tau=100.0;
DOUBLE           droptol=0.0;
INT              root,firsttag=0,error=-1;
INT              sym=SPOOLES_NONSYMMETRIC;
INT              nmycol;
DOUBLE          *recv;
ARRAY            recv_a;
INT              new;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("mlpcg_precond_spoolessolve");
#endif
/*----------------------------------------------------------------------*/
myrank      = actintra->intra_rank;
nproc       = actintra->intra_nprocs;
/*dsmemreport();fflush(stdout);*/
/*----------------------- INT the inprocs=1 case pivoting does not work */
if (nproc>1)   pivotingflag=1;
else           pivotingflag=0;
/*----------------------------------------------------------------------*/
numeq       = csr->numeq;
numeq_total = csr->numeq_total;
update      = csr->update.a.iv;
ia          = csr->ia.a.iv;
ja          = csr->ja.a.iv;
a           = csr->a.a.dv;
nnz         = ia[numeq];
recv        = amdef("recv",&recv_a,numeq_total,1,"DV");
/*----------------------------------------------------------------------*/
IVzero(20,stats);
DVzero(20,cpus);
/*----------------------------------------------------------------------*/
new=0;
if (csr->is_factored==0) new=1;
else if (csr->is_factored != mlprecond.ncall+1 && mlprecond.mod==0) new=1;
/*----------------------------------------------------------------------*/
if (new && mlprecond.ncall != 0)
{
      FrontMtx_free(sdata.frontmtx);
      InpMtx_free(sdata.newA);
      /*DenseMtx_free(sdata.newY);*/
      ETree_free(sdata.frontETree);
      SubMtxManager_free(sdata.mtxmanager);
      IV_free(sdata.newToOldIV);
      IV_free(sdata.oldToNewIV);
      IV_free(sdata.ownersIV);
      IV_free(sdata.vtxmapIV);
      IV_free(sdata.ownedColumnsIV);
      SolveMap_free(sdata.solvemap);
      IVL_free(sdata.symbfacIVL);
}
/*----------------------------------------------------------------------*/
if (new)
{
   sdata.mtxA = InpMtx_new();

   InpMtx_init(sdata.mtxA,INPMTX_BY_ROWS,1,nnz,0);

   for (i=0; i<numeq; i++)
   for (j=ia[i]; j<ia[i+1]; j++)
      InpMtx_inputRealEntry(sdata.mtxA,update[i],ja[j],a[j]);

   InpMtx_sortAndCompress(sdata.mtxA);

   InpMtx_changeStorageMode(sdata.mtxA,INPMTX_BY_VECTORS) ;
}
/*----------------------------------------------------------------------*/
sdata.mtxY = DenseMtx_new();

DenseMtx_init(sdata.mtxY,SPOOLES_REAL,0,0,numeq,1,1,numeq);

DenseMtx_zero(sdata.mtxY);

DenseMtx_rowIndices(sdata.mtxY, &nrow, &rowind1);

for (i=0; i<numeq; i++)
{
   rowind1[i] = update[i];
   DenseMtx_setRealEntry(sdata.mtxY,i,0,r[i]);
}
/*----------------------------------------------------------------------*/
if (new)
{
   sdata.graph  = Graph_new();

   sdata.adjIVL = InpMtx_MPI_fullAdjacency(sdata.mtxA,stats,msglvl,msgFile,actintra->MPI_INTRA_COMM);

   nedges = IVL_tsize(sdata.adjIVL);

   Graph_init2(sdata.graph,0,numeq_total,0,nedges,numeq_total,nedges,sdata.adjIVL,NULL,NULL);

   sdata.frontETree = orderViaMMD(sdata.graph,seed+myrank,msglvl,msgFile);

   Graph_free(sdata.graph);

   opcounts = DVinit(nproc,0.0);

   opcounts[myrank] = ETree_nFactorOps(sdata.frontETree,SPOOLES_REAL,sym);

   MPI_Allgather(&opcounts[myrank],1,MPI_DOUBLE,opcounts,1,MPI_DOUBLE,actintra->MPI_INTRA_COMM);

   minops = DVmin(nproc,opcounts,&root);

   DVfree(opcounts);

   sdata.frontETree = ETree_MPI_Bcast(sdata.frontETree,root,msglvl,msgFile,actintra->MPI_INTRA_COMM);
}
/*----------------------------------------------------------------------*/
if (new)
{
   sdata.oldToNewIV = ETree_oldToNewVtxPerm(sdata.frontETree);

   sdata.newToOldIV = ETree_newToOldVtxPerm(sdata.frontETree);

   ETree_permuteVertices(sdata.frontETree,sdata.oldToNewIV);

   InpMtx_permute(sdata.mtxA,IV_entries(sdata.oldToNewIV),IV_entries(sdata.oldToNewIV));

   if (sym==SPOOLES_SYMMETRIC || sym == SPOOLES_HERMITIAN)
   InpMtx_mapToUpperTriangle(sdata.mtxA);

   InpMtx_changeCoordType(sdata.mtxA,INPMTX_BY_CHEVRONS);

   InpMtx_changeStorageMode(sdata.mtxA,INPMTX_BY_VECTORS);
}
/*----------------------------------------------------------------------*/
DenseMtx_permuteRows(sdata.mtxY,sdata.oldToNewIV);
/*----------------------------------------------------------------------*/
if (new)
{
   cutoff   = 1./(2*nproc);

   sdata.cumopsDV = DV_new();

   DV_init(sdata.cumopsDV,nproc,NULL);

   sdata.ownersIV = ETree_ddMap(sdata.frontETree,SPOOLES_REAL,sym,sdata.cumopsDV,cutoff);

   DV_free(sdata.cumopsDV);

   sdata.vtxmapIV = IV_new();

   IV_init(sdata.vtxmapIV,numeq_total,NULL);

   IVgather(numeq_total,IV_entries(sdata.vtxmapIV),IV_entries(sdata.ownersIV),ETree_vtxToFront(sdata.frontETree));
}
/*----------------------------------------------------------------------*/
firsttag = 0;

if (new)
{
  sdata.newA = InpMtx_MPI_split(sdata.mtxA,sdata.vtxmapIV,stats,msglvl,msgFile,firsttag,actintra->MPI_INTRA_COMM);

  firsttag++;

  InpMtx_free(sdata.mtxA);

  sdata.mtxA = sdata.newA;

  InpMtx_changeStorageMode(sdata.mtxA,INPMTX_BY_VECTORS);
}

sdata.newY = DenseMtx_MPI_splitByRows(sdata.mtxY,sdata.vtxmapIV,stats,msglvl,msgFile,firsttag,actintra->MPI_INTRA_COMM);

DenseMtx_free(sdata.mtxY);

sdata.mtxY = sdata.newY;

firsttag += nproc;
/*----------------------------------------------------------------------*/
if (new)
{
   sdata.symbfacIVL = SymbFac_MPI_initFromInpMtx(sdata.frontETree,sdata.ownersIV,sdata.mtxA,stats,msglvl,msgFile,firsttag,actintra->MPI_INTRA_COMM);

   firsttag += sdata.frontETree->nfront;
}
/*----------------------------------------------------------------------*/
if (new)
{
   sdata.mtxmanager = SubMtxManager_new();

   SubMtxManager_init(sdata.mtxmanager, NO_LOCK, 0);

   sdata.frontmtx = FrontMtx_new();

   FrontMtx_init(sdata.frontmtx,sdata.frontETree,sdata.symbfacIVL,SPOOLES_REAL,sym,FRONTMTX_DENSE_FRONTS,pivotingflag,NO_LOCK,myrank,sdata.ownersIV,sdata.mtxmanager,msglvl,msgFile);
}
/*----------------------------------------------------------------------*/
if (new)
{
   sdata.chvmanager = ChvManager_new();

   ChvManager_init(sdata.chvmanager,NO_LOCK,0);

   sdata.rootchv = FrontMtx_MPI_factorInpMtx(sdata.frontmtx,sdata.mtxA,tau,droptol,sdata.chvmanager,sdata.ownersIV,0,&error,cpus,stats,msglvl,msgFile,firsttag,actintra->MPI_INTRA_COMM);

   ChvManager_free(sdata.chvmanager);

   firsttag += 3*(sdata.frontETree->nfront) + 2;
}
/*----------------------------------------------------------------------*/
if (new)
{
   FrontMtx_MPI_postProcess(sdata.frontmtx,sdata.ownersIV,stats,msglvl,msgFile,firsttag,actintra->MPI_INTRA_COMM);

   firsttag += 5*nproc;
}
/*----------------------------------------------------------------------*/
if (new)
{
   sdata.solvemap = SolveMap_new() ;

   SolveMap_ddMap(sdata.solvemap,sdata.frontmtx->symmetryflag,FrontMtx_upperBlockIVL(sdata.frontmtx),FrontMtx_lowerBlockIVL(sdata.frontmtx),nproc,sdata.ownersIV,FrontMtx_frontTree(sdata.frontmtx),seed,msglvl,msgFile);
}
/*----------------------------------------------------------------------*/
if (new)
{
   FrontMtx_MPI_split(sdata.frontmtx,sdata.solvemap,stats,msglvl,msgFile,firsttag,actintra->MPI_INTRA_COMM);
}
/*----------------------------------------------------------------------*/
if ( FRONTMTX_IS_PIVOTING(sdata.frontmtx) )
{
   IV   *rowmapIV ;

   rowmapIV = FrontMtx_MPI_rowmapIV(sdata.frontmtx,sdata.ownersIV,msglvl,msgFile,actintra->MPI_INTRA_COMM);

   sdata.newY = DenseMtx_MPI_splitByRows(sdata.mtxY,rowmapIV,stats,msglvl,msgFile,firsttag,actintra->MPI_INTRA_COMM);

   DenseMtx_free(sdata.mtxY);

   sdata.mtxY = sdata.newY;

   IV_free(rowmapIV);
}
/*----------------------------------------------------------------------*/
if (new)
   sdata.ownedColumnsIV = FrontMtx_ownedColumnsIV(sdata.frontmtx,myrank,sdata.ownersIV,msglvl,msgFile);

nmycol = IV_size(sdata.ownedColumnsIV);

sdata.mtxX = DenseMtx_new();

if ( nmycol > 0 )
{
   DenseMtx_init(sdata.mtxX,SPOOLES_REAL,0,0,nmycol,1,1,nmycol);

   DenseMtx_rowIndices(sdata.mtxX,&nrow,&rowind1);

   IVcopy(nmycol,rowind1,IV_entries(sdata.ownedColumnsIV));
}
/*----------------------------------------------------------------------*/
sdata.solvemanager = SubMtxManager_new();

SubMtxManager_init(sdata.solvemanager, NO_LOCK, 0);

FrontMtx_MPI_solve(sdata.frontmtx,sdata.mtxX,sdata.mtxY,sdata.solvemanager,sdata.solvemap,cpus,stats,msglvl,msgFile,firsttag,actintra->MPI_INTRA_COMM);

SubMtxManager_free(sdata.solvemanager);

DenseMtx_permuteRows(sdata.mtxX,sdata.newToOldIV);

IV_fill(sdata.vtxmapIV, 0);

firsttag++ ;

sdata.mtxX = DenseMtx_MPI_splitByRows(sdata.mtxX,sdata.vtxmapIV,stats,msglvl,msgFile,firsttag,actintra->MPI_INTRA_COMM);
/*----------------------------------------------------------------------*/
/* put complete solution to recv on proc 0 */
if (myrank==0)
   for (i=0; i<numeq_total; i++)
      recv[i] =  sdata.mtxX->entries[i];

/* broadcast solution */
MPI_Bcast(recv,numeq_total,MPI_DOUBLE,0,actintra->MPI_INTRA_COMM);

/* every proc puts his own piece to sol */
for (i=0; i<numeq; i++)
   z[i] = recv[update[i]];

/*----------------------------------------------------------------------*/
csr->is_factored=mlprecond.ncall+1;
DenseMtx_free(sdata.mtxY);
DenseMtx_free(sdata.mtxX);
amdel(&recv_a);
/*----------------------------------------------------------------------*/
exit:
#ifdef DEBUG
dstrc_exit();
#endif
#endif /* end of ifdef SPOOLES_PACKAGE */
return;
} /* end of mlpcg_precond_spoolessolve */


/*! @} (documentation module close)*/
#endif
#endif
