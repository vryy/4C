/*!----------------------------------------------------------------------
\file linalg_solver_spooles.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "linalg_solver.H"

extern "C" 
{
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
}

/*----------------------------------------------------------------------*
 |  solve (protected)                                        mwgee 02/07|
 *----------------------------------------------------------------------*/
void LINALG::Solver::Solve_spooles(const bool reset)
{
#ifdef PARALLEL
#ifdef SPOOLES_PACKAGE
  if (reset)
  {
    if (frontmtx_)       FrontMtx_free(frontmtx_);        frontmtx_      =NULL;
    if (newA_)           InpMtx_free(newA_);              newA_          =NULL;
    if (newY_)           DenseMtx_free(newY_);            newY_          =NULL;
    if (frontETree_)     ETree_free(frontETree_);         frontETree_    =NULL;
    if (mtxmanager_)     SubMtxManager_free(mtxmanager_); mtxmanager_    =NULL;
    if (newToOldIV_)     IV_free(newToOldIV_);            newToOldIV_    =NULL;
    if (oldToNewIV_)     IV_free(oldToNewIV_);            oldToNewIV_    =NULL;
    if (ownersIV_)       IV_free(ownersIV_);              ownersIV_      =NULL;
    if (vtxmapIV_)       IV_free(vtxmapIV_);              vtxmapIV_      =NULL;
    if (ownedColumnsIV_) IV_free(ownedColumnsIV_);        ownedColumnsIV_=NULL;
    if (solvemap_)       SolveMap_free(solvemap_);        solvemap_      =NULL;
    if (graph_)          Graph_free(graph_);              graph_         =NULL;
    if (mtxY_)           DenseMtx_free(mtxY_);            mtxY_          =NULL;
    if (mtxX_)           DenseMtx_free(mtxX_);            mtxX_          =NULL;
    if (mtxA_)           InpMtx_free(mtxA_);              mtxA_          =NULL;
    if (symbfacIVL_)     IVL_free(symbfacIVL_);           symbfacIVL_    =NULL;
  }
  

  // set some values
  int            nnz         = A_->NumGlobalNonzeros();
  int            numeq       = A_->NumMyRows();
  int            numeq_total = A_->NumGlobalRows();  
  int            nrow        = 0;
  int*           rowind1;
  int            root,firsttag=0,error=-1,stats[20];
  IVzero(20,stats);
  int            msglvl=0;
  FILE           *msgFile=NULL;
  const          Epetra_MpiComm& comm = dynamic_cast<const Epetra_MpiComm&>(A_->Comm());
  int            nedges;
  double         *opcounts,minops,cutoff,cpus[20],tau=100.;
  DVzero(20,cpus);
  int            seed = 10101; 
  int            sym=SPOOLES_NONSYMMETRIC;
  double         droptol=0.0;
  vector<double> recv(numeq_total);

/*
   --------------------------------------------
   STEP 1: read the entries
           and create the InpMtx object
   --------------------------------------------
*/
  if (!IsFactored())
  {
    mtxA_ = InpMtx_new();
    InpMtx_init(mtxA_,INPMTX_BY_ROWS,1,nnz,0);
    for (int i=0; i<A_->NumMyRows(); ++i)
    {
      int     numentries;
      double* values;
      int*    indices;
      int err = A_->ExtractMyRowView(i,numentries,values,indices);
      if (err) dserror("Epetra_CrsMatrix::ExtractMyRowView returned an error");
      int I = A_->RowMap().GID(i);
      for (int j=0; j<numentries; ++j)
      {
        int J = A_->ColMap().GID(indices[j]);
        InpMtx_inputRealEntry(mtxA_,I,J,values[j]);
      }
    }
    InpMtx_sortAndCompress(mtxA_);
    InpMtx_changeStorageMode(mtxA_,INPMTX_BY_VECTORS) ;
  }
  
/*
   -----------------------------------------
   STEP 2: read the right hand side matrix B
   -----------------------------------------
*/
  if (mtxY_) DenseMtx_free(mtxY_); 
  mtxY_ = DenseMtx_new();
  DenseMtx_init(mtxY_,SPOOLES_REAL,0,0,numeq,1,1,numeq);
  DenseMtx_zero(mtxY_);
  DenseMtx_rowIndices(mtxY_,&nrow,&rowind1);
  for (int i=0; i<numeq; i++)
  {
    rowind1[i] = b_->Map().GID(i);
    DenseMtx_setRealEntry(mtxY_,i,0,(*b_)[i]);
  }

/*
   -------------------------------------------------
   STEP 3 : find a low-fill ordering
   (1) create the Graph object
   (2) order the graph using multiple minimum degree
   -------------------------------------------------
*/
  if (!IsFactored())
  {
    graph_ = Graph_new();
    _IVL* adjIVL = InpMtx_MPI_fullAdjacency(mtxA_,stats,msglvl,msgFile,comm.GetMpiComm());
    nedges = IVL_tsize(adjIVL);
    Graph_init2(graph_,0,numeq_total,0,nedges,numeq_total,nedges,adjIVL,NULL,NULL);
    
    frontETree_ = orderViaMMD(graph_,seed+comm.MyPID(),msglvl,msgFile);
    Graph_free(graph_);
    graph_ = NULL;
    
    opcounts               = DVinit(comm.NumProc(),0.0);
    opcounts[comm.MyPID()] = ETree_nFactorOps(frontETree_,SPOOLES_REAL,sym);
    MPI_Allgather(&opcounts[comm.MyPID()],1,MPI_DOUBLE,
                  opcounts,1,MPI_DOUBLE,comm.GetMpiComm());
    minops = DVmin(comm.NumProc(),opcounts,&root);
    DVfree(opcounts);
    frontETree_ = ETree_MPI_Bcast(frontETree_,root,msglvl,msgFile,comm.GetMpiComm());
  }

/*
   -----------------------------------------------------
   STEP 4: get the permutation, permute the matrix and
           front tree and get the symbolic factorization,
           permute the right hand side.
   -----------------------------------------------------
*/
  if (!IsFactored())
  {
    oldToNewIV_ = ETree_oldToNewVtxPerm(frontETree_);
    newToOldIV_ = ETree_newToOldVtxPerm(frontETree_);
    ETree_permuteVertices(frontETree_,oldToNewIV_);
    InpMtx_permute(mtxA_,IV_entries(oldToNewIV_),IV_entries(oldToNewIV_));
    InpMtx_changeCoordType(mtxA_,INPMTX_BY_CHEVRONS);
    InpMtx_changeStorageMode(mtxA_,INPMTX_BY_VECTORS);
  }
  
/*
   -------------------------------------------
   STEP 5: generate the owners map IV object
           and the map from vertices to owners
   -------------------------------------------
*/
  if (!IsFactored())
  {
    cutoff    = 1./(2*comm.NumProc());
    cumopsDV_ = DV_new();
    DV_init(cumopsDV_,comm.NumProc(),NULL);
    ownersIV_ = ETree_ddMap(frontETree_,SPOOLES_REAL,sym,cumopsDV_,cutoff);
    DV_free(cumopsDV_);
    cumopsDV_ = NULL;
    vtxmapIV_ = IV_new();
    IV_init(vtxmapIV_,numeq_total,NULL);
    IVgather(numeq_total,IV_entries(vtxmapIV_),IV_entries(ownersIV_),
             ETree_vtxToFront(frontETree_));
  }
  
/*
   ---------------------------------------------------
   STEP 6: redistribute the matrix and right hand side
   ---------------------------------------------------
*/
  if (!IsFactored())
  {
    newA_ = InpMtx_MPI_split(mtxA_,vtxmapIV_,stats,msglvl,msgFile,
                             firsttag,comm.GetMpiComm());
    firsttag++;
    InpMtx_free(mtxA_);
    mtxA_ = newA_;
    newA_ = NULL;
    InpMtx_changeStorageMode(mtxA_,INPMTX_BY_VECTORS);
  }
  newY_ = DenseMtx_MPI_splitByRows(mtxY_,vtxmapIV_,stats,msglvl,
                                   msgFile,firsttag,comm.GetMpiComm());
  DenseMtx_free(mtxY_);
  mtxY_ = newY_;
  newY_ = NULL;
  firsttag += comm.NumProc();
  
/*
   ------------------------------------------
   STEP 7: compute the symbolic factorization
   ------------------------------------------
*/
  if (!IsFactored())
  {
    symbfacIVL_ = SymbFac_MPI_initFromInpMtx(frontETree_,ownersIV_,mtxA_,stats,
                                             msglvl,msgFile,firsttag,
                                             comm.GetMpiComm());
    firsttag += frontETree_->nfront;
    
  }
  
/*
   -----------------------------------
   STEP 8: initialize the front matrix
   -----------------------------------
*/
  int pivotingflag;
  if (comm.NumProc()>1) pivotingflag = 1;
  else                  pivotingflag = 0;
  if (!IsFactored())
  {
    mtxmanager_ = SubMtxManager_new();
    SubMtxManager_init(mtxmanager_,NO_LOCK,0);
    frontmtx_ = FrontMtx_new();
    FrontMtx_init(frontmtx_,frontETree_,symbfacIVL_,SPOOLES_REAL,sym,
                  FRONTMTX_DENSE_FRONTS,pivotingflag,NO_LOCK,comm.MyPID(),
                  ownersIV_,mtxmanager_,msglvl,msgFile);
  }

/*
   -----------------------------------------
   STEP 9: compute the numeric factorization
   -----------------------------------------
*/
  if (!IsFactored())
  {
    chvmanager_ = ChvManager_new();
    ChvManager_init(chvmanager_,NO_LOCK,0);
    FrontMtx_MPI_factorInpMtx(frontmtx_,mtxA_,tau,droptol,chvmanager_,ownersIV_,
                              0,&error,cpus,stats,msglvl,msgFile,firsttag,
                              comm.GetMpiComm());
    ChvManager_free(chvmanager_);
    chvmanager_ = NULL;
    firsttag += 3*(frontETree_->nfront) + 2;
   if ( error >= 0 )
      dserror("Error in spooles numeric factorization");
  }
  
/*
   --------------------------------------
   STEP 10: post-process the factorization
   --------------------------------------
*/
  if (!IsFactored())
  {
    FrontMtx_MPI_postProcess(frontmtx_,ownersIV_,stats,msglvl,msgFile,firsttag,
                             comm.GetMpiComm());
    firsttag += 5*comm.NumProc();
  }  

/*
   -----------------------------------
   STEP 11: create the solve map object
   -----------------------------------
*/
  if (!IsFactored())
  {
    solvemap_ = SolveMap_new();
    SolveMap_ddMap(solvemap_,frontmtx_->symmetryflag,
                   FrontMtx_upperBlockIVL(frontmtx_),
                   FrontMtx_lowerBlockIVL(frontmtx_),
                   comm.NumProc(),ownersIV_,FrontMtx_frontTree(frontmtx_),
                   seed,msglvl,msgFile);
  }

/*
   ----------------------------------------------------
   STEP 12: redistribute the submatrices of the factors
   ----------------------------------------------------
*/
  if (!IsFactored())
    FrontMtx_MPI_split(frontmtx_,solvemap_,stats,msglvl,msgFile,firsttag,
                       comm.GetMpiComm());
  
/*
   ------------------------------------------------
   STEP 13: permute and redistribute Y if necessary
   ------------------------------------------------
*/
  if ( FRONTMTX_IS_PIVOTING(frontmtx_))
  {
    _IV *rowmapIV;
    rowmapIV = FrontMtx_MPI_rowmapIV(frontmtx_,ownersIV_,msglvl,msgFile,
                                     comm.GetMpiComm());
    newY_ = DenseMtx_MPI_splitByRows(mtxY_,rowmapIV,stats,msglvl,msgFile,
                                     firsttag,comm.GetMpiComm());
    DenseMtx_free(mtxY_);
    mtxY_ = newY_;
    newY_ = NULL;
    IV_free(rowmapIV);
  }

/*
   ------------------------------------------
   STEP 14: create a solution DenseMtx object
   ------------------------------------------
*/
  if (!IsFactored())
  {
    ownedColumnsIV_ = FrontMtx_ownedColumnsIV(frontmtx_,comm.MyPID(),ownersIV_,
                                              msglvl,msgFile);
  }
  int nmycol = IV_size(ownedColumnsIV_);
  mtxX_  = DenseMtx_new();
  if (nmycol > 0)
  {
    DenseMtx_init(mtxX_,SPOOLES_REAL,0,0,nmycol,1,1,nmycol);
    DenseMtx_rowIndices(mtxX_,&nrow,&rowind1);
    IVcopy(nmycol,rowind1,IV_entries(ownedColumnsIV_));
  }

/*
   --------------------------------
   STEP 15: solve the linear system
   --------------------------------
*/
  _SubMtxManager *solvemanager = SubMtxManager_new();
  SubMtxManager_init(solvemanager, NO_LOCK, 0);
  FrontMtx_MPI_solve(frontmtx_,mtxX_,mtxY_,solvemanager,
                     solvemap_,cpus,stats,msglvl,msgFile,firsttag,
                     comm.GetMpiComm());
  SubMtxManager_free(solvemanager);

/*
   --------------------------------------------------------
   STEP 16: permute the solution into the original ordering
            and assemble the solution onto processor zero
   --------------------------------------------------------
*/
  DenseMtx_permuteRows(mtxX_,newToOldIV_);
  IV_fill(vtxmapIV_,0);
  firsttag++;
  {
    _DenseMtx* mtxX;
    mtxX = DenseMtx_MPI_splitByRows(mtxX_,vtxmapIV_,stats,msglvl,msgFile,
                                    firsttag,comm.GetMpiComm());
    DenseMtx_free(mtxX_);
    mtxX_ = mtxX;
  }
  if (comm.MyPID()==0)
    for (int i=0; i<numeq_total; ++i)
      recv[i] = mtxX_->entries[i];
  MPI_Bcast(&recv[0],numeq_total,MPI_DOUBLE,0,comm.GetMpiComm());
  for (int i=0; i<numeq; ++i)
    (*x_)[i] = recv[x_->Map().GID(i)];
  
/*
   -----------
   free memory
   -----------
*/
  DenseMtx_free(mtxX_);
  mtxX_ = NULL;

#endif
#endif
  return;
}

#endif  // #ifdef CCADISCRET
