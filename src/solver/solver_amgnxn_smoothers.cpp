/*!----------------------------------------------------------------------
\file solver_amgnxn_smoothers.cpp

<pre>
Maintainer: Francesc Verdugo
            verdugo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Feb 27, 2014
</pre>
*----------------------------------------------------------------------*/


#ifdef HAVE_MueLu

#include <iostream>

#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <MueLu_MLParameterListInterpreter_decl.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_EpetraOperator.hpp>
#include "EpetraExt_RowMatrixOut.h"
#include "../drt_lib/drt_dserror.H"
#include "solver_amgnxn_smoothers.H"
#include "solver_amgnxn_hierarchies.H"
#include "solver_amgnxn_preconditioner.H"


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/


Teuchos::RCP<LINALG::SparseOperator> LINALG::Multiply(
    Teuchos::RCP<SparseOperator> A, bool transA,
    Teuchos::RCP<SparseOperator> B, bool transB,
    bool complete)
{

  // This is only a prototype. Things can be done better.

  // If both matrices are sparse matrices, perform the standard multiply.
  // Then cast the result to SparseOperator. return.

  // If at least one is blocked
  // 1) If one is none blocked, wrap it into a blocked
  // 2) Perform standard multiplication between two blocked matrices
  // 3) If the result is 1x1, recover the SparseMatrix and cast to SparseOperator. return.
  //    Else. Cast the result to SparseOperator. return.

  Teuchos::RCP<SparseOperator> C = Teuchos::null;

  Teuchos::RCP<BlockSparseMatrixBase> Abl = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(A);
  Teuchos::RCP<BlockSparseMatrixBase> Bbl = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(B);
  Teuchos::RCP<SparseMatrix> Asp = Teuchos::rcp_dynamic_cast<SparseMatrix>(A);
  Teuchos::RCP<SparseMatrix> Bsp = Teuchos::rcp_dynamic_cast<SparseMatrix>(B);

  bool A_is_sp = (Abl==Teuchos::null);
  bool B_is_sp = (Bbl==Teuchos::null);

  if ( A_is_sp and B_is_sp )
  {
    if (Asp==Teuchos::null or Bsp==Teuchos::null)
      dserror("Something went wrong!");
    Teuchos::RCP<LINALG::SparseMatrix> Csp = LINALG::Multiply(*Asp,transA,*Bsp,transB,complete);
    C = Teuchos::rcp_dynamic_cast<SparseOperator>(Csp);
  }
  else
  {
    SOLVER::BlockSparseMatrix_Creator myBlockMatrixCreator;
    if(A_is_sp)
    {
      if(Asp==Teuchos::null)
        dserror("Something went wrong");
      std::vector<Teuchos::RCP<SparseMatrix> > blocks(1,Asp);
      Abl = myBlockMatrixCreator.CreateBlockSparseMatrix(blocks,1,1,View);
    }
    if(B_is_sp)
    {
      if(Bsp==Teuchos::null)
        dserror("Something went wrong");
      std::vector<Teuchos::RCP<SparseMatrix> > blocks(1,Bsp);
      Bbl = myBlockMatrixCreator.CreateBlockSparseMatrix(blocks,1,1,View);
    }
    Teuchos::RCP<BlockSparseMatrixBase> Cbl = Multiply(*Abl,transA,*Bbl,transB,true,false,complete);
    if( (Cbl->Rows()==1) and (Cbl->Cols()==1) )
    {
      const SparseMatrix& Caux = Cbl->Matrix(0,0);
      Teuchos::RCP<SparseMatrix> Csp = Teuchos::rcp(new SparseMatrix(Caux,Copy));
      C = Teuchos::rcp_dynamic_cast<SparseOperator>(Csp);
    }
    else
    {
      C = Teuchos::rcp_dynamic_cast<SparseOperator>(Cbl);
    }
  }


  return C;
}



/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::BlockSparseMatrixBase>
  LINALG::SOLVER::BlockSparseMatrix_Creator::CreateBlockSparseMatrix(
  std::vector< Teuchos::RCP<SparseMatrix> > blocks,
  int rows,
  int cols,
  Epetra_DataAccess access,
  bool explicitdirichlet,
  bool savegraph)
{

  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::BlockSparseMatrix_Creator::CreateBlockSparseMatrix");

  // Check if the number of given of blocks is consistent with rows and cols
  int NumBlocks = blocks.size();
  bool flag_all_blocks_are_given = false;
  if(NumBlocks==rows*cols)
    flag_all_blocks_are_given = true;
  else if(NumBlocks==rows and NumBlocks==cols)
    flag_all_blocks_are_given = false;
  else
    dserror("The number of given blocks is not consistent with the given number of block rows and columns");

  // Determine the estimated number of non zero entries per row
  int npr = 0;
  //for(int i=0; i<(int)blocks.size();i++)
  //{
  //  if(blocks[i]==Teuchos::null)
  //    dserror("The given blocks cannot be null pointers");
  //  if(blocks[i]->MaxNumEntries()>npr)
  //    npr = blocks[i]->MaxNumEntries();
  //}

  // Some checks
  if(flag_all_blocks_are_given)
  {
    // check that all the rows have the same range map
    for(int row=0;row<rows;row++)
      for(int col=1;col<cols;col++)
        if(!((blocks[row*cols+col]->RangeMap()).SameAs(blocks[row*cols+0]->RangeMap())))
          dserror("The range map must be the same for all blocks in the same row");
    // check that all the cols have the same domain map
    for(int col=0;col<cols;col++)
      for(int row=1;row<rows;row++)
        if(!((blocks[row*cols+col]->DomainMap()).SameAs(blocks[0*cols+col]->DomainMap())))
          dserror("The domain map must be the same for all blocks in the same col");
  }

  // build the partial and full domain maps
  std::vector< Teuchos::RCP<const Epetra_Map> > domain_maps(cols,Teuchos::null);
  for(int i=0;i<cols;i++)
    if(flag_all_blocks_are_given)
      domain_maps[i]= Teuchos::rcp(new Epetra_Map(blocks[0*cols+i]->DomainMap()));
  // we assume the rest of rows are consistent.
    else
      domain_maps[i]= Teuchos::rcp(new Epetra_Map(blocks[i]->DomainMap()));
  Teuchos::RCP<Epetra_Map> fullmap_domain = MultiMapExtractor::MergeMaps(domain_maps);
  Teuchos::RCP<MultiMapExtractor> domainmaps  =
    Teuchos::rcp(new MultiMapExtractor(*fullmap_domain ,domain_maps ));

  // build the partial and full range maps
  std::vector< Teuchos::RCP<const Epetra_Map> > range_maps (rows,Teuchos::null);
  for(int i=0;i<rows;i++)
    if(flag_all_blocks_are_given)
      range_maps[i] = Teuchos::rcp(new Epetra_Map(blocks[i*cols+0]->RangeMap()));
  // we assume the rest of cols are consistent.
    else
      range_maps[i] = Teuchos::rcp(new Epetra_Map(blocks[i]->RangeMap()));
  Teuchos::RCP<Epetra_Map> fullmap_range  = MultiMapExtractor::MergeMaps(range_maps);
  Teuchos::RCP<MultiMapExtractor> rangemaps  =
    Teuchos::rcp(new MultiMapExtractor(*fullmap_range ,range_maps ));

  // Create the concrete matrix
  Teuchos::RCP< LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > the_matrix =
    Teuchos::rcp( new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *domainmaps,
          *rangemaps,
          npr,
          explicitdirichlet,
          savegraph));

  // Assign the blocks
  if(flag_all_blocks_are_given)
  {
    for(int i=0;i<rows;i++)
      for(int j=0;j<cols;j++)
        the_matrix->Assign(i,j,access,*(blocks[i*cols+j]));
  }
  else
  {
    for(int i=0;i<rows;i++)
      the_matrix->Assign(i,i,access,*(blocks[i]));
    // Do not forget to zero out the off-diagonal blocks!!!
    for(int i=0;i<rows;i++)
      for(int j=0;j<cols;j++)
        if(i!=j)
        {
          the_matrix->Matrix(i,j).Zero();
          the_matrix->Matrix(i,j).Scale(0.0);
        }
  }

  // Call complete
  the_matrix->Complete();

  // Return
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> the_matrix_base =
    Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(the_matrix);

  return the_matrix_base;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

LINALG::SOLVER::BlockAggrupator::BlockAggrupator( Teuchos::RCP<BlockSparseMatrixBase> A,
    std::vector<std::vector<int> > superblocks_to_blocks):
  superblocks_to_blocks_(superblocks_to_blocks),
  A_(A),
  num_superblocks_(superblocks_to_blocks.size()) { Setup();}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector< Teuchos::RCP<LINALG::SparseOperator> >
LINALG::SOLVER::BlockAggrupator::GetSuperBlocksRowOrder()
{return superblocks_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::MultiMapExtractor> LINALG::SOLVER::BlockAggrupator::GetRangeMapExtractor()
{return range_ex_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::MultiMapExtractor> LINALG::SOLVER::BlockAggrupator::GetDomainMapExtractor()
{return domain_ex_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::BlockAggrupator::GetNumSuperBlocks(){return num_superblocks_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<std::vector<int> >
LINALG::SOLVER::BlockAggrupator::GetSuperBlocksToBlocks(){return superblocks_to_blocks_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SparseOperator>
LINALG::SOLVER::BlockAggrupator::GetSuperBlock(int srow, int scol)
{return superblocks_[srow*GetNumSuperBlocks()+scol];}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::BlockAggrupator::Setup()
{

  // =============================================================
  // Construct  the operators for the superblocks
  // =============================================================

  BlockSparseMatrix_Creator myBlockMatrixCreator;

  int num_all_blocks = (GetNumSuperBlocks())*(GetNumSuperBlocks());
  superblocks_.assign(num_all_blocks,Teuchos::null);
  for(int srow=0;srow<GetNumSuperBlocks();srow++)
  {
    std::vector<int> Rows = GetSuperBlocksToBlocks()[srow];
    int NumRows = Rows.size();
    for(int scol=0;scol<GetNumSuperBlocks();scol++)
    {
      std::vector<int> Cols = GetSuperBlocksToBlocks()[scol];
      int NumCols = Cols.size();
      Teuchos::RCP<SparseOperator> SuperBlock = Teuchos::null;
      if (NumRows == 1 and NumCols == 1)
      {
        // TODO does "false" avoids segfaults?
        SuperBlock = Teuchos::rcp(&(A_->Matrix(Rows[0],Cols[0])),false);
      }
      else
      {
        std::vector< Teuchos::RCP<SparseMatrix> > Blocks(NumRows*NumCols,Teuchos::null);
        for(int r=0;r<NumRows;r++)
          for(int c=0;c<NumCols;c++)
          {
            Blocks[r*NumCols+c] = Teuchos::rcp(&(A_->Matrix(Rows[r],Cols[c])),false);
          }
        //TODO view?
        SuperBlock = myBlockMatrixCreator.CreateBlockSparseMatrix(Blocks,NumRows,NumCols,View);
      }
      superblocks_[srow*GetNumSuperBlocks()+scol] = SuperBlock;
    }
  }

  // =============================================================
  // Construct the map extractors
  // =============================================================

  std::vector<Teuchos::RCP<const Epetra_Map> > DomainMaps(GetNumSuperBlocks(),Teuchos::null);
  std::vector<Teuchos::RCP<const Epetra_Map> > RangeMaps (GetNumSuperBlocks(),Teuchos::null);
  Teuchos::RCP<Epetra_Map> FullRangeMap  = Teuchos::null;
  Teuchos::RCP<Epetra_Map> FullDomainMap = Teuchos::null;

  for(int scol=0;scol<GetNumSuperBlocks();scol++)
    DomainMaps[scol] =
      Teuchos::rcp(&(superblocks_[0*GetNumSuperBlocks()+scol]->OperatorDomainMap()),false);
  FullDomainMap = MultiMapExtractor::MergeMaps(DomainMaps);
    domain_ex_  = Teuchos::rcp(new MultiMapExtractor(*FullDomainMap,DomainMaps));

  for(int srow=0;srow<GetNumSuperBlocks();srow++)
    RangeMaps[srow] =
      Teuchos::rcp(&(superblocks_[srow*GetNumSuperBlocks()+0]->OperatorRangeMap()),false);
  FullRangeMap = MultiMapExtractor::MergeMaps(RangeMaps);
    range_ex_  = Teuchos::rcp(new MultiMapExtractor(*FullRangeMap,RangeMaps));

  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::MueluSmootherWrapper::Apply
(const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero) const
{

  // Convert to Tpetra
  Teuchos::RCP<Epetra_MultiVector> X_rcp =
    Teuchos::rcp(new Epetra_MultiVector(X));
  Teuchos::RCP<Xpetra::EpetraMultiVector> Xex =
    Teuchos::rcp(new Xpetra::EpetraMultiVector(X_rcp));
  Teuchos::RCP<MultiVector> Xx =
    Teuchos::rcp_dynamic_cast<MultiVector>(Xex);
  Teuchos::RCP<Epetra_MultiVector> Y_rcp =
    Teuchos::rcp(new Epetra_MultiVector(Y));
  Teuchos::RCP<Xpetra::EpetraMultiVector> Yex =
    Teuchos::rcp(new Xpetra::EpetraMultiVector(Y_rcp));
  Teuchos::RCP<MultiVector> Yx = Teuchos::rcp_dynamic_cast<MultiVector>(Yex);

  // Apply underlying smoother
  S_->Apply(*Yx,*Xx,InitialGuessIsZero);

  // Convert to Epetra
  const Teuchos::RCP<Epetra_MultiVector>& Ye =
    MueLu::Utils<double,int,int,Node>::MV2NonConstEpetraMV(Yx);
  Y.Update(1.0,*Ye,0.0);

  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::MueluHierarchyWrapper::MueluHierarchyWrapper
(Teuchos::RCP<Hierarchy> H): H_(H)
{
      P_ = Teuchos::rcp(new MueLu::EpetraOperator(H_));
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::MueluHierarchyWrapper::Apply(
    const Epetra_MultiVector& X,
    Epetra_MultiVector& Y,
    bool InitialGuessIsZero) const
{
  P_->ApplyInverse(X,Y);
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::MueluAMGWrapper::MueluAMGWrapper(
    Teuchos::RCP<SparseMatrix> A,
    int num_pde,
    int null_space_dim,
    Teuchos::RCP<std::vector<double> > null_space_data,
    Teuchos::ParameterList muelu_list):
A_                (A                ),
num_pde_          (num_pde          ),
null_space_dim_   (null_space_dim   ),
null_space_data_  (null_space_data  ),
muelu_list_       (muelu_list       )
{Setup();}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::MueluAMGWrapper::Setup()
{

  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::MueluAMGWrapper::Setup");

  //Prepare operator for MueLu
  Teuchos::RCP<Epetra_CrsMatrix> A_crs
    = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(A_->EpetraOperator());
  if(A_crs==Teuchos::null)
    dserror("Make sure that the input matrix is a Epetra_CrsMatrix (or derived)");
  Teuchos::RCP<CrsMatrix> mueluA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A_crs));
  Teuchos::RCP<CrsMatrixWrap> mueluA_wrap = Teuchos::rcp(new CrsMatrixWrap(mueluA));
  Teuchos::RCP<Matrix> mueluOp = Teuchos::rcp_dynamic_cast<Matrix>(mueluA_wrap);

  // Prepare null space vector for MueLu
  Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > rowMap = mueluA->getRowMap();
  Teuchos::RCP<MultiVector> nspVector =
    Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
        rowMap,null_space_dim_,true);
  for ( size_t i=0; i < Teuchos::as<size_t>(null_space_dim_); i++) {
    Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
    const size_t myLength = nspVector->getLocalLength();
    for(size_t j=0; j<myLength; j++) {
      nspVectori[j] = (*null_space_data_)[i*myLength+j];
    }
  }


  // Input num eq and offset in the final level.
  // The amalgamation factory needs this info!
  int offsetFineLevel =  A_->RowMap().MinAllGID();
  mueluOp->SetFixedBlockSize(num_pde_,offsetFineLevel);
  Teuchos::ParameterList& MatrixList = muelu_list_.sublist("Matrix");
  MatrixList.set<int>("DOF offset",offsetFineLevel);
  MatrixList.set<int>("number of equations",num_pde_);



  // Build up hierarchy
  ParameterListInterpreter mueLuFactory(muelu_list_);
  H_ = mueLuFactory.CreateHierarchy();
  H_->SetDefaultVerbLevel(MueLu::Extreme); // TODO sure?
  H_->GetLevel(0)->Set("A", mueluOp);
  H_->GetLevel(0)->Set("Nullspace", nspVector);
  H_->GetLevel(0)->setlib(Xpetra::UseEpetra);
  H_->setlib(Xpetra::UseEpetra);
  mueLuFactory.SetupHierarchy(*H_);

  // Create the V-cycle
  P_ = Teuchos::rcp(new MueLu::EpetraOperator(H_));

  return;
}



/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::MueluAMGWrapper::Apply(
    const Epetra_MultiVector& X,
    Epetra_MultiVector& Y,
    bool InitialGuessIsZero) const
{
  P_->ApplyInverse(X,Y);
  return;
}




/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::HierarchyRemainderWrapper::HierarchyRemainderWrapper
  (Teuchos::RCP<Richardson_Vcycle_Operator> S,int start_level):
  start_level_(start_level),S_(S) {}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::HierarchyRemainderWrapper::HierarchyRemainderWrapper::Apply
  (const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero) const
{
  S_->Apply(X,Y,start_level_);
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::IfpackWrapper::IfpackWrapper(
    Teuchos::RCP<SparseMatrixBase> A,
    Teuchos::ParameterList& list)
  : A_(A)
{

  // Determine the preconditioner type
  type_ = list.get<std::string>("type","none");
  if(type_=="none")
    dserror("The type of preconditioner has to be provided.");

  // Extract list of parameters
  if(not(list.isSublist("ParameterList")))
    dserror("The parameter list has to be provided");
  list_ = list.sublist("ParameterList");

  // Create smoother
  Ifpack Factory;
  Arow_ = Teuchos::rcp_dynamic_cast<Epetra_RowMatrix>(A_->EpetraMatrix());
  if (Arow_==Teuchos::null)
    dserror("Something wrong. Be sure that the given matrix is not a block matrix");
  prec_ = Factory.Create(type_,Arow_.get());

  // Set parameter list and setup
  prec_->SetParameters(list_);
  prec_->Initialize();
  prec_->Compute();

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::IfpackWrapper::Apply(
    const Epetra_MultiVector& X,
    Epetra_MultiVector& Y,
    bool InitialGuessIsZero) const
{
  prec_->ApplyInverse(X,Y);
  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::DirectSolverWrapper::DirectSolverWrapper():
  solver_(Teuchos::null),
  A_(Teuchos::null),
  x_(Teuchos::null),
  b_(Teuchos::null),
  isSetUp_(false)
{}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::DirectSolverWrapper::Setup(
    Teuchos::RCP<LINALG::SparseMatrix>     matrix)
{

  //Set matrix
  A_ = Teuchos::rcp_dynamic_cast<Epetra_Operator>(matrix->EpetraMatrix());
  Teuchos::RCP<Epetra_CrsMatrix> crsA =
    Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(A_);
  if(crsA == Teuchos::null) dserror("Something wrong");

  // Set sol vector and rhs
  x_ = Teuchos::rcp(new Epetra_MultiVector(A_->OperatorDomainMap(),1));
  b_ = Teuchos::rcp(new Epetra_MultiVector(A_->OperatorRangeMap(),1));

  // Create linear solver
  Teuchos::RCP<Teuchos::ParameterList> solvparams = Teuchos::rcp(new Teuchos::ParameterList);
  solvparams->set("solver","klu");
  solver_ = Teuchos::rcp(new LINALG::Solver(solvparams,A_->Comm(),NULL));

  // Set up solver
  solver_->Setup(A_,x_,b_,true,true);

  isSetUp_=true;
  return;

}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::DirectSolverWrapper::Apply (
    const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero) const
{
  if (not(isSetUp_))
    dserror("The DirectSolverWrapper class should be set up before calling Apply");

  b_->Update( 1., X, 0. );
  solver_->Solve(A_,x_,b_,false,false);
  Y.Update( 1., *x_, 0. );

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::MergeAndSolve::MergeAndSolve():
  solver_(Teuchos::null),
  sparse_matrix_(Teuchos::null),
  A_(Teuchos::null),
  x_(Teuchos::null),
  b_(Teuchos::null),
  isSetUp_(false)
{}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::MergeAndSolve::Setup(
    Teuchos::RCP<LINALG::BlockSparseMatrixBase>     matrix)
{

  //Set matrix
  sparse_matrix_ = matrix->Merge();
  A_ = Teuchos::rcp_dynamic_cast<Epetra_Operator>(sparse_matrix_->EpetraMatrix());
  Teuchos::RCP<Epetra_CrsMatrix> crsA =
    Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(A_);
  if(crsA == Teuchos::null) dserror("Houston, something went wrong in merging the matrix");

  // Set sol vector and rhs
  x_ = Teuchos::rcp(new Epetra_MultiVector(A_->OperatorDomainMap(),1));
  b_ = Teuchos::rcp(new Epetra_MultiVector(A_->OperatorRangeMap(),1));

  // Create linear solver
  solver_ = Teuchos::rcp(new LINALG::Solver(A_->Comm(),NULL));

  // Set up solver
  solver_->Setup(A_,x_,b_,true,true);

  isSetUp_=true;
  return;

}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::MergeAndSolve::Apply (
    const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero) const
{
  if (not(isSetUp_))
    dserror("The MergeAndSolve class should be set up before calling Apply");

  b_->Update( 1., X, 0. );
  solver_->Solve(A_,x_,b_,false,false);
  Y.Update( 1., *x_, 0. );

  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::SIMPLE_BlockSmoother:: SIMPLE_BlockSmoother(
    Teuchos::RCP<LINALG::SparseOperator>  App,
    Teuchos::RCP<LINALG::SparseOperator>  Ass,
    Teuchos::RCP<LINALG::SparseOperator>  Aps,
    Teuchos::RCP<LINALG::SparseOperator>  Asp,
    Teuchos::RCP<LINALG::SparseOperator>  invApp,
    Teuchos::RCP<LINALG::SparseOperator>  S,
    Teuchos::RCP<AMGnxn_SmootherBase>  Smoother_App,
    Teuchos::RCP<AMGnxn_SmootherBase>  Smoother_S,
    Teuchos::RCP<MultiMapExtractor> range_ex,
    Teuchos::RCP<MultiMapExtractor> domain_ex,
    int p,
    int s,
    int iter,
    double omega,
    double alpha,
    std::string algorithm,
    std::string correction):
App_           (App         ),
Ass_           (Ass         ),
Aps_           (Aps         ),
Asp_           (Asp         ),
invApp_        (invApp      ),
S_             (S           ),
Smoother_App_  (Smoother_App),
Smoother_S_    (Smoother_S  ),
range_ex_      (range_ex    ),
domain_ex_     (domain_ex   ),
p_             (p           ),
s_             (s           ),
iter_          (iter        ),
omega_         (omega       ),
alpha_         (alpha       ),
algorithm_     (algorithm   ),
correction_    (correction  ),
DXp_      (Teuchos::null),
DXs_      (Teuchos::null),
DYp_      (Teuchos::null),
DYs_      (Teuchos::null),
DXp_aux_  (Teuchos::null),
DXs_aux_  (Teuchos::null),
DYp_aux_  (Teuchos::null)
{}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::SIMPLE_BlockSmoother::Apply (
    const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero) const
{


  //Extract blocks
  Teuchos::RCP<Epetra_MultiVector> Xp  = range_ex_->ExtractVector(X,p_);
  Teuchos::RCP<Epetra_MultiVector> Xs  = range_ex_->ExtractVector(X,s_);
  Teuchos::RCP<Epetra_MultiVector> Yp  = domain_ex_->ExtractVector(Y,p_);
  Teuchos::RCP<Epetra_MultiVector> Ys  = domain_ex_->ExtractVector(Y,s_);

  // Allocate working vectors only if necessary
  int NV = X.NumVectors();
  if (DXp_ == Teuchos::null)
  {
    const Teuchos::RCP<const Epetra_Map>& range_p  = range_ex_->Map(p_);
    const Teuchos::RCP<const Epetra_Map>& domain_p = domain_ex_->Map(p_);
    const Teuchos::RCP<const Epetra_Map>& range_s  = range_ex_->Map(s_);
    const Teuchos::RCP<const Epetra_Map>& domain_s = domain_ex_->Map(s_);
    DXp_      = Teuchos::rcp(new Epetra_MultiVector( *range_p,NV));
    DXs_      = Teuchos::rcp(new Epetra_MultiVector( *range_s,NV));
    DYp_      = Teuchos::rcp(new Epetra_MultiVector(*domain_p,NV));
    DYs_      = Teuchos::rcp(new Epetra_MultiVector(*domain_s,NV));
    DXp_aux_  = Teuchos::rcp(new Epetra_MultiVector( *range_p,NV));
    DXs_aux_  = Teuchos::rcp(new Epetra_MultiVector( *range_s,NV));
    DYp_aux_  = Teuchos::rcp(new Epetra_MultiVector(*domain_p,NV));
  }
  else if (DXp_->NumVectors() != NV)
  {
    const Teuchos::RCP<const Epetra_Map>& range_p  = range_ex_->Map(p_);
    const Teuchos::RCP<const Epetra_Map>& domain_p = domain_ex_->Map(p_);
    const Teuchos::RCP<const Epetra_Map>& range_s  = range_ex_->Map(s_);
    const Teuchos::RCP<const Epetra_Map>& domain_s = domain_ex_->Map(s_);
    DXp_      = Teuchos::rcp(new Epetra_MultiVector( *range_p,NV));
    DXs_      = Teuchos::rcp(new Epetra_MultiVector( *range_s,NV));
    DYp_      = Teuchos::rcp(new Epetra_MultiVector(*domain_p,NV));
    DYs_      = Teuchos::rcp(new Epetra_MultiVector(*domain_s,NV));
    DXp_aux_  = Teuchos::rcp(new Epetra_MultiVector( *range_p,NV));
    DXs_aux_  = Teuchos::rcp(new Epetra_MultiVector( *range_s,NV));
    DYp_aux_  = Teuchos::rcp(new Epetra_MultiVector(*domain_p,NV));
  }

  // Run several sweeps
  for(int k=0;k<iter_;k++)
  {

    if (algorithm_ == "simple")
    {
      // Compute tentative increment
      DXp_->Update(1.0,*Xp,0.0);  // DXp = Xp;
      App_->Apply(*Yp,*DXp_aux_);
      DXp_->Update(-1.0,*DXp_aux_,1.0); // DXp = Xp - App*Yp;
      Aps_->Apply(*Ys,*DXp_aux_);
      DXp_->Update(-1.0,*DXp_aux_,1.0); // DXp = Xp - App*Yp - Aps*Ys;
      DYp_->Scale(0.0);
      Smoother_App_->Apply(*DXp_,*DYp_); //DYp  = App^-1(Xp - App*Yp - Aps*Ys)

      // Compute Schur complement equation
      DXs_->Update(1.0,*Xs,0.0); // DXs_ = Xs
      Asp_->Apply(*Yp,*DXs_aux_);
      DXs_->Update(-1.0,*DXs_aux_,1.0); // DXs = Xs - Asp*Yp
      Ass_->Apply(*Ys,*DXs_aux_);
      DXs_->Update(-1.0,*DXs_aux_,1.0);// DXs = Xs - Asp*Yp - Ass*Ys
      Asp_->Apply(*DYp_,*DXs_aux_);
      DXs_->Update(-1.0,*DXs_aux_,1.0); // DXs = Xs - Asp*Yp - Ass*Ys  - Asp*DYp
      DYs_->Scale(0.0);
      Smoother_S_->Apply(*DXs_,*DYs_); // DYs = S^-1(Xs - Asp*Yp - Ass*Ys  - Asp*DYp)
      DYs_->Scale(alpha_);

      // Correction
      Aps_->Apply(*DYs_,*DXp_aux_); // DXp_aux = Aps*DYs
      if (correction_ == "approximated inverse")
        invApp_->Apply(*DXp_aux_,*DYp_aux_);// DYp_aux = App^-1*Aps*DYs
      else if (correction_ == "smoother")
      {
        DYp_aux_->Scale(0.0);
        Smoother_App_->Apply(*DXp_aux_,*DYp_aux_);// DYp_aux = App^-1*Aps*DYs
      }
      else dserror("Invalid strategy for computing the correction. Given correction = %s",
          correction_.c_str());
      DYp_->Update(-1.0,*DYp_aux_,1.0); // DYp = DYp - App^-1*Aps*DYs
    }
    else if (algorithm_ == "schur FSIAMG")
    {
      // Compute Schur complement equation
      DXs_->Update(1.0,*Xs,0.0); // DXs = Xs
      Asp_->Apply(*Yp,*DXs_aux_);
      DXs_->Update(-1.0,*DXs_aux_,1.0); // DXs = Xs - Asp*Yp
      Ass_->Apply(*Ys,*DXs_aux_);
      DXs_->Update(-1.0,*DXs_aux_,1.0);// DXs_ = Xs - Asp*Yp - Ass*Ys
      DYs_->Scale(0.0);
      Smoother_S_->Apply(*DXs_,*DYs_); // DYs = S^-1(Xs - Asp*Yp - Ass*Ys)
      DYs_->Scale(alpha_);

      // Correction
      DXp_->Update(1.0,*Xp,0.0); // DXp = Xp
      App_->Apply(*Yp,*DXp_aux_);
      DXp_->Update(-1.0,*DXp_aux_,1.0); // DXp = Xp - App*Yp
      Aps_->Apply(*Ys,*DXp_aux_);
      DXp_->Update(-1.0,*DXp_aux_,1.0);//DXp_ = Xp - App*Yp - Aps*Ys
      Aps_->Apply(*DYs_,*DXp_aux_);
      DXp_->Update(-1.0,*DXp_aux_,1.0); //DXp_ = Xp - App*Yp - Aps*Ys - Aps*DYs
      DYp_->Scale(0.0);
      Smoother_App_->Apply(*DXp_,*DYp_); // DYp = App^-1*(Xp - App*Yp - Aps*Ys - Aps*DYs)
    }
    else
      dserror ("Invalid value for algorithm. Given value = %s",algorithm_.c_str());

    // Update
    Yp->Update(omega_,*DYp_,1.0);
    Ys->Update(omega_,*DYs_,1.0);

  }

  domain_ex_->InsertVector(*Yp,p_,Y);
  domain_ex_->InsertVector(*Ys,s_,Y);
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::BGS_BlockSmoother::BGS_BlockSmoother(
    std::vector< Teuchos::RCP<LINALG::SparseOperator> > blocks,
    std::vector< Teuchos::RCP<AMGnxn_SmootherBase> > smoothers,
    std::vector<int> indices,
    Teuchos::RCP<MultiMapExtractor> range_ex,
    Teuchos::RCP<MultiMapExtractor> domain_ex,
    int size,
    int iter,
    double omega
    ):
blocks_     (blocks   ),
smoothers_  (smoothers),
indices_    (indices  ),
range_ex_   (range_ex ),
domain_ex_  (domain_ex),
size_       (size     ),
iter_       (iter     ),
omega_      (omega    )
{}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::BGS_BlockSmoother::Apply(
    const Epetra_MultiVector& X,
    Epetra_MultiVector& Y,
    bool InitialGuessIsZero) const
{

  int NV = X.NumVectors();

  // Auxiliary vectors
  std::vector< Teuchos::RCP<Epetra_MultiVector> > Y_blocks(size_);
  for(int i=0;i<size_;i++)
    Y_blocks[i] = domain_ex_->ExtractVector(Y,indices_[i]);


  // Several sweeps
  for(int k=0;k<iter_;k++)
  {
    // loop in blocks
    for(int i=0;i<size_;i++)
    {

      // Compute residual
      Teuchos::RCP<Epetra_MultiVector> DXi = range_ex_->ExtractVector(X,indices_[i]);
      Teuchos::RCP<Epetra_MultiVector> DXi_tmp = Teuchos::rcp(new Epetra_MultiVector(DXi->Map(),NV));
      for(int j=0;j<size_;j++)
      {
        blocks_[i*size_+j]->Apply(*Y_blocks[j],*DXi_tmp);
        DXi->Update(-1.0,*DXi_tmp,1.0);
      }

      // Solve Diagonal block
      Teuchos::RCP<Epetra_MultiVector> DYi =
        Teuchos::rcp(new Epetra_MultiVector(Y_blocks[i]->Map(),NV));
      DYi->PutScalar(0.0);
      smoothers_[i]->Apply(*DXi,*DYi,true);

      // Update
      Y_blocks[i]->Update(omega_,*DYi,1.0);

    }// loop in blocks

  }// Several sweeps

  // Insert vectors in the right place
  for(int i=0;i<size_;i++)
    domain_ex_->InsertVector(*(Y_blocks[i]),indices_[i],Y);

  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::AMG_BlockSmoother::AMG_BlockSmoother(
    Teuchos::RCP<BlockSparseMatrixBase> A,
    std::vector<int> num_pdes,
    std::vector<int> null_spaces_dim,
    std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data,
    const Teuchos::ParameterList& amgnxn_params,
    const Teuchos::ParameterList& smoothers_params)
{
  P_ = Teuchos::rcp(new AMGnxn_Operator(
        A,num_pdes,null_spaces_dim,null_spaces_data,amgnxn_params,smoothers_params,smoothers_params
        ));
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
void LINALG::SOLVER::AMG_BlockSmoother::Apply(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y, bool InitialGuessIsZero) const
{
  P_->ApplyInverse(X,Y);
  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

LINALG::SOLVER::AMGnxn_SmootherFactoryBase::AMGnxn_SmootherFactoryBase():
set_operator_              (false),
set_params_                (false),
set_params_subsolver_      (false),
set_hierarchies_           (false),
set_level_                 (false),
set_block_                 (false),
set_blocks_                (false),
set_type_                  (false),
set_verbosity_             (false),
set_null_space_            (false),
set_null_space_all_blocks_ (false)
{}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SparseOperator>
LINALG::SOLVER::AMGnxn_SmootherFactoryBase::GetOperator(){return operator_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::ParameterList
LINALG::SOLVER::AMGnxn_SmootherFactoryBase::GetParams() {return params_;};

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::ParameterList
LINALG::SOLVER::AMGnxn_SmootherFactoryBase::GetParamsSmoother(){return params_subsolver_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_Hierarchies>
LINALG::SOLVER::AMGnxn_SmootherFactoryBase::GetHierarchies(){return hierarchies_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGnxn_SmootherFactoryBase::GetLevel(){return level_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGnxn_SmootherFactoryBase::GetBlock(){return block_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<int> LINALG::SOLVER::AMGnxn_SmootherFactoryBase::GetBlocks(){return blocks_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::string  LINALG::SOLVER::AMGnxn_SmootherFactoryBase::GetSmootherName(){return subsolver_name_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::string  LINALG::SOLVER::AMGnxn_SmootherFactoryBase::GetType(){return type_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::string  LINALG::SOLVER::AMGnxn_SmootherFactoryBase::GetVerbosity(){return verbosity_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

LINALG::SOLVER::NullSpaceInfo  LINALG::SOLVER::AMGnxn_SmootherFactoryBase::GetNullSpace()
{return null_space_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<LINALG::SOLVER::NullSpaceInfo>
LINALG::SOLVER::AMGnxn_SmootherFactoryBase::GetNullSpaceAllBlocks(){return null_space_all_blocks_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactoryBase::SetOperator(Teuchos::RCP<SparseOperator> in)
{ set_operator_ = true; operator_ = in; return;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactoryBase::SetParams(const Teuchos::ParameterList& in)
{ set_params_ = true; params_ = in; return;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactoryBase::SetParamsSmoother(const Teuchos::ParameterList& in)
{ set_params_subsolver_= true; params_subsolver_= in; return;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactoryBase::SetHierarchies(Teuchos::RCP<AMGnxn_Hierarchies> in)
{ set_hierarchies_ = true; hierarchies_ = in; return;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactoryBase::SetLevel(int in)
{ set_level_ = true; level_ = in; return;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactoryBase::SetBlock(int in)
{ set_block_ = true; block_ = in; return;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactoryBase::SetBlocks(std::vector<int> in)
{ set_blocks_= true; blocks_= in; return;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactoryBase::SetSmootherName(std::string in)
{ set_subsolver_name_ = true; subsolver_name_ = in; return;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactoryBase::SetType(std::string in)
{ set_type_ = true; type_ = in; return;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactoryBase::SetVerbosity(std::string in)
{ set_verbosity_ = true; verbosity_ = in; return;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactoryBase::SetNullSpace(const NullSpaceInfo& in)
{ set_null_space_ = true; null_space_ = in; return;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactoryBase::SetNullSpaceAllBlocks
(const std::vector<NullSpaceInfo>& in)
{ set_null_space_all_blocks_ = true; null_space_all_blocks_ = in; return;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

bool LINALG::SOLVER::AMGnxn_SmootherFactoryBase::IsSetOperator(){return set_operator_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

bool LINALG::SOLVER::AMGnxn_SmootherFactoryBase::IsSetParams(){return set_params_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

bool LINALG::SOLVER::AMGnxn_SmootherFactoryBase::IsSetParamsSmoother()
{return set_params_subsolver_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

bool LINALG::SOLVER::AMGnxn_SmootherFactoryBase::IsSetHierarchies(){return set_hierarchies_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

bool LINALG::SOLVER::AMGnxn_SmootherFactoryBase::IsSetLevel(){return set_level_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

bool LINALG::SOLVER::AMGnxn_SmootherFactoryBase::IsSetBlock(){return set_block_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

bool LINALG::SOLVER::AMGnxn_SmootherFactoryBase::IsSetBlocks(){return set_blocks_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

bool LINALG::SOLVER::AMGnxn_SmootherFactoryBase::IsSetSmootherName(){return set_subsolver_name_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

bool LINALG::SOLVER::AMGnxn_SmootherFactoryBase::IsSetType(){return set_type_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

bool LINALG::SOLVER::AMGnxn_SmootherFactoryBase::IsSetVerbosity(){return set_verbosity_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

bool LINALG::SOLVER::AMGnxn_SmootherFactoryBase::IsSetNullSpace(){return set_null_space_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

bool LINALG::SOLVER::AMGnxn_SmootherFactoryBase::IsSetNullSpaceAllBlocks()
{return set_null_space_all_blocks_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/
LINALG::SOLVER::NullSpaceInfo::NullSpaceInfo(
    int num_pdes,
    int null_space_dim,
    Teuchos::RCP<std::vector<double> >  null_space_data):
  num_pdes_ (num_pdes),
  null_space_dim_ (null_space_dim),
  null_space_data_ (null_space_data)
{}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_SmootherFactory::SetTypeAndParams()
{

  // Valid types
  std::vector<std::string> valid_types;
  valid_types.push_back("BGS");
  valid_types.push_back("IFPACK");
  valid_types.push_back("REUSE_MUELU_SMOOTHER");
  valid_types.push_back("REUSE_MUELU_AMG");
  valid_types.push_back("NEW_MUELU_AMG");
  valid_types.push_back("DIRECT_SOLVER");
  valid_types.push_back("MERGE_AND_SOLVE");
  valid_types.push_back("BLOCK_AMG");

  std::string smoother_type;
  Teuchos::ParameterList smoother_params;
  if(GetParamsSmoother().isSublist(GetSmootherName()))
  {
    smoother_type  = GetParamsSmoother().sublist(GetSmootherName()).get<std::string>("type","none");
    smoother_params = GetParamsSmoother().sublist(GetSmootherName()).sublist("parameters");
  }
  else if (std::find(valid_types.begin(),valid_types.end(),GetSmootherName())!=valid_types.end())
    smoother_type = GetSmootherName();
  else
    smoother_type = "none";

  SetType(smoother_type);
  SetParams(smoother_params);
  return;

}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_SmootherBase> LINALG::SOLVER::AMGnxn_SmootherFactory::Create()
{

  // Expected parameters in GetParamsSmoother()
  //
  // <ParameterList name="mySmoother">
  //   <Parameter name="type"   type="string"  value="..."/>
  //   <ParameterList name="parameters">
  //
  //    ...    ...   ...   ...   ...   ...
  //
  //   </ParameterList>
  // </ParameterList>
  //
  // Available input?

  if (not IsSetParamsSmoother())
    dserror("IsSetParamsSmoother() returns false");
  if (not IsSetSmootherName())
    dserror("IsSetSmootherName() returns false");

  // Determine the type of smoother to be constructed and its parameters
  SetTypeAndParams();

  // Create the corresponding factory
  Teuchos::RCP<AMGnxn_SmootherFactoryBase> mySmootherFactory = Teuchos::null;

  if (GetType() == "BGS")
  {
    mySmootherFactory = Teuchos::rcp(new BGS_BlockSmootherFactory());
    mySmootherFactory->SetOperator(GetOperator());
    mySmootherFactory->SetParams(GetParams());
    mySmootherFactory->SetParamsSmoother(GetParamsSmoother());
    mySmootherFactory->SetHierarchies(GetHierarchies());
    mySmootherFactory->SetBlocks(GetBlocks());
    mySmootherFactory->SetNullSpaceAllBlocks(GetNullSpaceAllBlocks());
  }
  else if (GetType() == "BLOCK_AMG")
  {
    mySmootherFactory = Teuchos::rcp(new AMG_BlockSmootherFactory());
    mySmootherFactory->SetOperator(GetOperator());
    mySmootherFactory->SetParams(GetParams());
    mySmootherFactory->SetParamsSmoother(GetParamsSmoother());
    mySmootherFactory->SetBlocks(GetBlocks());
    mySmootherFactory->SetNullSpaceAllBlocks(GetNullSpaceAllBlocks());
  }
  else if (GetType() == "SIMPLE")
  {
    mySmootherFactory = Teuchos::rcp(new SIMPLE_BlockSmootherFactory());
    mySmootherFactory->SetOperator(GetOperator());
    mySmootherFactory->SetParams(GetParams());
    mySmootherFactory->SetParamsSmoother(GetParamsSmoother());
    mySmootherFactory->SetHierarchies(GetHierarchies());
    mySmootherFactory->SetBlocks(GetBlocks());
    mySmootherFactory->SetNullSpaceAllBlocks(GetNullSpaceAllBlocks());
  }
  else if (GetType() == "MERGE_AND_SOLVE")
  {
    mySmootherFactory = Teuchos::rcp(new MergeAndSolveFactory());
    mySmootherFactory->SetOperator(GetOperator());
    mySmootherFactory->SetBlocks(GetBlocks());
    mySmootherFactory->SetLevel(GetLevel());
  }
  else if (GetType() == "DIRECT_SOLVER")
  {
    mySmootherFactory = Teuchos::rcp(new DirectSolverWrapperFactory());
    mySmootherFactory->SetOperator(GetOperator());
    mySmootherFactory->SetParams(GetParams());
    mySmootherFactory->SetBlock(GetBlock());
  }
  else if (GetType() == "IFPACK")
  {
    mySmootherFactory = Teuchos::rcp(new IfpackWrapperFactory());
    mySmootherFactory->SetOperator(GetOperator());
    mySmootherFactory->SetParams(GetParams());
    mySmootherFactory->SetBlock(GetBlock());
  }
  else if (GetType() == "REUSE_MUELU_SMOOTHER")
  {
    mySmootherFactory = Teuchos::rcp(new MueluSmootherWrapperFactory());
    mySmootherFactory->SetHierarchies(GetHierarchies());
    mySmootherFactory->SetBlock(GetBlock());
  }
  else if (GetType() == "REUSE_MUELU_AMG")
  {
    mySmootherFactory = Teuchos::rcp(new HierarchyRemainderWrapperFactory());
    mySmootherFactory->SetHierarchies(GetHierarchies());
    mySmootherFactory->SetBlock(GetBlock());
  }
  else if (GetType() == "NEW_MUELU_AMG")
  {
    mySmootherFactory = Teuchos::rcp(new MueluAMGWrapperFactory());
    mySmootherFactory->SetOperator(GetOperator());
    mySmootherFactory->SetHierarchies(GetHierarchies());
    mySmootherFactory->SetBlock(GetBlock());
    mySmootherFactory->SetParams(GetParams());
    mySmootherFactory->SetParamsSmoother(GetParamsSmoother());
    mySmootherFactory->SetNullSpace(GetNullSpace());
  }
  else
    dserror("Unknown smoother type. Fix your xml file");

  //Build the smoother
  mySmootherFactory->SetVerbosity(GetVerbosity());
  mySmootherFactory->SetLevel(GetLevel());
  return mySmootherFactory->Create();

}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_SmootherBase>
LINALG::SOLVER::IfpackWrapperFactory::Create()
{
  // Expected parameters with default values
  //
  // <ParameterList name="parameters">
  //   <Parameter name="type"                           type="string"  value="point relaxation"/>
  //   <ParameterList name="ParameterList">
  //     <Parameter name="relaxation: type"             type="string"  value="Gauss-Seidel"/>
  //     <Parameter name="relaxation: backward mode"    type="bool"    value="false"/>
  //     <Parameter name="relaxation: sweeps"           type="int"     value="2"/>
  //     <Parameter name="relaxation: damping factor"   type="double"  value="1.0"/>
  //   </ParameterList>
  // </ParameterList>
  //

  // Check input
  if (not IsSetParams())
    dserror("IsSetParams() returns false");
  if (not IsSetOperator())
    dserror("IsSetOperator() returns false");

  // Fill myparams with default values where required
  Teuchos::ParameterList defaults;
  defaults.set<std::string>("type","point relaxation");
  defaults.sublist("ParameterList").set<std::string>("relaxation: type", "Gauss-Seidel");
  defaults.sublist("ParameterList").set<bool>("relaxation: backward mode",false);
  defaults.sublist("ParameterList").set<int>("relaxation: sweeps",2);
  defaults.sublist("ParameterList").set<double>("relaxation: damping factor",1.0);
  Teuchos::ParameterList myParamsAux = GetParams();
  myParamsAux.setParametersNotAlreadySet(defaults);
  SetParams(myParamsAux);

  // Some output
  if (GetVerbosity()=="on")
  {
    std::cout << std::endl;
    std::cout << "Creating an IFPACK smoother for block " << GetBlock()
      << " at level " << GetLevel() << std::endl;
    std::cout << "The Ifpack type is: " << GetParams().get<std::string>("type") << std::endl;
    std::cout << "The parameters are: " << std::endl;
    std::cout << GetParams().sublist("ParameterList");
    //std::cout << std::endl;
  }

  // Build the smoother
  Teuchos::RCP<SparseMatrixBase> Op = Teuchos::rcp_dynamic_cast<SparseMatrixBase>(GetOperator());
  if(Op==Teuchos::null)
    dserror("Ifpack smoother only works for non blocked matrices");
  Teuchos::ParameterList myParams = GetParams();
  return  Teuchos::rcp(new IfpackWrapper(Op,myParams));
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_SmootherBase>
LINALG::SOLVER::MueluSmootherWrapperFactory::Create()
{
  // Check input
  if (not IsSetLevel())
    dserror("IsSetLevel() returns false");
  if (not IsSetBlock())
    dserror("IsSetBlock() returns false");
  if (not IsSetHierarchies())
    dserror("IsSetHierarchies() returns false");

  if (GetVerbosity()=="on")
  {
    std::cout << std::endl;
    std::cout << "Creating a REUSE_MUELU_SMOOTHER smoother for block " << GetBlock() ;
    std::cout << " at level " << GetLevel() << std::endl;
  }
  return GetHierarchies()->GetSPre(GetBlock(),GetLevel());
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_SmootherBase>
LINALG::SOLVER::MueluAMGWrapperFactory::Create()
{


  // Expected parameters (example)
  // <ParameterList name="parameters">
  //   <Parameter name="xml file"      type="string"  value="myfile.xml"/>
  // </ParameterList>
  //
  // TODO or
  //
  // <ParameterList name="parameters">
  //   <Parameter name="parameter list"      type="string"  value="NameOfTheParameterList"/>
  // </ParameterList>
  //  In that case we look in GetParamsSmoother() for a list called "NameOfTheParameterList"
  //  which has to contain all the parameters defining a muelu hierarchy
  //
  // TODO or
  //
  // <ParameterList name="parameters">
  //  ... ... list defining the muelue hierarcy (i.e.) the contents of the xml file
  // </ParameterList>
  //
  //
  // Priority: first "xml file", then "parameter list", then other:
  // If the parameter "xml file" is found, then all other parameters are ignored
  // else, if "parameter list is found", then other parameters are ignored

  // Check input
  if (not IsSetLevel())
    dserror("IsSetLevel() returns false");
  if (not IsSetOperator())
    dserror("IsSetOperator() returns false");
  if (not IsSetBlock())
    dserror("IsSetBlock() returns false");
  if (not IsSetHierarchies())
    dserror("IsSetHierarchies() returns false");
  if (not IsSetParams())
    dserror("IsSetParams() returns false");
  if (not IsSetNullSpace())
    dserror("IsSetNullSpace() returns false");
  if (not IsSetParamsSmoother())
    dserror("IsSetSmoothersParams() returns false");

  if (GetVerbosity()=="on")
  {
    std::cout << std::endl;
    std::cout << "Creating a NEW_MUELU_AMG smoother for block " << GetBlock() ;
    std::cout << " at level " << GetLevel() << std::endl;
  }

  Teuchos::ParameterList myList;
  std::string xml_filename = GetParams().get<std::string>("xml file","none");
  std::string list_name = GetParams().get<std::string>("parameter list","none");
  if(xml_filename != "none")
  {
    Teuchos::updateParametersFromXmlFile(
        xml_filename,Teuchos::Ptr<Teuchos::ParameterList>(&myList));
    if (GetVerbosity()=="on")
    {
      std::cout << "The chosen parameters are:" << std::endl;
      std::cout << "xml file = : " << xml_filename << std::endl;
    }
  }
  else if (list_name != "none")
  {
    myList = GetParamsSmoother().sublist(list_name);
    if (GetVerbosity()=="on")
    {
      std::cout << "The chosen parameters are:" << std::endl;
      std::cout << "parameter list = : " << list_name << std::endl;
    }
  }
  else
    myList = GetParams();



  // TODO now we use the null space generated by baci, which only makes sense for the finest level.
  // We can obtain null spaces for other levels from inside the muelu hierarchies.
  if(GetLevel()!=0)
    dserror("Trying to create a NEW_MUELU_AMG smoother at a level > 0. Sorry, but this is not possible yet.");

  // Recover info
  Teuchos::RCP<SparseMatrix> A = Teuchos::rcp_dynamic_cast<SparseMatrix>(GetOperator());
  int num_pde = GetNullSpace().GetNumPDEs();
  int null_space_dim = GetNullSpace().GetNullSpaceDim();
  Teuchos::RCP<std::vector<double> > null_space_data = GetNullSpace().GetNullSpaceData();


  return Teuchos::rcp(new MueluAMGWrapper(A,num_pde,null_space_dim,null_space_data,myList));

}






/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_SmootherBase>
LINALG::SOLVER::HierarchyRemainderWrapperFactory::Create()
{
  // Check input
  if (not IsSetLevel())
    dserror("IsSetLevel() returns false");
  if (not IsSetBlock())
    dserror("IsSetBlock() returns false");
  if (not IsSetHierarchies())
    dserror("IsSetHierarchies() returns false");

  if (GetVerbosity()=="on")
  {
    std::cout << std::endl;
    std::cout << "Creating a REUSE_MUELU_AMG smoother for block " << GetBlock() ;
    std::cout << " at level " << GetLevel() << std::endl;
  }

  // Pick up info and cast
  int NumLevels = GetHierarchies()->GetNumLevels(GetBlock());
  std::vector< Teuchos::RCP<Epetra_Operator> > Avec(NumLevels,Teuchos::null);
  std::vector< Teuchos::RCP<Epetra_Operator> > Pvec(NumLevels-1,Teuchos::null);
  std::vector< Teuchos::RCP<Epetra_Operator> > Rvec(NumLevels-1,Teuchos::null);
  std::vector< Teuchos::RCP<AMGnxn_SmootherBase> > SvecPre(NumLevels,Teuchos::null);
  std::vector< Teuchos::RCP<AMGnxn_SmootherBase> > SvecPos(NumLevels-1,Teuchos::null);
  for(int level=0;level<NumLevels;level++)
  {
    Avec[level] =
      Teuchos::rcp_dynamic_cast<Epetra_Operator>(GetHierarchies()->GetA(GetBlock(),level));
    SvecPre[level] =
      Teuchos::rcp_dynamic_cast<AMGnxn_SmootherBase>(GetHierarchies()->GetSPre(GetBlock(),level));
  }
  for(int level=0;level<(NumLevels-1);level++)
  {
    Pvec[level] =
      Teuchos::rcp_dynamic_cast<Epetra_Operator>(GetHierarchies()->GetP(GetBlock(),level));
    Rvec[level] =
      Teuchos::rcp_dynamic_cast<Epetra_Operator>(GetHierarchies()->GetR(GetBlock(),level));
    SvecPos[level] =
      Teuchos::rcp_dynamic_cast<AMGnxn_SmootherBase>(GetHierarchies()->GetSPos(GetBlock(),level));
  }

  // Construct the V cycle
  Teuchos::RCP<Richardson_Vcycle_Operator> V =
    Teuchos::rcp( new Richardson_Vcycle_Operator(NumLevels));
  V->SetOperators  (Avec);
  V->SetProjectors (Pvec);
  V->SetRestrictors(Rvec);
  V->SetPreSmoothers (SvecPre);
  V->SetPosSmoothers (SvecPos);

  Teuchos::RCP<AMGnxn_SmootherBase> S = Teuchos::rcp( new HierarchyRemainderWrapper(V,GetLevel()));

  return S;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_SmootherBase>
LINALG::SOLVER::MergeAndSolveFactory::Create()
{
  // Check input
  if (not IsSetOperator())
    dserror("IsSetOperator() returns false");

  if (GetVerbosity()=="on")
  {
    std::cout << std::endl;
    std::cout << "Creating a MERGE_AND_SOLVE smoother for blocks (" ;
    for(size_t i=0;i<GetBlocks().size();i++)
    {
      std::cout << GetBlocks()[i] ;
      if (i<GetBlocks().size()-1)
        std::cout << ", ";
      else
        std::cout << ") ";
    }
    std::cout << " at level " << GetLevel() << std::endl;
  }

  Teuchos::RCP<MergeAndSolve> S = Teuchos::rcp(new MergeAndSolve);
  Teuchos::RCP<BlockSparseMatrixBase>  matrix
    = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(GetOperator());
  if(matrix==Teuchos::null)
    dserror("We expect here a block sparse matrix");
  S->Setup(matrix);

  Teuchos::RCP<AMGnxn_SmootherBase> SBase = Teuchos::rcp_dynamic_cast<AMGnxn_SmootherBase>(S);
  if (SBase==Teuchos::null)
    dserror("Something went wrong");

  return S;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_SmootherBase> LINALG::SOLVER::AMG_BlockSmootherFactory::Create()
{

  //<ParameterList name="parameters">
  //
  //  <Parameter name="number of levels"                 type="int"  value="..."/>
  //
  //  <Parameter name="smoother: all but coarsest level" type="string"  value="myFinestSmoother"/>
  //
  //  <Parameter name="smoother: coarsest level"         type="string"  value="myCoarsestSmoother"/>
  //
  //  <Parameter name="verbosity"                        type="string"  value="on"/>
  //
  //  <Parameter name="muelu parameters for block 0"       type="string"  value="myMuelu0"/>
  //
  //  <Parameter name="muelu parameters for block 1"       type="string"  value="myMuelu1"/>
  //
  //   ....
  //
  //  <Parameter name="muelu parameters for block N"       type="string"  value="myMueluN"/>
  //
  //</ParameterList>
  // WARNING: here the blocks are in local numeration of the submatrix passed to this factory

  // Recover the operator
  Teuchos::RCP<BlockSparseMatrixBase> A =
    Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(GetOperator());
  if (A==Teuchos::null)
    dserror("You have to provide a block sparse matrix in AMG_BlockSmootherFactory");

  // Recover the null space info
  int nBlocks = GetBlocks().size();
  const std::vector<int>& Blocks = GetBlocks();
  int b = 0;
  std::vector<int> num_pdes(nBlocks,0);
  std::vector<int> null_spaces_dim(nBlocks,0);
  std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data(nBlocks,Teuchos::null);
  for (int i=0;i<nBlocks;i++)
  {
    b = Blocks[i];
    num_pdes[i] = GetNullSpaceAllBlocks()[b].GetNumPDEs();
    null_spaces_dim[i] = GetNullSpaceAllBlocks()[b].GetNullSpaceDim();
    null_spaces_data[i] = GetNullSpaceAllBlocks()[b].GetNullSpaceData();
  }

  // Recover the lists
  const Teuchos::ParameterList& amgnxn_params = GetParams();
  const Teuchos::ParameterList& smoothers_params = GetParamsSmoother();

  return Teuchos::rcp(new AMG_BlockSmoother(
        A, num_pdes, null_spaces_dim, null_spaces_data, amgnxn_params, smoothers_params));

}




/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_SmootherBase>
LINALG::SOLVER::BGS_BlockSmootherFactory::Create()
{

  // Expected parameters (example)
  // <ParameterList name="parameters">
  //   <Parameter name="blocks"      type="string"  value="(1,2),(3,4),(5)"/>
  //   <Parameter name="smoothers"   type="string"  value="myBGS,mySIMPLE,IFPACK"/>
  //   <Parameter name="sweeps"      type="int"     value="3"/>
  //   <Parameter name="omega"       type="double"  value="1.0"/>
  // </ParameterList>

  // TODO Check that all required data is set

  // =============================================================
  // Parse parameters
  // =============================================================

  // determine how the blocks are grouped
  std::string blocks_string = GetParams().get<std::string>("blocks","none");
  std::vector<std::vector<int> > SuperBlocks2Blocks;
  std::vector<std::vector<int> > SuperBlocks2BlocksLocal;
  ParseBlocks(blocks_string,GetBlocks(),SuperBlocks2Blocks,SuperBlocks2BlocksLocal);

  // std::cout << "======================" << std::endl;
  // for(size_t i=0;i<SuperBlocks2Blocks.size();i++)
  // {
  //   for(size_t j=0;j<SuperBlocks2Blocks[i].size();j++)
  //     std::cout << SuperBlocks2Blocks[i][j] << ", ";
  //   std::cout << std::endl;
  // }


  // Determine the subsolver names
  std::string smoothers_string = GetParams().get<std::string>("smoothers","none");
  std::vector<std::string> SubSolverNames;
  ParseSmootherNames(smoothers_string,SubSolverNames,SuperBlocks2Blocks);

  // sweeps and damping
  int iter = GetParams().get<int>("sweeps",3);
  double omega = GetParams().get<double>("omega",1.0);

  // =============================================================
  // Some output
  // =============================================================
  if (GetVerbosity()=="on")
  {
    std::cout << std::endl;
    std::cout << "Creating a BGS smoother for blocks (" ;
    for(size_t i=0;i<GetBlocks().size();i++)
    {
      std::cout << GetBlocks()[i] ;
      if (i<GetBlocks().size()-1)
        std::cout << ", ";
      else
        std::cout << ") ";
    }
    std::cout << " at level " << GetLevel() << std::endl;
    std::cout << "The chosen parameters are" << std::endl;
    std::cout << "blocks = " ;
    for(size_t k=0;k<SuperBlocks2Blocks.size();k++)
    {
      std::cout << "(";
      for(size_t j=0;j<SuperBlocks2Blocks[k].size();j++)
      {
        std::cout << SuperBlocks2Blocks[k][j];
        if(j<(SuperBlocks2Blocks[k].size()-1))
          std::cout << ",";
      }
      if(k<(SuperBlocks2Blocks.size()-1))
        std::cout << "),";
      else
        std::cout << ")" << std::endl;
    }
    std::cout << "smoothers = ";
    for(size_t k=0;k<SubSolverNames.size();k++)
    {
      std::cout << SubSolverNames[k];
      if(k<(SubSolverNames.size()-1))
        std::cout << ",";
      else
        std::cout << std::endl;

    }
    std::cout << "sweeps = " << iter << std::endl;
    std::cout << "omega = " << omega << std::endl;
    //std::cout << std::endl;
  }

  // =============================================================
  // Rearrange blocks
  // =============================================================

  Teuchos::RCP<BlockSparseMatrixBase> OpBlocked = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(GetOperator());
  if (OpBlocked == Teuchos::null)
    dserror("I expect a block matrix here");
  BlockAggrupator myBlockAggrupator(OpBlocked,SuperBlocks2BlocksLocal);
  int NumSuperBlocks = SuperBlocks2Blocks.size();
  std::vector<int> indices(NumSuperBlocks,0);
  for(int i=0;i<(int)indices.size();i++)
    indices[i]=i;

  // =============================================================
  // Construct smoothers for diagonal superblocks
  // =============================================================

  std::vector< Teuchos::RCP<AMGnxn_SmootherBase> > SubSmoothers(NumSuperBlocks,Teuchos::null);
  for(int scol=0;scol<NumSuperBlocks;scol++)
  {
    AMGnxn_SmootherFactory mySmootherCreator;
    mySmootherCreator.SetSmootherName   (SubSolverNames[scol]);
    mySmootherCreator.SetParamsSmoother (GetParamsSmoother());
    mySmootherCreator.SetHierarchies    (GetHierarchies());
    mySmootherCreator.SetLevel          (GetLevel());
    mySmootherCreator.SetOperator       (myBlockAggrupator.GetSuperBlock(scol,scol));
    mySmootherCreator.SetVerbosity      (GetVerbosity());
    if (SuperBlocks2Blocks[scol].size()==1)
    {
      int thisblock = SuperBlocks2Blocks[scol][0];
      mySmootherCreator.SetBlock          (thisblock);
      mySmootherCreator.SetNullSpace      (GetNullSpaceAllBlocks()[thisblock]);
    }
    else
    {
      mySmootherCreator.SetBlocks (SuperBlocks2Blocks[scol]);
      mySmootherCreator.SetNullSpaceAllBlocks(GetNullSpaceAllBlocks());
    }

    SubSmoothers[scol] = mySmootherCreator.Create();
  }

  // =============================================================
  // Construct BGS smoother
  // =============================================================

  return Teuchos::rcp(new BGS_BlockSmoother(
        myBlockAggrupator.GetSuperBlocksRowOrder(),
        SubSmoothers,
        indices,
        myBlockAggrupator.GetRangeMapExtractor(),
        myBlockAggrupator.GetDomainMapExtractor(),
        myBlockAggrupator.GetNumSuperBlocks(),
        iter,
        omega));

}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_ParticularSmootherFactory::ParseBlocks(
    const std::string& block_string,
    const std::vector<int>& blocks,
    std::vector<std::vector<int> >& superblocks_to_blocks,
    std::vector<std::vector<int> >& superblocks_to_blocks_local)
{
  if (block_string == "none")
  {

    superblocks_to_blocks.resize(0);
    for(int i=0;i<(int)blocks.size();i++)
    {
      std::vector<int> vaux;
      vaux.push_back(blocks[i]);
      superblocks_to_blocks.push_back(vaux);
    }

    superblocks_to_blocks_local = superblocks_to_blocks;
    for(int i=0;i<(int)blocks.size();i++)
    {
      std::vector<int> vaux;
      vaux.push_back(i);
      superblocks_to_blocks_local[i] = vaux;
    }

    return;
  }

  // TODO add some checks

  // Parse the string
  superblocks_to_blocks.resize(0);
  bool brace_opened = false;
  std::vector<int> sb;
  std::string sint ="";
  for(int i=0;i<(int)block_string.size();i++)
  {
    std::string ch;
    ch.assign(1,block_string[i]);
    if (ch== "(")
    {
      if (brace_opened)
        dserror("Something wrong. Make sure you are setting correctly the blocks in your xml file");
      brace_opened = true;
      sb.resize(0);
    }
    else if (ch==")")
    {
      if (not(brace_opened))
        dserror("Something wrong. Make sure you are setting correctly the blocks in your xml file");
      sb.push_back(atoi(sint.c_str()));
      sint ="";
      superblocks_to_blocks.push_back(sb);
      brace_opened = false;
    }
    else if (ch==",")
    {
      if (brace_opened)
      {
        sb.push_back(atoi(sint.c_str()));
        sint ="";
      }
    }
    else
      sint += ch;
  }

  // Recover the local numeration
  superblocks_to_blocks_local = superblocks_to_blocks;
  for(int i=0;i<(int)superblocks_to_blocks.size();i++)
  {
    for(int j=0;j<(int)superblocks_to_blocks[i].size();j++)
    {
      if(std::find(blocks.begin(),blocks.end(),superblocks_to_blocks[i][j])==blocks.end() )
        dserror("Something wrong. Make sure in your xml file you are counting the blocks starting with 0 or you are not given to much blocks ");
      int pos =  std::find(blocks.begin(),blocks.end(),superblocks_to_blocks[i][j])
        - blocks.begin();
      superblocks_to_blocks_local[i][j] = pos;
    }
  }

  return;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::BGS_BlockSmootherFactory::ParseSmootherNames(
    const std::string& smoothers_string,
    std::vector<std::string>& smoothers_vector,std::vector<std::vector<int> > superblocks)
{

  if (smoothers_string=="none")
  {
    int NumSuperBlocks = superblocks.size();
    smoothers_vector.resize(0);
    for(int i=0;i<NumSuperBlocks;i++)
    {
      if(0 == (superblocks[i].size()) )
        dserror("Something wrong related with how the blocks are set in your xml file");
      else if(1 == (superblocks[i].size()) )
        smoothers_vector.push_back("IFPACK");
      else
        smoothers_vector.push_back("BGS");
    }
  }
  else
  {
    smoothers_vector.resize(0);
    std::string buf="";
    for(int i=0;i<(int)smoothers_string.size();i++)
    {
      std::string ch(1,smoothers_string[i]);
      if(ch==",")
      {
        smoothers_vector.push_back(buf);
        buf="";
      }
      else
        buf += ch;
    }
    if(not(buf==""))
      smoothers_vector.push_back(buf);
    buf="";
  }

  if (smoothers_vector.size() != superblocks.size() )
    dserror("Not given enough subsmoothers! Fix your xml file.");

  return ;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_SmootherBase>
LINALG::SOLVER::SIMPLE_BlockSmootherFactory::Create()
{

  // Expected parameters (example)
  // <ParameterList name="parameters">
  //   <Parameter name="algorithm"           type="string"  value="simple"/>
  //     <!-- Default: The standard simple algorithm  -->
  //   <Parameter name="algorithm"           type="string"  value="schur FSIAMG"/>
  //     <!-- The schur algorithm as it is implemented in the FSIAMG preconditioner  -->
  //   <Parameter name="predictor block"     type="string"  value="(1,2)"/>
  //   <Parameter name="predictor smoother"  type="string"  value="BGS"/>
  //   <Parameter name="predictor inverse"   type="string"  value="diagonal"/>
  //   <Parameter name="predictor inverse"   type="string"  value="row sums"/>
  //   <Parameter name="predictor inverse"   type="string"  value="row sums diagonal blocks"/>
  //   <Parameter name="schur block"         type="string"  value="(3)"/>
  //   <Parameter name="schur smoother"      type="string"  value="IFPACK"/>
  //   <Parameter name="correction"          type="string"  value="smoother"/>
  //   <Parameter name="correction"          type="string"  value="approximated inverse"/>
  //   <Parameter name="sweeps"              type="int"     value="3"/>
  //   <Parameter name="omega"               type="double"  value="1.0"/>
  //     <!-- Damping of the Richadson iteration -->
  //   <Parameter name="alpha"               type="double"  value="1.0"/>
  //     <!-- Damping of the "pressure" correction-->
  //   <Parameter name="beta"                type="double"  value="1.0"/>
  //     <!-- Coefficient that multiplies the approximate inverse of the predictor block
  //     in the Schur complement matrixCoefficient that multiplies the approximate inverse
  //     of the predictor block in the Schur complement matrix, i.e.
  //     S = A_22 - beta*A21*A11inv*A12-->
  // </ParameterList>

  // TODO Check that all required data is set

  // =============================================================
  // Parse parameters
  // =============================================================

  // determine how the blocks are grouped
  std::string predictor_block_string = GetParams().get<std::string>("predictor block","none");
  std::string schur_block_string     = GetParams().get<std::string>("schur block","none");
  if (predictor_block_string=="none")
    dserror("The field \"predictor block\" is mandatory for the SIMPLE smoother.Fix your xml file");
  if (schur_block_string=="none")
    dserror("The field \"schur block\" is mandatory for the SIMPLE smoother. Fix your xml file");
  std::string blocks_string = predictor_block_string + "," + schur_block_string;
  std::vector<std::vector<int> > SuperBlocks2Blocks;
  std::vector<std::vector<int> > SuperBlocks2BlocksLocal;
  ParseBlocks(blocks_string,GetBlocks(),SuperBlocks2Blocks,SuperBlocks2BlocksLocal);
  int pred=0;
  int schur=1;


  // Smoother names
  std::string predictor_smoother = GetParams().get<std::string>("predictor smoother","none");
  std::string schur_smoother     = GetParams().get<std::string>("schur smoother","none");
  if (predictor_smoother=="none")
    dserror("The field \"predictor smoother\" is mandatory for the SIMPLE smoother. Fix your xml file");
  if (schur_smoother=="none")
    dserror("The field \"schur smoother\" is mandatory for the SIMPLE smoother. Fix your xml file");
  if (schur_smoother=="REUSE_MUELU_SMOOTHER" or schur_smoother=="REUSE_MUELU_AMG")
    dserror("Invalid smoother for the schur block. We cannot reuse the smoothers generated by Muelu.");

  //other params
  int iter = GetParams().get<int>("sweeps",3);
  double omega = GetParams().get<double>("omega",1.0);
  double alpha = GetParams().get<double>("alpha",1.0);
  double beta = GetParams().get<double>("beta",1.0);
  std::string inverse_method = GetParams().get<std::string>
    ("predictor inverse","row sums diagonal blocks");
  std::string correction = GetParams().get<std::string>("correction","approximated inverse");
  std::string algorithm = GetParams().get<std::string>("algorithm","simple");

  // =============================================================
  // Some output
  // =============================================================

  if (GetVerbosity()=="on")
  {
    std::cout << std::endl;
    std::cout << "Creating a SIMPLE smoother for blocks (" ;
    for(size_t i=0;i<GetBlocks().size();i++)
    {
      std::cout << GetBlocks()[i] ;
      if (i<GetBlocks().size()-1)
        std::cout << ", ";
      else
        std::cout << ") ";
    }
    std::cout << " at level " << GetLevel() << std::endl;
    std::cout << "The chosen parameters are" << std::endl;
    std::cout << "algorithm = " << algorithm << std::endl;
    std::cout << "predictor block = " ;
    std::cout << "(";
    for(size_t j=0;j<SuperBlocks2Blocks[pred].size();j++)
    {
      std::cout << SuperBlocks2Blocks[pred][j];
      if(j<(SuperBlocks2Blocks[pred].size()-1))
        std::cout << ",";
    }
    std::cout << ")" << std::endl;
    std::cout << "predictor smoother = " << predictor_smoother << std::endl;
    std::cout << "predictor inverse = " << inverse_method << std::endl;
    std::cout << "schur block = " ;
    std::cout << "(";
    for(size_t j=0;j<SuperBlocks2Blocks[schur].size();j++)
    {
      std::cout << SuperBlocks2Blocks[schur][j];
      if(j<(SuperBlocks2Blocks[schur].size()-1))
        std::cout << ",";
    }
    std::cout << ")" << std::endl;
    std::cout << "schur smoother = " << schur_smoother << std::endl;
    std::cout << "sweeps = " << iter << std::endl;
    std::cout << "omega = " << omega << std::endl;
    std::cout << "alpha = " << alpha << std::endl;
    std::cout << "beta = " << beta << std::endl;
    if (algorithm == "simple")
      std::cout << "correction = " << correction << std::endl;
    //std::cout << std::endl;
  }

  // =============================================================
  // Rearrange blocks
  // =============================================================

  Teuchos::RCP<BlockSparseMatrixBase> OpBlocked = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(GetOperator());
  if (OpBlocked == Teuchos::null)
    dserror("I expect a block matrix here");
  BlockAggrupator myBlockAggrupator(OpBlocked,SuperBlocks2BlocksLocal);
  Teuchos::RCP<SparseOperator> App = myBlockAggrupator.GetSuperBlock(pred,pred);
  Teuchos::RCP<SparseOperator> Aps = myBlockAggrupator.GetSuperBlock(pred,schur);
  Teuchos::RCP<SparseOperator> Asp = myBlockAggrupator.GetSuperBlock(schur,pred);
  Teuchos::RCP<SparseOperator> Ass = myBlockAggrupator.GetSuperBlock(schur,schur);

  //{
  //  Epetra_Map myRange  = App->OperatorRangeMap();
  //  Epetra_Map myDomain = App->OperatorDomainMap();
  //  std::cout << "Matrix App" << std::endl;
  //  std::cout << "   Range   MinAllGID = " << myRange.MinAllGID()  << std::endl;
  //  std::cout << "   Range   MaxAllGID = " << myRange.MaxAllGID()  << std::endl;
  //  std::cout << "   Domain  MinAllGID = " << myDomain.MinAllGID()  << std::endl;
  //  std::cout << "   Domain  MaxAllGID = " << myDomain.MaxAllGID()  << std::endl;
  //}
  //{
  //  Epetra_Map myRange  = Ass->OperatorRangeMap();
  //  Epetra_Map myDomain = Ass->OperatorDomainMap();
  //  std::cout << "Matrix App" << std::endl;
  //  std::cout << "   Range   MinAllGID = " << myRange.MinAllGID()  << std::endl;
  //  std::cout << "   Range   MaxAllGID = " << myRange.MaxAllGID()  << std::endl;
  //  std::cout << "   Domain  MinAllGID = " << myDomain.MinAllGID()  << std::endl;
  //  std::cout << "   Domain  MaxAllGID = " << myDomain.MaxAllGID()  << std::endl;
  //}

  // =============================================================
  // Approximate the schur complement
  // =============================================================

  // Approximate the inverse of App
  Teuchos::RCP<SparseOperator>  invApp = Teuchos::null;
  Teuchos::RCP<SparseMatrix> App_sp = Teuchos::rcp_dynamic_cast<SparseMatrix>(App);
  Teuchos::RCP<BlockSparseMatrixBase> App_bl
    = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(App);
  if(App_bl==Teuchos::null)
  {
    if(App_sp==Teuchos::null)
      dserror("Something wrong");
    invApp = ApproximateInverse(*App_sp,inverse_method);
  }
  else
  {
    if(App_bl==Teuchos::null)
      dserror("Something wrong");
    if(inverse_method=="row sums")
      dserror("The \"row sums\" option is not implemented yet for Block matrices");
    if(App_bl->Rows()!=App_bl->Cols())
      dserror("The number of block rows and columns has to be the same!");
    std::vector<Teuchos::RCP<SparseMatrix> > DiagBlocks(App_bl->Rows(),Teuchos::null);
    for(int row=0;row<App_bl->Rows();row++)
      DiagBlocks[row] = ApproximateInverse(App_bl->Matrix(row,row),inverse_method);
    BlockSparseMatrix_Creator myBlockMatrixCreator;
    invApp
      = myBlockMatrixCreator.CreateBlockSparseMatrix(DiagBlocks,App_bl->Rows(),App_bl->Rows(),View);
  }

  // Compute the schur complement
  Teuchos::RCP<SparseOperator> S = ComputeSchurComplement(invApp,Aps,Asp,Ass);

  // =============================================================
  // Construct smoothers for diagonal superblocks
  // =============================================================

  AMGnxn_SmootherFactory mySmootherCreator;
  mySmootherCreator.SetParamsSmoother (GetParamsSmoother());
  mySmootherCreator.SetHierarchies    (GetHierarchies());
  mySmootherCreator.SetLevel          (GetLevel());
  mySmootherCreator.SetVerbosity      (GetVerbosity());

  // For predictor
  if (SuperBlocks2Blocks[pred].size()==1)
  {
      int thisblock = SuperBlocks2Blocks[pred][0];
      mySmootherCreator.SetBlock          (thisblock);
      mySmootherCreator.SetNullSpace      (GetNullSpaceAllBlocks()[thisblock]);
  }
  else
  {
      mySmootherCreator.SetBlocks (SuperBlocks2Blocks[pred]);
      mySmootherCreator.SetNullSpaceAllBlocks(GetNullSpaceAllBlocks());
  }
  mySmootherCreator.SetSmootherName   (predictor_smoother);
  mySmootherCreator.SetOperator       (App);
  Teuchos::RCP<AMGnxn_SmootherBase>  Smoother_App = mySmootherCreator.Create();

  // For schur
  if (SuperBlocks2Blocks[schur].size()==1)
  {
      int thisblock = SuperBlocks2Blocks[schur][0];
      mySmootherCreator.SetBlock          (thisblock);
      mySmootherCreator.SetNullSpace      (GetNullSpaceAllBlocks()[thisblock]);
  }
  else
  {
    mySmootherCreator.SetBlocks         (SuperBlocks2Blocks[schur]);
    mySmootherCreator.SetNullSpaceAllBlocks(GetNullSpaceAllBlocks());
  }

  mySmootherCreator.SetSmootherName   (schur_smoother);
  mySmootherCreator.SetOperator       (S);
  Teuchos::RCP<AMGnxn_SmootherBase>  Smoother_S = mySmootherCreator.Create();




  // =============================================================
  // Construct SIMPLE smoother
  // =============================================================

   return Teuchos::rcp(new SIMPLE_BlockSmoother(
         App,
         Ass,
         Aps,
         Asp,
         invApp,
         S,
         Smoother_App,
         Smoother_S,
         myBlockAggrupator.GetRangeMapExtractor(),
         myBlockAggrupator.GetDomainMapExtractor(),
         pred,
         schur,
         iter,
         omega,
         alpha,
         algorithm,
         correction));

  //Teuchos::RCP<MergeAndSolve> Skk = Teuchos::rcp(new MergeAndSolve());
  //Skk->Setup(OpBlocked);
  //Teuchos::RCP<AMGnxn_SmootherBase> Skkbase = Teuchos::rcp_dynamic_cast<AMGnxn_SmootherBase>(Skk);
  //return Skkbase;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SparseOperator>
LINALG::SOLVER::SIMPLE_BlockSmootherFactory::ComputeSchurComplement(
    Teuchos::RCP<SparseOperator> invApp,
    Teuchos::RCP<SparseOperator> Aps,
    Teuchos::RCP<SparseOperator> Asp,
    Teuchos::RCP<SparseOperator> Ass)
{

  Teuchos::RCP<SparseOperator> S = Teuchos::null;

  Teuchos::RCP<SparseMatrix> invApp_sp = Teuchos::rcp_dynamic_cast<SparseMatrix>(invApp);
  Teuchos::RCP<BlockSparseMatrixBase> invApp_bl = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(invApp);
  Teuchos::RCP<SparseMatrix> Ass_sp = Teuchos::rcp_dynamic_cast<SparseMatrix>(Ass);
  Teuchos::RCP<BlockSparseMatrixBase> Ass_bl = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Ass);


  if (invApp_sp != Teuchos::null and invApp_bl == Teuchos::null)
  {
    if (Ass_sp != Teuchos::null and Ass_bl == Teuchos::null)
    {
      Teuchos::RCP<SparseMatrix> Asp_sp = Teuchos::rcp_dynamic_cast<SparseMatrix>(Asp);
      Teuchos::RCP<SparseMatrix> Aps_sp = Teuchos::rcp_dynamic_cast<SparseMatrix>(Aps);
      Teuchos::RCP<SparseMatrix> temp = LINALG::MLMultiply(*Asp_sp,*invApp_sp, true);
      Teuchos::RCP<SparseMatrix> S_sp = LINALG::MLMultiply(*temp,*Aps_sp, true);
      S_sp->Add(*Ass_sp,false,1.0,-1.0);
      S_sp->Complete();
      S = Teuchos::rcp_dynamic_cast<SparseOperator>(S_sp);
    }
    else if (Ass_sp == Teuchos::null and Ass_bl != Teuchos::null)
    {
      dserror("TODO: Branch not implemented yet");
    }
    else dserror("Something went wrong");
  }
  else if (invApp_sp == Teuchos::null and invApp_bl != Teuchos::null)
  {
    if (Ass_sp != Teuchos::null and Ass_bl == Teuchos::null)
    {
      Teuchos::RCP<BlockSparseMatrixBase> Asp_bl = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Asp);
      Teuchos::RCP<BlockSparseMatrixBase> Aps_bl = Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Aps);
      int NumBlocks_pp = invApp_bl->Rows();
      Teuchos::RCP<SparseMatrix> S_sp = Teuchos::null;
      for(int b=0;b<NumBlocks_pp;b++)
      {
        Teuchos::RCP<SparseMatrix> temp = LINALG::MLMultiply(Asp_bl->Matrix(0,b),invApp_bl->Matrix(b,b),true);
        Teuchos::RCP<SparseMatrix> S_sp_tmp = LINALG::MLMultiply(*temp,Aps_bl->Matrix(b,0), true);
        if (b==0) S_sp = S_sp_tmp;
        else S_sp->Add(*S_sp_tmp,false,1.0,1.0);
      }
      S_sp->Add(*Ass_sp,false,1.0,-1.0);
      S_sp->Complete();
      S = Teuchos::rcp_dynamic_cast<SparseOperator>(S_sp);
    }
    else if (Ass_sp == Teuchos::null and Ass_bl != Teuchos::null)
    {
      dserror("TODO: Branch not implemented yet");
    }
    else dserror("Something went wrong");
  }
  else dserror("Something went wrong");


  return S;

}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SparseMatrix>
LINALG::SOLVER::SIMPLE_BlockSmootherFactory::ApproximateInverse(
    const SparseMatrix& A, std::string method)
{

  Teuchos::RCP<Epetra_Vector> invAVector = Teuchos::rcp(new Epetra_Vector(A.RowMap()));
  if(method=="diagonal")
  {
    A.ExtractDiagonalCopy(*invAVector);
    int err = invAVector->Reciprocal(*invAVector);
    if (err) dserror("Epetra_MultiVector::Reciprocal returned %d, are we dividing by 0?",err);
  }
  else if(method=="row sums" or method=="row sums diagonal blocks")
  {
    int err = A.EpetraMatrix()->InvRowSums(*invAVector);
    if (err) dserror("Epetra_CrsMatrix::InvRowSums returned %d, are we dividing by 0?",err);
  }
  else
    dserror("Invalid value for \"predictor inverse\". Fix your xml file.");
  Teuchos::RCP<SparseMatrix> S = Teuchos::rcp(new SparseMatrix(*invAVector));
  S->Complete(A.RowMap(),A.RowMap());
  return S;

}



/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_SmootherBase>
LINALG::SOLVER::DirectSolverWrapperFactory::Create()
{
  // Check input
  if (not IsSetOperator())
    dserror("IsSetOperator() returns false");

  if (GetVerbosity()=="on")
  {
    std::cout << std::endl;
    std::cout << "Creating a DIRECT_SOLVER for block " << GetBlock();
    std::cout << " at level " << GetLevel() << std::endl;
  }

  Teuchos::RCP<DirectSolverWrapper> S = Teuchos::rcp(new DirectSolverWrapper);
  Teuchos::RCP<SparseMatrix>  matrix
    = Teuchos::rcp_dynamic_cast<SparseMatrix>(GetOperator());
  if(matrix==Teuchos::null)
    dserror("We expect here a sparse matrix");
  S->Setup(matrix);

  Teuchos::RCP<AMGnxn_SmootherBase> SBase = Teuchos::rcp_dynamic_cast<AMGnxn_SmootherBase>(S);
  if (SBase==Teuchos::null)
    dserror("Something went wrong");

  return S;
}


#endif // HAVE_MueLu
