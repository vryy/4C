/*!----------------------------------------------------------------------
\file solver_amgnxn_hierarchies.cpp

<pre>
Maintainer: Francesc Verdugo
            verdugo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Aug 11, 2014
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
#include "solver_amgnxn_hierarchies.H"


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

LINALG::SOLVER::AMGnxn_Hierarchies::AMGnxn_Hierarchies(
    Teuchos::RCP<BlockSparseMatrixBase> A,
    std::vector<Teuchos::ParameterList> muelu_params,
    std::vector<int> num_pdes,
    std::vector<int> null_spaces_dim,
    std::vector<Teuchos::RCP<std::vector<double> > > null_spaces_data,
    int NumLevelAMG,
    std::string verbosity):
A_                (A               ),
muelu_params_     (muelu_params    ),
num_pdes_         (num_pdes        ),
null_spaces_dim_  (null_spaces_dim ),
null_spaces_data_ (null_spaces_data),
NumBlocks_        (A->Rows()       ),
NumLevelAMG_      (NumLevelAMG     ),
verbosity_        (verbosity       )
{
  // Plausibility checks
  if (A_->Rows() != A_->Cols() or
      (int)(muelu_params_.size()) != NumBlocks_ or
      (int)( num_pdes_.size()) != NumBlocks_ or
      (int)( null_spaces_dim_.size()) != NumBlocks_ or
      (int)( null_spaces_data_.size()) != NumBlocks_ )
    dserror("Something wrong");

  // Setput
  Setup();
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGnxn_Hierarchies::GetNumLevels(int block)
{
  return H_block_[block]->GetNumLevels();
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<LINALG::SparseMatrix> >
LINALG::SOLVER::AMGnxn_Hierarchies::GetA(int block)
{
  return A_block_level_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<LINALG::SparseMatrix> >
LINALG::SOLVER::AMGnxn_Hierarchies::GetP(int block)
{
  return P_block_level_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<LINALG::SparseMatrix> >
LINALG::SOLVER::AMGnxn_Hierarchies::GetR(int block)
{
  return R_block_level_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<LINALG::SOLVER::MueluSmootherWrapper> >
LINALG::SOLVER::AMGnxn_Hierarchies::GetSPre(int block)
{
  return SPre_block_level_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<LINALG::SOLVER::MueluSmootherWrapper> >
LINALG::SOLVER::AMGnxn_Hierarchies::GetSPos(int block)
{
  return SPos_block_level_[block];
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGnxn_Hierarchies::GetNumPDEs(int block)
{
  return num_pdes_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGnxn_Hierarchies::GetNullSpaceDim(int block)
{
  return null_spaces_dim_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<std::vector<double> >
LINALG::SOLVER::AMGnxn_Hierarchies::GetNullSpaceData(int block)
{
  return null_spaces_data_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_Hierarchies::Setup()
{

  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::AMGnxn_Hierarchies::Setup");

  // ===========================================================
  // Build up MueLu Hierarchies of each one of the blocks
  // ===========================================================

  H_block_.assign(NumBlocks_,Teuchos::null);
  std::vector<int> offsets(NumLevelAMG_-1,0);
  for(int block=0;block<NumBlocks_;block++)
  {
    int offsetFineLevel =  A_->Matrix(block,block).RowMap().MinAllGID();
    Teuchos::RCP<Epetra_Operator> A_eop = A_->Matrix(block,block).EpetraOperator();
    H_block_[block]=BuildMueLuHierarchy(
        muelu_params_[block],
        num_pdes_[block],
        null_spaces_dim_[block],
        null_spaces_data_[block],
        A_eop,
        block,
        NumBlocks_,
        offsets,
        offsetFineLevel);
  }

  // ===========================================================
  // Determine number of levels
  // ===========================================================

  NumLevelMax_ = -10000000;
  NumLevelMin_ =  10000000;
  for(int block=0;block<NumBlocks_;block++)
  {
    int NumLevel_this_block = H_block_[block]->GetNumLevels();
    if (NumLevel_this_block > NumLevelMax_)
      NumLevelMax_ = NumLevel_this_block;
    if (NumLevel_this_block < NumLevelMin_)
      NumLevelMin_ = NumLevel_this_block;
  }


  // ===========================================================
  // Extract matrices, transfer operators and smoothers from the hierarchies
  // ===========================================================

  for(int block=0;block<NumBlocks_;block++)
  {

    int NumLevel_this_block = H_block_[block]->GetNumLevels();
    std::vector<Teuchos::RCP<SparseMatrix> > A_level(NumLevel_this_block,Teuchos::null);
    std::vector<Teuchos::RCP<SparseMatrix> > P_level(NumLevel_this_block-1,Teuchos::null);
    std::vector<Teuchos::RCP<SparseMatrix> > R_level(NumLevel_this_block-1,Teuchos::null);
    std::vector<Teuchos::RCP<MueluSmootherWrapper> > SPre_level(NumLevel_this_block,Teuchos::null);
    std::vector<Teuchos::RCP<MueluSmootherWrapper> > SPos_level(NumLevel_this_block,Teuchos::null);

    Teuchos::RCP<Matrix>              myA = Teuchos::null;
    Teuchos::RCP<Epetra_CrsMatrix> myAcrs = Teuchos::null;
    Teuchos::RCP<SparseMatrix>     myAspa = Teuchos::null;
    Teuchos::RCP<SmootherBase>        myS = Teuchos::null;
    Teuchos::RCP<LINALG::SOLVER::MueluSmootherWrapper> mySWrap = Teuchos::null;

    bool explicitdirichlet = A_->Matrix(0,0).ExplicitDirichlet();
    bool savegraph         = A_->Matrix(0,0).SaveGraph();

    for(int level=0;level<NumLevel_this_block;level++)
    {
      Teuchos::RCP<Level> this_level = H_block_[block]->GetLevel(level);
      if (this_level->IsAvailable("A"))
      {
        myA    = this_level->Get< Teuchos::RCP<Matrix> >("A");
        myAcrs = MueLu::Utils<double,int,int,Node>::Op2NonConstEpetraCrs(myA);
        myAspa = Teuchos::rcp(new SparseMatrix(myAcrs,Copy,explicitdirichlet,savegraph));
        A_level[level] = myAspa;
      }
      else
        dserror("Error in extracting A");

      if (this_level->IsAvailable("PreSmoother"))
      {
        myS     = this_level->Get< Teuchos::RCP<SmootherBase> >("PreSmoother");
        mySWrap = Teuchos::rcp(new LINALG::SOLVER::MueluSmootherWrapper(myS));
        SPre_level[level]=mySWrap;
      }
      else
        dserror("Error in extracting PreSmoother");

      if(level<NumLevel_this_block-1)
      {
        if (this_level->IsAvailable("PostSmoother"))
        {
          myS     = this_level->Get< Teuchos::RCP<SmootherBase> >("PostSmoother");
          mySWrap = Teuchos::rcp(new LINALG::SOLVER::MueluSmootherWrapper(myS));
          SPos_level[level]=mySWrap;
        }
        else
          dserror("Error in extracting PostSmoother");
      }

      if(level!=0)
      {

        if (this_level->IsAvailable("P"))
        {
          myA    = this_level->Get< Teuchos::RCP<Matrix> >("P");
          myAcrs = MueLu::Utils<double,int,int,Node>::Op2NonConstEpetraCrs(myA);
          myAspa = Teuchos::rcp(new SparseMatrix(myAcrs,Copy,explicitdirichlet,savegraph));
          P_level[level-1]=myAspa;
        }
        else
          dserror("Error in extracting P");

        if (this_level->IsAvailable("R"))
        {
          myA = this_level->Get< Teuchos::RCP<Matrix> >("R");
          myAcrs =MueLu::Utils<double,int,int,Node>::Op2NonConstEpetraCrs(myA);
          myAspa = Teuchos::rcp(new SparseMatrix(myAcrs,Copy,explicitdirichlet,savegraph));
          R_level[level-1]=myAspa;
        }
        else
          dserror("Error in extracting R");
      }

    } // loop in levels
    A_block_level_.push_back(A_level);
    P_block_level_.push_back(P_level);
    R_block_level_.push_back(R_level);
    SPre_block_level_.push_back(SPre_level);
    SPos_block_level_.push_back(SPos_level);

  } // loop in blocks

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<Hierarchy> LINALG::SOLVER::AMGnxn_Hierarchies::BuildMueLuHierarchy(
    Teuchos::ParameterList paramListFromXml,
    int numdf,
    int dimns,
    Teuchos::RCP<std::vector<double> > nsdata,
    Teuchos::RCP<Epetra_Operator> A_eop,
    int block,
    int NumBlocks,
    std::vector<int>& offsets,
    int offsetFineLevel)
{

  //Some cheks
  if(numdf<1 or dimns<1)
    dserror("Error: PDE equations or null space dimension wrong.");
  if(nsdata==Teuchos::null)
    dserror("Error: null space data is empty");

  //Prepare operator for MueLu
  Teuchos::RCP<Epetra_CrsMatrix> A_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(A_eop);
  if(A_crs==Teuchos::null)
    dserror("Make sure that the input matrix is a Epetra_CrsMatrix (or derived)");
  Teuchos::RCP<CrsMatrix> mueluA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A_crs));
  Teuchos::RCP<CrsMatrixWrap> mueluA_wrap = Teuchos::rcp(new CrsMatrixWrap(mueluA));
  Teuchos::RCP<Matrix> mueluOp = Teuchos::rcp_dynamic_cast<Matrix>(mueluA_wrap);
  mueluOp->SetFixedBlockSize(numdf,offsetFineLevel);

  // Prepare null space vector for MueLu
  Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > rowMap = mueluA->getRowMap();
  Teuchos::RCP<MultiVector> nspVector =
    Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap,dimns,true);
  for ( size_t i=0; i < Teuchos::as<size_t>(dimns); i++) {
    Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
    const size_t myLength = nspVector->getLocalLength();
    for(size_t j=0; j<myLength; j++) {
      nspVectori[j] = (*nsdata)[i*myLength+j];
    }
  }

  // Build up hierarchy
  Teuchos::RCP<Hierarchy> H = Teuchos::null;


  // Add information about maps offsets
  std::string offsets_str("{");
  for(int i=0;i<(int)offsets.size();i++)
  {
    offsets_str= offsets_str + ConvertInt(offsets[i]);
    if(i<(int)offsets.size()-1)
      offsets_str= offsets_str + ", ";
  }
  offsets_str= offsets_str + "}";
  if(paramListFromXml.sublist("Factories",true).isSublist("myCoarseMapFactory123"))
    dserror("We are going to overwrite the factory 'myCoarseMapFactory123'. Please use an other name");
  Teuchos::ParameterList& myCoarseMapFactoryList =
    paramListFromXml.sublist("Factories",true).sublist("myCoarseMapFactory123");
  myCoarseMapFactoryList.set("factory","CoarseMapFactory");
  myCoarseMapFactoryList.set("Domain GID offsets",offsets_str);
  if(paramListFromXml.sublist("Hierarchy",true).sublist("All").isParameter("CoarseMap"))
    dserror("We are going to overwrite 'CoarseMap'. Don't use 'CoarseMap' here.");
  Teuchos::ParameterList& AllList =
    paramListFromXml.sublist("Hierarchy").sublist("All");
  AllList.set("CoarseMap","myCoarseMapFactory123");

  // Add offset for the finest level
  Teuchos::ParameterList& MatrixList = paramListFromXml.sublist("Matrix");
  MatrixList.set<int>("DOF offset",offsetFineLevel);
  MatrixList.set<int>("number of equations",numdf);

  // Build up hierarchy
  ParameterListInterpreter mueLuFactory(paramListFromXml);
  H = mueLuFactory.CreateHierarchy();
  H->SetDefaultVerbLevel(MueLu::Extreme); // TODO sure?
  H->GetLevel(0)->Set("A", mueluOp);
  H->GetLevel(0)->Set("Nullspace", nspVector);
  H->GetLevel(0)->setlib(Xpetra::UseEpetra);
  H->setlib(Xpetra::UseEpetra);
  mueLuFactory.SetupHierarchy(*H);

  // Recover information about the maps
  int NumLevel_block = H->GetNumLevels();
  Teuchos::RCP<Level>        this_level = Teuchos::null;
  Teuchos::RCP<Matrix>              myA = Teuchos::null;
  Teuchos::RCP<Epetra_CrsMatrix> myAcrs = Teuchos::null;
  for(int level=1;(level<NumLevel_block) and (level<(int)offsets.size()+1);level++)
  {
    this_level=H->GetLevel(level);
    if (this_level->IsAvailable("A"))
    {
      myA    = this_level->Get< Teuchos::RCP<Matrix> >("A"); //Matrix
      myAcrs = MueLu::Utils<double,int,int,Node>::Op2NonConstEpetraCrs(myA);
    }
    else
      dserror("Error in extracting A");

    offsets[level-1] =  offsets[level-1] + myAcrs->RangeMap().MaxAllGID() + 1;
  }

  return H;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGnxn_Hierarchies::GetNumLevelMin()
{
  return NumLevelMin_;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGnxn_Hierarchies::GetNumBlocks()
{
  return NumBlocks_;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::BlockSparseMatrixBase> LINALG::SOLVER::AMGnxn_Hierarchies::GetBlockMatrix()
{
  if (A_== Teuchos::null)
    dserror("No data");
  return A_;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<Hierarchy> LINALG::SOLVER::AMGnxn_Hierarchies::GetH(int block)
{
  Teuchos::RCP<Hierarchy> H = H_block_[block];
  if(H == Teuchos::null)
    dserror("No data");
  return H;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SparseMatrix> LINALG::SOLVER::AMGnxn_Hierarchies::GetA(int block, int level)
{
  Teuchos::RCP<SparseMatrix> A = A_block_level_[block][level];
  if(A == Teuchos::null)
    dserror("No data");
  return A;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SparseMatrix> LINALG::SOLVER::AMGnxn_Hierarchies::GetP(int block, int level)
{
  Teuchos::RCP<SparseMatrix> P = P_block_level_[block][level];
  if(P == Teuchos::null)
    dserror("No data");
  return P;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SparseMatrix> LINALG::SOLVER::AMGnxn_Hierarchies::GetR(int block, int level)
{
  Teuchos::RCP<SparseMatrix> R = R_block_level_[block][level];
  if(R == Teuchos::null)
    dserror("No data");
  return R;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::MueluSmootherWrapper> LINALG::SOLVER::AMGnxn_Hierarchies::GetSPre
(int block, int level)
{
  Teuchos::RCP<MueluSmootherWrapper> SPre = SPre_block_level_[block][level];
  if(SPre == Teuchos::null)
    dserror("No data");
  return SPre;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::MueluSmootherWrapper> LINALG::SOLVER::AMGnxn_Hierarchies::GetSPos
(int block, int level)
{
  Teuchos::RCP<MueluSmootherWrapper> SPos = SPos_block_level_[block][level];
  if(SPos == Teuchos::null)
    dserror("No data");
  return SPos;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

LINALG::SOLVER::AMGnxn_MonolithicHierarchy::AMGnxn_MonolithicHierarchy(
    Teuchos::RCP<AMGnxn_Hierarchies> H,
    const Teuchos::ParameterList& params,
    const Teuchos::ParameterList& params_smoothers):
  H_(H),
  params_(params),
  params_smoothers_(params_smoothers)
{

  // Expected parameters in params (example)
  //<ParameterList name="params">
  //
  //  <Parameter name="number of levels"                 type="int"  value="..."/>
  //
  //  <Parameter name="smoother: all but coarsest level" type="string"  value="myFinestSmoother"/>
  //
  //  <Parameter name="smoother: coarsest level"         type="string"  value="myCoarsestSmoother"/>
  //
  //  <Parameter name="verbosity"                        type="string"  value="on"/>
  //
  //</ParameterList>

  // Expected parameters in params_smoothers (example)
  //<ParameterList name="params_smoothers">
  //
  //  <ParameterList name="myFinestSmoother">
  //
  //   ...    ...    ...    ...    ...
  //
  //  </ParameterList>
  //
  //  <ParameterList name="myCoarsestSmoother">
  //
  //   ...    ...    ...    ...    ...
  //
  //  </ParameterList>
  //
  //</ParameterList>

  Setup();
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGnxn_MonolithicHierarchy::GetNumLevels()
{
  return NumLevels_;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGnxn_MonolithicHierarchy::Setup()
{

  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::AMGnxn_MonolithicHierarchy");

  // ====================================================
  // Create block transfer operators
  // ====================================================

  NumLevels_ = params_.get<int>("number of levels",-1);
  if(NumLevels_==-1)
    dserror("Missing \"number of levels\" in your xml file");
  NumLevels_ = std::min(NumLevels_,H_->GetNumLevelMin());
  NumBlocks_ = H_->GetNumBlocks();


  // ====================================================
  // Some output
  // ====================================================
  std::string verbosity = params_.get<std::string>("verbosity","off");
  if (verbosity=="on")
  {

    //std::cout << "===============================================" << std::endl;
    //std::cout << "AMGnxn preconditioner: debug info  (begin)" << std::endl;
    //std::cout << std::endl;
    //std::cout << "===============================================" << std::endl;
    std::cout << "number of blocks = " << NumBlocks_ << std::endl;
    std::cout << "number of levels = " << NumLevels_ << std::endl;
    for(int i=0;i<NumBlocks_;i++)
      std::cout << "block " << i << ": number of levels = "
        << GetHierarchies()->GetH(i)->GetNumLevels() << std::endl;
  }



  // ====================================================
  // Create block transfer operators
  // ====================================================

  BlockSparseMatrix_Creator myBlockMatrixCreator;
  P_.assign(NumLevels_-1,Teuchos::null);
  R_.assign(NumLevels_-1,Teuchos::null);
  for(int level=0;level<NumLevels_-1;level++)
  {
    std::vector< Teuchos::RCP<SparseMatrix> > Pblocks(NumBlocks_,Teuchos::null);
    std::vector< Teuchos::RCP<SparseMatrix> > Rblocks(NumBlocks_,Teuchos::null);
    for(int block=0;block<NumBlocks_;block++)
    {
      Pblocks[block]=H_->GetP(block,level);
      Rblocks[block]=H_->GetR(block,level);
    }
    P_[level] = myBlockMatrixCreator.CreateBlockSparseMatrix(Pblocks,NumBlocks_,NumBlocks_,Copy);
    R_[level] = myBlockMatrixCreator.CreateBlockSparseMatrix(Rblocks,NumBlocks_,NumBlocks_,Copy);
  }

  // ====================================================
  // Create coarser matrices
  // ====================================================

  A_.assign(NumLevels_,Teuchos::null);
  for(int level=0;level<NumLevels_;level++)
  {
    if(level==0)
      A_[level] = H_->GetBlockMatrix();
    else
    {
      std::vector< Teuchos::RCP<SparseMatrix> > Ablocks(NumBlocks_*NumBlocks_,Teuchos::null);
      for(int block=0;block<NumBlocks_;block++)
        Ablocks[block*NumBlocks_+block] = H_->GetA(block,level);
      for(int row=0;row<NumBlocks_;row++)
        for(int col=0;col<NumBlocks_;col++)
        {
          if(row!=col) // The diagonal blocks have been already assigned
          {
            const SparseMatrix& A_spa   = A_[level-1]->Matrix(row,col);
            const SparseMatrix& P_spa   = P_[level-1]->Matrix(col,col);
            const SparseMatrix& R_spa   = R_[level-1]->Matrix(row,row);
            Teuchos::RCP<SparseMatrix>  AP_spa  = Teuchos::null;
            AP_spa  = LINALG::Multiply(A_spa,false,P_spa,false,true);
            if(AP_spa==Teuchos::null)
              dserror("Error in AP");
            Teuchos::RCP<SparseMatrix>  RAP_spa = Teuchos::null;
            RAP_spa = LINALG::Multiply(R_spa,false,*AP_spa,false,true);
            if(RAP_spa==Teuchos::null)
              dserror("Error in RAP");
            Ablocks[row*NumBlocks_+col] = RAP_spa;
          }
        }
      A_[level] = myBlockMatrixCreator.CreateBlockSparseMatrix(Ablocks,NumBlocks_,NumBlocks_,Copy);
    }
  }

  // ====================================================
  // Create smoothers
  // ====================================================

  Spre_.assign(NumLevels_,Teuchos::null);
  Spos_.assign(NumLevels_-1,Teuchos::null);
  for(int level=0;level<NumLevels_;level++)
  {

    if (level<NumLevels_-1)
    {
      Spre_[level] = BuildSmoother(level);
      Spos_[level] = Spre_[level]; //BuildSmoother(level);
    }
    else
    {
      Spre_[level] = BuildSmoother(level);
    }
  }

  //if (verbosity=="on")
  //{
  //  //std::cout << "===============================================" << std::endl;
  //  std::cout << std::endl;
  //  std::cout << "AMGnxn preconditioner: debug info  (end)" << std::endl;
  //  std::cout << "===============================================" << std::endl;
  //}

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::BlockSparseMatrixBase>
LINALG::SOLVER::AMGnxn_MonolithicHierarchy::GetA(int level)
{return A_[level];}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_Hierarchies>
LINALG::SOLVER::AMGnxn_MonolithicHierarchy::GetHierarchies()
{return H_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGnxn_BlockSmootherBase>
LINALG::SOLVER::AMGnxn_MonolithicHierarchy::BuildSmoother(int level)
{

  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::AMGnxn_MonolithicHierarchy::BuildSmoother");

  std::string smother_name;
  if (level<NumLevels_-1)
  {
    smother_name = params_.get<std::string>("smoother: all but coarsest level","BGS");
    //if (smother_name=="BGS")
    //  std::cout << "WARNING file: " << __FILE__ << " line: "<< __LINE__ << " :: Using default values for \"smoother: all but coarsest level\" = BGS" << std::endl;
  }
  else
  {
    smother_name = params_.get<std::string>("smoother: coarsest level","BGS");
    //if (smother_name=="BGS")
    //  std::cout << "WARNING file: " << __FILE__ << " line: "<< __LINE__ << " :: Using default values for \"smoother: all but coarsest level\" = BGS" << std::endl;
  }

  std::string verbosity = params_.get<std::string>("verbosity","off");

  std::vector<int> blocks(H_->GetNumBlocks(),0);
  for(int i=0;i<H_->GetNumBlocks();i++)
    blocks[i]=i;

  AMGnxn_SmootherFactory mySmootherCreator;
  mySmootherCreator.SetOperator(GetA(level));
  mySmootherCreator.SetParamsSmoother(params_smoothers_);
  mySmootherCreator.SetHierarchies(GetHierarchies());
  mySmootherCreator.SetLevel(level);
  mySmootherCreator.SetBlocks(blocks);
  mySmootherCreator.SetSmootherName(smother_name);
  mySmootherCreator.SetVerbosity(verbosity);

  // Recover null spaces from the hierarchies and give it to the smoother creator
  std::vector<NullSpaceInfo> null_space_blocks;
  for(int i=0;i<NumBlocks_;i++)
  {
    NullSpaceInfo myNS(H_->GetNumPDEs(i),H_->GetNullSpaceDim(i),H_->GetNullSpaceData(i));
    null_space_blocks.push_back(myNS);
  }
  mySmootherCreator.SetNullSpaceAllBlocks(null_space_blocks);

  Teuchos::RCP<AMGnxn_SmootherBase> Sbase = mySmootherCreator.Create();
  Teuchos::RCP<AMGnxn_BlockSmootherBase>
    Sblock = Teuchos::rcp_dynamic_cast<AMGnxn_BlockSmootherBase>(Sbase);
  if(Sblock == Teuchos::null)
    dserror("Something wrong. Fix the xml file defining the smoother");
  return Sblock;

}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::Richardson_Vcycle_Operator> LINALG::SOLVER::AMGnxn_MonolithicHierarchy::BuildVCycle()
{
  int NumSweeps = 1; // Hard coded
  double omega = 1.0;// Hard coded
  Teuchos::RCP<Richardson_Vcycle_Operator> V = Teuchos::rcp(
      new  Richardson_Vcycle_Operator (NumLevels_,NumSweeps,omega));

  std::vector<Teuchos::RCP<Epetra_Operator> > Aop(A_.size(),Teuchos::null);
  for(int i=0;i< (int) Aop.size();i++)
    Aop[i] = Teuchos::rcp_dynamic_cast<Epetra_Operator>(A_[i]);
  V->SetOperators(Aop);

  std::vector<Teuchos::RCP<Epetra_Operator> > Pop(P_.size(),Teuchos::null);
  for(int i=0;i< (int) Pop.size();i++)
    Pop[i] = Teuchos::rcp_dynamic_cast<Epetra_Operator>(P_[i]);
  V->SetProjectors(Pop);

  std::vector<Teuchos::RCP<Epetra_Operator> > Rop(R_.size(),Teuchos::null);
  for(int i=0;i< (int) Rop.size();i++)
    Rop[i] = Teuchos::rcp_dynamic_cast<Epetra_Operator>(R_[i]);
  V->SetRestrictors(Rop);

  std::vector<Teuchos::RCP<AMGnxn_SmootherBase> > SpreOp(Spre_.size(),Teuchos::null);
  for(int i=0;i< (int) SpreOp.size();i++)
    SpreOp[i] = Teuchos::rcp_dynamic_cast<AMGnxn_SmootherBase>(Spre_[i]);
  V->SetPreSmoothers(SpreOp);

  std::vector<Teuchos::RCP<AMGnxn_SmootherBase> > SposOp(Spos_.size(),Teuchos::null);
  for(int i=0;i< (int) SposOp.size();i++)
    SposOp[i] = Teuchos::rcp_dynamic_cast<AMGnxn_SmootherBase>(Spos_[i]);
  V->SetPosSmoothers(SposOp);
  return V;
}




#endif // HAVE_MueLu
