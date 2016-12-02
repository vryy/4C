/*!----------------------------------------------------------------------
\file solver_amgnxn_hierarchies.cpp

<pre>
\brief Declaration
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15262
Created on: Aug 11, 2014
</pre>
*----------------------------------------------------------------------*/


#ifdef HAVE_MueLu
#include <iostream>

#include <Teuchos_PtrDecl.hpp>
#include <Epetra_Time.h>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <MueLu_MLParameterListInterpreter_decl.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_EpetraOperator.hpp>
#include "EpetraExt_RowMatrixOut.h"
#include "../drt_lib/drt_dserror.H"
#include "solver_amgnxn_hierarchies.H"
#include "solver_amgnxn_vcycle.H"
#include "../linalg/linalg_multiply.H"


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

LINALG::SOLVER::AMGNXN::Hierarchies::Hierarchies(
    Teuchos::RCP<AMGNXN::BlockedMatrix> A,
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
NumBlocks_        (A->GetNumRows() ),
NumLevelAMG_      (NumLevelAMG     ),
verbosity_        (verbosity       )
{
  // Plausibility checks
  if (A_->GetNumRows() != A_->GetNumCols() or
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

int LINALG::SOLVER::AMGNXN::Hierarchies::GetNumLevels(int block)
{
  if (H_block_[block] == Teuchos::null)
    return NumLevelMax_;
  else
    return H_block_[block]->GetNumLevels();
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<LINALG::SparseMatrix> >
LINALG::SOLVER::AMGNXN::Hierarchies::GetA(int block)
{
  return A_block_level_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<LINALG::SparseMatrix> >
LINALG::SOLVER::AMGNXN::Hierarchies::GetP(int block)
{
  return P_block_level_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<LINALG::SparseMatrix> >
LINALG::SOLVER::AMGNXN::Hierarchies::GetR(int block)
{
  return R_block_level_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<LINALG::SOLVER::AMGNXN::MueluSmootherWrapper> >
LINALG::SOLVER::AMGNXN::Hierarchies::GetSPre(int block)
{
  return SPre_block_level_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<LINALG::SOLVER::AMGNXN::MueluSmootherWrapper> >
LINALG::SOLVER::AMGNXN::Hierarchies::GetSPos(int block)
{
  return SPos_block_level_[block];
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGNXN::Hierarchies::GetNumPDEs(int block)
{
  return num_pdes_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGNXN::Hierarchies::GetNullSpaceDim(int block)
{
  return null_spaces_dim_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<std::vector<double> >
LINALG::SOLVER::AMGNXN::Hierarchies::GetNullSpaceData(int block)
{
  return null_spaces_data_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGNXN::Hierarchies::Setup()
{

  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::AMGNXN::Hierarchies::Setup");

  // ===========================================================
  // Build up MueLu Hierarchies of each one of the blocks
  // ===========================================================

  H_block_.assign(NumBlocks_,Teuchos::null);
  std::vector<int> offsets(NumLevelAMG_-1,0);
  for(int block=0;block<NumBlocks_;block++)
  {
    int offsetFineLevel =  A_->GetMatrix(block,block)->RowMap().MinAllGID();
    Teuchos::RCP<Epetra_Operator> A_eop = A_->GetMatrix(block,block)->EpetraOperator();
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
    if (H_block_[block] == Teuchos::null)
      continue;
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

    // create a dummy hierarchy by repeating the same matrix
    if (H_block_[block] == Teuchos::null)
    {

      std::vector<Teuchos::RCP<SparseMatrix> > A_level(NumLevelMax_,Teuchos::null);
      std::vector<Teuchos::RCP<SparseMatrix> > P_level(NumLevelMax_-1,Teuchos::null);
      std::vector<Teuchos::RCP<SparseMatrix> > R_level(NumLevelMax_-1,Teuchos::null);
      std::vector<Teuchos::RCP<AMGNXN::MueluSmootherWrapper> > SPre_level(NumLevelMax_,Teuchos::null);
      std::vector<Teuchos::RCP<AMGNXN::MueluSmootherWrapper> > SPos_level(NumLevelMax_-1,Teuchos::null);

      Teuchos::RCP<SparseMatrix> Abb = A_->GetMatrix(block,block);
      Teuchos::RCP<SparseMatrix> Peye = LINALG::Eye(Abb->DomainMap());
      Teuchos::RCP<SparseMatrix> Reye = LINALG::Eye(Abb->RangeMap());

      for(int level=0;level<NumLevelMax_;level++)
        A_level[level] = Abb;

      for(int level=0;level<NumLevelMax_-1;level++)
      {
        P_level[level] = Peye;
        R_level[level] = Reye;
      }

      A_block_level_.push_back(A_level);
      P_block_level_.push_back(P_level);
      R_block_level_.push_back(R_level);
      SPre_block_level_.push_back(SPre_level);
      SPos_block_level_.push_back(SPos_level);

    }
    else // Recover objects created by muelu
    {

      int NumLevel_this_block = H_block_[block]->GetNumLevels();
      std::vector<Teuchos::RCP<SparseMatrix> > A_level(NumLevel_this_block,Teuchos::null);
      std::vector<Teuchos::RCP<SparseMatrix> > P_level(NumLevel_this_block-1,Teuchos::null);
      std::vector<Teuchos::RCP<SparseMatrix> > R_level(NumLevel_this_block-1,Teuchos::null);
      std::vector<Teuchos::RCP<AMGNXN::MueluSmootherWrapper> > SPre_level(NumLevel_this_block,Teuchos::null);
      std::vector<Teuchos::RCP<AMGNXN::MueluSmootherWrapper> > SPos_level(NumLevel_this_block,Teuchos::null);

      // some local typedefs
      typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
      typedef MueLu::SmootherBase<Scalar,LocalOrdinal,GlobalOrdinal,Node> SmootherBase;

      Teuchos::RCP<Matrix>              myA = Teuchos::null;
      Teuchos::RCP<Epetra_CrsMatrix> myAcrs = Teuchos::null;
      Teuchos::RCP<SparseMatrix>     myAspa = Teuchos::null;
      Teuchos::RCP<SmootherBase>        myS = Teuchos::null;
      Teuchos::RCP<LINALG::SOLVER::AMGNXN::MueluSmootherWrapper> mySWrap = Teuchos::null;

      bool explicitdirichlet = A_->GetMatrix(0,0)->ExplicitDirichlet();
      bool savegraph         = A_->GetMatrix(0,0)->SaveGraph();

      for(int level=0;level<NumLevel_this_block;level++)
      {
        Teuchos::RCP<MueLu::Level> this_level = H_block_[block]->GetLevel(level);
        if (this_level->IsAvailable("A"))
        {
          myA    = this_level->Get< Teuchos::RCP<Matrix> >("A");
          myAcrs = MueLu::Utils<double,int,int,Node>::Op2NonConstEpetraCrs(myA);
          myAspa = Teuchos::rcp(new SparseMatrix(myAcrs,LINALG::Copy,explicitdirichlet,savegraph));
          A_level[level] = myAspa;
        }
        else
          dserror("Error in extracting A");

        if (this_level->IsAvailable("PreSmoother"))
        {
          myS     = this_level->Get< Teuchos::RCP<SmootherBase> >("PreSmoother");
          mySWrap = Teuchos::rcp(new LINALG::SOLVER::AMGNXN::MueluSmootherWrapper(myS));
          SPre_level[level]=mySWrap;
        }
        else
          dserror("Error in extracting PreSmoother");

        if(level<NumLevel_this_block-1)
        {
          if (this_level->IsAvailable("PostSmoother"))
          {
            myS     = this_level->Get< Teuchos::RCP<SmootherBase> >("PostSmoother");
            mySWrap = Teuchos::rcp(new LINALG::SOLVER::AMGNXN::MueluSmootherWrapper(myS));
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
            myAspa = Teuchos::rcp(new SparseMatrix(myAcrs,LINALG::Copy,explicitdirichlet,savegraph));
            P_level[level-1]=myAspa;
          }
          else
            dserror("Error in extracting P");

          if (this_level->IsAvailable("R"))
          {
            myA = this_level->Get< Teuchos::RCP<Matrix> >("R");
            myAcrs =MueLu::Utils<double,int,int,Node>::Op2NonConstEpetraCrs(myA);
            myAspa = Teuchos::rcp(new SparseMatrix(myAcrs,LINALG::Copy,explicitdirichlet,savegraph));
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

    } // else


  } // loop in blocks

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > LINALG::SOLVER::AMGNXN::Hierarchies::BuildMueLuHierarchy(
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

  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::AMGNXN::Hierarchies::BuildMueLuHierarchy");

  Epetra_Time timer(A_eop->Comm());
  timer.ResetStartTime();

  Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > H = Teuchos::null;
  bool create_uncoarsened_hierarchy = paramListFromXml.get<bool>("create un-coarsened hierarchy",false);
  bool fix_coarse_maps = paramListFromXml.get<bool>("fix coarse maps",false); // this is required in all the fields if we want to merge and solve a coarse level block matrix

  if (not create_uncoarsened_hierarchy)
  {


    //Some cheks
    if(numdf<1 or dimns<1)
      dserror("Error: PDE equations or null space dimension wrong.");
    if(nsdata==Teuchos::null)
      dserror("Error: null space data is empty");

    // Hack making TSI work with the trilinos Q1_2015. The Q1_2014 version worked without this
    if(numdf==1) offsetFineLevel = 0;

    //Prepare operator for MueLu
    Teuchos::RCP<Epetra_CrsMatrix> A_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(A_eop);
    if(A_crs==Teuchos::null)
      dserror("Make sure that the input matrix is a Epetra_CrsMatrix (or derived)");
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mueluA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(A_crs));
    Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mueluA_wrap = Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(mueluA));
    Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > mueluOp = Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >(mueluA_wrap);
    mueluOp->SetFixedBlockSize(numdf,offsetFineLevel);

    // Prepare null space vector for MueLu
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > rowMap = mueluA->getRowMap();
    Teuchos::RCP<Xpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > nspVector =
      Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(rowMap,dimns,true);
    for ( size_t i=0; i < Teuchos::as<size_t>(dimns); i++) {
      Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
      const size_t myLength = nspVector->getLocalLength();
      for(size_t j=0; j<myLength; j++) {
        nspVectori[j] = (*nsdata)[i*myLength+j];
      }
    }



    if (fix_coarse_maps)
    {
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

      if(A_eop->Comm().MyPID()==0)
        std::cout << "offsets_str " << offsets_str << std::endl;

    }

    // Add offset for the finest level
    Teuchos::ParameterList& MatrixList = paramListFromXml.sublist("Matrix");
    MatrixList.set<int>("DOF offset",offsetFineLevel);
    MatrixList.set<int>("number of equations",numdf);

    if(A_eop->Comm().MyPID()==0)
    {
      std::cout << "offsetFineLevel " << offsetFineLevel << std::endl;
    }

    // Build up hierarchy
    MueLu::ParameterListInterpreter<Scalar,LocalOrdinal,GlobalOrdinal,Node> mueLuFactory(paramListFromXml);
    H = mueLuFactory.CreateHierarchy();
    H->SetDefaultVerbLevel(MueLu::Extreme); // TODO sure?
    H->GetLevel(0)->Set("A", mueluOp);
    H->GetLevel(0)->Set("Nullspace", nspVector);
    H->GetLevel(0)->setlib(Xpetra::UseEpetra);
    H->setlib(Xpetra::UseEpetra);
    mueLuFactory.SetupHierarchy(*H);

    // Recover information about the maps
    if (fix_coarse_maps)
    {
      int NumLevel_block = H->GetNumLevels();
      Teuchos::RCP<MueLu::Level> this_level = Teuchos::null;
      Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > myA = Teuchos::null;
      Teuchos::RCP<Epetra_CrsMatrix> myAcrs = Teuchos::null;
      for(int level=1;(level<NumLevel_block) and (level<(int)offsets.size()+1);level++)
      {
        this_level=H->GetLevel(level);
        if (this_level->IsAvailable("A"))
        {
          myA    = this_level->Get<Teuchos::RCP<Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > >("A"); //Matrix
          myAcrs = MueLu::Utils<double,int,int,Node>::Op2NonConstEpetraCrs(myA);
        }
        else
          dserror("Error in extracting A");

        offsets[level-1] =  offsets[level-1] + myAcrs->RangeMap().MaxAllGID() + 1; // TODO I think we don't have to overwrite previous result
      }
    }


  }
  else // when we use a dummy hierarchy we sill have to compute the offsets
  {
    if (fix_coarse_maps)
    {
      for(int level=0; level<(int)offsets.size();level++)
        offsets[level] =  offsets[level] + A_eop->OperatorRangeMap().MaxAllGID() + 1; // TODO I think we don't have to overwrite previous result
    }
  }

  double elaptime =  timer.ElapsedTime();
  if(A_eop->Comm().MyPID()==0)
    std::cout <<  "       Calling LINALG::SOLVER::AMGNXN::Hierarchies::BuildMueLuHierarchy takes " << std::setw(16) << std::setprecision(6) << elaptime << " s" << std::endl ;
  return H;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGNXN::Hierarchies::GetNumLevelMin()
{
  return NumLevelMin_;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int LINALG::SOLVER::AMGNXN::Hierarchies::GetNumBlocks()
{
  return NumBlocks_;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGNXN::BlockedMatrix>
LINALG::SOLVER::AMGNXN::Hierarchies::GetBlockMatrix()
{
  if (A_== Teuchos::null)
    dserror("No data");
  return A_;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > LINALG::SOLVER::AMGNXN::Hierarchies::GetH(int block)
{
  Teuchos::RCP<MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> > H = H_block_[block];
  if(H == Teuchos::null)
    dserror("No data");
  return H;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SparseMatrix> LINALG::SOLVER::AMGNXN::Hierarchies::GetA(int block, int level)
{
  Teuchos::RCP<SparseMatrix> A = A_block_level_[block][level];
  if(A == Teuchos::null)
    dserror("No data");
  return A;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SparseMatrix> LINALG::SOLVER::AMGNXN::Hierarchies::GetP(int block, int level)
{
  Teuchos::RCP<SparseMatrix> P = P_block_level_[block][level];
  if(P == Teuchos::null)
    dserror("No data");
  return P;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SparseMatrix> LINALG::SOLVER::AMGNXN::Hierarchies::GetR(int block, int level)
{
  Teuchos::RCP<SparseMatrix> R = R_block_level_[block][level];
  if(R == Teuchos::null)
    dserror("No data");
  return R;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGNXN::MueluSmootherWrapper> LINALG::SOLVER::AMGNXN::Hierarchies::GetSPre
(int block, int level)
{
  Teuchos::RCP<AMGNXN::MueluSmootherWrapper> SPre = SPre_block_level_[block][level];
  if(SPre == Teuchos::null)
    dserror("No data");
  return SPre;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGNXN::MueluSmootherWrapper> LINALG::SOLVER::AMGNXN::Hierarchies::GetSPos
(int block, int level)
{
  Teuchos::RCP<AMGNXN::MueluSmootherWrapper> SPos = SPos_block_level_[block][level];
  if(SPos == Teuchos::null)
    dserror("No data");
  return SPos;
}



/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

LINALG::SOLVER::AMGNXN::MonolithicHierarchy::MonolithicHierarchy(
    Teuchos::RCP<AMGNXN::Hierarchies> H,
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

int LINALG::SOLVER::AMGNXN::MonolithicHierarchy::GetNumLevels()
{
  return NumLevels_;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void LINALG::SOLVER::AMGNXN::MonolithicHierarchy::Setup()
{

  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::AMGNXN::MonolithicHierarchy::Setup()");

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

  if(H_->GetBlockMatrix()->GetMatrix(0,0)->Comm().MyPID() != 0)
    verbosity = "off";

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
        << GetHierarchies()->GetNumLevels(i) << std::endl;
  }



  // ====================================================
  // Create block transfer operators
  // ====================================================

  P_.assign(NumLevels_-1,Teuchos::null);
  R_.assign(NumLevels_-1,Teuchos::null);
  for(int level=0;level<NumLevels_-1;level++)
  {
    P_[level] = Teuchos::rcp(new AMGNXN::DiagonalBlockedMatrix(NumBlocks_));
    R_[level] = Teuchos::rcp(new AMGNXN::DiagonalBlockedMatrix(NumBlocks_));
    for(int block=0;block<NumBlocks_;block++)
    {
      P_[level]->SetMatrix(H_->GetP(block,level),block,block);
      R_[level]->SetMatrix(H_->GetR(block,level),block,block);
    }
  }

  // ====================================================
  // Create coarser matrices
  // ====================================================
  //

  A_.assign(NumLevels_,Teuchos::null);
  for(int level=0;level<NumLevels_;level++)
  {
    if(level==0)
      A_[level] = H_->GetBlockMatrix();
    else
    {

      A_[level] = Teuchos::rcp(new AMGNXN::BlockedMatrix(NumBlocks_,NumBlocks_));

      for(int block=0;block<NumBlocks_;block++)
        A_[level]->SetMatrix(H_->GetA(block,level),block,block);

      for(int row=0;row<NumBlocks_;row++)
        for(int col=0;col<NumBlocks_;col++)
        {
          if(row!=col) // The diagonal blocks have been already assigned
          {

            Teuchos::RCP<SparseMatrix>  A_spa   = A_[level-1]->GetMatrix(row,col);
            Teuchos::RCP<SparseMatrix>  P_spa   = P_[level-1]->GetMatrix(col,col);
            Teuchos::RCP<SparseMatrix>  R_spa   = R_[level-1]->GetMatrix(row,row);

            Teuchos::RCP<SparseMatrix>  AP_spa  = Teuchos::null;
            AP_spa  = LINALG::MLMultiply(*A_spa,*P_spa,true);
            if(AP_spa==Teuchos::null)
              dserror("Error in AP");

            Teuchos::RCP<SparseMatrix>  RAP_spa = Teuchos::null;
            RAP_spa = LINALG::MLMultiply(*R_spa,*AP_spa,true);
            if(RAP_spa==Teuchos::null)
              dserror("Error in RAP");

            A_[level]->SetMatrix(RAP_spa,row,col);

          }
        }
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

Teuchos::RCP<LINALG::SOLVER::AMGNXN::BlockedMatrix>
LINALG::SOLVER::AMGNXN::MonolithicHierarchy::GetA(int level)
{return A_[level];}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGNXN::Hierarchies>
LINALG::SOLVER::AMGNXN::MonolithicHierarchy::GetHierarchies()
{return H_;}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGNXN::GenericSmoother>
LINALG::SOLVER::AMGNXN::MonolithicHierarchy::BuildSmoother(int level)
{

  TEUCHOS_FUNC_TIME_MONITOR("LINALG::SOLVER::AMGNXN::MonolithicHierarchy::BuildSmoother");

  std::string smother_name;
  if (level<NumLevels_-1)
  {
    smother_name = params_.get<std::string>("smoother: all but coarsest level","none");
    if (smother_name=="none")
      dserror("You have to set the fine level smoother. Fix your xml file.");
  }
  else
  {
    smother_name = params_.get<std::string>("smoother: coarsest level","none");
    if (smother_name=="none")
      dserror("You have to set the coarse level smoother. Fix your xml file.");
  }

  std::string verbosity = params_.get<std::string>("verbosity","off");

  if (GetA(0)->GetMatrix(0,0)->Comm().MyPID() != 0)
    verbosity = "off";

  std::vector<int> blocks(H_->GetNumBlocks(),0);
  for(int i=0;i<H_->GetNumBlocks();i++)
    blocks[i]=i;

  AMGNXN::SmootherFactory mySmootherCreator;
  mySmootherCreator.SetOperator(GetA(level));
  mySmootherCreator.SetParamsSmoother(params_smoothers_);
  mySmootherCreator.SetHierarchies(GetHierarchies());
  mySmootherCreator.SetLevel(level);
  mySmootherCreator.SetBlocks(blocks);
  mySmootherCreator.SetSmootherName(smother_name);
  mySmootherCreator.SetVerbosity(verbosity);

  // Recover null spaces from the hierarchies and give it to the smoother creator
  std::vector<AMGNXN::NullSpaceInfo> null_space_blocks;
  for(int i=0;i<NumBlocks_;i++)
  {
    AMGNXN::NullSpaceInfo myNS(H_->GetNumPDEs(i),H_->GetNullSpaceDim(i),H_->GetNullSpaceData(i));
    null_space_blocks.push_back(myNS);
  }
  mySmootherCreator.SetNullSpaceAllBlocks(null_space_blocks);

  Teuchos::RCP<AMGNXN::GenericSmoother> Sbase = mySmootherCreator.Create();
  Teuchos::RCP<AMGNXN::BlockedSmoother>
    Sblock = Teuchos::rcp_dynamic_cast<AMGNXN::BlockedSmoother>(Sbase);
  if(Sblock == Teuchos::null)
    dserror("We expect a block smoother. Fix the xml file defining the smoother");
  return Sbase;

}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<LINALG::SOLVER::AMGNXN::Vcycle>
LINALG::SOLVER::AMGNXN::MonolithicHierarchy::BuildVCycle()
{
  int NumSweeps = 1; // Hard coded
  int FirstLevel = 0;
  Teuchos::RCP<Vcycle> V = Teuchos::rcp( new Vcycle(NumLevels_,NumSweeps,FirstLevel) );

  V->SetOperators(A_);
  V->SetProjectors(P_);
  V->SetRestrictors(R_);
  V->SetPreSmoothers(Spre_);
  V->SetPosSmoothers(Spos_);

  return V;
}




#endif // HAVE_MueLu
