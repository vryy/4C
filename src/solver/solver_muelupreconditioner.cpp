/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for MueLu preconditioner

\level 1

*/

// Baci
#include "solver_muelupreconditioner.H"

#include "dserror.H"
#include "linalg_utils_sparse_algebra_manipulation.H"

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>

// MueLu
#include <MueLu_RAPFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_MLParameterListInterpreter_decl.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_AggregationExportFactory.hpp>
#include <MueLu_EpetraOperator.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_CreateXpetraPreconditioner.hpp>

// EpetraExt
#include <EpetraExt_BlockMapOut.h>

// Xpetra
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_EpetraMap.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_StridedMap.hpp>

// define some trillinos shortcuts
using SC = Scalar;
using LO = LocalOrdinal;
using GO = GlobalOrdinal;
using NO = Node;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuPreconditioner::MueLuPreconditioner(
    FILE* outfile, Teuchos::ParameterList& muelulist)
    : PreconditionerType(outfile), muelulist_(muelulist)
{
  P_ = Teuchos::null;
  Pmatrix_ = Teuchos::null;
  H_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  // check whether A is a Epetra_CrsMatrix i.e. no block matrix
  Teuchos::RCP<Epetra_CrsMatrix> A =
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Teuchos::rcp(matrix, false));
  if (A == Teuchos::null) dserror("Matrix is not a SparseMatrix");

  // store operator
  Pmatrix_ = A;

  // wrap Epetra_CrsMatrix to Xpetra::Matrix for use in MueLu
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mueluA =
      Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>(Pmatrix_));

  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mueluOp =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(mueluA));

  // remove unsupported flags, keep for now ... is used in many tests but everywhere 0
  muelulist_.remove("aggregation: threshold", false);  // no support for aggregation: threshold TODO

  if (muelulist_.get<bool>("MUELU_XML_ENFORCE"))
  {
    if (create)
    {
      // prepare nullspace vector for MueLu
      int numdf = muelulist_.get<int>("PDE equations", -1);
      int dimns = muelulist_.get<int>("null space: dimension", -1);
      if (dimns == -1 || numdf == -1)
        dserror("Error: PDE equations or null space dimension wrong.");

      Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap =
          mueluA->getRowMap();
      Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> nullspace;
      nullspace =
          LINALG::SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(rowMap, muelulist_);

      mueluOp->SetFixedBlockSize(numdf);

      if (!muelulist_.isParameter("MUELU_XML_FILE"))
        dserror(
            "XML-file w/ MueLu preconditioner configuration is missing in solver parameter list. "
            "Please set it as entry 'MUELU_XML_FILE'.");
      std::string xmlFileName = muelulist_.get<std::string>("MUELU_XML_FILE");

      MueLu::ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node> mueLuFactory(
          xmlFileName, *(mueluOp->getRowMap()->getComm()));

      Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H =
          mueLuFactory.CreateHierarchy();
      H->GetLevel(0)->Set("A", mueluOp);
      H->GetLevel(0)->Set("Nullspace", nullspace);
      H->setlib(Xpetra::UseEpetra);

      mueLuFactory.SetupHierarchy(*H);

      // set preconditioner
      P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));

      // store multigrid hierarchy
      H_ = H;

    }  // if (create)
    else
    {
      H_->setlib(Xpetra::UseEpetra);  // not very nice, but safe.
      H_->GetLevel(0)->Set("A", mueluOp);

      P_ = Teuchos::rcp(new MueLu::EpetraOperator(H_));

    }  // else (create)
  }    // if (xmlfile)
  else
  {
    // Standard case: use ML parameters from dat file
    // Setup MueLu Hierarchy
    MueLu::MLParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node> mueLuFactory(
        muelulist_);
    Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H =
        mueLuFactory.CreateHierarchy();
    H->GetLevel(0)->Set("A", mueluOp);
    H->GetLevel(0)->setlib(Xpetra::UseEpetra);
    H->setlib(Xpetra::UseEpetra);
    mueLuFactory.SetupHierarchy(*H);

    // set preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));

  }  // else (xml file)
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuFluidBlockPreconditioner::MueLuFluidBlockPreconditioner(
    FILE* outfile, Teuchos::ParameterList& muelulist)
    : MueLuPreconditioner(outfile, muelulist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuFluidBlockPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;

  SetupLinearProblem(matrix, x, b);

  // adapt nullspace for splitted pure fluid problem
  int nv = 0;     // number of velocity dofs
  int np = 0;     // number of pressure dofs
  int numdf = 0;  // dofs per node

  Teuchos::RCP<BlockSparseMatrixBase> A =
      Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
  if (A == Teuchos::null) dserror("Matrix is not a BlockSparseMatrix");

  // store operator
  Pmatrix_ = A;

  if (muelulist_.get<bool>("MUELU_XML_ENFORCE"))
  {
    if (create)
    {
      // fix null space for ML inverses
      numdf = muelulist_.sublist("NodalBlockInformation").get<int>("number of dofs per node", 0);
      nv = muelulist_.sublist("NodalBlockInformation").get<int>("number of momentum dofs", 0);
      np = muelulist_.sublist("NodalBlockInformation").get<int>("number of constraint dofs", 0);

      // build fluid null space in MueLu format
      if (numdf == 0 || nv == 0 || np == 0)
        dserror("Error: PDE equations or null space dimension wrong.");

      // define strided maps
      std::vector<size_t> stridingInfo;
      stridingInfo.push_back(nv);
      stridingInfo.push_back(np);

      // create maps
      Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> epetra_fullrangemap =
          Teuchos::rcp(new Xpetra::EpetraMapT<GO, NO>(Teuchos::rcpFromRef(A->FullRangeMap())));

      Teuchos::RCP<const Xpetra::StridedMap<LO, GO, NO>> fullrangemap =
          Xpetra::StridedMapFactory<LO, GO, NO>::Build(epetra_fullrangemap, stridingInfo, -1, 0);
      Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> strMap1 =
          Xpetra::StridedMapFactory<LO, GO, NO>::Build(fullrangemap, 0);
      Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> strMap2 =
          Xpetra::StridedMapFactory<LO, GO, NO>::Build(fullrangemap, 1);

      // split matrix into components
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA11 =
          Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0, 0).EpetraMatrix()));
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA12 =
          Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0, 1).EpetraMatrix()));
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA21 =
          Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1, 0).EpetraMatrix()));
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA22 =
          Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1, 1).EpetraMatrix()));

      // build map extractor
      std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>> xmaps;
      xmaps.push_back(strMap1);
      xmaps.push_back(strMap2);

      Teuchos::RCP<const Xpetra::MapExtractor<SC, LO, GO, NO>> map_extractor =
          Xpetra::MapExtractorFactory<SC, LO, GO, NO>::Build(fullrangemap->getMap(), xmaps);

      // build blocked Xpetra operator
      Teuchos::RCP<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>> bOp = Teuchos::rcp(
          new Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>(map_extractor, map_extractor, 10));

      bOp->setMatrix(0, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA11)));
      bOp->setMatrix(0, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA12)));
      bOp->setMatrix(1, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA21)));
      bOp->setMatrix(1, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA22)));

      bOp->fillComplete();

      // create velocity null space
      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nspVector1 =
          Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(strMap1, nv, true);
      for (int i = 0; i < numdf - 1; ++i)
      {
        Teuchos::ArrayRCP<SC> nsValues = nspVector1->getDataNonConst(i);
        int numBlocks = nsValues.size() / (numdf - 1);
        for (int j = 0; j < numBlocks; ++j)
        {
          nsValues[j * (numdf - 1) + i] = 1.0;
        }
      }

      // create pressure null space
      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nspVector2 =
          Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(strMap2, np, true);
      Teuchos::ArrayRCP<SC> nsValues2 = nspVector2->getDataNonConst(0);
      for (int j = 0; j < nsValues2.size(); ++j)
      {
        nsValues2[j] = 1.0;
      }

      if (!muelulist_.isParameter("MUELU_XML_FILE"))
        dserror(
            "XML-file w/ MueLu preconditioner configuration is missing in solver parameter list. "
            "Please set it as entry 'MUELU_XML_FILE'.");
      std::string xmlFileName = muelulist_.get<std::string>("MUELU_XML_FILE");

      MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(
          xmlFileName, *(bOp->getRangeMap()->getComm()));

      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = mueLuFactory.CreateHierarchy();

      Teuchos::RCP<MueLu::Level> Finest = H->GetLevel(0);
      H->GetLevel(0)->Set("A", Teuchos::rcp_dynamic_cast<Xpetra::Matrix<SC, LO, GO, NO>>(bOp));
      H->GetLevel(0)->Set("Nullspace1", nspVector1);
      H->GetLevel(0)->Set("Nullspace2", nspVector2);

      mueLuFactory.SetupHierarchy(*H);

      // set multigrid preconditioner
      P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));

    }  // else (create)
  }    // if (xmlfile)
  else
  {
    dserror("Only works with .xml file!");
  }  // else (xml file)
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuTsiBlockPreconditioner::MueLuTsiBlockPreconditioner(
    FILE* outfile, Teuchos::ParameterList& muelulist)
    : MueLuPreconditioner(outfile, muelulist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuTsiBlockPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;

  SetupLinearProblem(matrix, x, b);

  Teuchos::RCP<BlockSparseMatrixBase> A =
      Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
  if (A == Teuchos::null) dserror("matrix is not a BlockSparseMatrix");

  Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA11 =
      Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0, 0).EpetraMatrix()));
  Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA12 =
      Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0, 1).EpetraMatrix()));
  Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA21 =
      Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1, 0).EpetraMatrix()));
  Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA22 =
      Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1, 1).EpetraMatrix()));

  Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> fullrangemap =
      Teuchos::rcp(new Xpetra::EpetraMapT<GO, NO>(Teuchos::rcpFromRef(A->FullRangeMap())));

  Teuchos::ParameterList solidList = muelulist_.sublist("Inverse1");
  Teuchos::ParameterList thermoList = muelulist_.sublist("Inverse2");

  int solidDofs = solidList.get<int>("PDE equations");
  int thermoDofs = thermoList.get<int>("PDE equations");

  // define strided maps
  std::vector<size_t> solidStriding;
  std::vector<size_t> thermoStriding;
  solidStriding.push_back(solidDofs);
  thermoStriding.push_back(thermoDofs);

  Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> solidmap =
      Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(
          xA11->getRowMap(), solidStriding, xA11->getRowMap()->getIndexBase(), -1, 0));
  Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> thermomap =
      Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(
          xA22->getRowMap(), thermoStriding, xA22->getRowMap()->getIndexBase(), -1, 0));

  // build map extractor
  std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>> maps;
  maps.push_back(solidmap);
  maps.push_back(thermomap);

  Teuchos::RCP<const Xpetra::MapExtractor<SC, LO, GO, NO>> map_extractor =
      Xpetra::MapExtractorFactory<SC, LO, GO, NO>::Build(fullrangemap, maps);

  // build blocked Xpetra operator
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>> bOp =
      Teuchos::rcp(new Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>(map_extractor, map_extractor, 10));

  bOp->setMatrix(0, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA11)));
  bOp->setMatrix(0, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA12)));
  bOp->setMatrix(1, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA21)));
  bOp->setMatrix(1, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA22)));

  bOp->fillComplete();

  Teuchos::ParameterList mueluParameters = muelulist_.sublist("MueLu (TSI) Parameters");
  if (mueluParameters.get<bool>("MUELU_XML_ENFORCE"))
  {
    if (create)
    {
      if (!mueluParameters.isParameter("MUELU_XML_FILE"))
        dserror(
            "XML-file w/ MueLu preconditioner configuration is missing in solver parameter list. "
            "Please set it as entry 'MUELU_XML_FILE'.");
      std::string xmlFileName = mueluParameters.get<std::string>("MUELU_XML_FILE");

      int solidDimns = solidList.get<int>("null space: dimension", -1);
      if (solidDimns == -1 || solidDofs == -1)
        dserror("Error: PDE equations of solid or null space dimension wrong.");

      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace11 =
          LINALG::SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(
              solidmap, muelulist_.sublist("Inverse1"));

      int thermoDimns = thermoList.get<int>("null space: dimension", -1);
      if (thermoDimns == -1 || solidDofs == -1)
        dserror("Error: PDE equations of solid or null space dimension wrong.");

      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace22 =
          LINALG::SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(
              thermomap, muelulist_.sublist("Inverse2"));

      MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(
          xmlFileName, *(bOp->getRangeMap()->getComm()));

      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = mueLuFactory.CreateHierarchy();
      H->GetLevel(0)->Set("A", Teuchos::rcp_dynamic_cast<Xpetra::Matrix<SC, LO, GO, NO>>(bOp));
      H->GetLevel(0)->Set("Nullspace1", nullspace11);
      H->GetLevel(0)->Set("Nullspace2", nullspace22);

      mueLuFactory.SetupHierarchy(*H);

      // set multigrid preconditioner
      P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
    }
  }
  else
  {
    dserror("The MueLu preconditioner for TSI problems only works with an appropriate .xml file");
  }
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuContactSpPreconditioner::MueLuContactSpPreconditioner(
    FILE* outfile, Teuchos::ParameterList& muelulist)
    : MueLuPreconditioner(outfile, muelulist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuContactSpPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  using EpetraMap = Xpetra::EpetraMapT<int, Xpetra::EpetraNode>;
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;

  SetupLinearProblem(matrix, x, b);

  // Check whether input matrix is an actual blocked operator
  Teuchos::RCP<BlockSparseMatrixBase> A =
      Teuchos::rcp_dynamic_cast<BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
  if (A == Teuchos::null) dserror("Matrix is not a BlockSparseMatrix");

  // store blocked operator
  Pmatrix_ = A;

  ///////////////////////////////////////////////////////////////////////
  // Some preparations common for re-creating and re-using the preconditioner
  ///////////////////////////////////////////////////////////////////////

  // prepare nullspace vector for MueLu (block A11 only)
  if (!muelulist_.isParameter("PDE equations"))
    dserror("Multigrid parameter 'PDE equations' missing in solver parameter list.");
  if (!muelulist_.isParameter("null space: dimension"))
    dserror("Multigrid parameter 'null space: dimension' missing  in solver parameter list.");
  const int numdf = muelulist_.get<int>("PDE equations", -1);
  const int dimns = muelulist_.get<int>("null space: dimension", -1);
  if (numdf == -1) dserror("Multigrid parameter 'PDE equations' wrong. It has to be > 0.");
  if (dimns == -1) dserror("Multigrid parameter 'null space: dimension' wrong. It has to be > 0.");

  // create a Teuchos::Comm from EpetraComm
  Teuchos::RCP<const Teuchos::Comm<int>> comm = Xpetra::toXpetra(A->RangeMap(0).Comm());

  // Extract additional maps from parameter list
  //
  // These maps are provided by the STR::TimInt::PrepareContactMeshtying routine, that has access
  // to the contact manager class
  //
  // Note: Baci provides Epetra_Map objects. We will transform them to Xpetra::Map later.
  //
  Teuchos::RCP<Epetra_Map> epMasterDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epSlaveDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epActiveDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epInnerDofMap = Teuchos::null;
  Teuchos::RCP<Xpetra::Map<LO, GO, NO>> xSingleNodeAggMap = Teuchos::null;
  Teuchos::RCP<Xpetra::Map<LO, GO, NO>> xNearZeroDiagMap = Teuchos::null;
  if (muelulist_.isSublist("Linear System properties"))
  {
    const Teuchos::ParameterList& linSystemProps = muelulist_.sublist("Linear System properties");
    epMasterDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact masterDofMap");
    epSlaveDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact slaveDofMap");
    epActiveDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact activeDofMap");
    epInnerDofMap = linSystemProps.get<Teuchos::RCP<Epetra_Map>>("contact innerDofMap");
  }

  ///////////////////////////////////////////////////////////////////////
  // Transform from epetra to xpetra
  ///////////////////////////////////////////////////////////////////////

  // Transform maps
  Teuchos::RCP<EpetraMap> xSlaveDofMap = Teuchos::rcp(new EpetraMap(epSlaveDofMap));

  // Get maps and matrix blocks for blocked operator
  Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> fullrangemap = Teuchos::null;
  Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xCrsA11 = Teuchos::null;
  Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xCrsA12 = Teuchos::null;
  Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xCrsA21 = Teuchos::null;
  Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xCrsA22 = Teuchos::null;
  if (create)
  {
    // (Re-)create the preconditioner, so extract everything from the input matrix 'A'.

    // Prepare maps for blocked operator
    fullrangemap = Teuchos::rcp(new EpetraMap(Teuchos::rcpFromRef(A->FullRangeMap())));

    // Transform matrix blocks
    xCrsA11 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0, 0).EpetraMatrix()));
    xCrsA12 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0, 1).EpetraMatrix()));
    xCrsA21 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1, 0).EpetraMatrix()));
    xCrsA22 = Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1, 1).EpetraMatrix()));
  }
  else
  {
    // Re-use the preconditioner, so extract everything from the existing preconditioner 'Pmatrix_'.

    // create maps
    fullrangemap = Teuchos::rcp(new EpetraMap(Teuchos::rcpFromRef(Pmatrix_->FullRangeMap())));

    xCrsA11 = Teuchos::rcp(new EpetraCrsMatrix(Pmatrix_->Matrix(0, 0).EpetraMatrix()));
    xCrsA12 = Teuchos::rcp(new EpetraCrsMatrix(Pmatrix_->Matrix(0, 1).EpetraMatrix()));
    xCrsA21 = Teuchos::rcp(new EpetraCrsMatrix(Pmatrix_->Matrix(1, 0).EpetraMatrix()));
    xCrsA22 = Teuchos::rcp(new EpetraCrsMatrix(Pmatrix_->Matrix(1, 1).EpetraMatrix()));
  }

  // Define strided maps
  //
  // We have 'numdf' Lagrange multipliers per node at the contact interface,
  // so we also add 'numdf' entries to map of the Lagrange multipliers.
  //
  // Warning: we assume a vector valued Lagrange multiplier here.
  //

  std::vector<size_t> stridingInfoPrimal;
  stridingInfoPrimal.push_back(numdf);
  Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> stridedRangeMapPrimal =
      Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(
          xCrsA11->getRowMap(), stridingInfoPrimal, xCrsA11->getRowMap()->getIndexBase(), -1, 0));
  Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> stridedDomainMapPrimal =
      Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xCrsA11->getDomainMap(), stridingInfoPrimal,
          xCrsA11->getDomainMap()->getIndexBase(), -1, 0));

  std::vector<size_t> stridingInfoDual;
  stridingInfoDual.push_back(numdf);
  Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> stridedRangeMapDual =
      Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(
          xCrsA22->getRowMap(), stridingInfoDual, xCrsA22->getRowMap()->getIndexBase(), -1, 0));
  Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> stridedDomainMapDual =
      Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xCrsA22->getDomainMap(), stridingInfoDual,
          xCrsA22->getDomainMap()->getIndexBase(), -1, 0));

  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> xA11 =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xCrsA11));
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> xA12 =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xCrsA12));
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> xA21 =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xCrsA21));
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> xA22 =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xCrsA22));

  Xpetra::MatrixUtils<SC, LO, GO, NO>::convertMatrixToStridedMaps(
      xA11, stridingInfoPrimal, stridingInfoPrimal);
  Xpetra::MatrixUtils<SC, LO, GO, NO>::convertMatrixToStridedMaps(
      xA12, stridingInfoPrimal, stridingInfoDual);
  Xpetra::MatrixUtils<SC, LO, GO, NO>::convertMatrixToStridedMaps(
      xA21, stridingInfoDual, stridingInfoPrimal);
  Xpetra::MatrixUtils<SC, LO, GO, NO>::convertMatrixToStridedMaps(
      xA22, stridingInfoDual, stridingInfoDual);

  // build map extractor
  std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>> stridedMaps;
  stridedMaps.push_back(stridedRangeMapPrimal);
  stridedMaps.push_back(stridedRangeMapDual);

  Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> stridedFullMap =
      Xpetra::MapUtils<LO, GO, NO>::concatenateMaps(stridedMaps);

  Teuchos::RCP<const Xpetra::MapExtractor<SC, LO, GO, NO>> map_extractor =
      Xpetra::MapExtractorFactory<SC, LO, GO, NO>::Build(stridedFullMap, stridedMaps);

  // build blocked Xpetra operator
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>> bOp =
      Teuchos::rcp(new Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>(map_extractor, map_extractor, 81));

  bOp->setMatrix(0, 0, xA11);
  bOp->setMatrix(0, 1, xA12);
  bOp->setMatrix(1, 0, xA21);
  bOp->setMatrix(1, 1, xA22);

  Teuchos::RCP<const Xpetra::StridedMap<LO, GO, NO>> testMap = Teuchos::null;
  Teuchos::RCP<const Xpetra::Matrix<SC, LO, GO, NO>> xA11FromBOp = bOp->getMatrix(0, 0);
  testMap = Teuchos::rcp_dynamic_cast<const Xpetra::StridedMap<LO, GO, NO>>(
      xA11FromBOp->getRowMap("stridedMaps"));
  if (testMap.is_null()) dserror("Row map of A00 is no StridedMap.");

  Teuchos::RCP<const Xpetra::Matrix<SC, LO, GO, NO>> xA12FromBOp = bOp->getMatrix(0, 1);
  testMap = Teuchos::rcp_dynamic_cast<const Xpetra::StridedMap<LO, GO, NO>>(
      xA12FromBOp->getRowMap("stridedMaps"));
  if (testMap.is_null()) dserror("Row map of A01 is no StridedMap.");

  Teuchos::RCP<const Xpetra::Matrix<SC, LO, GO, NO>> xA21FromBOp = bOp->getMatrix(1, 0);
  testMap = Teuchos::rcp_dynamic_cast<const Xpetra::StridedMap<LO, GO, NO>>(
      xA21FromBOp->getRowMap("stridedMaps"));
  if (testMap.is_null()) dserror("Row map of A00 is no StridedMap.");

  bOp->SetFixedBlockSize(numdf);
  bOp->fillComplete();

  // Check for proper striding information
  {
    Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> myA11 = bOp->getMatrix(0, 0);
    Teuchos::rcp_dynamic_cast<const Xpetra::StridedMap<LO, GO, NO>>(
        myA11->getRowMap("stridedMaps"), true);
    Teuchos::rcp_dynamic_cast<const Xpetra::StridedMap<LO, GO, NO>>(
        myA11->getColMap("stridedMaps"), true);

    Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> myA12 = bOp->getMatrix(0, 1);
    Teuchos::rcp_dynamic_cast<const Xpetra::StridedMap<LO, GO, NO>>(
        myA12->getRowMap("stridedMaps"), true);
    Teuchos::rcp_dynamic_cast<const Xpetra::StridedMap<LO, GO, NO>>(
        myA12->getColMap("stridedMaps"), true);

    Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> myA21 = bOp->getMatrix(1, 0);
    Teuchos::rcp_dynamic_cast<const Xpetra::StridedMap<LO, GO, NO>>(
        myA21->getRowMap("stridedMaps"), true);
    Teuchos::rcp_dynamic_cast<const Xpetra::StridedMap<LO, GO, NO>>(
        myA21->getColMap("stridedMaps"), true);
  }

  // Re-create or re-use the preconditioner?
  if (create)
  {
    // (Re-)create the preconditioner

    // free old matrix first
    P_ = Teuchos::null;

    if (!muelulist_.isParameter("MUELU_XML_FILE"))
      dserror(
          "XML-file w/ MueLu preconditioner configuration is missing in solver parameter list. "
          "Please set it as entry 'MUELU_XML_FILE'.");
    std::string xml_file = muelulist_.get<std::string>("MUELU_XML_FILE");

    MueLu::ParameterList mueluParams;
    Teuchos::updateParametersFromXmlFileAndBroadcast(
        xml_file, Teuchos::Ptr<MueLu::ParameterList>(&mueluParams), *comm);

    // Get/compute nullspace vectors
    Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace11 = Teuchos::null;
    Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace22 = Teuchos::null;
    {
      // Extract pre-computed nullspace for block (0,0) from Baci's ML parameter list
      nullspace11 = LINALG::SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(
          stridedRangeMapPrimal, muelulist_);

      // Compute default nullspace for block (1,1)
      {
        const int dimNS2 = numdf;
        nullspace22 =
            Xpetra::MultiVectorFactory<SC, LO, GO, NO>::Build(stridedRangeMapDual, dimNS2);

        for (int i = 0; i < dimNS2; ++i)
        {
          Teuchos::ArrayRCP<Scalar> nsValues22 = nullspace22->getDataNonConst(i);
          int numBlocks = nsValues22.size() / dimNS2;
          for (int j = 0; j < numBlocks; ++j)
          {
            nsValues22[j * dimNS2 + i] = 1.0;
          }
        }
      }
    }

    // ToDo (mayr.mt) Switch to CreateXpetraPreconditioner. Pass nullspace via "user data" sublist
    // and use xml-entry "Fine level nullspace" in NullspaceFactory.
    //

    // ParameterListInterpreter mueLuFactory(xml_file, *comm);
    MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(mueluParams, comm);

    Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = mueLuFactory.CreateHierarchy();
    H->GetLevel(0)->Set("A", Teuchos::rcp_dynamic_cast<Xpetra::Matrix<SC, LO, GO, NO>>(bOp));
    H->GetLevel(0)->Set("Nullspace1", nullspace11);
    H->GetLevel(0)->Set("Nullspace2", nullspace22);
    H->GetLevel(0)->Set("Primal interface DOF map",
        Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO, GO, NO>>(xSlaveDofMap, true));

    mueLuFactory.SetupHierarchy(*H);

    // set multigrid preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));

    // store multigrid hierarchy
    H_ = H;
  }
  else
  {
    // Reuse existing multigrid hierarchy
    //
    // We use the same hierarchy, but just set the new fine level matrix.
    //

    H_->setlib(Xpetra::UseEpetra);  // not very nice, but safe.
    H_->GetLevel(0)->Set("A", Teuchos::rcp_dynamic_cast<Xpetra::Matrix<SC, LO, GO, NO>>(bOp, true));

    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H_));
  }

  return;
}

#ifdef TRILINOS_DEVELOP

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuBeamSolidBlockPreconditioner::MueLuBeamSolidBlockPreconditioner(
    FILE* outfile, Teuchos::ParameterList& muelulist)
    : MueLuPreconditioner(outfile, muelulist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuBeamSolidBlockPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  using EpetraMap = Xpetra::EpetraMapT<int, Xpetra::EpetraNode>;
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;

  SetupLinearProblem(matrix, x, b);

  // first we check if the input matrix is of normal CRS type
  Teuchos::RCP<Epetra_CrsMatrix> A =
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Teuchos::rcp(matrix, false));
  if (A == Teuchos::null) dserror("Matrix is not a SparseMatrix");

  ///////////////////////////////////////////////////////////////////////
  // First get the solid and beam rowmaps and split the matrix
  // before constructing the actual preconditioner
  ///////////////////////////////////////////////////////////////////////

  Teuchos::ParameterList solidList = muelulist_.sublist("Inverse1");
  Teuchos::RCP<Epetra_Map> solidDofRowmap =
      solidList.get<Teuchos::RCP<Epetra_Map>>("null space: map", Teuchos::null);
  if (solidDofRowmap == Teuchos::null) dserror("Solid row map is zero!");

  Teuchos::ParameterList beamList = muelulist_.sublist("Inverse2");
  Teuchos::RCP<Epetra_Map> beamDofRowmap =
      beamList.get<Teuchos::RCP<Epetra_Map>>("null space: map", Teuchos::null);
  if (beamDofRowmap == Teuchos::null) dserror("Beam row map is zero!");

  Teuchos::RCP<BlockSparseMatrix<DefaultBlockMatrixStrategy>> Ablock = Teuchos::null;
  LINALG::SplitMatrix2x2(A, Ablock, solidDofRowmap, beamDofRowmap);

  ///////////////////////////////////////////////////////////////////////
  // Here the actual construction of the preconditioner starts
  ///////////////////////////////////////////////////////////////////////

  Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA11 =
      Teuchos::rcp(new EpetraCrsMatrix(Ablock->Matrix(0, 0).EpetraMatrix()));
  Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA12 =
      Teuchos::rcp(new EpetraCrsMatrix(Ablock->Matrix(0, 1).EpetraMatrix()));
  Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA21 =
      Teuchos::rcp(new EpetraCrsMatrix(Ablock->Matrix(1, 0).EpetraMatrix()));
  Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA22 =
      Teuchos::rcp(new EpetraCrsMatrix(Ablock->Matrix(1, 1).EpetraMatrix()));

  Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> fullrangemap =
      Teuchos::rcp(new EpetraMap(Teuchos::rcpFromRef(A->RangeMap())));

  // define strided maps
  int solidDofs = solidList.get<int>("PDE equations");
  int beamDofs = beamList.get<int>("PDE equations");

  std::vector<size_t> solidStriding;
  std::vector<size_t> beamStriding;
  solidStriding.push_back(solidDofs);
  beamStriding.push_back(beamDofs);

  Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> solidmap =
      Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(
          xA11->getRowMap(), solidStriding, xA11->getRowMap()->getIndexBase(), -1, 0));
  Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> beammap =
      Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(
          xA22->getRowMap(), beamStriding, xA22->getRowMap()->getIndexBase(), -1, 0));

  // build map extractor
  std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>> maps;
  maps.push_back(solidmap);
  maps.push_back(beammap);

  Teuchos::RCP<const Xpetra::MapExtractor<SC, LO, GO, NO>> map_extractor =
      Xpetra::MapExtractorFactory<SC, LO, GO, NO>::Build(fullrangemap, maps);

  // build blocked Xpetra operator
  Teuchos::RCP<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>> bOp =
      Teuchos::rcp(new Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>(map_extractor, map_extractor, 81));

  bOp->setMatrix(0, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA11)));
  bOp->setMatrix(0, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA12)));
  bOp->setMatrix(1, 0, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA21)));
  bOp->setMatrix(1, 1, Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA22)));

  bOp->fillComplete();

  Teuchos::ParameterList mueluParameters = muelulist_.sublist("MueLu (BeamSolid) Parameters");
  if (mueluParameters.get<bool>("MUELU_XML_ENFORCE"))
  {
    if (create)
    {
      if (!mueluParameters.isParameter("MUELU_XML_FILE"))
        dserror(
            "XML-file w/ MueLu preconditioner configuration is missing in solver parameter list. "
            "Please set it as entry 'MUELU_XML_FILE'.");
      std::string xmlFileName = mueluParameters.get<std::string>("MUELU_XML_FILE");

      int solidDimns = solidList.get<int>("null space: dimension", -1);
      if (solidDimns == -1 || solidDofs == -1)
        dserror("Error: PDE equations of solid or null space dimension wrong.");

      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace11 =
          LINALG::SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(solidmap, solidList);

      int beamDimns = beamList.get<int>("null space: dimension", -1);
      if (beamDimns == -1 || beamDofs == -1)
        dserror("Error: PDE equations of beam or null space dimension wrong.");

      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace22 =
          LINALG::SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(beammap, beamList);

      MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(
          xmlFileName, *(bOp->getRangeMap()->getComm()));

      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = mueLuFactory.CreateHierarchy();
      H->GetLevel(0)->Set("A", Teuchos::rcp_dynamic_cast<Xpetra::Matrix<SC, LO, GO, NO>>(bOp));
      H->GetLevel(0)->Set("Nullspace1", nullspace11);
      H->GetLevel(0)->Set("Nullspace2", nullspace22);
      H->setlib(Xpetra::UseEpetra);

      mueLuFactory.SetupHierarchy(*H);

      // set preconditioner
      P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));
    }
  }
  else
  {
    dserror("Only works with .xml file");
  }
}

#endif

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>>
LINALG::SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(
    const Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>& rowMap, Teuchos::ParameterList& muelulist)
{
  // Extract info about nullspace dimension
  if (!muelulist.isParameter("null space: dimension"))
    dserror("Multigrid parameter 'null space: dimension' missing  in solver parameter list.");
  const int nspDimension = muelulist.get<int>("null space: dimension");
  if (nspDimension < 1)
    dserror("Multigrid parameter 'null space: dimension' wrong. It has to be > 0.");

  Teuchos::RCP<Epetra_MultiVector> nullspaceData =
      muelulist.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace", Teuchos::null);

  Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace =
      Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GO, NO>(nullspaceData));

  nullspace->replaceMap(rowMap);

  return nullspace;
}
