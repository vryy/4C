/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for MueLu preconditioner

\level 1

*/

// Baci
#include "baci_linear_solver_preconditioner_muelu.hpp"

#include "baci_linalg_utils_sparse_algebra_manipulation.hpp"
#include "baci_utils_exceptions.hpp"

#include <EpetraExt_BlockMapOut.h>
#include <MueLu_AggregationExportFactory.hpp>
#include <MueLu_CoalesceDropFactory.hpp>
#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_EpetraOperator.hpp>
#include <MueLu_ML2MueLuParameterTranslator.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_VerbosityLevel.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_EpetraMap.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_StridedMap.hpp>

FOUR_C_NAMESPACE_OPEN

// define some trillinos shortcuts
using SC = Scalar;
using LO = LocalOrdinal;
using GO = GlobalOrdinal;
using NO = Node;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
CORE::LINEAR_SOLVER::MueLuPreconditioner::MueLuPreconditioner(Teuchos::ParameterList& muelulist)
    : muelulist_(muelulist)
{
  P_ = Teuchos::null;
  pmatrix_ = Teuchos::null;
  H_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::MueLuPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  SetupLinearProblem(matrix, x, b);

  // check whether A is a Epetra_CrsMatrix i.e. no block matrix
  Teuchos::RCP<Epetra_CrsMatrix> A =
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Teuchos::rcp(matrix, false));
  if (A == Teuchos::null) FOUR_C_THROW("Matrix is not a SparseMatrix");

  // store operator
  pmatrix_ = A;

  // wrap Epetra_CrsMatrix to Xpetra::Matrix for use in MueLu
  Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mueluA =
      Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>(pmatrix_));

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
        FOUR_C_THROW("Error: PDE equations or null space dimension wrong.");

      Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap =
          mueluA->getRowMap();
      Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> nullspace;
      nullspace =
          CORE::LINEAR_SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(rowMap, muelulist_);

      mueluOp->SetFixedBlockSize(numdf);

      if (!muelulist_.isParameter("MUELU_XML_FILE"))
        FOUR_C_THROW(
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

      if (muelulist_.isParameter("Coordinates"))
      {
        Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> coordinates =
            Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GO, NO>(
                muelulist_.get<Teuchos::RCP<Epetra_MultiVector>>("Coordinates")));
        if (!coordinates.is_null()) H->GetLevel(0)->Set("Coordinates", coordinates);
      }

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
    // Default: setup MueLu hierarchy from an ML-style parameter list from the dat-file
    Teuchos::RCP<Teuchos::ParameterList> muelu_params_from_ml = Teuchos::getParametersFromXmlString(
        MueLu::ML2MueLuParameterTranslator::translate(muelulist_));
    MueLu::ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node> mueLuFactory(
        *muelu_params_from_ml);
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
CORE::LINEAR_SOLVER::MueLuFluidBlockPreconditioner::MueLuFluidBlockPreconditioner(
    Teuchos::ParameterList& muelulist)
    : MueLuPreconditioner(muelulist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::MueLuFluidBlockPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;

  SetupLinearProblem(matrix, x, b);

  // adapt nullspace for splitted pure fluid problem
  int nv = 0;     // number of velocity dofs
  int np = 0;     // number of pressure dofs
  int numdf = 0;  // dofs per node

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A =
      Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
  if (A == Teuchos::null) FOUR_C_THROW("Matrix is not a BlockSparseMatrix");

  // store operator
  pmatrix_ = A;

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
        FOUR_C_THROW("Error: PDE equations or null space dimension wrong.");

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
        FOUR_C_THROW(
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
    FOUR_C_THROW("Only works with .xml file!");
  }  // else (xml file)
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
CORE::LINEAR_SOLVER::MueLuTsiBlockPreconditioner::MueLuTsiBlockPreconditioner(
    Teuchos::ParameterList& muelulist)
    : MueLuPreconditioner(muelulist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::MueLuTsiBlockPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;

  SetupLinearProblem(matrix, x, b);

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A =
      Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
  if (A == Teuchos::null) FOUR_C_THROW("matrix is not a BlockSparseMatrix");

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

  Teuchos::RCP<const Xpetra::StridedMap<LO, GO, NO>> solidmap =
      Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xA11->getRowMap()->lib(),
          xA11->getRowMap()->getGlobalNumElements(), xA11->getRowMap()->getLocalElementList(),
          xA11->getRowMap()->getIndexBase(), solidStriding, xA11->getRowMap()->getComm(), -1));
  Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> thermomap =
      Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xA22->getRowMap()->lib(),
          xA22->getRowMap()->getGlobalNumElements(), xA22->getRowMap()->getLocalElementList(),
          xA22->getRowMap()->getIndexBase(), thermoStriding, xA22->getRowMap()->getComm(), -1));

  // build map extractor
  std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>> maps;
  maps.emplace_back(solidmap);
  maps.emplace_back(thermomap);

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
        FOUR_C_THROW(
            "XML-file w/ MueLu preconditioner configuration is missing in solver parameter list. "
            "Please set it as entry 'MUELU_XML_FILE'.");
      std::string xmlFileName = mueluParameters.get<std::string>("MUELU_XML_FILE");

      int solidDimns = solidList.get<int>("null space: dimension", -1);
      if (solidDimns == -1 || solidDofs == -1)
        FOUR_C_THROW("Error: PDE equations of solid or null space dimension wrong.");

      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace11 =
          CORE::LINEAR_SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(
              solidmap, muelulist_.sublist("Inverse1"));

      int thermoDimns = thermoList.get<int>("null space: dimension", -1);
      if (thermoDimns == -1 || thermoDofs == -1)
        FOUR_C_THROW("Error: PDE equations of solid or null space dimension wrong.");

      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace22 =
          CORE::LINEAR_SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(
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
    FOUR_C_THROW(
        "The MueLu preconditioner for TSI problems only works with an appropriate .xml file");
  }
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
CORE::LINEAR_SOLVER::MueLuContactSpPreconditioner::MueLuContactSpPreconditioner(
    Teuchos::ParameterList& muelulist)
    : MueLuPreconditioner(muelulist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::MueLuContactSpPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  using EpetraMap = Xpetra::EpetraMapT<int, Xpetra::EpetraNode>;
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;

  SetupLinearProblem(matrix, x, b);

  // Check whether input matrix is an actual blocked operator
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A =
      Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
  if (A == Teuchos::null) FOUR_C_THROW("Matrix is not a BlockSparseMatrix");

  // store blocked operator
  pmatrix_ = A;

  ///////////////////////////////////////////////////////////////////////
  // Some preparations common for re-creating and re-using the preconditioner
  ///////////////////////////////////////////////////////////////////////

  Teuchos::ParameterList contactList = muelulist_.sublist("MueLu (Contact) Parameters");

  // prepare nullspace vector for MueLu (block A11 only)
  if (!contactList.isParameter("PDE equations"))
    FOUR_C_THROW("Multigrid parameter 'PDE equations' missing in solver parameter list.");
  if (!contactList.isParameter("null space: dimension"))
    FOUR_C_THROW("Multigrid parameter 'null space: dimension' missing  in solver parameter list.");
  const int numdf = contactList.get<int>("PDE equations", -1);
  const int dimns = contactList.get<int>("null space: dimension", -1);
  if (numdf == -1) FOUR_C_THROW("Multigrid parameter 'PDE equations' wrong. It has to be > 0.");
  if (dimns == -1)
    FOUR_C_THROW("Multigrid parameter 'null space: dimension' wrong. It has to be > 0.");

  // create a Teuchos::Comm from EpetraComm
  Teuchos::RCP<const Teuchos::Comm<int>> comm = Xpetra::toXpetra(A->RangeMap(0).Comm());

  // Extract additional maps from parameter list
  //
  // These maps are provided by the STR::TimInt::PrepareContactMeshtying routine, that has access
  // to the contact manager class
  //
  // Note: Baci provides Epetra_Map objects. We will transform them to Xpetra::Map later.
  //
  Teuchos::RCP<Epetra_Map> epSlaveDofMap =
      muelulist_.sublist("Belos Parameters").get<Teuchos::RCP<Epetra_Map>>("contact slaveDofMap");

  if (epSlaveDofMap.is_null())
    FOUR_C_THROW(
        "CORE::LINEAR_SOLVER::MueLuContactSpPreconditioner::MueLuContactSpPreconditioner: "
        "Interface contact map is not available!");

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
    fullrangemap = Teuchos::rcp(new EpetraMap(Teuchos::rcpFromRef(pmatrix_->FullRangeMap())));

    xCrsA11 = Teuchos::rcp(new EpetraCrsMatrix(pmatrix_->Matrix(0, 0).EpetraMatrix()));
    xCrsA12 = Teuchos::rcp(new EpetraCrsMatrix(pmatrix_->Matrix(0, 1).EpetraMatrix()));
    xCrsA21 = Teuchos::rcp(new EpetraCrsMatrix(pmatrix_->Matrix(1, 0).EpetraMatrix()));
    xCrsA22 = Teuchos::rcp(new EpetraCrsMatrix(pmatrix_->Matrix(1, 1).EpetraMatrix()));
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
  if (testMap.is_null()) FOUR_C_THROW("Row map of A00 is no StridedMap.");

  Teuchos::RCP<const Xpetra::Matrix<SC, LO, GO, NO>> xA12FromBOp = bOp->getMatrix(0, 1);
  testMap = Teuchos::rcp_dynamic_cast<const Xpetra::StridedMap<LO, GO, NO>>(
      xA12FromBOp->getRowMap("stridedMaps"));
  if (testMap.is_null()) FOUR_C_THROW("Row map of A01 is no StridedMap.");

  Teuchos::RCP<const Xpetra::Matrix<SC, LO, GO, NO>> xA21FromBOp = bOp->getMatrix(1, 0);
  testMap = Teuchos::rcp_dynamic_cast<const Xpetra::StridedMap<LO, GO, NO>>(
      xA21FromBOp->getRowMap("stridedMaps"));
  if (testMap.is_null()) FOUR_C_THROW("Row map of A00 is no StridedMap.");

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

    if (!contactList.isParameter("MUELU_XML_FILE"))
      FOUR_C_THROW(
          "XML-file w/ MueLu preconditioner configuration is missing in solver parameter list. "
          "Please set it as entry 'MUELU_XML_FILE'.");
    std::string xml_file = contactList.get<std::string>("MUELU_XML_FILE");

    MueLu::ParameterList mueluParams;
    Teuchos::updateParametersFromXmlFileAndBroadcast(
        xml_file, Teuchos::Ptr<MueLu::ParameterList>(&mueluParams), *comm);

    // Get/compute nullspace vectors
    Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace11 = Teuchos::null;
    Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace22 = Teuchos::null;
    {
      // Extract pre-computed nullspace for block (0,0) from Baci's ML parameter list
      nullspace11 = CORE::LINEAR_SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(
          stridedRangeMapPrimal, contactList);

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


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
CORE::LINEAR_SOLVER::MueLuBeamSolidBlockPreconditioner::MueLuBeamSolidBlockPreconditioner(
    Teuchos::ParameterList& muelulist)
    : MueLuPreconditioner(muelulist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::MueLuBeamSolidBlockPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  using EpetraMap = Xpetra::EpetraMapT<int, Xpetra::EpetraNode>;
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;

  SetupLinearProblem(matrix, x, b);

  // first we check if the input matrix is of normal CRS type
  Teuchos::RCP<Epetra_CrsMatrix> A =
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Teuchos::rcp(matrix, false));
  if (A == Teuchos::null) FOUR_C_THROW("Matrix is not a SparseMatrix");

  ///////////////////////////////////////////////////////////////////////
  // First get the solid and beam rowmaps and split the matrix
  // before constructing the actual preconditioner
  ///////////////////////////////////////////////////////////////////////

  Teuchos::ParameterList solidList = muelulist_.sublist("Inverse1");
  Teuchos::RCP<Epetra_Map> solidDofRowmap =
      solidList.get<Teuchos::RCP<Epetra_Map>>("null space: map", Teuchos::null);
  if (solidDofRowmap == Teuchos::null) FOUR_C_THROW("Solid row map is zero!");

  Teuchos::ParameterList beamList = muelulist_.sublist("Inverse2");
  Teuchos::RCP<Epetra_Map> beamDofRowmap =
      beamList.get<Teuchos::RCP<Epetra_Map>>("null space: map", Teuchos::null);
  if (beamDofRowmap == Teuchos::null) FOUR_C_THROW("Beam row map is zero!");

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>> Ablock =
      Teuchos::null;
  CORE::LINALG::SplitMatrix2x2(A, Ablock, solidDofRowmap, beamDofRowmap);

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

  Teuchos::RCP<const Xpetra::StridedMap<LO, GO, NO>> solidmap =
      Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xA11->getRowMap()->lib(),
          xA11->getRowMap()->getGlobalNumElements(), xA11->getRowMap()->getLocalElementList(),
          xA11->getRowMap()->getIndexBase(), solidStriding, xA11->getRowMap()->getComm(), -1));
  Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> beammap =
      Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xA22->getRowMap()->lib(),
          xA22->getRowMap()->getGlobalNumElements(), xA22->getRowMap()->getLocalElementList(),
          xA22->getRowMap()->getIndexBase(), beamStriding, xA22->getRowMap()->getComm(), -1));

  // build map extractor
  std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>> maps;
  maps.emplace_back(solidmap);
  maps.emplace_back(beammap);

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
        FOUR_C_THROW(
            "XML-file w/ MueLu preconditioner configuration is missing in solver parameter list. "
            "Please set it as entry 'MUELU_XML_FILE'.");
      std::string xmlFileName = mueluParameters.get<std::string>("MUELU_XML_FILE");

      int solidDimns = solidList.get<int>("null space: dimension", -1);
      if (solidDimns == -1 || solidDofs == -1)
        FOUR_C_THROW("Error: PDE equations of solid or null space dimension wrong.");

      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace11 =
          CORE::LINEAR_SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(solidmap, solidList);

      int beamDimns = beamList.get<int>("null space: dimension", -1);
      if (beamDimns == -1 || beamDofs == -1)
        FOUR_C_THROW("Error: PDE equations of beam or null space dimension wrong.");

      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace22 =
          CORE::LINEAR_SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(beammap, beamList);

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
    FOUR_C_THROW("Only works with .xml file");
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
CORE::LINEAR_SOLVER::MueLuFsiBlockPreconditioner::MueLuFsiBlockPreconditioner(
    Teuchos::ParameterList& muelulist)
    : MueLuPreconditioner(muelulist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void CORE::LINEAR_SOLVER::MueLuFsiBlockPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;

  SetupLinearProblem(matrix, x, b);

  if (create)
  {
    // check wheter input matrix is an actual blocked operator
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> A =
        Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
    if (A == Teuchos::null) FOUR_C_THROW("matrix is not a BlockSparseMatrix");

    // split matrix into components
    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA11 =
        Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0, 0).EpetraMatrix()));
    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA12 =
        Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0, 1).EpetraMatrix()));
    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA13 =
        Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(0, 2).EpetraMatrix()));

    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA21 =
        Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1, 0).EpetraMatrix()));
    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA22 =
        Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1, 1).EpetraMatrix()));
    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA23 =
        Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(1, 2).EpetraMatrix()));

    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA31 =
        Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(2, 0).EpetraMatrix()));
    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA32 =
        Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(2, 1).EpetraMatrix()));
    Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> xA33 =
        Teuchos::rcp(new EpetraCrsMatrix(A->Matrix(2, 2).EpetraMatrix()));

    Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> fullrangemap =
        Teuchos::rcp(new Xpetra::EpetraMapT<GO, NO>(Teuchos::rcpFromRef(A->FullRangeMap())));

    Teuchos::ParameterList& solidList = muelulist_.sublist("Inverse1");
    Teuchos::ParameterList& fluidList = muelulist_.sublist("Inverse2");
    Teuchos::ParameterList& aleList = muelulist_.sublist("Inverse3");

    // information about the fsi subproblems
    const int numSolidDofsPerNode = solidList.get<int>("PDE equations");
    const int numFluidDofsPerNode = fluidList.get<int>("PDE equations");
    const int numAleDofsPerNode = aleList.get<int>("PDE equations");

    // explicit information about the fluid
    const int numVelocityDofsPerNode =
        fluidList.sublist("NodalBlockInformation").get<int>("number of momentum dofs", 0);
    const int numPressureDofsPerNode =
        fluidList.sublist("NodalBlockInformation").get<int>("number of constraint dofs", 0);

    if (numFluidDofsPerNode != (numVelocityDofsPerNode + numPressureDofsPerNode))
      FOUR_C_THROW("Number of fluid dofs does not match the number of velocity and pressure dofs");

    // define strided maps
    std::vector<size_t> solidStriding;
    std::vector<size_t> fluidStriding;
    std::vector<size_t> aleStriding;
    solidStriding.push_back(numSolidDofsPerNode);
    fluidStriding.push_back(numVelocityDofsPerNode);
    fluidStriding.push_back(numPressureDofsPerNode);
    aleStriding.push_back(numAleDofsPerNode);

    Teuchos::RCP<const Xpetra::StridedMap<LO, GO, NO>> solidmap =
        Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xA11->getRowMap()->lib(),
            xA11->getRowMap()->getGlobalNumElements(), xA11->getRowMap()->getLocalElementList(),
            xA11->getRowMap()->getIndexBase(), solidStriding, xA11->getRowMap()->getComm(), -1));
    Teuchos::RCP<const Xpetra::StridedMap<LO, GO, NO>> fluidmap =
        Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xA22->getRowMap()->lib(),
            xA22->getRowMap()->getGlobalNumElements(), xA22->getRowMap()->getLocalElementList(),
            xA22->getRowMap()->getIndexBase(), fluidStriding, xA22->getRowMap()->getComm(), -1));
    Teuchos::RCP<const Xpetra::StridedMap<LO, GO, NO>> alemap =
        Teuchos::rcp(new Xpetra::StridedMap<LO, GO, NO>(xA33->getRowMap()->lib(),
            xA33->getRowMap()->getGlobalNumElements(), xA33->getRowMap()->getLocalElementList(),
            xA33->getRowMap()->getIndexBase(), aleStriding, xA33->getRowMap()->getComm(), -1));

    // build map extractor for fsi problem
    std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>> maps;
    maps.emplace_back(solidmap);
    maps.emplace_back(fluidmap);
    maps.emplace_back(alemap);

    Teuchos::RCP<const Xpetra::MapExtractor<SC, LO, GO, NO>> map_extractor =
        Xpetra::MapExtractorFactory<SC, LO, GO, NO>::Build(fullrangemap, maps);

    Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> velocitymap =
        Xpetra::StridedMapFactory<LO, GO, NO>::Build(fluidmap, 0);
    Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> pressuremap =
        Xpetra::StridedMapFactory<LO, GO, NO>::Build(fluidmap, 1);

    // build blocked Xpetra operator
    pmatrix_ = Teuchos::rcp(
        new Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>(map_extractor, map_extractor, 42));

    pmatrix_->setMatrix(0, 0,
        Xpetra::MatrixFactory<SC, LO, GO, NO>::BuildCopy(
            Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA11))));
    pmatrix_->setMatrix(0, 1,
        Xpetra::MatrixFactory<SC, LO, GO, NO>::BuildCopy(
            Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA12))));
    pmatrix_->setMatrix(0, 2,
        Xpetra::MatrixFactory<SC, LO, GO, NO>::BuildCopy(
            Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA13))));
    pmatrix_->setMatrix(1, 0,
        Xpetra::MatrixFactory<SC, LO, GO, NO>::BuildCopy(
            Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA21))));
    pmatrix_->setMatrix(1, 1,
        Xpetra::MatrixFactory<SC, LO, GO, NO>::BuildCopy(
            Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA22))));
    pmatrix_->setMatrix(1, 2,
        Xpetra::MatrixFactory<SC, LO, GO, NO>::BuildCopy(
            Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA23))));
    pmatrix_->setMatrix(2, 0,
        Xpetra::MatrixFactory<SC, LO, GO, NO>::BuildCopy(
            Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA31))));
    pmatrix_->setMatrix(2, 1,
        Xpetra::MatrixFactory<SC, LO, GO, NO>::BuildCopy(
            Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA32))));
    pmatrix_->setMatrix(2, 2,
        Xpetra::MatrixFactory<SC, LO, GO, NO>::BuildCopy(
            Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(xA33))));

    pmatrix_->fillComplete();

    Teuchos::ParameterList mueluParameters = muelulist_.sublist("MueLu (FSI) Parameters");

    if (!mueluParameters.isParameter("MUELU_XML_FILE"))
      FOUR_C_THROW(
          "XML-file w/ MueLu preconditioner configuration is missing in solver parameter list. "
          "Please set it as entry 'MUELU_XML_FILE'.");
    std::string xmlFileName = mueluParameters.get<std::string>("MUELU_XML_FILE");

    const int solidDimns = solidList.get<int>("null space: dimension", -1);
    if (solidDimns == -1 || numSolidDofsPerNode == -1)
      FOUR_C_THROW("Error: PDE equations of solid or null space dimension wrong.");

    Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace11 =
        CORE::LINEAR_SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(
            solidmap, muelulist_.sublist("Inverse1"));

    const int fluidDimns = fluidList.get<int>("null space: dimension", -1);
    if (fluidDimns == -1 || numFluidDofsPerNode == -1)
      FOUR_C_THROW("Error: PDE equations of fluid or null space dimension wrong.");

    Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace22 =
        CORE::LINEAR_SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(
            fluidmap, muelulist_.sublist("Inverse2"));

    const int aleDimns = aleList.get<int>("null space: dimension", -1);
    if (aleDimns == -1 || numAleDofsPerNode == -1)
      FOUR_C_THROW("Error: PDE equations of ale or null space dimension wrong.");

    Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace33 =
        CORE::LINEAR_SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(
            alemap, muelulist_.sublist("Inverse3"));

    MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(
        xmlFileName, *(pmatrix_->getRangeMap()->getComm()));

    H_ = mueLuFactory.CreateHierarchy();
    H_->GetLevel(0)->Set("A", Teuchos::rcp_dynamic_cast<Xpetra::Matrix<SC, LO, GO, NO>>(pmatrix_));
    H_->GetLevel(0)->Set("Nullspace1", nullspace11);
    H_->GetLevel(0)->Set("Nullspace2", nullspace22);
    H_->GetLevel(0)->Set("Nullspace3", nullspace33);
    H_->GetLevel(0)->Set(
        "Map1", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO, GO, NO>>(velocitymap));
    H_->GetLevel(0)->Set(
        "Map2", Teuchos::rcp_dynamic_cast<const Xpetra::Map<LO, GO, NO>>(pressuremap));

    mueLuFactory.SetupHierarchy(*H_);

    // set multigrid preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H_));
  }
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>>
CORE::LINEAR_SOLVER::MUELU::UTILS::ExtractNullspaceFromParameterlist(
    const Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>& rowMap, Teuchos::ParameterList& muelulist)
{
  // Extract info about nullspace dimension
  if (!muelulist.isParameter("null space: dimension"))
    FOUR_C_THROW("Multigrid parameter 'null space: dimension' missing  in solver parameter list.");
  const int nspDimension = muelulist.get<int>("null space: dimension");
  if (nspDimension < 1)
    FOUR_C_THROW("Multigrid parameter 'null space: dimension' wrong. It has to be > 0.");

  Teuchos::RCP<Epetra_MultiVector> nullspaceData =
      muelulist.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace", Teuchos::null);

  Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace =
      Teuchos::rcp(new Xpetra::EpetraMultiVectorT<GO, NO>(nullspaceData));

  nullspace->replaceMap(rowMap);

  return nullspace;
}

FOUR_C_NAMESPACE_CLOSE
