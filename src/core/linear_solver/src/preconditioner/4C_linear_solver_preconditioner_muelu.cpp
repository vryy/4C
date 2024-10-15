/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for MueLu preconditioner

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_linear_solver_preconditioner_muelu.hpp"

#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method_parameters.hpp"
#include "4C_utils_exceptions.hpp"

#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_EpetraOperator.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_UseDefaultTypes.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_EpetraMap.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapExtractor.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_StridedMap.hpp>

FOUR_C_NAMESPACE_OPEN

using SC = Scalar;
using LO = LocalOrdinal;
using GO = GlobalOrdinal;
using NO = Node;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Core::LinearSolver::MueLuPreconditioner::MueLuPreconditioner(Teuchos::ParameterList& muelulist)
    : muelulist_(muelulist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void Core::LinearSolver::MueLuPreconditioner::setup(bool create, Epetra_Operator* matrix,
    Core::LinAlg::MultiVector<double>* x, Core::LinAlg::MultiVector<double>* b)
{
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<GO, NO>;
  using EpetraMap = Xpetra::EpetraMapT<GO, NO>;
  using EpetraMultiVector = Xpetra::EpetraMultiVectorT<GO, NO>;

  if (create)
  {
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> A =
        Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(
            Teuchos::rcpFromRef(*matrix));

    if (A.is_null())
    {
      Teuchos::RCP<Epetra_CrsMatrix> crsA =
          Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Teuchos::rcpFromRef(*matrix));

      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> mueluA =
          Teuchos::make_rcp<EpetraCrsMatrix>(crsA);
      pmatrix_ = Xpetra::MatrixFactory<SC, LO, GO, NO>::BuildCopy(
          Teuchos::make_rcp<Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(mueluA));

      Teuchos::ParameterList& inverseList = muelulist_.sublist("MueLu Parameters");

      std::string xmlFileName = inverseList.get<std::string>("MUELU_XML_FILE");
      if (xmlFileName == "none") FOUR_C_THROW("MUELU_XML_FILE parameter not set!");

      Teuchos::RCP<Teuchos::ParameterList> muelu_params =
          Teuchos::make_rcp<Teuchos::ParameterList>();
      auto comm = pmatrix_->getRowMap()->getComm();
      Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, muelu_params.ptr(), *comm);

      const int number_of_equations = inverseList.get<int>("PDE equations");
      pmatrix_->SetFixedBlockSize(number_of_equations);

      Teuchos::RCP<const Xpetra::Map<LO, GO, NO>> row_map = mueluA->getRowMap();
      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace =
          Core::LinearSolver::Parameters::extract_nullspace_from_parameterlist(
              row_map, inverseList);

      Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> coordinates =
          Teuchos::make_rcp<EpetraMultiVector>(
              inverseList.get<Teuchos::RCP<Core::LinAlg::MultiVector<double>>>("Coordinates")
                  ->get_ptr_of_Epetra_MultiVector());

      muelu_params->set("number of equations", number_of_equations);
      Teuchos::ParameterList& user_param_list = muelu_params->sublist("user data");
      user_param_list.set("Nullspace", nullspace);
      user_param_list.set("Coordinates", coordinates);

      H_ = MueLu::CreateXpetraPreconditioner(pmatrix_, *muelu_params);
      P_ = Teuchos::make_rcp<MueLu::EpetraOperator>(H_);
    }
    else
    {
      std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>> maps;

      for (int block = 0; block < A->rows(); block++)
      {
        EpetraCrsMatrix crsA(A->matrix(block, block).epetra_matrix());

        std::string inverse = "Inverse" + std::to_string(block + 1);
        Teuchos::ParameterList& inverseList =
            muelulist_.sublist(inverse).sublist("MueLu Parameters");
        const int number_of_equations = inverseList.get<int>("PDE equations");

        std::vector<size_t> striding;
        striding.emplace_back(number_of_equations);

        Teuchos::RCP<const Xpetra::StridedMap<LO, GO, NO>> map =
            Teuchos::make_rcp<Xpetra::StridedMap<LO, GO, NO>>(crsA.getRowMap()->lib(),
                crsA.getRowMap()->getGlobalNumElements(), crsA.getRowMap()->getLocalElementList(),
                crsA.getRowMap()->getIndexBase(), striding, crsA.getRowMap()->getComm(), -1);

        maps.emplace_back(map);
      }

      Teuchos::RCP<const EpetraMap> fullrangemap =
          Teuchos::make_rcp<EpetraMap>(Teuchos::rcpFromRef(A->full_range_map()));
      Teuchos::RCP<const Xpetra::MapExtractor<SC, LO, GO, NO>> map_extractor =
          Xpetra::MapExtractorFactory<SC, LO, GO, NO>::Build(fullrangemap, maps);

      auto bOp = Teuchos::make_rcp<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>>(
          map_extractor, map_extractor, 42);

      for (int row = 0; row < A->rows(); row++)
      {
        for (int col = 0; col < A->cols(); col++)
        {
          Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> crsA =
              Teuchos::make_rcp<EpetraCrsMatrix>(A->matrix(row, col).epetra_matrix());
          Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> mat =
              Xpetra::MatrixFactory<SC, LO, GO, NO>::BuildCopy(
                  Teuchos::make_rcp<Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(crsA));
          bOp->setMatrix(row, col, mat);
        }
      }
      bOp->fillComplete();
      pmatrix_ = bOp;

      if (!muelulist_.sublist("MueLu Parameters").isParameter("MUELU_XML_FILE"))
        FOUR_C_THROW("MUELU_XML_FILE parameter not set!");

      std::string xmlFileName =
          muelulist_.sublist("MueLu Parameters").get<std::string>("MUELU_XML_FILE");
      Teuchos::RCP<Teuchos::ParameterList> mueluParams =
          Teuchos::make_rcp<Teuchos::ParameterList>();
      auto comm = pmatrix_->getRowMap()->getComm();
      Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, mueluParams.ptr(), *comm);

      MueLu::ParameterListInterpreter<SC, LO, GO, NO> mueLuFactory(xmlFileName, *comm);
      Teuchos::RCP<MueLu::Hierarchy<SC, LO, GO, NO>> H = mueLuFactory.CreateHierarchy();
      H->GetLevel(0)->Set("A", Teuchos::rcp_dynamic_cast<Xpetra::Matrix<SC, LO, GO, NO>>(pmatrix_));

      for (int block = 0; block < A->rows(); block++)
      {
        std::string inverse = "Inverse" + std::to_string(block + 1);
        Teuchos::ParameterList& inverseList =
            muelulist_.sublist(inverse).sublist("MueLu Parameters");

        Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace =
            Teuchos::make_rcp<EpetraMultiVector>(
                inverseList.get<Teuchos::RCP<Core::LinAlg::MultiVector<double>>>("nullspace")
                    ->get_ptr_of_Epetra_MultiVector());

        H->GetLevel(0)->Set("Nullspace" + std::to_string(block + 1), nullspace);
      }

      mueLuFactory.SetupHierarchy(*H);

      P_ = Teuchos::make_rcp<MueLu::EpetraOperator>(H);
    }
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Core::LinearSolver::MueLuContactSpPreconditioner::MueLuContactSpPreconditioner(
    Teuchos::ParameterList& muelulist)
    : MueLuPreconditioner(muelulist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void Core::LinearSolver::MueLuContactSpPreconditioner::setup(bool create, Epetra_Operator* matrix,
    Core::LinAlg::MultiVector<double>* x, Core::LinAlg::MultiVector<double>* b)
{
  using EpetraMap = Xpetra::EpetraMapT<int, Xpetra::EpetraNode>;
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>;

  // Check whether input matrix is an actual blocked operator
  Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> A =
      Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(Teuchos::rcpFromRef(*matrix));
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
  Teuchos::RCP<const Teuchos::Comm<int>> comm = Xpetra::toXpetra(A->range_map(0).Comm());

  // Extract additional maps from parameter list
  //
  // These maps are provided by the Solid::TimInt::prepare_contact_meshtying routine, that has
  // access to the contact manager class
  //
  // Note: 4C provides Epetra_Map objects. We will transform them to Xpetra::Map later.
  //
  Teuchos::RCP<Epetra_Map> epSlaveDofMap =
      muelulist_.sublist("Belos Parameters").get<Teuchos::RCP<Epetra_Map>>("contact slaveDofMap");

  if (epSlaveDofMap.is_null())
    FOUR_C_THROW(
        "Core::LinearSolver::MueLuContactSpPreconditioner::MueLuContactSpPreconditioner: "
        "Interface contact map is not available!");

  ///////////////////////////////////////////////////////////////////////
  // Transform from epetra to xpetra
  ///////////////////////////////////////////////////////////////////////

  // Transform maps
  Teuchos::RCP<EpetraMap> xSlaveDofMap = Teuchos::make_rcp<EpetraMap>(epSlaveDofMap);

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
    fullrangemap = Teuchos::make_rcp<EpetraMap>(Teuchos::rcpFromRef(A->full_range_map()));

    // Transform matrix blocks
    xCrsA11 = Teuchos::make_rcp<EpetraCrsMatrix>(A->matrix(0, 0).epetra_matrix());
    xCrsA12 = Teuchos::make_rcp<EpetraCrsMatrix>(A->matrix(0, 1).epetra_matrix());
    xCrsA21 = Teuchos::make_rcp<EpetraCrsMatrix>(A->matrix(1, 0).epetra_matrix());
    xCrsA22 = Teuchos::make_rcp<EpetraCrsMatrix>(A->matrix(1, 1).epetra_matrix());
  }
  else
  {
    // Re-use the preconditioner, so extract everything from the existing preconditioner 'Pmatrix_'.

    // create maps
    fullrangemap = Teuchos::make_rcp<EpetraMap>(Teuchos::rcpFromRef(pmatrix_->full_range_map()));

    xCrsA11 = Teuchos::make_rcp<EpetraCrsMatrix>(pmatrix_->matrix(0, 0).epetra_matrix());
    xCrsA12 = Teuchos::make_rcp<EpetraCrsMatrix>(pmatrix_->matrix(0, 1).epetra_matrix());
    xCrsA21 = Teuchos::make_rcp<EpetraCrsMatrix>(pmatrix_->matrix(1, 0).epetra_matrix());
    xCrsA22 = Teuchos::make_rcp<EpetraCrsMatrix>(pmatrix_->matrix(1, 1).epetra_matrix());
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
      Teuchos::make_rcp<Xpetra::StridedMap<LO, GO, NO>>(
          xCrsA11->getRowMap(), stridingInfoPrimal, xCrsA11->getRowMap()->getIndexBase(), -1, 0);
  Xpetra::StridedMap<LO, GO, NO> stridedDomainMapPrimal(
      xCrsA11->getDomainMap(), stridingInfoPrimal, xCrsA11->getDomainMap()->getIndexBase(), -1, 0);

  std::vector<size_t> stridingInfoDual;
  stridingInfoDual.push_back(numdf);
  Teuchos::RCP<Xpetra::StridedMap<LO, GO, NO>> stridedRangeMapDual =
      Teuchos::make_rcp<Xpetra::StridedMap<LO, GO, NO>>(
          xCrsA22->getRowMap(), stridingInfoDual, xCrsA22->getRowMap()->getIndexBase(), -1, 0);
  Xpetra::StridedMap<LO, GO, NO> stridedDomainMapDual(
      xCrsA22->getDomainMap(), stridingInfoDual, xCrsA22->getDomainMap()->getIndexBase(), -1, 0);

  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> xA11 =
      Teuchos::make_rcp<Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(xCrsA11);
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> xA12 =
      Teuchos::make_rcp<Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(xCrsA12);
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> xA21 =
      Teuchos::make_rcp<Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(xCrsA21);
  Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> xA22 =
      Teuchos::make_rcp<Xpetra::CrsMatrixWrap<SC, LO, GO, NO>>(xCrsA22);

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
      Teuchos::make_rcp<Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>>(map_extractor, map_extractor, 81);

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
      // Extract pre-computed nullspace for block (0,0) from 4C's MueLu parameter list
      nullspace11 = Core::LinearSolver::Parameters::extract_nullspace_from_parameterlist(
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
    P_ = Teuchos::make_rcp<MueLu::EpetraOperator>(H);

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

    P_ = Teuchos::make_rcp<MueLu::EpetraOperator>(H_);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
