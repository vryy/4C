/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of MueLu-base contact saddle-point preconditioners

\level 2

*----------------------------------------------------------------------*/

#ifdef TRILINOS_DEVELOP

// Baci
#include "solver_muelucontactsppreconditioner.H"
#include "muelu/muelu_utils.H"

#include "../drt_lib/drt_dserror.H"

// MueLu
#include <MueLu_CreateXpetraPreconditioner.hpp>
#include <MueLu_EpetraOperator.hpp>
#include <MueLu_UseShortNames.hpp>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_DefaultComm.hpp>

// EpetraExt
#include <EpetraExt_BlockMapOut.h>

// Xpetra
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_EpetraMap.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_MatrixUtils.hpp>
#include <Xpetra_UseShortNamesOrdinal.hpp>
#include <Xpetra_UseShortNamesScalar.hpp>

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuContactSpPreconditioner::MueLuContactSpPreconditioner(
    FILE* outfile, Teuchos::ParameterList& mllist)
    : PreconditionerType(outfile), mllist_(mllist)
{
  P_ = Teuchos::null;
  Pmatrix_ = Teuchos::null;
  H_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::MueLuContactSpPreconditioner::Setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  using EpetraMap = Xpetra::EpetraMapT<int, Xpetra::EpetraNode>;

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
  if (!mllist_.isParameter("PDE equations"))
    dserror("Multigrid parameter 'PDE equations' missing in solver parameter list.");
  if (!mllist_.isParameter("null space: dimension"))
    dserror("Multigrid parameter 'null space: dimension' missing  in solver parameter list.");
  const int numdf = mllist_.get<int>("PDE equations", -1);
  const int dimns = mllist_.get<int>("null space: dimension", -1);
  if (numdf == -1) dserror("Multigrid parameter 'PDE equations' wrong. It has to be > 0.");
  if (dimns == -1) dserror("Multigrid parameter 'null space: dimension' wrong. It has to be > 0.");

  // create a Teuchos::Comm from EpetraComm
  Teuchos::RCP<const Teuchos::Comm<int>> comm = Xpetra::toXpetra(A->RangeMap(0).Comm());

  /* Extract additional maps from parameter list
   *
   * These maps are provided by the STR::TimInt::PrepareContactMeshtying routine, that has access
   * to the contact manager class
   *
   * Note: Baci provides Epetra_Map objects. We will transform them to Xpetra::Map later.
   */
  Teuchos::RCP<Epetra_Map> epMasterDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epSlaveDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epActiveDofMap = Teuchos::null;
  Teuchos::RCP<Epetra_Map> epInnerDofMap = Teuchos::null;
  Teuchos::RCP<Map> xSingleNodeAggMap = Teuchos::null;
  Teuchos::RCP<Map> xNearZeroDiagMap = Teuchos::null;
  if (mllist_.isSublist("Linear System properties"))
  {
    const Teuchos::ParameterList& linSystemProps = mllist_.sublist("Linear System properties");
    // extract information provided by solver (e.g. PermutedAztecSolver)
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
  Teuchos::RCP<const Map> fullrangemap = Teuchos::null;
  Teuchos::RCP<CrsMatrix> xCrsA11 = Teuchos::null;
  Teuchos::RCP<CrsMatrix> xCrsA12 = Teuchos::null;
  Teuchos::RCP<CrsMatrix> xCrsA21 = Teuchos::null;
  Teuchos::RCP<CrsMatrix> xCrsA22 = Teuchos::null;
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

  /* Define strided maps
   *
   * We have 'numdf' Lagrange multipliers per node at the contact interface,
   * so we also add 'numdf' entries to map of the Lagrange multipliers.
   *
   * Warning: we assume a vector valued Lagrange multiplier here.
   */
  std::vector<size_t> stridingInfoPrimal;
  stridingInfoPrimal.push_back(numdf);
  Teuchos::RCP<StridedMap> stridedRangeMapPrimal = Teuchos::rcp(new StridedMap(
      xCrsA11->getRowMap(), stridingInfoPrimal, xCrsA11->getRowMap()->getIndexBase(), -1, 0));
  Teuchos::RCP<StridedMap> stridedDomainMapPrimal = Teuchos::rcp(new StridedMap(
      xCrsA11->getDomainMap(), stridingInfoPrimal, xCrsA11->getDomainMap()->getIndexBase(), -1, 0));

  std::vector<size_t> stridingInfoDual;
  stridingInfoDual.push_back(numdf);
  Teuchos::RCP<StridedMap> stridedRangeMapDual = Teuchos::rcp(new StridedMap(
      xCrsA22->getRowMap(), stridingInfoDual, xCrsA22->getRowMap()->getIndexBase(), -1, 0));
  Teuchos::RCP<StridedMap> stridedDomainMapDual = Teuchos::rcp(new StridedMap(
      xCrsA22->getDomainMap(), stridingInfoDual, xCrsA22->getDomainMap()->getIndexBase(), -1, 0));

  RCP<Matrix> xA11 = Teuchos::rcp(new CrsMatrixWrap(xCrsA11));
  RCP<Matrix> xA12 = Teuchos::rcp(new CrsMatrixWrap(xCrsA12));
  RCP<Matrix> xA21 = Teuchos::rcp(new CrsMatrixWrap(xCrsA21));
  RCP<Matrix> xA22 = Teuchos::rcp(new CrsMatrixWrap(xCrsA22));

  MatrixUtils::convertMatrixToStridedMaps(xA11, stridingInfoPrimal, stridingInfoPrimal);
  MatrixUtils::convertMatrixToStridedMaps(xA12, stridingInfoPrimal, stridingInfoDual);
  MatrixUtils::convertMatrixToStridedMaps(xA21, stridingInfoDual, stridingInfoPrimal);
  MatrixUtils::convertMatrixToStridedMaps(xA22, stridingInfoDual, stridingInfoDual);

  // build map extractor
  std::vector<Teuchos::RCP<const Map>> stridedMaps;
  stridedMaps.push_back(stridedRangeMapPrimal);
  stridedMaps.push_back(stridedRangeMapDual);

  RCP<const Map> stridedFullMap = MapUtils::concatenateMaps(stridedMaps);

  Teuchos::RCP<const ::MapExtractor> map_extractor =
      MapExtractorFactory::Build(stridedFullMap, stridedMaps);

  // build blocked Xpetra operator
  Teuchos::RCP<BlockedCrsMatrix> bOp =
      Teuchos::rcp(new BlockedCrsMatrix(map_extractor, map_extractor, 81));

  bOp->setMatrix(0, 0, xA11);
  bOp->setMatrix(0, 1, xA12);
  bOp->setMatrix(1, 0, xA21);
  bOp->setMatrix(1, 1, xA22);

  RCP<const StridedMap> testMap = Teuchos::null;
  RCP<const Matrix> xA11FromBOp = bOp->getMatrix(0, 0);
  testMap = Teuchos::rcp_dynamic_cast<const StridedMap>(xA11FromBOp->getRowMap("stridedMaps"));
  if (testMap.is_null()) dserror("Row map of A00 is no StridedMap.");

  RCP<const Matrix> xA12FromBOp = bOp->getMatrix(0, 1);
  testMap = Teuchos::rcp_dynamic_cast<const StridedMap>(xA12FromBOp->getRowMap("stridedMaps"));
  if (testMap.is_null()) dserror("Row map of A01 is no StridedMap.");

  RCP<const Matrix> xA21FromBOp = bOp->getMatrix(1, 0);
  testMap = Teuchos::rcp_dynamic_cast<const StridedMap>(xA21FromBOp->getRowMap("stridedMaps"));
  if (testMap.is_null()) dserror("Row map of A00 is no StridedMap.");

  bOp->SetFixedBlockSize(numdf);
  bOp->fillComplete();

  // Check for proper striding information
  {
    RCP<Matrix> myA11 = bOp->getMatrix(0, 0);
    Teuchos::rcp_dynamic_cast<const StridedMap>(myA11->getRowMap("stridedMaps"), true);
    Teuchos::rcp_dynamic_cast<const StridedMap>(myA11->getColMap("stridedMaps"), true);

    RCP<Matrix> myA12 = bOp->getMatrix(0, 1);
    Teuchos::rcp_dynamic_cast<const StridedMap>(myA12->getRowMap("stridedMaps"), true);
    Teuchos::rcp_dynamic_cast<const StridedMap>(myA12->getColMap("stridedMaps"), true);

    RCP<Matrix> myA21 = bOp->getMatrix(1, 0);
    Teuchos::rcp_dynamic_cast<const StridedMap>(myA21->getRowMap("stridedMaps"), true);
    Teuchos::rcp_dynamic_cast<const StridedMap>(myA21->getColMap("stridedMaps"), true);
  }

  // Re-create or re-use the preconditioner?
  if (create)
  {
    // (Re-)create the preconditioner

    // free old matrix first
    P_ = Teuchos::null;

    if (!mllist_.isParameter("MUELU_XML_FILE"))
      dserror(
          "XML-file w/ MueLu preconditioner configuration is missing in solver parameter list. "
          "Please set it as entry 'MUELU_XML_FILE'.");
    std::string xml_file = mllist_.get<std::string>("MUELU_XML_FILE");

    ParameterList mueluParams;
    Teuchos::updateParametersFromXmlFileAndBroadcast(
        xml_file, Teuchos::Ptr<ParameterList>(&mueluParams), *comm);

    // Get/compute nullspace vectors
    RCP<MultiVector> nullspace11 = Teuchos::null;
    RCP<MultiVector> nullspace22 = Teuchos::null;
    {
      // Extract pre-computed nullspace for block (0,0) from Baci's ML parameter list
      nullspace11 =
          LINALG::SOLVER::MUELU::UTILS::ExtractNullspaceFromMLList(stridedRangeMapPrimal, mllist_);

      // Compute default nullspace for block (1,1)
      {
        const int dimNS2 = numdf;
        nullspace22 = MultiVectorFactory::Build(stridedRangeMapDual, dimNS2);

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

    /* ToDo (mayr.mt) Switch to CreateXpetraPreconditioner. Pass nullspace via "user data" sublist
     * and use xml-entry "Fine level nullspace" in NullspaceFactory.
     */

    // ParameterListInterpreter mueLuFactory(xml_file, *comm);
    ParameterListInterpreter mueLuFactory(mueluParams, comm);

    RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
    H->GetLevel(0)->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp));
    H->GetLevel(0)->Set("Nullspace1", nullspace11);
    H->GetLevel(0)->Set("Nullspace2", nullspace22);
    H->GetLevel(0)->Set(
        "Primal interface DOF map", Teuchos::rcp_dynamic_cast<const Map>(xSlaveDofMap, true));

    mueLuFactory.SetupHierarchy(*H);

    // set multigrid preconditioner
    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H));

    // store multigrid hierarchy
    H_ = H;
  }
  else
  {
    /* Reuse existing multigrid hierarchy
     *
     * We use the same hierarchy, but just set the new fine level matrix.
     */

    H_->setlib(Xpetra::UseEpetra);  // not very nice, but safe.
    H_->GetLevel(0)->Set("A", Teuchos::rcp_dynamic_cast<Matrix>(bOp, true));

    P_ = Teuchos::rcp(new MueLu::EpetraOperator(H_));
  }

  return;
}

#endif  // #ifdef TRILINOS_DEVELOP
