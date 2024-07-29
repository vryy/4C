/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for Teko block preconditioner

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_linear_solver_preconditioner_teko.hpp"

#include "4C_comm_utils.hpp"
#include "4C_linear_solver_method_parameters.hpp"

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Stratimikos_MueLuHelpers.hpp>
#include <Teko_EpetraInverseOpWrapper.hpp>
#include <Teko_InverseLibrary.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_StridedMap.hpp>
#include <Xpetra_ThyraUtils.hpp>

FOUR_C_NAMESPACE_OPEN

using SC = Scalar;
using LO = LocalOrdinal;
using GO = GlobalOrdinal;
using NO = Node;

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
Core::LinearSolver::TekoPreconditioner::TekoPreconditioner(Teuchos::ParameterList& tekolist)
    : tekolist_(tekolist)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void Core::LinearSolver::TekoPreconditioner::setup(
    bool create, Epetra_Operator* matrix, Epetra_MultiVector* x, Epetra_MultiVector* b)
{
  using EpetraCrsMatrix = Xpetra::EpetraCrsMatrixT<GO, NO>;
  using EpetraMap = Xpetra::EpetraMapT<GO, NO>;
  using EpetraMultiVector = Xpetra::EpetraMultiVectorT<GO, NO>;

  if (create)
  {
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> A =
        Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(Teuchos::rcp(matrix, false));
    if (A == Teuchos::null) FOUR_C_THROW("Matrix is not a BlockSparseMatrix!");

    std::vector<Teuchos::RCP<const Xpetra::Map<LO, GO, NO>>> maps;

    for (int block = 0; block < A->rows(); block++)
    {
      Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> crsA =
          Teuchos::rcp(new EpetraCrsMatrix(A->matrix(block, block).epetra_matrix()));

      Teuchos::RCP<Xpetra::Map<LO, GO, NO>> map =
          Xpetra::MapFactory<LO, GO, NO>::Build(crsA->getRowMap()->lib(),
              crsA->getRowMap()->getGlobalNumElements(), crsA->getRowMap()->getLocalElementList(),
              crsA->getRowMap()->getIndexBase(), crsA->getRowMap()->getComm());

      maps.emplace_back(map);
    }

    Teuchos::RCP<const EpetraMap> fullrangemap =
        Teuchos::rcp(new EpetraMap(Teuchos::rcpFromRef(A->full_range_map())));
    Teuchos::RCP<const Xpetra::MapExtractor<SC, LO, GO, NO>> map_extractor =
        Xpetra::MapExtractorFactory<SC, LO, GO, NO>::Build(fullrangemap, maps);

    pmatrix_ = Teuchos::rcp(
        new Xpetra::BlockedCrsMatrix<SC, LO, GO, NO>(map_extractor, map_extractor, 42));

    for (int row = 0; row < A->rows(); row++)
    {
      for (int col = 0; col < A->cols(); col++)
      {
        Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO>> crsA =
            Teuchos::rcp(new EpetraCrsMatrix(A->matrix(row, col).epetra_matrix()));
        Teuchos::RCP<Xpetra::Matrix<SC, LO, GO, NO>> mat =
            Xpetra::MatrixFactory<SC, LO, GO, NO>::BuildCopy(
                Teuchos::rcp(new Xpetra::CrsMatrixWrap<SC, LO, GO, NO>(crsA)));

        pmatrix_->setMatrix(row, col, mat);
      }
    }
    pmatrix_->fillComplete();

    if (!tekolist_.sublist("Teko Parameters").isParameter("TEKO_XML_FILE"))
      FOUR_C_THROW("TEKO_XML_FILE parameter not set!");
    std::string xmlFileName =
        tekolist_.sublist("Teko Parameters").get<std::string>("TEKO_XML_FILE");

    Teuchos::RCP<Teuchos::ParameterList> tekoParams = Teuchos::rcp(new Teuchos::ParameterList());
    auto comm = pmatrix_->getRowMap()->getComm();
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, tekoParams.ptr(), *comm);

    // Check if multigrid is used as preconditioner for single field inverse approximation and
    // attach nullspace and coordinate information to the respective inverse parameter list.
    for (size_t block = 0; block < pmatrix_->Rows(); block++)
    {
      // Get the single field preconditioner sub-list of a matrix block hardwired under
      // "Inverse<1...n>".
      std::string inverse = "Inverse" + std::to_string(block + 1);
      Teuchos::ParameterList& inverseList = tekolist_.sublist(inverse);

      if (tekoParams->sublist(inverse).get<std::string>("Type") == "MueLu")
      {
        const int number_of_equations = inverseList.get<int>("PDE equations");

        Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> nullspace = Teuchos::rcp(
            new EpetraMultiVector(inverseList.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace")));

        Teuchos::RCP<Xpetra::MultiVector<SC, LO, GO, NO>> coordinates =
            Teuchos::rcp(new EpetraMultiVector(
                inverseList.get<Teuchos::RCP<Epetra_MultiVector>>("Coordinates")));

        tekoParams->sublist(inverse).set("number of equations", number_of_equations);
        Teuchos::ParameterList& userParamList = tekoParams->sublist(inverse).sublist("user data");
        userParamList.set("Nullspace", nullspace);
        userParamList.set("Coordinates", coordinates);
      }
    }

    Teuchos::RCP<Stratimikos::LinearSolverBuilder<double>> builder =
        Teuchos::rcp(new Stratimikos::DefaultLinearSolverBuilder);
    Stratimikos::enableMueLu<SC, LO, GO, NO>(*builder);

    Teuchos::RCP<Teko::InverseLibrary> invLib =
        Teko::InverseLibrary::buildFromParameterList(*tekoParams, builder);

    // Get the block preconditioner definition parameter sub-list hardwired under "Preconditioner".
    Teuchos::RCP<const Teko::InverseFactory> inverse = invLib->getInverseFactory("Preconditioner");
    Teuchos::RCP<const Thyra::LinearOpBase<double>> inverseOp =
        Teko::buildInverse(*inverse, pmatrix_->getThyraOperator());

    p_ = Teuchos::rcp(new Teko::Epetra::EpetraInverseOpWrapper(inverseOp));
  }
}

FOUR_C_NAMESPACE_CLOSE
