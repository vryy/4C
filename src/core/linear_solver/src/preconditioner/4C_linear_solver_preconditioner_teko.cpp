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
#include <Teko_StratimikosFactory.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

FOUR_C_NAMESPACE_OPEN

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
  using EpetraMultiVector = Xpetra::EpetraMultiVectorT<GlobalOrdinal, Node>;
  using XpetraMultiVector = Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  if (create)
  {
    if (!tekolist_.sublist("Teko Parameters").isParameter("TEKO_XML_FILE"))
      FOUR_C_THROW("TEKO_XML_FILE parameter not set!");
    std::string xmlFileName =
        tekolist_.sublist("Teko Parameters").get<std::string>("TEKO_XML_FILE");

    Teuchos::RCP<Teuchos::ParameterList> tekoParams = Teuchos::make_rcp<Teuchos::ParameterList>();
    auto comm = Core::Communication::to_teuchos_comm<int>(matrix->Comm());
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, tekoParams.ptr(), *comm);

    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> A =
        Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(Teuchos::RCP(matrix, false));

    // wrap linear operators
    if (A.is_null())
    {
      auto A_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Teuchos::RCP(matrix, false));
      pmatrix_ = Thyra::epetraLinearOp(A_crs);
    }
    else
    {
      pmatrix_ = Thyra::defaultBlockedLinearOp<double>();

      Teko::toBlockedLinearOp(pmatrix_)->beginBlockFill(A->rows(), A->cols());
      for (int row = 0; row < A->rows(); row++)
      {
        for (int col = 0; col < A->cols(); col++)
        {
          auto A_crs = Teuchos::make_rcp<Epetra_CrsMatrix>(*A->matrix(row, col).epetra_matrix());
          Teko::toBlockedLinearOp(pmatrix_)->setBlock(row, col, Thyra::epetraLinearOp(A_crs));
        }
      }
      Teko::toBlockedLinearOp(pmatrix_)->endBlockFill();

      // check if multigrid is used as preconditioner for single field inverse approximation and
      // attach nullspace and coordinate information to the respective inverse parameter list.
      for (int block = 0; block < A->rows(); block++)
      {
        std::string inverse = "Inverse" + std::to_string(block + 1);

        if (tekolist_.isSublist(inverse))
        {
          // get the single field preconditioner sub-list of a matrix block hardwired under
          // "Inverse<1...n>".
          Teuchos::ParameterList& inverseList = tekolist_.sublist(inverse);

          if (tekoParams->sublist("Inverse Factory Library")
                  .sublist(inverse)
                  .get<std::string>("Type") == "MueLu")
          {
            const int number_of_equations = inverseList.get<int>("PDE equations");

            Teuchos::RCP<XpetraMultiVector> nullspace = Teuchos::make_rcp<EpetraMultiVector>(
                inverseList.get<Teuchos::RCP<Epetra_MultiVector>>("nullspace"));

            Teuchos::RCP<XpetraMultiVector> coordinates = Teuchos::make_rcp<EpetraMultiVector>(
                inverseList.get<Teuchos::RCP<Epetra_MultiVector>>("Coordinates"));

            tekoParams->sublist("Inverse Factory Library")
                .sublist(inverse)
                .set("number of equations", number_of_equations);
            Teuchos::ParameterList& userParamList = tekoParams->sublist("Inverse Factory Library")
                                                        .sublist(inverse)
                                                        .sublist("user data");
            userParamList.set("Nullspace", nullspace);
            userParamList.set("Coordinates", coordinates);
          }
        }
      }
    }

    // setup preconditioner builder and enable relevant packages
    Teuchos::RCP<Stratimikos::LinearSolverBuilder<double>> builder =
        Teuchos::make_rcp<Stratimikos::DefaultLinearSolverBuilder>();

    Stratimikos::enableMueLu<Scalar, LocalOrdinal, GlobalOrdinal, Node>(*builder);
    Teko::addTekoToStratimikosBuilder(*builder);

    Teuchos::RCP<Teuchos::ParameterList> stratimikos_params =
        Teuchos::make_rcp<Teuchos::ParameterList>(*builder->getValidParameters());
    Teuchos::ParameterList& tekoList =
        stratimikos_params->sublist("Preconditioner Types").sublist("Teko");
    tekoList.setParameters(*tekoParams);
    builder->setParameterList(stratimikos_params);

    Teuchos::RCP<Thyra::PreconditionerFactoryBase<double>> precFactory =
        builder->createPreconditioningStrategy("Teko");
    Teuchos::RCP<Thyra::PreconditionerBase<double>> prec =
        Thyra::prec<double>(*precFactory, pmatrix_);
    Teko::LinearOp inverseOp = prec->getUnspecifiedPrecOp();

    p_ = Teuchos::make_rcp<Teko::Epetra::EpetraInverseOpWrapper>(inverseOp);
  }
}

FOUR_C_NAMESPACE_CLOSE
