/*----------------------------------------------------------------------*/
/*! \file

\brief Interface class for MueLu preconditioner

\level 1

*/

// Baci
#include "solver_muelupreconditioner.H"
#include "muelu/muelu_utils.H"

#include "../drt_lib/drt_dserror.H"

// Teuchos
#include <Teuchos_RCP.hpp>
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
#include <MueLu_EpetraOperator.hpp>   // Aztec interface
#include <MueLu_UseDefaultTypes.hpp>  // => Scalar=double, LocalOrdinal=GlobalOrdinal=int

// Xpetra
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_EpetraMap.hpp>
#include <Xpetra_IO.hpp>
#include <Xpetra_Map.hpp>

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::MueLuPreconditioner::MueLuPreconditioner(
    FILE* outfile, Teuchos::ParameterList& mllist)
    : PreconditionerType(outfile), mllist_(mllist)
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
#ifdef TRILINOS_DEVELOP
      Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>(Pmatrix_));
#else
      Teuchos::rcp(new Xpetra::EpetraCrsMatrix(Pmatrix_));
#endif
  Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mueluOp =
      Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(mueluA));

  // remove unsupported flags, keep for now ... is used in many tests but everywhere 0
  mllist_.remove("aggregation: threshold", false);  // no support for aggregation: threshold TODO

  if (mllist_.get<bool>("MUELU_XML_ENFORCE"))
  {
    if (create)
    {
      if (!mllist_.isParameter("MUELU_XML_FILE"))
        dserror(
            "XML-file w/ MueLu preconditioner configuration is missing in solver parameter list. "
            "Please set it as entry 'MUELU_XML_FILE'.");
      std::string xmlFileName = mllist_.get<std::string>("MUELU_XML_FILE");

      // prepare nullspace vector for MueLu
      int numdf = mllist_.get<int>("PDE equations", -1);
      int dimns = mllist_.get<int>("null space: dimension", -1);
      if (dimns == -1 || numdf == -1)
        dserror("Error: PDE equations or null space dimension wrong.");

      Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap =
          mueluA->getRowMap();
      Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> nullspace;
      nullspace = LINALG::SOLVER::MUELU::UTILS::ExtractNullspaceFromMLList(rowMap, mllist_);

      mueluOp->SetFixedBlockSize(numdf);

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
        mllist_ /*, vec*/);
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