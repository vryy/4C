/*----------------------------------------------------------------------*/
/*! \file

\brief Solver class using Trilinos STRATIMIKOS framework

\level 1


Created on: Jul 13, 2011
*----------------------------------------------------------------------*/

#ifdef HAVE_Stratimikos

#include <Epetra_Comm.h>
//#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
//#endif
#include <Epetra_SerialComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_LinearOpBase_decl.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_DefaultSpmdVectorSpaceFactory_def.hpp>

#include "Teuchos_DefaultSerialComm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"

#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "solver_stratimikossolver.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::StratimikosSolver::StratimikosSolver(
    const Epetra_Comm& comm, Teuchos::ParameterList& params, FILE* outfile)
    : comm_(comm), params_(params), outfile_(outfile), ncall_(0)
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::StratimikosSolver::~StratimikosSolver()
{
  A_ = Teuchos::null;
  x_ = Teuchos::null;
  b_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::StratimikosSolver::Setup(Teuchos::RCP<Epetra_Operator> matrix,
    Teuchos::RCP<Epetra_MultiVector> x, Teuchos::RCP<Epetra_MultiVector> b, bool refactor,
    bool reset, Teuchos::RCP<LINALG::KrylovProjector> projector)
{
  if (!Params().isSublist("Stratimikos Parameters"))
    dserror("Do not have stratimikos parameter list");

  x_ = x;
  b_ = b;
  A_ = matrix;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::StratimikosSolver::Solve()
{
  Teuchos::RCP<Teuchos::ParameterList> stratimikoslist =
      Teuchos::rcp(new Teuchos::ParameterList(Params().sublist("Stratimikos Parameters")));
  // Teuchos::rcp(&(Params().sublist("Stratimikos Parameters")));
  std::cout << "Stratimikos List from dat file" << std::endl;
  std::cout << *stratimikoslist << std::endl << std::endl;


  // std::string xmlfile = stratimikoslist.get<std::string>("xml file");

  // Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder("stratimikos_KLU.xml");
  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
  linearSolverBuilder.setParameterList(stratimikoslist);
  // Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder(xmlfile);
  linearSolverBuilder.readParameters(&std::cout);
  std::cout << *(linearSolverBuilder.getParameterList()) << std::endl;

  // see whether operator is a Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> epetra_A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(A_);

  // create a dummy one-column Thyra::VectorSpaceBase
  Teuchos::RCP<const Teuchos::Comm<Teuchos::Ordinal>> TeuchosComm =
      COMM_UTILS::toTeuchosComm<Teuchos::Ordinal>(A_->Comm());

  Teuchos::RCP<Thyra::DefaultSpmdVectorSpaceFactory<double>> dummyDomainSpaceFac =
      Teuchos::rcp(new Thyra::DefaultSpmdVectorSpaceFactory<double>(TeuchosComm));

  Teuchos::RCP<const Thyra::VectorSpaceBase<double>> dummyDomainSpace =
      dummyDomainSpaceFac->createVecSpc(x_->NumVectors());

  // wrap Epetra -> Thyra
  Teuchos::RCP<const Thyra::LinearOpBase<double>> A = Thyra::epetraLinearOp(epetra_A);
  Teuchos::RCP<Thyra::MultiVectorBase<double>> x =
      Thyra::create_MultiVector(x_, A->domain(), dummyDomainSpace);
  Teuchos::RCP<const Thyra::MultiVectorBase<double>> b =
      Thyra::create_MultiVector(b_, A->range(), dummyDomainSpace);

  // Create a linear solver factory given information read from the
  // parameter list.
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double>> lowsFactory =
      linearSolverBuilder.createLinearSolveStrategy("");

  // Setup output stream and the verbosity level
  // lowsFactory->setOStream(&out);
  // lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

  // Create a linear solver based on the forward operator A
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double>> lows =
      Thyra::linearOpWithSolve(*lowsFactory, A);

  // Solve the linear system (note: the initial guess in 'x' is critical)
  Thyra::SolveStatus<double> status = Thyra::solve<double>(*lows, Thyra::NOTRANS, *b, x.ptr());
  std::cout << "\nSolve status:\n" << status << std::endl;


  // Wipe out the Thyra wrapper for x to guarantee that the solution will be
  // written back to epetra_x!  At the time of this writting this is not
  // really needed but the behavior may change at some point so this is a
  // good idea.
  x = Teuchos::null;
  b = Teuchos::null;


  // print some output if desired
  /*if (comm_.MyPID()==0 && outfile_)
  {
    fprintf(outfile_,"AztecOO: unknowns/iterations/time %d  %d  %f\n",
            A_->OperatorRangeMap().NumGlobalElements(),(int)status[AZ_its],status[AZ_solve_time]);
    fflush(outfile_);
  }*/

  ncall_ += 1;
  return 0;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::StratimikosSolver::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y)
{
  dserror("ApplyInverse not implemented for StratimikosSolver");
  return -1;
}

#endif
