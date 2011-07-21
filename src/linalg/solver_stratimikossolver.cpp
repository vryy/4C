/*
 * solver_stratimikossolver.cpp
 *
 *  Created on: Jul 13, 2011
 *      Author: wiesner
 */

#ifdef TRILINOS_DEV

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_CrsMatrix.h>

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_LinearOpBase_decl.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>

#include "../drt_lib/drt_dserror.H"
#include "solver_stratimikossolver.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::StratimikosSolver::StratimikosSolver( const Epetra_Comm & comm,
                                            Teuchos::ParameterList & params,
                                            FILE * outfile )
  : comm_( comm ),
    params_( params ),
    outfile_( outfile ),
    ncall_( 0 )
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
void LINALG::SOLVER::StratimikosSolver::Setup( Teuchos::RCP<Epetra_Operator> matrix,
                                               Teuchos::RCP<Epetra_Vector> x,
                                               Teuchos::RCP<Epetra_Vector> b,
                                               bool refactor,
                                               bool reset,
                                               Teuchos::RCP<Epetra_MultiVector> weighted_basis_mean,
                                               Teuchos::RCP<Epetra_MultiVector> kernel_c,
                                               bool project)
{
  if (!Params().isSublist("Stratimikos Parameters"))
    dserror("Do not have stratimikos parameter list");

  x_ = x;
  b_ = b;
  A_ = matrix;

}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::StratimikosSolver::Solve()
{
  Teuchos::RCP<Teuchos::ParameterList> stratimikoslist = Teuchos::rcp(new Teuchos::ParameterList(Params().sublist("Stratimikos Parameters")));
      //Teuchos::rcp(&(Params().sublist("Stratimikos Parameters")));
  cout << "Stratimikos List from dat file" << endl;
  cout << *stratimikoslist << endl << endl;


  //std::string xmlfile = stratimikoslist.get<std::string>("xml file");

  //Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder("stratimikos_KLU.xml");
  Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
  linearSolverBuilder.setParameterList(stratimikoslist);
  //Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder(xmlfile);
  linearSolverBuilder.readParameters(&cout);
  cout << *(linearSolverBuilder.getParameterList()) << endl;

  // see whether operator is a Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> epetra_A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(A_);

  // wrap Epetra -> Thyra
  Teuchos::RCP<const Thyra::LinearOpBase<double> > A = Thyra::epetraLinearOp( epetra_A );
  Teuchos::RCP<Thyra::VectorBase<double> > x = Thyra::create_Vector( x_, A->domain());
  Teuchos::RCP<const Thyra::VectorBase<double> > b = Thyra::create_Vector( b_, A->range());

  // Create a linear solver factory given information read from the
  // parameter list.
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
    linearSolverBuilder.createLinearSolveStrategy("");

  // Setup output stream and the verbosity level
  //lowsFactory->setOStream(&out);
  //lowsFactory->setVerbLevel(Teuchos::VERB_LOW);

  // Create a linear solver based on the forward operator A
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > lows =
    Thyra::linearOpWithSolve(*lowsFactory, A);

  // Solve the linear system (note: the initial guess in 'x' is critical)
  Thyra::SolveStatus<double> status =
       Thyra::solve<double>(*lows, Thyra::NOTRANS, *b, x.ptr());
  cout << "\nSolve status:\n" << status << endl;


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
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::StratimikosSolver::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)
{
  dserror("ApplyInverse not implemented for StratimikosSolver");
  return -1;
}


#endif /* TRILINOS_DEV */
