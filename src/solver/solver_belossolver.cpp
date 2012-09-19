/*
 * solver_belossolver.cpp
 *
 *  Created on: Jul 5, 2011
 *      Author: wiesner
 */

// Belos headers
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"

// BACI headers
#include "solver_belossolver.H"
#include "solver_pointpreconditioner.H"
#include "solver_blockpreconditioners.H"
#include "solver_krylovprojectionpreconditioner.H"
#include "solver_ifpackpreconditioner.H"
#include "solver_mlpreconditioner.H"

#include "../linalg/linalg_solver.H"

void LINALG::Solver::BuildBelosSolver(const Epetra_Comm & comm, Teuchos::ParameterList & params, FILE * outfile)
{
  solver_ = Teuchos::rcp( new LINALG::SOLVER::BelosSolver( comm, params, outfile));
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::BelosSolver::BelosSolver( const Epetra_Comm & comm,
                                            Teuchos::ParameterList & params,
                                            FILE * outfile )
 : KrylovSolver(comm,params,outfile)
{
  ncall_ = 0;
  preconditioner_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::BelosSolver::~BelosSolver()
{
  preconditioner_ = Teuchos::null;
  A_ = Teuchos::null;
  x_ = Teuchos::null;
  b_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::BelosSolver::Setup(  Teuchos::RCP<Epetra_Operator>     matrix            ,
                                          Teuchos::RCP<Epetra_MultiVector>       x            ,
                                          Teuchos::RCP<Epetra_MultiVector>       b            ,
                                          bool                             refactor           ,
                                          bool                             reset              ,
                                          Teuchos::RCP<Epetra_MultiVector>  weighted_basis_mean,
                                          Teuchos::RCP<Epetra_MultiVector>  kernel_c           ,
                                          bool                             project)
{
  x_ = x;
  b_ = b;
  A_ = matrix;

  // see whether operator is a Epetra_CrsMatrix
  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>( A_.get() );

  if (!Params().isSublist("Belos Parameters"))
    dserror("Do not have belos parameter list");
  Teuchos::ParameterList& belist = Params().sublist("Belos Parameters");

  int reuse = belist.get("reuse",0); // TODO: fix me!
  bool create = not Ncall() or not reuse or (Ncall() % reuse ) == 0;
  if (create)
  {
    ncall_ = 0;
    CreatePreconditioner(belist, A!=NULL, weighted_basis_mean, kernel_c, project );
  }

  preconditioner_->Setup(create,&*A_, &*x_,&*b_);

}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::BelosSolver::Solve()
{
  Teuchos::ParameterList& belist = Params().sublist("Belos Parameters");

  typedef double                            ST;
  typedef Teuchos::ScalarTraits<ST>        SCT;
  typedef SCT::magnitudeType                MT;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  typedef Belos::MultiVecTraits<ST,MV>     MVT;
  typedef Belos::OperatorTraits<ST,MV,OP>  OPT;

  // build Belos linear problem
  Teuchos::RCP<Belos::LinearProblem<double, MV, OP> > problem = Teuchos::rcp(new Belos::LinearProblem<double,MV,OP>(A_, x_, b_) );
  // TODO support for left preconditioner?
  if (preconditioner_ != Teuchos::null)
  {
    // prepare preconditioner in preconditioner_->PrecOperator() for Belos
    Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp(new Belos::EpetraPrecOp(Teuchos::rcp(preconditioner_->PrecOperator(),false)));
    problem->setRightPrec(belosPrec);
  }
  bool set = problem->setProblem();
  if (set == false)
  {
    std::cout << std::endl << "ERROR: Belos::LinearProblem failed to set up correctly!" << std::endl;
  }

  // create iterative solver manager
  Teuchos::RCP<Belos::SolverManager<double,MV,OP> > newSolver;
  std::string solverType = belist.get<std::string>("Solver Type");
  if(solverType=="GMRES")
     newSolver = Teuchos::rcp(new Belos::BlockGmresSolMgr<double,MV,OP>(problem,Teuchos::rcp(&belist,false)));
  else if (solverType=="CG")
     newSolver = Teuchos::rcp(new Belos::BlockCGSolMgr<double,MV,OP>(problem,Teuchos::rcp(&belist,false)));
  else
    dserror("unknown solver type for Belos");

  //
  // Perform solve
  //
  Belos::ReturnType ret = newSolver->solve();

  // TODO: check me -> access solution x from linear problem???
  if(preconditioner_!=Teuchos::null)
    preconditioner_->Finish( &*A_, &*x_, &*b_ );

  if (ret!=Belos::Converged)
  {
    std::cout << std::endl << "WARNING: Belos did not converge!" << std::endl;
  }

  ncall_ += 1; // increment counter of solver calls
}
