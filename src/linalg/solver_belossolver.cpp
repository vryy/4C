/*
 * solver_belossolver.cpp
 *
 *  Created on: Jul 5, 2011
 *      Author: wiesner
 */

#include "solver_belossolver.H"

#include "solver_pointpreconditioner.H"
#include "solver_blockpreconditioners.H"
#include "solver_krylovprojectionpreconditioner.H"
#include "solver_ifpackpreconditioner.H"
#include "solver_mlpreconditioner.H"

#include "linalg_solver.H"

void LINALG::Solver::BuildBelosSolver(const Epetra_Comm & comm, Teuchos::ParameterList & params, FILE * outfile)
{
  cout << "build Belos solver" << endl;
  solver_ = Teuchos::rcp( new LINALG::SOLVER::BelosSolver( comm, params, outfile));
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::BelosSolver::BelosSolver( const Epetra_Comm & comm,
                                            Teuchos::ParameterList & params,
                                            FILE * outfile )
  : comm_( comm ),
    params_( params ),
    outfile_( outfile ),
    ncall_( 0 )
{
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
void LINALG::SOLVER::BelosSolver::Setup(  Teuchos::RCP<Epetra_Operator>     matrix             ,
                                          Teuchos::RCP<Epetra_Vector>       x                  ,
                                          Teuchos::RCP<Epetra_Vector>       b                  ,
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
  Teuchos::RCP<Belos::SolverManager<double,MV,OP> > newSolver = Teuchos::rcp(new Belos::BlockGmresSolMgr<double,MV,OP>(problem,Teuchos::rcp(&belist,false)));

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

}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::BelosSolver::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)
{
  return preconditioner_->ApplyInverse( X, Y );
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::BelosSolver::CreatePreconditioner( Teuchos::ParameterList & azlist,
    bool isCrsMatrix,
    Teuchos::RCP<Epetra_MultiVector> weighted_basis_mean,
    Teuchos::RCP<Epetra_MultiVector> kernel_c,
    bool project )
{
  preconditioner_ = Teuchos::null;

  if ( isCrsMatrix )
  {
    // get type of preconditioner and build either Ifpack or ML
    // if we have an ifpack parameter list, we do ifpack
    // if we have an ml parameter list we do ml
    // if we have a downwinding flag we downwind the linear problem
    if ( Params().isSublist("IFPACK Parameters") )
    {
      preconditioner_ = Teuchos::rcp( new IFPACKPreconditioner( outfile_, Params().sublist("IFPACK Parameters"), azlist ) );
    }
    else if ( Params().isSublist("ML Parameters") )
    {
      preconditioner_ = Teuchos::rcp( new MLPreconditioner( outfile_, Params().sublist("ML Parameters") ) );
    }
    else if (azlist.get<int>("AZ_precond") == AZ_none)
    {
      preconditioner_ = Teuchos::rcp( new NonePreconditioner( outfile_, Params() ) );
    }
    else
    {
      dserror( "unknown preconditioner" );
    }

    // decide whether we do what kind of scaling
    string scaling = azlist.get("scaling","none");
    if (scaling=="none")
    {
    }
    else if (scaling=="infnorm")
    {
      preconditioner_ = Teuchos::rcp( new InfNormPreconditioner( preconditioner_ ) );
    }
    else if (scaling=="symmetric")
    {
      preconditioner_ = Teuchos::rcp( new SymDiagPreconditioner( preconditioner_ ) );
    }
    else
      dserror("Unknown type of scaling found in parameter list");

    if ( azlist.get<bool>("downwinding",false) )
    {
      preconditioner_ = Teuchos::rcp( new DWindPreconditioner( outfile_, preconditioner_, azlist ) );
    }

    if ( project )
    {
      preconditioner_ = Teuchos::rcp( new KrylovProjectionPreconditioner( outfile_, preconditioner_, weighted_basis_mean, kernel_c ) );
    }
  }
  else
  {
    // assume block matrix

    if ( Params().isSublist("SIMPLER") )
    {
      preconditioner_ = Teuchos::rcp( new SimplePreconditioner( outfile_, Params(), Params().sublist("SIMPLER") ) );
    }
    else if ( Params().isSublist("BGS Parameters") )
    {
      preconditioner_ = Teuchos::rcp( new BGSPreconditioner( outfile_, Params(), Params().sublist("BGS Parameters") ) );
    }
    else if ( Params().isSublist("AMGBS Parameters") )
    {
      preconditioner_ = Teuchos::rcp( new AMGBSPreconditioner( outfile_, Params() ) );
    }
    else
    {
      dserror( "unknown preconditioner" );
    }
  }

#if 0
  preconditioner_->Print( std::cout );
  std::cout << "\n";
#endif
}
