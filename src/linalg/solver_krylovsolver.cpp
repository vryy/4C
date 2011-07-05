/*
 * solver_krylovsolver.cpp
 *
 *  Created on: Jul 4, 2011
 *      Author: wiesner
 */

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include "../drt_lib/drt_dserror.H"
#include "solver_krylovsolver.H"

#include "solver_pointpreconditioner.H"
#include "solver_blockpreconditioners.H"
#include "solver_krylovprojectionpreconditioner.H"
#include "solver_ifpackpreconditioner.H"
#include "solver_mlpreconditioner.H"


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::KrylovSolver::KrylovSolver( const Epetra_Comm & comm,
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
LINALG::SOLVER::KrylovSolver::~KrylovSolver()
{
  preconditioner_ = Teuchos::null;
  A_ = Teuchos::null;
  x_ = Teuchos::null;
  b_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::KrylovSolver::Setup( Teuchos::RCP<Epetra_Operator> matrix,
                                          Teuchos::RCP<Epetra_Vector> x,
                                          Teuchos::RCP<Epetra_Vector> b,
                                          bool refactor,
                                          bool reset,
                                          Teuchos::RCP<Epetra_MultiVector> weighted_basis_mean,
                                          Teuchos::RCP<Epetra_MultiVector> kernel_c,
                                          bool project)
{
#ifdef WRITEOUTSTATISTICS
  dtimeprecondsetup_ = 0.;
  Epetra_Time tttcreate(Comm()); // time measurement for creation of preconditioner
#endif

  if (!Params().isSublist("Aztec Parameters"))
    dserror("Do not have aztec parameter list");
  Teuchos::ParameterList& azlist = Params().sublist("Aztec Parameters");
  //int azoutput = azlist.get<int>("AZ_output",0);

  x_ = x;
  b_ = b;
  A_ = matrix;

  // see whether operator is a Epetra_CrsMatrix
  Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>( A_.get() );

#ifdef WRITEOUTSTATISTICS
  tttcreate.ResetStartTime();
#endif

  // decide whether we recreate preconditioners
  int  reuse  = azlist.get("reuse",0);
  bool create = reset or not Ncall() or not reuse or ( Ncall() % reuse )==0;
  if ( create )
  {
    ncall_ = 0;
    CreatePreconditioner( azlist, A!=NULL, weighted_basis_mean, kernel_c, project );
  }

  preconditioner_->Setup( create, &*A_, &*x_, &*b_ );

#ifdef WRITEOUTSTATISTICS
  dtimeprecondsetup = tttcreate.ElapsedTime();
#endif
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::KrylovSolver::Solve()
{
#ifdef WRITEOUTSTATISTICS
  Epetra_Time ttt(Comm());       // time measurement for whole routine
  ttt.ResetStartTime();
#endif

  Teuchos::ParameterList& azlist = Params().sublist("Aztec Parameters");

  // Allocate an aztec solver with default parameters
  // We do this every time because reusing the solver object
  // does lead to crashes that are not understood

  // create an aztec solver
  AztecOO aztec;
  aztec.SetAztecDefaults();

  // tell aztec to which stream to write
  aztec.SetOutputStream(std::cout);
  aztec.SetErrorStream(std::cerr);

  // Don't want linear problem to alter our aztec parameters (idiot feature!)
  // this is why we set our list here AFTER the linear problem has been set
  aztec.SetProblem( preconditioner_->LinearProblem() );

  {
    // We don't want to use Aztec's scaling capabilities as we prefer to do
    // the scaling ourselves (so we precisely know what happens)
    // Therefore set scaling parameter to none and reset it after aztec has made
    // its internal copy of the parameter list
    string scaling = azlist.get("scaling","none");
    azlist.set("scaling","none");
    aztec.SetParameters(azlist,false);
    azlist.set("scaling",scaling);
  }

  aztec.SetPrecOperator( preconditioner_->PrecOperator() );

  // iterate on the solution
  int iter = azlist.get("AZ_max_iter",500);
  double tol = azlist.get("AZ_tol",1.0e-6);

  // This hurts! It supresses error messages. This needs to be fixed.
#if 0
  // create an aztec convergence test as combination of
  // L2-norm and Inf-Norm to be both satisfied where we demand
  // L2 < tol and Linf < 10*tol
  {
    Epetra_Operator* op  = aztec.GetProblem()->GetOperator();
    Epetra_Vector*   rhs = static_cast<Epetra_Vector*>(aztec.GetProblem()->GetRHS());
    Epetra_Vector*   lhs = static_cast<Epetra_Vector*>(aztec.GetProblem()->GetLHS());
    // max iterations
    aztest_maxiter_ = rcp(new AztecOO_StatusTestMaxIters(iter));
    // L2 norm
    aztest_norm2_ = rcp(new AztecOO_StatusTestResNorm(*op,*lhs,*rhs,tol));
    aztest_norm2_->DefineResForm(AztecOO_StatusTestResNorm::Implicit,
                                 AztecOO_StatusTestResNorm::TwoNorm);
    aztest_norm2_->DefineScaleForm(AztecOO_StatusTestResNorm::NormOfInitRes,
                                   AztecOO_StatusTestResNorm::TwoNorm);
    // Linf norm (demanded to be 1.0 times L2-norm now, to become an input parameter?)
    aztest_norminf_ = rcp(new AztecOO_StatusTestResNorm(*op,*lhs,*rhs,1.0*tol));
    aztest_norminf_->DefineResForm(AztecOO_StatusTestResNorm::Implicit,
                                   AztecOO_StatusTestResNorm::InfNorm);
    aztest_norminf_->DefineScaleForm(AztecOO_StatusTestResNorm::NormOfInitRes,
                                     AztecOO_StatusTestResNorm::InfNorm);
    // L2 AND Linf
    aztest_combo1_ = rcp(new AztecOO_StatusTestCombo(AztecOO_StatusTestCombo::SEQ));
    // maxiters OR (L2 AND Linf)
    aztest_combo2_ = rcp(new AztecOO_StatusTestCombo(AztecOO_StatusTestCombo::OR));
    aztest_combo1_->AddStatusTest(*aztest_norm2_);
    aztest_combo1_->AddStatusTest(*aztest_norminf_);
    aztest_combo2_->AddStatusTest(*aztest_maxiter_);
    aztest_combo2_->AddStatusTest(*aztest_combo1_);
    // set status test
    aztec.SetStatusTest(aztest_combo2_.get());
  }
#endif

  // if you want to get some information on eigenvalues of the Hessenberg matrix/the
  // estimated condition number of the preconditioned system, uncomment the following
  // line and set AZOUTPUT>0 in your .dat-file
  // aztec_->SetAztecOption(AZ_solver,AZ_gmres_condnum);

  //------------------------------- just do it----------------------------------------
  aztec.Iterate(iter,tol);
  //----------------------------------------------------------------------------------

  preconditioner_->Finish( &*A_, &*x_, &*b_ );

  // check status of solution process
  const double* status = aztec.GetAztecStatus();
#if 0
  AztecOO_StatusType stat = aztest_combo2_->GetStatus();
  if (stat!=Converged)
  {
    bool resolve = false;
    if (stat==Unconverged)
    {
      if (comm_.MyPID()==0) printf("Max iterations reached in AztecOO\n");
    }
    else if (stat==Failed || stat==NaN || stat==PartialFailed)
    {
      if (comm_.MyPID()==0) printf("Numerical breakdown in AztecOO\n");
    }
    else dserror("Aztec returned unknown nonzero status %d",(int)stat);
  }
#else
  if (status[AZ_why] != AZ_normal)
  {
    if (status[AZ_why] == AZ_breakdown)
    {
      if (comm_.MyPID()==0) printf("Numerical breakdown in AztecOO\n");
    }
    else if (status[AZ_why] == AZ_ill_cond)
    {
      if (comm_.MyPID()==0) printf("Problem is near singular in AztecOO\n");
    }
    else if (status[AZ_why] == AZ_maxits)
    {
      if (comm_.MyPID()==0) printf("Max iterations reached in AztecOO\n");
    }
  } // if (status[AZ_why] != AZ_normal)
#endif

#ifdef WRITEOUTSTATISTICS
    if(outfile_)
    {
      fprintf(outfile_,"LinIter %i\tNumGlobalElements %i\tAZ_solve_time %f\tAztecSolveTime %f\tAztecPrecondSetup %f\t\n",
              (int)status[AZ_its],
              A_->OperatorRangeMap().NumGlobalElements(),
              status[AZ_solve_time],
              dtimeprecondsetup_ + ttt.ElapsedTime(),
              dtimeprecondsetup_);
      fflush(outfile_);
    }
#endif

  // print some output if desired
  if (comm_.MyPID()==0 && outfile_)
  {
    fprintf(outfile_,"AztecOO: unknowns/iterations/time %d  %d  %f\n",
            A_->OperatorRangeMap().NumGlobalElements(),(int)status[AZ_its],status[AZ_solve_time]);
    fflush(outfile_);
  }

  ncall_ += 1;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::KrylovSolver::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)
{
  return preconditioner_->ApplyInverse( X, Y );
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::KrylovSolver::CreatePreconditioner( Teuchos::ParameterList & azlist,
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
      preconditioner_ = Teuchos::rcp( new LINALG::SOLVER::IFPACKPreconditioner( outfile_, Params().sublist("IFPACK Parameters"), azlist ) );
    }
    else if ( Params().isSublist("ML Parameters") )
    {
      preconditioner_ = Teuchos::rcp( new LINALG::SOLVER::MLPreconditioner( outfile_, Params().sublist("ML Parameters") ) );
    }
    else if ( Params().isSublist("AMGBS Parameters") )
    {
      preconditioner_ = Teuchos::rcp( new LINALG::SOLVER::AMGBSPreconditioner( outfile_, Params() ) );
    }
    else if (azlist.get<int>("AZ_precond") == AZ_none)
    {
      preconditioner_ = Teuchos::rcp( new LINALG::SOLVER::NonePreconditioner( outfile_, Params() ) );
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
      preconditioner_ = Teuchos::rcp( new LINALG::SOLVER::InfNormPreconditioner( preconditioner_ ) );
    }
    else if (scaling=="symmetric")
    {
      preconditioner_ = Teuchos::rcp( new LINALG::SOLVER::SymDiagPreconditioner( preconditioner_ ) );
    }
    else
      dserror("Unknown type of scaling found in parameter list");

    if ( azlist.get<bool>("downwinding",false) )
    {
      preconditioner_ = Teuchos::rcp( new LINALG::SOLVER::DWindPreconditioner( outfile_, preconditioner_, azlist ) );
    }

    if ( project )
    {
      preconditioner_ = Teuchos::rcp( new LINALG::SOLVER::KrylovProjectionPreconditioner( outfile_, preconditioner_, weighted_basis_mean, kernel_c ) );
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
