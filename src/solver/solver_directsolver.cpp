/*
 * solver_directsolver.cpp
 *
 *  Created on: Jul 4, 2011
 *      Author: wiesner
 */

// Amesos headers
#include <Amesos_Klu.h>
#include <Amesos_Lapack.h>
#ifndef HAVENOT_UMFPACK
#include <Amesos_Umfpack.h>
#endif
#ifndef HAVENOT_SUPERLU
#include <Amesos_Superludist.h>
#endif

// EpetraExt headers
#include <EpetraExt_Reindex_LinearProblem.h>

// BACI headers
#include "solver_directsolver.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::DirectSolver::DirectSolver( std::string solvertype )
  : solvertype_( solvertype ),
    factored_( false )
{
  lp_ = Teuchos::rcp( new Epetra_LinearProblem() );
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::DirectSolver::~DirectSolver()
{
  amesos_    = Teuchos::null;
  reindexer_ = Teuchos::null;
  lp_        = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::DirectSolver::Setup( RCP<Epetra_Operator> matrix,
                                          RCP<Epetra_MultiVector> x,
                                          RCP<Epetra_MultiVector> b,
                                          bool refactor,
                                          bool reset,
                                          RCP<Epetra_MultiVector> weighted_basis_mean,
                                          RCP<Epetra_MultiVector> kernel_c,
                                          bool project)
{
  if ( project )
    dserror( "a projection of Krylov space basis vectors is possible only for aztec type iterative solvers" );

  x_ = x;
  b_ = b;
  A_ = matrix;

  // fill the linear problem
  lp_->SetRHS(b_.get());
  lp_->SetLHS(x_.get());
  lp_->SetOperator(A_.get());

  if (reset or refactor or not IsFactored())
  {
    reindexer_ = Teuchos::rcp(new EpetraExt::LinearProblem_Reindex(NULL));

    if ( solvertype_=="klu" )
    {
      amesos_ = Teuchos::rcp(new Amesos_Klu((*reindexer_)(*lp_)));
    }
    else if ( solvertype_=="umfpack" )
    {
#ifndef HAVENOT_UMFPACK
      amesos_ = Teuchos::rcp(new Amesos_Umfpack((*reindexer_)(*lp_)));
#else
      dserror( "no umfpack here" );
#endif
    }
    else if ( solvertype_=="superlu" )
    {
#ifndef HAVENOT_SUPERLU
#ifdef PARALLEL
      amesos_ = Teuchos::rcp(new Amesos_Superludist((*reindexer_)(*lp_)));
#else
      dserror("Distributed SuperLU only in parallel");
#endif
#else
      dserror( "no superlu_dist here" );
#endif
    }
    else if ( solvertype_=="lapack" )
    {
      amesos_ = Teuchos::rcp(new Amesos_Lapack((*reindexer_)(*lp_)));
    }

    factored_ = false;
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::DirectSolver::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y)
{
  x_->Update( 1., X, 0. );
  Solve();
  Y.Update( 1., *b_, 0. );

  return 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::DirectSolver::Solve()
{
  if (amesos_==Teuchos::null) dserror("No solver allocated");

  // Problem has not been factorized before
  if (not IsFactored())
  {
    int err = amesos_->SymbolicFactorization();
    if (err) dserror("Amesos::SymbolicFactorization returned an err");
    err = amesos_->NumericFactorization();
    if (err) dserror("Amesos::NumericFactorization returned an err");

    factored_ = true;
  }

  int err = amesos_->Solve();
  if (err) dserror("Amesos::Solve returned an err");
}




