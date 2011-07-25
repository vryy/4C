/*
 * solver_tekopreconditioner.cpp
 *
 *  Created on: Jul 11, 2011
 *      Author: wiesner
 */

#ifdef TRILINOS_DEV

#include "../drt_lib/drt_dserror.H"
#include "solver_tekopreconditioner.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::TekoPreconditioner::TekoPreconditioner( FILE * outfile,
                                                        Teuchos::ParameterList & tekolist)
  : PreconditionerType( outfile ),
    tekolist_( tekolist )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::TekoPreconditioner::Setup( bool create,
                                                Epetra_Operator * matrix,
                                                Epetra_MultiVector * x,
                                                Epetra_MultiVector * b )
{
  SetupLinearProblem( matrix, x, b );

  if ( create )
  {
    Epetra_CrsMatrix* A = dynamic_cast<Epetra_CrsMatrix*>( matrix );
    if ( A==NULL )
      dserror( "CrsMatrix expected" );

    // free old matrix first
    prec_    = Teuchos::null;
    Pmatrix_ = Teuchos::null;

    // create a copy of the scaled matrix
    // so we can reuse the preconditioner
    Pmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*A));


  }
}


#endif /* TRILINOS_DEV */
