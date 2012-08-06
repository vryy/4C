/*
 * solver_krylovsolver.cpp
 *
 *  Created on: Jul 4, 2011
 *      Author: wiesner
 */

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <az_aztec_defs.h> // for AZ_none (provokes compiler warning due to redeclaration of HAVE_SYS_TIME_H in mpi.h and AztecOO_config.h -> AztecOO problem)

#include "../drt_lib/drt_dserror.H"
#include "solver_krylovsolver.H"

#include "solver_pointpreconditioner.H"
#include "solver_blockpreconditioners.H"
#include "solver_krylovprojectionpreconditioner.H"
#include "solver_ifpackpreconditioner.H"
#include "solver_mlpreconditioner.H"
#include "solver_muelupreconditioner.H"
#include "solver_muelucontactpreconditioner.H"
#ifdef HAVE_TEKO
#include "solver_tekopreconditioner.H"
#endif

#include <Teuchos_TimeMonitor.hpp>

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
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::Solver:  1.1)   CreatePreconditioner");

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
    else if ( Params().isSublist("MueLu Parameters") )
    {
#ifdef HAVE_MueLu
      preconditioner_ = Teuchos::rcp( new LINALG::SOLVER::MueLuPreconditioner( outfile_, Params().sublist("MueLu Parameters") ) );
#else
      dserror("MueLu only available in DEV version of BACI with Trilinos Q1/2012 or newer.");
#endif
    }
    else if ( Params().isSublist("MueLu (Contact) Parameters") )
    {
#ifdef HAVE_MueLu
      preconditioner_ = Teuchos::rcp( new LINALG::SOLVER::MueLuContactPreconditioner( outfile_, Params().sublist("MueLu (Contact) Parameters") ) );
#else
      dserror("MueLu (Contact) preconditioner only available in DEV version of BACI with Trilinos Q1/2012 or newer.");
#endif
    }
    else if ( Params().isSublist("AMGBS Parameters") )
    {
      preconditioner_ = Teuchos::rcp( new LINALG::SOLVER::AMGBSPreconditioner( outfile_, Params() ) );
    }
    else if (azlist.get<int>("AZ_precond") == AZ_none)  // FIXME Attention: this is dangerous.
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

    if ( Params().isSublist("SIMPLER") ) // old BACI::(Cheap)SIMPLER preconditioner TODO: remove/replace me
    {
      dserror("SIMPLER sublist not supported any more.");
      preconditioner_ = Teuchos::rcp( new SimplePreconditioner( outfile_, Params() ) ); // Michael's SIMPLE for Fluid
    }
    else if( Params().isSublist("CheapSIMPLE Parameters") )
    {
      preconditioner_ = Teuchos::rcp( new SimplePreconditioner( outfile_, Params() ) );
    }
    else if ( Params().isSublist("BGS Parameters") )
    {
      preconditioner_ = Teuchos::rcp( new BGSPreconditioner( outfile_, Params(), Params().sublist("BGS Parameters") ) );
    }
    else if ( Params().isSublist("AMGBS Parameters") )
    {
      preconditioner_ = Teuchos::rcp( new AMGBSPreconditioner( outfile_, Params() ) );
    }
    else if ( Params().isSublist("Teko Parameters") )
    {
#ifdef HAVE_TEKO
      preconditioner_ = Teuchos::rcp( new TekoPreconditioner( outfile_, Params() ));
#else
      dserror("You need the HAVE_TEKO define flag set. Works only for TRILINOS_DEV Q4/2011 or newer.");
#endif
    }
    else
    {
      dserror( "unknown preconditioner for block matrix solver" );
    }
  }

#if 0
  preconditioner_->Print( std::cout );
  std::cout << "\n";
#endif
}
