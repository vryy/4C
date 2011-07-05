/*
 * solver_blockpreconditioners.cpp
 *
 *  Created on: Jul 4, 2011
 *      Author: wiesner
 */

#include "../drt_lib/drt_dserror.H"

#include "solver_blockpreconditioners.H"


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::SimplePreconditioner::SimplePreconditioner( FILE * outfile,
                                                            Teuchos::ParameterList & params,
                                                            Teuchos::ParameterList & simpleparams )
  : LINALG::SOLVER::PreconditionerType( outfile ),
    params_( params ),
    simpleparams_( simpleparams )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::SimplePreconditioner::Setup( bool create,
                                                  Epetra_Operator * matrix,
                                                  Epetra_MultiVector * x,
                                                  Epetra_MultiVector * b )
{
  SetupLinearProblem( matrix, x, b );

  if ( create )
  {
    // SIMPLER does not need copy of preconditioning matrix to live
    // SIMPLER does not use the downwinding installed here, it does
    // its own downwinding inside if desired

    // free old matrix first
    P_ = Teuchos::null;

    // temporary hack: distinguish between "old" SIMPLER_Operator (for fluid
    // only) and "new" more general test implementation
    bool mt = simpleparams_.get<bool>("MESHTYING",false);
    bool co = simpleparams_.get<bool>("CONTACT",false);
    bool cstr = simpleparams_.get<bool>("CONSTRAINT",false);
    if (mt || co || cstr)
    {
      P_ = Teuchos::rcp(new LINALG::SIMPLER_BlockPreconditioner(Teuchos::rcp( matrix, false ),params_,simpleparams_,outfile_));
    }
    else
    {
      P_ = Teuchos::rcp(new LINALG::SIMPLER_Operator(Teuchos::rcp( matrix, false ),params_,simpleparams_,outfile_));
    }
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::AMGBSPreconditioner::AMGBSPreconditioner( FILE * outfile,
    Teuchos::ParameterList & params )
  : LINALG::SOLVER::PreconditionerType( outfile ),
    params_( params )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::AMGBSPreconditioner::Setup( bool create,
                                                 Epetra_Operator * matrix,
                                                 Epetra_MultiVector * x,
                                                 Epetra_MultiVector * b )
{
  SetupLinearProblem( matrix, x, b );

  if ( create )
  {
    // free old matrix first
    P_ = Teuchos::null;

    // Params().sublist("AMGBS") just contains the Fluid Pressure Solver block
    // from the dat file (not needed anymore)
    P_ = Teuchos::rcp(new LINALG::SaddlePointPreconditioner(Teuchos::rcp( matrix, false ),params_,outfile_));
  }
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::BGSPreconditioner::BGSPreconditioner( FILE * outfile,
    Teuchos::ParameterList & params, Teuchos::ParameterList & bgslist )
  : LINALG::SOLVER::PreconditionerType( outfile ),
    params_( params ),
    bgslist_( bgslist )
{
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::BGSPreconditioner::Setup( bool create,
                                               Epetra_Operator * matrix,
                                               Epetra_MultiVector * x,
                                               Epetra_MultiVector * b )
{
  SetupLinearProblem( matrix, x, b );

  if ( create )
  {
    P_ = Teuchos::null;

    int numblocks = bgslist_.get<int>("numblocks");

    if (numblocks == 2) // BGS2x2
    {
      // check whether sublists for individual block solvers are present
      bool haveprec1 = params_.isSublist("PREC1");
      bool haveprec2 = params_.isSublist("PREC2");
      if (!haveprec1 or !haveprec2)
        dserror("individual block solvers for BGS2x2 need to be specified");

      int global_iter = bgslist_.get<int>("global_iter");
      double global_omega = bgslist_.get<double>("global_omega");
      int block1_iter = bgslist_.get<int>("block1_iter");
      double block1_omega = bgslist_.get<double>("block1_omega");
      int block2_iter = bgslist_.get<int>("block2_iter");
      double block2_omega = bgslist_.get<double>("block2_omega");
      bool fliporder = bgslist_.get<bool>("fliporder");

      P_ = Teuchos::rcp(new LINALG::BGS2x2_Operator(Teuchos::rcp( matrix, false ),
                                           params_.sublist("PREC1"),
                                           params_.sublist("PREC2"),
                                           global_iter,
                                           global_omega,
                                           block1_iter,
                                           block1_omega,
                                           block2_iter,
                                           block2_omega,
                                           fliporder,
                                           outfile_));
    }
    else
      dserror("Block Gauss-Seidel is currently only implemented for a 2x2 system");
  }
}
