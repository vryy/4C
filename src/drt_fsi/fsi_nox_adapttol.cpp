#ifdef CCADISCRET

#include <string>

#include "fsi_nox_adapttol.H"

#include <NOX_Solver_Generic.H>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
NOX::FSI::AdaptiveTolerance::AdaptiveTolerance(Teuchos::RCP<NOX::Utils> utils,
                                               Teuchos::ParameterList& lsParams)
  : utils_(utils),
    lsParams_(lsParams),
    desirednlnres_(0.),
    currentnlnres_(0.),
    plaintol_(1e-4),
    better_(0.1)
{
  lsParams_.set<std::string>("Convergence Test","r0");
  lsParams_.set<double>("Tolerance",plaintol_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void NOX::FSI::AdaptiveTolerance::runPreIterate(const NOX::Solver::Generic& solver)
{
  // do adaptive linear solver tolerance (not in first solve)
  int itnum = solver.getNumIterations();
  if (itnum>0 and currentnlnres_!=0)
  {
    if (currentnlnres_*plaintol_ < desirednlnres_)
    {
      double tol = desirednlnres_*better_/currentnlnres_;
      if (tol>plaintol_)
      {
        lsParams_.set<double>("Tolerance",tol);
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void NOX::FSI::AdaptiveTolerance::runPostIterate(const NOX::Solver::Generic& solver)
{
#if 0
  NOX::StatusTest::StatusType status = solver.getStatus();
  if (status==NOX::StatusTest::Converged)
  {
    currentnlnres_ = 0.;
    return;
  }
#endif
}


#endif
