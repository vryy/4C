
#ifdef TRILINOS_PACKAGE
#ifdef D_FSI

#include "fsi_aitken_nox.H"
#include <NOX.H>

extern "C"
{
#include "../headers/standardtypes.h"
}

#include "../discret/dstrc.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Aitken::reset(const Teuchos::RefCountPtr< NOX::GlobalData > &gd, Teuchos::ParameterList &params)
{
  nu_ = 0;
  if (!is_null(del_))
  {
    del_->init(1e20);
  }
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Aitken::compute(NOX::Abstract::Vector &dir,
                     NOX::Abstract::Group &soln,
                     const NOX::Solver::Generic &solver)
{
  DSTraceHelper("Aitken::compute");
//   NOX::Abstract::Group::ReturnType status;

  // Compute F at current solution.
//   status = soln.computeF();
//   if (status != NOX::Abstract::Group::Ok)
//     dserror("Unable to compute F");

  //const NOX::Abstract::Vector& x = soln.getX();
  const NOX::Abstract::Vector& F = soln.getF();

  if (is_null(del_))
  {
    del_  = F.clone();
    del2_ = F.clone();
    del_->init(1e20);
  }

  del2_->update(1,*del_,1,F);
  del_ ->update(-1,F);

  double top = del2_->innerProduct(*del_);
  double den = del2_->innerProduct(*del2_);

  nu_ = nu_ + (nu_ - 1.)*top/den;
  double omega = 1. - nu_;

  std::cout << "omega = " << omega << "\n";

  dir = F;
  dir.scale(omega);

  return true;
}

#endif
#endif
