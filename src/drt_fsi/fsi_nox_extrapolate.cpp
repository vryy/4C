
#ifdef TRILINOS_PACKAGE

#include "fsi_nox_extrapolate.H"

#include "NOX_Common.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_Group.H"
#include "NOX_Solver_Generic.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_GlobalData.H"

#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseSolver.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C"
{
#include "../headers/standardtypes.h"
}

using namespace NOX;
using namespace NOX::LineSearch;

NOX::FSI::Extrapolate::Extrapolate(const Teuchos::RefCountPtr<NOX::Utils>& utils,
                                   Teuchos::ParameterList& params)
  : utils_(utils)
{
}

NOX::FSI::Extrapolate::~Extrapolate()
{
}

bool NOX::FSI::Extrapolate::reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
                                  Teuchos::ParameterList& params)
{
  r_.clear();
  x_.clear();

  utils_ = gd->getUtils();
  //utils_->out() << RED_LIGHT "Extrapolate::reset" END_COLOR "\n";
  return true;
}

bool NOX::FSI::Extrapolate::compute(Abstract::Group& grp, double& step,
                                    const Abstract::Vector& dir,
                                    const Solver::Generic& s)
{
  //utils_->out() << YELLOW_LIGHT "Extrapolate::compute" END_COLOR "\n";
  const Abstract::Group& oldGrp = s.getPreviousSolutionGroup();

  r_.push_back(grp.getF().clone(DeepCopy));
  x_.push_back(grp.getX().clone(DeepCopy));
  x_.back()->update(1.0, dir, 1.0);

  //x_.back()->print(utils_->out());

  if (r_.size()==1)
  {
    grp.computeX(oldGrp, dir, 1.0);
  }
  else
  {
    Epetra_SerialDenseMatrix A(r_.size()-1,r_.size()-1);
    Epetra_SerialDenseVector rhs(r_.size()-1);
    Epetra_SerialDenseVector sol(r_.size()-1);

    Teuchos::RefCountPtr<NOX::Abstract::Vector> rk = r_.back();
    Teuchos::RefCountPtr<NOX::Abstract::Vector> xk = x_.back();
    std::vector<Teuchos::RefCountPtr<NOX::Abstract::Vector> > dr;
    std::vector<Teuchos::RefCountPtr<NOX::Abstract::Vector> > dx;
    for (unsigned i=0; i<r_.size()-1; ++i)
    {
      dr.push_back(rk->clone(ShapeCopy));
      dr.back()->update(1.0, *r_[i], -1.0, *rk, 0.0);

      dx.push_back(xk->clone(ShapeCopy));
      dx.back()->update(1.0, *x_[i], -1.0, *xk, 0.0);

      rhs[i] = - rk->innerProduct(*dr.back());
    }

    for (unsigned i=0; i<r_.size()-1; ++i)
      for (unsigned j=0; j<r_.size()-1; ++j)
      {
        A(i,j) = dr[i]->innerProduct(*dr[j]);
      }

    Epetra_SerialDenseSolver solver;
    solver.SetMatrix(A);
    solver.SetVectors(sol,rhs);

    int err = solver.Solve();
    if (err==0)
    {
      // We should check if the solution really helps.

      Teuchos::RefCountPtr<NOX::Abstract::Vector> d = rk->clone(ShapeCopy);
      d->init(0);
      for (unsigned i=0; i<r_.size()-1; ++i)
      {
        d->update(sol[i], *dx[i], 1.0);
      }

      // take the latest x and add the extrapolation to it
      grp.setX(*xk);
      grp.computeX(grp, *d, 1.0);
    }
    else
    {
      utils_->out() << RED_LIGHT "failed to solve linear system in Extrapolate::compute" END_COLOR "\n";

      // fallback
      grp.computeX(oldGrp, dir, 1.0);
    }
  }
  return true;
}

#endif
