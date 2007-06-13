
#ifdef TRILINOS_PACKAGE

#include "fsi_nox_michler.H"

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
using namespace NOX::Direction;

Michler::Michler(const Teuchos::RefCountPtr<NOX::Utils>& utils,
                 Teuchos::ParameterList& params)
  : utils_(utils), params_(params)
{
}

Michler::~Michler()
{
}

bool Michler::reset(const Teuchos::RefCountPtr<NOX::GlobalData>& gd,
                    Teuchos::ParameterList& params)
{
  r_.clear();
  x_.clear();

  utils_ = gd->getUtils();
  //utils_->out() << RED_LIGHT "Michler::reset" END_COLOR "\n";

  nu_ = params_.sublist("Michler").get("nu", 1e-2);

  return true;
}


bool Michler::compute(NOX::Abstract::Vector& dir,
                      NOX::Abstract::Group& soln,
                      const NOX::Solver::Generic& solver)
{
  NOX::Abstract::Group::ReturnType status;

  Teuchos::RefCountPtr<NOX::Abstract::Vector> x0 = soln.getX().clone();

  // Compute F at current solution
  status = soln.computeF();
  if (status != NOX::Abstract::Group::Ok)
    throwError("compute", "Unable to compute F");

  // So we got the residuum we start with.
  Teuchos::RefCountPtr<NOX::Abstract::Vector> ri = soln.getF().clone();
  double eta = ri->norm();
  double epsilon = 0.1*eta;
  double nu = nu_*eta;

  std::vector<Teuchos::RefCountPtr<NOX::Abstract::Vector> > dr;
  std::vector<Teuchos::RefCountPtr<NOX::Abstract::Vector> > dx;

  // save current residuum
  r_.push_back(soln.getF().clone());

  // save current interface displacement
  x_.push_back(soln.getX().clone());
  x_.back()->update(1.0, *r_.back(), 1.0);

  while (eta > epsilon)
  {
    Teuchos::RefCountPtr<NOX::Abstract::Vector> rj = r_.back();
    Teuchos::RefCountPtr<NOX::Abstract::Vector> xj = x_.back();

    // calculate interface displacement difference
    dx.push_back(xj->clone(ShapeCopy));
    Teuchos::RefCountPtr<NOX::Abstract::Vector> dxj = dx.back();
    dxj->update(1.0, *xj, -1.0, *x0);

    // calculate orthonormal (relaxed) basis of interface displacements
    for (unsigned k=0; k<dx.size()-1; ++k)
    {
      Teuchos::RefCountPtr<NOX::Abstract::Vector> dxk = dx[k];
      dxj->update(-dxj->innerProduct(*dxk)/dxk->innerProduct(*dxk), *dxk, 1.0);
    }

    // That's bad!
    dxj->scale(nu/dxj->norm());
    //dxj->scale(nu_);

    // modify current displacement
    xj->update(1.0, *x0, 1.0, *dxj);

    // compute new residual
    soln.setX(*xj);
    status = soln.computeF();
    if (status != NOX::Abstract::Group::Ok)
      throwError("compute", "Unable to compute F");

    // save current residuum
    r_.push_back(soln.getF().clone());
    Teuchos::RefCountPtr<NOX::Abstract::Vector> rj1 = r_.back();

    // save current interface displacement
    x_.push_back(soln.getX().clone());
    Teuchos::RefCountPtr<NOX::Abstract::Vector> xj1 = x_.back();
    xj1->update(1.0, *rj1, 1.0);

    // calculate residuum difference
    dr.push_back(rj1->clone(ShapeCopy));
    Teuchos::RefCountPtr<NOX::Abstract::Vector> drj = dr.back();
    drj->update(1.0, *rj1, -1.0, *ri);

    // now find the minimum
    Epetra_SerialDenseMatrix A(dr.size(),dr.size());
    Epetra_SerialDenseVector rhs(dr.size());
    Epetra_SerialDenseVector sol(dr.size());

    for (unsigned i=0; i<dr.size(); ++i)
    {
      rhs[i] = - ri->innerProduct(*dr[i]);
    }

    for (unsigned i=0; i<dr.size(); ++i)
      for (unsigned j=0; j<dr.size(); ++j)
      {
        A(i,j) = dr[i]->innerProduct(*dr[j]);
      }

    Epetra_SerialDenseSolver solver;
    solver.SetMatrix(A);
    solver.SetVectors(sol,rhs);

    int err = solver.Solve();
    if (err==0)
    {
      Teuchos::RefCountPtr<NOX::Abstract::Vector> r = ri->clone();
      for (unsigned i=0; i<dr.size(); ++i)
      {
        r->update(sol[i], *dr[i], 1.0);
      }
      eta = r->norm();

      // in case we are done: calculate new direction
      dir.init(0);
      for (unsigned i=0; i<dx.size(); ++i)
      {
        dir.update(sol[i], *dx[i], 1.0);
      }
      dir.update(1.0, soln.getF());
    }
    else
    {
      dserror("failed to solve linear system: err=%d",err);
    }
  }

  // reset interface displacement
  soln.setX(*x0);

  return true;
}

bool Michler::compute(NOX::Abstract::Vector& dir,
                      NOX::Abstract::Group& soln,
                      const NOX::Solver::LineSearchBased& solver)
{
  return NOX::Direction::Generic::compute( dir, soln, solver );
}


void Michler::throwError(const string& functionName, const string& errorMsg)
{
    if (utils_->isPrintType(NOX::Utils::Error))
      utils_->err() << "Michler::" << functionName
                    << " - " << errorMsg << endl;
    throw "NOX Error";
}

////////////////////////////////////////////////////////////

#if 0

bool Michler::compute(Abstract::Group& grp, double& step,
                          const Abstract::Vector& dir,
                          const Solver::Generic& s)
{
  //utils_->out() << YELLOW_LIGHT "Michler::compute" END_COLOR "\n";
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
      utils_->out() << RED_LIGHT "failed to solve linear system in Michler::compute" END_COLOR "\n";

      // fallback
      grp.computeX(oldGrp, dir, 1.0);
    }
  }
  return true;
}

#endif

#endif
