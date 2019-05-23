/*----------------------------------------------------------------------*/
/*!

\brief Solve FSI problems using minimal polynomial vector extrapolation

\maintainer Matthias Mayr

\level 1

*/
/*----------------------------------------------------------------------*/

#include "fsi_nox_mpe.H"
#include <NOX_GlobalData.H>
#include <NOX_Abstract_Group.H>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <NOX_Epetra_Group.H>
#include <NOX_Epetra_Vector.H>

#include <vector>

#include <Epetra_Comm.h>
#include <Epetra_Time.h>

#include <Epetra_Vector.h>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"


NOX::FSI::MinimalPolynomial::MinimalPolynomial(
    const Teuchos::RCP<NOX::Utils>& utils, Teuchos::ParameterList& params)
    : utils_(utils)
{
  Teuchos::ParameterList& mpeparams = params.sublist("Extrapolation");
  kmax_ = mpeparams.get("kmax", 10);
  omega_ = mpeparams.get("omega", 0.01);
  eps_ = mpeparams.get("Tolerance", 1e-1);
  mpe_ = mpeparams.get("Method", "RRE") == "MPE";
}


NOX::FSI::MinimalPolynomial::~MinimalPolynomial() {}


bool NOX::FSI::MinimalPolynomial::reset(
    const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params)
{
  utils_ = gd->getUtils();
  return true;
}


bool NOX::FSI::MinimalPolynomial::compute(
    NOX::Abstract::Vector& dir, NOX::Abstract::Group& soln, const NOX::Solver::Generic& solver)
{
  // We work in a local copy of the group so that we do not spoil the
  // current state.
  NOX::Epetra::Group grp(dynamic_cast<NOX::Epetra::Group&>(soln));

  const NOX::Abstract::Vector& x = soln.getX();
  const NOX::Epetra::Vector& ex = dynamic_cast<const NOX::Epetra::Vector&>(x);

  std::vector<Teuchos::RCP<NOX::Epetra::Vector>> q;
  LINALG::SerialDenseMatrix r(kmax_ + 1, kmax_ + 1, true);
  LINALG::SerialDenseVector c(kmax_ + 1, true);
  LINALG::SerialDenseVector gamma(kmax_ + 1, true);

  int k;
  for (k = 0; k < kmax_; ++k)
  {
    Epetra_Time t(ex.getEpetraVector().Comm());
    NOX::Abstract::Group::ReturnType status;

    // Compute F at current solution
    status = grp.computeF();
    if (status != NOX::Abstract::Group::Ok) throwError("compute", "Unable to compute F");

    // get f = d(k+1) - d(k)
    const NOX::Epetra::Vector& f = dynamic_cast<const NOX::Epetra::Vector&>(grp.getF());

    // We have to work on the scaled residual here.
    Teuchos::RCP<NOX::Epetra::Vector> y = Teuchos::rcp(new NOX::Epetra::Vector(f));
    y->scale(omega_);

    // modified Gram-Schmidt
    for (int j = 0; j < k; ++j)
    {
      r(j, k) = y->innerProduct(*q[j]);
      y->update(-r(j, k), *q[j], 1.);
    }
    r(k, k) = sqrt(y->innerProduct(*y));

    // store new direction
    if (r(k, k) > 1e-32 * r(0, 0) and k < kmax_)
    {
      y->scale(1. / r(k, k));
      q.push_back(y);
    }
    else if (r(k, k) <= 1e-32 * r(0, 0))
    {
      if (utils_->isPrintType(NOX::Utils::Error))
        utils_->err() << "r(" << k << "," << k << ") <= " << 1e-32 << "*r(0,0)\n";
      break;
    }

    double res = 0;

    if (mpe_)
    {
      // MPE gamma calculation
      for (int i = k - 1; i >= 0; --i)
      {
        double ci = -r(i, k);
        for (int j = i + 1; j < k; ++j)
        {
          ci -= r(i, j) * c(j);
        }
        c(i) = ci / r(i, i);
      }
      c(k) = 1.;

      // sum over all entries of c
      double sc = 0.0;
      for (int l = 0; l < c.Length(); l++) sc += c(l);

      if (fabs(sc) < 1e-16)
      {
        throwError("compute", "sum(c) equals zero");
      }

      gamma.Update(1 / sc, c, 0.0);
      res = r(k, k) * fabs(gamma(k));

      if (utils_->isPrintType(NOX::Utils::InnerIteration))
      {
        utils_->out() << "MPE:  k=" << k << "  res=" << res << "  eps*r(0,0)=" << eps_ * r(0, 0)
                      << "  r(k,k)=" << r(k, k) << std::endl;
      }
    }
    else
    {
      // RRE gamma calculation
      for (int i = 0; i <= k; ++i)
      {
        double ci = 1.;
        for (int j = 0; j < i; ++j)
        {
          ci -= r(j, i) * c(j);
        }
        c(i) = ci / r(i, i);
      }
      for (int i = k; i >= 0; --i)
      {
        double ci = c(i);
        for (int j = i + 1; j <= k; ++j)
        {
          ci -= r(i, j) * gamma(j);
        }
        gamma(i) = ci / r(i, i);
      }

      // sum over all entries of gamma
      double sc = 0.0;
      for (int l = 0; l < gamma.Length(); l++) sc += gamma(l);

      gamma.Scale(1 / sc);
      res = 1. / sqrt(fabs(sc));

      if (utils_->isPrintType(NOX::Utils::InnerIteration))
      {
        utils_->out() << "RRE:  k=" << std::setw(2) << k << "  res=" << std::scientific << res
                      << "  eps*r(0,0)=" << std::scientific << eps_ * r(0, 0)
                      << "  r(k,k)=" << std::scientific << r(k, k) << "  time=" << std::scientific
                      << t.ElapsedTime() << std::endl;
      }
    }

    // leave if we are close enough to the solution
    if (res <= eps_ * r(0, 0) or r(k, k) <= 1e-32 * r(0, 0))
    {
      k += 1;
      break;
    }

    // Update the group to go another round
    // Note: We do not use any extrapolated vector here but simply go
    // on with the series of vectors. The fixed relaxation is needed
    // to keep the iteration from diverging.
    grp.computeX(grp, f, omega_);
  }

  // calc extrapolated vector
  LINALG::SerialDenseVector xi(kmax_, true);
  xi(0) = 1. - gamma(0);
  for (int j = 1; j < k; ++j)
  {
    xi(j) = xi(j - 1) - gamma(j);
  }

  NOX::Epetra::Vector s(dynamic_cast<const NOX::Epetra::Vector&>(x));
  for (int j = 0; j < k; ++j)
  {
    double hp = 0.;
    for (int i = j; i < k; ++i)
    {
      hp += r(j, i) * xi(i);
    }
    s.update(hp, *q[j], 1.);
  }

  // set direction from original position
  dir.update(1., s, -1., x, 0.);

  return true;
}


bool NOX::FSI::MinimalPolynomial::compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& soln,
    const NOX::Solver::LineSearchBased& solver)
{
  return NOX::Direction::Generic::compute(dir, soln, solver);
}


void NOX::FSI::MinimalPolynomial::throwError(
    const std::string& functionName, const std::string& errorMsg)
{
  if (utils_->isPrintType(NOX::Utils::Error))
    utils_->err() << "MinimalPolynomial::" << functionName << " - " << errorMsg << std::endl;
  throw "NOX Error";
}
