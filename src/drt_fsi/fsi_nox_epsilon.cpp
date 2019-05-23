/*!----------------------------------------------------------------------

\brief Calculates vector extrapolation via epsilon algorithm.

\maintainer Matthias Mayr

\level 1
*/
/*---------------------------------------------------------------------*/

#include "fsi_nox_epsilon.H"

#include <NOX_GlobalData.H>
#include <NOX_Abstract_Group.H>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <NOX_Epetra_Group.H>
#include <NOX_Epetra_Vector.H>

#include <vector>

NOX::FSI::EpsilonExtrapolation::EpsilonExtrapolation(
    const Teuchos::RCP<NOX::Utils>& utils, Teuchos::ParameterList& params)
    : utils_(utils)
{
  Teuchos::ParameterList& mpeparams = params.sublist("Epsilon");
  kmax_ = mpeparams.get("kmax", 10);
  omega_ = mpeparams.get("omega", 0.1);
  // eps_ = mpeparams.get("Tolerance", 1e-8);
  maxcol_ = mpeparams.get("max col", 10);
}


NOX::FSI::EpsilonExtrapolation::~EpsilonExtrapolation() {}


bool NOX::FSI::EpsilonExtrapolation::reset(
    const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params)
{
  utils_ = gd->getUtils();
  return true;
}


bool NOX::FSI::EpsilonExtrapolation::compute(
    NOX::Abstract::Vector& dir, NOX::Abstract::Group& soln, const NOX::Solver::Generic& solver)
{
  // We work in a local copy of the group so that we do not spoil the
  // current state.
  NOX::Epetra::Group grp(dynamic_cast<NOX::Epetra::Group&>(soln));

  // Compute F at current solution
  NOX::Abstract::Group::ReturnType status = grp.computeF();
  if (status != NOX::Abstract::Group::Ok) throwError("compute", "Unable to compute F");

  // get f = d(k+1) - d(k)
  const NOX::Epetra::Vector& f = dynamic_cast<const NOX::Epetra::Vector&>(grp.getF());
  grp.computeX(grp, f, omega_);

  const NOX::Abstract::Vector& x = grp.getX();

  std::vector<Teuchos::RCP<NOX::Epetra::Vector>> epslist(maxcol_ + 1);
  epslist[0] = Teuchos::rcp(new NOX::Epetra::Vector(dynamic_cast<const NOX::Epetra::Vector&>(x)));

  Teuchos::RCP<NOX::Epetra::Vector> wg1;
  Teuchos::RCP<NOX::Epetra::Vector> wg2;
  Teuchos::RCP<NOX::Epetra::Vector> wg3;
  // Teuchos::RCP<NOX::Epetra::Vector> wg4;

  int indm = 1;
  for (int k = 0; k < kmax_; ++k)
  {
    // Compute F at current solution
    status = grp.computeF();
    if (status != NOX::Abstract::Group::Ok) throwError("compute", "Unable to compute F");

    // get f = d(k+1) - d(k)
    const NOX::Epetra::Vector& f = dynamic_cast<const NOX::Epetra::Vector&>(grp.getF());

    // We have to work on the scaled residual here.
    Teuchos::RCP<NOX::Epetra::Vector> y = Teuchos::rcp(new NOX::Epetra::Vector(f));
    y->scale(omega_);

    indm = k + 1;
    if (indm > maxcol_) indm = maxcol_;

    wg2 = Teuchos::rcp(new NOX::Epetra::Vector(*epslist[0]));
    wg2->update(1., *y, 1.);

    for (int i = 0; i < indm; ++i)
    {
      // epsilon extrapolation without care for instabilities
      wg1 = Teuchos::rcp(new NOX::Epetra::Vector(*wg2));
      wg1->update(-1., *epslist[i], 1.);

      double rd = wg1->norm();
      double rdd = epslist[i]->norm();
      if (rd < 1e-16 or rdd < 1e-16)
      {
        throwError("compute", "division by zero");
      }

      // vector inverse
      wg1->scale(1. / rd);

      if (i != 0)
      {
        wg1->update(1., *epslist[i - 1], 1.);
      }

      // shift the rhombus
      if (i != 0)
      {
        epslist[i - 1] = wg3;
      }
      wg3 = wg2;
      wg2 = wg1;
    }

    // save for the final values of the diagonal
    epslist[indm - 1] = wg3;
    epslist[indm] = wg2;

    // Update the group to go another round
    // Note: We do not use any extrapolated vector here but simply go
    // on with the series of vectors. The fixed relaxation is needed
    // to keep the iteration from diverging.
    grp.computeX(grp, f, omega_);
  }

  // set direction from original position
  if (indm % 2 == 0)
  {
    dir.update(1., *epslist[indm], -1., x, 0.);
  }
  else
  {
    dir.update(1., *epslist[indm - 1], -1., x, 0.);
  }

  return true;
}


bool NOX::FSI::EpsilonExtrapolation::compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& soln,
    const NOX::Solver::LineSearchBased& solver)
{
  return NOX::Direction::Generic::compute(dir, soln, solver);
}


void NOX::FSI::EpsilonExtrapolation::throwError(
    const std::string& functionName, const std::string& errorMsg)
{
  if (utils_->isPrintType(NOX::Utils::Error))
    utils_->err() << "EpsilonExtrapolation::" << functionName << " - " << errorMsg << std::endl;
  throw "NOX Error";
}
