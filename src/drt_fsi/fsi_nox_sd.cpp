/*----------------------------------------------------------------------*/
/*!
\file fsi_nox_sd.cpp

\brief Solve FSI problems using NonlinearCG

\maintainer Andreas Rauch

\level 1

*/
/*----------------------------------------------------------------------*/

#include "fsi_nox_sd.H"

#include <NOX_Epetra_Interface_Required.H>
#include <NOX_Common.H>
#include <NOX_Abstract_Vector.H>
#include <NOX_Abstract_Group.H>
#include <NOX_Epetra_Vector.H>
#include <NOX_Epetra_Group.H>
#include <NOX_Solver_Generic.H>
#include <Teuchos_ParameterList.hpp>
#include <NOX_Utils.H>
#include <NOX_GlobalData.H>

#include <Epetra_Vector.h>
#include <Epetra_Comm.h>
#include <NOX_Epetra_Vector.H>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"


NOX::FSI::SDRelaxation::SDRelaxation(const Teuchos::RCP<NOX::Utils>& utils,
                                     Teuchos::ParameterList& params)
  : utils_(utils)
{
}


NOX::FSI::SDRelaxation::~SDRelaxation()
{
}


bool NOX::FSI::SDRelaxation::reset(const Teuchos::RCP<NOX::GlobalData>& gd,
                                   Teuchos::ParameterList& params)
{
  utils_ = gd->getUtils();
  //Teuchos::ParameterList& p = params.sublist("SDRelaxation");
  return true;
}


bool NOX::FSI::SDRelaxation::compute(NOX::Abstract::Group& newgrp,
                                     double& step,
                                     const NOX::Abstract::Vector& dir,
                                     const NOX::Solver::Generic& s)
{
  if (utils_->isPrintType(NOX::Utils::InnerIteration))
  {
    utils_->out() << "\n" << NOX::Utils::fill(72) << "\n"
                 << "-- SDRelaxation Line Search -- \n";
  }

  const NOX::Abstract::Group& oldgrp = s.getPreviousSolutionGroup();
  NOX::Epetra::Group& egrp = dynamic_cast<NOX::Epetra::Group&>(newgrp);

  // Perform single-step linesearch

  // Note that the following could be wrapped with a while loop to allow
  // iterations to be attempted

  double numerator = oldgrp.getF().innerProduct(dir);
  double denominator = computeDirectionalDerivative(dir, *egrp.getRequiredInterface())
                       .innerProduct(dir);

  step = - numerator / denominator;
  utils_->out() << "          RELAX = " << std::setw(5) << step << "\n";

  newgrp.computeX(oldgrp, dir, step);
  newgrp.computeF();

  double checkOrthogonality = fabs( newgrp.getF().innerProduct(dir) );

  if (utils_->isPrintType(NOX::Utils::InnerIteration)) {
    utils_->out() << std::setw(3) << "1" << ":";
    utils_->out() << " step = " << utils_->sciformat(step);
    utils_->out() << " orth = " << utils_->sciformat(checkOrthogonality);
    utils_->out() << "\n" << NOX::Utils::fill(72) << "\n" << std::endl;
  }

  // write omega
  double fnorm = oldgrp.getF().norm();
  if (dynamic_cast<const NOX::Epetra::Vector&>(oldgrp.getF()).getEpetraVector().Comm().MyPID()==0)
  {
    static int count;
    static std::ofstream* out;
    if (out==NULL)
    {
      std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
      s.append(".omega");
      out = new std::ofstream(s.c_str());
    }
    (*out) << count << " "
           << step << " "
           << fnorm
           << "\n";
    count += 1;
    out->flush();
  }

  return true;
}


NOX::Abstract::Vector&
NOX::FSI::SDRelaxation::computeDirectionalDerivative(const NOX::Abstract::Vector& dir,
                                                     NOX::Epetra::Interface::Required& interface)
{
  // Allocate space for vecPtr and grpPtr if necessary
  if (Teuchos::is_null(vecPtr_))
    vecPtr_ = dir.clone(NOX::ShapeCopy);

  const NOX::Epetra::Vector& edir = dynamic_cast<const NOX::Epetra::Vector&>(dir);
  NOX::Epetra::Vector& evec = dynamic_cast<NOX::Epetra::Vector&>(*vecPtr_);

  // we do not want the group to remember this solution
  // and we want to set our own flag
  // this tells computeF to do a SD relaxation calculation
  interface.computeF(edir.getEpetraVector(),
                     evec.getEpetraVector(),
                     NOX::Epetra::Interface::Required::User);

  return *vecPtr_;
}

