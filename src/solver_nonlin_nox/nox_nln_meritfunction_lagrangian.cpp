/*-----------------------------------------------------------*/
/*!
\file nox_nln_meritfunction_lagrangian.cpp

\brief Implementation of the Lagrangian merit function for
       constrained problems.

\maintainer Michael Hiermeier

\date Jun 9, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_meritfunction_lagrangian.H"   // class header
#include "nox_nln_constraint_group.H"

#include <NOX_Abstract_Vector.H>
#include <NOX_Utils.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::MeritFunction::Lagrangian::Lagrangian(const Teuchos::RCP<NOX::Utils>& u) :
  meritFunctionName_("Lagrangian"),
  meritFunctionEnum_(mrtfct_lagrangian)
{
  utils_ = u;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computef(const NOX::Abstract::Group& grp) const
{
  double mrtFctVal = 0.0;

  if (!grp.isF())
  {
    throwError("computef()","The current function value was not computed yet. Please call "
        "computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if (constr_grp_ptr == NULL)
  {
    throwError("computef()","Dynamic cast to NOX::NLN::Constraint::Group failed!");
  }

  // Get the structural energy
  mrtFctVal = constr_grp_ptr->GetObjectiveModelValue("energy");
  // Get the part of the merit function, which is based on the constraint equations.
  // Get the constraint interfaces map
  const std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >&
      constr_interfaces = constr_grp_ptr->GetConstrInterfaces();
  std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >::const_iterator
      constr_iter;
  // loop over the second entries of the stl_map
  for (constr_iter=constr_interfaces.begin();constr_iter!=constr_interfaces.end();++constr_iter)
    mrtFctVal += constr_iter->second->GetConstrObjectiveModelValue(nameAsEnum());

  return mrtFctVal;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computeSlope(const NOX::Abstract::Vector& dir,
               const NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throwError("computeSlope()","The current function value was not computed yet. Please call "
      "computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if (constr_grp_ptr==NULL)
    throwError("computeSlope()","Dynamic cast to NOX::NLN::Constraint::Group failed!");

  // Get the constraint interfaces map
  const std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >&
      constr_interfaces = constr_grp_ptr->GetConstrInterfaces();
  std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >::const_iterator
      constr_iter;

  // compute the slope
  double slope = 0.0;
  for (constr_iter=constr_interfaces.begin();constr_iter!=constr_interfaces.end();++constr_iter)
  {
    // get the first order linearization terms of the Lagrangian objective model
    Teuchos::RCP<const std::vector<double> > firstOrderTerms =
        constr_iter->second->GetLinearizedObjectiveModelTerms(nameAsEnum(),
            linorder_first);

    for (std::size_t i=0;i<firstOrderTerms->size();++i)
      slope += firstOrderTerms->at(i);
  }
  return slope;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computeSaddlePointModel(
       const double& stepPV, const double& stepLM,
       const NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throwError("computeModel()","The current function value was not computed yet. Please call "
        "computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if (constr_grp_ptr==NULL)
  {
    throwError("computeModel()","Dynamic cast to NOX::NLN::Constraint::Group failed!");
  }

  double model = 0.0;
  // Get the constraint interfaces map
    const std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >&
        constr_interfaces = constr_grp_ptr->GetConstrInterfaces();
    std::map<NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> >::const_iterator
        constr_iter;
  for (constr_iter=constr_interfaces.begin();constr_iter!=constr_interfaces.end();++constr_iter)
  {
    // --------------------------------------------
    // Get the 1st order linearization terms of the Lagrangian objective model
    // w.r.t the primary degrees of freedom
    // --------------------------------------------
    Teuchos::RCP<const std::vector<double> > pFirstOrderTerms =
        constr_iter->second->GetLinearizedObjectiveModelTerms(nameAsEnum(),
        linorder_first,lin_wrt_primary_dofs);

    for (std::size_t i=0;i<pFirstOrderTerms->size();++i)
      model += stepPV * pFirstOrderTerms->at(i);

    // --------------------------------------------
    // Get the 1st order linearization terms of the Lagrangian objective model
    // w.r.t the Lagrange multiplier degrees of freedom
    // --------------------------------------------
    Teuchos::RCP<const std::vector<double> > lmFirstOrderTerms =
        constr_iter->second->GetLinearizedObjectiveModelTerms(nameAsEnum(),
        linorder_first,lin_wrt_lagrange_multiplier_dofs);

    for (std::size_t i=0;i<lmFirstOrderTerms->size();++i)
      model += stepLM * lmFirstOrderTerms->at(i);

    // --------------------------------------------
    // Get the 2nd order linearization terms of the Lagrangian objective model
    // w.r.t the Lagrange multiplier AND primary degrees of freedom
    // --------------------------------------------
    Teuchos::RCP<const std::vector<double> > pLmSecondOrderTerms =
        constr_iter->second->GetLinearizedObjectiveModelTerms(nameAsEnum(),
        linorder_second,lin_wrt_mixed_dofs);

    for (std::size_t i=0;i< pLmSecondOrderTerms->size();++i)
      model += stepLM * stepPV * pLmSecondOrderTerms->at(i);
  }

  return model;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computeSaddlePointModel(
       const double& step, const NOX::Abstract::Group& grp) const
{
  return computeSaddlePointModel(step,step,grp);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::string& NOX::NLN::MeritFunction::Lagrangian::name() const
{
  return meritFunctionName_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const NOX::NLN::MeritFunction::MeritFctName&
    NOX::NLN::MeritFunction::Lagrangian::nameAsEnum() const
{
  return meritFunctionEnum_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::MeritFunction::Lagrangian::throwError(
    const std::string& functionName,
    const std::string& errorMsg) const
{
  if (utils_->isPrintType(NOX::Utils::Error)) {
    utils_->out() << "ERROR - NOX::NLN::MeritFunction::Lagrangian::" << functionName
     << " - " << errorMsg << std::endl;
  }
  throw "NOX Error";
}
