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

#include "../linalg/linalg_serialdensevector.H"

#include <NOX_Abstract_Vector.H>
#include <NOX_Utils.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::MeritFunction::Lagrangian::Lagrangian(const Teuchos::RCP<NOX::Utils>& u)
    : meritFunctionName_("Lagrangian")
{
  utils_ = u;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computef(
    const NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throwError("computef()","The current function value was not computed yet. "
        "Please call computeF() on the group passed into this function.");
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if ( not constr_grp_ptr )
  {
    throwError("computef()","Dynamic cast to NOX::NLN::Constraint::Group failed!");
  }

  // Get the primary contribution and constraint contributions
  return constr_grp_ptr->GetModelValue( Type() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computeSlope(
    const NOX::Abstract::Vector& dir,
    const NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throwError( "computeSlope()", "The current function value was not computed yet. "
        "Please call computeF() on the group passed into this function." );
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if ( not constr_grp_ptr )
    throwError( "computeSlope()", "Dynamic cast to NOX::NLN::Constraint::Group failed!" );

  // compute the slope
  return constr_grp_ptr->GetLinearizedModelTerms( dir, Type(),
      linorder_first, lin_wrt_all_dofs );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computeMixed2ndOrderTerms(
    const NOX::Abstract::Vector& dir,
    const NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throwError( "computeSlope()", "The current function value was not computed yet. "
        "Please call computeF() on the group passed into this function." );
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if ( not constr_grp_ptr )
    throwError( "computeSlope()", "Dynamic cast to NOX::NLN::Constraint::Group failed!" );

  // compute the slope
  return constr_grp_ptr->GetLinearizedModelTerms( dir, Type(),
      linorder_second, lin_wrt_mixed_dofs );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computeSaddlePointModel(
    const double& stepPV,
    const double& stepLM,
    const NOX::Abstract::Vector& dir,
    const NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    throwError( "computeSaddlePointModel()", "The current function value was not "
        "computed yet. Please call computeF() on the group passed into this function." );
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if ( not constr_grp_ptr )
  {
    throwError("computeModel()","Dynamic cast to NOX::NLN::Constraint::Group failed!");
  }

  const NOX::NLN::CONSTRAINT::Group& constr_grp = *constr_grp_ptr;

  // compute the function value
  double model = 0.0;

  // --------------------------------------------
  // Get the 1-st order linearization terms of the Lagrangian objective model
  // w.r.t the primary degrees of freedom
  // --------------------------------------------
  model += stepPV * constr_grp.GetLinearizedModelTerms( dir, Type(),
      linorder_first, lin_wrt_primary_dofs );

  // --------------------------------------------
  // Get the 1-st order linearization terms of the Lagrangian objective model
  // w.r.t the Lagrange multiplier degrees of freedom
  // --------------------------------------------
  model += stepLM * constr_grp.GetLinearizedModelTerms( dir, Type(),
      linorder_first, lin_wrt_lagrange_multiplier_dofs );

  // --------------------------------------------
  // Get the 2-nd order linearization terms of the Lagrangian objective model
  // w.r.t the Lagrange multiplier AND primary degrees of freedom
  // --------------------------------------------
  model += stepLM * stepPV * constr_grp.GetLinearizedModelTerms( dir, Type(),
      linorder_second, lin_wrt_mixed_dofs );

  return model;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Lagrangian::computeSaddlePointModel(
    const double& step,
    const NOX::Abstract::Vector& dir,
    const NOX::Abstract::Group& grp) const
{
  return computeSaddlePointModel( step, step, dir, grp );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::string& NOX::NLN::MeritFunction::Lagrangian::name() const
{
  return meritFunctionName_;
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
