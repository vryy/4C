/*-----------------------------------------------------------*/
/*!
\file nox_nln_meritfunction_infeasibility.cpp

\brief Implementation of the infeasibility merit function for
       constrained problems. Especially useful for the filter method.

\maintainer Michael Hiermeier

\date Apr 20, 2017

\level 3

*/
/*-----------------------------------------------------------*/


#include "nox_nln_meritfunction_infeasibility.H"
#include "nox_nln_constraint_group.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::MeritFunction::Infeasibility::Infeasibility(
    const Teuchos::ParameterList& params,
    const NOX::Utils& u )
    : /* utils_( u ), */
      infeasibility_type_( type_vague )
{
  const std::string& type_name = params.get<std::string>( "Type" );
  SetType( type_name );

  meritFunctionName_ = "Infeasibility measure (" + type_name +")";
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::MeritFunction::Infeasibility::SetType( const std::string& type_name )
{
  if ( type_name == "Two Norm" )
    infeasibility_type_ = type_two_norm;
  else
    dserror( "Unknown type name: \"%s\"", type_name.c_str() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Infeasibility::computef(
    const NOX::Abstract::Group& grp) const
{
  if ( not grp.isF() )
    dserror( "The current function value was not computed yet. "
        "Please call computeF() on the group passed into this function." );

  // cast the nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if ( not constr_grp_ptr )
    dserror( "Dynamic cast to NOX::NLN::Constraint::Group failed!" );

  return constr_grp_ptr->GetModelValue( Type() );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::MeritFunction::Infeasibility::computeGradient(
    const NOX::Abstract::Group& group,
    NOX::Abstract::Vector& result ) const
{
  dserror( "Currently unsupported." );
  exit( EXIT_FAILURE );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Infeasibility::computeSlope(
    const NOX::Abstract::Vector& dir,
    const NOX::Abstract::Group& grp) const
{
  if (!grp.isF())
  {
    dserror( "The current function value was not computed yet. Please call "
        "computeF() on the group passed into this function." );
  }

  // cast the underlying nox-group to the constraint group
  const NOX::NLN::CONSTRAINT::Group* constr_grp_ptr =
      dynamic_cast<const NOX::NLN::CONSTRAINT::Group*>(&grp);
  if ( not constr_grp_ptr )
    dserror( "Dynamic cast to NOX::NLN::Constraint::Group failed!" );

  // compute the slope
  return constr_grp_ptr->GetLinearizedModelTerms( dir, Type(),
      linorder_first, lin_wrt_all_dofs );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double NOX::NLN::MeritFunction::Infeasibility::computeQuadraticModel(
    const NOX::Abstract::Vector& dir,
    const NOX::Abstract::Group& grp) const
{
  dserror( "Currently unsupported." );
  exit( EXIT_FAILURE );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::MeritFunction::Infeasibility::computeQuadraticMinimizer(
    const NOX::Abstract::Group& grp,
    NOX::Abstract::Vector& result) const
{
  dserror( "Currently unsupported." );
  exit( EXIT_FAILURE );
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const std::string& NOX::NLN::MeritFunction::Infeasibility::name() const
{
  return meritFunctionName_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum NOX::NLN::MeritFunction::MeritFctName
NOX::NLN::MeritFunction::Infeasibility::Type() const
{
  switch ( infeasibility_type_ )
  {
    case type_two_norm:
      return mrtfct_infeasibility_two_norm;
    default:
    {
      dserror( "Unknown infeasibility type!" );
      exit( EXIT_FAILURE );
    }
  }

}
