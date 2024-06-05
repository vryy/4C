/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for biomedical simulations

\level 3


*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_INPAR_BIO_HPP
#define FOUR_C_INPAR_BIO_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
/*----------------------------------------------------------------------*/

// ToDo: move these enums to namespace INPAR::ARTNET etc.
//       is the typedef really needed?

/*!----------------------------------------------------------------------
\brief enum of arterial network dynamic types

\author ismail \date 11/08
This is the enumeration of all types of different integration schemes

*-----------------------------------------------------------------------*/


/*!----------------------------------------------------------------------
\brief enum of reduced dimensional airways dynamic types

\author ismail \date 01/10
This is the enumeration of all types of different integration schemes

*-----------------------------------------------------------------------*/
typedef enum RedAirwaysDyntype
{
  one_step_theta,
  linear,
  nonlinear
} _RED_AIRWAYS_DYNTYPE;



namespace INPAR
{
  namespace ARTDYN
  {
    enum TimeIntegrationScheme
    {
      tay_gal,
      stationary
    };

    /// initial field for artery problem
    enum InitialField
    {
      initfield_zero_field,
      initfield_field_by_function,
      initfield_field_by_condition
    };

    //! element implementation type
    enum ImplType
    {
      impltype_undefined,
      impltype_lin_exp,
      impltype_pressure_based
    };

    /// set the arterial dynamic parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);
  }  // namespace ARTDYN

  namespace ARTNET
  {
    /*!----------------------------------------------------------------------
    \brief enum of reduced dimensional relaxation type

    \author roth \date 10/13
    This is the enumeration of all types of different relaxation types

    *-----------------------------------------------------------------------*/
    enum Relaxtype3D0D
    {
      norelaxation,
      fixedrelaxation,
      Aitken,
      SD
    };

    //! type of coupling between artery network and poromultiphasescatra-framework
    enum ArteryPoroMultiphaseScatraCouplingMethod
    {
      none,   // none
      nodal,  // nodal
      gpts,   // Gauss-Point-To-Segment
      mp,     // Mortar-Penalty
      ntp     // 1Dnode-to-point in 2D/3D
    };

    /// set the artnet parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific artnet conditions
    void SetValidConditions(
        std::vector<Teuchos::RCP<CORE::Conditions::ConditionDefinition>>& condlist);

  }  // namespace ARTNET

  namespace BIOFILM
  {
    /// set the biofilm parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific biofilm conditions
    void SetValidConditions(
        std::vector<Teuchos::RCP<CORE::Conditions::ConditionDefinition>>& condlist);
  }  // namespace BIOFILM

  namespace REDAIRWAYS
  {
    /// set the reduced airways parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific reduced airways conditions
    void SetValidConditions(
        std::vector<Teuchos::RCP<CORE::Conditions::ConditionDefinition>>& condlist);
  }  // namespace REDAIRWAYS
}  // namespace INPAR

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
