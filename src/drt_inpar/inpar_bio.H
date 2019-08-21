/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for biomedical simulations

\level 3

\maintainer Johannes Kremheller

*/
/*----------------------------------------------------------------------*/
#ifndef INPAR_ARTNET_H
#define INPAR_ARTNET_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

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
typedef enum _RED_AIRWAYS_DYNTYPE
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
    enum RELAXTYPE_3D_0D
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
      mp      // Mortar-Penalty
    };

    /// set the artnet parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific artnet conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace ARTNET

  namespace BIOFILM
  {
    /// set the biofilm parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific biofilm conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);
  }  // namespace BIOFILM

  namespace PATSPEC
  {
    /// set the patient specific parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific patspec conditions (AAA)
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);
  }  // namespace PATSPEC

  namespace REDAIRWAYS
  {
    /// set the reduced airways parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific reduced airways conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);
  }  // namespace REDAIRWAYS
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_ARTNET_H
