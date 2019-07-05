/*-----------------------------------------------------------*/
/*! \file
\file inpar_fbi.H

\brief input parameter for Fluid-Beam Interaction

\maintainer Nora Hagmeyer

\level 3

*/
/*-----------------------------------------------------------*/
#ifndef INPAR_FBI_H
#define INPAR_FBI_H

#include "inpar_parameterlist_utils.H"
#include <boost/algorithm/string/predicate.hpp>

#include "../drt_lib/drt_dserror.H"

#include "../drt_fem_general/drt_utils_integration.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace FBI
  {
    /// Coupling of the Fluid and the beam problems
    enum class BeamToFluidCoupling
    {
      fluid,  //< Coupling on the fluid partition, while the beam is not influenced
      solid,  //< Coupling on the structure partition, while the fluid is not influenced
      twoway  //< Full two-way FBI coupling
    };
    /// Constraint enforcement for beam to fluid meshtying.
    enum class BeamToFluidConstraintEnforcement
    {
      //! Default value.
      none,
      //! Penalty method.
      penalty
    };

    /// Discretization approach for beam to fluid meshtying.
    enum class BeamToFluidDiscretization
    {
      none,                    //< Default value
      gauss_point_to_segment,  //< Gauss point to segment approach
      mortar                   //< mortar-type segment to segment approach
    };

    /**
     * \brief Shape function for the mortar Lagrange-multiplicators
     */
    enum class BeamToFluidMeshtingMortarShapefunctions
    {
      //! Default value.
      none,
      //! Linear.
      line2,
      //! Quadratic.
      line3,
      //! Cubic.
      line4
    };

    /// set the beam interaction parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set beam interaction specific conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace FBI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif
