/*-----------------------------------------------------------*/
/*! \file
\file inpar_fbi.H

\brief input parameter for Fluid-Beam Interaction


\level 3

*/
/*-----------------------------------------------------------*/
#ifndef FOUR_C_INPAR_FBI_HPP
#define FOUR_C_INPAR_FBI_HPP

#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_integration.hpp"
#include "baci_utils_exceptions.hpp"
#include "baci_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace INPUT
{
  class ConditionDefinition;
}

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

    /// Parallel presorting strategy to be used for the beam mesh
    enum class BeamToFluidPreSortStrategy
    {
      bruteforce,  //< each processor searches for each beam if it is near one of its fluid elements
      binning  //< each processor only searches for beam elements which lie in or around its bins
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
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

  }  // namespace FBI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
