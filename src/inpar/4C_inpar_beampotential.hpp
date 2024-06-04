/*----------------------------------------------------------------------*/
/*! \file

\brief input parameter definitions for beam potential-based interactions

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_BEAMPOTENTIAL_HPP
#define FOUR_C_INPAR_BEAMPOTENTIAL_HPP

#include "4C_config.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace BEAMPOTENTIAL
  {
    /// type of potential interaction
    /// (this enum represents the input file parameter BEAMPOTENTIAL_TYPE)
    enum BeamPotentialType
    {
      beampot_surf,  ///< surface potential
      beampot_vol,   ///< volume potential
      beampot_vague
    };

    /// available strategies/methods to evaluate potential interaction
    /// (this enum represents the input file parameter STRATEGY)
    enum BeamPotentialStrategy
    {
      strategy_doublelengthspec_largesepapprox,         ///< double length specific potential, large
                                                        ///< separations
      strategy_doublelengthspec_smallsepapprox,         ///< double length specific potential, small
                                                        ///< separations
      strategy_singlelengthspec_smallsepapprox,         ///< single length specific potential, small
                                                        ///< separations
      strategy_singlelengthspec_smallsepapprox_simple,  ///< reduced variant of the previous one
      strategy_vague
    };

    /// available types to regularize the force law for separations smaller than
    /// the specified regularization parameter
    enum BeamPotentialRegularizationType
    {
      regularization_linear,    ///< linear extrapolation
      regularization_constant,  ///< constant extrapolation, i.e. f(r)=f(r_reg) for all r<r_reg
      regularization_none       ///< no regularization
    };

    /**
     * \brief rule for how to assign the role of slave and master to beam elements
     */
    enum class MasterSlaveChoice
    {
      smaller_eleGID_is_slave,
      higher_eleGID_is_slave,
      choice_master_slave_vague
    };

    /// set the beam potential parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set beam potential specific conditions
    void SetValidConditions(
        std::vector<Teuchos::RCP<CORE::Conditions::ConditionDefinition>>& condlist);

  }  // namespace BEAMPOTENTIAL

}  // namespace INPAR

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
