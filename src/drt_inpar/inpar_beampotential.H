/*----------------------------------------------------------------------*/
/*! \file

\brief input parameter definitions for beam potential-based interactions

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------*/

#ifndef INPAR_BEAMPOTENTIAL_H
#define INPAR_BEAMPOTENTIAL_H

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

    /// rule for how to assign the role of slave and master to beam elements
    enum MasterSlaveChoice
    {
      smaller_eleGID_is_slave,
      higher_eleGID_is_slave,
      choice_master_slave_vague
    };

    /// data format for written numeric data via vtk
    enum OutputDataFormat
    {
      binary,
      ascii,
      vague
    };

    /// set the beam potential parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set beam potential specific conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace BEAMPOTENTIAL

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_BEAMPOTENTIAL_H
