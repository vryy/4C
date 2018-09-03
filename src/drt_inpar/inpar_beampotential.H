/*----------------------------------------------------------------------*/
/*!
\file inpar_beampotential.H

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
  /// Input parameters for potential-based beam interaction
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
      strategy_doublelengthspec_largesepapprox,  ///< double length specific potential, large
                                                 ///< separations
      strategy_doublelengthspec_smallsepapprox,  ///< double length specific potential, small
                                                 ///< separations
      strategy_singlelengthspec_smallsepapprox,  ///< single length specific potential, small
                                                 ///< separations
      strategy_vague
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
