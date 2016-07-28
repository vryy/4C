/*----------------------------------------------------------------------*/
/*!
\file inpar_cavitation.H

\brief Input parameters for cavitation problems

\level 3

\maintainer Georg Hammerl
*-----------------------------------------------------------------------*/
#ifndef INPAR_CAVITATION_H
#define INPAR_CAVITATION_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

/*----------------------------------------------------------------------*/

namespace INPAR
{
  /// Input parameters for cavitation
  namespace CAVITATION
  {

    /// Type of coupling strategy for cavitation problems
    enum CouplingStrategyOverFields
    {
      OneWay,
      TwoWayMomentum,
      VoidFracOnly,
      TwoWayFull_weak,
      TwoWayFull_strong
    };

    /// Strategy for calculation of void fraction of fluid elements
    enum VoidFractionCalculation
    {
      analytical_constpoly,
      analytical_quadraticpoly,
      analytical_quarticpoly,
      gaussian_integration
    };

    //! reconstruction type of fluid fraction
    enum FluidFracReconstructionMethod
    {
      fluidfracreco_spr,
      fluidfracreco_l2
    };

    /// set the cavitation parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace CAVITATION

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CAVITATION_H
