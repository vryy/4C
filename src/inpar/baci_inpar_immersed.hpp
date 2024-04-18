/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for immersed

\level 1


*/

/*----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_IMMERSED_HPP
#define FOUR_C_INPAR_IMMERSED_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace INPUT
{
  class ConditionDefinition;
}

namespace INPAR
{
  namespace IMMERSED
  {
    /*----------------------------------------------------------------------*
     | Coupling Methods                                                     |
     *----------------------------------------------------------------------*/
    enum ImmersedCoupling
    {
      partitioned,
      monolithic
    };

    typedef enum ParitionedScheme
    {
      cell_coupling_undefined = 0,
      cell_basic_sequ_stagg = 1,
      cell_iter_stagg_fixed_rel_param = 2,
      cell_iter_stagg_AITKEN_rel_param = 3,
    } PARITIONED_SCHEME;

    enum ImmersedCouplingScheme
    {
      neumannneumann,
      dirichletneumann
    };

    enum ImmersedProjection
    {
      shapefunctions,
      mortar
    };

    enum ImmersedRelaxation
    {
      globally,
      selectively
    };

    enum ImmersedNlnsolver
    {
      nlnsolver_stop,
      nlnsolver_continue
    };

    enum ImmersedRelaxationparam
    {
      fixed,
      aitken
    };


    /// set the immersed parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific immersed conditions
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

  }  // namespace IMMERSED
}  // namespace INPAR
/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
