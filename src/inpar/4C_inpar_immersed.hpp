#ifndef FOUR_C_INPAR_IMMERSED_HPP
#define FOUR_C_INPAR_IMMERSED_HPP

#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace Core::Conditions
{
  class ConditionDefinition;
}
namespace Inpar
{
  namespace Immersed
  {
    /*----------------------------------------------------------------------*
     | Coupling Methods                                                     |
     *----------------------------------------------------------------------*/
    enum ImmersedCoupling
    {
      partitioned,
      monolithic
    };

    typedef enum PartitionedScheme
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
    void set_valid_parameters(Teuchos::ParameterList& list);

    /// set specific immersed conditions
    void set_valid_conditions(
        std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist);

  }  // namespace Immersed
}  // namespace Inpar
/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
