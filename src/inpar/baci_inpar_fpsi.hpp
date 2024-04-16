/*----------------------------------------------------------------------*/
/*! \file
\brief fpsi parameters
\level 1
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_INPAR_FPSI_HPP
#define FOUR_C_INPAR_FPSI_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

// forward declaration
namespace INPUT
{
  class ConditionDefinition;
}

/*----------------------------------------------------------------------*
 | Coupling Methods                                                                  |
 *----------------------------------------------------------------------*/
enum _FPSI_COUPLING
{
  fpsi_monolithic_plain,
  partitioned
};

namespace INPAR
{
  namespace FPSI
  {
    // type of partitioned coupling for FPSI problems
    enum PartitionedCouplingMethod
    {
      RobinNeumann,
      monolithic,
      nocoupling
    };

    // type of norm to check for convergence
    enum ConvergenceNorm
    {
      absoluteconvergencenorm,            // compare absolute value with single tolerance
      absoluteconvergencenorm_sys_split,  // compare absolute value with correction of systemsize
                                          // with different tolerances for each field
      relativconvergencenorm_sys  // compare relative value with correction of systemsize with
                                  // single tolerance
    };

    // type of norm to check for convergence
    enum BinaryOp
    {
      bop_and,
      bop_or
    };

    enum FluidFieldHierachy
    {
      fluid,
      porofluid
    };

    /// set the fpsi parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set specific fpsi conditions
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

  }  // namespace FPSI

}  // namespace INPAR

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
