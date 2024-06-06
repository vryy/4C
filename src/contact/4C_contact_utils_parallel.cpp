/*----------------------------------------------------------------------------*/
/*! \file
\brief Contact utility functions related to parallel runs


\level 1

*/
/*----------------------------------------------------------------------------*/

#include "4C_contact_utils_parallel.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_mortar.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_legacy_enum_definitions_problem_type.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::UTILS::UseSafeRedistributeAndGhosting(const Teuchos::ParameterList& contactParams)
{
  /* Limit the use of the new safe "redistribute & ghosting" branch to our core contact
   * capabilities. If your case of interest is missing here, feel free to migrate your scenario to
   * the new safe branch.
   */
  bool use_safe_ghosting_branch = false;
  {
    const Teuchos::ParameterList& sdyn = GLOBAL::Problem::Instance()->structural_dynamic_params();
    const enum INPAR::STR::IntegrationStrategy intstrat =
        CORE::UTILS::IntegralValue<INPAR::STR::IntegrationStrategy>(sdyn, "INT_STRATEGY");

    if (intstrat == INPAR::STR::int_old)
    {
      /* Enable new safe ghosting only for interface discretization type "mortar"
       *
       * There's a conflict with create_volume_ghosting(). This affects all Nitsche-type algorithms
       * and also classical Penalty with Gauss-point-to-segment (GPTS).
       *
       * In theory, penalty with GPTS should work just fine, because it should never need a volume
       * ghosting. However, penalty with GPTS is implemented in the NitscheStrategy, which
       * always requires volume ghosting.
       *
       * Other cases require volume ghosting as well and, thus, have to stick to the old code
       * branch. They are:
       * - Everything porous media related has to stick to the old code branch as well.
       * - "Large" wear, i.e. using Structure-ALE.
       */
      if (CORE::UTILS::IntegralValue<INPAR::MORTAR::AlgorithmType>(contactParams, "ALGORITHM") ==
              INPAR::MORTAR::algorithm_mortar &&
          (GLOBAL::Problem::Instance()->GetProblemType() != CORE::ProblemType::poroelast &&
              GLOBAL::Problem::Instance()->GetProblemType() != CORE::ProblemType::poroscatra &&
              GLOBAL::Problem::Instance()->GetProblemType() != CORE::ProblemType::struct_ale))
        use_safe_ghosting_branch = true;
    }
    else
    {
      // Use old code path, if structure uses the new time integration.
    }
  }

  return use_safe_ghosting_branch;
}

FOUR_C_NAMESPACE_CLOSE
