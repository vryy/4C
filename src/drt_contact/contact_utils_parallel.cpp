/*----------------------------------------------------------------------------*/
/*! \file
\brief Contact utility functions related to parallel runs

\maintainer Matthias Mayr

\level 1

*/
/*----------------------------------------------------------------------------*/

#include "contact_utils_parallel.H"

#include "../drt_inpar/inpar_contact.H"
#include "../drt_inpar/inpar_mortar.H"
#include "../drt_inpar/inpar_structure.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_globalproblem_enums.H"

#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool CONTACT::UTILS::UseSafeRedistributeAndGhosting(const Teuchos::ParameterList& contactParams)
{
  using namespace INPAR::CONTACT;

  /* Limit the use of the new safe "redistribute & ghosting" branch to our core contact
   * capabilities. If your case of interest is missing here, feel free to migrate your scenario to
   * the new safe branch.
   */
  bool use_safe_ghosting_branch = false;
  {
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    const enum INPAR::STR::IntegrationStrategy intstrat =
        DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(sdyn, "INT_STRATEGY");

    if (intstrat == INPAR::STR::int_old)
    {
      /* Enable new safe ghosting only for interface discretization type "mortar"
       *
       * There's a conflict with CreateVolumeGhosting(). This affects all Nitsche-type algorithms
       * and also classical Penalty with Gauss-point-to-segment (GPTS).
       *
       * In theory, penalty with GPTS should work just fine, because it should never need a volume
       * ghosting. However, penalty with GPTS is implemented in the CoNitscheStrategy, which
       * always requires volume ghosting.
       *
       * Furthermore, everything porous media related has to stick to the old code branch as well.
       */
      if (DRT::INPUT::IntegralValue<INPAR::MORTAR::AlgorithmType>(contactParams, "ALGORITHM") ==
              INPAR::MORTAR::algorithm_mortar &&
          DRT::Problem::Instance()->GetProblemType() != prb_poroelast)
        use_safe_ghosting_branch = true;
    }
    else
    {
      // Use old code path, if structure uses the new time integration.
    }
  }

  return use_safe_ghosting_branch;
}
