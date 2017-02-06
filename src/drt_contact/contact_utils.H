/*----------------------------------------------------------------------------*/
/**
\file contact_utils.H

\brief Contains a summary of contact utility functions

\maintainer Michael Hiermeier

\date Jun 15, 2016

\level 2

*/
/*----------------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_CONTACT_UTILS_H_
#define SRC_DRT_CONTACT_CONTACT_UTILS_H_

#include <Teuchos_RCP.hpp>

namespace DRT {
  class Condition;
  class DiscretizationInterface;
} // namespace DRT
namespace CONTACT {
namespace UTILS {
  /// Get the solid to solid contact conditions
  int GetContactConditions(std::vector<DRT::Condition*>& cond_group,
      const std::vector<DRT::Condition*> beamandsolidcontactconditions,
      const bool & throw_error = true);

  /// Find the solid to solid contact conditions and combine them to contact condition groups
  int GetContactConditionGroups(
      std::vector<std::vector<DRT::Condition*> >& ccond_grps,
      const DRT::DiscretizationInterface& discret_wrapper,
      const bool & throw_error = true);

  /// Combine the solid to solid contact conditions to contact condition groups
  void GetContactConditionGroups(
      std::vector<std::vector<DRT::Condition*> >& ccond_grps,
      const std::vector<DRT::Condition*>& cconds);

  /// Gather information which side is master and which side is slave
  void GetMasterSlaveSideInfo(
      std::vector<bool>& isslave,
      std::vector<bool>& isself,
      const std::vector<DRT::Condition*> cond_grp);

} // namespace UTILS
} // namespace CONTACT


#endif /* SRC_DRT_CONTACT_CONTACT_UTILS_H_ */
