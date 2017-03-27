/*---------------------------------------------------------------------*/
/*!
\file contact_aug_utils.H

\brief Global utility file for the augmented contact formulation

\level 3

\maintainer Michael Hiermeier

\date Mar 23, 2017

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_CONTACT_AUG_CONTACT_AUG_UTILS_H_
#define SRC_DRT_CONTACT_AUG_CONTACT_AUG_UTILS_H_

#include "../drt_inpar/inpar_contact.H"

namespace CONTACT {
  class CoInterface;
  class CoAbstractStrategy;
namespace AUG {

  typedef std::vector< Teuchos::RCP<CONTACT::CoAbstractStrategy> > plain_strategy_set;
  typedef std::vector<Teuchos::RCP<CONTACT::CoInterface> > plain_interface_set;
  typedef std::vector< plain_interface_set > plain_interface_sets;
  typedef std::vector<INPAR::CONTACT::SolvingStrategy> plain_strattype_set;
} // namespace AUG
} // namespace CONTACT


#endif /* SRC_DRT_CONTACT_AUG_CONTACT_AUG_UTILS_H_ */
