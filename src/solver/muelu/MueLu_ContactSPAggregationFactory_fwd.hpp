/*!----------------------------------------------------------------------
\file MueLu_ContactSPAggregationFactory_fwd.hpp

\brief MueLu contact aggregation factory class for saddle point formulations
\level 2
\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTSPAGGREGATIONFACTORY_FWD_HPP_
#define MUELU_CONTACTSPAGGREGATIONFACTORY_FWD_HPP_

#ifdef HAVE_MueLu

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ContactSPAggregationFactory;
}

#ifndef MUELU_CONTACTSPAGGREGATIONFACTORY_SHORT
#define MUELU_CONTACTSPAGGREGATIONFACTORY_SHORT
#endif

#endif  // HAVE_MueLu

#endif /* MUELU_CONTACTSPAGGREGATIONFACTORY_FWD_HPP_ */
