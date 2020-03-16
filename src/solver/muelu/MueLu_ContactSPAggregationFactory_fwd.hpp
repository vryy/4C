/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu contact aggregation factory class for saddle point formulations
\level 2
\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTSPAGGREGATIONFACTORY_FWD_HPP_
#define MUELU_CONTACTSPAGGREGATIONFACTORY_FWD_HPP_

#include <Trilinos_version.h>
#if !(TRILINOS_MAJOR_MINOR_VERSION >= 121400) || defined(HAVE_MueLuContact)

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ContactSPAggregationFactory;
}

#ifndef MUELU_CONTACTSPAGGREGATIONFACTORY_SHORT
#define MUELU_CONTACTSPAGGREGATIONFACTORY_SHORT
#endif

#endif  // HAVE_MueLuContact

#endif /* MUELU_CONTACTSPAGGREGATIONFACTORY_FWD_HPP_ */
