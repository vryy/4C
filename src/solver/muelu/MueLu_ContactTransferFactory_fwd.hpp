/*!----------------------------------------------------------------------
\file MueLu_ContactTransferFactory_fwd.hpp

\brief MueLu transfer factory class for contact
\level 2
\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTTRANSFERFACTORY_FWD_HPP_
#define MUELU_CONTACTTRANSFERFACTORY_FWD_HPP_

#ifdef HAVE_MueLu

namespace MueLu {
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ContactTransferFactory;
}

#ifndef MUELU_CONTACTTRANSFERFACTORY_SHORT
#define MUELU_CONTACTTRANSFERFACTORY_SHORT
#endif

#endif // HAVE_MueLu

#endif /* MUELU_CONTACTTRANSFERFACTORY_FWD_HPP_ */
