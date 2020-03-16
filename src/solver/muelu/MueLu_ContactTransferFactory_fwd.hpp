/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu transfer factory class for contact
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_CONTACTTRANSFERFACTORY_FWD_HPP_
#define MUELU_CONTACTTRANSFERFACTORY_FWD_HPP_

#include <Trilinos_version.h>
#if !(TRILINOS_MAJOR_MINOR_VERSION >= 121400) || defined(HAVE_MueLuContact)

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ContactTransferFactory;
}

#ifndef MUELU_CONTACTTRANSFERFACTORY_SHORT
#define MUELU_CONTACTTRANSFERFACTORY_SHORT
#endif

#endif  // HAVE_MueLuContact

#endif /* MUELU_CONTACTTRANSFERFACTORY_FWD_HPP_ */
