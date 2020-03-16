/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu factory class for BACI
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/
#ifndef MUELU_CONTACTAFILTERFACTORY_FWD_HPP
#define MUELU_CONTACTAFILTERFACTORY_FWD_HPP

#include <Trilinos_version.h>
#if !(TRILINOS_MAJOR_MINOR_VERSION >= 121400) || defined(HAVE_MueLuContact)

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ContactAFilterFactory;
}

#ifndef MUELU_CONTACTAFILTERFACTORY_SHORT
#define MUELU_CONTACTAFILTERFACTORY_SHORT
#endif

#endif  // HAVE_MueLuContact

#endif  // MUELU_CONTACTAFILTERFACTORY_FWD_HPP
