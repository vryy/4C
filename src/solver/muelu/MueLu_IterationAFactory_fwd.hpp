/*!----------------------------------------------------------------------
\file MueLu_IterationAFactory_fwd.hpp

\brief MueLu iteration factory class
\level 2
\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#ifndef MUELU_ITERATIONAFACTORY_FWD_HPP_
#define MUELU_ITERATIONAFACTORY_FWD_HPP_

#ifdef HAVE_MueLu

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class IterationAFactory;
}

#ifndef MUELU_ITERATIONAFACTORY_SHORT
#define MUELU_ITERATIONAFACTORY_SHORT
#endif

#endif  // HAVE_MueLu

#endif /* MUELU_ITERATIONAFACTORY_FWD_HPP_ */
