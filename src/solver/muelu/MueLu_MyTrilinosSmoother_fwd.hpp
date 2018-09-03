/*!----------------------------------------------------------------------
\file MueLu_MyTrilinosSmoother_fwd.hpp

\brief MueLu smoother interface
\level 2
\maintainer Martin Kronbichler

*----------------------------------------------------------------------*/

#ifndef MUELU_MYTRILINOSSMOOTHER_FWD_HPP
#define MUELU_MYTRILINOSSMOOTHER_FWD_HPP

#ifdef HAVE_MueLu

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class MyTrilinosSmoother;
}

#ifndef MUELU_MYTRILINOSSMOOTHER_SHORT
#define MUELU_MYTRILINOSSMOOTHER_SHORT
#endif

#endif  // HAVE_MueLu

#endif  // MUELU_MYTRILINOSSMOOTHER_FWD_HPP
