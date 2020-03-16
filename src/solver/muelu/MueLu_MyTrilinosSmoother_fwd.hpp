/*----------------------------------------------------------------------*/
/*! \file

\brief MueLu smoother interface
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_MYTRILINOSSMOOTHER_FWD_HPP
#define MUELU_MYTRILINOSSMOOTHER_FWD_HPP

#ifdef HAVE_MueLuContact

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class MyTrilinosSmoother;
}

#ifndef MUELU_MYTRILINOSSMOOTHER_SHORT
#define MUELU_MYTRILINOSSMOOTHER_SHORT
#endif

#endif  // HAVE_MueLuContact

#endif  // MUELU_MYTRILINOSSMOOTHER_FWD_HPP
