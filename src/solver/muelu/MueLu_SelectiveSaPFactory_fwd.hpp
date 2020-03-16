/*----------------------------------------------------------------------*/
/*! \file

\brief Factory class for MueLu
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_SELECTIVESAPFACTORY_FWD_HPP_
#define MUELU_SELECTIVESAPFACTORY_FWD_HPP_


#ifdef HAVE_MueLuContact

namespace MueLu
{
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class SelectiveSaPFactory;
}

#ifndef MUELU_SELECTIVESAPFACTORY_SHORT
#define MUELU_SELECTIVESAPFACTORY_SHORT
#endif

#endif  // HAVE_MueLuContact

#endif /* MUELU_SELECTIVESAPFACTORY_FWD_HPP_ */
