/*
 * MueLu_SelectiveSaPFactory_fwd.hpp
 *
 *  Created on: Jan 17, 2013
 *      Author: wiesner
 */

#ifndef MUELU_SELECTIVESAPFACTORY_FWD_HPP_
#define MUELU_SELECTIVESAPFACTORY_FWD_HPP_


#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q1_2013

namespace MueLu {
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class SelectiveSaPFactory;
}

#ifndef MUELU_SELECTIVESAPFACTORY_SHORT
#define MUELU_SELECTIVESAPFACTORY_SHORT
#endif

#endif // #ifdef HAVE_Trilinos_Q1_2013
#endif // HAVE_MueLu

#endif /* MUELU_SELECTIVESAPFACTORY_FWD_HPP_ */
