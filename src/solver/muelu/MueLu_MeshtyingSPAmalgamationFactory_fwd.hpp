/*
 * MueLu_MeshtyingSPAmalgamationFactory_fwd.hpp
 *
 *  Created on: 18.03.2013
 *      Author: wiesner
 */

#ifndef MUELU_MESHTYINGSPAMALGAMATIONFACTORY_FWD_HPP_
#define MUELU_MESHTYINGSPAMALGAMATIONFACTORY_FWD_HPP_

#ifdef HAVE_MueLu
#ifdef HAVE_Trilinos_Q1_2013

namespace MueLu {
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class MeshtyingSPAmalgamationFactory;
}

#ifndef MUELU_MESHTYINGSPAMALGAMATIONFACTORY_SHORT
#define MUELU_MESHTYINGSPAMALGAMATIONFACTORY_SHORT
#endif

#endif // Q1/2013
#endif // HAVE_MueLu

#endif /* MUELU_MESHTYINGSPAMALGAMATIONFACTORY_FWD_HPP_ */
