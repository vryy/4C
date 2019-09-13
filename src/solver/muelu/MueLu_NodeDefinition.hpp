/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of the serial node type for MueLu
\level 2
\maintainer Matthias Mayr

*----------------------------------------------------------------------*/

#ifndef MUELU_NODEDEFINITION_HPP_

#ifdef HAVE_MueLu

#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> KokkosSerialNode;


#endif  // HAVE_MueLu

#endif /* MUELU_ITERATIONAFACTORY_DECL_HPP_ */
