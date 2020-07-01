/*----------------------------------------------------------------------*/
/*! \file

\brief Definition of the serial node type for MueLu
\level 2

*----------------------------------------------------------------------*/

#ifndef MUELU_NODEDEFINITION_HPP_

#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> KokkosSerialNode;

#endif /* MUELU_ITERATIONAFACTORY_DECL_HPP_ */
