/*
 * MueLu_NodeDefinition.hpp
 *
 * Maintainer: Martin Kronbichler
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15235
 */

#ifndef MUELU_NODEDEFINITION_HPP_

#ifdef HAVE_MueLu

#include <KokkosCompat_ClassicNodeAPI_Wrapper.hpp>

typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> KokkosSerialNode;


#endif // HAVE_MueLu

#endif /* MUELU_ITERATIONAFACTORY_DECL_HPP_ */
