/*----------------------------------------------------------------------*/
/*! \file
\brief Main control routine for all arterial network solvers.


\level 3

*----------------------------------------------------------------------*/
/************/
#ifndef BACI_ART_NET_DYN_DRT_H
#define BACI_ART_NET_DYN_DRT_H


#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

// forward declaration
namespace ADAPTER
{
  class ArtNet;
}

void dyn_art_net_drt();

Teuchos::RCP<ADAPTER::ArtNet> dyn_art_net_drt(bool CoupledTo3D);


BACI_NAMESPACE_CLOSE

#endif  // ART_NET_DYN_DRT_H
