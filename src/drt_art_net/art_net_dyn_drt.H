/*!----------------------------------------------------------------------
\brief Main control routine for all arterial network solvers.

\maintainer Johannes Kremheller

\level 3

*----------------------------------------------------------------------*/
/************/
#ifndef ART_NET_DYN_DRT_H
#define ART_NET_DYN_DRT_H


#include "Teuchos_RCP.hpp"

// forward declaration
namespace ADAPTER
{
  class ArtNet;
}

void dyn_art_net_drt();

Teuchos::RCP<ADAPTER::ArtNet> dyn_art_net_drt(bool CoupledTo3D);


#endif /* #ifndef ART_NET_DYN_DRT_H */
