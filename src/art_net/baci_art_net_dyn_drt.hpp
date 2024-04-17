/*----------------------------------------------------------------------*/
/*! \file
\brief Main control routine for all arterial network solvers.


\level 3

*----------------------------------------------------------------------*/
/************/
#ifndef FOUR_C_ART_NET_DYN_DRT_HPP
#define FOUR_C_ART_NET_DYN_DRT_HPP


#include "baci_config.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declaration
namespace ADAPTER
{
  class ArtNet;
}

void dyn_art_net_drt();

Teuchos::RCP<ADAPTER::ArtNet> dyn_art_net_drt(bool CoupledTo3D);


FOUR_C_NAMESPACE_CLOSE

#endif
