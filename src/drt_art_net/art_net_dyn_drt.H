/*!----------------------------------------------------------------------
\file art_net_dyn_drt.H
\brief Control routine for arterial network (time) integration. Calls


<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/
/************/
#ifndef ART_NET_DYN_DRT_H
#define ART_NET_DYN_DRT_H


#include "Teuchos_RCP.hpp"

// forward declaration
namespace ART{
 class ArtNetExplicitTimeInt;
}

void dyn_art_net_drt();

Teuchos::RCP<ART::ArtNetExplicitTimeInt>  dyn_art_net_drt(bool CoupledTo3D);


#endif  /* #ifndef ART_NET_DYN_DRT_H */

