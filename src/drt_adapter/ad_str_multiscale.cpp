/*----------------------------------------------------------------------*/
/*!
\file ad_str_multiscale.cpp

\brief Wrapper for the structural time loop for multiscale problems

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/

#include "ad_str_multiscale.H"
#include "../drt_stru_multi/microstatic.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureMultiScale::PostTimeLoop()
{
  // stop supporting processors in multi scale simulations
  STRUMULTI::stop_np_multiscale();

  return;
}

