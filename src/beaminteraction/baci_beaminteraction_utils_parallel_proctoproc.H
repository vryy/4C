/*-----------------------------------------------------------*/
/*! \file

\brief utils for proc to proc communication


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_BEAMINTERACTION_UTILS_PARALLEL_PROCTOPROC_H
#define BACI_BEAMINTERACTION_UTILS_PARALLEL_PROCTOPROC_H


#include "baci_config.H"

#include <Teuchos_RCP.hpp>

#include <map>
#include <vector>

BACI_NAMESPACE_OPEN

namespace DRT
{
  // forward declarations
  class Discretization;

  namespace UTILS
  {
    //! send data to rank map key and recv data
    void ISendReceiveAny(Teuchos::RCP<DRT::Discretization> const& discret,
        std::map<int, std::vector<std::pair<int, std::vector<int>>>> const& toranktosenddata,
        std::vector<std::pair<int, std::vector<int>>>& recvdata);

  }  // namespace UTILS
}  // namespace DRT

BACI_NAMESPACE_CLOSE

#endif
