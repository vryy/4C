/*-----------------------------------------------------------*/
/*! \file

\brief utils for proc to proc communication


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_BEAMINTERACTION_UTILS_PARALLEL_PROCTOPROC_HPP
#define FOUR_C_BEAMINTERACTION_UTILS_PARALLEL_PROCTOPROC_HPP


#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

#include <map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

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

FOUR_C_NAMESPACE_CLOSE

#endif
