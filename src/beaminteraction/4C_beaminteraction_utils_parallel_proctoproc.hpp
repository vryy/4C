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

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Discret
{
  namespace Utils
  {
    //! send data to rank map key and recv data
    void i_send_receive_any(Core::FE::Discretization& discret,
        std::map<int, std::vector<std::pair<int, std::vector<int>>>> const& toranktosenddata,
        std::vector<std::pair<int, std::vector<int>>>& recvdata);

  }  // namespace Utils
}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
