/*-----------------------------------------------------------*/
/*! \file

\brief utils for proc to proc communication


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef UTILS_PARALLEL_PROCTOPROC_H
#define UTILS_PARALLEL_PROCTOPROC_H


#include <map>
#include <vector>

#include "Teuchos_RCP.hpp"


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

#endif
