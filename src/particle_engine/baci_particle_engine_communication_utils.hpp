/*---------------------------------------------------------------------------*/
/*! \file
\brief communication utils for particle problem
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_ENGINE_COMMUNICATION_UTILS_HPP
#define BACI_PARTICLE_ENGINE_COMMUNICATION_UTILS_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.hpp"

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>

#include <map>
#include <set>
#include <vector>

BACI_NAMESPACE_OPEN

namespace PARTICLEENGINE
{
  namespace COMMUNICATION
  {
    /*!
     * \brief communicate data via non-buffered send from processor to processor
     *
     * Communicate data via non-buffered send from processor to processor via point-to-point
     * communication. Collective communication is required to inform the target processors about the
     * size of data to be received. The send buffer only relates data to be send to the target
     * processor.
     *
     * \note This method has to be called by all processors of the communicator as it contains
     * collective communication.
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in]  comm  communicator
     * \param[in]  sdata send buffers related to corresponding target processors
     * \param[out] rdata receive buffers related to corresponding source processors
     */
    void ImmediateRecvBlockingSend(const Epetra_Comm& comm, std::map<int, std::vector<char>>& sdata,
        std::map<int, std::vector<char>>& rdata);

  }  // namespace COMMUNICATION

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
