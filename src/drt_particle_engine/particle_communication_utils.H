/*---------------------------------------------------------------------------*/
/*! \file
\brief communication utils for particle problem

\level 1

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_COMMUNICATION_UTILS_H
#define PARTICLE_COMMUNICATION_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <Epetra_Comm.h>
#include "Epetra_MpiComm.h"

#include <map>
#include <set>
#include <vector>

namespace PARTICLEENGINE
{
  namespace COMMUNICATION
  {
    /*!
     * \brief communicate data via non-buffered send from processor to processor
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
#endif
