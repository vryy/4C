/*---------------------------------------------------------------------------*/
/*!
\file particle_communication_utils.H

\brief communication utils for particle engine

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_COMMUNICATION_UTILS_H
#define PARTICLE_COMMUNICATION_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
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
    //! communicate data via non-buffered send from proc to proc
    void ImmediateRecvBlockingSend(const Epetra_Comm& comm, std::map<int, std::vector<char>>& sdata,
        std::map<int, std::vector<char>>& rdata);

  }  // namespace COMMUNICATION

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
