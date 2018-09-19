/*---------------------------------------------------------------------------*/
/*!
\file particle_algorithm_utils.H

\brief utils for particle algorithm

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
#ifndef PARTICLE_ALGORITHM_UTILS_H
#define PARTICLE_ALGORITHM_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include <map>

namespace PARTICLEALGORITHM
{
  namespace UTILS
  {
    //! read parameters relating particle types to IDs
    void ReadParamsTypesRelatedToIDs(const Teuchos::ParameterList& params, const std::string& name,
        std::map<PARTICLEENGINE::TypeEnum, int>& typetoidmap);

  }  // namespace UTILS

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
