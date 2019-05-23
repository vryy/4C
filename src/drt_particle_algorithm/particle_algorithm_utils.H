/*---------------------------------------------------------------------------*/
/*!
\brief utils for particle algorithm

\level 3

\maintainer  Sebastian Fuchs
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
    //! read parameters relating particle types to values
    template <typename valtype>
    void ReadParamsTypesRelatedToValues(const Teuchos::ParameterList& params,
        const std::string& name, std::map<PARTICLEENGINE::TypeEnum, valtype>& typetovalmap);

  }  // namespace UTILS

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
