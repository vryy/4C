/*---------------------------------------------------------------------------*/
/*! \file
\brief utils for particle algorithm
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_ALGORITHM_UTILS_H
#define PARTICLE_ALGORITHM_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_inpar/inpar_parameterlist_utils.H"

#include <map>

namespace PARTICLEALGORITHM
{
  namespace UTILS
  {
    /*!
     * \brief read parameters relating particle types to values
     *
     * Read parameters relating particle types to specific values from the parameter list.
     *
     * \author Sebastian Fuchs \date 07/2018
     *
     * \tparam valtype type of value
     *
     * \param[in]  params       particle simulation parameter list
     * \param[in]  name         parameter name
     * \param[out] typetovalmap map relating particle types to specific values
     */
    template <typename valtype>
    void ReadParamsTypesRelatedToValues(const Teuchos::ParameterList& params,
        const std::string& name, std::map<PARTICLEENGINE::TypeEnum, valtype>& typetovalmap);

  }  // namespace UTILS

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
