/*---------------------------------------------------------------------------*/
/*! \file
\brief utils for particle algorithm
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_ALGORITHM_UTILS_H
#define BACI_PARTICLE_ALGORITHM_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.H"

#include "baci_inpar_parameterlist_utils.H"
#include "baci_particle_engine_typedefs.H"

#include <map>

BACI_NAMESPACE_OPEN

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
BACI_NAMESPACE_CLOSE

#endif
