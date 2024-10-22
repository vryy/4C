// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_PARTICLE_ALGORITHM_UTILS_HPP
#define FOUR_C_PARTICLE_ALGORITHM_UTILS_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_particle_engine_typedefs.hpp"
#include "4C_utils_parameter_list.hpp"

#include <map>

FOUR_C_NAMESPACE_OPEN

namespace PARTICLEALGORITHM
{
  namespace Utils
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
    template <typename Valtype>
    void read_params_types_related_to_values(const Teuchos::ParameterList& params,
        const std::string& name, std::map<PARTICLEENGINE::TypeEnum, Valtype>& typetovalmap);

  }  // namespace Utils

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
