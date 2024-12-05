// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_PASI_UTILS_HPP
#define FOUR_C_PASI_UTILS_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_utils_parameter_list.fwd.hpp"

#include <mpi.h>


FOUR_C_NAMESPACE_OPEN

namespace PaSI
{
  namespace Utils
  {
    /*!
     * \brief modification of time parameter list
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void change_time_parameter(MPI_Comm comm, const Teuchos::ParameterList& pasi_params,
        Teuchos::ParameterList& particle_params, Teuchos::ParameterList& struct_params);

    /*!
     * \brief print particle structure interaction logo
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void logo();

  }  // namespace Utils

}  // namespace PaSI

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
