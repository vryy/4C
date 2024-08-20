/*---------------------------------------------------------------------------*/
/*! \file
\brief utility methods for particle structure interaction
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PASI_UTILS_HPP
#define FOUR_C_PASI_UTILS_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

namespace PaSI
{
  namespace UTILS
  {
    /*!
     * \brief modification of time parameter list
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void change_time_parameter(const Epetra_Comm& comm, const Teuchos::ParameterList& pasi_params,
        Teuchos::ParameterList& particle_params, Teuchos::ParameterList& struct_params);

    /*!
     * \brief print particle structure interaction logo
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void logo();

  }  // namespace UTILS

}  // namespace PaSI

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
