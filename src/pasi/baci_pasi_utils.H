/*---------------------------------------------------------------------------*/
/*! \file
\brief utility methods for particle structure interaction
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PASI_UTILS_H
#define BACI_PASI_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_config.H"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN

namespace PASI
{
  namespace UTILS
  {
    /*!
     * \brief modification of time parameter list
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void ChangeTimeParameter(const Epetra_Comm& comm, const Teuchos::ParameterList& pasi_params,
        Teuchos::ParameterList& particle_params, Teuchos::ParameterList& struct_params);

    /*!
     * \brief print particle structure interaction logo
     *
     * \author Sebastian Fuchs \date 02/2017
     */
    void Logo();

  }  // namespace UTILS

}  // namespace PASI

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
