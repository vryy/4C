/*---------------------------------------------------------------------------*/
/*! \file
\brief utility methods for particle structure interaction

\level 3

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PASI_UTILS_H
#define PASI_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

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
#endif
