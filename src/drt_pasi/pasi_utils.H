/*---------------------------------------------------------------------------*/
/*!
\brief utility methods for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
#ifndef PASI_UTILS_H
#define PASI_UTILS_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2017 |
 *---------------------------------------------------------------------------*/
#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

namespace PASI
{
  namespace UTILS
  {
    //! modification of time parameter list
    void ChangeTimeParameter(const Epetra_Comm& comm, const Teuchos::ParameterList& pasi_params,
        Teuchos::ParameterList& particle_params, Teuchos::ParameterList& struct_params);

    //! print particle structure interaction logo
    void Logo();

  }  // namespace UTILS

}  // namespace PASI

/*---------------------------------------------------------------------------*/
#endif
