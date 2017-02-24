/*!----------------------------------------------------------------------
\file pasi_utils.H

\brief utility methods for particle structure interaction

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                           sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
#ifndef PASI_UTILS_H
#define PASI_UTILS_H

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*
 | pasi utils                                            sfuchs 02/2017 |
 *----------------------------------------------------------------------*/
namespace PASI
{

namespace UTILS
{

//! modification of time parameter list
void ChangeTimeParameter(
    const Epetra_Comm& comm,                    //! local epetra communicator
    const Teuchos::ParameterList& pasi_params,  //! input parameters for particle structure interaction
    Teuchos::ParameterList& particle_params,    //! input parameters for particle dynamics
    Teuchos::ParameterList& struct_params       //! input parameters for structural dynamics
    );

//! prints the particle structure interaction logo on the screen
void Logo();

} // namespace UTILS

} // namespace PASI

/*----------------------------------------------------------------------*/
#endif // PASI_UTILS_H
