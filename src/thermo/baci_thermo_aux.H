/*----------------------------------------------------------------------*/
/*! \file
\brief Various auxiliar methods needed in thermal analysis
\level 1
*/


/*----------------------------------------------------------------------*
 | definitions                                              bborn 08/09 |
 *----------------------------------------------------------------------*/
#ifndef BACI_THERMO_AUX_H
#define BACI_THERMO_AUX_H


/*----------------------------------------------------------------------*
 | headers                                                  bborn 08/09 |
 *----------------------------------------------------------------------*/
#include "baci_inpar_thermo.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | belongs to thermal dynamics namespace                    bborn 08/09 |
 *----------------------------------------------------------------------*/
namespace THR
{
  /*====================================================================*/
  namespace AUX
  {
    //! Determine norm of force residual
    double CalculateVectorNorm(const enum INPAR::THR::VectorNorm norm,  //!< norm to use
        const Teuchos::RCP<Epetra_Vector> vect                          //!< the vector of interest
    );

  }  // namespace AUX

}  // namespace THR

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // THERMO_AUX_H
