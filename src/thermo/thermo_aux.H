/*----------------------------------------------------------------------*/
/*! \file
\brief Various auxiliar methods needed in thermal analysis
\level 1
*/


/*----------------------------------------------------------------------*
 | definitions                                              bborn 08/09 |
 *----------------------------------------------------------------------*/
#ifndef THERMO_AUX_H
#define THERMO_AUX_H


/*----------------------------------------------------------------------*
 | headers                                                  bborn 08/09 |
 *----------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>

#include "inpar_thermo.H"


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
#endif  // THERMO_AUX_H
