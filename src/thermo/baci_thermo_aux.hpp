/*----------------------------------------------------------------------*/
/*! \file
\brief Various auxiliar methods needed in thermal analysis
\level 1
*/


/*----------------------------------------------------------------------*
 | definitions                                              bborn 08/09 |
 *----------------------------------------------------------------------*/
#ifndef FOUR_C_THERMO_AUX_HPP
#define FOUR_C_THERMO_AUX_HPP


/*----------------------------------------------------------------------*
 | headers                                                  bborn 08/09 |
 *----------------------------------------------------------------------*/
#include "baci_config.hpp"

#include "baci_inpar_thermo.hpp"

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

#endif
