/*----------------------------------------------------------------------*/
/*!
\file thr_aux.H
\brief Various auxiliar methods needed in thermal analysis
\level 1
<pre>
\maintainer Christoph Meier

</pre>
*/


/*----------------------------------------------------------------------*
 | definitions                                              bborn 08/09 |
 *----------------------------------------------------------------------*/
#ifndef THR_AUX_H
#define THR_AUX_H


/*----------------------------------------------------------------------*
 | headers                                                  bborn 08/09 |
 *----------------------------------------------------------------------*/
#include "Teuchos_RCP.hpp"
#include "Epetra_Vector.h"

#include "../drt_inpar/inpar_thermo.H"


/*----------------------------------------------------------------------*
 | belongs to thermal dynamics namespace                    bborn 08/09 |
 *----------------------------------------------------------------------*/
namespace THR
{
  /*====================================================================*/
  //! A namespace for the various unsorted auxiliar kind of
  //! methods needed in thermal analysis
  //!
  //! \author bborn
  //! \date 11/08
  namespace AUX
  {
    //! Determine norm of force residual
    double CalculateVectorNorm(const enum INPAR::THR::VectorNorm norm,  //!< norm to use
        const Teuchos::RCP<Epetra_Vector> vect                          //!< the vector of interest
    );

  }  // namespace AUX

}  // namespace THR

/*----------------------------------------------------------------------*/
#endif  // #ifndef THR_AUX_H
