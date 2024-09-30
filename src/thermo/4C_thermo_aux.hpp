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
#include "4C_config.hpp"

#include "4C_inpar_thermo.hpp"
#include "4C_linalg_vector.hpp"

#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | belongs to thermal dynamics namespace                    bborn 08/09 |
 *----------------------------------------------------------------------*/
namespace Thermo
{
  /*====================================================================*/
  namespace Aux
  {
    //! Determine norm of force residual
    double calculate_vector_norm(const enum Inpar::Thermo::VectorNorm norm,  //!< norm to use
        const Teuchos::RCP<Core::LinAlg::Vector> vect  //!< the vector of interest
    );

  }  // namespace Aux

}  // namespace Thermo

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
