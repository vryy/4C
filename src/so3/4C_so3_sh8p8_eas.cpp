/*----------------------------------------------------------------------*/
/*! \file
\brief 8-node solid shell element
\level 2
*/

/*----------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_so3_sh8p8.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Discret::ELEMENTS::SoSh8p8::eas_init()
{
  // ordinary matrices and vectors
  soh8_easinit();

  // specials
  // EAS matrix K_{alpha pres}
  Core::LinAlg::SerialDenseMatrix Kap(neas_, NUMPRES_);
  easdata_.Kap = Kap;
}

/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
