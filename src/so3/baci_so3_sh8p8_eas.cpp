/*----------------------------------------------------------------------*/
/*! \file
\brief 8-node solid shell element
\level 2
*/

/*----------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_so3_sh8p8.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoSh8p8::EasInit()
{
  // ordinary matrices and vectors
  soh8_easinit();

  // specials
  // EAS matrix K_{alpha pres}
  CORE::LINALG::SerialDenseMatrix Kap(neas_, NUMPRES_);
  easdata_.Kap = Kap;
}

/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
