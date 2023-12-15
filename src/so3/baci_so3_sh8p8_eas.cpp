/*----------------------------------------------------------------------*/
/*! \file
\brief 8-node solid shell element
\level 2
*/

/*----------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------*/
/* headers */
#include "baci_so3_sh8p8.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::EasInit()
{
  // ordinary matrices and vectors
  soh8_easinit();

  // specials
  // EAS matrix K_{alpha pres}
  CORE::LINALG::SerialDenseMatrix Kap(neas_, NUMPRES_);
  data_.Add("Kap", Kap);

  // quit
  return;
}

/*----------------------------------------------------------------------*/

BACI_NAMESPACE_CLOSE
