/*----------------------------------------------------------------------*/
/*!
\file so_sh8p8_eas.cpp
\brief
\level 2
\maintainer Alexander Seitz
*/

/*----------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------*/
/* headers */
#include "so_sh8p8.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh8p8::EasInit()
{
  // ordinary matrices and vectors
  soh8_easinit();

  // specials
  // EAS matrix K_{alpha pres}
  Epetra_SerialDenseMatrix Kap(neas_, NUMPRES_);
  data_.Add("Kap", Kap);

  // quit
  return;
}

/*----------------------------------------------------------------------*/
