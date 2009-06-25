/*----------------------------------------------------------------------*/
/*!
\file so_sh8p8_eas.cpp
\brief

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* defintions */
#ifdef D_SOLID3
#ifdef CCADISCRET

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
  Epetra_SerialDenseMatrix Kap(neas_,NUMPRES_);
  data_.Add("Kap",Kap);

  // quit
  return;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
