/*----------------------------------------------------------------------*/
/*!
\file scatra_utils.cpp

\brief utility functions for scalar transport problems

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_utils.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/// does the given scatratype belong to the ELCH group of problems?
bool SCATRA::IsElch(const enum INPAR::SCATRA::ScaTraType scatratype)
{
  return ((scatratype==INPAR::SCATRA::scatratype_elch_enc)
      or (scatratype==INPAR::SCATRA::scatratype_elch_enc_pde)
      or (scatratype==INPAR::SCATRA::scatratype_elch_enc_pde_elim)
      or (scatratype==INPAR::SCATRA::scatratype_elch_poisson)
      or (scatratype==INPAR::SCATRA::scatratype_elch_laplace));
};

