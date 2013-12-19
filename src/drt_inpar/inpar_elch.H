/*----------------------------------------------------------------------*/
/*!
\file inpar_elch.H

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_ELCH_H
#define INPAR_ELCH_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

  /// Input parameters for ELCH
  namespace ELCH
  {

  /// type of elch problem
  enum ElchType
  {
    elchtype_undefined,
    elchtype_enc,
    elchtype_enc_pde,
    elchtype_enc_pde_elim,
    elchtype_poisson,
    elchtype_laplace,
    elchtype_diffcond
  };

  /// possible types of moving boundary simulation
  enum ElchMovingBoundary
  {
    elch_mov_bndry_no,
    elch_mov_bndry_pseudo_transient,
    elch_mov_bndry_fully_transient,
  };

  /// type of scalar transport problem
  enum EquPot
  {
    equpot_enc,
    equpot_divi
  };

  }  // namespace ELCH

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_ELCH_H
