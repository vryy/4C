/*----------------------------------------------------------------------*/
/*!
\file inpar_elch.H

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
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
    equpot_undefined,
    equpot_enc,
    equpot_enc_pde,
    equpot_enc_pde_elim,
    equpot_poisson,
    equpot_laplace,
    equpot_divi
  };

  enum DiffCondMat
  {
    diffcondmat_newman,
    diffcondmat_ion,
    diffcondmat_undefined
  };

  enum ApproxElectResist
  {
    approxelctresist_relpotcur,
    approxelctresist_effleninitcond,
    approxelctresist_efflenintegcond
  };

  // faraday constant (general definition)
  const double faraday_const = 96485.3399;

  // gas constant (general definition)
  const double gas_const = 8.314472;

  // permittivität (not a general definition, but only once used)
  const double epsilon_const = 1.e-4;

  }  // namespace ELCH

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_ELCH_H
