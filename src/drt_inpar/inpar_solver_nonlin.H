/*----------------------------------------------------------------------*/
/*!
\file inpar_solver_nonlin.H

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289 10362
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef INPAR_SOLVER_NONLIN_H
#define INPAR_SOLVER_NONLIN_H

#include "inpar_parameterlist_utils.H"

/*----------------------------------------------------------------------*/
namespace INPAR
{

//! Input parameters for nonlinear solver package
namespace NLNSOL
{

namespace LINESEARCH
{
  //! type of sufficient decrease condition
  enum SufficientDeacreaseType
  {
    suffdecr_none, ///< no type chosen
    suffdecr_armijo, ///< Armijo's rule involving gradient information
    suffdecr_aredpred, ///< ration of actual reduction to predicted reduction
    suffdecr_loose ///< loose conditione: \f$||F(x_k+\lambda p)|| \leq \kappa||F(x_k)||\f$ with \f$\kappa\geq0\f$
  };
} // namespace LINESARCH

namespace NONLINCG
{
  //! possible formulas to compute beta within a nonlinear CG
  enum BetaType
  {
    beta_none, ///< no formulation chosen
    beta_fletcherreeves, ///< Fletcher-Reeves: \f$\beta^{FR} = \frac{x_k^Tf_k}{x_{k-1}^Tf_{k-1}}\f$
    beta_polakribiere, ///< Polak-Ribiere: \f$\beta^{PR} = \frac{F_k^TM^{-1}(F_k-F_{k-1})}{F_{k-1}^TM^{-1}F_{k-1}}\f$
    beta_hestenesstiefel ///< Hestenes-Stiefel:
  };
} // namespace NONLINCG
} // namespace NLNSOL
} // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_SOLVER_NONLIN_H
