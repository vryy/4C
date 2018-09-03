/*----------------------------------------------------------------------*/
/*!
\file inpar_solver_nonlin.H


\brief Input parameters for nonlinear solvers

\maintainer Michael Hiermeier

\level 1
*/

/*----------------------------------------------------------------------*/

#ifndef INPAR_SOLVER_NONLIN_H
#define INPAR_SOLVER_NONLIN_H

/*----------------------------------------------------------------------*/
namespace INPAR
{
  //! Input parameters for nonlinear solver package
  namespace NLNSOL
  {
    namespace FAS
    {
      //! type of multigrid cycle
      enum CycleType
      {
        cycle_none,           ///< no type chosen
        cycle_v,              ///< V-cycle (recursive definition)
        cycle_v_two_level,    ///< V-cycle with 2 levels (explicit definition)
        cycle_v_three_level,  ///< V-cycle with 3 levels (explicit definition)
        cycle_w               ///< W-cycle
      };
    }  // namespace FAS

    namespace LINESEARCH
    {
      //! type of sufficient decrease condition
      enum SufficientDeacreaseType
      {
        suffdecr_none,      ///< no type chosen
        suffdecr_armijo,    ///< Armijo's rule involving gradient information
        suffdecr_aredpred,  ///< ration of actual reduction to predicted reduction
        suffdecr_loose  ///< loose conditione: \f$||F(x_k+\lambda p)|| \leq \kappa||F(x_k)||\f$ with
                        ///< \f$\kappa\geq0\f$
      };
    }  // namespace LINESEARCH

    namespace NONLINCG
    {
      //! possible formulas to compute beta within a nonlinear CG
      enum BetaType
      {
        beta_none,            ///< no formulation chosen
        beta_fletcherreeves,  ///< Fletcher-Reeves: \f$\beta^{FR} =
                              ///< \frac{x_k^Tf_k}{x_{k-1}^Tf_{k-1}}\f$
        beta_polakribiere,    ///< Polak-Ribiere: \f$\beta^{PR} =
                              ///< \frac{F_k^TM^{-1}(F_k-F_{k-1})}{F_{k-1}^TM^{-1}F_{k-1}}\f$
        beta_hestenesstiefel  ///< Hestenes-Stiefel:
      };
    }  // namespace NONLINCG


    /// set the nonlinear solver parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace NLNSOL
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_SOLVER_NONLIN_H
