/*-----------------------------------------------------------*/
/*! \file

\brief Stationary problem driver for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_PORO_STAT_H
#define FLD_TIMINT_PORO_STAT_H


#include "fluid_timint_poro.H"
#include "fluid_timint_stat.H"

#include "linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntPoroStat : public TimIntStationary, public TimIntPoro
  {
   public:
    //! Standard Constructor
    TimIntPoroStat(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /*!
    \brief initialization

    */
    void Init() override;

    /*!
    \brief return scheme-specific time integration parameter
    */
    double TimIntParam() const override { return 0.0; };

    /*!
    \brief read restart data
    */
    void ReadRestart(int step) override;
  };

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_PORO_STAT_H
