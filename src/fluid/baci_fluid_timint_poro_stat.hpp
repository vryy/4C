/*-----------------------------------------------------------*/
/*! \file

\brief Stationary problem driver for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_PORO_STAT_HPP
#define FOUR_C_FLUID_TIMINT_PORO_STAT_HPP


#include "baci_config.hpp"

#include "baci_fluid_timint_poro.hpp"
#include "baci_fluid_timint_stat.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"

BACI_NAMESPACE_OPEN


namespace FLD
{
  class TimIntPoroStat : public TimIntStationary, public TimIntPoro
  {
   public:
    //! Standard Constructor
    TimIntPoroStat(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
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


BACI_NAMESPACE_CLOSE

#endif
