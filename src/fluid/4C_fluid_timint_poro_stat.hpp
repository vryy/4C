/*-----------------------------------------------------------*/
/*! \file

\brief Stationary problem driver for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_PORO_STAT_HPP
#define FOUR_C_FLUID_TIMINT_PORO_STAT_HPP


#include "4C_config.hpp"

#include "4C_fluid_timint_poro.hpp"
#include "4C_fluid_timint_stat.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntPoroStat : public TimIntStationary, public TimIntPoro
  {
   public:
    //! Standard Constructor
    TimIntPoroStat(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<CORE::IO::DiscretizationWriter>& output, bool alefluid = false);

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
    void read_restart(int step) override;
  };

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
