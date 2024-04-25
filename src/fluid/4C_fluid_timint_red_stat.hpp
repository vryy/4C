/*-----------------------------------------------------------*/
/*! \file

\brief Stationary driver for reduced models


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_RED_STAT_HPP
#define FOUR_C_FLUID_TIMINT_RED_STAT_HPP


#include "4C_config.hpp"

#include "4C_fluid_timint_red.hpp"
#include "4C_fluid_timint_stat.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntRedModelsStat : public TimIntStationary, public TimIntRedModels
  {
   public:
    /// Standard Constructor
    TimIntRedModelsStat(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void Init() override;

    /*!
    \brief read restart data

    */
    void ReadRestart(int step) override;


   protected:
   private:
  };  // class TimIntRedModelsStat

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
