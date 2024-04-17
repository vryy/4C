/*-----------------------------------------------------------*/
/*! \file

\brief Generalized-alpha time integration for reduced models


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_RED_GENALPHA_HPP
#define FOUR_C_FLUID_TIMINT_RED_GENALPHA_HPP


#include "baci_config.hpp"

#include "baci_fluid_timint_genalpha.hpp"
#include "baci_fluid_timint_red.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntRedModelsGenAlpha : public TimIntGenAlpha, public TimIntRedModels
  {
   public:
    /// Standard Constructor
    TimIntRedModelsGenAlpha(const Teuchos::RCP<DRT::Discretization>& actdis,
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
  };  // class TimIntRedModelsGenAlpha

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
