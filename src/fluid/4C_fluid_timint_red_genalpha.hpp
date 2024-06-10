/*-----------------------------------------------------------*/
/*! \file

\brief Generalized-alpha time integration for reduced models


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_RED_GENALPHA_HPP
#define FOUR_C_FLUID_TIMINT_RED_GENALPHA_HPP


#include "4C_config.hpp"

#include "4C_fluid_timint_genalpha.hpp"
#include "4C_fluid_timint_red.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntRedModelsGenAlpha : public TimIntGenAlpha, public TimIntRedModels
  {
   public:
    /// Standard Constructor
    TimIntRedModelsGenAlpha(const Teuchos::RCP<Core::FE::Discretization>& actdis,
        const Teuchos::RCP<Core::LinAlg::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void Init() override;

    /*!
    \brief read restart data

    */
    void read_restart(int step) override;


   protected:
   private:
  };  // class TimIntRedModelsGenAlpha

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
