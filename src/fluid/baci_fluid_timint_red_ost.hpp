/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration for reduced models


\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_FLUID_TIMINT_RED_OST_HPP
#define FOUR_C_FLUID_TIMINT_RED_OST_HPP


#include "baci_config.hpp"

#include "baci_fluid_timint_ost.hpp"
#include "baci_fluid_timint_red.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"

BACI_NAMESPACE_OPEN


namespace FLD
{
  class TimIntRedModelsOst : public TimIntOneStepTheta, public TimIntRedModels
  {
   public:
    /// Standard Constructor
    TimIntRedModelsOst(const Teuchos::RCP<DRT::Discretization>& actdis,
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
  };  // class TimIntRedModelsOst

}  // namespace FLD


BACI_NAMESPACE_CLOSE

#endif
