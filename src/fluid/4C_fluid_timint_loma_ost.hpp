/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integrator for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_LOMA_OST_HPP
#define FOUR_C_FLUID_TIMINT_LOMA_OST_HPP


#include "4C_config.hpp"

#include "4C_fluid_timint_loma.hpp"
#include "4C_fluid_timint_ost.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntLomaOst : public TimIntOneStepTheta, public TimIntLoma
  {
   public:
    /// Standard Constructor
    TimIntLomaOst(const Teuchos::RCP<Discret::Discretization>& actdis,
        const Teuchos::RCP<Core::LinAlg::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void Init() override;


   protected:
   private:
  };  // class TimIntLomaOst

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
