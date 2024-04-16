/*-----------------------------------------------------------*/
/*! \file

\brief Gen-alpha time integrator for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_LOMA_GENALPHA_HPP
#define FOUR_C_FLUID_TIMINT_LOMA_GENALPHA_HPP


#include "baci_config.hpp"

#include "baci_fluid_timint_genalpha.hpp"
#include "baci_fluid_timint_loma.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"

BACI_NAMESPACE_OPEN


namespace FLD
{
  class TimIntLomaGenAlpha : public TimIntGenAlpha, public TimIntLoma
  {
   public:
    /// Standard Constructor
    TimIntLomaGenAlpha(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void Init() override;


   protected:
   private:
  };  // class TimIntLomaGenAlpha

}  // namespace FLD


BACI_NAMESPACE_CLOSE

#endif
