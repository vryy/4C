/*-----------------------------------------------------------*/
/*! \file

\brief Gen-alpha time integrator for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_LOMA_GENALPHA_HPP
#define FOUR_C_FLUID_TIMINT_LOMA_GENALPHA_HPP


#include "4C_config.hpp"

#include "4C_fluid_timint_genalpha.hpp"
#include "4C_fluid_timint_loma.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntLomaGenAlpha : public TimIntGenAlpha, public TimIntLoma
  {
   public:
    /// Standard Constructor
    TimIntLomaGenAlpha(const Teuchos::RCP<Core::FE::Discretization>& actdis,
        const Teuchos::RCP<Core::LinAlg::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void init() override;


   protected:
   private:
  };  // class TimIntLomaGenAlpha

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
