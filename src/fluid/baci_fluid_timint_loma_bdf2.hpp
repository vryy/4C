/*-----------------------------------------------------------*/
/*! \file

\brief BDF2 time integrator for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_LOMA_BDF2_HPP
#define FOUR_C_FLUID_TIMINT_LOMA_BDF2_HPP


#include "baci_config.hpp"

#include "baci_fluid_timint_bdf2.hpp"
#include "baci_fluid_timint_loma.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"

BACI_NAMESPACE_OPEN


namespace FLD
{
  class TimIntLomaBDF2 : public TimIntBDF2, public TimIntLoma
  {
   public:
    /// Standard Constructor
    TimIntLomaBDF2(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void Init() override;


   protected:
   private:
  };  // class TimIntLomaBDF2

}  // namespace FLD


BACI_NAMESPACE_CLOSE

#endif  // FLUID_TIMINT_LOMA_BDF2_H
