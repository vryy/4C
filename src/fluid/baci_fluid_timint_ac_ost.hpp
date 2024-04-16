/*-----------------------------------------------------------*/
/*! \file

\brief OST time integrator for FS3I-AC problems. Is in diamond inhertance with TimIntOneStepTheta
       and TimIntAC.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_AC_OST_HPP
#define FOUR_C_FLUID_TIMINT_AC_OST_HPP


#include "baci_config.hpp"

#include "baci_fluid_timint_ac.hpp"
#include "baci_fluid_timint_ost.hpp"

BACI_NAMESPACE_OPEN

namespace FLD
{
  class TimIntACOst : public TimIntOneStepTheta, public TimIntAC
  {
   public:
    /// Standard Constructor
    TimIntACOst(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);


    /// read restart data
    void ReadRestart(int step) override;

   protected:
   private:
  };  // class TimIntACOst

}  // namespace FLD


BACI_NAMESPACE_CLOSE

#endif
