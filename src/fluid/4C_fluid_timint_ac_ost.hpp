/*-----------------------------------------------------------*/
/*! \file

\brief OST time integrator for FS3I-AC problems. Is in diamond inhertance with TimIntOneStepTheta
       and TimIntAC.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_AC_OST_HPP
#define FOUR_C_FLUID_TIMINT_AC_OST_HPP


#include "4C_config.hpp"

#include "4C_fluid_timint_ac.hpp"
#include "4C_fluid_timint_ost.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  class TimIntACOst : public TimIntOneStepTheta, public TimIntAC
  {
   public:
    /// Standard Constructor
    TimIntACOst(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<CORE::IO::DiscretizationWriter>& output, bool alefluid = false);


    /// read restart data
    void read_restart(int step) override;

   protected:
   private:
  };  // class TimIntACOst

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
