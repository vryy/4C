/*-----------------------------------------------------------*/
/*! \file

\brief OST time integrator for FS3I-AC problems. Is in diamond inhertance with TimIntOneStepTheta
       and TimIntAC.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_FLUID_TIMINT_AC_OST_H
#define BACI_FLUID_TIMINT_AC_OST_H


#include "baci_fluid_timint_ac.H"
#include "baci_fluid_timint_ost.H"

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

    /// Destructor
    virtual ~TimIntACOst();

    /// read restart data
    virtual void ReadRestart(int step);

   protected:
   private:
  };  // class TimIntACOst

}  // namespace FLD


#endif  // FLUID_TIMINT_AC_OST_H
