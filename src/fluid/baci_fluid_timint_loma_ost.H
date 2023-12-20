/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integrator for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef BACI_FLUID_TIMINT_LOMA_OST_H
#define BACI_FLUID_TIMINT_LOMA_OST_H


#include "baci_config.H"

#include "baci_fluid_timint_loma.H"
#include "baci_fluid_timint_ost.H"
#include "baci_linalg_utils_sparse_algebra_math.H"

BACI_NAMESPACE_OPEN


namespace FLD
{
  class TimIntLomaOst : public TimIntOneStepTheta, public TimIntLoma
  {
   public:
    /// Standard Constructor
    TimIntLomaOst(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void Init() override;


   protected:
   private:
  };  // class TimIntLomaOst

}  // namespace FLD


BACI_NAMESPACE_CLOSE

#endif  // FLUID_TIMINT_LOMA_OST_H
