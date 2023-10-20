/*-----------------------------------------------------------*/
/*! \file

\brief Gen-alpha time integrator for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef BACI_FLUID_TIMINT_LOMA_GENALPHA_H
#define BACI_FLUID_TIMINT_LOMA_GENALPHA_H


#include "baci_fluid_timint_genalpha.H"
#include "baci_fluid_timint_loma.H"
#include "baci_linalg_utils_sparse_algebra_math.H"


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


#endif  // FLUID_TIMINT_LOMA_GENALPHA_H
