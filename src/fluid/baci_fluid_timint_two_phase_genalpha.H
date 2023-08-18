/*-----------------------------------------------------------*/
/*! \file

\brief Gen-alpha time integration for two-phase flow


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef BACI_FLUID_TIMINT_TWO_PHASE_GENALPHA_H
#define BACI_FLUID_TIMINT_TWO_PHASE_GENALPHA_H


#include "baci_fluid_timint_genalpha.H"
#include "baci_fluid_timint_two_phase.H"
#include "baci_linalg_utils_sparse_algebra_math.H"


namespace FLD
{
  class TimIntTwoPhaseGenAlpha : public TimIntGenAlpha, public TimIntTwoPhase
  {
   public:
    /// Standard Constructor
    TimIntTwoPhaseGenAlpha(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntTwoPhaseGenAlpha();

    /*!
    \brief initialization

    */
    virtual void Init();


   protected:
   private:
  };  // class TimIntTwoPhaseGenAlpha

}  // namespace FLD


#endif  // FLUID_TIMINT_TWO_PHASE_GENALPHA_H
