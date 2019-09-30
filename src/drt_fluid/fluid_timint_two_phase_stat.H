/*-----------------------------------------------------------*/
/*! \file

\brief stationary two-phase flow problems

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_TWO_PHASE_STAT_H
#define FLD_TIMINT_TWO_PHASE_STAT_H


#include "fluid_timint_two_phase.H"
#include "fluid_timint_stat.H"

#include "../linalg/linalg_utils.H"


namespace FLD
{
  class TimIntTwoPhaseStat : public TimIntStationary, public TimIntTwoPhase
  {
   public:
    /// Standard Constructor
    TimIntTwoPhaseStat(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntTwoPhaseStat();

    /*!
    \brief initialization

     */
    virtual void Init();


   protected:
   private:
  };  // class TimIntTwoPhaseStat

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_TWO_PHASE_STAT_H
