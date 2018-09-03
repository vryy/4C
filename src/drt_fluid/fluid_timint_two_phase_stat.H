/*!----------------------------------------------------------------------
\file fluid_timint_two_phase_stat.H
\brief

\level 2

<pre>
\maintainer Christoph Ager & Martin Kronbichler
             {winter,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>

*----------------------------------------------------------------------*/

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
