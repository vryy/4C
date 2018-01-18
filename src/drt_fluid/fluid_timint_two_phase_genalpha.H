/*!----------------------------------------------------------------------
\file fluid_timint_two_phase_genalpha.H
\brief

\level 2

<pre>
\maintainer  Christoph Ager & Martin Kronbichler
             {ager,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_TWO_PHASE_GENALPHA_H
#define FLD_TIMINT_TWO_PHASE_GENALPHA_H


#include "fluid_timint_two_phase.H"
#include "fluid_timint_genalpha.H"

#include "../linalg/linalg_utils.H"


namespace FLD
{

class TimIntTwoPhaseGenAlpha : public TimIntGenAlpha, public TimIntTwoPhase
{
public:
  /// Standard Constructor
  TimIntTwoPhaseGenAlpha(const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

  /// Destructor
  virtual ~TimIntTwoPhaseGenAlpha();

  /*!
  \brief initialization

  */
  virtual void Init();


protected:

private:


}; // class TimIntTwoPhaseGenAlpha

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_TWO_PHASE_GENALPHA_H
