/*!----------------------------------------------------------------------
\file fluid_timint_loma_genalpha.H
\brief

<pre>
Maintainers: Ursula Rasthofer & Martin Kronbichler
             {rasthofer,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_LOMA_GENALPHA_H
#define FLD_TIMINT_LOMA_GENALPHA_H


#include "fluid_timint_loma.H"
#include "fluid_timint_genalpha.H"

#include "../linalg/linalg_utils.H"


namespace FLD
{

class TimIntLomaGenAlpha : public TimIntGenAlpha, public TimIntLoma
{
public:
  /// Standard Constructor
  TimIntLomaGenAlpha(const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

  /// Destructor
  virtual ~TimIntLomaGenAlpha();




protected:

private:


}; // class TimIntLomaGenAlpha

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_LOMA_GENALPHA_H
