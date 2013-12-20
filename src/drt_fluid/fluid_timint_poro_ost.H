/*!----------------------------------------------------------------------
\file fluid_timint_poro_ost.H
\brief

<pre>
Maintainers: Ursula Rasthofer & Martin Kronbichler
             {rasthofer,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_PORO_OST_H
#define FLD_TIMINT_PORO_OST_H


#include "fluid_timint_poro.H"
#include "fluid_timint_ost.H"

#include "../linalg/linalg_utils.H"


namespace FLD
{

class TimIntPoroOst : public TimIntOneStepTheta, public TimIntPoro
{
public:
  /// Standard Constructor
  TimIntPoroOst(const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

  /// Destructor
  virtual ~TimIntPoroOst();




protected:

private:


}; // class TimIntPoroOst

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_PORO_OST_H
