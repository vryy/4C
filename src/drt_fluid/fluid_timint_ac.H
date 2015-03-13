/*!----------------------------------------------------------------------
\file fluid_timint_ac.H
\brief

<pre>
   Maintainer: Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10364
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_AC_H
#define FLD_TIMINT_AC_H


#include "fluidimplicitintegration.H"

namespace FLD
{

class TimIntAC : public virtual FluidImplicitTimeInt
{
public:
  /// Standard Constructor
  TimIntAC(const Teuchos::RCP<DRT::Discretization>&     actdis,
      const Teuchos::RCP<LINALG::Solver>&               solver,
      const Teuchos::RCP<Teuchos::ParameterList>&       params,
      const Teuchos::RCP<IO::DiscretizationWriter>&     output,
      bool                                              alefluid = false);

  /// Destructor
  virtual ~TimIntAC();

  /// read restart from step
  void ReadRestart(int step);

  /// write output
  virtual void Output();

  /// time update of stresses
  virtual void TimeUpdateStresses(){ stressmanager_->GetStresses(trueresidual_,dta_); }; //we do the time update of wss from inside the ac-fsi framework!

protected:

private:


}; // class TimIntAC

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_AC_H
