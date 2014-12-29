/*!----------------------------------------------------------------------
\file fluid_timint_ac_ost.H
\brief

<pre>
   Maintainer: Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.mhpc.mw.tum.de
               089 - 289-10364
</pre>
*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_AC_OST_H
#define FLD_TIMINT_AC_OST_H


#include "fluid_timint_ac.H"
#include "fluid_timint_ost.H"

namespace FLD
{

class TimIntACOst : public TimIntOneStepTheta, public TimIntAC
{
public:
  /// Standard Constructor
  TimIntACOst(const Teuchos::RCP<DRT::Discretization>&    actdis,
      const Teuchos::RCP<LINALG::Solver>&                 solver,
      const Teuchos::RCP<Teuchos::ParameterList>&         params,
      const Teuchos::RCP<IO::DiscretizationWriter>&       output,
      bool                                                alefluid = false);

  /// Destructor
  virtual ~TimIntACOst();

  /*!
  \brief read restart data

  */
  virtual void ReadRestart(int step);

protected:

private:


}; // class TimIntACOst

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_AC_OST_H
