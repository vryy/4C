/*!----------------------------------------------------------------------
\file fluid_timint_topopt_stat.H
\brief

\level 3

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235/-235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_TOPOPT_STAT_H
#define FLD_TIMINT_TOPOPT_STAT_H


#include "fluid_timint_topopt.H"
#include "fluid_timint_stat.H"

#include "../linalg/linalg_utils.H"


namespace FLD
{

class TimIntTopOptStat : public TimIntStationary, public TimIntTopOpt
{
public:
  /// Standard Constructor
  TimIntTopOptStat(const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

  /// Destructor
  virtual ~TimIntTopOptStat();

  /*!
  \brief initialization

  */
  virtual void Init();


protected:

private:


}; // class TimIntTopOptStat

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_TOPOPT_STAT_H
