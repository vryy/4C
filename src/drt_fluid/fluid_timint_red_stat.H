/*!----------------------------------------------------------------------
\file fluid_timint_red_stat.H
\brief

\level 2

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_RED_STAT_H
#define FLD_TIMINT_RED_STAT_H


#include "fluid_timint_red.H"
#include "fluid_timint_stat.H"

#include "../linalg/linalg_utils.H"


namespace FLD
{

class TimIntRedModelsStat : public TimIntStationary, public TimIntRedModels
{
public:
  /// Standard Constructor
  TimIntRedModelsStat(const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

  /// Destructor
  virtual ~TimIntRedModelsStat();

  /*!
  \brief initialization

  */
  virtual void Init();

  /*!
  \brief read restart data

  */
  virtual void ReadRestart(int step);


protected:

private:


}; // class TimIntRedModelsStat

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_RED_STAT_H
