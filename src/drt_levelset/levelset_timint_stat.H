/*!----------------------------------------------------------------------
\file levelset_timint_stat.H
\brief

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifndef LEVELSET_TIMINT_STAT_H
#define LEVELSET_TIMINT_STAT_H

#include "../drt_scatra/scatra_timint_stat.H"
#include "levelset_algorithm.H"


namespace SCATRA
{

class LevelSetTimIntStationary : public LevelSetAlgorithm, public TimIntStationary
{
public:

  /// Standard Constructor
  LevelSetTimIntStationary(Teuchos::RCP<DRT::Discretization>       dis,
      Teuchos::RCP<LINALG::Solver>            solver,
      Teuchos::RCP<Teuchos::ParameterList>    params,
      Teuchos::RCP<Teuchos::ParameterList>    sctratimintparams,
      Teuchos::RCP<Teuchos::ParameterList>    extraparams,
      Teuchos::RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~LevelSetTimIntStationary();

  /// initialize time-integration scheme
  virtual void Init();

protected:

  /// update phi within the reinitialization loop
  void UpdateReinit(){dserror("You should not need this function!"); return;};

private:

}; // class TimIntStationary

}// namespace SCATRA

#endif  // #ifndef LEVELSET_TIMINT_STAT_H
