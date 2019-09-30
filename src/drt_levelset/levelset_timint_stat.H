/*----------------------------------------------------------------------*/
/*! \file

\brief stationary time integration scheme for level-set problems (for coupled problems only)
       just a dummy

\level 2

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

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
    LevelSetTimIntStationary(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// Destructor
    virtual ~LevelSetTimIntStationary();

    /// initialize time-integration scheme
    virtual void Init();

    /// setup time-integration scheme
    virtual void Setup();

    /// read restart data
    virtual void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null)
    {
      dserror("You should not need this function!");
      return;
    };

    /// redistribute the scatra discretization and vectors according to nodegraph
    void Redistribute(const Teuchos::RCP<Epetra_CrsGraph>& nodegraph)
    {
      dserror("You should not need this function!");
      return;
    };

   protected:
    /// update state vectors
    /// current solution becomes old solution of next time step
    void UpdateState()
    {
      dserror("You should not need this function!");
      return;
    };

    /// update the solution after Solve()
    /// extended version for coupled level-set problems including reinitialization
    void Update(const int num = 0)
    {
      dserror("You should not need this function!");
      return;
    };

    /// update phi within the reinitialization loop
    void UpdateReinit()
    {
      dserror("You should not need this function!");
      return;
    };

   private:
  };  // class TimIntStationary

}  // namespace SCATRA

#endif  // #ifndef LEVELSET_TIMINT_STAT_H
