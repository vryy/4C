/*----------------------------------------------------------------------*/
/*! \file

\brief stationary time integration scheme for level-set problems (for coupled problems only)
       just a dummy

\level 2


*----------------------------------------------------------------------*/

#ifndef BACI_LEVELSET_TIMINT_STAT_H
#define BACI_LEVELSET_TIMINT_STAT_H

#include "baci_levelset_algorithm.H"
#include "baci_scatra_timint_stat.H"


namespace SCATRA
{
  class LevelSetTimIntStationary : public LevelSetAlgorithm, public TimIntStationary
  {
   public:
    /// Standard Constructor
    LevelSetTimIntStationary(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);


    /// initialize time-integration scheme
    void Init() override;

    /// setup time-integration scheme
    void Setup() override;

    /// read restart data
    void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null) override
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
    void UpdateState() override
    {
      dserror("You should not need this function!");
      return;
    };

    /// update the solution after Solve()
    /// extended version for coupled level-set problems including reinitialization
    void Update(const int num = 0) override
    {
      dserror("You should not need this function!");
      return;
    };

    /// update phi within the reinitialization loop
    void UpdateReinit() override
    {
      dserror("You should not need this function!");
      return;
    };

   private:
  };  // class TimIntStationary

}  // namespace SCATRA

#endif  // LEVELSET_TIMINT_STAT_H
