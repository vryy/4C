/*----------------------------------------------------------------------*/
/*! \file

\brief one-step theta time integration scheme for level-set problems

\level 2

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236

*----------------------------------------------------------------------*/

#ifndef LEVELSET_TIMINT_OST_H
#define LEVELSET_TIMINT_OST_H

#include "../drt_scatra/scatra_timint_ost.H"
#include "levelset_algorithm.H"


namespace SCATRA
{
  class LevelSetTimIntOneStepTheta : public LevelSetAlgorithm, public TimIntOneStepTheta
  {
   public:
    /// standard Constructor
    LevelSetTimIntOneStepTheta(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    /// destructor
    virtual ~LevelSetTimIntOneStepTheta();

    /// initialize time-integration scheme
    virtual void Init();

    /// setup time-integration scheme
    virtual void Setup();

    /// read restart data
    virtual void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null);

    /// redistribute the scatra discretization and vectors according to nodegraph
    void Redistribute(const Teuchos::RCP<Epetra_CrsGraph>& nodegraph);

    /// interpolate phi to intermediate time n+theta with 0<theta<1
    Teuchos::RCP<Epetra_Vector> Phinptheta(const double theta_inter);

    /// interpolate phidt to intermediate time n+theta with 0<theta<1
    Teuchos::RCP<Epetra_Vector> Phidtnptheta(const double theta_inter);

   protected:
    /// Print information about current time step to screen (reimplementation for OST)
    void PrintTimeStepInfo();

    /// calculate consistent initial scalar time derivatives in compliance with initial scalar field
    void CalcInitialTimeDerivative();

    /// additional predictor not intended for level-set methods
    void ExplicitPredictor() const { return; };

    /// Set the part of the righthandside belonging to the last timestep.
    void SetOldPartOfRighthandside();

    /// update state vectors
    /// current solution becomes old solution of next time step
    void UpdateState();

    /// update the solution after Solve()
    /// extended version for coupled level-set problems including reinitialization
    void Update(const int num = 0);

    /// update phi within the reinitialization loop
    void UpdateReinit();

   private:
  };  // class LevelSetTimIntOneStepTheta

}  // namespace SCATRA

#endif  // #ifndef LEVELSET_TIMINT_OST_H
