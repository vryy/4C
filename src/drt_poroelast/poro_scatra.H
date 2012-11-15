/*!----------------------------------------------------------------------
 \file poro_scatra.H
 \brief

 <pre>
 Maintainer: Anh-Tu Vuong
 vuong@lnm.mw.tum.de
 http://www.lnm.mw.tum.de
 089 - 289-15264
 </pre>

 *----------------------------------------------------------------------*/

#ifndef PORO_SCATRA_H
#define PORO_SCATRA_H

/*----------------------------------------------------------------------*
 | header inclusions                                                     |
 *----------------------------------------------------------------------*/
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------*
 | forward declarations                                                  |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
}

namespace POROELAST
{
  class PoroBase;
}

/*----------------------------------------------------------------------*
 |                                                                       |
 *----------------------------------------------------------------------*/
namespace PORO_SCATRA
{
  /// partitioned algorithm for scalar transport in porous media
  class PartPORO_SCATRA
  {
  public:

    ///Constructor
    explicit PartPORO_SCATRA(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);

    ///destructor
    virtual ~PartPORO_SCATRA()
    {
    }
    ;

    //Methods

    //! solve one time step of porous media problem
    void DoPoroStep();
    //! solve one time step of scalar transport problem
    void DoScatraStep();
    //! update current time and time step
    void IncrementTimeAndStep();

    //! flag whether time loop should be finished
    bool NotFinished()
    {
      return step_ < numstep_ and time_ <= timemax_;
    }
    ;

    //! apply solution of poro-problem to scatra
    void SetPoroSolution();
    //! apply displacement fields to scatra
    void SetMeshDisp();
    //! setup discretization, includes cloning the structure discretization
    void SetupDiscretizations(const Epetra_Comm& comm);
    //! apply velocity fields to scatra
    void SetVelocityFields();
    //! Main time loop.
    void Timeloop();

    //! setup for single fields
    virtual void SetupSystem();
    //! read and set fields needed for restart
    virtual void ReadRestart();

    //! perform result test
    virtual void TestResults(const Epetra_Comm& comm);

  protected:
    //Members

    //! Pointer to the poroelastic problem. (drt_poroelastic)
    Teuchos::RCP<POROELAST::PoroBase> poroelast_subproblem_;
    //! Pointer to the ScaTra problem.     (drt_scatra)
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_subproblem_;

    //! current time step number
    int step_;
    //! total number of time steps
    int numstep_;

    //! time step length
    double dt_;
    //! current time
    double time_;
    //! total simulation time
    double timemax_;
  };
}

#endif
