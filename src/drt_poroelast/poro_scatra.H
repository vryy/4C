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

#include "poroelast_monolithic.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include <Teuchos_RefCountPtr.hpp>
#include <Epetra_Vector.h>

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{
  class Coupling;
  class FSIStructureWrapper;
}

namespace PORO_SCATRA
{
  class PartPORO_SCATRA
  {
  public:

    //Constructors/destructors

    PartPORO_SCATRA(const Epetra_Comm& comm,
        const Teuchos::ParameterList& timeparams);// Problem builder

    virtual ~PartPORO_SCATRA()
    {
    }
    ; // Problem destructor

    //Methods


    void DoPoroStep();
    void DoScatraStep();
    void IncrementTimeAndStep();
    bool NotFinished()
    {
      return step_ < numstep_ and time_ <= timemax_;
    }
    ; //! flag whether time loop should be finished
    void SetPoroSolution();
    void SetMeshDisp();
    void SetupDiscretizations(const Epetra_Comm& comm);
    void SetVelocityFields();
    void Timeloop(); // Main time loop.

    virtual void SetupSystem();
    //! read and set fields needed for restart
    virtual void ReadRestart();

    virtual void TestResults(const Epetra_Comm& comm);

  protected:
    //Members

    Teuchos::RCP<POROELAST::Monolithic> poroelast_subproblem_; // Pointer to the poroelastic problem. (drt_poroelastic)
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_subproblem_; // Pointer to the ScaTra problem.     (drt_scatra)

    int step_;
    int numstep_;

    double dt_;
    double time_;
    double timemax_;
  };
}

#endif
