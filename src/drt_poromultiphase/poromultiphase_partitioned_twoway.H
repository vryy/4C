/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_partitioned_twoway.H

 \brief two-way coupled solution algorithm
        for porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_PARTITIONED_TWOWAY_H_
#define SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_PARTITIONED_TWOWAY_H_

#include "poromultiphase_partitioned.H"

namespace POROMULTIPHASE
{

  //! Base class of all solid-scatra algorithms
  class PoroMultiPhasePartitionedTwoWay : public PoroMultiPhasePartitioned
  {
  public:

    /// create using a Epetra_Comm
    explicit PoroMultiPhasePartitionedTwoWay(
        const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);

    /// virtual destructor
    virtual ~PoroMultiPhasePartitionedTwoWay(){};

    /// initialization
    virtual void Init(
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams,
        const Teuchos::ParameterList& structparams,
        const Teuchos::ParameterList& fluidparams,
        const std::string& struct_disname,
        const std::string& fluid_disname,
        bool isale,
        int nds_disp,
        int nds_vel,
        int nds_solidpressure);

    /// setup
    virtual void SetupSystem();

    /// timeloop of coupled problem
    virtual void Timeloop();

    /// time step of coupled problem
    virtual void TimeStep() {OuterLoop();};

    /// time step of coupled problem
    virtual void PrepareTimeStep() {PrepareTimeStep(false);};

    /// time step of coupled problem
    virtual void PrepareTimeStep(bool printheader);

    //! update time step and print to screen
    virtual void UpdateAndOutput();

  private:

    //! perform iteration loop between fields
    virtual void OuterLoop();

    /// prepare timeloop of coupled problem
    virtual void PrepareTimeLoop();

    //! perform iteration step of structure field and set the new disp and vel states in the fluid field
    virtual void DoStructStep();

    //! perform iteration step of scatra field and set the new phi state in the structure field
    virtual void DoFluidStep();

    //! update the current states in every iteration
    //! states are set to the last solutions obtained
    virtual void IterUpdateStates();

    //! convergence check of outer loop
    virtual bool ConvergenceCheck(int itnum);

    //! pressure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> presincnp_;
    //! displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnp_;

    //! convergence tolerance
    double ittol_;
    //! maximum iteration steps
    int itmax_;
    //! current iteration step
    int itnum_;

  }; //PoroMultiPhasePartitioned


}



#endif /* SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_PARTITIONED_TWOWAY_H_ */
