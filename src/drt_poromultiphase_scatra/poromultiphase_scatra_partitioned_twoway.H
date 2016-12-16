/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_partitioned_twoway.H

 \brief two-way coupled partitioned algorithm for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_PARTITIONED_TWOWAY_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_PARTITIONED_TWOWAY_H_


#include "poromultiphase_scatra_partitioned.H"

namespace POROMULTIPHASESCATRA
{

  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseScaTraPartitionedTwoWay : public PoroMultiPhaseScaTraPartitioned
  {
  public:

    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraPartitionedTwoWay(
        const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);

    /// virtual destructor
    virtual ~PoroMultiPhaseScaTraPartitionedTwoWay(){};

    /// initialization
    virtual void Init(
        const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams,
        const Teuchos::ParameterList& poroparams,
        const Teuchos::ParameterList& structparams,
        const Teuchos::ParameterList& fluidparams,
        const Teuchos::ParameterList& scatraparams,
        const std::string& struct_disname,
        const std::string& fluid_disname,
        const std::string& scatra_disname,
        bool isale,
        int nds_disp,
        int nds_vel,
        int nds_solidpressure,
        int ndsporofluid_scatra);

    /// setup
    virtual void SetupSystem();

    /// prepare timeloop of coupled problem
    virtual void PrepareTimeLoop();

    /// timeloop of coupled problem
    virtual void Timeloop();

    /// time step of coupled problem
    virtual void TimeStep() {return;};

    /// time step of coupled problem
    virtual void PrepareTimeStep() {PrepareTimeStep(false);};

    /// time step of coupled problem
    virtual void PrepareTimeStep(bool printheader);

    //! update time step and print to screen
    virtual void UpdateAndOutput();

    //! perform iteration loop between fields
    void Solve();


  private:

    //! perform iteration step of structure field
    void DoPoroStep();

    //! perform iteration step of scatra field
    void DoScatraStep();

    //! convergence check of outer loop
    bool ConvergenceCheck(int itnum);

    //! scalar increment of the outer loop
    Teuchos::RCP<Epetra_Vector> scaincnp_;
    //! structure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> structincnp_;
    //! fluid increment of the outer loop
    Teuchos::RCP<Epetra_Vector> fluidincnp_;

    //! maximum iteration steps
    int itmax_;
    //! convergence tolerance
    double ittol_;


  }; //PoroMultiPhaseScatraPartitioned


}



#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_PARTITIONED_TWOWAY_H_ */
