/*----------------------------------------------------------------------*/
/*! \file
 \brief two-way coupled partitioned solution algorithm
        for porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_PARTITIONED_TWOWAY_H_
#define SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_PARTITIONED_TWOWAY_H_

#include "poromultiphase_partitioned.H"

#include "../drt_inpar/inpar_poromultiphase.H"

namespace POROMULTIPHASE
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhasePartitionedTwoWay : public PoroMultiPhasePartitioned
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhasePartitionedTwoWay(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    /// virtual destructor
    virtual ~PoroMultiPhasePartitionedTwoWay(){};

    /// initialization
    virtual void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& structparams,
        const Teuchos::ParameterList& fluidparams, const std::string& struct_disname,
        const std::string& fluid_disname, bool isale, int nds_disp, int nds_vel,
        int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs);

    /// setup
    virtual void SetupSystem();

    /// time step of coupled problem
    virtual void TimeStep() { OuterLoop(); };

    /// read restart
    virtual void ReadRestart(int restart);

    // update
    virtual void UpdateAndOutput();

    // update
    virtual Teuchos::RCP<const Epetra_Vector> RelaxedFluidPhinp() const { return fluidphinp_; }

   private:
    //! perform iteration loop between fields
    virtual void OuterLoop();

    //! perform iteration step of structure field and set the new disp and vel states in the fluid
    //! field
    virtual void DoStructStep();

    //! perform iteration step of scatra field and set the new phi state in the structure field
    virtual void DoFluidStep();

    //! update the current states in every iteration
    //! states are set to the last solutions obtained
    virtual void IterUpdateStates();

    //! convergence check of outer loop
    virtual bool ConvergenceCheck(int itnum);

    //! perform relaxation
    virtual void PerformRelaxation(Teuchos::RCP<const Epetra_Vector> phi, const int itnum);

    /// set (relaxed) fluid solution on structure field
    void SetRelaxedFluidSolution();

    /// perform aitken
    void AitkenRelaxation(double& omega, const int itnum);

    //! pressure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> phiincnp_;
    //! displacement increment of the outer loop
    Teuchos::RCP<Epetra_Vector> dispincnp_;

    //! fluid primary variable at time n+1, iteration i+1
    Teuchos::RCP<Epetra_Vector> fluidphinp_;
    //! fluid primary variable at time n+1, iteration i
    Teuchos::RCP<Epetra_Vector> fluidphioldnp_;
    //! fluid primary variable increment: phi,n+1^i+1 - phi,n+1^i
    Teuchos::RCP<Epetra_Vector> fluidphiincnp_;
    //! old fluid primary variablee increment: phi,n+1^i+1 - phi,n+1^i
    Teuchos::RCP<Epetra_Vector> fluidphiincnpold_;

    //! convergence tolerance
    double ittol_;
    //! relaxation parameter
    double omega_;
    //! start value for relaxation parameter (or constant value)
    double startomega_;
    //! minimum value for relaxation parameter
    double omegamin_;
    //! maximum value for relaxation parameter
    double omegamax_;
    //! maximum iteration steps
    int itmax_;
    //! current iteration step
    int itnum_;
    //! write restart every n steps
    int writerestartevery_;

    //! relaxation method
    INPAR::POROMULTIPHASE::RelaxationMethods relaxationmethod_;

  };  // PoroMultiPhasePartitioned


}  // namespace POROMULTIPHASE



#endif /* SRC_DRT_POROMULTIPHASE_POROMULTIPHASE_PARTITIONED_TWOWAY_H_ */
