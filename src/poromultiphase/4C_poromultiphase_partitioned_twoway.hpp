/*----------------------------------------------------------------------*/
/*! \file
 \brief two-way coupled partitioned solution algorithm
        for porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_PARTITIONED_TWOWAY_HPP
#define FOUR_C_POROMULTIPHASE_PARTITIONED_TWOWAY_HPP

#include "4C_config.hpp"

#include "4C_inpar_poromultiphase.hpp"
#include "4C_poromultiphase_partitioned.hpp"

FOUR_C_NAMESPACE_OPEN

namespace POROMULTIPHASE
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhasePartitionedTwoWay : public PoroMultiPhasePartitioned
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhasePartitionedTwoWay(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams);

    /// initialization
    void init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& structparams,
        const Teuchos::ParameterList& fluidparams, const std::string& struct_disname,
        const std::string& fluid_disname, bool isale, int nds_disp, int nds_vel,
        int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs) override;

    /// setup
    void SetupSystem() override;

    /// time step of coupled problem
    void TimeStep() override { outer_loop(); };

    /// read restart
    void read_restart(int restart) override;

    // update
    void update_and_output() override;

    // update
    Teuchos::RCP<const Epetra_Vector> RelaxedFluidPhinp() const override { return fluidphinp_; }

   private:
    //! perform iteration loop between fields
    virtual void outer_loop();

    //! perform iteration step of structure field and set the new disp and vel states in the fluid
    //! field
    virtual void do_struct_step();

    //! perform iteration step of scatra field and set the new phi state in the structure field
    virtual void do_fluid_step();

    //! update the current states in every iteration
    //! states are set to the last solutions obtained
    virtual void iter_update_states();

    //! convergence check of outer loop
    virtual bool convergence_check(int itnum);

    //! perform relaxation
    void PerformRelaxation(Teuchos::RCP<const Epetra_Vector> phi, const int itnum) override;

    /// set (relaxed) fluid solution on structure field
    void set_relaxed_fluid_solution() override;

    /// perform aitken
    void aitken_relaxation(double& omega, const int itnum);

    //! pressure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> phiincnp_;
    //! artery pressure increment of the outer loop
    Teuchos::RCP<Epetra_Vector> arterypressincnp_;
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
    //! is artery coupling active
    bool artery_coupling_active_;

    //! relaxation method
    Inpar::POROMULTIPHASE::RelaxationMethods relaxationmethod_;

  };  // PoroMultiPhasePartitioned


}  // namespace POROMULTIPHASE



FOUR_C_NAMESPACE_CLOSE

#endif
