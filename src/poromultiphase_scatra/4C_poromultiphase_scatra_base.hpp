/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_SCATRA_BASE_HPP
#define FOUR_C_POROMULTIPHASE_SCATRA_BASE_HPP

#include "4C_config.hpp"

#include "4C_adapter_algorithmbase.hpp"
#include "4C_inpar_poromultiphase_scatra.hpp"
#include "4C_porofluidmultiphase_utils.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>

#include <set>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Adapter
{
  class ScaTraBaseAlgorithm;
  class PoroMultiPhase;
}  // namespace Adapter

namespace Discret
{
  class Discretization;
}

namespace ScaTra
{
  class MeshtyingStrategyArtery;
}

namespace PoroMultiPhaseScaTra
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseScaTraBase : public Adapter::AlgorithmBase
  {
   public:
    //! create using a Epetra_Comm
    PoroMultiPhaseScaTraBase(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);  // Problem builder

    //! initialization
    virtual void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& poroparams,
        const Teuchos::ParameterList& structparams, const Teuchos::ParameterList& fluidparams,
        const Teuchos::ParameterList& scatraparams, const std::string& struct_disname,
        const std::string& fluid_disname, const std::string& scatra_disname, bool isale,
        int nds_disp, int nds_vel, int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs) = 0;

    //! read restart
    void read_restart(int restart) override;

    //! create result test for subproblems
    void CreateFieldTest();

    //! setup
    virtual void SetupSystem() = 0;

    //! setup solver (only needed in monolithic case)
    virtual void SetupSolver() = 0;

    //! prepare timeloop of coupled problem
    void prepare_time_loop();

    //! timeloop of coupled problem
    void Timeloop();

    //! time step of coupled problem --> here the actual action happens (overwritten by sub-classes)
    virtual void TimeStep() = 0;

    //! time step of coupled problem
    void prepare_time_step() override { prepare_time_step(false); };

    //! time step of coupled problem
    void prepare_time_step(bool printheader);

    //! update time step and print to screen
    void update_and_output();

    //! apply solution of poro-problem to scatra
    void SetPoroSolution();

    //! apply solution of scatra to poro
    void SetScatraSolution();

    //! apply the additional Dirichlet boundary condition for volume fraction species
    void apply_additional_dbc_for_vol_frac_species();

    //! access to poro field
    const Teuchos::RCP<Adapter::PoroMultiPhase>& poro_field() { return poromulti_; }

    //! access to fluid field
    const Teuchos::RCP<Adapter::ScaTraBaseAlgorithm>& ScatraAlgo() { return scatra_; }

    //! dof map of vector of unknowns of scatra field
    virtual Teuchos::RCP<const Epetra_Map> ScatraDofRowMap() const;

    //! handle divergence of solver
    void HandleDivergence() const;

   private:
    //! underlying poroelast multi phase
    Teuchos::RCP<Adapter::PoroMultiPhase> poromulti_;

    //! underlying scatra problem
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra_;

    //! flux-reconstruction method
    Inpar::POROFLUIDMULTIPHASE::FluxReconstructionMethod fluxreconmethod_;

    //! dofset of scatra field on fluid dis
    //! TODO: find a better way to do this. Perhaps this should be moved to the adapter?
    int ndsporofluid_scatra_;

    Teuchos::Time timertimestep_;  //!< timer for measurement of duration of one time-step
    double dttimestep_;            //!< duration of one time step

   protected:
    //! what to do when nonlinear solution fails
    enum Inpar::PoroMultiPhaseScaTra::DivContAct divcontype_;
    //! do we perform coupling with 1D artery
    const bool artery_coupl_;

    //! additional volume-fraction species Dirichlet conditions
    Teuchos::RCP<Epetra_Map> add_dirichmaps_volfrac_spec_;

    Teuchos::RCP<ScaTra::MeshtyingStrategyArtery> scatramsht_;

  };  // PoroMultiPhaseScaTraBase


}  // namespace PoroMultiPhaseScaTra



FOUR_C_NAMESPACE_CLOSE

#endif
