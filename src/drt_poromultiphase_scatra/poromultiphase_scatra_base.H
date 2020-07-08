/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_BASE_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_BASE_H_

#include <set>

#include "../drt_adapter/adapter_algorithmbase.H"
#include "../drt_porofluidmultiphase/porofluidmultiphase_utils.H"
#include "../drt_inpar/inpar_poromultiphase_scatra.H"

#include "Teuchos_RCP.hpp"
#include <Epetra_Vector.h>
#include <Epetra_Time.h>

// forward declarations
namespace ADAPTER
{
  class ScaTraBaseAlgorithm;
  class PoroMultiPhase;
}  // namespace ADAPTER

namespace DRT
{
  class Discretization;
}

namespace SCATRA
{
  class MeshtyingStrategyArtery;
}

namespace POROMULTIPHASESCATRA
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseScaTraBase : public ADAPTER::AlgorithmBase
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseScaTraBase(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);  // Problem builder

    /// virtual destructor
    virtual ~PoroMultiPhaseScaTraBase(){};

    /// initialization
    virtual void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& poroparams,
        const Teuchos::ParameterList& structparams, const Teuchos::ParameterList& fluidparams,
        const Teuchos::ParameterList& scatraparams, const std::string& struct_disname,
        const std::string& fluid_disname, const std::string& scatra_disname, bool isale,
        int nds_disp, int nds_vel, int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs) = 0;

    /// read restart
    void ReadRestart(int restart);

    /// create result test for subproblems
    void CreateFieldTest();

    /// setup
    virtual void SetupSystem() = 0;

    /// setup solver (only needed in monolithic case)
    virtual void SetupSolver() = 0;

    /// prepare timeloop of coupled problem
    void PrepareTimeLoop();

    /// timeloop of coupled problem
    void Timeloop();

    /// time step of coupled problem --> here the actual action happens (overwritten by sub-classes)
    virtual void TimeStep() = 0;

    /// time step of coupled problem
    void PrepareTimeStep() { PrepareTimeStep(false); };

    /// time step of coupled problem
    void PrepareTimeStep(bool printheader);

    //! update time step and print to screen
    void UpdateAndOutput();

    //! apply solution of poro-problem to scatra
    void SetPoroSolution();

    //! apply solution of scatra to poro
    void SetScatraSolution();

    //! apply the additional Dirichlet boundary condition for volume fraction species
    void ApplyAdditionalDBCForVolFracSpecies();

    //! access to poro field
    const Teuchos::RCP<ADAPTER::PoroMultiPhase>& PoroField() { return poromulti_; }

    //! access to fluid field
    const Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm>& ScatraAlgo() { return scatra_; }

    /// dof map of vector of unknowns of scatra field
    virtual Teuchos::RCP<const Epetra_Map> ScatraDofRowMap() const;

    //! handle divergence of solver
    void HandleDivergence() const;

   private:
    /// underlying poroelast multi phase
    Teuchos::RCP<ADAPTER::PoroMultiPhase> poromulti_;

    /// underlying scatra problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra_;

    // flux-reconstruction method
    INPAR::POROFLUIDMULTIPHASE::FluxReconstructionMethod fluxreconmethod_;

    //! dofset of scatra field on fluid dis
    // TODO: find a better way to do this. Perhaps this should be moved to the adapter?
    int ndsporofluid_scatra_;

    Epetra_Time timertimestep_;  //!< timer for measurement of duration of one time-step
    double dttimestep_;          //!< duration of one time step

   protected:
    //! what to do when nonlinear solution fails
    enum INPAR::POROMULTIPHASESCATRA::DivContAct divcontype_;
    // do we perform coupling with 1D artery
    const bool artery_coupl_;

    // additional volume-fraction species Dirichlet conditions
    Teuchos::RCP<Epetra_Map> add_dirichmaps_volfrac_spec_;

    Teuchos::RCP<SCATRA::MeshtyingStrategyArtery> scatramsht_;

  };  // PoroMultiPhaseScaTraBase


}  // namespace POROMULTIPHASESCATRA



#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_BASE_H_ */
