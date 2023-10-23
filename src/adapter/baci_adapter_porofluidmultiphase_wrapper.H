/*----------------------------------------------------------------------*/
/*! \file
 \brief a wrapper for porous multiphase flow algorithms

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_ADAPTER_POROFLUIDMULTIPHASE_WRAPPER_H
#define BACI_ADAPTER_POROFLUIDMULTIPHASE_WRAPPER_H

#include "baci_adapter_porofluidmultiphase.H"

#include <Teuchos_RCP.hpp>

// forward declarations
namespace DRT
{
  class ResultTest;
  class Discretization;
}  // namespace DRT

namespace POROFLUIDMULTIPHASE
{
  class TimIntImpl;
}

namespace Teuchos
{
  class ParameterList;
}

namespace ADAPTER
{
  /// basic multiphase porous flow adapter
  class PoroFluidMultiphaseWrapper : public PoroFluidMultiphase
  {
   public:
    /// constructor
    explicit PoroFluidMultiphaseWrapper(Teuchos::RCP<PoroFluidMultiphase> porofluid);

    /// initialization
    void Init(const bool isale,         ///< ALE flag
        const int nds_disp,             ///< number of dofset associated with displacements
        const int nds_vel,              ///< number of dofset associated with fluid velocities
        const int nds_solidpressure,    ///< number of dofset associated with solid pressure
        const int ndsporofluid_scatra,  ///< number of dofset associated with scalar on fluid
                                        ///< discretization
        const std::map<int, std::set<int>>*
            nearbyelepairs  ///< possible interaction partners between porofluid and artery
                            ///< discretization
        ) override;

    /// create result test for multiphase porous fluid field
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest() override;

    /// read restart
    void ReadRestart(int restart) override;

    /// access dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds = 0) const override;

    /// access dof row map
    Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const override;

    /// access coupled system matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> ArteryPorofluidSysmat() const override;

    /// direct access to discretization
    Teuchos::RCP<DRT::Discretization> Discretization() const override;

    //! apply moving mesh data
    void ApplyMeshMovement(Teuchos::RCP<const Epetra_Vector> dispnp  //!< displacement vector
        ) override;

    //! set state on discretization
    void SetState(
        unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state) override;

    //! set convective velocity field (+ pressure and acceleration field as
    //! well as fine-scale velocity field, if required)
    void SetVelocityField(Teuchos::RCP<const Epetra_Vector> vel  //!< velocity vector
        ) override;

    //! set solution of scatra problem
    void SetScatraSolution(unsigned nds, Teuchos::RCP<const Epetra_Vector> scalars);

    //! return primary field at time n+1
    Teuchos::RCP<const Epetra_Vector> Phinp() const override;

    //! return primary field at time n
    Teuchos::RCP<const Epetra_Vector> Phin() const override;

    //! return solid pressure field at time n+1
    Teuchos::RCP<const Epetra_Vector> SolidPressure() const override;

    //! return pressure field at time n+1
    Teuchos::RCP<const Epetra_Vector> Pressure() const override;

    //! return saturation field at time n+1
    Teuchos::RCP<const Epetra_Vector> Saturation() const override;

    //! return valid volume fraction species dof vector
    Teuchos::RCP<const Epetra_Vector> ValidVolFracSpecDofs() const override;

    //! return phase flux field at time n+1
    Teuchos::RCP<const Epetra_MultiVector> Flux() const override;

    //! return number of dof set associated with solid pressure
    int GetDofSetNumberOfSolidPressure() const override;

    //! do time integration (time loop)
    void TimeLoop() override;

    //! initialization procedure prior to evaluation of a time step
    void PrepareTimeStep() override;

    //! output solution and restart data to file
    void Output() override;

    //! update the solution after convergence of the nonlinear iteration.
    void Update() override;

    //! calculate error compared to analytical solution
    void EvaluateErrorComparedToAnalyticalSol() override;

    //! general solver call for coupled algorithms
    void Solve() override;

    /// prepare timeloop of coupled problem
    void PrepareTimeLoop() override;

    //! Return MapExtractor for Dirichlet boundary conditions
    const Teuchos::RCP<const CORE::LINALG::MapExtractor> GetDBCMapExtractor() const override;

    //! right-hand side alias the dynamic force residual
    Teuchos::RCP<const Epetra_Vector> RHS() const override;

    //! right-hand side alias the dynamic force residual for coupled system
    Teuchos::RCP<const Epetra_Vector> ArteryPorofluidRHS() const override;

    //! iterative update of phinp
    void UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc) override;

    //! reconstruct pressures and saturation from current solution
    void ReconstructPressuresAndSaturations() override;

    //! reconstruct flux from current solution
    void ReconstructFlux() override;

    //! build linear system tangent matrix, rhs/force residual
    void Evaluate() override;

    // Assemble Off-Diagonal Fluid-Structure Coupling matrix
    void AssembleFluidStructCouplingMat(Teuchos::RCP<CORE::LINALG::SparseOperator> k_fs) override;

    // Assemble Off-Diagonal Fluid-Scatra Coupling matrix
    void AssembleFluidScatraCouplingMat(Teuchos::RCP<CORE::LINALG::SparseOperator> k_pfs) override;

    /// direct access to system matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix() override;

    // return arterial network time integrator
    Teuchos::RCP<ADAPTER::ArtNet> ArtNetTimInt() override;

   private:
    /// multiphase porous flow time integrator
    Teuchos::RCP<PoroFluidMultiphase> porofluid_;

  };  // class PoroFluidMultiphaseWrapper

}  // namespace ADAPTER



#endif  // ADAPTER_POROFLUIDMULTIPHASE_WRAPPER_H
