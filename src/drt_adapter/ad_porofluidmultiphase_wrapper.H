/*----------------------------------------------------------------------*/
/*!
 \file ad_porofluidmultiphase_wrapper.H

 \brief a wrapper for porous multiphase flow algorithms

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_ADAPTER_AD_POROFLUIDMULTIPHASE_WRAPPER_H_
#define SRC_DRT_ADAPTER_AD_POROFLUIDMULTIPHASE_WRAPPER_H_

#include "ad_porofluidmultiphase.H"

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

    /// virtual destructor to support polymorph destruction
    virtual ~PoroFluidMultiphaseWrapper(){};

    /// initialization
    void Init(const bool isale,        ///< ALE flag
        const int nds_disp,            ///< number of dofset associated with displacements
        const int nds_vel,             ///< number of dofset associated with fluid velocities
        const int nds_solidpressure,   ///< number of dofset associated with solid pressure
        const int ndsporofluid_scatra  ///< number of dofset associated with scalar on fluid
                                       ///< discretization
    );

    /// create result test for multiphase porous fluid field
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    /// read restart
    void ReadRestart(int restart);

    /// access dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds = 0) const;

    /// access dof row map
    Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const;

    /// access coupled system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> ArteryPorofluidSysmat() const;

    /// direct access to discretization
    Teuchos::RCP<DRT::Discretization> Discretization() const;

    //! apply moving mesh data
    void ApplyMeshMovement(Teuchos::RCP<const Epetra_Vector> dispnp  //!< displacement vector
    );

    //! set state on discretization
    void SetState(unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state);

    //! set convective velocity field (+ pressure and acceleration field as
    //! well as fine-scale velocity field, if required)
    virtual void SetVelocityField(Teuchos::RCP<const Epetra_Vector> vel  //!< velocity vector
    );

    //! set solution of scatra problem
    void SetScatraSolution(unsigned nds, Teuchos::RCP<const Epetra_Vector> scalars);

    //! return primary field at time n+1
    Teuchos::RCP<const Epetra_Vector> Phinp() const;

    //! return solid pressure field at time n+1
    Teuchos::RCP<const Epetra_Vector> SolidPressure() const;

    //! return pressure field at time n+1
    Teuchos::RCP<const Epetra_Vector> Pressure() const;

    //! return saturation field at time n+1
    Teuchos::RCP<const Epetra_Vector> Saturation() const;

    //! return saturation field at time n+1
    Teuchos::RCP<const Epetra_Vector> ValidVolFracPressDofs() const;

    //! return phase flux field at time n+1
    Teuchos::RCP<const Epetra_MultiVector> Flux() const;

    //! return number of dof set associated with solid pressure
    virtual int GetDofSetNumberOfSolidPressure() const;

    //! do time integration (time loop)
    virtual void TimeLoop();

    //! initialization procedure prior to evaluation of a time step
    virtual void PrepareTimeStep();

    //! output solution and restart data to file
    virtual void Output();

    //! update the solution after convergence of the nonlinear iteration.
    virtual void Update();

    //! calculate error compared to analytical solution
    virtual void EvaluateErrorComparedToAnalyticalSol();

    //! general solver call for coupled algorithms
    virtual void Solve();

    /// prepare timeloop of coupled problem
    virtual void PrepareTimeLoop();

    //! Return MapExtractor for Dirichlet boundary conditions
    const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() const;

    //! right-hand side alias the dynamic force residual
    Teuchos::RCP<const Epetra_Vector> RHS() const;

    //! right-hand side alias the dynamic force residual for coupled system
    Teuchos::RCP<const Epetra_Vector> ArteryPorofluidRHS() const;

    //! iterative update of phinp
    virtual void UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc);

    //! reconstruct pressures and saturation from current solution
    virtual void ReconstructPressuresAndSaturations();

    //! reconstruct flux from current solution
    virtual void ReconstructFlux();

    //! build linear system tangent matrix, rhs/force residual
    virtual void Evaluate();

    // Assemble Off-Diagonal Fluid-Structure Coupling matrix
    virtual void AssembleFluidStructCouplingMat(Teuchos::RCP<LINALG::SparseOperator> k_fs);

    // Assemble Off-Diagonal Fluid-Scatra Coupling matrix
    virtual void AssembleFluidScatraCouplingMat(Teuchos::RCP<LINALG::SparseOperator> k_pfs);

    /// direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    // return arterial network time integrator
    virtual Teuchos::RCP<ADAPTER::ArtNet> ArtNetTimInt();

   private:
    /// multiphase porous flow time integrator
    Teuchos::RCP<PoroFluidMultiphase> porofluid_;

  };  // class PoroFluidMultiphaseWrapper

}  // namespace ADAPTER



#endif /* SRC_DRT_ADAPTER_AD_POROFLUIDMULTIPHASE_WRAPPER_H_ */
