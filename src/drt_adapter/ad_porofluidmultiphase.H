/*----------------------------------------------------------------------*/
/*!
 \file ad_porofluidmultiphase.H

 \brief adapter for multiphase porous flow

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_ADAPTER_AD_POROFLUIDMULTIPHASE_H_
#define SRC_DRT_ADAPTER_AD_POROFLUIDMULTIPHASE_H_


#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include "../linalg/linalg_utils.H"

// forward declarations
namespace DRT
{
  class ResultTest;
}

namespace POROFLUIDMULTIPHASE
{
  class TimIntImpl;
}

namespace Teuchos
{
  class ParameterList;
}

// Algorithm adapter classes for use in multi-field problems
namespace ADAPTER {

  /// basic multiphase porous flow adapter
  class PoroFluidMultiphase
  {
  public:

    /// constructor
    PoroFluidMultiphase(){};

    /// virtual destructor to support polymorph destruction
    virtual ~PoroFluidMultiphase(){};

    /// initialization
    virtual void Init(
        const bool                      isale,                ///< ALE flag
        const int                       nds_disp,             ///< number of dofset associated with displacements
        const int                       nds_vel,              ///< number of dofset associated with fluid velocities
        const int                       nds_solidpressure,    ///< number of dofset associated with solid pressure
        const int                       ndsporofluid_scatra   ///< number of dofset associated with scalar on fluid discretization
        )=0;

    /// create result test for multiphase porous fluid field
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest()=0;

    /// read restart
    virtual void ReadRestart(int restart)=0;

    /// access dof row map
    virtual Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds=0) const=0;

    //! apply moving mesh data
    virtual void ApplyMeshMovement(
      Teuchos::RCP<const Epetra_Vector> dispnp //!< displacement vector
    )=0;

    //! set convective velocity field (+ pressure and acceleration field as
    //! well as fine-scale velocity field, if required)
    virtual void SetVelocityField(
      Teuchos::RCP<const Epetra_Vector>   vel                       //!< velocity vector
     )=0;

    //! set state on discretization
    virtual void SetState(
       unsigned nds,
       const std::string& name,
       Teuchos::RCP<const Epetra_Vector> state
    )=0;

    //! return primary field at time n+1
    virtual Teuchos::RCP<const Epetra_Vector> Phinp() const=0;

    //! return solid pressure field at time n+1
    virtual Teuchos::RCP<const Epetra_Vector> SolidPressure() const=0;

    //! return pressure field at time n+1
    virtual Teuchos::RCP<const Epetra_Vector> Pressure() const=0;

    //! return saturation field at time n+1
    virtual Teuchos::RCP<const Epetra_Vector> Saturation() const=0;

    //! return phase flux field at time n+1
    virtual Teuchos::RCP<const Epetra_MultiVector> Flux() const=0;

    //! do time integration (time loop)
    virtual void TimeLoop()=0;

    //! initialization procedure prior to evaluation of a time step
    virtual void PrepareTimeStep()=0;

    //! output solution and restart data to file
    virtual void Output()=0;

    //! update the solution after convergence of the nonlinear iteration.
    virtual void Update()=0;

    //! calculate error compared to analytical solution
    virtual void EvaluateErrorComparedToAnalyticalSol()=0;

    //! general solver call for coupled algorithms
    virtual void Solve()=0;

    /// prepare timeloop of coupled problem
    virtual void PrepareTimeLoop()=0;

    //! return number of dof set associated with solid pressure
    virtual int GetDofSetNumberOfSolidPressure() const =0;

    //! Return MapExtractor for Dirichlet boundary conditions
    virtual const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() const=0;

    //! right-hand side alias the dynamic force residual
    virtual Teuchos::RCP<const Epetra_Vector> RHS() const=0;

    //! iterative update of phinp
    virtual void UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc)=0;

    //! reconstruct pressures and saturation from current solution
    virtual void ReconstructPressuresAndSaturations()=0;

    //! reconstruct flux from current solution
    virtual void ReconstructFlux()=0;

    //! build linear system tangent matrix, rhs/force residual
    virtual void Evaluate()=0;

    // Assemble Off-Diagonal Fluid-Structure Coupling matrix
    virtual void AssembleFluidCouplingMat(Teuchos::RCP< LINALG::SparseOperator> k_fs)=0;

    //! direct access to system matrix
    virtual Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() = 0;


  }; // class PoroFluidMultiphase

} // namespace ADAPTER


#endif /* SRC_DRT_ADAPTER_AD_POROFLUIDMULTIPHASE_H_ */
