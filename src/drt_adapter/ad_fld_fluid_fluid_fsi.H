/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for embedded (ALE-)fluid-fluid problems using XFEM

\level 2

\maintainer  Martin Kronbichler

*/
/*----------------------------------------------------------------------*/
#ifndef ADAPTER_FLUID_FLUID_FSI_H
#define ADAPTER_FLUID_FLUID_FSI_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

#include "ad_fld_fluid_fsi.H"
#include "../drt_inpar/inpar_xfem.H"

namespace DRT
{
  class Condition;
}

namespace LINALG
{
  class Solver;
  class MapExtractor;
}  // namespace LINALG

namespace IO
{
  class DiscretizationWriter;
}

namespace FLD
{
  class XFluidFluid;
  namespace UTILS
  {
    class MapExtractor;
    class XFluidFluidMapExtractor;
  }  // namespace UTILS
}  // namespace FLD

namespace ADAPTER
{
  class FluidFluidFSI : public FluidFSI
  {
   public:
    /// constructor
    FluidFluidFSI(Teuchos::RCP<Fluid> xfluidfluid, Teuchos::RCP<Fluid> embfluid,
        Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        bool isale, bool dirichletcond);

    /// initialize and prepare maps
    void Init();

    /// prepare time step
    void PrepareTimeStep();

    /// save results of current time step, do XFEM cut and refresh the
    /// merged fluid map extractor
    void Update();

    /// solve for pure fluid-fluid-ale problem
    void Solve();

    Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
    {
      dserror("Do not call RexationSolve for XFFSI.");
      return Teuchos::null;
    }

    Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces()
    {
      dserror("Do not call ExtractInterfaceForces for XFFSI.");
      return Teuchos::null;
    }

    /// @name Accessors
    //@{

    // get merged xfluid-fluid dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap();

    /// communication object at the interface
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& Interface() const
    {
      return mergedfluidinterface_;
    }

    Teuchos::RCP<const Epetra_Vector> GridVel();
    Teuchos::RCP<Epetra_Vector> WriteAccessGridVel();

    Teuchos::RCP<const Epetra_Vector> Dispnp();
    Teuchos::RCP<Epetra_Vector> WriteAccessDispnp();
    Teuchos::RCP<const Epetra_Vector> Dispn();

    /// get the velocity row map of the embedded fluid
    Teuchos::RCP<const Epetra_Map> VelocityRowMap();

    /// get block system matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

    // access to embedded discretization
    const Teuchos::RCP<DRT::Discretization>& Discretization();

    // return discretization writer of embedded fluid discretization (for special purpose output)
    const Teuchos::RCP<IO::DiscretizationWriter>& DiscWriter() { return output_; }

    /// get map extractor for background/embedded fluid
    Teuchos::RCP<FLD::UTILS::XFluidFluidMapExtractor> const& XFluidFluidMapExtractor();

    //@}

    /// Apply initial mesh displacement
    void ApplyInitialMeshDisplacement(Teuchos::RCP<const Epetra_Vector> initfluiddisp)
    {
      dserror("Not implemented, yet!");
    }

    // apply ALE-mesh displacements to embedded fluid
    void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp);

    /// evaluate the fluid and update the merged fluid/FSI DOF-map extractor in case of a change in
    /// the DOF-maps
    void Evaluate(Teuchos::RCP<const Epetra_Vector>
            stepinc  ///< solution increment between time step n and n+1
    );

    /// request fluid system matrix & shapederivatives as blockmatrices when called with true
    /// (indicated monolithic XFFSI with fluidsplit)
    void UseBlockMatrix(bool split_fluidsysmat = false);

    /// determine, whether the ALE-mesh should be relaxed at current time step
    bool IsAleRelaxationStep(int step) const;

    /// get type of monolithic XFFSI approach
    INPAR::XFEM::Monolithic_xffsi_Approach MonolithicXffsiApproach() const
    {
      return monolithic_approach_;
    }

    Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives();

   private:
    /// setup of map extractor to distinguish between FSI DOF-map and
    /// merged inner embedded fluid and background fluid DOF-map
    void SetupInterface(const int nds_master = 0);

    /// prepare underlying extended shape derivatives matrix, that is based
    /// on the merged fluid dof-map (with background fluid dof set to zero),
    /// as it may change
    void PrepareShapeDerivatives();

    /// type cast pointer to XFluidFluid
    Teuchos::RCP<FLD::XFluidFluid> xfluidfluid_;

    /// fsi map extractor for merged fluid maps (to keep fsi interface-DOF apart from
    /// merged inner DOF (inner embedded fluid together with background fluid)
    Teuchos::RCP<FLD::UTILS::MapExtractor> mergedfluidinterface_;

    /// type of monolithic XFluid-Fluid approach (decides whether ALE-mesh is fixed during
    /// Newton iteration)
    enum INPAR::XFEM::Monolithic_xffsi_Approach monolithic_approach_;

    /// flag, that indicates, whether ALE-relaxation is activated
    bool relaxing_ale_;

    /// no. of timesteps, after which ALE-mesh should be relaxed
    int relaxing_ale_every_;
  };
}  // namespace ADAPTER

#endif
