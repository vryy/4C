/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for embedded (ALE-)fluid-fluid problems using XFEM

\level 2


*/
/*----------------------------------------------------------------------*/
#ifndef FOUR_C_ADAPTER_FLD_FLUID_FLUID_FSI_HPP
#define FOUR_C_ADAPTER_FLD_FLUID_FLUID_FSI_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_inpar_xfem.hpp"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class Solver;
  class MapExtractor;
}  // namespace CORE::LINALG

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
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        bool isale, bool dirichletcond);

    /// initialize and prepare maps
    void Init() override;

    /// prepare time step
    void prepare_time_step() override;

    /// save results of current time step, do XFEM cut and refresh the
    /// merged fluid map extractor
    void Update() override;

    /// solve for pure fluid-fluid-ale problem
    void Solve() override;

    Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel) override
    {
      FOUR_C_THROW("Do not call RexationSolve for XFFSI.");
      return Teuchos::null;
    }

    Teuchos::RCP<Epetra_Vector> extract_interface_forces() override
    {
      FOUR_C_THROW("Do not call extract_interface_forces for XFFSI.");
      return Teuchos::null;
    }

    /// @name Accessors
    //@{

    // get merged xfluid-fluid dof row map
    Teuchos::RCP<const Epetra_Map> dof_row_map() override;

    /// communication object at the interface
    Teuchos::RCP<FLD::UTILS::MapExtractor> const& Interface() const override
    {
      return mergedfluidinterface_;
    }

    Teuchos::RCP<const Epetra_Vector> GridVel() override;
    Teuchos::RCP<Epetra_Vector> WriteAccessGridVel();

    Teuchos::RCP<const Epetra_Vector> Dispnp() override;
    Teuchos::RCP<Epetra_Vector> WriteAccessDispnp();
    Teuchos::RCP<const Epetra_Vector> Dispn() override;

    /// get the velocity row map of the embedded fluid
    Teuchos::RCP<const Epetra_Map> VelocityRowMap() override;

    /// get block system matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() override;

    // access to embedded discretization
    const Teuchos::RCP<DRT::Discretization>& Discretization() override;

    // return discretization writer of embedded fluid discretization (for special purpose output)
    const Teuchos::RCP<IO::DiscretizationWriter>& DiscWriter() override { return output_; }

    /// get map extractor for background/embedded fluid
    Teuchos::RCP<FLD::UTILS::XFluidFluidMapExtractor> const& x_fluid_fluid_map_extractor();

    //@}

    /// Apply initial mesh displacement
    void apply_initial_mesh_displacement(Teuchos::RCP<const Epetra_Vector> initfluiddisp) override
    {
      FOUR_C_THROW("Not implemented, yet!");
    }

    // apply ALE-mesh displacements to embedded fluid
    void apply_mesh_displacement(Teuchos::RCP<const Epetra_Vector> fluiddisp) override;

    /// evaluate the fluid and update the merged fluid/FSI DOF-map extractor in case of a change in
    /// the DOF-maps
    void Evaluate(Teuchos::RCP<const Epetra_Vector>
            stepinc  ///< solution increment between time step n and n+1
        ) override;

    /// request fluid system matrix & shapederivatives as blockmatrices when called with true
    /// (indicated monolithic XFFSI with fluidsplit)
    void UseBlockMatrix(bool split_fluidsysmat = false) override;

    /// determine, whether the ALE-mesh should be relaxed at current time step
    bool IsAleRelaxationStep(int step) const;

    /// get type of monolithic XFFSI approach
    INPAR::XFEM::MonolithicXffsiApproach monolithic_xffsi_approach() const
    {
      return monolithic_approach_;
    }

    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> ShapeDerivatives() override;

   private:
    /// setup of map extractor to distinguish between FSI DOF-map and
    /// merged inner embedded fluid and background fluid DOF-map
    void SetupInterface(const int nds_master = 0) override;

    /// prepare underlying extended shape derivatives matrix, that is based
    /// on the merged fluid dof-map (with background fluid dof set to zero),
    /// as it may change
    void prepare_shape_derivatives();

    /// type cast pointer to XFluidFluid
    Teuchos::RCP<FLD::XFluidFluid> xfluidfluid_;

    /// fsi map extractor for merged fluid maps (to keep fsi interface-DOF apart from
    /// merged inner DOF (inner embedded fluid together with background fluid)
    Teuchos::RCP<FLD::UTILS::MapExtractor> mergedfluidinterface_;

    /// type of monolithic XFluid-Fluid approach (decides whether ALE-mesh is fixed during
    /// Newton iteration)
    enum INPAR::XFEM::MonolithicXffsiApproach monolithic_approach_;

    /// flag, that indicates, whether ALE-relaxation is activated
    bool relaxing_ale_;

    /// no. of timesteps, after which ALE-mesh should be relaxed
    int relaxing_ale_every_;
  };
}  // namespace ADAPTER

FOUR_C_NAMESPACE_CLOSE

#endif
