#ifndef FOUR_C_ADAPTER_FLD_FLUID_FLUID_FSI_HPP
#define FOUR_C_ADAPTER_FLD_FLUID_FLUID_FSI_HPP

#include "4C_config.hpp"

#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_inpar_xfem.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  class Solver;
  class MapExtractor;
}  // namespace Core::LinAlg

namespace Core::IO
{
  class DiscretizationWriter;
}

namespace FLD
{
  class XFluidFluid;
  namespace Utils
  {
    class MapExtractor;
    class XFluidFluidMapExtractor;
  }  // namespace Utils
}  // namespace FLD

namespace Adapter
{
  class FluidFluidFSI : public FluidFSI
  {
   public:
    /// constructor
    FluidFluidFSI(Teuchos::RCP<Fluid> xfluidfluid, Teuchos::RCP<Fluid> embfluid,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        bool isale, bool dirichletcond);

    /// initialize and prepare maps
    void init() override;

    /// prepare time step
    void prepare_time_step() override;

    /// save results of current time step, do XFEM cut and refresh the
    /// merged fluid map extractor
    void update() override;

    /// solve for pure fluid-fluid-ale problem
    void solve() override;

    Teuchos::RCP<Core::LinAlg::Vector<double>> relaxation_solve(
        Teuchos::RCP<Core::LinAlg::Vector<double>> ivel) override
    {
      FOUR_C_THROW("Do not call RexationSolve for XFFSI.");
      return Teuchos::null;
    }

    Teuchos::RCP<Core::LinAlg::Vector<double>> extract_interface_forces() override
    {
      FOUR_C_THROW("Do not call extract_interface_forces for XFFSI.");
      return Teuchos::null;
    }

    /// @name Accessors
    //@{

    // get merged xfluid-fluid dof row map
    Teuchos::RCP<const Epetra_Map> dof_row_map() override;

    /// communication object at the interface
    Teuchos::RCP<FLD::Utils::MapExtractor> const& interface() const override
    {
      return mergedfluidinterface_;
    }

    Teuchos::RCP<const Core::LinAlg::Vector<double>> grid_vel() override;
    Teuchos::RCP<Core::LinAlg::Vector<double>> write_access_grid_vel();

    Teuchos::RCP<const Core::LinAlg::Vector<double>> dispnp() override;
    Teuchos::RCP<Core::LinAlg::Vector<double>> write_access_dispnp();
    Teuchos::RCP<const Core::LinAlg::Vector<double>> dispn() override;

    /// get the velocity row map of the embedded fluid
    Teuchos::RCP<const Epetra_Map> velocity_row_map() override;

    /// get block system matrix
    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> block_system_matrix() override;

    // access to embedded discretization
    const Teuchos::RCP<Core::FE::Discretization>& discretization() override;

    // return discretization writer of embedded fluid discretization (for special purpose output)
    const Teuchos::RCP<Core::IO::DiscretizationWriter>& disc_writer() override { return output_; }

    /// get map extractor for background/embedded fluid
    Teuchos::RCP<FLD::Utils::XFluidFluidMapExtractor> const& x_fluid_fluid_map_extractor();

    //@}

    /// Apply initial mesh displacement
    void apply_initial_mesh_displacement(
        Teuchos::RCP<const Core::LinAlg::Vector<double>> initfluiddisp) override
    {
      FOUR_C_THROW("Not implemented, yet!");
    }

    // apply ALE-mesh displacements to embedded fluid
    void apply_mesh_displacement(
        Teuchos::RCP<const Core::LinAlg::Vector<double>> fluiddisp) override;

    /// evaluate the fluid and update the merged fluid/FSI DOF-map extractor in case of a change in
    /// the DOF-maps
    void evaluate(Teuchos::RCP<const Core::LinAlg::Vector<double>>
            stepinc  ///< solution increment between time step n and n+1
        ) override;

    /// request fluid system matrix & shapederivatives as blockmatrices when called with true
    /// (indicated monolithic XFFSI with fluidsplit)
    void use_block_matrix(bool split_fluidsysmat = false) override;

    /// determine, whether the ALE-mesh should be relaxed at current time step
    bool is_ale_relaxation_step(int step) const;

    /// get type of monolithic XFFSI approach
    Inpar::XFEM::MonolithicXffsiApproach monolithic_xffsi_approach() const
    {
      return monolithic_approach_;
    }

    Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> shape_derivatives() override;

   private:
    /// setup of map extractor to distinguish between FSI DOF-map and
    /// merged inner embedded fluid and background fluid DOF-map
    void setup_interface(const int nds_master = 0) override;

    /// prepare underlying extended shape derivatives matrix, that is based
    /// on the merged fluid dof-map (with background fluid dof set to zero),
    /// as it may change
    void prepare_shape_derivatives();

    /// type cast pointer to XFluidFluid
    Teuchos::RCP<FLD::XFluidFluid> xfluidfluid_;

    /// fsi map extractor for merged fluid maps (to keep fsi interface-DOF apart from
    /// merged inner DOF (inner embedded fluid together with background fluid)
    Teuchos::RCP<FLD::Utils::MapExtractor> mergedfluidinterface_;

    /// type of monolithic XFluid-Fluid approach (decides whether ALE-mesh is fixed during
    /// Newton iteration)
    enum Inpar::XFEM::MonolithicXffsiApproach monolithic_approach_;

    /// flag, that indicates, whether ALE-relaxation is activated
    bool relaxing_ale_;

    /// no. of timesteps, after which ALE-mesh should be relaxed
    int relaxing_ale_every_;
  };
}  // namespace Adapter

FOUR_C_NAMESPACE_CLOSE

#endif
