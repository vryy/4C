/*----------------------------------------------------------------------*/
/*! \file
 \brief base class for all porous multiphase flow through elastic medium problems

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROMULTIPHASE_BASE_HPP
#define FOUR_C_POROMULTIPHASE_BASE_HPP


#include "4C_config.hpp"

#include "4C_adapter_algorithmbase.hpp"
#include "4C_adapter_poromultiphase.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

namespace ADAPTER
{
  class PoroFluidMultiphaseWrapper;
  class Structure;
}  // namespace ADAPTER

namespace DRT
{
  class Discretization;
}

namespace POROMULTIPHASE
{
  //! Base class of all solid-scatra algorithms
  class PoroMultiPhaseBase : public ADAPTER::AlgorithmBase, public ADAPTER::PoroMultiPhase
  {
   public:
    /// create using a Epetra_Comm
    PoroMultiPhaseBase(const Epetra_Comm& comm,
        const Teuchos::ParameterList& globaltimeparams);  // Problem builder

    /// initialization
    void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& structparams,
        const Teuchos::ParameterList& fluidparams, const std::string& struct_disname,
        const std::string& fluid_disname, bool isale, int nds_disp, int nds_vel,
        int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs) override = 0;

    /// read restart
    void read_restart(int restart) override;

    /// test results (if necessary)
    void CreateFieldTest() override;

    /// setup
    void SetupSystem() override = 0;

    /// prepare timeloop of coupled problem
    void prepare_time_loop() override;

    /// timeloop of coupled problem
    void Timeloop() override;

    /// time step of coupled problem
    void TimeStep() override = 0;

    /// prepare time step of coupled problem
    void prepare_time_step() override;

    //! update fields after convergence
    void UpdateAndOutput() override;

    /// dof map of vector of unknowns of structure field
    Teuchos::RCP<const Epetra_Map> StructDofRowMap() const override;

    /// dof map of vector of unknowns of fluid field
    Teuchos::RCP<const Epetra_Map> FluidDofRowMap() const override;

    /// dof map of vector of unknowns of artery field
    Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const override;

    /// system matrix of coupled artery porofluid problem
    virtual Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> artery_porofluid_sysmat() const;

    //! access to structural field
    const Teuchos::RCP<ADAPTER::Structure>& StructureField() override { return structure_; }

    //! access to fluid field
    const Teuchos::RCP<ADAPTER::PoroFluidMultiphaseWrapper>& fluid_field() override
    {
      return fluid_;
    }

    /// set structure solution on scatra field
    void SetStructSolution(
        Teuchos::RCP<const Epetra_Vector> disp, Teuchos::RCP<const Epetra_Vector> vel) override;

    /// set scatra solution on fluid field
    void SetScatraSolution(unsigned nds, Teuchos::RCP<const Epetra_Vector> scalars) override;

    //! setup solver (for monolithic only)
    bool SetupSolver() override { return false; };

    /// unknown displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> StructDispnp() const override;

    /// unknown velocity at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> StructVelnp() const override;

    /// return fluid flux
    Teuchos::RCP<const Epetra_MultiVector> FluidFlux() const override;

    /// return fluid solution variable
    Teuchos::RCP<const Epetra_Vector> FluidPhinp() const override;

    /// return relaxed fluid solution variable (partitioned coupling will overwrite this method)
    Teuchos::RCP<const Epetra_Vector> RelaxedFluidPhinp() const override { return FluidPhinp(); };

    /// set (relaxed) fluid solution on structure field (partitioned coupling will overwrite this
    /// method)
    void set_relaxed_fluid_solution() override
    {
      FOUR_C_THROW("set_relaxed_fluid_solution() only available for partitioned schemes!");
      return;
    };

    /// return fluid solution variable
    Teuchos::RCP<const Epetra_Vector> FluidSaturation() const override;

    /// return fluid solution variable
    Teuchos::RCP<const Epetra_Vector> FluidPressure() const override;

    /// return fluid solution variable
    Teuchos::RCP<const Epetra_Vector> SolidPressure() const override;

    //! unique map of all dofs that should be constrained with DBC
    Teuchos::RCP<const Epetra_Map> combined_dbc_map() const override
    {
      FOUR_C_THROW("combined_dbc_map() only available for monolithic schemes!");
      return Teuchos::null;
    };

    //! build the block null spaces
    void build_block_null_spaces(Teuchos::RCP<CORE::LINALG::Solver>& solver) override
    {
      FOUR_C_THROW("build_block_null_spaces() only available for monolithic schemes!");
      return;
    };

    //! build the block null spaces
    void build_artery_block_null_space(
        Teuchos::RCP<CORE::LINALG::Solver>& solver, const int& arteryblocknum) override
    {
      FOUR_C_THROW("build_artery_block_null_space() only available for monolithic schemes!");
      return;
    };

    //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    void Evaluate(Teuchos::RCP<const Epetra_Vector> sx, Teuchos::RCP<const Epetra_Vector> fx,
        const bool firstcall) override
    {
      FOUR_C_THROW("Evaluate() only available for monolithic schemes!");
      return;
    };

    //! update all fields after convergence (add increment on displacements and fluid primary
    //! variables)
    void update_fields_after_convergence(
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx) override
    {
      FOUR_C_THROW("update_fields_after_convergence() only available for monolithic schemes!");
      return;
    };

    /// perform relaxaton (only for partitioned schemes)
    void PerformRelaxation(Teuchos::RCP<const Epetra_Vector> phi, const int itnum) override
    {
      FOUR_C_THROW("PerformRelaxation() only available for partitioned schemes!");
      return;
    };

    //! get monolithic rhs vector
    Teuchos::RCP<const Epetra_Vector> RHS() const override
    {
      FOUR_C_THROW("RHS() only available for monolithic schemes!");
      return Teuchos::null;
    };

    //! get extractor
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> Extractor() const override
    {
      FOUR_C_THROW("Extractor() only available for monolithic schemes!");
      return Teuchos::null;
    };

    //! get monolithic block system matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix() const override
    {
      FOUR_C_THROW("BlockSystemMatrix() only available for monolithic schemes!");
      return Teuchos::null;
    };

   private:
    /// set structure mesh displacement on fluid field
    void set_mesh_disp(Teuchos::RCP<const Epetra_Vector> disp);

    /// set structure velocity field on fluid field
    void set_velocity_fields(Teuchos::RCP<const Epetra_Vector> vel);

    /// underlying structure of the PoroMultiPhase problem
    Teuchos::RCP<ADAPTER::Structure> structure_;

    /// underlying fluid problem of the PoroMultiPhase problem
    Teuchos::RCP<ADAPTER::PoroFluidMultiphaseWrapper> fluid_;

   protected:
    /// a zero vector of full length of structure dofs
    Teuchos::RCP<Epetra_Vector> struct_zeros_;
    //! here the computation of the structure can be skipped, this is helpful if only fluid-scatra
    //! coupling should be calculated
    bool solve_structure_;

    /// coupling with 1D artery network
    const bool artery_coupl_;

    /// Print user output that structure field is disabled
    void print_structure_disabled_info();

  };  // PoroMultiPhaseBase


}  // namespace POROMULTIPHASE



FOUR_C_NAMESPACE_CLOSE

#endif
