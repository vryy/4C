/*----------------------------------------------------------------------*/
/*! \file
 \brief

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_ADAPTER_POROMULTIPHASE_HPP
#define FOUR_C_ADAPTER_POROMULTIPHASE_HPP

#include "4C_config.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Adapter
{
  class PoroFluidMultiphaseWrapper;
  class Structure;
}  // namespace Adapter

namespace Core::LinAlg
{
  class Solver;
}

namespace Adapter
{
  class PoroMultiPhase
  {
   public:
    /// constructor
    PoroMultiPhase(){};

    /// virtual destructor
    virtual ~PoroMultiPhase() = default;

    /// initialization
    virtual void Init(const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& algoparams, const Teuchos::ParameterList& structparams,
        const Teuchos::ParameterList& fluidparams, const std::string& struct_disname,
        const std::string& fluid_disname, bool isale, int nds_disp, int nds_vel,
        int nds_solidpressure, int ndsporofluid_scatra,
        const std::map<int, std::set<int>>* nearbyelepairs) = 0;

    /// read restart
    virtual void read_restart(int restart) = 0;

    /// test results (if necessary)
    virtual void CreateFieldTest() = 0;

    /// setup
    virtual void SetupSystem() = 0;

    /// setup the solver (only for monolithic system)
    virtual bool SetupSolver() = 0;

    /// perform relaxation (only for partitioned system)
    virtual void PerformRelaxation(Teuchos::RCP<const Epetra_Vector> phi, const int itnum) = 0;

    /// get relaxed fluid solution (only for partitioned system)
    virtual Teuchos::RCP<const Epetra_Vector> RelaxedFluidPhinp() const = 0;

    /// set relaxed fluid solution on structure (only for partitioned system)
    virtual void set_relaxed_fluid_solution() = 0;

    /// prepare timeloop of coupled problem
    virtual void prepare_time_loop() = 0;

    /// timeloop of coupled problem
    virtual void Timeloop() = 0;

    /// time step of coupled problem
    virtual void TimeStep() = 0;

    /// time step of coupled problem
    virtual void prepare_time_step() = 0;

    //! update time step and print to screen
    virtual void update_and_output() = 0;

    /// set structure solution on scatra field
    virtual void set_struct_solution(
        Teuchos::RCP<const Epetra_Vector> disp, Teuchos::RCP<const Epetra_Vector> vel) = 0;

    /// set scatra solution on fluid field
    virtual void SetScatraSolution(unsigned nds, Teuchos::RCP<const Epetra_Vector> scalars) = 0;

    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> StructDofRowMap() const = 0;

    /// unknown displacements at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> StructDispnp() const = 0;

    /// unknown velocity at \f$t_{n+1}\f$
    virtual Teuchos::RCP<const Epetra_Vector> StructVelnp() const = 0;

    /// dof map of vector of unknowns
    virtual Teuchos::RCP<const Epetra_Map> FluidDofRowMap() const = 0;

    /// dof map of vector of unknowns of artery field
    virtual Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const = 0;

    /// return fluid flux
    virtual Teuchos::RCP<const Epetra_MultiVector> FluidFlux() const = 0;

    /// return fluid solution variable
    virtual Teuchos::RCP<const Epetra_Vector> FluidPhinp() const = 0;

    /// return fluid solution variable
    virtual Teuchos::RCP<const Epetra_Vector> FluidSaturation() const = 0;

    /// return fluid solution variable
    virtual Teuchos::RCP<const Epetra_Vector> FluidPressure() const = 0;

    /// return fluid solution variable
    virtual Teuchos::RCP<const Epetra_Vector> SolidPressure() const = 0;

    //! unique map of all dofs that should be constrained with DBC
    virtual Teuchos::RCP<const Epetra_Map> combined_dbc_map() const = 0;

    //! evaluate all fields at x^n+1 with x^n+1 = x_n + stepinc
    virtual void Evaluate(Teuchos::RCP<const Epetra_Vector> sx,
        Teuchos::RCP<const Epetra_Vector> fx, const bool firstcall) = 0;

    //! access to monolithic right-hand side vector
    virtual Teuchos::RCP<const Epetra_Vector> RHS() const = 0;

    //! update all fields after convergence (add increment on displacements and fluid primary
    //! variables)
    virtual void update_fields_after_convergence(
        Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx) = 0;

    //! get the extractor
    virtual Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> Extractor() const = 0;

    //! get the monolithic system matrix
    virtual Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> BlockSystemMatrix() const = 0;

    //! get structure field
    virtual const Teuchos::RCP<Adapter::Structure>& structure_field() = 0;

    //! get fluid field
    virtual const Teuchos::RCP<Adapter::PoroFluidMultiphaseWrapper>& fluid_field() = 0;

    //! build the block null spaces
    virtual void build_block_null_spaces(Teuchos::RCP<Core::LinAlg::Solver>& solver) = 0;

    //! build the block null spaces
    virtual void build_artery_block_null_space(
        Teuchos::RCP<Core::LinAlg::Solver>& solver, const int& arteryblocknum) = 0;
  };
}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
