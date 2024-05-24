/*----------------------------------------------------------------------*/
/*! \file

\brief scatra time integration for elch

\level 2

 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_ELCH_SCL_HPP
#define FOUR_C_SCATRA_TIMINT_ELCH_SCL_HPP

#include "4C_config.hpp"

#include "4C_coupling_adapter.hpp"
#include "4C_inpar_elch.hpp"
#include "4C_io_runtime_csv_writer.hpp"
#include "4C_scatra_timint_elch.hpp"

FOUR_C_NAMESPACE_OPEN

/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace FLD
{
  class Meshtying;
}

namespace ADAPTER
{
  class Coupling;
  class ScaTraBaseAlgorithm;
}  // namespace ADAPTER

namespace SCATRA
{

  class ScaTraTimIntElchSCL : public virtual ScaTraTimIntElch
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    ScaTraTimIntElchSCL(Teuchos::RCP<DRT::Discretization> dis,
        Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! create result test for micro field
    Teuchos::RCP<CORE::UTILS::ResultTest> create_micro_field_test();

    //! get time integration of micro problem
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> MicroScaTraField();

    void NonlinearSolve() override;

    void check_and_write_output_and_restart() override;

    void prepare_time_loop() override;

    void prepare_time_step() override;

    void read_restart_problem_specific(int step, IO::DiscretizationReader& reader) override;

    void Setup() override;

    void TestResults() override;

    void Update() override;

   protected:
    void add_problem_specific_parameters_and_vectors(Teuchos::ParameterList& params) override;

    void calc_initial_potential_field() override;

    void create_meshtying_strategy() override;

   private:
    //! assemble micro and macro and apply mesh tying between micro and macro model
    void assemble_and_apply_mesh_tying();

    //! stop Netwon loop on convergence and print L2-Norm of increments and residuals
    bool break_newton_loop_and_print_convergence();

    //! copy solution from coupling nodes from macro discretization to micro discretization
    void copy_solution_to_micro_field();

    //! redistribute micro discretization to minimize processor interfaces
    void redistribute_micro_discretization();

    //! scale micro problem with associated area of macro field
    void scale_micro_problem();

    //! setup coupling between micro and macro field
    void setup_coupling();

    //! update increments in micro and macro field
    void update_iter_micro_macro();

    //! write coupled nodes and node coordinates to csv file
    //! \param glob_micro_macro_coupled_node_gids      coupled micro macro nodes of all procs
    //! \param glob_macro_slave_node_master_node_gids  macro slave/master nodes of all procs
    void write_coupling_to_csv(const std::map<int, int>& glob_micro_macro_coupled_node_gids,
        const std::map<int, int>& glob_macro_slave_node_master_node_gids);

    //! the micro problem is split into sub discretiations. This map relates all nodes in the sub
    //! problem (key) to the coupled node of each sub problem (value)
    std::map<int, int> coupled_micro_nodes_;

    //! DBC maps for coupled problem
    Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps_elch_scl_;

    //! block map of coupled ELCH-SCL problem
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> full_block_map_elch_scl_;

    //! map of coupled ELCH-SCL problem
    Teuchos::RCP<const Epetra_Map> full_map_elch_scl_;

    //! increment of coupled ELCH-SCL problem
    Teuchos::RCP<Epetra_Vector> increment_elch_scl_;

    //! map extractor to get the coupled dofs from macro discretization (CondMap) out of all macro
    //! dofs
    Teuchos::RCP<CORE::LINALG::MapExtractor> macro_coupling_dofs_;

    //! coupling adapter between micro (slave) and macro discretization (master).
    Teuchos::RCP<CORE::ADAPTER::Coupling> macro_micro_coupling_adapter_;

    //! map extractor to get micro and macro dofs from global vector
    //! cond. map: micro, other map: macro
    Teuchos::RCP<CORE::LINALG::MapExtractor> macro_micro_dofs_;

    //! type of system matrix of coupled ELCH-SCL problem
    const CORE::LINALG::MatrixType matrixtype_elch_scl_;

    //! map extractor to get the coupled dofs from micro discretization (CondMap) out of all micro
    //! dofs
    Teuchos::RCP<CORE::LINALG::MapExtractor> micro_coupling_dofs_;

    //! time integrator for micro problem
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> micro_timint_;

    //! residual of coupled ELCH-SCL problem
    Teuchos::RCP<Epetra_Vector> residual_elch_scl_;

    //! solver for coupled ELCH-SCL problem
    Teuchos::RCP<CORE::LINALG::Solver> solver_elch_scl_;

    //! system matrix of coupled ELCH-SCL problem
    Teuchos::RCP<CORE::LINALG::SparseOperator> system_matrix_elch_scl_;
  };
}  // namespace SCATRA
FOUR_C_NAMESPACE_CLOSE

#endif