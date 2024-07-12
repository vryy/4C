/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2


*/
/*--------------------------------------------------------------------------*/

#ifndef FOUR_C_SSI_MONOLITHIC_HPP
#define FOUR_C_SSI_MONOLITHIC_HPP

#include "4C_config.hpp"

#include "4C_ssi_base.hpp"

#include <Teuchos_Time.hpp>

FOUR_C_NAMESPACE_OPEN

// forward declarations
namespace Inpar::SSI
{
  enum class ScaTraTimIntType;
}  // namespace Inpar::SSI

namespace Core::LinAlg
{
  class Solver;
  class Equilibration;
  enum class EquilibrationMethod;
  enum class MatrixType;
}  // namespace Core::LinAlg

namespace SSI
{
  namespace UTILS
  {
    class SSIMaps;
    class SSIMatrices;
    class SSIVectors;
  }  // namespace UTILS

  class AssembleStrategyBase;
  class ContactStrategyBase;
  class DBCHandlerBase;
  class ManifoldMeshTyingStrategyBase;
  class MeshtyingStrategyBase;
  class ScatraStructureOffDiagCoupling;
  class ScaTraManifoldScaTraFluxEvaluator;

  //! equilibration methods applied to system matrix
  struct SSIMonoEquilibrationMethod
  {
    const Core::LinAlg::EquilibrationMethod global;     //! unique equilibration
    const Core::LinAlg::EquilibrationMethod scatra;     //! equilibration for scatra block
    const Core::LinAlg::EquilibrationMethod structure;  //! equilibration for structure block
  };

  enum class Subproblem : int
  {
    scalar_transport,
    structure,
    manifold
  };

  class SsiMono : public SSIBase
  {
   public:
    //! constructor
    explicit SsiMono(const Epetra_Comm& comm,           //!< communicator
        const Teuchos::ParameterList& globaltimeparams  //!< parameter list for time integration
    );

    //! return global map of degrees of freedom
    const Teuchos::RCP<const Epetra_Map>& dof_row_map() const;

    void init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
        const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
        const std::string& struct_disname, const std::string& scatra_disname, bool isAle) override;

    //! return global map extractor (0: scalar transport, 1: structure, [2: scatra manifold])
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> maps_sub_problems() const;

    //! return map extractor associated with all degrees of freedom inside scatra field
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_scatra() const;

    //! return map extractor associated with all degrees of freedom inside scatra manifold field
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_scatra_manifold() const;

    //! return map extractor associated with all degrees of freedom inside structural field
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_structure() const;

    //! return map extractor associated with blocks of global system matrix
    Teuchos::RCP<const Core::LinAlg::MultiMapExtractor> block_map_system_matrix() const;

    //! Return matrix type of global system matrix
    Core::LinAlg::MatrixType matrix_type() const { return matrixtype_; };

    void read_restart(int restart) override;

    void setup() override;

    void setup_system() override;

    /*!
     * @brief solves the linear system
     *
     * @note in case an equilibration method (scaling of rows and columns) is defined this is also
     * performed within this call
     */
    void solve_linear_system();

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps() const { return ssi_maps_; }

    //! return algebraic solver for global system of equations
    const Core::LinAlg::Solver& solver() const { return *solver_; };

    void timeloop() override;

   private:
    //! strategies for Newton-Raphson convergence check
    class ConvCheckStrategyBase;
    class ConvCheckStrategyElch;
    class ConvCheckStrategyElchScaTraManifold;
    class ConvCheckStrategyStd;

    //! apply the contact contributions to matrices and residuals of the sub problems
    void apply_contact_to_sub_problems();

    //! apply the Dirichlet boundary conditions to the ssi system, i.e. matrices and residuals
    void apply_dbc_to_system();

    //! apply mesh tying between manifold domains on matrices and residuals
    void apply_manifold_meshtying();

    //! perform mesh tying on matrices and residuals as obtained from sub problems
    void apply_meshtying_to_sub_problems();

    //! assemble global system of equations
    void assemble_mat_and_rhs();

    //! assemble linearization of scatra residuals to system matrix
    void assemble_mat_scatra();

    //! assemble linearization of scatra on manifold residuals to system matrix
    void assemble_mat_scatra_manifold();

    //! assemble linearization of structural residuals to system matrix
    void assemble_mat_structure();

    //! build null spaces associated with blocks of global system matrix
    void build_null_spaces() const;

    //! calc initial potential field for monolithic SSI problem including scatra and scatra manifold
    //! fields
    void calc_initial_potential_field();

    //! calc initial time derivative of transported scalars for monolithic SSI problem including
    //! scatra and scatra manifold fields
    void calc_initial_time_derivative();

    //! call complete on the sub problem matrices
    void complete_subproblem_matrices();

    //! distribute solution to all other fields
    //! \param restore_velocity   restore velocity when structure_field()->set_state() is called
    void distribute_solution_all_fields(bool restore_velocity = false);

    //! evaluate all off-diagonal matrix contributions
    void evaluate_off_diag_contributions();

    //! Evaluate ScaTra including copy to corresponding ssi matrix
    void evaluate_scatra();

    //! Evaluate ScaTra on manifold incl. coupling with scatra
    void evaluate_scatra_manifold();

    //! get matrix and right-hand-side for all subproblems incl. coupling
    void evaluate_subproblems();

    //! build and return vector of equilibration methods for each block of system matrix
    std::vector<Core::LinAlg::EquilibrationMethod> get_block_equilibration();

    /*!
     * @note This is only necessary in the first iteration of the simulation, since only there the
     * graph of the matrix changes
     *
     * @return flag indicating if we need to uncomplete the matrices before adding the mesh tying
     * contributions.
     */
    bool is_uncomplete_of_matrices_necessary_for_mesh_tying() const;

    void output() override;

    //! do everything, that has to be done once before first time step
    void prepare_time_loop();

    void prepare_time_step() override;

    //! prepare output for subproblems if needed
    void prepare_output();

    //! print system matrix, rhs, and map of system matrix to file
    void print_system_matrix_rhs_to_mat_lab_format();

    //! print time step size, time, and number of time step
    void print_time_step_info();

    //! set scatra manifold solution on scatra field
    void set_scatra_manifold_solution(Teuchos::RCP<const Epetra_Vector> phi);

    //! evaluate time step using Newton-Raphson iteration
    void newton_loop();

    void update() override;

    //! update ScaTra state within Newton iteration
    void update_iter_scatra();

    //! update structure state within Newton iteration
    void update_iter_structure();

    //! Dirichlet boundary condition handler
    Teuchos::RCP<SSI::DBCHandlerBase> dbc_handler_;

    //! time for element evaluation and assembly of global system of equations
    double dt_eval_ = 0.0;

    //! time for solution of global system of equations
    double dt_solve_ = 0.0;

    //! equilibration method applied to system matrix
    const struct SSIMonoEquilibrationMethod equilibration_method_;

    //! Evaluation of coupling flux between scatra and manifold on scatra
    Teuchos::RCP<SSI::ScaTraManifoldScaTraFluxEvaluator> manifoldscatraflux_;

    //! type of global system matrix in global system of equations
    const Core::LinAlg::MatrixType matrixtype_;

    //! print system matrix, rhs, and map of system matrix to file
    const bool print_matlab_;

    //! relax the tolerance of the linear solver in case it is an iterative solver by scaling the
    //! convergence tolerance with factor @p relax_lin_solver_tolerance_
    const double relax_lin_solver_tolerance_;

    //! relax the tolerance of the linear solver within the first @p relax_lin_solver_step_ steps
    const int relax_lin_solver_iter_step_;

    //! all OD evaluation is in here
    Teuchos::RCP<SSI::ScatraStructureOffDiagCoupling> scatrastructure_off_diagcoupling_;

    //! algebraic solver for global system of equations
    Teuchos::RCP<Core::LinAlg::Solver> solver_;

    //! this object holds all maps relevant to monolithic scalar transport - structure interaction
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps_;

    //! this object holds the system matrix and all sub blocks
    Teuchos::RCP<SSI::UTILS::SSIMatrices> ssi_matrices_;

    //! this object holds the system residuals and increment
    Teuchos::RCP<SSI::UTILS::SSIVectors> ssi_vectors_;

    //! strategy how to assembly system matrix and rhs
    Teuchos::RCP<SSI::AssembleStrategyBase> strategy_assemble_;

    //! strategy how to apply contact contributions to sub matrices and rhs
    Teuchos::RCP<SSI::ContactStrategyBase> strategy_contact_;

    //! strategy for Newton-Raphson convergence check
    Teuchos::RCP<SSI::SsiMono::ConvCheckStrategyBase> strategy_convcheck_;

    //! all equilibration of global system matrix and RHS is done in here
    Teuchos::RCP<Core::LinAlg::Equilibration> strategy_equilibration_;

    //! strategy how to apply mesh tying on manifold domains
    Teuchos::RCP<SSI::ManifoldMeshTyingStrategyBase> strategy_manifold_meshtying_;

    //! strategy how to apply mesh tying to system matrix and rhs
    Teuchos::RCP<SSI::MeshtyingStrategyBase> strategy_meshtying_;

    //! timer for Newton-Raphson iteration
    Teuchos::RCP<Teuchos::Time> timer_;
  };
}  // namespace SSI
FOUR_C_NAMESPACE_CLOSE

#endif
