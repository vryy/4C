/*-----------------------------------------------------------*/
/*! \file

\brief Manager of the model evaluator calls.


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_MANAGER_HPP
#define FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_MANAGER_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"  // necessary due to enums

#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_Vector;
namespace NOX
{
  namespace Solver
  {
    class Generic;
  }  // namespace Solver
}  // namespace NOX

FOUR_C_NAMESPACE_OPEN

namespace Core::IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace Core::IO

namespace Core::LinAlg
{
  class SparseOperator;
  class SparseMatrix;
}  // namespace Core::LinAlg

namespace NOX
{
  namespace Nln
  {
    enum class CorrectionType : int;
    class Group;
  }  // namespace Nln
}  // namespace NOX

namespace Solid
{
  class Integrator;

  namespace TimeInt
  {
    class Base;
    class BaseDataSDyn;
    class BaseDataGlobalState;
    class BaseDataIO;
  }  // namespace TimeInt

  namespace ModelEvaluator
  {
    class Data;
    class Generic;
  }  // namespace ModelEvaluator

  /*! \brief Wrapper class for the Solid::ModelEvaluator::Generic derived objects.
   *
   * Manages the access to the different distinct model evaluators. Calling any routine on this
   * evaluator basically will trigger a loop over all active model evaluators that also implement
   * this routine.
   */
  class ModelEvaluatorManager
  {
   public:
    typedef std::map<enum Inpar::Solid::ModelType, Teuchos::RCP<Solid::ModelEvaluator::Generic>>
        Map;
    typedef std::vector<Teuchos::RCP<Solid::ModelEvaluator::Generic>> Vector;

    //! constructor
    ModelEvaluatorManager();

    //! destructor
    virtual ~ModelEvaluatorManager() = default;

    /*! \brief Initialize
     *
     * \todo Document remaining input arguments.
     *
     * \param[in] eval_data_ptr ??
     * \param[in] sdyn_ptr Pointer to the structural dynamic data container
     * \param[in] gstate_ptr Pointer to the global state data container
     * \param[in] gio_ptr Pointer to the input/output data container
     * \param[in] int_ptr ??
     * \param[in] timint_ptr Pointer to the underlying time integrator (read-only)
     */
    void init(const Teuchos::RCP<Solid::ModelEvaluator::Data>& eval_data_ptr,
        const Teuchos::RCP<Solid::TimeInt::BaseDataSDyn>& sdyn_ptr,
        const Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState>& gstate_ptr,
        const Teuchos::RCP<Solid::TimeInt::BaseDataIO>& gio_ptr,
        const Teuchos::RCP<Solid::Integrator>& int_ptr,
        const Teuchos::RCP<const Solid::TimeInt::Base>& timint_ptr);

    //! setup
    void setup();

    //! setup the MultiMapExtractor in the global state
    void setup_multi_map_extractor();

    //! @name General evaluate routines
    //!@{

    bool initialize_inertia_and_damping(const Epetra_Vector& x, Core::LinAlg::SparseOperator& jac);

    bool apply_initial_force(const Epetra_Vector& x, Epetra_Vector& f);

    /*! \brie Apply force
     *
     * @param[in] x Current solution
     * @param[in/out] f Residual vector (empty on input, filled on output)
     * @param[in] timefac_np Time integration factor for the current contribution at \f$t_{n+1}\f$
     *
     * @return Boolean flag to indicate success (true) or failure (false)
     */
    bool apply_force(const Epetra_Vector& x, Epetra_Vector& f, const double& timefac_np) const;

    /*! \brief Apply stiffness
     *
     * @param[in] x Current solution
     * @param[in/out] jac Jacobian matrix (empty on input, filled on output)
     * @param[in] timefac_np Time integration factor for the current contribution at \f$t_{n+1}\f$
     *
     * @return Boolean flag to indicate success (true) or failure (false)
     */
    bool apply_stiff(
        const Epetra_Vector& x, Core::LinAlg::SparseOperator& jac, const double& timefac_np) const;

    /*! \brief Apply model specific stiff
     *
     * @param[in] mt Type of model to be evaluated
     * @param[in] x Current solution
     * @param[in/out] jac Jacobian matrix (empty on input, filled on output)
     * @param[in] timefac_np Time integration factor for the current contribution at \f$t_{n+1}\f$
     *
     * @return Boolean flag to indicate success (true) or failure (false)
     */
    bool apply_stiff(const Inpar::Solid::ModelType& mt, const Epetra_Vector& x,
        Core::LinAlg::SparseOperator& jac, const double& timefac_np) const;

    /*! \brief Apply force and stiffness
     *
     * @param[in] x Current solution
     * @param[in/out] f Residual vector (empty on input, filled on output)
     * @param[in/out] jac Jacobian matrix (empty on input, filled on output)
     * @param[in] timefac_np Time integration factor for the current contribution at \f$t_{n+1}\f$
     *
     * @return Boolean flag to indicate success (true) or failure (false)
     */
    bool apply_force_stiff(const Epetra_Vector& x, Epetra_Vector& f,
        Core::LinAlg::SparseOperator& jac, const double& timefac_np) const;

    /*! \brief Compute cheap second order correction right hand side
     *
     * \todo Document remaining input arguments.
     *
     * @param type ??
     * @param constraint_models ??
     * @param[in] x Current solution
     * @param[in/out] f Residual vector (empty on input, filled on output)
     * @param[in] timefac_np Time integration factor for the current contribution at \f$t_{n+1}\f$
     *
     * @return Boolean flag to indicate success (true) or failure (false)
     */
    bool apply_cheap_soc_rhs(const enum NOX::Nln::CorrectionType type,
        const std::vector<Inpar::Solid::ModelType>& constraint_models, const Epetra_Vector& x,
        Epetra_Vector& f, const double& timefac_np) const;

    bool correct_parameters(const enum NOX::Nln::CorrectionType type) const;

    /*! \brief Remove any condensed contributions from the structural right-hand side
     *
     *  Recover the original right hand side vector by removing any
     *  contributions stemming from the condensation algorithm. A typical
     *  example is the condensed contact.
     *
     *  \author hiermeier \date 03/18 */
    void remove_condensed_contributions_from_rhs(Epetra_Vector& rhs) const;

    /*! \brief predict all internal variables in model evaluators
     *
     * @param[in] pred_type Type of predictor to be applied
     */
    void predict(const Inpar::Solid::PredEnum& pred_type) const;

    /** \brief Assembly of all force contributions
     *
     *  \param timefac_np           (in) : time integration factor for the current contribution
     *                                     \f$t_{n+1}\f$
     *  \param f                   (out) : force vector which is going to be assembled.
     *  \param without_these_models (in) : Assemble all models, except the models in this
     *                                     vector (optional)
     *  \return Boolean flag to indicate success (true) or failure (false)
     *
     *  \author hiermeier \date 03/17 */
    bool assemble_force(const double timefac_np, Epetra_Vector& f,
        const std::vector<Inpar::Solid::ModelType>* without_these_models) const;


    /** \brief Assembly of all jacobian contributions
     *
     *  \param timefac_np           (in) : time integration factor for the current contribution
     *                                     \f$t_{n+1}\f$
     *  \param jac                 (out) : jacobian which is going to be assembled.
     *  \param without_these_models (in) : Assemble all models, except the models in this
     *                                     vector (optional)
     *
     *  \return Boolean flag to indicate success (true) or failure (false)
     *
     *  \author farah \date 07/17 */
    bool assemble_jacobian(const double timefac_np, Core::LinAlg::SparseOperator& jac,
        const std::vector<Inpar::Solid::ModelType>* without_these_models) const;

    /** \brief Assembly of all force contributions
     *
     *  \param timefac_np (in) : time integration factor for the current contribution \f$t_{n+1}\f$
     *  \param f         (out) : force vector which is going to be assembled.
     *
     *  \return Boolean flag to indicate success (true) or failure (false)
     *
     *  \author hiermeier \date 03/17 */
    inline bool assemble_force(const double timefac_np, Epetra_Vector& f) const
    {
      return assemble_force(*me_vec_ptr_, timefac_np, f);
    }

    /** \brief Assembly of all Jacobian contributions
     *
     *  \param timefac_np (in) : time integration factor for the current contribution \f$t_{n+1}\f$
     *  \param jac        (out) : jacobian which is going to be assembled.
     *
     *  \return Boolean flag to indicate success (true) or failure (false)
     */
    inline bool assemble_jacobian(const double timefac_np, Core::LinAlg::SparseOperator& jac) const
    {
      return assemble_jacobian(*me_vec_ptr_, timefac_np, jac);
    }

    /** \brief Assembly of a sub-set of force contributions
     *
     *  \param me_vec     (in) : user provided vector of all models which are going to be assembled
     *  \param timefac_np (in) : time integration factor for the current contribution \f$t_{n+1}\f$
     *  \param f         (out) : force vector which is going to be assembled.
     *
     *  \return Boolean flag to indicate success (true) or failure (false)
     *
     *  \author hiermeier \date 03/17 */
    inline bool assemble_force(
        const Vector& me_vec, const double timefac_np, Epetra_Vector& f) const
    {
      bool ok = true;
      assemble_force(ok, me_vec, timefac_np, f);
      return ok;
    }

    /*! brief Assemble element contributions to global PTC scaling operator
     *
     * \todo Document remaining input arguments.
     *
     * @param me_vec User-provided vector of all models which are going to be assembled
     * @param timefac_np Time integration factor for the current contribution at \f$t_{n+1}\f$
     * @param modjac ??
     */
    void assemble_jacobian_contributions_from_element_level_for_ptc(const Vector& me_vec,
        const double timefac_np, Teuchos::RCP<Core::LinAlg::SparseMatrix>& modjac);

    /** \brief Assembly of a sub-set of stiffness contributions
     *
     *  \param me_vec     (in) : user provided vector of all models which are going to be assembled
     *  \param timefac_np (in) : time integration factor for the current contribution \f$t_{n+1}\f$
     *  \param jac       (out) : jacobian which is going to be assembled.
     *
     *  \author farah \date 07/17 */
    inline bool assemble_jacobian(
        const Vector& me_vec, const double timefac_np, Core::LinAlg::SparseOperator& jac) const
    {
      bool ok = true;
      assemble_jacobian(ok, me_vec, timefac_np, jac);
      return ok;
    }

    //!@}

    void create_backup_state(const Epetra_Vector& dir);

    void recover_from_backup_state();

    /*! \brief reset all model states (incl. the structural dynamic state)
     *
     *  \param x (in) : current full state vector */
    void reset_states(const Epetra_Vector& x) const;

    /*! \brief reset all model states (optional even. the structural dynamic
     *  state)
     *
     *  \param x (in) : current full state vector
     *  \param setstate (in) : flag to set state */
    void reset_states(const Epetra_Vector& x, bool setstate) const;

    /*! \brief reset a sub-set of all model states (optional even the structural
     *  dynamic state)
     *
     *  \param x (in) : current full state vector
     *  \param setstate (in) : flag to set state
     *  \param me_vec (in) : vector containing the sub-set of model evaluators */
    void reset_states(const Epetra_Vector& x, bool setstate, Vector& me_vec) const;

    //! Write current restart
    void write_restart(
        Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const;

    //! Read restart information
    void read_restart(Core::IO::DiscretizationReader& ioreader);

    //! @name Accessors
    //!@{

    //! return global state (read-only)
    const Solid::TimeInt::BaseDataGlobalState& get_global_state() const;

    //! return global state pointer (read and write access of the data)
    const Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState>& global_state_ptr();

    //! return pointer to the underlying time integrator (read-only)
    const Teuchos::RCP<const Solid::TimeInt::Base>& get_tim_int_ptr() const;

    /*! \brief Access one specific model evaluator
     *
     * \param[in] mt Type of model evaluator to be accessed
     */
    Solid::ModelEvaluator::Generic& evaluator(const enum Inpar::Solid::ModelType& mt);
    const Solid::ModelEvaluator::Generic& evaluator(const enum Inpar::Solid::ModelType& mt) const;

    //!@}

    //! @name Monolithic update routines
    //!@{

    //! Update configuration after time step
    void update_step_state(const double& timefac_n);

    //! Update everything on element level after time step and after output
    void update_step_element();

    //! Compute the residual by difference of {n+1} and {n} state
    void update_residual();

    //! calculation of stresses and strains
    void determine_stress_strain();

    //! calculation of engery
    void determine_energy();

    //! calculation of an optional quantitiy
    void determine_optional_quantity();

    //! Write the current step state
    void output_step_state(Core::IO::DiscretizationWriter& iowriter) const;

    /**
     * \brief Do stuff that has to be done before the runtime output is written.
     */
    void runtime_pre_output_step_state();

    //! Write the current step state during runtime
    void runtime_output_step_state() const;

    //! Do things after writing output
    void post_output();

    /*!
     * \brief Reset the current state variables to the ones of the previous timestep
     *
     * This is used for example to output the last successfull timestep.
     */
    void reset_step_state();

    /*! \brief Recover the current state
     *
     * Necessary for condensed systems, e.g. EAS, dual mortar, etc.*/
    void run_recover() const;

    /*! \brief Recover the current state
     *
     * Necessary for condensed systems, e.g. EAS, dual mortar, etc.*/
    void run_post_compute_x(const Epetra_Vector& xold, const Epetra_Vector& dir, const double& step,
        const Epetra_Vector& xnew, const bool isdefaultstep) const;

    /*! \brief Executed before the solution vector is going to be updated
     *
     *  \author hiermeier \date 03/17 */
    void run_pre_compute_x(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
        const double& step, const NOX::Nln::Group& curr_grp, const bool isdefaultstep) const;

    /*! \brief Executed at the end of the ::NOX::Solver::Step() (f.k.a. Iterate()) method
     *
     *  \author hiermeier \date 03/17 */
    void run_post_iterate(const ::NOX::Solver::Generic& solver, const double step,
        const bool isdefaultstep, const int num_corrs) const;

    /*! \brief Executed at the beginning of the ::NOX::Solver::solve() method
     *
     *  \author hiermeier */
    void run_pre_solve(
        const ::NOX::Solver::Generic& solver, const double step, const bool isdefaultstep) const;

    /*! \brief Executed at the end of the NOX::Nln::Group::applyJacobianInverse
     *  method
     *
     *  \author hiermeier \date 12/17 */
    void run_post_apply_jacobian_inverse(const Epetra_Vector& rhs, Epetra_Vector& result,
        const Epetra_Vector& xold, const NOX::Nln::Group& grp) const;

    /*! \brief Executed before the solution of the linear system
     *
     *  \author seitz \date 04/17 */
    void run_pre_apply_jacobian_inverse(const Epetra_Vector& rhs, Epetra_Vector& result,
        const Epetra_Vector& xold, const NOX::Nln::Group& grp) const;

    //!@}

    //! computes element based scaling contributions for PTC
    void compute_jacobian_contributions_from_element_level_for_ptc(
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& scalingMatrixOpPtr);

   protected:
    //! Returns the init flag.
    inline const bool& is_init() const { return isinit_; };

    //! Returns the setup flag.
    inline const bool& is_setup() const { return issetup_; };

    //! Check the init and setup state.
    void check_init_setup() const;

    //! Check the init state
    void check_init() const;

   private:
    Teuchos::RCP<Solid::ModelEvaluatorManager::Vector> transform_to_vector(
        const Solid::ModelEvaluatorManager::Map& model_map) const;

    /** \brief Assembly of all force contributions
     *
     *  \note There is also a PUBLIC alternative.
     *
     *  \param ok     (in/out) : flag indicating whether anything went wrong or not
     *  \param me_vec     (in) : vector of all models which are going to be assembled
     *  \param timefac_np (in) : time integration factor for the current contribution \f$t_{n+1}\f$
     *  \param f         (out) : force vector which is going to be assembled.
     *
     *  \author hiermeier \date 03/17 */
    void assemble_force(
        bool& ok, const Vector& me_vec, const double timefac_np, Epetra_Vector& f) const;

    /** \brief Assembly of all jacobian contributions
     *
     *  \param ok     (in/out) : flag indicating whether anything went wrong or not
     *  \param me_vec     (in) : vector of all models which are going to be assembled
     *  \param timefac_np (in) : time integration factor for the current contribution \f$t_{n+1}\f$
     *  \param f         (out) : jacobian matrix object which is going to be assembled.
     *
     *  \author hiermeier \date 03/17 */
    void assemble_jacobian(bool& ok, const Vector& me_vec, const double timefac_np,
        Core::LinAlg::SparseOperator& jac) const;

    void assemble_cheap_soc_rhs(
        bool& ok, const Vector& me_vec, const double timefac_np, Epetra_Vector& f) const;

    void evaluate_force(bool& ok, const Vector& me_vec) const;

    void evaluate_stiff(bool& ok, const Vector& me_vec) const;

    void evaluate_force_stiff(bool& ok, const Vector& me_vec) const;

    void evaluate_cheap_soc_rhs(bool& ok, const Vector& me_vec) const;

    void post_evaluate(bool ok, const Vector& me_vec) const;

    void pre_evaluate(bool ok, const Vector& me_vec) const;

    /** \brief split the internally stored model vector and get the set without
     *  the specified models */
    void split_model_vector(Vector& partial_me_vec,
        const std::vector<Inpar::Solid::ModelType>& without_these_models) const;

    /** \brief Extract from the internally stored model vector all models
     *  with the desired types */
    void extract_model_vector(Solid::ModelEvaluatorManager::Vector& partial_me_vec,
        const std::vector<Inpar::Solid::ModelType>& only_these_models) const;

   private:
    //! Flag to indicate whether init() has been called
    bool isinit_;

    //! Flag to indicate whether setup() has been called
    bool issetup_;

    Teuchos::RCP<Solid::ModelEvaluatorManager::Map> me_map_ptr_;

    Teuchos::RCP<Solid::ModelEvaluatorManager::Vector> me_vec_ptr_;

    Teuchos::RCP<Solid::ModelEvaluator::Data> eval_data_ptr_;

    //! Pointer to the structural dynamic data container
    Teuchos::RCP<Solid::TimeInt::BaseDataSDyn> sdyn_ptr_;

    //! Pointer to the global state data container
    Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState> gstate_ptr_;

    //! Pointer to the input/output data container
    Teuchos::RCP<Solid::TimeInt::BaseDataIO> gio_ptr_;

    Teuchos::RCP<Solid::Integrator> int_ptr_;

    //! Pointer to the underlying time integrator (read-only)
    Teuchos::RCP<const Solid::TimeInt::Base> timint_ptr_;

  };  // class ModelEvaluatorManager
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
