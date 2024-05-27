/*-----------------------------------------------------------*/
/*! \file

\brief Manager of the model evaluator calls.


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_HPP
#define FOUR_C_STRUCTURE_NEW_MODEL_EVALUATOR_HPP

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

namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace IO

namespace CORE::LINALG
{
  class SparseOperator;
  class SparseMatrix;
}  // namespace CORE::LINALG

namespace NOX
{
  namespace NLN
  {
    enum class CorrectionType : int;
    class Group;
  }  // namespace NLN
}  // namespace NOX

namespace STR
{
  class Integrator;

  namespace TIMINT
  {
    class Base;
    class BaseDataSDyn;
    class BaseDataGlobalState;
    class BaseDataIO;
  }  // namespace TIMINT

  namespace MODELEVALUATOR
  {
    class Data;
    class Generic;
  }  // namespace MODELEVALUATOR

  /*! \brief Wrapper class for the STR::MODELEVALUATOR::Generic derived objects.
   *
   * Manages the access to the different distinct model evaluators. Calling any routine on this
   * evaluator basically will trigger a loop over all active model evaluators that also implement
   * this routine.
   */
  class ModelEvaluator
  {
   public:
    typedef std::map<enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic>> Map;
    typedef std::vector<Teuchos::RCP<STR::MODELEVALUATOR::Generic>> Vector;

    //! constructor
    ModelEvaluator();

    //! destructor
    virtual ~ModelEvaluator() = default;

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
    void Init(const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataSDyn>& sdyn_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
        const Teuchos::RCP<STR::Integrator>& int_ptr,
        const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr);

    //! setup
    void Setup();

    //! setup the MultiMapExtractor in the global state
    void setup_multi_map_extractor();

    //! @name General evaluate routines
    //!@{

    bool initialize_inertia_and_damping(const Epetra_Vector& x, CORE::LINALG::SparseOperator& jac);

    bool ApplyInitialForce(const Epetra_Vector& x, Epetra_Vector& f);

    /*! \brie Apply force
     *
     * @param[in] x Current solution
     * @param[in/out] f Residual vector (empty on input, filled on output)
     * @param[in] timefac_np Time integration factor for the current contribution at \f$t_{n+1}\f$
     *
     * @return Boolean flag to indicate success (true) or failure (false)
     */
    bool ApplyForce(const Epetra_Vector& x, Epetra_Vector& f, const double& timefac_np) const;

    /*! \brief Apply stiffness
     *
     * @param[in] x Current solution
     * @param[in/out] jac Jacobian matrix (empty on input, filled on output)
     * @param[in] timefac_np Time integration factor for the current contribution at \f$t_{n+1}\f$
     *
     * @return Boolean flag to indicate success (true) or failure (false)
     */
    bool ApplyStiff(
        const Epetra_Vector& x, CORE::LINALG::SparseOperator& jac, const double& timefac_np) const;

    /*! \brief Apply model specific stiff
     *
     * @param[in] mt Type of model to be evaluated
     * @param[in] x Current solution
     * @param[in/out] jac Jacobian matrix (empty on input, filled on output)
     * @param[in] timefac_np Time integration factor for the current contribution at \f$t_{n+1}\f$
     *
     * @return Boolean flag to indicate success (true) or failure (false)
     */
    bool ApplyStiff(const INPAR::STR::ModelType& mt, const Epetra_Vector& x,
        CORE::LINALG::SparseOperator& jac, const double& timefac_np) const;

    /*! \brief Apply force and stiffness
     *
     * @param[in] x Current solution
     * @param[in/out] f Residual vector (empty on input, filled on output)
     * @param[in/out] jac Jacobian matrix (empty on input, filled on output)
     * @param[in] timefac_np Time integration factor for the current contribution at \f$t_{n+1}\f$
     *
     * @return Boolean flag to indicate success (true) or failure (false)
     */
    bool ApplyForceStiff(const Epetra_Vector& x, Epetra_Vector& f,
        CORE::LINALG::SparseOperator& jac, const double& timefac_np) const;

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
    bool ApplyCheapSOCRhs(const enum NOX::NLN::CorrectionType type,
        const std::vector<INPAR::STR::ModelType>& constraint_models, const Epetra_Vector& x,
        Epetra_Vector& f, const double& timefac_np) const;

    bool correct_parameters(const enum NOX::NLN::CorrectionType type) const;

    /*! \brief Remove any condensed contributions from the structural right-hand side
     *
     *  Recover the original right hand side vector by removing any
     *  contributions stemming from the condensation algorithm. A typical
     *  example is the condensed contact.
     *
     *  \author hiermeier \date 03/18 */
    void remove_condensed_contributions_from_rhs(Epetra_Vector& rhs) const;

    /*! \brief Predict all internal variables in model evaluators
     *
     * @param[in] pred_type Type of predictor to be applied
     */
    void Predict(const INPAR::STR::PredEnum& pred_type) const;

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
        const std::vector<INPAR::STR::ModelType>* without_these_models) const;


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
    bool assemble_jacobian(const double timefac_np, CORE::LINALG::SparseOperator& jac,
        const std::vector<INPAR::STR::ModelType>* without_these_models) const;

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
    inline bool assemble_jacobian(const double timefac_np, CORE::LINALG::SparseOperator& jac) const
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
        const double timefac_np, Teuchos::RCP<CORE::LINALG::SparseMatrix>& modjac);

    /** \brief Assembly of a sub-set of stiffness contributions
     *
     *  \param me_vec     (in) : user provided vector of all models which are going to be assembled
     *  \param timefac_np (in) : time integration factor for the current contribution \f$t_{n+1}\f$
     *  \param jac       (out) : jacobian which is going to be assembled.
     *
     *  \author farah \date 07/17 */
    inline bool assemble_jacobian(
        const Vector& me_vec, const double timefac_np, CORE::LINALG::SparseOperator& jac) const
    {
      bool ok = true;
      assemble_jacobian(ok, me_vec, timefac_np, jac);
      return ok;
    }

    //!@}

    void CreateBackupState(const Epetra_Vector& dir);

    void recover_from_backup_state();

    /*! \brief reset all model states (incl. the structural dynamic state)
     *
     *  \param x (in) : current full state vector */
    void ResetStates(const Epetra_Vector& x) const;

    /*! \brief reset all model states (optional even. the structural dynamic
     *  state)
     *
     *  \param x (in) : current full state vector
     *  \param setstate (in) : flag to set state */
    void ResetStates(const Epetra_Vector& x, bool setstate) const;

    /*! \brief reset a sub-set of all model states (optional even the structural
     *  dynamic state)
     *
     *  \param x (in) : current full state vector
     *  \param setstate (in) : flag to set state
     *  \param me_vec (in) : vector containing the sub-set of model evaluators */
    void ResetStates(const Epetra_Vector& x, bool setstate, Vector& me_vec) const;

    //! Write current restart
    void write_restart(IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const;

    //! Read restart information
    void read_restart(IO::DiscretizationReader& ioreader);

    //! @name Accessors
    //!@{

    //! return global state (read-only)
    const STR::TIMINT::BaseDataGlobalState& GetGlobalState() const;

    //! return global state pointer (read and write access of the data)
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& global_state_ptr();

    //! return pointer to the underlying time integrator (read-only)
    const Teuchos::RCP<const STR::TIMINT::Base>& GetTimIntPtr() const;

    /*! \brief Access one specific model evaluator
     *
     * \param[in] mt Type of model evaluator to be accessed
     */
    STR::MODELEVALUATOR::Generic& Evaluator(const enum INPAR::STR::ModelType& mt);
    const STR::MODELEVALUATOR::Generic& Evaluator(const enum INPAR::STR::ModelType& mt) const;

    //!@}

    //! @name Monolithic update routines
    //!@{

    //! Update configuration after time step
    void UpdateStepState(const double& timefac_n);

    //! Update everything on element level after time step and after output
    void UpdateStepElement();

    //! Compute the residual by difference of {n+1} and {n} state
    void update_residual();

    //! calculation of stresses and strains
    void determine_stress_strain();

    //! calculation of engery
    void DetermineEnergy();

    //! calculation of an optional quantitiy
    void determine_optional_quantity();

    //! Write the current step state
    void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    /**
     * \brief Do stuff that has to be done before the runtime output is written.
     */
    void runtime_pre_output_step_state();

    //! Write the current step state during runtime
    void runtime_output_step_state() const;

    //! Do things after writing output
    void PostOutput();

    /*!
     * \brief Reset the current state variables to the ones of the previous timestep
     *
     * This is used for example to output the last successfull timestep.
     */
    void ResetStepState();

    /*! \brief Recover the current state
     *
     * Necessary for condensed systems, e.g. EAS, dual mortar, etc.*/
    void RunRecover() const;

    /*! \brief Recover the current state
     *
     * Necessary for condensed systems, e.g. EAS, dual mortar, etc.*/
    void run_post_compute_x(const Epetra_Vector& xold, const Epetra_Vector& dir, const double& step,
        const Epetra_Vector& xnew, const bool isdefaultstep) const;

    /*! \brief Executed before the solution vector is going to be updated
     *
     *  \author hiermeier \date 03/17 */
    void run_pre_compute_x(const Epetra_Vector& xold, Epetra_Vector& dir_mutable,
        const double& step, const NOX::NLN::Group& curr_grp, const bool isdefaultstep) const;

    /*! \brief Executed at the end of the ::NOX::Solver::Step() (f.k.a. Iterate()) method
     *
     *  \author hiermeier \date 03/17 */
    void run_post_iterate(const ::NOX::Solver::Generic& solver, const double step,
        const bool isdefaultstep, const int num_corrs) const;

    /*! \brief Executed at the beginning of the ::NOX::Solver::solve() method
     *
     *  \author hiermeier */
    void RunPreSolve(
        const ::NOX::Solver::Generic& solver, const double step, const bool isdefaultstep) const;

    /*! \brief Executed at the end of the NOX::NLN::Group::applyJacobianInverse
     *  method
     *
     *  \author hiermeier \date 12/17 */
    void run_post_apply_jacobian_inverse(const Epetra_Vector& rhs, Epetra_Vector& result,
        const Epetra_Vector& xold, const NOX::NLN::Group& grp) const;

    /*! \brief Executed before the solution of the linear system
     *
     *  \author seitz \date 04/17 */
    void run_pre_apply_jacobian_inverse(const Epetra_Vector& rhs, Epetra_Vector& result,
        const Epetra_Vector& xold, const NOX::NLN::Group& grp) const;

    //!@}

    //! computes element based scaling contributions for PTC
    void compute_jacobian_contributions_from_element_level_for_ptc(
        Teuchos::RCP<CORE::LINALG::SparseMatrix>& scalingMatrixOpPtr);

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
    Teuchos::RCP<STR::ModelEvaluator::Vector> transform_to_vector(
        const STR::ModelEvaluator::Map& model_map) const;

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
        CORE::LINALG::SparseOperator& jac) const;

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
        const std::vector<INPAR::STR::ModelType>& without_these_models) const;

    /** \brief Extract from the internally stored model vector all models
     *  with the desired types */
    void extract_model_vector(STR::ModelEvaluator::Vector& partial_me_vec,
        const std::vector<INPAR::STR::ModelType>& only_these_models) const;

   private:
    //! Flag to indicate whether Init() has been called
    bool isinit_;

    //! Flag to indicate whether Setup() has been called
    bool issetup_;

    Teuchos::RCP<STR::ModelEvaluator::Map> me_map_ptr_;

    Teuchos::RCP<STR::ModelEvaluator::Vector> me_vec_ptr_;

    Teuchos::RCP<STR::MODELEVALUATOR::Data> eval_data_ptr_;

    //! Pointer to the structural dynamic data container
    Teuchos::RCP<STR::TIMINT::BaseDataSDyn> sdyn_ptr_;

    //! Pointer to the global state data container
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

    //! Pointer to the input/output data container
    Teuchos::RCP<STR::TIMINT::BaseDataIO> gio_ptr_;

    Teuchos::RCP<STR::Integrator> int_ptr_;

    //! Pointer to the underlying time integrator (read-only)
    Teuchos::RCP<const STR::TIMINT::Base> timint_ptr_;

  };  // class ModelEvaluator
}  // namespace STR


FOUR_C_NAMESPACE_CLOSE

#endif
