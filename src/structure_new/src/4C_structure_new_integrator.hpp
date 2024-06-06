/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all explicit/implicit time integrators.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_INTEGRATOR_HPP
#define FOUR_C_STRUCTURE_NEW_INTEGRATOR_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"               // enumerators
#include "4C_solver_nonlin_nox_enum_lists.hpp"  // enumerators

#include <NOX_Abstract_Vector.H>  // enumerators
#include <Teuchos_RCP.hpp>

// forward declarations
class Epetra_Vector;
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos

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

namespace Core::IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace Core::IO

namespace STR
{
  class ModelEvaluator;
  class Dbc;
  class MonitorDbc;

  namespace MODELEVALUATOR
  {
    class Data;
    class Generic;
  }  // namespace MODELEVALUATOR

  namespace TimeInt
  {
    class Base;
    class BaseDataSDyn;
    class BaseDataGlobalState;
    class BaseDataIO;
  }  // namespace TimeInt

  /*! \brief Base class for all structural time integrators
   *
   */
  class Integrator
  {
   public:
    //! constructor
    Integrator();

    //! destructor
    virtual ~Integrator() = default;

    /*! \brief Initialization
     *
     * \param[in] sdyn_ptr Pointer to the structural dynamic data container
     * \param[in] gstate_ptr Pointer to the global state data container
     * \param[in] gio_ptr Pointer to the input/output data container
     * \param[in] dpc_ptr Pointer to the dirichlet boundary condition object
     * \param[in] timint_ptr Pointer to the underlying time integrator (read-only)
     */
    virtual void Init(const Teuchos::RCP<STR::TimeInt::BaseDataSDyn>& sdyn_ptr,
        const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& gstate_ptr,
        const Teuchos::RCP<STR::TimeInt::BaseDataIO>& gio_ptr,
        const Teuchos::RCP<STR::Dbc>& dbc_ptr,
        const Teuchos::RCP<const STR::TimeInt::Base>& timint_ptr);

    //! Setup (has to be implemented by the derived classes)
    virtual void Setup();

    //! Post setup operation (compute initial equilibrium state), should be run directly after the
    //! setup routine has been finished
    virtual void post_setup() = 0;

    //! Set state variables
    virtual void set_state(const Epetra_Vector& x) = 0;

    //! Set initial displacement field
    virtual void set_initial_displacement(
        const Inpar::STR::InitialDisp init, const int startfuncno);

    /*! \brief Reset state variables of all models
     *  (incl. the structural dynamic state variables)
     *
     *  \param x (in) : current full state vector */
    void reset_model_states(const Epetra_Vector& x);

    /*! \brief Add the viscous and mass contributions to the right hand side (TR-rule)
     *
     * \remark The remaining contributions have been considered in the corresponding model
     *         evaluators. This is due to the fact, that some models use a different
     *         time integration scheme for their terms (e.g. GenAlpha for the structure
     *         and OST for the remaining things). */
    virtual void add_visco_mass_contributions(Epetra_Vector& f) const = 0;

    /*! \brief Add the viscous and mass contributions to the jacobian (TR-rule)
     *
     * \remark The remaining blocks have been considered in the corresponding model
     *         evaluators. This is due to the fact, that some models use a different
     *         time integration scheme for their terms (e.g. GenAlpha for the structure
     *         and OST for the remaining things). Furthermore, constraint/Lagrange
     *         multiplier blocks need no scaling anyway. */
    virtual void add_visco_mass_contributions(Core::LinAlg::SparseOperator& jac) const = 0;

    //! Apply the right hand side only
    virtual bool apply_force(const Epetra_Vector& x, Epetra_Vector& f) = 0;

    /*! \brief Apply the stiffness only
     *
     * Normally this one is unnecessary, since it makes more sense
     * to evaluate the stiffness and right hand side at once, because of
     * the lower computational overhead. */
    virtual bool apply_stiff(const Epetra_Vector& x, Core::LinAlg::SparseOperator& jac) = 0;

    /*! \brief Apply force and stiff at once
     *
     *  Only one loop over all elements. Especially in the contact case,
     *  the difference between this call and first call apply_force and
     *  then apply_stiff is mentionable because of the projection operations. */
    virtual bool apply_force_stiff(
        const Epetra_Vector& x, Epetra_Vector& f, Core::LinAlg::SparseOperator& jac) = 0;

    /*! \brief Modify the right hand side and Jacobian corresponding to the requested correction
     * action of one (or several) second order constraint (SOC) model(s)
     */
    virtual bool apply_correction_system(const enum NOX::Nln::CorrectionType type,
        const std::vector<Inpar::STR::ModelType>& constraint_models, const Epetra_Vector& x,
        Epetra_Vector& f, Core::LinAlg::SparseOperator& jac) = 0;

    /*! \brief Remove contributions from the structural right-hand side stemming
     *  from any condensation operations (typical example is contact)
     */
    virtual void remove_condensed_contributions_from_rhs(Epetra_Vector& rhs) const = 0;

    /*! \brief Calculate characteristic/reference norms for forces
     *
     *  The reference norms are used to scale the calculated iterative
     *  displacement norm and/or the residual force norm. For this
     *  purpose we only need the right order of magnitude, so we don't
     *  mind evaluating the corresponding norms at possibly different
     *  points within the time-step (end point, generalized midpoint). */
    virtual double calc_ref_norm_force(
        const enum ::NOX::Abstract::Vector::NormType& type) const = 0;

    //! compute the scaling operator for element based scaling using PTC
    virtual void compute_jacobian_contributions_from_element_level_for_ptc(
        Teuchos::RCP<Core::LinAlg::SparseMatrix>& scalingMatrixOpPtr) = 0;

    //! Assemble the right hand side
    virtual bool assemble_force(Epetra_Vector& f,
        const std::vector<Inpar::STR::ModelType>* without_these_models = nullptr) const = 0;

    //! Assemble Jacobian
    virtual bool assemble_jac(Core::LinAlg::SparseOperator& jac,
        const std::vector<Inpar::STR::ModelType>* without_these_models = nullptr) const
    {
      return false;
    };

    //! Create backup state
    void create_backup_state(const Epetra_Vector& dir);

    //! recover state from the stored backup
    void recover_from_backup_state();

    //! return integration factor
    virtual double get_int_param() const = 0;
    virtual double get_acc_int_param() const { return get_int_param(); }

    /*!
     * \brief Allows to stop the simulation before the max. time or max timestep is reached
     *
     * \return true In case of an early stop is desired
     * \return false In case of no early stop (default)
     */
    virtual bool early_stopping() const { return false; }

    //! @name Restart and output related functions
    //!@{

    /*! write restart information of the different time integration schemes
     *  and model evaluators */
    void write_restart(Core::IO::DiscretizationWriter& iowriter) const
    {
      write_restart(iowriter, false);
    };
    virtual void write_restart(
        Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const = 0;

    /*! read restart information of the different time integration schemes
     *  and model evaluators */
    virtual void read_restart(Core::IO::DiscretizationReader& ioreader) = 0;

    //!@}

    //! @name Monolithic update routines
    //!@{

    //! things that should be done before updating
    virtual void pre_update() = 0;

    /*! \brief Update configuration after time step
     *
     *  Thus the 'last' converged is lost and a reset of the time step
     *  becomes impossible. We are ready and keen awaiting the next
     *  time step. */
    virtual void update_step_state() = 0;

    /*! \brief Update everything on element level after time step and after output
     *
     *  Thus the 'last' converged is lost and a reset of the time step
     *  becomes impossible. We are ready and keen awaiting the next
     *  time step. */
    virtual void update_step_element() = 0;

    //! calculate stresses and strains in the different model evaluators
    void determine_stress_strain();

    //! calculate the energy in the different model evaluators
    void determine_energy();

    //! get the model value in accordance with the currently active time integration
    virtual double get_model_value(const Epetra_Vector& x);

    /*! return the total structural energy evaluated at the actual mid-time
     *  in accordance to the used time integration scheme */
    double get_total_mid_time_str_energy(const Epetra_Vector& x);

    /// update the structural energy variable in the end of a successful time step
    void update_structural_energy();

    //! calculate an optional quantity in the different model evaluators
    void determine_optional_quantity();

    /// compute the current volumes of all elements
    bool determine_element_volumes(const Epetra_Vector& x, Teuchos::RCP<Epetra_Vector>& ele_vols);

    /*! \brief Output to file
     *
     *  This routine prints always the last converged state, i.e.
     *  \f$D_{n}, V_{n}, A_{n}\f$. So, #UpdateIncrement should be called
     *  upon object prior to writing stuff here.
     *
     *  \author mwgee (originally)                         \date 03/07 */
    void output_step_state(Core::IO::DiscretizationWriter& iowriter) const;

    /**
     * \brief Do stuff that has to be done before runtime output is written.
     */
    void runtime_pre_output_step_state();

    /*! \brief Output to file during runtime (no filter necessary)
     *
     *  This routine prints always the last converged state, i.e.
     *  \f$D_{n}, V_{n}, A_{n}\f$. So, #UpdateIncrement should be called
     *  upon object prior to writing stuff here.
     *
     *                                                     \date 04/17 */
    void runtime_output_step_state() const;

    /** \brief reset step configuration after time step
     *
     *  This function is supposed to reset all variables which are directly related
     *  to the current new step n+1. To be more precise all variables ending with "Np"
     *  have to be reseted. */
    virtual void reset_step_state();

    /// things that should be done after updating
    virtual void post_update() = 0;

    /// things that should be done after the timeloop
    virtual void post_time_loop(){};

    //! update constant contributions of the current state for the new time step \f$t_{n+1}\f$
    virtual void update_constant_state_contributions() = 0;

    //! things that should be done after output
    virtual void post_output();

    void monitor_dbc(Core::IO::DiscretizationWriter& writer) const;
    //!@}

    //! @name Accessors
    //!@{

    double get_condensed_update_norm(const enum NOX::Nln::StatusTest::QuantityType& qtype) const;

    double get_condensed_previous_sol_norm(
        const enum NOX::Nln::StatusTest::QuantityType& qtype) const;

    double get_condensed_solution_update_rms(
        const enum NOX::Nln::StatusTest::QuantityType& qtype) const;

    int get_condensed_dof_number(const enum NOX::Nln::StatusTest::QuantityType& qtype) const;

    //! Return the model evaluator control object (read and write)
    STR::ModelEvaluator& model_eval();

    //! Return the model evaluator control object (read-only)
    const STR::ModelEvaluator& model_eval() const;
    Teuchos::RCP<const STR::ModelEvaluator> model_eval_ptr() const;

    //! Return the model evaluator object for the given model type
    STR::MODELEVALUATOR::Generic& evaluator(const Inpar::STR::ModelType& mt);

    //! Return the model evaluator object for the given model type
    const STR::MODELEVALUATOR::Generic& evaluator(const Inpar::STR::ModelType& mt) const;

    //! Return the model evaluator data object (read-only)
    const STR::MODELEVALUATOR::Data& eval_data() const;

    //! Return the model evaluator data object (read and write access)
    STR::MODELEVALUATOR::Data& eval_data();

    //! Return the Dirichlet boundary condition object (read-only)
    const STR::Dbc& get_dbc() const;

    //!@}

   protected:
    //! returns init state
    inline const bool& is_init() const { return isinit_; };

    //! returns setup state
    inline const bool& is_setup() const { return issetup_; };

    //! Check the init state
    void check_init() const;

    //! Check the setup state
    void check_init_setup() const;

    /*! \brief Equilibriate system at initial state and identify consistent accelerations
     */
    void equilibrate_initial_state();

    /*! \brief Check if current state is equilibrium (with respect to
     *  a given tolerance of the inf-norm)
     *
     *  This is a sanity check in case of nonlinear mass
     *  and non-additive rotation pseudo-vector DoFs where
     *  determination of consistent accelerations is more intricate
     *  and not supported yet */
    bool current_state_is_equilibrium(const double& tol);

    //! Return the structural dynamic data container
    STR::TimeInt::BaseDataSDyn& sdyn();

    //! Return the structural dynamic data container (read-only)
    const STR::TimeInt::BaseDataSDyn& s_dyn() const;

    //! Return the global state data container
    STR::TimeInt::BaseDataGlobalState& global_state();

    //! Return the global state data container (read-only)
    const STR::TimeInt::BaseDataGlobalState& global_state() const;

    //! Return the Dirichlet boundary condition object
    STR::Dbc& dbc();

    //! Return the time integration strategy object (read-only)
    const STR::TimeInt::Base& tim_int() const;

    //! reset the time step dependent parameters for the element evaluation
    virtual void reset_eval_params(){};

    double get_condensed_global_norm(const enum NOX::Nln::StatusTest::QuantityType& qtype,
        const enum ::NOX::Abstract::Vector::NormType& normtype, double& mynorm) const;

   protected:
    //! indicates if the Init() function has been called
    bool isinit_;

    //! indicates if the Setup() function has been called
    bool issetup_;

    //! Mid-time energy container
    struct MidTimeEnergy
    {
      /// constructor
      MidTimeEnergy(const Integrator& integrator);

      /// setup
      void Setup();

      /// can this container be used?
      bool is_correctly_configured() const;

      bool store_energy_n() const;

      /// print energy info to output stream
      void Print(std::ostream& os) const;

      /// Get total energy measure in accordance to the surrounding time integrator
      double get_total() const;

      /// average quantities based on the used averaging type
      Teuchos::RCP<const Epetra_Vector> Average(
          const Epetra_Vector& state_np, const Epetra_Vector& state_n, const double fac_n) const;

      /// copy current state to old state (during update)
      void CopyNpToN();

      /// internal energy at \f$t_{n+1}\f$
      double int_energy_np_ = 0.0;

      /// external energy at \f$t_{n+1}\f$
      double ext_energy_np_ = 0.0;

      /// kinetic energy at \f$t_{n+1}\f$
      double kin_energy_np_ = 0.0;

      /// internal energy at \f$t_{n}\f$
      double int_energy_n_ = 0.0;

      /// external energy at \f$t_{n}\f$
      double ext_energy_n_ = 0.0;

      /// kinetic energy at \f$t_{n}\f$
      double kin_energy_n_ = 0.0;

     private:
      /// call-back
      const Integrator& integrator_;

      /// mid-time energy averaging type
      enum Inpar::STR::MidAverageEnum avg_type_;

      /// setup flag
      bool issetup_ = false;
    };
    MidTimeEnergy mt_energy_ = MidTimeEnergy(*this);

   private:
    //! pointer to the model evaluator
    Teuchos::RCP<STR::ModelEvaluator> modelevaluator_ptr_;

    //! pointer to model evaluator data
    Teuchos::RCP<STR::MODELEVALUATOR::Data> eval_data_ptr_;

    //! pointer to the structural dynamic data container
    Teuchos::RCP<STR::TimeInt::BaseDataSDyn> sdyn_ptr_;

    //! pointer to the global state data container
    Teuchos::RCP<STR::TimeInt::BaseDataGlobalState> gstate_ptr_;

    //! pointer to the input/output data container
    Teuchos::RCP<STR::TimeInt::BaseDataIO> io_ptr_;

    //! pointer to the dirichlet boundary condition object
    Teuchos::RCP<STR::Dbc> dbc_ptr_;

    //! pointer to the underlying time integrator (read-only)
    Teuchos::RCP<const STR::TimeInt::Base> timint_ptr_;

    //! pointer to the dirichlet boundary condition monitor
    Teuchos::RCP<STR::MonitorDbc> monitor_dbc_ptr_ = Teuchos::null;
  };  // namespace STR
}  // namespace STR

FOUR_C_NAMESPACE_CLOSE

#endif
