/*---------------------------------------------------------------------*/
/*! \file
\brief This strategy allows the combination of an arbitrary number of
       augmented contact solving strategies.

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef FOUR_C_CONTACT_AUG_COMBO_STRATEGY_HPP
#define FOUR_C_CONTACT_AUG_COMBO_STRATEGY_HPP

#include "4C_config.hpp"

#include "4C_contact_abstract_strategy.hpp"
#include "4C_contact_aug_utils.hpp"

FOUR_C_NAMESPACE_OPEN

// #define DEBUG_COMBO_STRATEGY

namespace CONTACT
{
  namespace Aug
  {
    class Strategy;
    class Interface;
    class DataContainer;

    /** \brief Allows the combination of different CONTACT::AUG strategies
     *
     *  The switch between the combined strategies is controlled by an internal
     *  switching strategy, which can be chosen by the user.
     *
     *  \author hiermeier */
    class ComboStrategy : public AbstractStrategy
    {
     public:
      /// create a new combo strategy object
      static Teuchos::RCP<AbstractStrategy> Create(
          const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& data,
          const Epetra_Map* dof_row_map, const Epetra_Map* node_row_map,
          const Teuchos::ParameterList& params, const plain_interface_set& ref_interfaces,
          const int dim, const Teuchos::RCP<const Epetra_Comm>& comm, const int maxdof,
          CONTACT::ParamsInterface* cparams_interface);

     private:
      /// create the interface objects for the combined strategies
      static void create_strategy_interfaces(const enum Inpar::CONTACT::SolvingStrategy strat_type,
          const plain_interface_set& ref_interfaces, plain_interface_set& strat_interfaces);

      /// create the linear solver objects for the different combined contact strategies
      static void create_strategy_linear_solvers(const CONTACT::AbstractStrategy& strategy,
          const std::string& lin_solver_id_str, const Teuchos::ParameterList& params,
          CONTACT::ParamsInterface* cparams_interface, plain_lin_solver_set& strat_lin_solvers);

     public:
      /// constructor
      ComboStrategy(const Teuchos::RCP<CONTACT::AbstractStratDataContainer>& stratData,
          const Epetra_Map* dof_row_map, const Epetra_Map* node_row_map,
          const Teuchos::ParameterList& params, const plain_strategy_set& strategies,
          const plain_lin_solver_set& lin_solvers, const int dim,
          const Teuchos::RCP<const Epetra_Comm>& comm, const int maxdof);

      /* Hide the public derived methods as well, since we do not want to access
       * any of these methods directly.                         hiermeier 03/17 */
     protected:
      /// @name Derived public methods from the base class ( a.k.a. AbstractStrategy )
      /// @{

      /// return the type of currently active wrapped contact strategy
      Inpar::CONTACT::SolvingStrategy Type() const override;

      /// return the linear solver for the currently active strategy
      Core::LinAlg::Solver* GetLinearSolver() const override;

      /// function wrapper
      void reset(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& dispnp,
          const Epetra_Vector& xnew) override;

      /// function wrapper
      bool IsSaddlePointSystem() const override;

      /// function wrapper
      void ResetActiveSet() override;

      /// function wrapper
      double ConstraintNorm() const override;

      /// function wrapper
      void SaveReferenceState(Teuchos::RCP<const Epetra_Vector> dis) override;

      /// function wrapper
      bool ActiveSetConverged() override;

      /// function wrapper
      int ActiveSetSteps() override;

      /// function wrapper
      void UpdateActiveSet() override;

      /// function wrapper
      void evaluate_rel_mov_predict() override;

      /// function wrapper
      bool active_set_semi_smooth_converged() const override;

      /// function wrapper
      Teuchos::RCP<const Epetra_Map> get_old_active_row_nodes() const override;

      /// function wrapper
      Teuchos::RCP<const Epetra_Map> GetOldSlipRowNodes() const override;

      /// function wrapper
      Teuchos::RCP<const Epetra_Map> sl_normal_do_f_row_map_ptr(const bool& redist) const override;

      /// function wrapper
      const Epetra_Map& SlNormalDoFRowMap(const bool& redist) const override;

      /// function wrapper
      Teuchos::RCP<const Epetra_Map> sl_tangential_do_f_row_map_ptr(
          const bool& redist) const override;

      /// function wrapper
      const Epetra_Map& sl_tangential_do_f_row_map(const bool& redist) const override;

      /// function wrapper
      Teuchos::RCP<const Epetra_Vector> GetRhsBlockPtr(
          const enum CONTACT::VecBlockType& bt) const override;

      /// function wrapper
      Teuchos::RCP<const Epetra_Vector> get_rhs_block_ptr_for_norm_check(
          const enum CONTACT::VecBlockType& bt) const override;

      /// function wrapper
      Teuchos::RCP<const Epetra_Vector> GetCondensedRhsPtr(
          Epetra_Vector& f, const double& timefac_np) const override;

      /// function wrapper
      Teuchos::RCP<Core::LinAlg::SparseMatrix> GetMatrixBlockPtr(
          const enum CONTACT::MatBlockType& bt,
          const CONTACT::ParamsInterface* cparams = nullptr) const override;

      /// function wrapper
      Teuchos::RCP<Core::LinAlg::SparseMatrix> get_condensed_matrix_block_ptr(
          Teuchos::RCP<Core::LinAlg::SparseMatrix>& kteff, const double& timefac_np) const override;

      /// function wrapper
      Teuchos::RCP<Epetra_Vector> ConstrRhs() override;

      /// function wrapper
      void initialize() override;

      /// function wrapper
      void EvalConstrRHS() override;

      /// function wrapper
      void update_active_set_semi_smooth(const bool firstStepPredictor = false) override;

      /// function wrapper
      void DoReadRestart(Core::IO::DiscretizationReader& reader,
          Teuchos::RCP<const Epetra_Vector> dis,
          Teuchos::RCP<CONTACT::ParamsInterface> cparams_ptr) override;

      /// function wrapper
      void Update(Teuchos::RCP<const Epetra_Vector> dis) override;

      /// function wrapper
      double GetPotentialValue(
          const enum NOX::Nln::MeritFunction::MeritFctName mrt_type) const override;

      /// function wrapper
      double get_linearized_potential_value_terms(const Epetra_Vector& dir,
          const enum NOX::Nln::MeritFunction::MeritFctName mrt_type,
          const enum NOX::Nln::MeritFunction::LinOrder linorder,
          const enum NOX::Nln::MeritFunction::LinType lintype) const override;

      /// function wrapper
      void WriteOutput(Core::IO::DiscretizationWriter& writer) const override;

      /// function wrapper
      void evaluate_reference_state() override;

      /// function wrapper
      bool dyn_redistribute_contact(const Teuchos::RCP<const Epetra_Vector>& dis,
          Teuchos::RCP<const Epetra_Vector> vel, const int nlniter) override;

      /// @}

      /// @{

      /// @name accessors to public augmented Lagrangian methods
      /// @{

      /// function wrapper
      bool was_in_contact_last_iter() const;

      /// @}

     protected:
      /// @name Derived protected methods from the base class ( a.k.a. AbstractStrategy )
      /// @{

      /// Get the set of currently active interfaces
      std::vector<Teuchos::RCP<CONTACT::Interface>>& interfaces() override;

      /// Get the set of currently active interfaces
      const std::vector<Teuchos::RCP<CONTACT::Interface>>& interfaces() const override;

      /// function wrapper
      void compute_contact_stresses() final;

      /** \brief function wrapper: Redistribute and setup augmented Lagrangian members
       *
       *  \author hiermeier \date 03/17 */
      void post_setup(bool redistributed, bool init) override;

      /** \brief function wrapper: Compute force terms
       *
       *  \author hiermeier \date 03/17 */
      void eval_force(CONTACT::ParamsInterface& cparams) override;

      /** \brief function wrapper: Compute force and stiffness terms
       *
       *  \author hiermeier \date 03/17 */
      void eval_force_stiff(CONTACT::ParamsInterface& cparams) override;

      /// function wrapper
      void eval_static_constraint_rhs(CONTACT::ParamsInterface& cparams) override;

      /// function wrapper
      void run_pre_evaluate(CONTACT::ParamsInterface& cparams) override;

      /// function wrapper
      void run_post_evaluate(CONTACT::ParamsInterface& cparams) override;

      /// function wrapper
      void run_pre_solve(const Teuchos::RCP<const Epetra_Vector>& curr_disp,
          const CONTACT::ParamsInterface& cparams) override;

      /// function wrapper
      void run_post_iterate(const CONTACT::ParamsInterface& cparams) override;

      /// function wrapper
      void run_pre_compute_x(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
          Epetra_Vector& dir_mutable) override;

      /// function wrapper
      void run_post_compute_x(const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xold,
          const Epetra_Vector& dir, const Epetra_Vector& xnew) override;

      /// function wrapper
      void run_post_apply_jacobian_inverse(const CONTACT::ParamsInterface& cparams,
          const Epetra_Vector& rhs, Epetra_Vector& result, const Epetra_Vector& xold,
          const NOX::Nln::Group& grp) override;

      /// function wrapper
      void remove_condensed_contributions_from_rhs(Epetra_Vector& str_rhs) const override;

      /// function wrapper
      void correct_parameters(
          CONTACT::ParamsInterface& cparams, const NOX::Nln::CorrectionType type) override;

      /// function wrapper
      void reset_lagrange_multipliers(
          const CONTACT::ParamsInterface& cparams, const Epetra_Vector& xnew) override;

      void post_store_dirichlet_status(
          Teuchos::RCP<const Core::LinAlg::MapExtractor> dbcmaps) override;

      /// @}
     private:
      /// @name class member functions
      /// @{

      /// access the currently active contact strategy
      CONTACT::Aug::Strategy& get();

      /// access the currently active contact strategy (read-only)
      const CONTACT::Aug::Strategy& get() const;

      /// update the strategy switching conditions
      void switch_update(CONTACT::ParamsInterface& cparams);

      /// run after EvalForce
      void run_post_eval_force(CONTACT::ParamsInterface& cparams);

      /// run after eval_force_stiff
      void run_post_eval_force_stiff(CONTACT::ParamsInterface& cparams);

      /// run after EvalStaticConstraontRHS
      void run_post_eval_static_constraint_rhs(CONTACT::ParamsInterface& cparams);

      /// @}

      //! @name Unsupported derived routines (dead-end)
      //! @{

      //! @name Deprecated methods
      //! @{
      void EvaluateContact(Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff,
          Teuchos::RCP<Epetra_Vector>& feff) override
      {
        FOUR_C_THROW("Deprecated function call!");
      };
      void EvaluateFriction(Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff,
          Teuchos::RCP<Epetra_Vector>& feff) override
      {
        FOUR_C_THROW("Deprecated function call!");
      };
      void build_saddle_point_system(Teuchos::RCP<Core::LinAlg::SparseOperator> kdd,
          Teuchos::RCP<Epetra_Vector> fd, Teuchos::RCP<Epetra_Vector> sold,
          Teuchos::RCP<Core::LinAlg::MapExtractor> dbcmaps, Teuchos::RCP<Epetra_Operator>& blockMat,
          Teuchos::RCP<Epetra_Vector>& blocksol, Teuchos::RCP<Epetra_Vector>& blockrhs) override
      {
        FOUR_C_THROW("Deprecated function call!");
      };
      void update_displacements_and_l_mincrements(
          Teuchos::RCP<Epetra_Vector> sold, Teuchos::RCP<const Epetra_Vector> blocksol) override
      {
        FOUR_C_THROW("Deprecated function call!");
      };
      void Recover(Teuchos::RCP<Epetra_Vector> disi) override
      {
        FOUR_C_THROW("Deprecated function call! Replaced by run_post_compute_x().");
      };
      //! @}

      /*! @name Dead-end for penalty and Uzawa methods (wrong strategy)
       *
       * Please note, that the definition of these functions seems completely unnecessary here.
       * Actually it would be a much better idea to cast the object to the right strategy at the
       * place where it is needed.                                            hiermeier 05/16 */
      //! @{
      double InitialPenalty() override
      {
        FOUR_C_THROW("Wrong strategy!");
        exit(EXIT_FAILURE);
      };
      void InitializeUzawa(Teuchos::RCP<Core::LinAlg::SparseOperator>& kteff,
          Teuchos::RCP<Epetra_Vector>& feff) override
      {
        FOUR_C_THROW("Wrong strategy!");
      };
      void ResetPenalty() override { FOUR_C_THROW("Wrong strategy!"); };
      void ModifyPenalty() override { FOUR_C_THROW("Wrong strategy!"); };
      void update_uzawa_augmented_lagrange() override { FOUR_C_THROW("Wrong strategy!"); };
      void update_constraint_norm(int uzawaiter = 0) override { FOUR_C_THROW("Wrong strategy!"); };
      bool IsPenalty() const override { return false; };
      //! @}
      //! @}

     private:
      // forward declaration of nested classes
      class Switching;
      class PreAsymptoticSwitching;

     private:
      /// wrapped contact solution strategies
      plain_strategy_set strategies_;

      /// lin solvers for the wrapped contact solution strategies
      plain_lin_solver_set lin_solvers_;

      /// interface of the wrapped contact solution strategies
      plain_interface_sets interface_sets_;

      CONTACT::Aug::DataContainer& data_;

      /*----------------------------------------------------------------------*/
      /** \brief nested container for the non-dbc dof handling
       *
       *  The class members contain all entries of the counterparts without any
       *  DBC dofs. */
      struct GlobalNoDbc
      {
        /// constructor
        GlobalNoDbc()
            : slMaMap_(Teuchos::null),
              slMap_(Teuchos::null),
              maMap_(Teuchos::null),
              slMaForce_(Teuchos::null),
              slForce_(Teuchos::null),
              maForce_(Teuchos::null)
        { /*do nothing*/
        }

        /// assemble the maps
        void Assemble(const Epetra_Map& dbcmap, const CONTACT::Aug::DataContainer& data);

        /// handle a parallel redistribution
        void Redistribute(const CONTACT::Aug::DataContainer& data);

        /// reset container members
        void reset(const Epetra_Map& slMaMap, const CONTACT::Aug::DataContainer& data);

        /// slave/master DOF map without DBC
        Teuchos::RCP<Epetra_Map> slMaMap_;

        /// slave DOF map without DBC
        Teuchos::RCP<Epetra_Map> slMap_;

        /// master DOF map without DBC
        Teuchos::RCP<Epetra_Map> maMap_;

        /// slave/master force vector without DBC
        Teuchos::RCP<Epetra_Vector> slMaForce_;

        /// slave force vector without DBC
        Teuchos::RCP<Epetra_Vector> slForce_;

        /// master force vector without DBC
        Teuchos::RCP<Epetra_Vector> maForce_;
      };

      /*----------------------------------------------------------------------*/
      /// nested class which controls the screen output of the surrounding class
      class Output
      {
       public:
        /// init this class
        void initScreenOutput(bool print2screen);

        /// get output screen stream
        std::ostream& oscreen() const;

       private:
        /// output screen stream
        std::ostream* oscreen_ = nullptr;

        /// blackhole stream. This kind of stream will prevent any output if set.
        const Teuchos::RCP<Teuchos::oblackholestream> blackhole_ =
            Teuchos::rcp(new Teuchos::oblackholestream);
      };

      /// non-DBC DOF container
      GlobalNoDbc no_dbc_;

      /// output object for this class
      Output output_;

      /// pointer to the switching strategy
      Teuchos::RCP<Switching> switch_;
    };

    /*--------------------------------------------------------------------------*/
    /// Generic base class for all switching strategies
    class ComboStrategy::Switching
    {
     public:
      /// create and return the desired switching strategy
      static Teuchos::RCP<Switching> Create(ComboStrategy& combo);

      /// constructor
      Switching(ComboStrategy& combo, const Teuchos::ParameterList& p_combo);

      /// destructor
      virtual ~Switching() = default;

      /// return the ID of the given solving strategy
      unsigned Id(enum Inpar::CONTACT::SolvingStrategy sol_type) const;

      /// return the active ID
      virtual unsigned Id() const = 0;

      /// check the switching condition and update the member variables
      virtual void Update(CONTACT::ParamsInterface& cparams, std::ostream& os) = 0;

     protected:
      /// find the id corresponding to \c sol_type
      unsigned find_id(Inpar::CONTACT::SolvingStrategy sol_type) const;

      /// detect strategy types of the wrapped strategies
      void get_strategy_types(
          const plain_strategy_set& strategies, plain_strattype_set& strat_types) const;

     protected:
      /// reference to the combo strategy
      ComboStrategy& combo_;

      /// set of all combined contact strategy types
      plain_strattype_set strat_types_;
    };  // class SwitchingStrategy

    /*--------------------------------------------------------------------------*/
    /** \brief Switching strategy for a preasymptotic/asymptotic switch
     *
     *  This implementation allows the consideration of two different
     *  contact solving strategies. One for the preasymptotic phase and another
     *  one for the asymptotic phase. In this way improved robustness in the
     *  far field and quadratic convergence near the solution can be combined.
     *
     *  \author hiermeier */
    class ComboStrategy::PreAsymptoticSwitching : public ComboStrategy::Switching
    {
     public:
      /// constructor
      PreAsymptoticSwitching(ComboStrategy& combo, const Teuchos::ParameterList& p_combo);

     private:
      /// derived
      unsigned Id() const override;

      /// derived
      void Update(CONTACT::ParamsInterface& cparams, std::ostream& os) override;

      /// check the current penetration (first rough test)
      bool check_penetration(std::ostream& os);

      /** \brief perform checks based on the structure/contact force residual
       *
       *  All the here considered tests are based on the fact that at the solution
       *  the gradient of all active constraints scaled by the active Lagrange multipliers
       *  and the gradient of the structural potential must be collinear and of equal
       *  length at all involved slave/master DOFs. In a slightly weaker form, this can
       *  be used as a reliable indicator to initiate a switch. Therefore, the angle
       *  between the gradients as well as the difference in magnitude must be under a
       *  predefined relative bound.
       *
       *  \author hiermeier */
      bool check_residual(CONTACT::ParamsInterface& cparams, std::ostream& os);

      /// get the penetration bound based on the element edge length
      double get_penetration_bound() const;

      /// return the structural force without DBC DOFs
      Teuchos::RCP<Epetra_Vector> get_structural_force_without_dbc_dofs(
          const CONTACT::ParamsInterface& cparams);

      /** check the relative difference between the Lagrange multiplier contact force
       *  and the structural force at all active contact DOFs. */
      bool check_contact_residual_norm(const Epetra_Vector& str_slmaforce,
          const Epetra_Vector& constr_slmaforce, std::ostream& os) const;

      /** check the angle between the Lagrange multiplier contact force
       *  and the structural force at all active contact DOFs. */
      bool check_angle_between_str_force_and_contact_force(const Epetra_Vector& str_slmaforce,
          const Epetra_Vector& constr_slmaforce, std::ostream& os) const;

      /// check cn-bound (fall-back strategy, which is usually not activated)
      bool check_cn_bound(std::ostream& os) const;

      /// return SlMa forces over all active interfaces
      void get_active_sl_ma_forces(const Epetra_Vector& str_force,
          Teuchos::RCP<Epetra_Vector>& str_slmaforce,
          Teuchos::RCP<Epetra_Vector>& constr_slmaforce) const;

      /// get global dof maps for the active sl/ma force vector
      void get_global_sl_ma_active_force_maps(const Epetra_Vector& slforce,
          const Epetra_Vector& maforce, Teuchos::RCP<Epetra_Map>& gSlActiveForceMap,
          Teuchos::RCP<Epetra_Map>& gMaActiveForceMap) const;

      /// print header
      void print_update_head(std::ostream& os) const;

     private:
      /*----------------------------------------------------------------------*/
      /** container for the max. detected absolute value of all involved
       *  active contact nodes during the preasymptotic and asymptotic phase */
      struct MaxAbsAWGap
      {
        MaxAbsAWGap(PreAsymptoticSwitching& switching)
            : switch_(switching), asymptotic_(1.0e12), pre_asymptotic_(1.0e12){};

        /// update the container variables based on the current solution phase
        inline void Update(const double max_awgap)
        {
          if (switch_.is_asymptotic_)
            asymptotic_ = max_awgap;
          else
            pre_asymptotic_ = max_awgap;
        }

        /// call-back
        const PreAsymptoticSwitching& switch_;

        /// max. detected absolute nodal averaged weighted gap value in the asymptotic phase
        double asymptotic_;

        /// max. detected absolute nodal averaged weighted gap value in the asymptotic phase
        double pre_asymptotic_;
      };

      /// id of the pre-asymptotic solution strategy
      unsigned preasymptotic_id_;

      /// id of the asymptotic solution strategy
      unsigned asymptotic_id_;

      /// which phase is currently active?
      bool is_asymptotic_;

      /// internal container (see container description for more info)
      MaxAbsAWGap maxabsawgap_;
    };
  }  // namespace Aug
}  // namespace CONTACT


FOUR_C_NAMESPACE_CLOSE

#endif
