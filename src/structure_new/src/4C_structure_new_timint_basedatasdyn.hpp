/*-----------------------------------------------------------*/
/*! \file

\brief Structural dynamics data container for the structural (time)
       integration


\level 3

*/
/*-----------------------------------------------------------*/


#ifndef FOUR_C_STRUCTURE_NEW_TIMINT_BASEDATASDYN_HPP
#define FOUR_C_STRUCTURE_NEW_TIMINT_BASEDATASDYN_HPP

#include "4C_config.hpp"

#include "4C_inpar_structure.hpp"
#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_utils_exceptions.hpp"

#include <NOX_Abstract_Vector.H>

#include <set>

// forward declaration
namespace Teuchos
{
  class Time;
};

FOUR_C_NAMESPACE_OPEN

namespace CORE::LINALG
{
  class Solver;
}  // namespace CORE::LINALG
namespace DRT
{
  class Discretization;
}  // namespace DRT
namespace CORE::GEO
{
  namespace MESHFREE
  {
    class BoundingBox;
  }  // namespace MESHFREE
}  // namespace CORE::GEO
namespace STR
{
  namespace MODELEVALUATOR
  {
    class Generic;
  }  // namespace MODELEVALUATOR
  namespace TIMINT
  {
    /** \brief Structural dynamics data container for the structural (time) integration
     *
     * This data container holds everything, which refers directly to the
     * control mechanisms of the structural dynamics algorithms. The input parameters
     * are transformed into enumerators. Furthermore, the maximum time and step number
     * can be found here.
     *
     * \author Michael Hiermeier */
    class BaseDataSDyn
    {
     public:
      /// constructor
      BaseDataSDyn();

      /// destructor
      virtual ~BaseDataSDyn() = default;

      /// initialize class variables (already existing)
      virtual void Init(const Teuchos::RCP<DRT::Discretization> discret,
          const Teuchos::ParameterList& sDynParams, const Teuchos::ParameterList& xparams,
          const Teuchos::RCP<std::set<enum INPAR::STR::ModelType>> modeltypes,
          const Teuchos::RCP<std::set<enum INPAR::STR::EleTech>> eletechs,
          const Teuchos::RCP<
              std::map<enum INPAR::STR::ModelType, Teuchos::RCP<CORE::LINALG::Solver>>>
              linsolvers);

      /// setup new class variables
      virtual void Setup();

     protected:
      /// get the indicator state
      inline const bool& is_init() const { return isinit_; };

      /// get the indicator state
      inline const bool& is_setup() const { return issetup_; };

      /// Check if Init() has been called, yet.
      inline void check_init() const
      {
        FOUR_C_ASSERT(is_init(), "Init() has not been called, yet!");
      }

      /// Check if Init() and Setup() have been called, yet.
      inline void check_init_setup() const
      {
        FOUR_C_ASSERT(is_init() and is_setup(), "Call Init() and Setup() first!");
      }

     public:
      /// @name Get general control parameters (read only access)
      ///@{

      /// Returns final time \f$t_\text{fin}\f$
      const double& GetTimeMax() const
      {
        check_init_setup();
        return timemax_;
      };

      /// Returns final time step \f$N\f$
      const int& GetStepMax() const
      {
        check_init_setup();
        return stepmax_;
      };

      bool is_restarting_initial_state() const
      {
        check_init_setup();
        return isrestarting_initial_state_;
      }

      /// Returns timer for solution technique
      Teuchos::RCP<const Teuchos::Time> GetTimer() const
      {
        check_init_setup();
        return timer_;
      };

      /// Returns dynamic type
      const enum INPAR::STR::DynamicType& GetDynamicType() const
      {
        check_init_setup();
        return dyntype_;
      };

      /// Get type of shell thickness scaling
      const enum INPAR::STR::StcScale& GetSTCAlgoType() const
      {
        check_init_setup();
        return stcscale_;
      };

      /// Get number of layers for multilayered shell thickness scaling case
      const int& GetSTCLayer() const
      {
        check_init_setup();
        return stclayer_;
      };

      /// Returns minimal non-linear iteration number
      const int& GetIterMin() const
      {
        check_init_setup();
        return itermin_;
      };

      /// Returns maximal non-linear iteration number
      const int& GetIterMax() const
      {
        check_init_setup();
        return itermax_;
      };

      /// Returns true if the external load should be linearized
      const bool& GetLoadLin() const
      {
        check_init_setup();
        return loadlin_;
      }

      // Return time until the prestressing algorthm should be applied
      double GetPreStressTime() const
      {
        check_init_setup();
        return prestresstime_;
      }

      /// Returns the tolerance for the displacements during prestressing
      double get_pre_stress_displacement_tolerance() const
      {
        check_init_setup();
        return prestress_displacement_tolerance_;
      }

      /// Returns the minimum number of load steps during prestressing
      [[nodiscard]] int get_pre_stress_minimum_number_of_load_steps() const
      {
        check_init_setup();
        return prestress_min_number_of_load_steps_;
      }

      /// Returns prestress type
      const enum INPAR::STR::PreStress& GetPreStressType() const
      {
        check_init_setup();
        return prestresstype_;
      };

      /// Returns predictor type
      const enum INPAR::STR::PredEnum& GetPredictorType() const
      {
        check_init_setup();
        return predtype_;
      };

      /// Returns nonlinear solver type
      const enum INPAR::STR::NonlinSolTech& GetNlnSolverType() const
      {
        check_init_setup();
        return nlnsolvertype_;
      };

      /// Returns the divergence action
      /// Short: What to do if the non-linear solver fails.
      const enum INPAR::STR::DivContAct& GetDivergenceAction() const
      {
        check_init_setup();
        return divergenceaction_;
      };

      /// Returns mid-time energy type
      enum INPAR::STR::MidAverageEnum get_mid_time_energy_type() const
      {
        check_init_setup();
        return mid_time_energy_type_;
      }

      /// Returns number of times you want to halve your timestep in case nonlinear solver diverges
      const int& get_max_div_con_refine_level() const
      {
        check_init_setup();
        return maxdivconrefinementlevel_;
      };

      /// Returns nox parameters
      const Teuchos::ParameterList& GetNoxParams() const
      {
        check_init_setup();
        return *noxparams_;
      }

      /// Returns loca parameters
      const Teuchos::ParameterList& GetLocaParams() const
      {
        check_init_setup();
        return *locaparams_;
      }

      /// Returns the inital pseudo time step for the PTC method
      const double& get_initial_ptc_pseudo_time_step() const { return ptc_delta_init_; }
      ///@}

      /// @name Get mutable linear solver variables (read only access)
      ///@{

      /// Returns linear solvers pointer
      const std::map<enum INPAR::STR::ModelType, Teuchos::RCP<CORE::LINALG::Solver>>&
      GetLinSolvers() const
      {
        check_init_setup();
        return *linsolvers_;
      }
      ///@}


      /// @name Get damping control parameters (read only access)
      ///@{
      /// Returns damping type
      const enum INPAR::STR::DampKind& GetDampingType() const
      {
        check_init_setup();
        return damptype_;
      };

      /// Returns damping factor for stiffness \f$c_\text{K}\f$
      const double& get_damping_stiffness_factor() const
      {
        check_init_setup();
        return dampk_;
      };

      /// Returns damping factor for mass \f$c_\text{M}\f$
      const double& get_damping_mass_factor() const
      {
        check_init_setup();
        return dampm_;
      };
      ///@}

      /// @name Get mass and inertia control parameters (read only access)
      ///@{
      /// Returns mass linearization type
      const enum INPAR::STR::MassLin& GetMassLinType() const
      {
        check_init_setup();
        return masslintype_;
      };

      const bool& IsMassLumping() const
      {
        check_init_setup();
        return lumpmass_;
      }

      const bool& NeglectInertia() const
      {
        check_init_setup();
        return neglectinertia_;
      }
      ///@}

      /// @name Get model evaluator control parameters (read only access)
      ///@{
      /// Returns types of the current models
      const std::set<enum INPAR::STR::ModelType>& GetModelTypes() const
      {
        check_init_setup();
        return *modeltypes_;
      };

      /// Returns current active element technologies
      const std::set<enum INPAR::STR::EleTech>& get_element_technologies() const
      {
        check_init_setup();
        return *eletechs_;
      };

      /// check if the given model type is active.
      bool HaveModelType(const INPAR::STR::ModelType& modeltype) const;

      /// check if the given element technology is active.
      bool HaveEleTech(const INPAR::STR::EleTech& eletech) const;
      ///@}

      /// @name Get model specific data container
      ///@{
      /// Returns periodic bounding box object (read access)
      Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> const& get_periodic_bounding_box() const
      {
        check_init_setup();
        return periodic_boundingbox_;
      }
      ///@}

      /// @name Get the different status test control parameters (read only)
      ///@{
      /// Returns the STR vector norm type
      const enum INPAR::STR::VectorNorm& GetNormType() const
      {
        check_init_setup();
        return normtype_;
      }

      /// Returns the NOX normtype
      const enum ::NOX::Abstract::Vector::NormType& GetNoxNormType() const
      {
        check_init_setup();
        return nox_normtype_;
      }

      /// Return random factor for timestep in case nonlinear solver diverges
      const double& get_random_time_step_factor() const
      {
        check_init_setup();
        return rand_tsfac_;
      }


      /// Return level of refinement in case of divercontype_ == adapt_step
      const int& get_div_con_refine_level() const
      {
        check_init_setup();
        return divconrefinementlevel_;
      }


      /// Return number of fine steps in case of in case of divercontype_ == adapt_step
      const int& get_div_con_num_fine_step() const
      {
        check_init_setup();
        return divconnumfinestep_;
      }

      /// @name Get residual and increment related parameters
      ///@{
      /// Returns the combination type of the two quantities
      enum INPAR::STR::BinaryOp GetResIncrComboType(
          const enum NOX::NLN::StatusTest::QuantityType& qtype_res,
          const enum NOX::NLN::StatusTest::QuantityType& qtype_incr) const;
      ///@}

      /// @name Get residual related parameters
      ///@{
      /// Returns the tolerance values for the different quantities
      double GetResTolerance(const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

      /// Returns the tolerance type of the different quantities
      enum INPAR::STR::ConvNorm GetResToleranceType(
          const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

      /// Returns the combination type of the different quantities
      enum INPAR::STR::BinaryOp GetResComboType(
          const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

      enum INPAR::STR::BinaryOp GetResComboType(
          const enum NOX::NLN::StatusTest::QuantityType& qtype_1,
          const enum NOX::NLN::StatusTest::QuantityType& qtype_2) const;
      ///@}

      /// @name Get increment related parameters
      ///@{
      /// Returns the tolerance values for the different quantities
      double GetIncrTolerance(const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

      /// Returns the tolerance type of the different quantities
      enum INPAR::STR::ConvNorm get_incr_tolerance_type(
          const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

      enum INPAR::STR::ConvNorm get_incr_tolerance_type(
          const enum NOX::NLN::StatusTest::QuantityType& qtype_1,
          const enum NOX::NLN::StatusTest::QuantityType& qtype_2) const;

      /// Returns the combination type of the different quantities
      enum INPAR::STR::BinaryOp GetIncrComboType(
          const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

      enum INPAR::STR::BinaryOp GetIncrComboType(
          const enum NOX::NLN::StatusTest::QuantityType& qtype_1,
          const enum NOX::NLN::StatusTest::QuantityType& qtype_2) const;
      ///@}
      ///@}

      /// @name Get mutable general control parameters (read and write access)
      ///@{

      /// Returns final time \f$t_\text{fin}\f$
      double& GetTimeMax()
      {
        check_init_setup();
        return timemax_;
      };

      /// Returns final time step \f$N\f$
      int& GetStepMax()
      {
        check_init_setup();
        return stepmax_;
      };

      /// Returns timer for solution technique
      Teuchos::RCP<Teuchos::Time>& GetTimer()
      {
        check_init_setup();
        return timer_;
      };

      /// Returns minimal non-linear iteration number
      int& GetIterMin()
      {
        check_init_setup();
        return itermin_;
      };

      /// Returns maximal non-linear iteration number
      int& GetIterMax()
      {
        check_init_setup();
        return itermax_;
      };

      double& GetPreStressTime()
      {
        check_init_setup();
        return prestresstime_;
      }

      double& get_pre_stress_displacement_tolerance()
      {
        check_init_setup();
        return prestress_displacement_tolerance_;
      }

      [[nodiscard]] int& get_pre_stress_minimum_number_of_load_steps()
      {
        check_init_setup();
        return prestress_min_number_of_load_steps_;
      }

      /// Returns prestress type
      enum INPAR::STR::PreStress& GetPreStressType()
      {
        check_init_setup();
        return prestresstype_;
      };

      /// Returns predictor type
      enum INPAR::STR::PredEnum& GetPredictorType()
      {
        check_init_setup();
        return predtype_;
      };

      /// Returns dynamic type
      enum INPAR::STR::DynamicType& GetDynamicType()
      {
        check_init_setup();
        return dyntype_;
      };

      /// Returns nonlinear solver type
      enum INPAR::STR::NonlinSolTech& GetNlnSolverType()
      {
        check_init_setup();
        return nlnsolvertype_;
      };

      /// Returns mid-time energy type
      enum INPAR::STR::MidAverageEnum& get_mid_time_energy_type()
      {
        check_init_setup();
        return mid_time_energy_type_;
      }

      /// Returns the divergence action
      /// Short: What to do if the non-linear solver fails.
      enum INPAR::STR::DivContAct& GetDivergenceAction()
      {
        check_init_setup();
        return divergenceaction_;
      };

      /// Returns nox parameters
      Teuchos::ParameterList& GetNoxParams()
      {
        check_init_setup();
        return *noxparams_;
      }

      /// Returns nox parameters
      Teuchos::RCP<Teuchos::ParameterList> GetNoxParamsPtr()
      {
        check_init_setup();
        return noxparams_;
      }

      /// Returns loca parameters
      Teuchos::ParameterList& GetLocaParams()
      {
        check_init_setup();
        return *locaparams_;
      }

      /// Return random factor for time step in case nonlinear solver diverges
      double& get_random_time_step_factor()
      {
        check_init_setup();
        return rand_tsfac_;
      }

      /// Return level of refinement in case of divercontype_ == adapt_step
      int& get_div_con_refine_level()
      {
        check_init_setup();
        return divconrefinementlevel_;
      }

      /// Return number of fine steps in case of in case of divercontype_ == adapt_step
      int& get_div_con_num_fine_step()
      {
        check_init_setup();
        return divconnumfinestep_;
      }

      ///@}

      /// @name Get mutable linear solver variables (read and write access)
      ///@{
      /// Returns linear solvers pointer
      Teuchos::RCP<std::map<enum INPAR::STR::ModelType, Teuchos::RCP<CORE::LINALG::Solver>>>
      GetLinSolversPtr()
      {
        check_init_setup();
        return linsolvers_;
      }

      /// Returns linear solvers pointer
      std::map<enum INPAR::STR::ModelType, Teuchos::RCP<CORE::LINALG::Solver>>& GetLinSolvers()
      {
        check_init_setup();
        return *linsolvers_;
      }
      ///@}


      /// @name Get mutable damping control parameters (read and write access)
      ///@{
      /// Returns damping type
      enum INPAR::STR::DampKind& GetDampingType()
      {
        check_init_setup();
        return damptype_;
      };

      /// Returns damping factor for stiffness \f$c_\text{K}\f$
      double& get_damping_stiffness_factor()
      {
        check_init_setup();
        return dampk_;
      };

      /// Returns damping factor for mass \f$c_\text{M}\f$
      double& get_damping_mass_factor()
      {
        check_init_setup();
        return dampm_;
      };
      ///@}

      /// @name Get mutable mass and inertia control parameters (read and write access)
      ///@{
      /// Returns mass linearization type
      enum INPAR::STR::MassLin& GetMassLinType()
      {
        check_init_setup();
        return masslintype_;
      };
      ///@}

      /// @name Get model evaluator control parameters (read and write access)
      ///@{
      /// Returns types of the current active models
      std::set<enum INPAR::STR::ModelType>& GetModelTypes()
      {
        check_init_setup();
        return *modeltypes_;
      };

      /// Returns the current active element technologies
      std::set<enum INPAR::STR::EleTech>& get_element_technologies()
      {
        check_init_setup();
        return *eletechs_;
      };

      /// Return a pointer to the coupling model evaluator
      const Teuchos::RCP<STR::MODELEVALUATOR::Generic>& CouplingModelPtr()
      {
        return coupling_model_ptr_;
      }
      ///@}

      /// @name Get the initial displacement method
      ///@{
      /// Returns type of initial displacement
      enum INPAR::STR::InitialDisp GetInitialDisp() const
      {
        check_init_setup();
        return initial_disp_;
      };

      /// Returns the function index to initialize the displacement
      int StartFuncNo() const
      {
        check_init_setup();
        return start_func_no_;
      }
      ///@}

     protected:
      /** \brief returns the structural dynamics parameter list (internal use only!)
       *
       * Do not touch this. It should be used only in derived Setup routines.
       * Do not call it from outside! */
      const Teuchos::ParameterList& get_s_dyn_params() const
      {
        check_init();
        return *sdynparams_ptr_;
      }

     protected:
      /// @name variables for internal use only
      ///@{
      /// init flag
      bool isinit_;

      /// setup flag
      bool issetup_;
      ///@}

     private:
      /// @name General control parameters
      ///@{

      /// final time \f$t_\text{fin}\f$
      double timemax_;

      /// final time step \f$N\f$
      int stepmax_;

      /// flag indicating that if the initial state needs to be re-calculated
      bool isrestarting_initial_state_;

      ///@}

      /// @name Time measurement
      ///@{
      /// timer for solution technique
      Teuchos::RCP<Teuchos::Time> timer_;

      ///@}

      /// @name Damping control parameters
      /// Rayleigh damping means \f${C} = c_\text{K} {K} + c_\text{M} {M}\f$
      ///@{
      /// damping type
      enum INPAR::STR::DampKind damptype_;

      /// damping factor for stiffness \f$c_\text{K}\f$
      double dampk_;

      /// damping factor for mass \f$c_\text{M}\f$
      double dampm_;
      ///@}

      /// @name Mass and inertia control parameters
      ///@{
      /// have inertia forces to be linearized?
      enum INPAR::STR::MassLin masslintype_;

      /// lump the mass matrix?
      bool lumpmass_;

      /// neglect inertia?
      bool neglectinertia_;
      ///@}

      /// @name Model evaluator control parameters
      ///@{

      /// current active model types for the model evaluator
      Teuchos::RCP<std::set<enum INPAR::STR::ModelType>> modeltypes_;

      /// current active element technologies
      Teuchos::RCP<std::set<enum INPAR::STR::EleTech>> eletechs_;

      /// pointer to the coupling model evaluator object
      Teuchos::RCP<STR::MODELEVALUATOR::Generic> coupling_model_ptr_;

      ///@}

      /// @name implicit and explicit time integrator parameters
      ///@{
      /// dynamic type
      enum INPAR::STR::DynamicType dyntype_;

      /// scaled thickness conditioning type (necessary for shells)
      enum INPAR::STR::StcScale stcscale_;

      /// number of layers for multilayered case
      int stclayer_;

      ///@}

      /// @name implicit time integrator parameters
      ///@{
      /// minimal non-linear iteration number
      int itermin_;

      /// maximal non-linear iteration number
      int itermax_;

      /// linearization of external follower load in Newton
      bool loadlin_;

      /// Time until the prestressing algorithm should be applied
      double prestresstime_;

      /// Tolerance for the residual displacements during prestressing
      double prestress_displacement_tolerance_;

      /// Minimum number of load steps during prestressing
      int prestress_min_number_of_load_steps_;

      /// type of pre-stressing
      enum INPAR::STR::PreStress prestresstype_;

      /// type of the predictor
      enum INPAR::STR::PredEnum predtype_;

      /// type of nonlinear solver
      enum INPAR::STR::NonlinSolTech nlnsolvertype_;

      /// action to be performed when the non-linear solver diverges
      enum INPAR::STR::DivContAct divergenceaction_;

      /// mid-time energy type
      enum INPAR::STR::MidAverageEnum mid_time_energy_type_;

      /// how often you want to half your timestep until you give up
      int maxdivconrefinementlevel_;

      /// nox parameters list
      Teuchos::RCP<Teuchos::ParameterList> noxparams_;

      /// loca parameter list
      Teuchos::RCP<Teuchos::ParameterList> locaparams_;

      /// initial pseudo time step for the pseudo transient continuation (PTC) method
      double ptc_delta_init_;
      ///@}

      /// @name linear solver variables
      ///@{

      /// pointer to the linear solvers map
      Teuchos::RCP<std::map<enum INPAR::STR::ModelType, Teuchos::RCP<CORE::LINALG::Solver>>>
          linsolvers_;
      ///@}

      /// @name status test control parameters
      ///@{

      /// vector norm type for the status/convergence test
      enum INPAR::STR::VectorNorm normtype_;

      /// vector norm type for the status/convergence test
      enum ::NOX::Abstract::Vector::NormType nox_normtype_;

      /// @name primary variables
      ///@{
      /// tolerance residual displacements
      double tol_disp_incr_;

      /// tolerance type for the convergence check of the displacement vector
      enum INPAR::STR::ConvNorm toltype_disp_incr_;

      /// tolerance force residual
      double tol_fres_;

      /// tolerance type for the convergence check of the force residual
      enum INPAR::STR::ConvNorm toltype_fres_;

      /// tolerance pressure residual
      double tol_pres_;

      /// tolerance type for the convergence check of the pressure residual
      enum INPAR::STR::ConvNorm toltype_pres_;

      /// tolerance incompressible residual/ residual pressure forces
      double tol_inco_;

      /// tolerance type for the convergence check if the incompressible residual/ residual pressure
      /// forces
      enum INPAR::STR::ConvNorm toltype_inco_;

      /// tolerance plasticity residual
      double tol_plast_res_;

      /// tolerance type for the convergence check of the plasticity residual
      enum INPAR::STR::ConvNorm toltype_plast_res_;

      /// tolerance plasticity increment
      double tol_plast_incr_;

      /// tolerance type for the convergence check of the plasticity increment
      enum INPAR::STR::ConvNorm toltype_plast_incr_;

      /// tolerance EAS residual
      double tol_eas_res_;

      /// tolerance type for the convergence check of the EAS residual
      enum INPAR::STR::ConvNorm toltype_eas_res_;

      /// tolerance EAS increment
      double tol_eas_incr_;

      /// tolerance type for the convergence check of the EAS increment
      enum INPAR::STR::ConvNorm toltype_eas_incr_;

      /// type of combination of the displacement and the pressure test
      enum INPAR::STR::BinaryOp normcombo_disp_pres_;

      /// type of combination of the force and the pressure residual test
      enum INPAR::STR::BinaryOp normcombo_fres_inco_;

      /// type of combination of the force and the EAS residual
      enum INPAR::STR::BinaryOp normcombo_fres_eas_res_;

      /// type of combination of the displacement and the EAS increment
      enum INPAR::STR::BinaryOp normcombo_disp_eas_incr_;

      /// type of combination of the force and the plasticity residual
      enum INPAR::STR::BinaryOp normcombo_fres_plast_res_;

      /// type of combination of the displacement and the plasticity increment
      enum INPAR::STR::BinaryOp normcombo_disp_plast_incr_;

      /// type of combination of the force and the displacement test
      enum INPAR::STR::BinaryOp normcombo_fres_disp_;
      ///@}

      /// @name constraint variables
      ///@{

      /// tolerance type for the convergence check of the 0D cardiovascular residual
      enum INPAR::STR::ConvNorm toltype_cardvasc0d_res_;

      /// tolerance 0D cardiovascular residual
      double tol_cardvasc0d_res_;

      /// tolerance type for the convergence check of the 0D cardiovascular dof increment
      enum INPAR::STR::ConvNorm toltype_cardvasc0d_incr_;

      /// tolerance 0D cardiovascular dof increment
      double tol_cardvasc0d_incr_;

      /// tolerance type for the convergence check of the constraint residual
      enum INPAR::STR::ConvNorm toltype_constr_res_;

      /// tolerance constraint residual
      double tol_constr_res_;

      /// tolerance type for the convergence check of the constraint LM increment
      enum INPAR::STR::ConvNorm toltype_constr_incr_;

      /// tolerance constraint LM increment
      double tol_constr_incr_;

      /// tolerance type for the convergence check of the contact residual
      enum INPAR::STR::ConvNorm toltype_contact_res_;

      /// tolerance contact constraint residual
      double tol_contact_res_;

      /// tolerance type  for the convergence check of the contact lagrange increment
      enum INPAR::STR::ConvNorm toltype_contact_lm_incr_;

      /// tolerance contact lagrange increment
      double tol_contact_lm_incr_;

      /// type of combination of the force and the contact residual
      enum INPAR::STR::BinaryOp normcombo_fres_contact_res_;

      /// type of combination of the displacement and the contact lagrange multiplier increment test
      enum INPAR::STR::BinaryOp normcombo_disp_contact_lm_incr_;

      /// type of combination of the force and the 0D cardiovascular model residual
      enum INPAR::STR::BinaryOp normcombo_fres_cardvasc0d_res_;

      /// type of combination of the displacement and the 0D cardiovascular model dof increment test
      enum INPAR::STR::BinaryOp normcombo_disp_cardvasc0d_incr_;

      /// type of combination of the force and the constraint residual
      enum INPAR::STR::BinaryOp normcombo_fres_constr_res_;

      /// type of combination of the displacement and the 0D cardiovascular model dof increment test
      enum INPAR::STR::BinaryOp normcombo_disp_constr_incr_;

      /// random factor for modifying time-step size in case this way of continuing non-linear
      /// iteration was chosen
      double rand_tsfac_;

      /// number of refinement level in case of divercontype_ == adapt_step
      int divconrefinementlevel_;

      /// number of converged time steps on current refinement level in case of divercontype_ ==
      /// adapt_step
      int divconnumfinestep_;

      ///@}

      /// @name initial displacement variables
      ///@{

      /// type of initial displacement initialization
      enum INPAR::STR::InitialDisp initial_disp_;
      int start_func_no_;

      ///@}

      /// @name pointer to model specific data container
      ///@{

      /// pointer to the linear solvers map
      Teuchos::RCP<CORE::GEO::MESHFREE::BoundingBox> periodic_boundingbox_;
      ///@}

      ///@}
      /** \brief Ptr to the structural dynamics parameter list
       *
       * Do not touch this. It should be used only in derived Setup routines. Do not call it from
       * outside! */
      Teuchos::RCP<const Teuchos::ParameterList> sdynparams_ptr_;
    };  // class BaseDataSDyn

    /** \brief Generalized alpha structural dynamics data container
     *
     * This data container is derived from the standard structural dynamics data container
     * and contains some special GenAlpha input parameters.
     *
     * \author Michael Hiermeier */
    class GenAlphaDataSDyn : public BaseDataSDyn
    {
     public:
      //! contructor
      GenAlphaDataSDyn();

      //! Setup function [derived]
      void Setup() override;

     public:
      //! @name Read-only accessors
      //!@{

      //! returns the mid-average type (for more information see MidAverageEnum)
      const enum INPAR::STR::MidAverageEnum& GetMidAverageType() const
      {
        check_init_setup();
        return midavg_;
      };

      //! Return time integration parameter \f$\beta\f$
      const double& GetBeta() const
      {
        check_init_setup();
        return beta_;
      };

      //! Return time integration parameter \f$\gamma\f$
      const double& GetGamma() const
      {
        check_init_setup();
        return gamma_;
      };

      //! Return time integration parameter \f$\alpha_f\f$
      const double& GetAlphaF() const
      {
        check_init_setup();
        return alphaf_;
      };

      //! Return time integration parameter \f$\alpha_m\f$
      const double& GetAlphaM() const
      {
        check_init_setup();
        return alpham_;
      };

      //! Return spectral radius \f$\rho_\infty\f$
      const double& GetRhoInf() const
      {
        check_init_setup();
        return rhoinf_;
      };

      //!@}

     private:
      //! mid-average type more at MidAverageEnum
      enum INPAR::STR::MidAverageEnum midavg_;

      //! Time integration parameter \f$\beta \in (0,1/2]\f$
      double beta_;

      //! Time integration parameter \f$\gamma \in (0,1]\f$
      double gamma_;

      //! Time integation parameter \f$\alpha_f \in [0,1)\f$
      double alphaf_;

      //! Time integration parameter \f$\alpha_m \in [-1,1)\f$
      double alpham_;

      //! Spectral radius \f$\rho_\infty \in [0,1]\f$
      double rhoinf_;
    };

    /** \brief One-step Theta structural dynamics data container
     *
     * This data container is derived from the standard structural dynamics data container
     * and contains some special OneStepTheta input parameters.
     *
     * \author Michael Hiermeier */
    class OneStepThetaDataSDyn : public BaseDataSDyn
    {
     public:
      //! contructor
      OneStepThetaDataSDyn();

      //! Setup function [derived]
      void Setup() override;

     public:
      //! @name Read-only accessors
      //! @{

      //! Return time integration parameter \f$\theta\f$
      const double& get_theta() const
      {
        check_init_setup();
        return theta_;
      };

      //! @}

     private:
      //! Time integration parameter \f$\theta \in [0,1]\f$
      double theta_;
    };

    /** \brief Forward Euler structural dynamics data container
     *
     * This data container is derived from the standard structural dynamics data container
     * and contains some special Forward Euler input parameters.
     *
     */
    class ExplEulerDataSDyn : public BaseDataSDyn
    {
     public:
      //! contructor
      ExplEulerDataSDyn();

      //! Setup function [derived]
      void Setup() override;

     public:
      //! @name Read-only accessors
      //! @{

      //! Return time integration parameter \f$\theta\f$
      bool get_modified_forward_euler() const
      {
        check_init_setup();
        return modexpleuler_;
      };

      //! @}

     private:
      //! Modified forward Euler flag
      bool modexpleuler_;
    };

  }  // namespace TIMINT
}  // namespace STR


FOUR_C_NAMESPACE_CLOSE

#endif
