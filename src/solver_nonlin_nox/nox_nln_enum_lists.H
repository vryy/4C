/*---------------------------------------------------------------------*/
/*!
\file nox_nln_enum_lists.H

\brief Contains ONLY lists of enumerators and is supposed to be included
       into the header files, if necessary and till the C++11 standard is
       available in BACI.

\maintainer Michael Hiermeier

\date May 1, 2016

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_ENUM_LISTS_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_ENUM_LISTS_H_

#include <string>
#include <boost/algorithm/string/predicate.hpp>
#include <NOX_StatusTest_Generic.H>

namespace NOX
{
  namespace NLN
  {
    //! Supported solution type names
    enum SolutionType : int
    {
      sol_unknown,
      sol_structure,           ///< structural problem
      sol_contact,             ///< contact problem
      sol_meshtying,           ///< meshtying problem
      sol_cardiovascular0d,    ///< 0D cardiovascular problem
      sol_lag_pen_constraint,  ///< Lagrange or/and penalty enforced constraint problem
      sol_scatra               ///< scalar transport problem
    };

    //! Map quantity enum to std::string
    inline std::string SolutionType2String(const enum SolutionType type)
    {
      switch (type)
      {
        case sol_structure:
          return "Structure";
        case sol_contact:
          return "Contact";
        case sol_meshtying:
          return "Meshtying";
        case sol_cardiovascular0d:
          return "Cardiovascular0D";
        case sol_lag_pen_constraint:
          return "Lag-Pen-Constraint";
        case sol_scatra:
          return "ScaTra";
        case sol_unknown:
        default:
          return "Unknown Solution Type";
      }
      exit(EXIT_FAILURE);
    };

    //! Map quantity std::string to enum
    inline enum SolutionType String2SolutionType(const std::string& name)
    {
      SolutionType type = sol_unknown;
      if (boost::iequals(name, "structure"))
        type = sol_structure;
      else if (boost::iequals(name, "contact"))
        type = sol_contact;
      else if (boost::iequals(name, "meshtying"))
        type = sol_meshtying;
      else if (boost::iequals(name, "cardiovascular0d"))
        type = sol_cardiovascular0d;
      else if (boost::iequals(name, "lag-pen-constraint"))
        type = sol_lag_pen_constraint;
      else if (boost::iequals(name, "scatra"))
        type = sol_scatra;

      return type;
    };

    //! type of the optimization problem
    enum OptimizationProblemType : int
    {
      opt_unconstrained,          ///< unconstrained optimization problem
      opt_equality_constrained,   ///< pure equality constrained optimization problem
      opt_inequality_constrained  ///< inequality or mixed optimization problem
    };

    /// Correction types
    enum class CorrectionType : int
    {
      vague,          ///< undefined correction type
      soc_automatic,  ///< let the algorithm choose which SOC strategy to follow
      soc_cheap,      ///< cheap SOC step, reuse the old system of equations
      soc_full        ///< do full SOC step in the flavor of the watch dog method
    };

    //! Map second order correction type enum to std::string
    inline std::string CorrectionType2String(const enum CorrectionType type)
    {
      switch (type)
      {
        case CorrectionType::vague:
          return "CorrectionType::vague";
        case CorrectionType::soc_cheap:
          return "CorrectionType::soc_cheap";
        case CorrectionType::soc_full:
          return "CorrectionType::soc_full";
        case CorrectionType::soc_automatic:
          return "CorrectionType::soc_automatic";
        default:
          return "Unknown second order correction Type";
      }
      exit(EXIT_FAILURE);
    };

    namespace LinSystem
    {
      enum class ConditionNumber
      {
        one_norm,         ///< 1-norm condition number
        inf_norm,         ///< inf-norm condition number
        max_min_ev_ratio  ///< condition number as ration of largest and smallest eigenvalue
      };

      inline std::string ConditionNumber2String(const enum ConditionNumber condtype)
      {
        switch (condtype)
        {
          case ConditionNumber::one_norm:
            return "one_norm";
          case ConditionNumber::inf_norm:
            return "inf_norm";
          case ConditionNumber::max_min_ev_ratio:
            return "max_min_ev_ratio";
          default:
            return "unknown_condition_number_type";
        }
      }

      //! supported LinearSystem types
      enum LinearSystemType : int
      {
        linear_system_structure,
        linear_system_structure_meshtying,
        linear_system_structure_contact,
        linear_system_structure_cardiovascular0d,
        linear_system_structure_lag_pen_constraint,
        linear_system_scatra,
        linear_system_undefined
      };

      /// Map linear system type to std::string
      inline std::string LinearSystemType2String(const enum LinearSystemType type)
      {
        switch (type)
        {
          case linear_system_structure:
            return "linear_system_structure";
          case linear_system_structure_contact:
            return "linear_system_structure_contact";
          case linear_system_structure_cardiovascular0d:
            return "linear_system_structure_cardiovascular0d";
          case linear_system_structure_lag_pen_constraint:
            return "linear_system_structure_lag_pen_constraint";
          case linear_system_undefined:
            return "linear_system_undefined";
          case linear_system_scatra:
            return "linear_system_scatra";
          default:
            return "unknown operator type";
        }
        exit(EXIT_FAILURE);
      };

      //! List of types of epetra objects that can be used for the Jacobian and/or Preconditioner.
      enum OperatorType : int
      {
        LinalgSparseOperator,     ///< A LINALG_SparseOperator object.
        LinalgSparseMatrixBase,   ///< A LINALG_SparseMatrixBase object.
        LinalgBlockSparseMatrix,  ///< A LINALG_BlockSparseMatrix object.
        LinalgSparseMatrix        ///< A LINALG_SparseMatrix object.
      };

      /// Map operator type to std::string
      inline std::string OperatorType2String(const enum OperatorType type)
      {
        switch (type)
        {
          case LinalgSparseOperator:
            return "LINALG_SparseOperator";
          case LinalgBlockSparseMatrix:
            return "LINALG_BlockSparseMatrix";
          case LinalgSparseMatrixBase:
            return "LINALG_SparseMatrixBase";
          case LinalgSparseMatrix:
            return "LINALG_SparseMatrix";
          default:
            return "unknown operator type";
        }
        exit(EXIT_FAILURE);
      };
    }  // namespace LinSystem


    namespace Direction
    {
      enum class DefaultStepTest : int
      {
        none,                  ///< unspecified
        volume_change_control  ///< use the volume change test
      };

      inline std::string DefaultStepTest2String(const enum DefaultStepTest test_type)
      {
        switch (test_type)
        {
          case DefaultStepTest::none:
            return "DefaultStepTest::none";
          case DefaultStepTest::volume_change_control:
            return "DefaultStepTest::volume_change_control";
          default:
            return "UNKNOWN DefaultStepTest enumerator";
        }
      }
    }  // namespace Direction

    namespace MeritFunction
    {
      //! order of the linearization term
      enum LinOrder : int
      {
        //! all orders
        linorder_all,
        //! first order linearization
        linorder_first,
        //! second order linearization terms
        linorder_second
      };

      //! type of the linearization term, i.e. with respect to which quantity the linearization was
      //! performed
      enum LinType : int
      {
        //! linearization with respect to primary and Lagrange multiplier degrees of freedom
        lin_wrt_all_dofs,
        //! linearization with repsect to primary degrees of freedom
        lin_wrt_primary_dofs,
        //! linearization with respect to Lagrange multiplier degrees of freedom
        lin_wrt_lagrange_multiplier_dofs,
        //! linearization with respect to mixed degrees of freedom
        lin_wrt_mixed_dofs
      };

      //! merit function names
      enum MeritFctName : int
      {
        mrtfct_sum_of_squares,          //!< sum of squares merit function
        mrtfct_lagrangian,              //!< lagrangian merit function
        mrtfct_lagrangian_active,       //!< lagrangian considering only the active contributions
        mrtfct_infeasibility_two_norm,  //!< infeasibility merit function
        mrtfct_infeasibility_two_norm_active, /*!< infeasibility merit function,
                                               *  considering only the active contributions */
        mrtfct_energy,                        //!< representative energy value
        mrtfct_vague                          //!< unspecified
      };

      inline std::string MeritFuncName2String(const enum MeritFctName mrt_func)
      {
        switch (mrt_func)
        {
          case mrtfct_sum_of_squares:
            return "mrtfct_sum_of_squares";
          case mrtfct_lagrangian:
            return "mrtfct_lagrangian";
          case mrtfct_lagrangian_active:
            return "mrtfct_lagrangian_active";
          case mrtfct_infeasibility_two_norm:
            return "mrtfct_infeasibility_two_norm";
          case mrtfct_infeasibility_two_norm_active:
            return "mrtfct_infeasibility_two_norm_active";
          case mrtfct_energy:
            return "mrtfct_energy";
          case mrtfct_vague:
            return "mrtfct_vague";
          default:
            return "INVALID_merit_function_name";
        }
      }

    }  // namespace MeritFunction

    namespace StatusTest
    {
      //! map status type (outer status test) to std::string
      inline std::string StatusType2String(const enum NOX::StatusTest::StatusType stype)
      {
        switch (stype)
        {
          case NOX::StatusTest::Unevaluated:
            return "Unevaluated";
          case NOX::StatusTest::Unconverged:
            return "Unconverged";
          case NOX::StatusTest::Converged:
            return "Converged";
          case NOX::StatusTest::Failed:
            return "Failed";
          default:
            return "Unknown StatusType";
        }
      };

      //! Supported quantity names for distinguished status tests
      enum QuantityType : int
      {
        quantity_unknown,         ///< unknown quantity (dummy)
        quantity_structure,       ///< check structural quantities
        quantity_contact_normal,  ///< check (semi-smooth) contact quantities (normal/frictionless)
        quantity_contact_friction,    ///< check (semi-smooth) contact quantities (frictionless)
        quantity_meshtying,           ///< check meshtying quantities
        quantity_cardiovascular0d,    ///< check 0d cardiovascular quantities
        quantity_lag_pen_constraint,  ///< check Lagrange/penalty enforced constraint quantites
        quantity_plasticity,          ///< check semi-smooth plasticity
        quantity_pressure,            ///< check pressure dofs
        quantity_eas,                 ///< check eas dofs
        quantity_levelset_reinit      ///< check levelset reinitialization
      };

      /// Map quantity name to std::string
      inline std::string QuantityType2String(const enum QuantityType type)
      {
        switch (type)
        {
          case quantity_structure:
            return "STRUCTURE";
          case quantity_contact_normal:
            return "CONTACT-NORMAL";
          case quantity_contact_friction:
            return "CONTACT-FRICTION";
          case quantity_meshtying:
            return "MESHTYING";
          case quantity_cardiovascular0d:
            return "CARDIOVASCULAR0D";
          case quantity_lag_pen_constraint:
            return "LAG-PEN-CONSTRAINT";
          case quantity_plasticity:
            return "PLASTICITY";
          case quantity_pressure:
            return "PRESSURE";
          case quantity_eas:
            return "EAS";
          case quantity_levelset_reinit:
            return "LEVELSET-REINIT";
          case quantity_unknown:
          default:
            return "unknown quantity type";
        }
        exit(EXIT_FAILURE);
      };

      /// Map std::string to quantity type
      inline QuantityType String2QuantityType(const std::string& name)
      {
        QuantityType type = quantity_unknown;

        // -- case insensitive comparison --
        if (boost::iequals(name, "structure"))
          type = quantity_structure;
        else if (boost::iequals(name, "contact-normal"))
          type = quantity_contact_normal;
        else if (boost::iequals(name, "contact-friction"))
          type = quantity_contact_friction;
        else if (boost::iequals(name, "meshtying"))
          type = quantity_meshtying;
        else if (boost::iequals(name, "cardiovascular0d"))
          type = quantity_cardiovascular0d;
        else if (boost::iequals(name, "lag-pen-constraint"))
          type = quantity_lag_pen_constraint;
        else if (boost::iequals(name, "plasticity"))
          type = quantity_plasticity;
        else if (boost::iequals(name, "pressure"))
          type = quantity_pressure;
        else if (boost::iequals(name, "eas"))
          type = quantity_eas;
        else if (boost::iequals(name, "levelset-reinit"))
          type = quantity_levelset_reinit;

        return type;
      };
    }  // namespace StatusTest
  }    // namespace NLN
}  // namespace NOX



#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_ENUM_LISTS_H_ */
