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

namespace NOX {
namespace NLN {
  //! Supported solution type names
  enum SolutionType {
    sol_unknown,
    // structural quantities
    sol_structure,
    // contact quantities
    sol_contact,
    // meshtying quantities
    sol_meshtying,
    // 0D cardiovascular quantities
    sol_cardiovascular0d,
    // Lagrange or/and penalty enforced constraint quantities
    sol_lag_pen_constraint
  };

  //! Map quantity enum to std::string
  static inline const std::string SolutionType2String(const enum SolutionType& type)
  {
    switch (type)
    {
    case sol_structure :
      return "Structure";
      break;
    case sol_contact :
      return "Contact";
      break;
    case sol_meshtying :
      return "Meshtying";
      break;
    case sol_cardiovascular0d:
      return "Cardiovascular0D";
      break;
    case sol_lag_pen_constraint:
      return "Lag-Pen-Constraint";
      break;
    case sol_unknown:
    default :
      return "Unknown Solution Type";
      break;
    }
    return "";
  };

  //! Map quantity std::string to enum
  inline enum SolutionType String2SolutionType(const std::string& name)
  {
    SolutionType type = sol_unknown;
    if (boost::iequals(name,"structure"))
      type = sol_structure;
    else if (boost::iequals(name,"contact"))
      type = sol_contact;
    else if (boost::iequals(name,"meshtying"))
      type = sol_meshtying;
    else if (boost::iequals(name,"cardiovascular0d"))
      type = sol_cardiovascular0d;
    else if (boost::iequals(name,"lag-pen-constraint"))
      type = sol_lag_pen_constraint;

    return type;
  };

  //! type of the optimization problem
  enum OptimizationProblemType
  {
    // unconstrained optimization problem
    opt_unconstrained,
    // pure equality constrained optimization problem
    opt_equality_constrained,
    // inequality or mixed optimization problem
    opt_inequality_constrained
  };
namespace LinSystem {
  //! supported LinearSystem types
  enum LinearSystemType {
    linear_system_structure,
    linear_system_structure_contact,
    linear_system_structure_cardiovascular0d,
    linear_system_structure_lag_pen_constraint,
    linear_system_undefined
  };

  /// Map quantity name to std::string
  static inline const std::string LinearSystemType2String(const enum LinearSystemType& type)
  {
    switch (type)
    {
    case linear_system_structure:
      return "linear_system_structure";
      break;
    case linear_system_structure_contact:
      return "linear_system_structure_contact";
      break;
    case linear_system_structure_cardiovascular0d:
      return "linear_system_structure_cardiovascular0d";
      break;
    case linear_system_structure_lag_pen_constraint:
      return "linear_system_structure_lag_pen_constraint";
      break;
    case linear_system_undefined:
      return "linear_system_undefined";
      break;
    default :
      return "unknown operator type";
      break;
    }
    return "";
  };
} // namespace LinSystem
namespace MeritFunction {
  //! order of the linearization term
  enum LinOrder
  {
    //! all orders
    linorder_all,
    //! first order linearization
    linorder_first,
    //! second order linearization terms
    linorder_second
  };

  //! type of the linearization term, i.e. with respect to which quantity the linearization was performed
  enum LinType
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
  enum MeritFctName
  {
    //! sum of squares merit function
    mrtfct_sum_of_squares,
    //! lagrangian merit function
    mrtfct_lagrangian
  };
} // namespace MeritFunction
namespace StatusTest {
  //! Supported quantity names for distinguished status tests
  enum QuantityType {
    quantity_unknown,
    // check structural quantities
    quantity_structure,
    // check (semi-smooth) contact quantities (normal/frictionless)
    quantity_contact_normal,
    // check (semi-smooth) contact quantities (frictionless)
    quantity_contact_friction,
    // check meshtying quantities
    quantity_meshtying,
    // check 0d cardiovascular quantities
    quantity_cardiovascular0d,
    // check Lagrange/penalty enforced constraint quantites
    quantity_lag_pen_constraint,
    // check semi-smooth contact
    quantity_plasticity,
    // check pressure dofs
    quantity_pressure,
    // check eas dofs
    quantity_eas
  };

  /// Map quantity name to std::string
  static inline const std::string QuantityType2String(const enum QuantityType type)
  {
    switch (type)
    {
    case quantity_structure:
      return "STRUCTURE";
      break;
    case quantity_contact_normal:
      return "CONTACT-NORMAL";
      break;
    case quantity_contact_friction:
      return "CONTACT-FRICTION";
      break;
    case quantity_meshtying:
      return "MESHTYING";
      break;
    case quantity_cardiovascular0d:
      return "CARDIOVASCULAR0D";
      break;
    case quantity_lag_pen_constraint:
      return "LAG-PEN-CONSTRAINT";
      break;
    case quantity_plasticity:
      return "PLASTICITY";
      break;
    case quantity_pressure:
      return "PRESSURE";
      break;
    case quantity_eas:
      return "EAS";
      break;
    case quantity_unknown:
    default :
      return "unknown quantity type";
      break;
    }
    return "";
  };

  inline QuantityType String2QuantityType(const std::string& name)
  {
    QuantityType type = quantity_unknown;

    // -- case insensitive comparison --
    if (boost::iequals(name,"structure"))
      type = quantity_structure;
    else if (boost::iequals(name,"contact-normal"))
      type = quantity_contact_normal;
    else if (boost::iequals(name,"contact-friction"))
      type = quantity_contact_friction;
    else if (boost::iequals(name,"meshtying"))
      type = quantity_meshtying;
    else if (boost::iequals(name,"cardiovascular0d"))
      type = quantity_cardiovascular0d;
    else if (boost::iequals(name,"lag-pen-constraint"))
      type = quantity_lag_pen_constraint;
    else if (boost::iequals(name,"plasticity"))
      type = quantity_plasticity;
    else if (boost::iequals(name,"pressure"))
      type = quantity_pressure;
    else if (boost::iequals(name,"eas"))
      type = quantity_eas;

    return type;
  };
} // namespace StatusTest
} // namespace NLN
} // namespace NOX



#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_ENUM_LISTS_H_ */
