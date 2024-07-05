/*-----------------------------------------------------------*/
/*! \file

\brief Required interface for constrained problems.
       (necessary for the NOX::Nln::CONSTRAINT::Group
        and the evaluation of special constraint status
        tests)

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_SOLVER_NONLIN_NOX_CONSTRAINT_INTERFACE_REQUIRED_HPP
#define FOUR_C_SOLVER_NONLIN_NOX_CONSTRAINT_INTERFACE_REQUIRED_HPP

#include "4C_config.hpp"

#include "4C_solver_nonlin_nox_enum_lists.hpp"
#include "4C_solver_nonlin_nox_forward_decl.hpp"
#include "4C_utils_exceptions.hpp"

#include <NOX_Abstract_Vector.H>
#include <NOX_StatusTest_Generic.H>

FOUR_C_NAMESPACE_OPEN

// forward declaration...
namespace Core::LinAlg
{
  class SparseOperator;
  class SerialDenseVector;
}  // namespace Core::LinAlg

namespace NOX
{
  namespace Nln
  {
    namespace CONSTRAINT
    {
      namespace Interface
      {
        class Required
        {
         public:
          //! Constructor
          Required(){};

          //! Destrcutor
          virtual ~Required() = default;

          /*! @name Merit function support functions
           *  These functions are optional. They become only necessary, if you want to use the full
           * functionality of the non-linear constraint solver framework (e.g. filter methods,
           * etc.). */
          //! @{

          /*! \brief Get the objective model
           *
           *  This value can be calculated as a combination of the objective function,
           *  which we try to minimize and the subjected constraint equations.
           *  Typical examples are the Lagrangian function value and the augmented Lagrangian
           *  function value. */
          virtual double get_model_value(NOX::Nln::MeritFunction::MeritFctName name) const
          {
            FOUR_C_THROW("get_model_value() is not implemented!");
          };

          //! Get the desired linearization terms of the objective model
          virtual double get_linearized_model_terms(const Epetra_Vector& dir,
              const enum NOX::Nln::MeritFunction::MeritFctName name,
              const enum NOX::Nln::MeritFunction::LinOrder order,
              const enum NOX::Nln::MeritFunction::LinType type) const
          {
            FOUR_C_THROW("get_linearized_model_terms() is not implemented!");
          };

          //! @}

          //! @name Status test support functions
          //! @{

          //! Returns the constraint right-hand-side norms
          double get_constraint_rhs_norms(
              const Epetra_Vector& F, NOX::Nln::StatusTest::QuantityType chQ) const
          {
            return get_constraint_rhs_norms(F, chQ, ::NOX::Abstract::Vector::TwoNorm, false);
          };
          double get_constraint_rhs_norms(const Epetra_Vector& F,
              NOX::Nln::StatusTest::QuantityType chQ, ::NOX::Abstract::Vector::NormType type) const
          {
            return get_constraint_rhs_norms(F, chQ, type, false);
          };
          virtual double get_constraint_rhs_norms(const Epetra_Vector& F,
              NOX::Nln::StatusTest::QuantityType chQ, ::NOX::Abstract::Vector::NormType type,
              bool isScaled) const = 0;

          //! Returns the Root Mean Square (abbr.: RMS) of the Lagrange multiplier updates
          double get_lagrange_multiplier_update_rms(const Epetra_Vector& xNew,
              const Epetra_Vector& xOld, double aTol, double rTol,
              NOX::Nln::StatusTest::QuantityType checkQuantity) const
          {
            return get_lagrange_multiplier_update_rms(xNew, xOld, aTol, rTol, checkQuantity, false);
          };
          virtual double get_lagrange_multiplier_update_rms(const Epetra_Vector& xNew,
              const Epetra_Vector& xOld, double aTol, double rTol,
              NOX::Nln::StatusTest::QuantityType checkQuantity,
              bool disable_implicit_weighting) const = 0;

          //! Returns the increment norm of the largange multiplier DoFs
          virtual double get_lagrange_multiplier_update_norms(const Epetra_Vector& xNew,
              const Epetra_Vector& xOld, NOX::Nln::StatusTest::QuantityType checkQuantity,
              ::NOX::Abstract::Vector::NormType type = ::NOX::Abstract::Vector::TwoNorm,
              bool isScaled = false) const = 0;

          //! Returns the previous solution norm of the largange multiplier DoFs
          virtual double get_previous_lagrange_multiplier_norms(const Epetra_Vector& xOld,
              NOX::Nln::StatusTest::QuantityType checkQuantity,
              ::NOX::Abstract::Vector::NormType type = ::NOX::Abstract::Vector::TwoNorm,
              bool isScaled = false) const = 0;

          /*! @name Handle active set changes.
           *  This is optional and only relevant for inequality constraint problems. */
          //! @{
          virtual enum ::NOX::StatusTest::StatusType get_active_set_info(
              enum NOX::Nln::StatusTest::QuantityType qt, int& activeset_size) const
          {
            activeset_size = -1;
            return ::NOX::StatusTest::Unevaluated;
          }

          virtual Teuchos::RCP<const Epetra_Map> get_current_active_set_map(
              enum NOX::Nln::StatusTest::QuantityType qt) const
          {
            return Teuchos::null;
          };

          virtual Teuchos::RCP<const Epetra_Map> get_old_active_set_map(
              enum NOX::Nln::StatusTest::QuantityType qt) const
          {
            return Teuchos::null;
          };
          //! @}
          //! @}
        };
      }  // end namespace Interface
      // typedef
      typedef std::map<NOX::Nln::SolutionType,
          Teuchos::RCP<NOX::Nln::CONSTRAINT::Interface::Required>>
          ReqInterfaceMap;
    }  // namespace CONSTRAINT
  }    // end namespace Nln
}  // end namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
