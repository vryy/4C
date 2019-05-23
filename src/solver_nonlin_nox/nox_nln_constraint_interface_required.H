/*-----------------------------------------------------------*/
/*!

\brief Required interface for constrained problems.
       (necessary for the NOX::NLN::CONSTRAINT::Group
        and the evaluation of special constraint status
        tests)

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_CONSTRAINT_INTERFACE_REQUIRED_H_
#define NOX_NLN_CONSTRAINT_INTERFACE_REQUIRED_H_

#include "nox_nln_enum_lists.H"
#include "../drt_lib/drt_dserror.H"

#include <NOX_StatusTest_Generic.H>
#include <NOX_Abstract_Vector.H>


// forward declaration...
namespace LINALG
{
  class SparseOperator;
  class SerialDenseVector;
}  // namespace LINALG

class Epetra_Operator;
class Epetra_Vector;
class Epetra_Map;

namespace NOX
{
  namespace Epetra
  {
    class Vector;
  }  // namespace Epetra
  namespace NLN
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
          virtual ~Required(){};

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
          virtual double GetModelValue(NOX::NLN::MeritFunction::MeritFctName name) const
          {
            dserror("GetObjectiveModelValue() is not implemented!");
            exit(EXIT_FAILURE);
          };

          //! Get the desired linearization terms of the objective model
          virtual double GetLinearizedModelTerms(const Epetra_Vector& dir,
              const enum NOX::NLN::MeritFunction::MeritFctName name,
              const enum NOX::NLN::MeritFunction::LinOrder order,
              const enum NOX::NLN::MeritFunction::LinType type) const
          {
            dserror("GetLinearizedObjectiveModelTerms() is not implemented!");
            exit(EXIT_FAILURE);
          };

          //! @}

          //! @name Status test support functions
          //! @{

          //! Returns the constraint right-hand-side norms
          double GetConstraintRHSNorms(
              const Epetra_Vector& F, NOX::NLN::StatusTest::QuantityType chQ) const
          {
            return GetConstraintRHSNorms(F, chQ, NOX::Abstract::Vector::TwoNorm, false);
          };
          double GetConstraintRHSNorms(const Epetra_Vector& F,
              NOX::NLN::StatusTest::QuantityType chQ, NOX::Abstract::Vector::NormType type) const
          {
            return GetConstraintRHSNorms(F, chQ, type, false);
          };
          virtual double GetConstraintRHSNorms(const Epetra_Vector& F,
              NOX::NLN::StatusTest::QuantityType chQ, NOX::Abstract::Vector::NormType type,
              bool isScaled) const = 0;

          //! Returns the Root Mean Square (abbr.: RMS) of the Lagrange multiplier updates
          double GetLagrangeMultiplierUpdateRMS(const Epetra_Vector& xNew,
              const Epetra_Vector& xOld, double aTol, double rTol,
              NOX::NLN::StatusTest::QuantityType checkQuantity) const
          {
            return GetLagrangeMultiplierUpdateRMS(xNew, xOld, aTol, rTol, checkQuantity, false);
          };
          virtual double GetLagrangeMultiplierUpdateRMS(const Epetra_Vector& xNew,
              const Epetra_Vector& xOld, double aTol, double rTol,
              NOX::NLN::StatusTest::QuantityType checkQuantity,
              bool disable_implicit_weighting) const = 0;

          //! Returns the increment norm of the largange multiplier DoFs
          virtual double GetLagrangeMultiplierUpdateNorms(const Epetra_Vector& xNew,
              const Epetra_Vector& xOld, NOX::NLN::StatusTest::QuantityType checkQuantity,
              NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TwoNorm,
              bool isScaled = false) const = 0;

          //! Returns the previous solution norm of the largange multiplier DoFs
          virtual double GetPreviousLagrangeMultiplierNorms(const Epetra_Vector& xOld,
              NOX::NLN::StatusTest::QuantityType checkQuantity,
              NOX::Abstract::Vector::NormType type = NOX::Abstract::Vector::TwoNorm,
              bool isScaled = false) const = 0;

          /*! @name Handle active set changes.
           *  This is optional and only relevant for inequality constraint problems. */
          //! @{
          virtual enum NOX::StatusTest::StatusType GetActiveSetInfo(
              enum NOX::NLN::StatusTest::QuantityType qt, int& activeset_size) const
          {
            activeset_size = -1;
            return NOX::StatusTest::Unevaluated;
          }

          virtual Teuchos::RCP<const Epetra_Map> GetCurrentActiveSetMap(
              enum NOX::NLN::StatusTest::QuantityType qt) const
          {
            return Teuchos::null;
          };

          virtual Teuchos::RCP<const Epetra_Map> GetOldActiveSetMap(
              enum NOX::NLN::StatusTest::QuantityType qt) const
          {
            return Teuchos::null;
          };
          //! @}
          //! @}
        };
      }  // end namespace Interface
      // typedef
      typedef std::map<::NOX::NLN::SolutionType,
          Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required>>
          ReqInterfaceMap;
    }  // namespace CONSTRAINT
  }    // end namespace NLN
}  // end namespace NOX


#endif /* NOX_NLN_CONSTRAINT_INTERFACE_REQUIRED_H_ */
