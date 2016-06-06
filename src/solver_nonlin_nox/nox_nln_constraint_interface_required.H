/*-----------------------------------------------------------*/
/*!
\file nox_nln_constraint_interface_required.H

\brief Required interface for constrained problems.
       (necessary for the NOX::NLN::CONSTRAINT::Group
        and the evaluation of special constraint status
        tests)

\maintainer Michael Hiermeier

\date Jun 9, 2015

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
namespace LINALG {
  class SparseOperator;
} // namespace LINALG

class Epetra_Operator;
class Epetra_Vector;
class Epetra_Map;

namespace NOX {
namespace Epetra {
  class Vector;
} // namespace Epetra
namespace NLN {
namespace CONSTRAINT {
namespace Interface {
  class Required
  {
  public:
    //! Constructor
    Required() {};

    //! Destrcutor
    virtual ~Required() {};

    /*! @name Merit function support functions
     *  These functions are optional. They become only necessary, if you want to use the full functionality
     *  of the non-linear constraint solver framework (e.g. filter methods, etc.). */
    //! @{

    /*! \brief Get the objective model
     *
     *  This value can be calculated as a combination of the objective function,
     *  which we try to minimize and the subjected constraint equations.
     *  Typical examples are the Lagrangian function value and the augmented Lagrangian
     *  function value. */
    virtual double GetConstrObjectiveModelValue(const NOX::NLN::MeritFunction::MeritFctName& name) const
    {
      dserror("GetConstrObjectiveModelValue() is not implemented!");
      return -1.0;
    };

    //! Get the desired linearization terms of the objective model
    Teuchos::RCP<const std::vector<double> > GetLinearizedObjectiveModelTerms(
          const NOX::NLN::MeritFunction::MeritFctName& name) const
      { return GetLinearizedObjectiveModelTerms(name,NOX::NLN::MeritFunction::linorder_first,
          NOX::NLN::MeritFunction::lin_wrt_all_dofs); };
    Teuchos::RCP<const std::vector<double> > GetLinearizedObjectiveModelTerms(
        const NOX::NLN::MeritFunction::MeritFctName& name,
        const NOX::NLN::MeritFunction::LinOrder& order) const
    { return GetLinearizedObjectiveModelTerms(name,order,NOX::NLN::MeritFunction::lin_wrt_all_dofs); };
    virtual Teuchos::RCP<const std::vector<double> >
    GetLinearizedObjectiveModelTerms(
        const NOX::NLN::MeritFunction::MeritFctName& name,
        const NOX::NLN::MeritFunction::LinOrder& order,
        const NOX::NLN::MeritFunction::LinType& type) const
    {
      dserror("GetLinearizedObjectiveModelTerms() is not implemented!");
      return Teuchos::null;
    };

    //! Get value of the infeasibility model
    virtual double GetInfeasibilityModelValue() const
    {
      dserror("GetInfeasibilityModelValue() is not implemented!");
      return -1.0;
    };

    //! Get the desired linearization terms of the infeasibility model
    Teuchos::RCP<const std::vector<double> >
      GetLinearizedInfeasibilityModelTerms() const
    { return GetLinearizedInfeasibilityModelTerms(NOX::NLN::MeritFunction::linorder_first,
        NOX::NLN::MeritFunction::lin_wrt_all_dofs); };
    Teuchos::RCP<const std::vector<double> >
      GetLinearizedInfeasibilityModelTerms(
          const NOX::NLN::MeritFunction::LinOrder& order) const
    { return GetLinearizedInfeasibilityModelTerms(order,NOX::NLN::MeritFunction::lin_wrt_all_dofs); };
    virtual Teuchos::RCP<const std::vector<double> >
    GetLinearizedInfeasibilityModelTerms(
        const NOX::NLN::MeritFunction::LinOrder& order,
        const NOX::NLN::MeritFunction::LinType& type) const
    {
      dserror("GetLinearizedInfeasibilityModelTerms() is not implemented!");
      return Teuchos::null;
    };
    //! @}

    //! @name Status test support functions
    //! @{

    //! Returns the constraint right-hand-side norms
    double GetConstraintRHSNorms(
        const NOX::NLN::StatusTest::QuantityType& chQ) const
    { return GetConstraintRHSNorms(chQ,NOX::Abstract::Vector::TwoNorm,false); };
    double GetConstraintRHSNorms(
        const NOX::NLN::StatusTest::QuantityType& chQ,
        const NOX::Abstract::Vector::NormType& type) const
    { return GetConstraintRHSNorms(chQ,type,false); };
    virtual double GetConstraintRHSNorms(
        const NOX::NLN::StatusTest::QuantityType& chQ,
        const NOX::Abstract::Vector::NormType& type,
        const bool& isScaled) const = 0;

    //! Returns the Root Mean Square (abbr.: RMS) of the Lagrange multiplier updates
    double GetLagrangeMultiplierUpdateRMS(
        const Epetra_Vector& xNew, const Epetra_Vector& xOld,
        const double& aTol, const double& rTol,
        const NOX::NLN::StatusTest::QuantityType& checkQuantity) const
    { return GetLagrangeMultiplierUpdateRMS(xNew,xOld,aTol,rTol,checkQuantity,false); };
    virtual double GetLagrangeMultiplierUpdateRMS(
        const Epetra_Vector& xNew, const Epetra_Vector& xOld,
        const double& aTol, const double& rTol,
        const NOX::NLN::StatusTest::QuantityType& checkQuantity,
        const bool& disable_implicit_weighting) const = 0;

    //! Returns the increment norm of the largange multiplier DoFs
    virtual double GetLagrangeMultiplierUpdateNorms(
        const Epetra_Vector& xNew, const Epetra_Vector& xOld,
        const NOX::NLN::StatusTest::QuantityType& checkQuantity,
        const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
        const bool& isScaled = false) const = 0;

    //! Returns the previous solution norm of the largange multiplier DoFs
    virtual double GetPreviousLagrangeMultiplierNorms(
        const Epetra_Vector& xOld,
        const NOX::NLN::StatusTest::QuantityType& checkQuantity,
        const NOX::Abstract::Vector::NormType& type = NOX::Abstract::Vector::TwoNorm,
        const bool& isScaled = false) const = 0;

    /*! @name Handle active set changes.
     *  This is optional and only relevant for inequality constraint problems. */
    //! @{
    virtual enum NOX::StatusTest::StatusType GetActiveSetInfo(
        const enum NOX::NLN::StatusTest::QuantityType& qt,
        int& activeset_size) const
    {
      activeset_size = -1;
      return NOX::StatusTest::Unevaluated;
    }

    virtual Teuchos::RCP<const Epetra_Map> GetCurrentActiveSetMap(
        const enum NOX::NLN::StatusTest::QuantityType& qt) const
    { return Teuchos::null; };

    virtual Teuchos::RCP<const Epetra_Map> GetOldActiveSetMap(
        const enum NOX::NLN::StatusTest::QuantityType& qt) const
    { return Teuchos::null; };
    //! @}
    //! @}
  };
} // end namespace Interface
  // typedef
  typedef std::map<::NOX::NLN::SolutionType,Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> > ReqInterfaceMap;
} // end namespace Constraint
} // end namespace NLN
} // end namespace NOX


#endif /* NOX_NLN_CONSTRAINT_INTERFACE_REQUIRED_H_ */
