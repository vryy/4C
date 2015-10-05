/*-----------------------------------------------------------*/
/*!
\file nox_nln_constraint_interface_required.H

\maintainer Michael Hiermeier

\date Jun 9, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_CONSTRAINT_INTERFACE_REQUIRED_H_
#define NOX_NLN_CONSTRAINT_INTERFACE_REQUIRED_H_

#include "../drt_lib/drt_dserror.H"

#include "nox_nln_globaldata.H"
#include "nox_nln_statustest_factory.H"

#include <NOX_StatusTest_Generic.H>
#include <NOX_Abstract_Vector.H>


// forward declaration...
namespace LINALG {
  class SparseOperator;
}

class Epetra_Operator;

namespace NOX {
namespace Epetra {
  class Vector;
}

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

  //! @name Merit function support functions

  //! Get the objective model
  /*! This value can be calculated as a combination of the objective function,
      which we try to minimize and the subjected constraint equations.
      Typical examples are the Lagrangian function value and the augmented Lagrangian
      function value. */
  virtual const double GetConstrObjectiveModelValue(const NOX::NLN::GlobalData::MeritFctName& name) const = 0;

  //! Get the desired linearization terms of the objective model
  virtual Teuchos::RCP<const std::vector<double> >
  GetLinearizedObjectiveModelTerms(const NOX::NLN::GlobalData::MeritFctName& name,
      const NOX::NLN::GlobalData::LinOrder& order=NOX::NLN::GlobalData::linorder_first,
      const NOX::NLN::GlobalData::LinType& type=NOX::NLN::GlobalData::lin_wrt_all_dofs) const = 0;

  //! Get value of the infeasibility model
  virtual const double GetInfeasibilityModelValue() const = 0;

  //! Get the desired linearization terms of the infeasibility model
  virtual Teuchos::RCP<const std::vector<double> >
  GetLinearizedInfeasibilityModelTerms(
      const NOX::NLN::GlobalData::LinOrder& order=NOX::NLN::GlobalData::linorder_first,
      const NOX::NLN::GlobalData::LinType& type=NOX::NLN::GlobalData::lin_wrt_all_dofs) const = 0;

  //! @}

  //! @name Basic functions

  //! Get the systemtype of the underlying linear system: condensed or saddlepoint
  virtual const bool IsCondensed() const = 0;

  //! @}

  //! @name Status test support functions

  //! Returns the constraint right-hand-side norms
  virtual double GetConstraintRHSNorms(
      const NOX::NLN::StatusTest::QuantityType& chQ,
      const NOX::Abstract::Vector::NormType& type=NOX::Abstract::Vector::TwoNorm,
      const bool& isScaled = false) const = 0;

  //! Returns the Root Mean Square (abbr.: RMS) of the Lagrange multiplier updates
  virtual double GetLagrangeMultiplierUpdateRMS(
      const double& aTol, const double& rTol,
      const NOX::NLN::StatusTest::QuantityType& checkQuantity,
      const bool& disable_implicit_weighting=false) const = 0;

  //! Returns the active set status.
  //! This is optional and only relevant for inequality constraint problems.
  virtual enum NOX::StatusTest::StatusType GetActiveSetStatus(int& activesetsize,
      int& cyclesize) const
  {
    return NOX::StatusTest::Unevaluated;
  }

  //! @}

  //! Sets the Lagrange multipliers which are saved internally by the
  //! interface classes.
  virtual void SetLagrangeMultiplier(
      const NOX::Epetra::Vector& xOld,
      const NOX::Epetra::Vector& d,
      const double& step)
  { /*will be deleted in a up-coming commit*/ };
};
} // end namespace Interface
} // end namespace Constraint
} // end namespace NLNSOL
} // end namespace NOX


#endif /* NOX_NLN_CONSTRAINT_INTERFACE_REQUIRED_H_ */
