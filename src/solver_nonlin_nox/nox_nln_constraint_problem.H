/*-----------------------------------------------------------*/
/*!
\file nox_nln_constraint_problem.H

\maintainer Michael Hiermeier

\date Jun 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_CONSTRAINT_PROBLEM_H_
#define NOX_NLN_CONSTRAINT_PROBLEM_H_

#include "nox_nln_problem.H"    // base class
#include "nox_nln_constraint_interface_required.H"

// forward declaration
class Epetra_Vector;

namespace NOX {
namespace Epetra {
  class Vector;
} // namespace Epetra

namespace NLN {
namespace CONSTRAINT {
  class NoxProblem : public NOX::NLN::NoxProblem
  {
  public:
    //! constructor
    NoxProblem(Teuchos::RCP<NOX::NLN::GlobalData>& nlnGlobalData,
        Teuchos::RCP<Epetra_Vector>& x,
        Teuchos::RCP<LINALG::SparseOperator> jac);

    //! destructor
    virtual ~NoxProblem() {};

    // derived
    virtual void Initialize(Teuchos::RCP<Epetra_Vector>& x1, Teuchos::RCP<LINALG::SparseOperator> A11);

    // create a nox group for constraint problems
    Teuchos::RCP<NOX::Abstract::Group> CreateNoxGroup(
        const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys) const;

    //! Returns the constraint interface pointer
    virtual const Teuchos::RCP<NOX::NLN::CONSTRAINT::Interface::Required> GetConstrPtr() const;
  };
} // namespace CONSTRAINT
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_CONSTRAINT_PROBLEM_H_ */
