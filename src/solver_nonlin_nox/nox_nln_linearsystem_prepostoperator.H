/*-----------------------------------------------------------*/
/*!
\file nox_nln_linearsystem_prepostoperator.H

\maintainer Michael Hiermeier

\date Oct 7, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_LINEARSYSTEM_PREPOSTOPERATOR_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_LINEARSYSTEM_PREPOSTOPERATOR_H_

#include <Teuchos_RCP.hpp>
#include "nox_nln_linearsystem_prepostoperator_generic.H"


// forward declarations
class Epetra_LinearProblem;
namespace LINALG {
  class SparseOperator;
} // namespace LINALG
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
namespace NOX {
  class Utils;
namespace Abstract {
  class Vector;
} // namespace Abstract
namespace NLN {
  class LinearSystem;
namespace LinSystem {
/*!
  @brief Functor to process the pre/post operator object in the parameter list for the linear system.

  This is a wrapper class for a user derived  NOX::NLN::LinSystem::PrePostOp::Generic (ppo) object.
  All NOX::NLN linear systems use this class so we don't have to repeat all parsing code in each
  NOX::NLN linear system class. This class searches the "Linear Solver" parameter list passed into
  the constructor and if a ppo is found will wrap the object.

  For instructions on how to implement a PrePostOperator, see
  NOX::NLN::LinSystem::PrePostOp::Generic.

  \author Michael Hiermeier
*/
  class PrePostOperator
  {
  private:
    //! Disallow default constructor.
    PrePostOperator();

    //! Disallow copy constructor.
    PrePostOperator(const PrePostOperator& p);

    //! Disallow assignment operator.
    PrePostOperator& operator=(const PrePostOperator& ppo);

  public:

    //! Allowed constructor.
    PrePostOperator(Teuchos::ParameterList& linearSolverSubList);

    //! Destructor.
    virtual ~PrePostOperator();

    //! Resets the pre/post operator.
    virtual void reset(Teuchos::ParameterList& linearSolverSublist);

    /*! User defined method that will be executed at the start
        of a call to NOX::NLN::LinearSystem::applyJacobianInverse().

       \param rhs    : full access to the rhs vector
       \param jac    : full access to the jacobian
       \param linsys : read only access to the linear system object
     */
    virtual void runPreApplyJacobianInverse(NOX::Abstract::Vector& rhs,
        LINALG::SparseOperator& jac,
        const NOX::NLN::LinearSystem& linsys);

    /*! User defined method that will be executed at the end
        of a call to NOX::NLN::LinearSystem::applyJacobianInverse().

       \param rhs    : full access to the rhs vector
       \param jac    : full access to the jacobian
       \param linsys : read only access to the linear system object
     */
    virtual void runPostApplyJacobianInverse(NOX::Abstract::Vector& rhs,
        LINALG::SparseOperator& jac,
        const NOX::NLN::LinearSystem& linsys);

  protected:
    //! Flag that determines if a pre/post operator has been supplied by user.
    bool havePrePostOperator_;

    //! Points to user defined pre/post operator for the linear system.
    Teuchos::RCP<NOX::NLN::LinSystem::PrePostOp::Generic> prePostOperatorPtr_;
  };  // class PrePostOperator
} // namespace LinSystem
} // namespace NLN
} // namespace NOX

inline void NOX::NLN::LinSystem::PrePostOperator::runPreApplyJacobianInverse(
    NOX::Abstract::Vector& rhs,
    LINALG::SparseOperator& jac,
    const NOX::NLN::LinearSystem& linsys)
{
  if (havePrePostOperator_)
    prePostOperatorPtr_->runPreApplyJacobianInverse(rhs,jac,linsys);
}

inline void NOX::NLN::LinSystem::PrePostOperator::runPostApplyJacobianInverse(
    NOX::Abstract::Vector& rhs,
    LINALG::SparseOperator& jac,
    const NOX::NLN::LinearSystem& linsys)
{
  if (havePrePostOperator_)
    prePostOperatorPtr_->runPostApplyJacobianInverse(rhs,jac,linsys);
}

#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_LINEARSYSTEM_PREPOSTOPERATOR_H_ */
