/*-----------------------------------------------------------*/
/*!
\file nox_nln_abstract_prepostoperator.H

\brief %NOX::NLN's pure virtual class to allow users to insert pre and post
       operations into different NOX::NLN classes.

\level 3

\maintainer Michael Hiermeier

\date Oct 19, 2015



*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_ABSTRACT_PREPOSTOPERATOR_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_ABSTRACT_PREPOSTOPERATOR_H_


#include <NOX_Abstract_PrePostOperator.H>     // base class

// forward declaration
class Epetra_Vector;
namespace LINALG {
  class SparseOperator;
} // namespace LINALG
namespace NOX {
namespace Abstract {
  class Vector;
} // namespace Abstract
namespace LineSearch {
  class Generic;
} // namespace LineSearch
namespace NLN {
  class Group;
  class LinearSystem;
namespace Abstract {
/*!
  \brief %NOX::NLN's pure virtual class to allow users to insert pre and post
  operations into different NOX::NLN classes.

  The user should implement their own concrete implementation of this class
  and register it as a Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator>
  in the corresponding sublist.

  For example: To create and use a user defined pre/post operator for the linear system you have to
  insert the object in the "Linear Solver" sub-sublist:

  <ol>

  <li> Create a pre/post operator that derives from
  NOX::NLN::Abstract::PrePostOperator. For example, the pre/post operator \c
  Foo might be defined as shown below.

  \code
  class Foo : public NOX::NLN::Abstract::PrePostOperator
  {
  // Insert class definition here
  }
  \endcode

  <li> Create the appropriate entries in the linear solver parameter list which belongs to
  the current direction method, as follows.

  \code
  Teuchos::RCP<Foo> foo = Teuchos::rcp(new Foo);
  const std::string& dir_str = paramsPtr->sublist("Direction").get<std::string>("Method");
  Teuchos::ParameterList& p_linsolver = paramsPtr->sublist("Direction").sublist(dir_str).
      sublist("Linear Solver").set<Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator> >
      ("User Defined Pre/Post Operator",foo);
  \endcode

  <li> See also the nox_nln_solver_ptc implementation for a short example.

  </ol>

  \author Michael Hiermeier
 */
  class PrePostOperator : public NOX::Abstract::PrePostOperator
  {
  public:
    //! constructor (does nothing)
    PrePostOperator() {};

    //! Copy constructor (does nothing)
    PrePostOperator(const NOX::NLN::Abstract::PrePostOperator& source) {};

    //! destructor (does nothing)
    virtual ~PrePostOperator() {};

    /** @name Solver Pre/Post Operator
     *  Non-linear solver pre/post-operator functions. See the  NOX::Solver::PrePostOperator
     *  class and its derived classes for more information. The virtual functions can be found in
     *  the base class.
     */
    ///@{
    /** User defined method that will be executed at the start of a call to NOX::Solver::Generic::iterate().
    virtual void runPreIterate(const NOX::Solver::Generic& solver); */

    /** User defined method that will be executed at the end of a call to NOX::Solver::Generic::iterate().
    virtual void runPostIterate(const NOX::Solver::Generic& solver); */

    /** User defined method that will be executed at the start of a call to NOX::Solver::Generic::solve().
    virtual void runPreSolve(const NOX::Solver::Generic& solver); */

    /** User defined method that will be executed at the end of a call to NOX::Solver::Generic::solve().
    virtual void runPostSolve(const NOX::Solver::Generic& solver); */
    ///@}

    /** @name NLN::LinearSystem Pre/Post Operator
     *  This pre/post operator is used in the NOX::NLN::LinearSystem class and its derived classes.
     */
    ///@{
    /** User defined method that will be executed at the start
     *  of a call to NOX::NLN::LinearSystem::applyJacobianInverse().
     *
     * \param rhs    : full access to the rhs vector
     * \param jac    : full access to the jacobian
     * \param linsys : read only access to the linear system object
     */
    virtual void runPreApplyJacobianInverse(
        NOX::Abstract::Vector& rhs,
        LINALG::SparseOperator& jac,
        const NOX::NLN::LinearSystem& linsys)
    { return; };

    /** User defined method that will be executed at the end
     *  of a call to NOX::NLN::LinearSystem::applyJacobianInverse().
     *
     * \param result : full access to the result vector
     * \param rhs    : full access to the rhs vector
     * \param jac    : full access to the jacobian
     * \param linsys : read only access to the linear system object
     */
    virtual void runPostApplyJacobianInverse(
        NOX::Abstract::Vector& result,
        NOX::Abstract::Vector& rhs,
        LINALG::SparseOperator& jac,
        const NOX::NLN::LinearSystem& linsys)
    { return; };

    /** User defined method that will be executed at the start of a call to NOX::NLN::LinearSystem::applyJacobianInverse().
     *
     * \param jac    : full access to the jacobian operator
     * \param x      : read only access to the current solution point
     * \param linsys : read only access to the linear system object
     */
    virtual void runPreComputeJacobian(
        LINALG::SparseOperator& jac,
        const Epetra_Vector& x,
        const NOX::NLN::LinearSystem& linsys)
    { return; };

    /** User defined method that will be executed at the end of a call to NOX::NLN::LinearSystem::applyJacobianInverse().
     *
     * \param jac    : full access to the jacobian operator
     * \param x      : read only access to the current solution point
     * \param linsys : read only access to the linear system object
     */
    virtual void runPostComputeJacobian(
        LINALG::SparseOperator& jac,
        const Epetra_Vector& x,
        const NOX::NLN::LinearSystem& linsys)
    { return; };

    /** User defined method that will be executed at the start of a call to NOX::NLN::LinearSystem::computeFandJacobian().
     *
     * \param rhs    : full access to the right-hand-side vector
     * \param jac    : full access to the jacobian operator
     * \param x      : read only access to the current solution point
     * \param linsys : read only access to the linear system object
     */
    virtual void runPreComputeFandJacobian(
        Epetra_Vector& rhs,
        LINALG::SparseOperator& jac,
        const Epetra_Vector& x,
        const NOX::NLN::LinearSystem& linsys)
    { return; };

    /** User defined method that will be executed at the end of a call to NOX::NLN::LinearSystem::computeFandJacobian().
     *
     * \param rhs    : full access to the right-hand-side vector
     * \param jac    : full access to the jacobian operator
     * \param x      : read only access to the current solution point
     * \param linsys : read only access to the linear system object
     */
    virtual void runPostComputeFandJacobian(
        Epetra_Vector& rhs,
        LINALG::SparseOperator& jac,
        const Epetra_Vector& x,
        const NOX::NLN::LinearSystem& linsys)
    { return; };

    ///@}

    /** @name NLN::Group Pre/Post Operator
     *  This pre/post operator is used in the NOX::NLN::Group class and its derived classes.
     */
    ///@{
    /** User defined method that will be executed at the start of a call to NOX::NLN::Group::computeF().
     *
     * \param F        : full access to the right hand side vector of the NOX::NLN::Group.
     * \param grp      : read only access to the NOX::NLN::Group object.
     */
    virtual void runPreComputeF(Epetra_Vector& F,
        const NOX::NLN::Group& grp)
    { return; };

    /** User defined method that will be executed at the end of a call to NOX::NLN::Group::computeF().
     *
     * \param F        : full access to the right hand side vector of the NOX::NLN::Group.
     * \param grp      : read only access to the NOX::NLN::Group object.
     */
    virtual void runPostComputeF( Epetra_Vector& F,
        const NOX::NLN::Group& grp)
    { return; };

    /** User defined method that will be executed at the start of a call to NOX::NLN::Group::computeX().
     *
     * \param input_grp: read only access to the input group (holds the old X).
     * \param dir      : read only access to the direction vector (step length equal 1.0).
     * \param step     : read only access to the current step length (line search).
     * \param curr_grp : read only access to the called/current group (will hold the new X).
     */
    virtual void runPreComputeX(
        const NOX::NLN::Group& input_grp,
        const Epetra_Vector& dir,
        const double& step,
        const NOX::NLN::Group& curr_grp)
    { return; };

    /** User defined method that will be executed at the end of a call to NOX::NLN::Group::computeX().
     *
     * \param input_grp: read only access to the input group (holds the old X).
     * \param dir      : read only access to the direction vector (step length equal 1.0).
     * \param step     : read only access to the current step length (line search).
     * \param curr_grp : read only access to the called/current group (holds the new X).
     */
    virtual void runPostComputeX(
        const NOX::NLN::Group& input_grp,
        const Epetra_Vector& dir,
        const double& step,
        const NOX::NLN::Group& curr_grp)
    { return; };

    /*! User defined method that will be executed at the beginning
     *  of a call to NOX::NLN::Group::applyJacobianInverse().
     *
     *  \param rhs    : read-only access to the rhs vector
     *  \param result : full access to the result vector
     *  \param xold   : read-only access to the jacobian
     *  \param grp    : read only access to the group object
     */
    virtual void runPreApplyJacobianInverse(
        const NOX::Abstract::Vector& rhs,
        NOX::Abstract::Vector& result,
        const NOX::Abstract::Vector& xold,
        const NOX::NLN::Group& grp )
    { return; };

    /*! User defined method that will be executed at the end
     *  of a call to NOX::NLN::Group::applyJacobianInverse().
     *
     *  \param rhs    : read-only access to the rhs vector
     *  \param result : full access to the result vector
     *  \param xold   : read-only access to the jacobian
     *  \param grp    : read only access to the group object
     */
    virtual void runPostApplyJacobianInverse(
        const NOX::Abstract::Vector& rhs,
        NOX::Abstract::Vector& result,
        const NOX::Abstract::Vector& xold,
        const NOX::NLN::Group& grp )
    {return; };

    ///@}

    /** @name NLN::LineSearch Pre/Post Operator
     *  This pre/post operator is used in the NOX::NLN::LineSearch classes.
     */
    ///@{
    /** User defined method that will be executed before the step is modified in
     *  the line search routine.
     *
     * \param solver     : Access to the underlying solver object.
     * \param linesearch : Access to the line search object. */
    virtual void runPreModifyStepLength(
        const NOX::Solver::Generic& solver,
        const NOX::LineSearch::Generic& linesearch )
    { return; }

    ///@}

  };  // class PrePostOperator
} // namespace Abstract
} // namespace NLN
} // namespace NOX


#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_ABSTRACT_PREPOSTOPERATOR_H_ */
