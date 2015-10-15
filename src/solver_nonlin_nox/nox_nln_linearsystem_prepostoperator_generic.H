/*-----------------------------------------------------------*/
/*!
\file nox_nln_linearsystem_prepostoperator_generic.H

\maintainer Michael Hiermeier

\date Oct 7, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_LINEARSYSTEM_PREPOSTOPERATOR_GENERIC_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_LINEARSYSTEM_PREPOSTOPERATOR_GENERIC_H_

// forward declaration
class Epetra_Problem;
namespace LINALG {
  class SparseOperator;
} // namespace LINALG
namespace NOX {
namespace Abstract {
  class Vector;
} // namespace Abstract
namespace NLN {
  class LinearSystem;
namespace LinSystem {
namespace PrePostOp {
/*!
  \brief %NOX::NLN's pure virtual class to allow users to insert pre and post
  operations into NOX::NLN linear systems (before and after the
  NOX::NLN::LinearSystem::applyJacobianInverse() method).

  The user should implement their own concrete implementation of this class
  and register it as a Teuchos::RCP<NOX::NLN::LinSystem::PrePostOp::Generic>
  in the "Linear Solver" sub-sublist.

  To create and use a user defined pre/post operator for the linearsystem:

  <ol>

  <li> Create a pre/post operator that derives from
  NOX::NLN::LinSystem::PrePostOp::Generic. For example, the pre/post operator \c
  Foo might be defined as shown below.

  \code
  class Foo : public NOX::NLN::LinSystem::PrePostOp::Generic
  {
  // Insert class definition here
  }
  \endcode

  <li> Create the appropriate entries in the linear solver parameter list which belongs to the current direction method, as follows.

  \code
  Teuchos::RCP<Foo> foo = Teuchos::rcp(new Foo);
  const std::string& dir_str = paramsPtr->sublist("Direction").get<std::string>("Method");
  Teuchos::ParameterList& p_linsolver = paramsPtr->sublist("Direction").sublist(dir_str).
      sublist("Linear Solver").set<Teuchos::RCP<NOX::NLN::LinSystem::PrePostOp::Generic> >
      ("User Defined Pre/Post Operator",foo);
  \endcode

  <li> See also the nox_nln_solver_ptc implementation for a short example.

  </ol>

  \author Michael Hiermeier
 */
  class Generic
  {
  public:
    //! constructor (does nothing)
    Generic() {};

    //! Copy constructor (does nothing)
    Generic(const NOX::NLN::LinSystem::PrePostOp::Generic& source) {};

    //! destructor (does nothing)
    virtual ~Generic() {};

    /*! User defined method that will be executed at the start
        of a call to NOX::NLN::LinearSystem::applyJacobianInverse().

       \param rhs    : full access to the rhs vector
       \param jac    : full access to the jacobian
       \param linsys : read only access to the linear system object
     */
    virtual void runPreApplyJacobianInverse(
        NOX::Abstract::Vector& rhs,
        LINALG::SparseOperator& jac,
        const NOX::NLN::LinearSystem& linsys);

    /*! User defined method that will be executed at the end
        of a call to NOX::NLN::LinearSystem::applyJacobianInverse().

       \param rhs    : full access to the rhs vector
       \param jac    : full access to the jacobian
       \param linsys : read only access to the linear system object
     */
    virtual void runPostApplyJacobianInverse(
        NOX::Abstract::Vector& rhs,
        LINALG::SparseOperator& jac,
        const NOX::NLN::LinearSystem& linsys);

  };  // class Generic
} // namespace PrePostOp
} // namespace LinSystem
} // namespace NLN
} // namespace NOX

inline void NOX::NLN::LinSystem::PrePostOp::Generic::runPreApplyJacobianInverse(
    NOX::Abstract::Vector& rhs,
    LINALG::SparseOperator& jac,
    const NOX::NLN::LinearSystem& linsys)
{
  return;
}

inline void NOX::NLN::LinSystem::PrePostOp::Generic::runPostApplyJacobianInverse(
    NOX::Abstract::Vector& rhs,
    LINALG::SparseOperator& jac,
    const NOX::NLN::LinearSystem& linsys)
{
  return;
}

#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_LINEARSYSTEM_PREPOSTOPERATOR_GENERIC_H_ */
