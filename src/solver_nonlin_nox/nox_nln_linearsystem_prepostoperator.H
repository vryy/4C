/*-----------------------------------------------------------*/
/*! \file

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_LINEARSYSTEM_PREPOSTOPERATOR_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_LINEARSYSTEM_PREPOSTOPERATOR_H_

#include <Teuchos_RCP.hpp>
#include "nox_nln_abstract_prepostoperator.H"  // wrapped abstract class


// forward declarations
class Epetra_LinearProblem;
namespace LINALG
{
  class SparseOperator;
}  // namespace LINALG
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace NOX
{
  class Utils;
  namespace Abstract
  {
    class Vector;
  }  // namespace Abstract
  namespace NLN
  {
    class LinearSystem;
    namespace LinSystem
    {
      //! currently supported pre/post operators
      enum PrePostOpType
      {
        prepost_ptc,  //!< pseudo transient continuation pre/post operator
        prepost_dbc   //!< dirichlet boundary condition pre/post operator
      };

      /*!
        @brief Functor to process the pre/post operator object in the parameter list for the linear
        system.

        This is a wrapper class for a user derived  NOX::NLN::Abstract::PrePostOperator (ppo)
        object. All NOX::NLN linear systems use this class so we don't have to repeat all parsing
        code in each NOX::NLN linear system class. This class searches the "Linear Solver" parameter
        list passed into the constructor and if a ppo is found will wrap the object.

        For instructions on how to implement a PrePostOperator, see
        NOX::NLN::Abstract::PrePostOperator.

        \author Michael Hiermeier
      */
      class PrePostOperator
      {
       public:
        typedef std::map<enum PrePostOpType, Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator>> Map;

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
        virtual ~PrePostOperator(){};

        //! Resets the pre/post operator.
        virtual void reset(Teuchos::ParameterList& linearSolverSublist);

        /*! User defined method that will be executed at the start
            of a call to NOX::NLN::LinearSystem::applyJacobianInverse().

           \param rhs    : full access to the rhs vector
           \param jac    : full access to the jacobian
           \param linsys : read only access to the linear system object
         */
        virtual void runPreApplyJacobianInverse(NOX::Abstract::Vector& rhs,
            LINALG::SparseOperator& jac, const NOX::NLN::LinearSystem& linsys);

        /*! User defined method that will be executed at the end
            of a call to NOX::NLN::LinearSystem::applyJacobianInverse().

           \param result : full access to the result vector
           \param rhs    : full access to the rhs vector
           \param jac    : full access to the jacobian
           \param linsys : read only access to the linear system object
         */
        virtual void runPostApplyJacobianInverse(NOX::Abstract::Vector& result,
            NOX::Abstract::Vector& rhs, LINALG::SparseOperator& jac,
            const NOX::NLN::LinearSystem& linsys);

        /** User defined method that will be executed at the start of a call to
         * NOX::NLN::LinearSystem::applyJacobianInverse().
         *
         * \param jac    : full access to the jacobian operator
         * \param x      : read only access to the current solution point
         * \param linsys : read only access to the linear system object
         */
        virtual void runPreComputeJacobian(LINALG::SparseOperator& jac, const Epetra_Vector& x,
            const NOX::NLN::LinearSystem& linsys);

        /** User defined method that will be executed at the end of a call to
         * NOX::NLN::LinearSystem::applyJacobianInverse().
         *
         * \param jac    : full access to the jacobian operator
         * \param x      : read only access to the current solution point
         * \param linsys : read only access to the linear system object
         */
        virtual void runPostComputeJacobian(LINALG::SparseOperator& jac, const Epetra_Vector& x,
            const NOX::NLN::LinearSystem& linsys);

        /** User defined method that will be executed at the start of a call to
         * NOX::NLN::LinearSystem::computeFandJacobian().
         *
         * \param rhs    : full access to the right-hand-side vector
         * \param jac    : full access to the jacobian operator
         * \param x      : read only access to the current solution point
         * \param linsys : read only access to the linear system object
         */
        virtual void runPreComputeFandJacobian(Epetra_Vector& rhs, LINALG::SparseOperator& jac,
            const Epetra_Vector& x, const NOX::NLN::LinearSystem& linsys);

        /** User defined method that will be executed at the end of a call to
         * NOX::NLN::LinearSystem::computeFandJacobian().
         *
         * \param rhs    : full access to the right-hand-side vector
         * \param jac    : full access to the jacobian operator
         * \param x      : read only access to the current solution point
         * \param linsys : read only access to the linear system object
         */
        virtual void runPostComputeFandJacobian(Epetra_Vector& rhs, LINALG::SparseOperator& jac,
            const Epetra_Vector& x, const NOX::NLN::LinearSystem& linsys);

       protected:
        //! Flag that determines if a pre/post operator has been supplied by user.
        bool havePrePostOperator_;

        //! Points to user defined pre/post operator for the linear system.
        Teuchos::RCP<Map> prePostOperatorMapPtr_;
      };  // class PrePostOperator
      namespace PrePostOp
      {
        // non-member function
        /*! Returns the inherent pre/post operator std::map of the "Linear Solver" sublist.
         *  If the corresponding parameter called "User Defined Pre/Post Operator" is not yet
         *  defined, a empty std::map is generated and set into the parameter list first. */
        NOX::NLN::LinSystem::PrePostOperator::Map& GetMutableMap(
            Teuchos::ParameterList& p_linsolver);
      }  // namespace PrePostOp
    }    // namespace LinSystem
  }      // namespace NLN
}  // namespace NOX

inline void NOX::NLN::LinSystem::PrePostOperator::runPreApplyJacobianInverse(
    NOX::Abstract::Vector& rhs, LINALG::SparseOperator& jac, const NOX::NLN::LinearSystem& linsys)
{
  if (havePrePostOperator_)
  {
    Map::iterator it;
    for (it = prePostOperatorMapPtr_->begin(); it != prePostOperatorMapPtr_->end(); ++it)
      it->second->runPreApplyJacobianInverse(rhs, jac, linsys);
  }
}

inline void NOX::NLN::LinSystem::PrePostOperator::runPostApplyJacobianInverse(
    NOX::Abstract::Vector& result, NOX::Abstract::Vector& rhs, LINALG::SparseOperator& jac,
    const NOX::NLN::LinearSystem& linsys)
{
  if (havePrePostOperator_)
  {
    Map::iterator it;
    for (it = prePostOperatorMapPtr_->begin(); it != prePostOperatorMapPtr_->end(); ++it)
      it->second->runPostApplyJacobianInverse(result, rhs, jac, linsys);
  }
}

inline void NOX::NLN::LinSystem::PrePostOperator::runPreComputeJacobian(
    LINALG::SparseOperator& jac, const Epetra_Vector& x, const NOX::NLN::LinearSystem& linsys)
{
  if (havePrePostOperator_)
  {
    Map::iterator it;
    for (it = prePostOperatorMapPtr_->begin(); it != prePostOperatorMapPtr_->end(); ++it)
      it->second->runPreComputeJacobian(jac, x, linsys);
  }
}

inline void NOX::NLN::LinSystem::PrePostOperator::runPostComputeJacobian(
    LINALG::SparseOperator& jac, const Epetra_Vector& x, const NOX::NLN::LinearSystem& linsys)
{
  if (havePrePostOperator_)
  {
    Map::iterator it;
    for (it = prePostOperatorMapPtr_->begin(); it != prePostOperatorMapPtr_->end(); ++it)
      it->second->runPostComputeJacobian(jac, x, linsys);
  }
}

inline void NOX::NLN::LinSystem::PrePostOperator::runPreComputeFandJacobian(Epetra_Vector& rhs,
    LINALG::SparseOperator& jac, const Epetra_Vector& x, const NOX::NLN::LinearSystem& linsys)
{
  if (havePrePostOperator_)
  {
    Map::iterator it;
    for (it = prePostOperatorMapPtr_->begin(); it != prePostOperatorMapPtr_->end(); ++it)
      it->second->runPreComputeFandJacobian(rhs, jac, x, linsys);
  }
}

inline void NOX::NLN::LinSystem::PrePostOperator::runPostComputeFandJacobian(Epetra_Vector& rhs,
    LINALG::SparseOperator& jac, const Epetra_Vector& x, const NOX::NLN::LinearSystem& linsys)
{
  if (havePrePostOperator_)
  {
    Map::iterator it;
    for (it = prePostOperatorMapPtr_->begin(); it != prePostOperatorMapPtr_->end(); ++it)
      it->second->runPostComputeFandJacobian(rhs, jac, x, linsys);
  }
}

#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_LINEARSYSTEM_PREPOSTOPERATOR_H_ */
