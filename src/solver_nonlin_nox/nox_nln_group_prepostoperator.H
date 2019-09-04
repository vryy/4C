/*-----------------------------------------------------------*/
/*! \file

\brief wrapper class for a user derived NOX PrePostOperator

\maintainer Anh-Tu Vuong


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_GROUP_PREPOSTOPERATOR_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_GROUP_PREPOSTOPERATOR_H_

#include <Teuchos_RCP.hpp>
#include "nox_nln_abstract_prepostoperator.H"


// forward declarations
class Epetra_Vector;
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace NOX
{
  namespace NLN
  {
    namespace GROUP
    {
      //! Currently supported group pre/post operators.
      enum PrePostOpType
      {
        prepost_ptc,           //!< pseudo transient continuation pre/post operator
        prepost_tangdis,       //!< tangdis (predictor) pre/post operator
        prepost_impl_generic,  //!< implicit generic (time integrator) pre/post operator
        prepost_rotvecupdate   //!< update of non-additive rotation vector DoFs
      };

      /*!
        @brief Functor to process the pre/post operator object in the parameter list for the
        NOX::NLN::Group objects.

        This is a wrapper class for a user derived  NOX::NLN::Abstract::PrePostOperator (ppo)
        object. All NOX::NLN groups use this class so we don't have to repeat all parsing code in
        each NOX::NLN group class. This class searches the "Group Options" parameter list passed
        into the constructor and if a ppo is found will wrap the object.

        For instructions on how to implement a PrePostOperator, see
        NOX::NLN::Abstract::PrePostOperator or one of the currently supported implementations (enum
        list).

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
        PrePostOperator(Teuchos::ParameterList& groupOptionsSublist);

        //! Destructor.
        virtual ~PrePostOperator(){};

        //! Resets the pre/post operator.
        virtual void reset(Teuchos::ParameterList& groupOptionsSublist);

        /** User defined method that will be executed at the start of a call to
         * NOX::NLN::Group::computeF().
         *
         * \param F        : full access to the right hand side vector of the NOX::NLN::Group.
         * \param grp      : read only access to the NOX::NLN::Group object.
         */
        virtual void runPreComputeF(Epetra_Vector& F, const NOX::NLN::Group& grp);

        /** User defined method that will be executed at the end of a call to
         * NOX::NLN::Group::computeF().
         *
         * \param F        : full access to the right hand side vector of the NOX::NLN::Group.
         * \param grp      : read only access to the NOX::NLN::Group object.
         */
        virtual void runPostComputeF(Epetra_Vector& F, const NOX::NLN::Group& grp);

        /** User defined method that will be executed at the start of a call to
         * NOX::NLN::Group::computeX().
         *
         * \param input_grp: read only access to the input group (holds the old X).
         * \param dir      : read only access to the direction vector (step length equal 1.0).
         * \param step     : read only access to the current step length (line search).
         * \param curr_grp : read only access to the called/current group (will hold the new X).
         */
        virtual void runPreComputeX(const NOX::NLN::Group& input_grp, const Epetra_Vector& dir,
            const double& step, const NOX::NLN::Group& curr_grp);

        /** User defined method that will be executed at the end of a call to
         * NOX::NLN::Group::computeX().
         *
         * \param input_grp: read only access to the input group (holds the old X).
         * \param dir      : read only access to the direction vector (step length equal 1.0).
         * \param step     : read only access to the current step length (line search).
         * \param curr_grp : read only access to the called/current group (holds the new X).
         */
        virtual void runPostComputeX(const NOX::NLN::Group& input_grp, const Epetra_Vector& dir,
            const double& step, const NOX::NLN::Group& curr_grp);

        /*! User defined method that will be executed at the beginning
         *  of a call to NOX::NLN::Group::applyJacobianInverse().
         *
         *  \param rhs    : read-only access to the rhs vector
         *  \param result : full access to the result vector
         *  \param xold   : read-only access to the state vector
         *  \param grp    : read only access to the group object
         */
        virtual void runPreApplyJacobianInverse(const NOX::Abstract::Vector& rhs,
            NOX::Abstract::Vector& result, const NOX::Abstract::Vector& xold,
            const NOX::NLN::Group& grp);

        /*! User defined method that will be executed at the end
         *  of a call to NOX::NLN::Group::applyJacobianInverse().
         *
         *  \param rhs    : read-only access to the rhs vector
         *  \param result : full access to the result vector
         *  \param xold   : read-only access to the old state vector
         *  \param grp    : read only access to the group object
         */
        virtual void runPostApplyJacobianInverse(const NOX::Abstract::Vector& rhs,
            NOX::Abstract::Vector& result, const NOX::Abstract::Vector& xold,
            const NOX::NLN::Group& grp);

       protected:
        //! Flag that determines if a pre/post operator has been supplied by user.
        bool havePrePostOperator_;

        //! Points to user defined pre/post operator for the linear system.
        Teuchos::RCP<Map> prePostOperatorMapPtr_;
      };  // class PrePostOperator
      namespace PrePostOp
      {
        // non-member function
        /*! Returns the inherent pre/post operator std::map of the "Group Options" sublist.
         *  If the corresponding parameter called "User Defined Pre/Post Operator" is not yet
         *  defined, a empty std::map is generated and set into the parameter list first. */
        NOX::NLN::GROUP::PrePostOperator::Map& GetMutableMap(Teuchos::ParameterList& p_grp_opt);
      }  // namespace PrePostOp
    }    // namespace GROUP
  }      // namespace NLN
}  // namespace NOX

inline void NOX::NLN::GROUP::PrePostOperator::runPreComputeF(
    Epetra_Vector& F, const NOX::NLN::Group& grp)
{
  if (havePrePostOperator_)
  {
    Map::iterator it;
    for (it = prePostOperatorMapPtr_->begin(); it != prePostOperatorMapPtr_->end(); ++it)
      it->second->runPreComputeF(F, grp);
  }
}

inline void NOX::NLN::GROUP::PrePostOperator::runPostComputeF(
    Epetra_Vector& F, const NOX::NLN::Group& grp)
{
  if (havePrePostOperator_)
  {
    Map::iterator it;
    for (it = prePostOperatorMapPtr_->begin(); it != prePostOperatorMapPtr_->end(); ++it)
      it->second->runPostComputeF(F, grp);
  }
}

inline void NOX::NLN::GROUP::PrePostOperator::runPreComputeX(const NOX::NLN::Group& input_grp,
    const Epetra_Vector& dir, const double& step, const NOX::NLN::Group& curr_grp)
{
  if (havePrePostOperator_)
  {
    Map::iterator it;
    for (it = prePostOperatorMapPtr_->begin(); it != prePostOperatorMapPtr_->end(); ++it)
      it->second->runPreComputeX(input_grp, dir, step, curr_grp);
  }
}

inline void NOX::NLN::GROUP::PrePostOperator::runPostComputeX(const NOX::NLN::Group& input_grp,
    const Epetra_Vector& dir, const double& step, const NOX::NLN::Group& curr_grp)
{
  if (havePrePostOperator_)
  {
    Map::iterator it;
    for (it = prePostOperatorMapPtr_->begin(); it != prePostOperatorMapPtr_->end(); ++it)
      it->second->runPostComputeX(input_grp, dir, step, curr_grp);
  }
}

inline void NOX::NLN::GROUP::PrePostOperator::runPreApplyJacobianInverse(
    const NOX::Abstract::Vector& rhs, NOX::Abstract::Vector& result,
    const NOX::Abstract::Vector& xold, const NOX::NLN::Group& grp)
{
  if (havePrePostOperator_)
  {
    Map::iterator it;
    for (it = prePostOperatorMapPtr_->begin(); it != prePostOperatorMapPtr_->end(); ++it)
      it->second->runPreApplyJacobianInverse(rhs, result, xold, grp);
  }
}

inline void NOX::NLN::GROUP::PrePostOperator::runPostApplyJacobianInverse(
    const NOX::Abstract::Vector& rhs, NOX::Abstract::Vector& result,
    const NOX::Abstract::Vector& xold, const NOX::NLN::Group& grp)
{
  if (havePrePostOperator_)
  {
    Map::iterator it;
    for (it = prePostOperatorMapPtr_->begin(); it != prePostOperatorMapPtr_->end(); ++it)
      it->second->runPostApplyJacobianInverse(rhs, result, xold, grp);
  }
}

#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_GROUP_PREPOSTOPERATOR_H_ */
