/*
 * nox_nln_group_prepostoperator.H
 *
 *  Created on: Oct 19, 2015
 *      Author: hiermeier
 */

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_GROUP_PREPOSTOPERATOR_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_GROUP_PREPOSTOPERATOR_H_

#include <Teuchos_RCP.hpp>
#include "nox_nln_abstract_prepostoperator.H"


// forward declarations
class Epetra_Vector;
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
namespace NOX {
namespace NLN {
namespace GROUP {
/*!
  @brief Functor to process the pre/post operator object in the parameter list for the NOX::NLN::Group objects.

  This is a wrapper class for a user derived  NOX::NLN::Abstract::PrePostOperator (ppo) object.
  All NOX::NLN groups use this class so we don't have to repeat all parsing code in each
  NOX::NLN group class. This class searches the "Group Options" parameter list passed into
  the constructor and if a ppo is found will wrap the object.

  For instructions on how to implement a PrePostOperator, see
  NOX::NLN::Abstract::PrePostOperator.

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
    PrePostOperator(Teuchos::ParameterList& groupOptionsSublist);

    //! Destructor.
    virtual ~PrePostOperator() {};

    //! Resets the pre/post operator.
    virtual void reset(Teuchos::ParameterList& groupOptionsSublist);

    /** User defined method that will be executed at the start of a call to NOX::NLN::Group::computeF().
     *
     * \param F        : full access to the right hand side vector of the NOX::NLN::Group.
     * \param grp      : read only access to the NOX::NLN::Group object.
     */
    virtual void runPreComputeF(Epetra_Vector& F,
        const NOX::NLN::Group& grp);

    /** User defined method that will be executed at the end of a call to NOX::NLN::Group::computeF().
     *
     * \param F        : full access to the right hand side vector of the NOX::NLN::Group.
     * \param grp      : read only access to the NOX::NLN::Group object.
     */
    virtual void runPostComputeF( Epetra_Vector& F,
        const NOX::NLN::Group& grp);

  protected:
    //! Flag that determines if a pre/post operator has been supplied by user.
    bool havePrePostOperator_;

    //! Points to user defined pre/post operator for the linear system.
    Teuchos::RCP<NOX::NLN::Abstract::PrePostOperator> prePostOperatorPtr_;
  };  // class PrePostOperator
} // namespace GROUP
} // namespace NLN
} // namespace NOX

inline void NOX::NLN::GROUP::PrePostOperator::runPreComputeF(
    Epetra_Vector& F,
    const NOX::NLN::Group& grp)
{
  if (havePrePostOperator_)
    prePostOperatorPtr_->runPreComputeF(F,grp);
}

inline void NOX::NLN::GROUP::PrePostOperator::runPostComputeF(
    Epetra_Vector& F,
    const NOX::NLN::Group& grp)
{
  if (havePrePostOperator_)
    prePostOperatorPtr_->runPostComputeF(F,grp);
}


#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_GROUP_PREPOSTOPERATOR_H_ */
