/*---------------------------------------------------------------------*/
/*!
\file str_model_evaluator_contact.H

\brief Evaluation and assembly of all contact terms

\maintainer Michael Hiermeier

\date Feb 3, 2016

\level 3

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_CONTACT_H_
#define SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_CONTACT_H_

#include "str_model_evaluator_generic.H"

// forward declarations...
namespace CONTACT {
  class Manager;
  class CoAbstractStrategy;
} // namespace CONTACT
namespace MORTAR {
  class StrategyBase;
} // namespace MORTAR

namespace STR {
namespace MODELEVALUATOR {
  class ContactData;

  class Contact : public Generic
  {
  public:
    //! constructor
    Contact();

    //! destructor
    virtual ~Contact() {};

    //! setup class variables [derived]
    virtual void Setup();

    //! @name Functions which are derived from the base generic class
    //! @{
    //! [derived]
    virtual bool ApplyForce(
        const Epetra_Vector& x,
        Epetra_Vector& f);

    //! [derived]
    virtual bool ApplyStiff(
        const Epetra_Vector& x,
        LINALG::SparseOperator& jac);

    //! [derived]
    virtual bool ApplyForceStiff(
        const Epetra_Vector& x,
        Epetra_Vector& f,
        LINALG::SparseOperator& jac);

    //! [derived]
    virtual void WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const;

    //! [derived]
    virtual void ReadRestart(IO::DiscretizationReader& ioreader);

    //! recover condensed Lagrange multipliers
    virtual void RecoverState(
        const Epetra_Vector& xold,
        const Epetra_Vector& dir,
        const Epetra_Vector& xnew);

    //! [derived]
    virtual void UpdateStepState();

    //! [derived]
    virtual void UpdateStepElement();

    //! [derived]
    virtual void DetermineStressStrain();

    //! [derived]
    virtual void DetermineEnergy();

    //! [derived]
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    //! [derived]
    virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const;

    //! [derived]
    virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const;

    //! [derived]
    virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const;
    //! @}

    //! @name Accessors
    //! @{

    const Teuchos::RCP<CONTACT::CoAbstractStrategy>& StrategyPtr();
    //! @}
  protected:
    //! Returns the underlying contact strategy object
    CONTACT::CoAbstractStrategy& Strategy();
    const CONTACT::CoAbstractStrategy& Strategy() const;

    //! reset class variables (without jacobian) [derived]
    virtual void Reset(const Epetra_Vector& x);

    //! reset jacobian and call Reset(x) [derived]
    virtual void Reset(const Epetra_Vector& x, LINALG::SparseOperator& jac);

    void AssembleRhs(Epetra_Vector& f) const;

    void AssembleJacobian(LINALG::SparseOperator& jac) const;

    STR::MODELEVALUATOR::ContactData& EvalContact();

  private:
    void PostUpdateStepState();

    void CheckPseudo2D() const;

    void PrintBanner() const;

  private:
    //! contact evaluation data container
    Teuchos::RCP<STR::MODELEVALUATOR::ContactData> eval_contact_ptr_;

    //! contact strategy
    Teuchos::RCP<CONTACT::CoAbstractStrategy> strategy_ptr_;
  }; // class Contact
} // namespace MODELEVALUATOR
} // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_CONTACT_H_ */
