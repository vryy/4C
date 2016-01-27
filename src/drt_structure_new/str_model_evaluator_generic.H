/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_generic.H

\maintainer Michael Hiermeier

\date Aug 11, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_MODEL_EVALUATOR_GENERIC_H_
#define STR_MODEL_EVALUATOR_GENERIC_H_

#include <Teuchos_RCP.hpp>

// forward declaration
class Epetra_Vector;
namespace LINALG {
  class SparseOperator;
} // namespace LINALG
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
namespace DRT {
  class Discretization;
} // namespace DRT
namespace STR {
namespace TIMINT {
  class BaseDataGlobalState;
  class BaseDataIO;
  class Base;
} // namespace TIMINT
namespace MODELEVALUATOR {
  class Generic
  {
  public:
    //! constructor
    Generic();

    //! destructor
    virtual ~Generic() {};

    //! initialize the class variables
    virtual void Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
        const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr);

    //! setup class variables
    virtual void Setup() = 0;

  protected:
    //! reset class variables (without jacobian)
    virtual void Reset(const Epetra_Vector& x) = 0;
    //! reset jacobian and call Reset(x)
    virtual void Reset(const Epetra_Vector& x,LINALG::SparseOperator& jac) = 0;

    //! Returns true, if Init() has been called
    inline const bool& IsInit() const { return isinit_; };

    //! Returns true, if Setup() has been called
    inline const bool& IsSetup() const { return issetup_; };

    //! Checks, if Init() and Setup() have been called
    void CheckInitSetup() const;

    void CheckInit() const;

    /*! \brief Check the evaluation procedures for possible errors.
     *
     *  In the standard case we check for exceptions like overflow, invalid results
     *  or divide by zero operations. Furthermore, we look in the parameter-list for
     *  an (optional) parameter named "eval_error". This is universal and should be
     *  usable by all model evaluators. */
    virtual bool EvalErrorCheck(const Teuchos::ParameterList& p) const;

  public:

    virtual bool ApplyForce(
        const Epetra_Vector& x,
        Epetra_Vector& f) = 0;

    virtual bool ApplyStiff(
        const Epetra_Vector& x,
        LINALG::SparseOperator& jac) = 0;

    virtual bool ApplyForceStiff(
        const Epetra_Vector& x,
        Epetra_Vector& f,
        LINALG::SparseOperator& jac) = 0;

    virtual void UpdateStepState() = 0;

    virtual void UpdateStepElement() = 0;

    //! output routine for model evlaluator
    virtual void OutputStepState() = 0;


  protected:
    //! init flag
    bool isinit_;

    //! setup flag
    bool issetup_;

    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

    Teuchos::RCP<STR::TIMINT::BaseDataIO> gio_ptr_;

    //! problem discretization
    Teuchos::RCP<DRT::Discretization> discret_ptr_;

    Teuchos::RCP<const STR::TIMINT::Base> timint_ptr_;
  }; // class Generic

} // namespace MODELEVALUATOR
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_GENERIC_H_ */
