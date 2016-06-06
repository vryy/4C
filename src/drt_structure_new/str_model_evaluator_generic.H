/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_generic.H

\brief Generic class for all model evaluators.

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
class Epetra_Map;
namespace LINALG {
  class SparseOperator;
} // namespace LINALG
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
namespace IO {
  class DiscretizationWriter;
  class DiscretizationReader;
} // namespace IO
namespace DRT {
  class Discretization;
} // namespace DRT
namespace STR {
  class Integrator;
namespace TIMINT {
  class BaseDataGlobalState;
  class BaseDataIO;
  class Base;
} // namespace TIMINT
namespace MODELEVALUATOR {
  class Data;
  class Generic
  {
  public:
    //! constructor
    Generic();

    //! destructor
    virtual ~Generic() {};

    //! initialize the class variables
    virtual void Init(
        const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
        const Teuchos::RCP<STR::Integrator>& int_ptr,
        const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr,
        const int& dof_offset);

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
     *  or divide by zero operations. Furthermore, we look for an (optional) parameter
     *  named ele_eval_error_flag_. This is universal and should be usable by all model
     *  evaluators.
     *
     *  \author hiermeier*/
    virtual bool EvalErrorCheck() const;

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

    virtual void WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const = 0;

    virtual void ReadRestart(IO::DiscretizationReader& ioreader) = 0;

    virtual void RecoverState(
        const Epetra_Vector& xold,
        const Epetra_Vector& dir,
        const Epetra_Vector& xnew) = 0;

    virtual void UpdateStepState() = 0;

    virtual void UpdateStepElement() = 0;

    virtual void DetermineStressStrain() = 0;

    virtual void DetermineEnergy() = 0;

    //! output routine for model evlaluator
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const = 0;

    //! @name Accessors to model specific things
    //! @{

    //! Returns a pointer to the model specific dof row map
    virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const = 0;

    //! Returns a pointer to the current model solution vector (usually the Lagrange multiplier vector)
    virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const = 0;

    //! Returns a pointer to the model solution vector of the last time step (usually the Lagrange multiplier vector)
    virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const = 0;
    //! @}
  protected:
    //! @name internal accessors
    //! @{
    //! Returns the model evaluator data container
    STR::MODELEVALUATOR::Data& EvalData();
    Teuchos::RCP<STR::MODELEVALUATOR::Data>& EvalDataPtr();

    //! Returns the global state data container
    STR::TIMINT::BaseDataGlobalState& GState();
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& GStatePtr();
    const STR::TIMINT::BaseDataGlobalState& GState() const;

    //! Returns the global input/output data container
    STR::TIMINT::BaseDataIO& GInOutput();
    const STR::TIMINT::BaseDataIO& GInOutput() const;

    //! Returns the (structural) discretization
    DRT::Discretization& Discret();
    Teuchos::RCP<DRT::Discretization>& DiscretPtr();
    const DRT::Discretization& Discret() const;

    //! Returns the underlying STR::Integrator object
    STR::Integrator& Int();

    //! Returns the underlying STR::TIMINT object
    const STR::TIMINT::Base& TimInt() const;

    const int& DofOffset() const;
    //! @}
  protected:
    //! init flag
    bool isinit_;

    //! setup flag
    bool issetup_;
  private:
    //! pointer to the model evaluator data container
    Teuchos::RCP<STR::MODELEVALUATOR::Data> eval_data_ptr_;

    //! pointer to the global state data container
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

    //! pointer to input/ouput data container
    Teuchos::RCP<STR::TIMINT::BaseDataIO> gio_ptr_;

    //! problem discretization
    Teuchos::RCP<DRT::Discretization> discret_ptr_;

    //! pointer to the structural (time) integrator
    Teuchos::RCP<STR::Integrator> int_ptr_;

    //! pointer to the time integrator strategy object
    Teuchos::RCP<const STR::TIMINT::Base> timint_ptr_;

    int dof_offset_;
  }; // class Generic

} // namespace MODELEVALUATOR
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_GENERIC_H_ */
