/*-----------------------------------------------------------*/
/*!
\file str_integrator.H

\maintainer Michael Hiermeier

\date Dec 7, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_INTEGRATOR_H_
#define SRC_DRT_STRUCTURE_NEW_STR_INTEGRATOR_H_

#include <Teuchos_RCP.hpp>

// forward declaration
class Epetra_Vector;
namespace STR {
  class ModelEvaluator;
namespace TIMINT {
  class Base;
  class BaseDataGlobalState;
} // namespace TIMINT
  class Integrator
  {
  public:
    //! constructor
    Integrator();

    //! destructor
    virtual ~Integrator() {};

    //! Initialization
    virtual void Init(const Teuchos::RCP<STR::ModelEvaluator>& me_ptr);

    //! Setup (has to be implemented by the derived classes)
    virtual void Setup() = 0;

    //! Set state variables
    virtual void SetState(const Epetra_Vector& x) = 0;

    //! Apply the right hand side only
    virtual bool ApplyForce(const Epetra_Vector& x,
        Epetra_Vector& f) = 0;

    //! return integration factor
    virtual double GetIntParam() const = 0;

    //! @name Monolithic update routines
    //!@{
    /// things that should be done before updating
    virtual void PreUpdate() = 0;

    /*! \brief Update configuration after time step
     *
     *  Thus the 'last' converged is lost and a reset of the time step
     *  becomes impossible. We are ready and keen awaiting the next
     *  time step. */
    virtual void UpdateStepState() = 0;

    /*! \brief Update everything on element level after time step and after output
     *
     *  Thus the 'last' converged is lost and a reset of the time step
     *  becomes impossible. We are ready and keen awaiting the next
     *  time step. */
    virtual void UpdateStepElement() = 0;

    /// things that should be done after updating
    virtual void PostUpdate() = 0;
    //!@}

  protected:
    //! returns init state
    inline const bool& IsInit() const
    { return isinit_; };

    //! returns setup state
    inline const bool& IsSetup() const
    { return issetup_; };

    //! Check the init state
    void CheckInit() const;

    //! Check the setup state
    void CheckInitSetup() const;

    STR::ModelEvaluator& ModelEval();

    const STR::ModelEvaluator& ModelEval() const;

    STR::TIMINT::BaseDataGlobalState& GlobalState();

    const STR::TIMINT::BaseDataGlobalState& GlobalState() const;

    const STR::TIMINT::Base& TimInt() const;

  protected:
    //! indicates if the Init() function has been called
    bool isinit_;

    //! indicates if the Setup() function has been called
    bool issetup_;

  private:
    //! pointer to the model evaluator
    Teuchos::RCP<STR::ModelEvaluator> modelevaluator_ptr_;

    //! pointer to the global state data
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

    //! pointer to the underlying time integrator (read-only)
    Teuchos::RCP<const STR::TIMINT::Base> timint_ptr_;
  };
} // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_INTEGRATOR_H_ */
