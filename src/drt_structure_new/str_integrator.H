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
#include "../drt_inpar/inpar_structure.H" // enumerators
#include "../solver_nonlin_nox/nox_nln_statustest_factory.H" // enumerators
#include <NOX_Abstract_Vector.H> // enumerators

// forward declaration
class Epetra_Vector;
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
namespace STR {
  class ModelEvaluator;
  class Dbc;
namespace MODELEVALUATOR {
  class Data;
} // namespace MODELEVALUATOR
namespace TIMINT {
  class Base;
  class BaseDataSDyn;
  class BaseDataGlobalState;
  class BaseDataIO;
} // namespace TIMINT
  class Integrator
  {
  public:
    //! constructor
    Integrator();

    //! destructor
    virtual ~Integrator() {};

    //! Initialization
    virtual void Init(
        const Teuchos::RCP<STR::TIMINT::BaseDataSDyn>& sdyn_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
        const Teuchos::RCP<STR::Dbc>& dbc_ptr,
        const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr
        );

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

    //! calculate stresses and strains in the different model evaluators
    void DetermineStressStrain();

    //! calculate the energy in the different model evaluators
    void DetermineEnergy();

    /*! \brief Output to file
     *
     *  This routine prints always the last converged state, i.e.
     *  \f$D_{n}, V_{n}, A_{n}\f$. So, #UpdateIncrement should be called
     *  upon object prior to writing stuff here.
     *
     *  \author mwgee (originally)                         \date 03/07 */
    void OutputStepState();

    /// things that should be done after updating
    virtual void PostUpdate() = 0;
    //!@}

    //! @name Accessors
    //! @{
    //! return the EquilibriateInitalState flag
    const bool& IsEquilibriateInitialState() const
    {
      CheckInit();
      return isequalibriate_initial_state_;
    }

    double GetCondensedUpdateNorm(
        const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

    double GetCondensedPreviousSolNorm(
        const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

    double GetCondensedSolutionUpdateRMS(
        const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

    int GetCondensedDofNumber(
        const enum NOX::NLN::StatusTest::QuantityType& qtype) const;

    //! Return the model evaluator control object (read-only)
    const STR::ModelEvaluator& ModelEval() const;

    //! Return the model evaluator data object (read-only)
    const STR::MODELEVALUATOR::Data& EvalData() const;
    //! @}

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

    /*! Equilibriate system at initial state and
     *  identify consistent accelerations */
    void EquilibriateInitialState();

    //! Return the model evaluator control object
    STR::ModelEvaluator& ModelEval();

    //! Return the model evaluator data object
    STR::MODELEVALUATOR::Data& EvalData();

    //! Return the structural dynamic data container
    STR::TIMINT::BaseDataSDyn& SDyn();

    //! Return the structural dynamic data container (read-only)
    const STR::TIMINT::BaseDataSDyn& SDyn() const;

    //! Return the global state data container
    STR::TIMINT::BaseDataGlobalState& GlobalState();

    //! Return the global state data container (read-only)
    const STR::TIMINT::BaseDataGlobalState& GlobalState() const;

    //! Return the Dirichlet boundary condition object
    STR::Dbc& Dbc();

    //! Return the Dirichlet boundary condition object (read-only)
    const STR::Dbc& Dbc() const;

    //! Return the time integration strategy object (read-only)
    const STR::TIMINT::Base& TimInt() const;

    //! reset the time step dependent parameters for the element evaluation
    virtual void ResetEvalParams(){ };

    double GetCondensedGlobalNorm(
        const enum NOX::NLN::StatusTest::QuantityType& qtype,
        const enum NOX::Abstract::Vector::NormType& normtype,
        double& mynorm) const;

  protected:
    //! indicates if the Init() function has been called
    bool isinit_;

    //! indicates if the Setup() function has been called
    bool issetup_;

  private:
    //! pointer to the model evaluator
    Teuchos::RCP<STR::ModelEvaluator> modelevaluator_ptr_;

    //! pointer to model evaluator data
    Teuchos::RCP<STR::MODELEVALUATOR::Data> eval_data_ptr_;

    //! pointer to the structural dynamic data container
    Teuchos::RCP<STR::TIMINT::BaseDataSDyn> sdyn_ptr_;

    //! pointer to the global state data container
    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

    //! pointer to the dirichlet boundary condition object
    Teuchos::RCP<STR::Dbc> dbc_ptr_;

    //! pointer to the underlying time integrator (read-only)
    Teuchos::RCP<const STR::TIMINT::Base> timint_ptr_;

    //! indicates that the EquilibriateInitialState() routine is active
    bool isequalibriate_initial_state_;
  };
} // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_INTEGRATOR_H_ */
