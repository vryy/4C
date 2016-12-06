/*-----------------------------------------------------------*/
/*!
\file str_integrator.H

\brief Generic class for all explicit/implicit time integrators.

\maintainer Michael Hiermeier

\date Dec 7, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_STR_INTEGRATOR_H_
#define SRC_DRT_STRUCTURE_NEW_STR_INTEGRATOR_H_

#include <Teuchos_RCP.hpp>
#include "../drt_inpar/inpar_structure.H" // enumerators
#include "../solver_nonlin_nox/nox_nln_enum_lists.H" // enumerators
#include <NOX_Abstract_Vector.H> // enumerators

// forward declaration
class Epetra_Vector;
namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
namespace IO {
  class DiscretizationWriter;
  class DiscretizationReader;
} // namespace IO
namespace STR {
  class ModelEvaluator;
  class Dbc;
namespace MODELEVALUATOR {
  class Data;
  class Generic;
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
    virtual void Setup();

    //! Set state variables
    virtual void SetState(const Epetra_Vector& x) = 0;

    /*! \brief Reset state variables of all models
     *  (incl. the structural dynamic state variables)
     *
     *  \param x (in) : current full state vector */
    void ResetModelStates(const Epetra_Vector& x);

    //! Apply the right hand side only
    virtual bool ApplyForce(const Epetra_Vector& x,
        Epetra_Vector& f) = 0;

    //! return integration factor
    virtual double GetIntParam() const = 0;

    //! @name Restart and output related functions
    //! @{

    /*! write restart information of the different time integration schemes
     *  and model evaluators */
    void WriteRestart(IO::DiscretizationWriter& iowriter) const
    { WriteRestart(iowriter,false); };
    virtual void WriteRestart(IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const = 0;

    /*! read restart information of the different time integration schemes
     *  and model evaluators */
    virtual void ReadRestart(IO::DiscretizationReader& ioreader) = 0;
    //! @}

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
    void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    /** \brief reset step configuration after time step
     *
     *  This function is supposed to reset all variables which are directly related
     *  to the current new step n+1. To be more precise all variables ending with "Np"
     *  have to be reseted. */
    virtual void ResetStepState();

    /// things that should be done after updating
    virtual void PostUpdate() = 0;

    /// things that should be done after output
    virtual void PostOutput();
    //!@}

    //! @name Accessors
    //! @{

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
    Teuchos::RCP<const STR::ModelEvaluator> ModelEvalPtr() const;

    //! Return the model evaluator control object (read-write)
    STR::ModelEvaluator& ModelEval();

    //! Return the model evaluator object for the given model type
    STR::MODELEVALUATOR::Generic& Evaluator(const INPAR::STR::ModelType& mt);

    //! Return the model evaluator object for the given model type
    const STR::MODELEVALUATOR::Generic& Evaluator(const INPAR::STR::ModelType& mt) const;

    //! Return the model evaluator data object (read-only)
    const STR::MODELEVALUATOR::Data& EvalData() const;

    //! Return the model evaluator data object (read and write access)
    STR::MODELEVALUATOR::Data& EvalData();

    //! Return the Dirichlet boundary condition object (read-only)
    const STR::Dbc& GetDbc() const;
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

    /*! \brief Equilibriate system at initial state and identify consistent accelerations
     *
     *  \author hiermeier */
    void EquilibriateInitialState();

    /*! check if current state is equilibrium (with respect to
     *  a given tolerance of the inf-norm)
     *  this is a sanity check in case of nonlinear mass
     *  and non-additive rotation pseudo-vector DoFs where
     *  determination of consistent accelerations is more intricate
     *  and not supported yet */
    bool CurrentStateIsEquilibrium(const double& tol);

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

    //! pointer to the input/output data container
    Teuchos::RCP<STR::TIMINT::BaseDataIO> io_ptr_;

    //! pointer to the dirichlet boundary condition object
    Teuchos::RCP<STR::Dbc> dbc_ptr_;

    //! pointer to the underlying time integrator (read-only)
    Teuchos::RCP<const STR::TIMINT::Base> timint_ptr_;
  };
} // namespace STR


#endif /* SRC_DRT_STRUCTURE_NEW_STR_INTEGRATOR_H_ */
