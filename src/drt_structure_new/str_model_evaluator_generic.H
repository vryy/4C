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
#include "../drt_inpar/inpar_structure.H"

// forward declaration
class Epetra_Vector;
class Epetra_Map;
namespace NOX {
namespace Solver {
  class Generic;
} // namespace Solver
namespace NLN {
  class Group;
} // namespace NLN
} // namespace NOX
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
  class DiscretizationInterface;
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
  /*! \brief This is the abstract base class of all model evaluators
   *
   *  This class summarizes the functionality which all model evaluators share
   *  and/or have to implement. Look in the derived classes for examples. A minimal
   *  example can be found at \ref STR::MODELEVALUATOR::PartitionedFSI.
   *
   *  \date 08/15
   *  \author hiermeier */
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
    //! Returns true, if Init() has been called
    inline const bool& IsInit() const { return isinit_; };

    //! Returns true, if Setup() has been called
    inline const bool& IsSetup() const { return issetup_; };

    //! Checks, if Init() and Setup() have been called
    virtual void CheckInitSetup() const;

    virtual void CheckInit() const;

  public:
    //! Returns the type of the current model evaluator
    virtual INPAR::STR::ModelType Type() const = 0;

    /*! \brief reset model specific variables (without jacobian)
     *
     *  This function is always called before the actual Evaluate() routine is going to start.
     *  You can use it to reset model specific stuff at the beginning of a new evaluation round.
     *
     *  \param x   (in) : current full state vector
     *
     *  \date 07/2016
     *  \author hiermeier */
    virtual void Reset(const Epetra_Vector& x) = 0;

    /*! \brief Evaluate the current right-hand-side at \f$t_{n+1}\f$
     *
     *  \date 07/2016
     *  \author hiermeier */
    virtual bool EvaluateForce() = 0;

    /*! \brief Evaluate the initial right hand side (overload if needed for specific model)
     *
     *  \date 07/2016
     *  \author hiermeier */
    virtual bool EvaluateInitialForce(){
      return EvaluateForce();
    };

    /*! \brief Evaluate the current tangential stiffness matrix at \f$t_{n+1}\f$
     *
     *  \date 07/2016
     *  \author hiermeier*/
    virtual bool EvaluateStiff() = 0;

    /*! \brief Evaluate the current right-hand-side vector and tangential stiffness matrix at \f$t_{n+1}\f$
     *
     *  \date 07/2016
     *  \author hiermeier */
    virtual bool EvaluateForceStiff() = 0;

    /*! \brief Called at the end of each call to one of the STR::ModelEvaluator::Evaluate routines,
     *         i.e. EvaluateForce, EvaluateStiff, EvaluateForceStiff.
     *
     *  \author hiermeier \date 03/17 */
    virtual void PostEvaluate() = 0;

    /*! \brief Assemble the force right-hand-side
     *
     *  After the evaluation of all models at the new state \f$t_{n+1}\f$ is finished we
     *  start to put everything together. At this point we build the right-hand-side at
     *  the desired mid-point. The reason why we do it here and not at time integrator level
     *  is that you have now the possibility to use a different time integration (not the underlying
     *  structural one). This makes it more flexible. However, the contributions of the old
     *  time step should be stored in the UpdateStepState() routine. There you can scale the
     *  old contributions with the time integration factor you like and save them once, thus you just
     *  have to add them to complete the mid-point right-hand-side.
     *
     *  \param f         (out) : Time factor scaled right hand side w/o viscous and/or mass effects
     *  \param timefac_np (in) : Time factor of the underlying structural time integrator for the
     *                           new state at \f$t_{n+1}\f$. (To scale the old state of the previous
     *                           time step see UpdateStepState.)
     *
     *  \date 07/2016
     *  \author hiermeier */
    virtual bool AssembleForce(Epetra_Vector& f,
        const double & timefac_np) const = 0;

    //! Assemble the jacobian
    virtual bool AssembleJacobian(
        LINALG::SparseOperator& jac,
        const double & timefac_np) const = 0;

    /*! \brief write model specific restart
     *
     *  \param iowriter            (in) : output writer
     *  \param forced_writerestart (in) : special treatment is necessary, if the restart is forced */
    virtual void WriteRestart(
        IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const = 0;

    /*! \brief read model specific restart information
     *
     *  \param ioreader (in) : input reader
     *
     *  \date 07/2016
     *  \author hiermeier */
    virtual void ReadRestart(IO::DiscretizationReader& ioreader) = 0;

    /*! \brief Predict the values for DoFs that are defined in
     *         the respective model evaluators, e.g. condensed variables.*/
    virtual void Predict(const INPAR::STR::PredEnum& pred_type) = 0;

    /*! \brief Recover condensed solution variables
     *
     *  This method is supposed to be used to recover condensed solution variables.
     *  Typical examples are the EAS degrees of freedom or the dual Lagrange multipliers.
     *  Do NOT use it to reset your model variables! Use the Reset() method instead.
     *
     *  \author hiermeier \date 07/2016 */
    virtual void RecoverState(
        const Epetra_Vector& xold,
        const Epetra_Vector& dir,
        const Epetra_Vector& xnew) = 0;

    /*! \brief Executed before the solution vector is going to be updated
     *
     *  \author hiermeier \date 03/17 */
    virtual void RunPreComputeX(
        const Epetra_Vector& xold,
        Epetra_Vector& dir_mutable,
        const NOX::NLN::Group& curr_grp ) = 0;

    /*! \brief Executed at the end of the NOX::Solver::Generic::Step() (f.k.a. Iterate()) method
     *
     *  \param solver (in) : reference to the non-linear nox solver object (read-only)
     *
     *  \author hiermeier \date 03/17 */
    virtual void RunPostIterate( const NOX::Solver::Generic& solver ) = 0;


    virtual void UpdateStepState(const double& timefac_n) = 0;

    virtual void UpdateStepElement() = 0;

    /*! \brief calculate the stress/strain contributions of each model evaluator
     *
     *  \remark This function is called from STR::TIMINT::Base::PrepareOutput() and calculates
     *  missing quantities, which were not evaluated during the standard evaluate call and are
     *  only calculated once per load/time step. You can not do the calculations during the
     *  OutputStepState() routine, because of the const status of the named function!
     *
     *  \author hiermeier*/
    virtual void DetermineStressStrain() = 0;

    /*! \brief calculate energy contributions of each model evaluator
     *
     *  \remark This function is called from STR::TIMINT::Base::PrepareOutput() and calculates
     *  missing quantities, which were not evaluated during the standard evaluate call and are
     *  only calculated once per load/time step. You can not do the calculations during the
     *  OutputStepState() routine, because of the const status of the named function!
     *
     *  \author hiermeier*/
    virtual void DetermineEnergy() = 0;

    //! output routine for model evaluator
    virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const = 0;

    //! runtime output routine for model evaluator
    virtual void RuntimeOutputStepState() const {return;};

    //! reset routine for model evlaluator
    virtual void ResetStepState() = 0;

    //! post output routine for model evlaluator
    virtual void PostOutput() = 0;

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
    /*! \brief Check the evaluation procedures for possible errors.
     *
     *  In the standard case we check for exceptions like overflow, invalid results
     *  or divide by zero operations. Furthermore, we look for an (optional) parameter
     *  named ele_eval_error_flag_. This is universal and should be usable by all model
     *  evaluators.
     *
     *  \author hiermeier*/
    virtual bool EvalErrorCheck() const;

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
    Teuchos::RCP<STR::TIMINT::BaseDataIO> GInOutputPtr();
    const STR::TIMINT::BaseDataIO& GInOutput() const;

    //! Returns the (structural) discretization
    DRT::DiscretizationInterface& Discret();
    Teuchos::RCP<DRT::DiscretizationInterface>& DiscretPtr();
    const DRT::DiscretizationInterface& Discret() const;

    //! Returns the underlying STR::Integrator object
    STR::Integrator& Int();
    const STR::Integrator& Int() const;
    Teuchos::RCP<STR::Integrator>& IntPtr();

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

    //! pointer to the problem discretization
    Teuchos::RCP<DRT::DiscretizationInterface> discret_ptr_;

    //! pointer to the structural (time) integrator
    Teuchos::RCP<STR::Integrator> int_ptr_;

    //! pointer to the time integrator strategy object
    Teuchos::RCP<const STR::TIMINT::Base> timint_ptr_;

    /*! \brief initial dof offset
     *
     *  This variable becomes important when you plan to create a
     *  block for the system of equations. E.g. saddle-point system
     *  in contact case. */
    int dof_offset_;

  }; // class Generic

} // namespace MODELEVALUATOR
} // namespace STR

#endif /* STR_MODEL_EVALUATOR_GENERIC_H_ */
