/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator.H

\brief Manager of the model evaluator calls.

\maintainer Michael Hiermeier

\date Nov 30, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_H_
#define SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_H_

#include <Teuchos_RCP.hpp>
// necessary due to enums
#include "../drt_inpar/inpar_structure.H"

// forward declaration ...
class Epetra_Vector;
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
namespace IO {
  class DiscretizationWriter;
  class DiscretizationReader;
} // namespace IO
namespace STR {
  class Integrator;
namespace TIMINT {
  class Base;
  class BaseDataSDyn;
  class BaseDataGlobalState;
  class BaseDataIO;
} // namespace TIMINT
namespace MODELEVALUATOR {
  class Data;
  class Generic;
} // namespace MODELEVALUATOR
  //! Wrapper class for the STR::MODELEVALUATOR::Generic derived objects.
  class ModelEvaluator
  {
  public:
    typedef std::map<enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > Map;
    typedef std::vector<Teuchos::RCP<STR::MODELEVALUATOR::Generic> > Vector;

    //! constructor
    ModelEvaluator();

    //! destructor
    virtual ~ModelEvaluator() {};

    //! initialize
    void Init(const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataSDyn>& sdyn_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
        const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
        const Teuchos::RCP<STR::Integrator>& int_ptr,
        const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr);

    //! setup
    void Setup();

    //! setup the MultiMapExtractor in the global state
    void SetupMultiMapExtractor();

    //! @name General evaluate routines
    //! @{
    bool InitializeInertiaAndDamping(const Epetra_Vector& x,
        LINALG::SparseOperator& jac);

    bool ApplyInitialForce(const Epetra_Vector& x,
        Epetra_Vector& f);

    //! apply force
    bool ApplyForce(const Epetra_Vector& x, Epetra_Vector& f,
        const double& timefac_np) const;

    //! apply stiffness
    bool ApplyStiff(const Epetra_Vector& x, LINALG::SparseOperator& jac,
        const double& timefac_np) const;

    //! apply model specific stiff
    bool ApplyStiff(
        const INPAR::STR::ModelType& mt,
        const Epetra_Vector& x,
        LINALG::SparseOperator& jac,
        const double& timefac_np) const;

    //! apply force and stiffness
    bool ApplyForceStiff(const Epetra_Vector& x,
        Epetra_Vector& f,
        LINALG::SparseOperator& jac,
        const double& timefac_np) const;

    //! predict all internal variables in model evaluators
    void Predict(const INPAR::STR::PredEnum& pred_type) const;

    /** \brief Assembly of all force contributions
     *
     *  \param timefac_np           (in) : time integration factor for the current contribution
     *                                     \f$t_{n+1}\f$
     *  \param f                   (out) : force vector which is going to be assembled.
     *  \param without_these_models (in) : Assemble all models, except the models in this
     *                                     vector (optional)
     *
     *  \author hiermeier \date 03/17 */
    bool AssembleForce(
        const double timefac_np,
        Epetra_Vector& f,
        const std::vector<INPAR::STR::ModelType>* without_these_models ) const;

    /** \brief Assembly of all force contributions
     *
     *  \param timefac_np (in) : time integration factor for the current contribution \f$t_{n+1}\f$
     *  \param f         (out) : force vector which is going to be assembled.
     *
     *  \author hiermeier \date 03/17 */
    inline bool AssembleForce(
        const double timefac_np,
        Epetra_Vector& f ) const
    {
      return AssembleForce( *me_vec_ptr_, timefac_np, f );
    }

    /** \brief Assembly of a sub-set of force contributions
     *
     *  \param me_vec     (in) : user provided vector of all models which are going to be assembled
     *  \param timefac_np (in) : time integration factor for the current contribution \f$t_{n+1}\f$
     *  \param f         (out) : force vector which is going to be assembled.
     *
     *  \author hiermeier \date 03/17 */
    inline bool AssembleForce(
        const Vector& me_vec,
        const double timefac_np,
        Epetra_Vector& f ) const
    {
      bool ok = true;
      AssembleForce( ok, me_vec, timefac_np, f );
      return ok;
    }

    //! @}

    /*! \brief reset all model states (incl. the structural dynamic state)
     *
     *  \param x (in) : current full state vector */
    void ResetStates(const Epetra_Vector& x) const;

    /*! \brief reset all model states (incl. the structural dynamic state)
     *
     *  \param x (in) : current full state vector
     *  \param setstate (in) : flag to set state */
    void ResetStates(const Epetra_Vector& x, bool setstate) const;

    //! Write current restart
    void WriteRestart(IO::DiscretizationWriter& iowriter,
        const bool& forced_writerestart) const;

    //! Read restart information
    void ReadRestart(IO::DiscretizationReader& ioreader);

    //! @name Accessors
    //! @{
    //! return global state (read-only)
    const STR::TIMINT::BaseDataGlobalState& GetGlobalState() const;

    //! return global state pointer (read and write access of the data)
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& GlobalStatePtr();

    //! return pointer to the underlying time integrator (read-only)
    const Teuchos::RCP<const STR::TIMINT::Base>& GetTimIntPtr() const;

    //! access one specific model evaluator
    STR::MODELEVALUATOR::Generic& Evaluator(const enum INPAR::STR::ModelType& mt);
    const STR::MODELEVALUATOR::Generic& Evaluator(const enum INPAR::STR::ModelType& mt) const;
    //! @}

    //! @name Monolithic update routines
    //! @{
    //! Update configuration after time step
    void UpdateStepState(const double& timefac_n);

    //! Update everything on element level after time step and after output
    void UpdateStepElement();

    //! calculation of stresses and strains
    void DetermineStressStrain();

    //! calculation of engery
    void DetermineEnergy();

    //! Write the current step state
    void OutputStepState(IO::DiscretizationWriter& iowriter) const;

    //! Write the current step state during runtime
    void RuntimeOutputStepState() const;

    //! Do things after writing output
    void PostOutput();

    //! reset all variables which are directly related to the current new step n+1.
    void ResetStepState();

    /*! \brief Recover the current state
     *
     * Necessary for condensed systems, e.g. EAS, dual mortar, etc.*/
    void RecoverState(
        const Epetra_Vector& xold,
        const Epetra_Vector& dir,
        const double& step,
        const Epetra_Vector& xnew,
        const bool isdefaultstep) const;

    /*! \brief Executed before the solution vector is going to be updated
     *
     *  \author hiermeier \date 03/17 */
    void RunPreComputeX(
        const Epetra_Vector& xold,
        Epetra_Vector& dir_mutable,
        const double& step,
        const NOX::NLN::Group& curr_grp,
        const bool isdefaultstep) const;

    /*! \brief Executed at the end of the NOX::Solver::Step() (f.k.a. Iterate()) method
     *
     *  \author hiermeier \date 03/17 */
    void RunPostIterate(
        const NOX::Solver::Generic& solver,
        const double step,
        const bool isdefaultstep ) const;

    //! @}

  protected:
    //! Returns the init flag.
    inline const bool& IsInit() const { return isinit_; };

    //! Returns the setup flag.
    inline const bool& IsSetup() const { return issetup_; };

    //! Check the init and setup state.
    void CheckInitSetup() const;

    //! Check the init state
    void CheckInit() const;

  private:
    Teuchos::RCP<STR::ModelEvaluator::Vector> TransformToVector(
        const STR::ModelEvaluator::Map& model_map ) const;

    /** \brief Assembly of all force contributions
     *
     *  \note There is also a PUBLIC alternative.
     *
     *  \param ok     (in/out) : flag indicating whether anything went wrong or not
     *  \param me_vec     (in) : vector of all models which are going to be assembled
     *  \param timefac_np (in) : time integration factor for the current contribution \f$t_{n+1}\f$
     *  \param f         (out) : force vector which is going to be assembled.
     *
     *  \author hiermeier \date 03/17 */
    void AssembleForce(
        bool& ok,
        const Vector& me_vec,
        const double timefac_np,
        Epetra_Vector& f ) const;

    /** \brief Assembly of all jacobian contributions
     *
     *  \param ok     (in/out) : flag indicating whether anything went wrong or not
     *  \param me_vec     (in) : vector of all models which are going to be assembled
     *  \param timefac_np (in) : time integration factor for the current contribution \f$t_{n+1}\f$
     *  \param f         (out) : jacobian matrix object which is going to be assembled.
     *
     *  \author hiermeier \date 03/17 */
    void AssembleJacobian(
        bool& ok,
        const Vector& me_vec,
        const double timefac_np,
        LINALG::SparseOperator& jac ) const;

    void EvaluateForce( bool& ok, const Vector& me_vec ) const;

    void EvaluateStiff( bool& ok, const Vector& me_vec ) const;

    void EvaluateForceStiff( bool& ok, const Vector& me_vec ) const;

    void PostEvaluate( bool ok, const Vector& me_vec ) const;

  private:
    bool isinit_;

    bool issetup_;

    Teuchos::RCP<STR::ModelEvaluator::Map> me_map_ptr_;

    Teuchos::RCP<STR::ModelEvaluator::Vector> me_vec_ptr_;

    Teuchos::RCP<STR::MODELEVALUATOR::Data> eval_data_ptr_;

    Teuchos::RCP<STR::TIMINT::BaseDataSDyn> sdyn_ptr_;

    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

    Teuchos::RCP<STR::TIMINT::BaseDataIO> gio_ptr_;

    Teuchos::RCP<STR::Integrator> int_ptr_;

    Teuchos::RCP<const STR::TIMINT::Base> timint_ptr_;
  };  // class ModelEvaluator
} // namespace STR



#endif /* SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_H_ */
