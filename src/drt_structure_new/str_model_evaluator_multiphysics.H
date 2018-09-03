/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator_multiphysics.H

\brief Base class for modelevaluators in partitioned algorithms.

\maintainer Andreas Rauch

\date Nov 28, 2016

\level 3

*/
/*-----------------------------------------------------------*/


#ifndef STR_MODEL_EVALUATOR_MULTIPHYSICS_H_
#define STR_MODEL_EVALUATOR_MULTIPHYSICS_H_


#include "str_model_evaluator_generic.H"


// forward declaration
class Epetra_Vector;
class Epetra_Map;
namespace LINALG
{
  class SparseOperator;
}  // namespace LINALG
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace IO
{
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace IO
namespace DRT
{
  class Discretization;
}  // namespace DRT
namespace STR
{
  class Integrator;
  namespace TIMINT
  {
    class BaseDataGlobalState;
    class BaseDataIO;
    class Base;
  }  // namespace TIMINT
  namespace MODELEVALUATOR
  {
    class Data;


    //! supported multiphysic problems
    enum MultiphysicType
    {
      mt_none = 0,  //!< none specific default value
      mt_fsi = 1,   //!< multiphysics type fluid-structure-interaction
      mt_ssi = 2    //!< multiphysics type structure-scalar-interaction
    };              // MultiphysicType


    /*! \brief This is the base class for all multiphysics models.
     *
     *  This class summarizes the functionality which all model multiphysics model
     *  evaluators share.
     *
     *  \date 11/16
     *  \author rauch */
    class Multiphysics : public Generic
    {
     public:
      //! constructor
      Multiphysics();

      //! destructor
      virtual ~Multiphysics(){};

      //! initialize the class variables
      virtual void Init(const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
          const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
          const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
          const Teuchos::RCP<STR::Integrator>& int_ptr,
          const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr, const int& dof_offset);

      //! setup class variables
      virtual void Setup();

      //! set the active model type wrapped in this class.
      //! only active model type is evaluated.
      //! e.g. mt_fsi in case fluid-structure interaction is to be evaluated
      void SetActiveModelType(enum STR::MODELEVALUATOR::MultiphysicType mtype)
      {
        active_mt_ = mtype;
      };

      void CheckActiveModelType() const
      {
        if (active_mt_ == mt_none) dserror("No active model evaluator set for Multiphysics");
      };

      //! @name Functions which are derived from the base generic class
      //! @{
      //! [derived]
      virtual INPAR::STR::ModelType Type() const { return INPAR::STR::model_partitioned_coupling; }

      //! reset class variables (without jacobian) [derived]
      virtual void Reset(const Epetra_Vector& x);

      //! [derived]
      virtual bool EvaluateForce();

      //! [derived]
      virtual bool EvaluateStiff();

      //! [derived] not needed in partitioned scheme
      virtual bool EvaluateForceStiff();

      //! derived
      virtual void PostEvaluate() { return; };

      //! derived
      virtual bool AssembleForce(Epetra_Vector& f, const double& timefac_np) const;

      //! Assemble the jacobian at \f$t_{n+1}\f$ not needed in partitioned scheme
      virtual bool AssembleJacobian(LINALG::SparseOperator& jac, const double& timefac_np) const;

      //! [derived]
      virtual void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const
      {
        return;
      };

      //! [derived]
      virtual void ReadRestart(IO::DiscretizationReader& ioreader) { return; };

      //! [derived]
      virtual void Predict(const INPAR::STR::PredEnum& pred_type) { return; };

      //! derived
      virtual void RunPreComputeX(
          const Epetra_Vector& xold, Epetra_Vector& dir_mutable, const NOX::NLN::Group& curr_grp)
      {
        return;
      };

      //! recover condensed Lagrange multipliers
      virtual void RunPostComputeX(
          const Epetra_Vector& xold, const Epetra_Vector& dir, const Epetra_Vector& xnew)
      {
        return;
      };

      //! derived
      virtual void RunPostIterate(const NOX::Solver::Generic& solver) { return; };

      //! [derived]
      virtual void UpdateStepState(const double& timefac_n);

      //! [derived]
      virtual void UpdateStepElement() { return; };

      //! [derived]
      virtual void DetermineStressStrain() { return; };

      //! [derived]
      virtual void DetermineEnergy() { return; };

      //! [derived]
      virtual void DetermineOptionalQuantity() override { return; };

      //! [derived]
      virtual void OutputStepState(IO::DiscretizationWriter& iowriter) const { return; };

      //! derived
      virtual void ResetStepState() { return; };

      //! [derived]
      virtual void PostOutput() { return; };

      //! @name Accessors to model specific things
      //! @{
      //! Returns a pointer to the model specific dof row map
      virtual Teuchos::RCP<const Epetra_Map> GetBlockDofRowMapPtr() const { return Teuchos::null; };

      //! Returns a pointer to the current model solution vector (usually the Lagrange multiplier
      //! vector)
      virtual Teuchos::RCP<const Epetra_Vector> GetCurrentSolutionPtr() const
      {
        return Teuchos::null;
      };

      //! Returns a pointer to the model solution vector of the last time step (usually the Lagrange
      //! multiplier vector)
      virtual Teuchos::RCP<const Epetra_Vector> GetLastTimeStepSolutionPtr() const
      {
        return Teuchos::null;
      };

      //! @}

     protected:
      //! map containing the model evaluators of the sub modules
      std::map<enum STR::MODELEVALUATOR::MultiphysicType,
          Teuchos::RCP<STR::MODELEVALUATOR::Generic>>
          me_map_;

      //! currently active model evaluator type
      STR::MODELEVALUATOR::MultiphysicType active_mt_;

      //! return reference to map containing the model evaluators
      std::map<enum STR::MODELEVALUATOR::MultiphysicType,
          Teuchos::RCP<STR::MODELEVALUATOR::Generic>>&
      GetModelEvalutaorMap()
      {
        return me_map_;
      };

     public:
      //! return RCP to model evaluator of specific MultiphysicType
      Teuchos::RCP<STR::MODELEVALUATOR::Generic> GetModelEvaluatorFromMap(
          enum STR::MODELEVALUATOR::MultiphysicType mtype) const
      {
        return me_map_.at(mtype);
      }


    };  // class Multiphysics

  }  // namespace MODELEVALUATOR
}  // namespace STR

#endif /* STR_MODEL_EVALUATOR_MULTIPHYSICS_H_ */
