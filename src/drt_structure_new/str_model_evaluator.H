/*-----------------------------------------------------------*/
/*!
\file str_model_evaluator.H

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
namespace LINALG {
  class SparseOperator;
} // namespace LINALG
namespace STR {
namespace TIMINT {
  class Base;
  class BaseDataGlobalState;
} // namespace TIMINT
namespace MODELEVALUATOR {
  class Generic;
} // namespace MODELEVALUATOR
  //! Wrapper class for the STR::MODELEVALUATOR::Generic derived objects.
  class ModelEvaluator
  {
  public:
    typedef std::map<enum INPAR::STR::ModelType, Teuchos::RCP<STR::MODELEVALUATOR::Generic> > Map;
    //! constructor
    ModelEvaluator();

    //! destructor
    virtual ~ModelEvaluator() {};

    //! initialize
    void Init(const std::set<enum INPAR::STR::ModelType>& mt,
        const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
        const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr);

    //! setup
    void Setup();

    //! apply force
    bool ApplyForce(const Epetra_Vector& x, Epetra_Vector& f) const;

    //! apply stiffness
    bool ApplyStiff(const Epetra_Vector& x, LINALG::SparseOperator& jac) const;

    //! apply force and stiffness
    bool ApplyForceStiff(const Epetra_Vector& x, Epetra_Vector& f, LINALG::SparseOperator& jac) const;

    //! @name Accessors
    //! @{
    //! return global state (read-only)
    const STR::TIMINT::BaseDataGlobalState& GetGlobalState() const;

    //! return global state pointer (read and write access of the data)
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& GlobalStatePtr();

    //! return pointer to the underlying time integrator (read-only)
    const Teuchos::RCP<const STR::TIMINT::Base>& GetTimIntPtr() const;

    //! @}

    //! @name Monolithic update routines
    //! @{
    //! Update configuration after time step
    void UpdateStepState();

    //! Update everything on element level after time step and after output
    void UpdateStepElement();
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

  protected:
    bool isinit_;

    bool issetup_;

    Teuchos::RCP<const std::set<enum INPAR::STR::ModelType> > modeltypes_ptr_;

    Teuchos::RCP<STR::ModelEvaluator::Map> modelevaluators_ptr_;

    Teuchos::RCP<STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

    Teuchos::RCP<const STR::TIMINT::Base> timint_ptr_;
  };  // class ModelEvaluator
} // namespace STR



#endif /* SRC_DRT_STRUCTURE_NEW_STR_MODEL_EVALUATOR_H_ */
