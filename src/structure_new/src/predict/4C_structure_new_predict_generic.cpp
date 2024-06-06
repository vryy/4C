/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all predictors.


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_predict_generic.hpp"

#include "4C_io_pstream.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_structure_new_dbc.hpp"
#include "4C_structure_new_impl_generic.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_timint_base.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Predict::Generic::Generic()
    : isinit_(false),
      issetup_(false),
      type_(Inpar::STR::pred_vague),
      implint_ptr_(Teuchos::null),
      dbc_ptr_(Teuchos::null),
      noxparams_ptr_(Teuchos::null)
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Predict::Generic::Init(const enum Inpar::STR::PredEnum& type,
    const Teuchos::RCP<STR::IMPLICIT::Generic>& implint_ptr, const Teuchos::RCP<STR::Dbc>& dbc_ptr,
    const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TimeInt::BaseDataIO>& iodata_ptr,
    const Teuchos::RCP<Teuchos::ParameterList>& noxparams_ptr)
{
  issetup_ = false;

  // initialize the predictor type
  type_ = type;
  implint_ptr_ = implint_ptr;
  dbc_ptr_ = dbc_ptr;
  gstate_ptr_ = gstate_ptr;
  iodata_ptr_ = iodata_ptr;
  noxparams_ptr_ = noxparams_ptr;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Predict::Generic::PrePredict(::NOX::Abstract::Group& grp)
{
  check_init_setup();
  Print();
  dbc_ptr_->UpdateLocSysManager();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Predict::Generic::Predict(::NOX::Abstract::Group& grp)
{
  check_init_setup();
  bool& ispredict = gstate_ptr_->IsPredict();
  ispredict = true;

  // pre-process the prediction step
  PrePredict(grp);

  // compute the actual prediction step
  Compute(grp);

  // post-process the prediction step
  PostPredict(grp);

  ispredict = false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Predict::Generic::PostPredict(::NOX::Abstract::Group& grp)
{
  check_init_setup();

  dbc().apply_dirichlet_bc(global_state().GetTimeNp(), global_state().GetDisNp(),
      global_state().GetVelNp(), global_state().GetAccNp(), false);

  // Create the new solution vector
  Teuchos::RCP<::NOX::Epetra::Vector> x_vec = global_state().CreateGlobalVector(
      TimeInt::BaseDataGlobalState::VecInitType::init_current_state, impl_int().ModelEvalPtr());
  // resets all isValid flags
  grp.setX(*x_vec);

  NOX::Nln::Group* nlngrp_ptr = dynamic_cast<NOX::Nln::Group*>(&grp);
  FOUR_C_ASSERT(nlngrp_ptr != nullptr, "Group cast failed!");
  // evaluate the right hand side and the jacobian
  implint_ptr_->SetIsPredictorState(true);
  nlngrp_ptr->computeFandJacobian();
  implint_ptr_->SetIsPredictorState(false);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::string STR::Predict::Generic::Name() const
{
  check_init();
  return Inpar::STR::PredEnumString(type_);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Predict::Generic::check_init() const { FOUR_C_ASSERT(is_init(), "Call Init() first!"); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Predict::Generic::check_init_setup() const
{
  FOUR_C_ASSERT(is_init() and is_setup(), "Call Init() and Setup() first!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::IMPLICIT::Generic>& STR::Predict::Generic::impl_int_ptr()
{
  check_init();
  return implint_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::Generic& STR::Predict::Generic::impl_int()
{
  check_init();
  return *implint_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::Dbc>& STR::Predict::Generic::dbc_ptr()
{
  check_init();
  return dbc_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::Dbc& STR::Predict::Generic::dbc()
{
  check_init();
  return *dbc_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& STR::Predict::Generic::global_state_ptr()
{
  check_init();
  return gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TimeInt::BaseDataGlobalState& STR::Predict::Generic::global_state()
{
  check_init();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<STR::TimeInt::BaseDataIO>& STR::Predict::Generic::io_data_ptr()
{
  check_init();
  return iodata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TimeInt::BaseDataIO& STR::Predict::Generic::io_data()
{
  check_init();
  return *iodata_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const STR::TimeInt::BaseDataGlobalState& STR::Predict::Generic::GlobalState() const
{
  check_init();
  return *gstate_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Teuchos::ParameterList>& STR::Predict::Generic::nox_params_ptr()
{
  check_init();
  return noxparams_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::ParameterList& STR::Predict::Generic::nox_params()
{
  check_init();
  return *noxparams_ptr_;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::Predict::Generic::Print() const
{
  check_init_setup();
  if (gstate_ptr_->GetMyRank() == 0 and iodata_ptr_->get_print2_screen_every_n_step() and
      gstate_ptr_->GetStepN() % iodata_ptr_->get_print2_screen_every_n_step() == 0)
  {
    Core::IO::cout << "=== Structural predictor: " << Name().c_str() << " ===" << Core::IO::endl;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::Predict::Generic::pre_apply_force_external(Epetra_Vector& fextnp) const
{
  // do nothing
  return false;
}

FOUR_C_NAMESPACE_CLOSE
