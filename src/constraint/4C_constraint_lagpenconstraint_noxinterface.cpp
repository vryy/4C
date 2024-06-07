/*---------------------------------------------------------------------*/
/*! \file

\brief Concrete mplementation of all the %NOX::Nln::CONSTRAINT::Interface::Required
       (pure) virtual routines.

\level 3


\date July 29, 2016

*/
/*---------------------------------------------------------------------*/

#include "4C_constraint_lagpenconstraint_noxinterface.hpp"

#include "4C_discretization_condition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_solver_nonlin_nox_aux.hpp"

#include <Epetra_Vector.h>
#include <NOX_Epetra_Vector.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LAGPENCONSTRAINT::NoxInterface::NoxInterface()
    : isinit_(false), issetup_(false), gstate_ptr_(Teuchos::null)
{
  // should stay empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LAGPENCONSTRAINT::NoxInterfacePrec::NoxInterfacePrec()
    : isinit_(false), issetup_(false), gstate_ptr_(Teuchos::null)
{
  // should stay empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LAGPENCONSTRAINT::NoxInterface::Init(
    const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& gstate_ptr)
{
  issetup_ = false;

  gstate_ptr_ = gstate_ptr;

  // set flag at the end
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LAGPENCONSTRAINT::NoxInterfacePrec::Init(
    const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& gstate_ptr)
{
  issetup_ = false;

  gstate_ptr_ = gstate_ptr;

  // set flag at the end
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LAGPENCONSTRAINT::NoxInterface::Setup()
{
  check_init();

  // set flag at the end
  issetup_ = true;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LAGPENCONSTRAINT::NoxInterfacePrec::Setup()
{
  check_init();

  // set flag at the end
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double LAGPENCONSTRAINT::NoxInterface::get_constraint_rhs_norms(const Epetra_Vector& F,
    NOX::Nln::StatusTest::QuantityType chQ, ::NOX::Abstract::Vector::NormType type,
    bool isScaled) const
{
  if (chQ != NOX::Nln::StatusTest::quantity_lag_pen_constraint) return -1.0;

  Teuchos::RCP<Epetra_Vector> constrRhs =
      gstate_ptr_->extract_model_entries(Inpar::STR::model_lag_pen_constraint, F);

  // no constraint contributions present
  if (constrRhs.is_null()) return 0.0;

  Teuchos::RCP<const ::NOX::Epetra::Vector> constrRhs_nox =
      Teuchos::rcp(new ::NOX::Epetra::Vector(constrRhs, ::NOX::Epetra::Vector::CreateView));

  double constrNorm = -1.0;
  constrNorm = constrRhs_nox->norm(type);
  if (isScaled) constrNorm /= static_cast<double>(constrRhs_nox->length());

  return constrNorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double LAGPENCONSTRAINT::NoxInterface::get_lagrange_multiplier_update_rms(const Epetra_Vector& xNew,
    const Epetra_Vector& xOld, double aTol, double rTol,
    NOX::Nln::StatusTest::QuantityType checkQuantity, bool disable_implicit_weighting) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_lag_pen_constraint) return -1.0;

  double rms = -1.0;

  // export the constraint solution
  Teuchos::RCP<Epetra_Vector> lagincr_ptr =
      gstate_ptr_->extract_model_entries(Inpar::STR::model_lag_pen_constraint, xOld);
  Teuchos::RCP<const Epetra_Vector> lagnew_ptr =
      gstate_ptr_->extract_model_entries(Inpar::STR::model_lag_pen_constraint, xNew);

  lagincr_ptr->Update(1.0, *lagnew_ptr, -1.0);
  Teuchos::RCP<const ::NOX::Epetra::Vector> lagincr_nox_ptr =
      Teuchos::rcp(new ::NOX::Epetra::Vector(lagincr_ptr, ::NOX::Epetra::Vector::CreateView));

  rms = NOX::Nln::Aux::RootMeanSquareNorm(
      aTol, rTol, lagnew_ptr, lagincr_ptr, disable_implicit_weighting);

  return rms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double LAGPENCONSTRAINT::NoxInterface::get_lagrange_multiplier_update_norms(
    const Epetra_Vector& xNew, const Epetra_Vector& xOld,
    NOX::Nln::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
    bool isScaled) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_lag_pen_constraint) return -1.0;

  // export the constraint solution
  Teuchos::RCP<Epetra_Vector> lagincr_ptr =
      gstate_ptr_->extract_model_entries(Inpar::STR::model_lag_pen_constraint, xOld);
  Teuchos::RCP<const Epetra_Vector> lagnew_ptr =
      gstate_ptr_->extract_model_entries(Inpar::STR::model_lag_pen_constraint, xNew);

  lagincr_ptr->Update(1.0, *lagnew_ptr, -1.0);
  Teuchos::RCP<const ::NOX::Epetra::Vector> lagincr_nox_ptr =
      Teuchos::rcp(new ::NOX::Epetra::Vector(lagincr_ptr, ::NOX::Epetra::Vector::CreateView));

  double updatenorm = -1.0;

  updatenorm = lagincr_nox_ptr->norm(type);
  // do scaling if desired
  if (isScaled) updatenorm /= static_cast<double>(lagincr_nox_ptr->length());

  return updatenorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double LAGPENCONSTRAINT::NoxInterface::get_previous_lagrange_multiplier_norms(
    const Epetra_Vector& xOld, NOX::Nln::StatusTest::QuantityType checkQuantity,
    ::NOX::Abstract::Vector::NormType type, bool isScaled) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_lag_pen_constraint) return -1.0;

  // export the constraint solution
  Teuchos::RCP<Epetra_Vector> lagold_ptr =
      gstate_ptr_->extract_model_entries(Inpar::STR::model_lag_pen_constraint, xOld);

  Teuchos::RCP<const ::NOX::Epetra::Vector> lagold_nox_ptr =
      Teuchos::rcp(new ::NOX::Epetra::Vector(lagold_ptr, ::NOX::Epetra::Vector::CreateView));

  double lagoldnorm = -1.0;

  lagoldnorm = lagold_nox_ptr->norm(type);
  // do scaling if desired
  if (isScaled) lagoldnorm /= static_cast<double>(lagold_nox_ptr->length());

  return lagoldnorm;
}



/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool LAGPENCONSTRAINT::NoxInterfacePrec::IsSaddlePointSystem() const
{
  Teuchos::RCP<const Discret::Discretization> dis = gstate_ptr_->get_discret();

  // ---------------------------------------------------------------------------
  // check type of constraint conditions (Lagrange multiplier vs. penalty)
  // ---------------------------------------------------------------------------
  bool have_lag_constraint = false;
  std::vector<Core::Conditions::Condition*> lagcond_volconstr3d(0);
  std::vector<Core::Conditions::Condition*> lagcond_areaconstr3d(0);
  std::vector<Core::Conditions::Condition*> lagcond_areaconstr2d(0);
  std::vector<Core::Conditions::Condition*> lagcond_mpconline2d(0);
  std::vector<Core::Conditions::Condition*> lagcond_mpconplane3d(0);
  std::vector<Core::Conditions::Condition*> lagcond_mpcnormcomp3d(0);
  dis->GetCondition("VolumeConstraint_3D", lagcond_volconstr3d);
  dis->GetCondition("AreaConstraint_3D", lagcond_areaconstr3d);
  dis->GetCondition("AreaConstraint_2D", lagcond_areaconstr2d);
  dis->GetCondition("MPC_NodeOnLine_2D", lagcond_mpconline2d);
  dis->GetCondition("MPC_NodeOnPlane_3D", lagcond_mpconplane3d);
  dis->GetCondition("MPC_NormalComponent_3D", lagcond_mpcnormcomp3d);
  if (lagcond_volconstr3d.size() or lagcond_areaconstr3d.size() or lagcond_areaconstr2d.size() or
      lagcond_mpconline2d.size() or lagcond_mpconplane3d.size() or lagcond_mpcnormcomp3d.size())
    have_lag_constraint = true;

  return have_lag_constraint;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool LAGPENCONSTRAINT::NoxInterfacePrec::IsCondensedSystem() const
{
  //  std::cout << "IsCondensedSystem" << std::endl;
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LAGPENCONSTRAINT::NoxInterfacePrec::fill_maps_for_preconditioner(
    std::vector<Teuchos::RCP<Epetra_Map>>& maps) const
{
  //  std::cout << "fill_maps_for_preconditioner" << std::endl;
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool LAGPENCONSTRAINT::NoxInterfacePrec::computePreconditioner(
    const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  //  std::cout << "computePreconditioner" << std::endl;
  check_init_setup();
  // currently not supported
  // ToDo add the scaled thickness conditioning (STC) approach here
  return false;
}

FOUR_C_NAMESPACE_CLOSE
