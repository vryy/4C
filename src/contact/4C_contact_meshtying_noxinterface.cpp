#include "4C_contact_meshtying_noxinterface.hpp"

#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_solver_nonlin_nox_aux.hpp"

#include <NOX_Epetra_Vector.H>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::MtNoxInterface::MtNoxInterface()
    : isinit_(false), issetup_(false), gstate_ptr_(Teuchos::null)
{
  // should stay empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::MtNoxInterface::init(
    const Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState>& gstate_ptr)
{
  issetup_ = false;

  gstate_ptr_ = gstate_ptr;

  // set flag at the end
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::MtNoxInterface::setup()
{
  check_init();

  // set flag at the end
  issetup_ = true;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::MtNoxInterface::get_constraint_rhs_norms(const Core::LinAlg::Vector<double>& F,
    NOX::Nln::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
    bool isScaled) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_meshtying) return -1.0;

  Teuchos::RCP<Core::LinAlg::Vector<double>> constrRhs =
      gstate_ptr_->extract_model_entries(Inpar::Solid::model_meshtying, F);

  // no constraint contributions present
  if (constrRhs.is_null()) return 0.0;

  const ::NOX::Epetra::Vector constrRhs_nox(
      constrRhs->get_ptr_of_Epetra_Vector(), ::NOX::Epetra::Vector::CreateView);

  double constrNorm = -1.0;
  constrNorm = constrRhs_nox.norm(type);
  if (isScaled) constrNorm /= static_cast<double>(constrRhs_nox.length());

  return constrNorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::MtNoxInterface::get_lagrange_multiplier_update_rms(
    const Core::LinAlg::Vector<double>& xNew, const Core::LinAlg::Vector<double>& xOld, double aTol,
    double rTol, NOX::Nln::StatusTest::QuantityType checkQuantity,
    bool disable_implicit_weighting) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_meshtying) return -1.0;

  double rms = -1.0;

  // export the constraint solution
  Teuchos::RCP<Core::LinAlg::Vector<double>> lagincr_ptr =
      gstate_ptr_->extract_model_entries(Inpar::Solid::model_meshtying, xOld);
  Teuchos::RCP<const Core::LinAlg::Vector<double>> lagnew_ptr =
      gstate_ptr_->extract_model_entries(Inpar::Solid::model_meshtying, xNew);

  lagincr_ptr->Update(1.0, *lagnew_ptr, -1.0);
  const ::NOX::Epetra::Vector lagincr_nox_ptr(
      lagincr_ptr->get_ptr_of_Epetra_Vector(), ::NOX::Epetra::Vector::CreateView);

  rms = NOX::Nln::Aux::root_mean_square_norm(
      aTol, rTol, *lagnew_ptr, *lagincr_ptr, disable_implicit_weighting);

  return rms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::MtNoxInterface::get_lagrange_multiplier_update_norms(
    const Core::LinAlg::Vector<double>& xNew, const Core::LinAlg::Vector<double>& xOld,
    NOX::Nln::StatusTest::QuantityType checkQuantity, ::NOX::Abstract::Vector::NormType type,
    bool isScaled) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_meshtying) return -1.0;

  // export the constraint solution
  Teuchos::RCP<Core::LinAlg::Vector<double>> lagincr_ptr =
      gstate_ptr_->extract_model_entries(Inpar::Solid::model_meshtying, xOld);
  Teuchos::RCP<const Core::LinAlg::Vector<double>> lagnew_ptr =
      gstate_ptr_->extract_model_entries(Inpar::Solid::model_meshtying, xNew);

  lagincr_ptr->Update(1.0, *lagnew_ptr, -1.0);
  const ::NOX::Epetra::Vector lagincr_nox_ptr(
      lagincr_ptr->get_ptr_of_Epetra_Vector(), ::NOX::Epetra::Vector::CreateView);

  double updatenorm = -1.0;

  updatenorm = lagincr_nox_ptr.norm(type);
  // do scaling if desired
  if (isScaled) updatenorm /= static_cast<double>(lagincr_nox_ptr.length());

  return updatenorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::MtNoxInterface::get_previous_lagrange_multiplier_norms(
    const Core::LinAlg::Vector<double>& xOld, NOX::Nln::StatusTest::QuantityType checkQuantity,
    ::NOX::Abstract::Vector::NormType type, bool isScaled) const
{
  if (checkQuantity != NOX::Nln::StatusTest::quantity_meshtying) return -1.0;

  // export the constraint solution
  Teuchos::RCP<Core::LinAlg::Vector<double>> lagold_ptr =
      gstate_ptr_->extract_model_entries(Inpar::Solid::model_meshtying, xOld);

  const ::NOX::Epetra::Vector lagold_nox_ptr(
      lagold_ptr->get_ptr_of_Epetra_Vector(), ::NOX::Epetra::Vector::CreateView);

  double lagoldnorm = -1.0;

  lagoldnorm = lagold_nox_ptr.norm(type);
  // do scaling if desired
  if (isScaled) lagoldnorm /= static_cast<double>(lagold_nox_ptr.length());

  return lagoldnorm;
}

FOUR_C_NAMESPACE_CLOSE
