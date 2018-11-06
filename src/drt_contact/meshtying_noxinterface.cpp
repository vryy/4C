/*---------------------------------------------------------------------*/
/*!
\file meshtying_noxinterface.cpp

\brief Concrete implementation of all the %NOX::NLN::CONSTRAINT::Interface::Required
       (pure) virtual routines.

\level 3

\maintainer Matthias Mayr

*/
/*---------------------------------------------------------------------*/

#include "../solver_nonlin_nox/nox_nln_aux.H"

#include "../linalg/linalg_utils.H"

#include <NOX_Epetra_Vector.H>
#include <Epetra_Vector.h>
#include "meshtying_noxinterface.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
CONTACT::MtNoxInterface::MtNoxInterface()
    : isinit_(false), issetup_(false), gstate_ptr_(Teuchos::null)
{
  // should stay empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::MtNoxInterface::Init(const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr)
{
  issetup_ = false;

  gstate_ptr_ = gstate_ptr;

  // set flag at the end
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void CONTACT::MtNoxInterface::Setup()
{
  CheckInit();

  // set flag at the end
  issetup_ = true;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::MtNoxInterface::GetConstraintRHSNorms(const Epetra_Vector& F,
    NOX::NLN::StatusTest::QuantityType chQ, NOX::Abstract::Vector::NormType type,
    bool isScaled) const
{
  if (chQ != NOX::NLN::StatusTest::quantity_meshtying) return -1.0;

  Teuchos::RCP<Epetra_Vector> constrRhs =
      gstate_ptr_->ExtractModelEntries(INPAR::STR::model_meshtying, F);

  // no constraint contributions present
  if (constrRhs.is_null()) return 0.0;

  Teuchos::RCP<const NOX::Epetra::Vector> constrRhs_nox =
      Teuchos::rcp(new NOX::Epetra::Vector(constrRhs, NOX::Epetra::Vector::CreateView));

  double constrNorm = -1.0;
  constrNorm = constrRhs_nox->norm(type);
  if (isScaled) constrNorm /= static_cast<double>(constrRhs_nox->length());

  return constrNorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::MtNoxInterface::GetLagrangeMultiplierUpdateRMS(const Epetra_Vector& xNew,
    const Epetra_Vector& xOld, double aTol, double rTol, NOX::NLN::StatusTest::QuantityType chQ,
    bool disable_implicit_weighting) const
{
  if (chQ != NOX::NLN::StatusTest::quantity_meshtying) return -1.0;

  double rms = -1.0;

  // export the constraint solution
  Teuchos::RCP<Epetra_Vector> lagincr_ptr =
      gstate_ptr_->ExtractModelEntries(INPAR::STR::model_meshtying, xOld);
  Teuchos::RCP<const Epetra_Vector> lagnew_ptr =
      gstate_ptr_->ExtractModelEntries(INPAR::STR::model_meshtying, xNew);

  lagincr_ptr->Update(1.0, *lagnew_ptr, -1.0);
  Teuchos::RCP<const NOX::Epetra::Vector> lagincr_nox_ptr =
      Teuchos::rcp(new NOX::Epetra::Vector(lagincr_ptr, NOX::Epetra::Vector::CreateView));

  rms = NOX::NLN::AUX::RootMeanSquareNorm(
      aTol, rTol, lagnew_ptr, lagincr_ptr, disable_implicit_weighting);

  return rms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::MtNoxInterface::GetLagrangeMultiplierUpdateNorms(const Epetra_Vector& xNew,
    const Epetra_Vector& xOld, NOX::NLN::StatusTest::QuantityType chQ,
    NOX::Abstract::Vector::NormType type, bool isScaled) const
{
  if (chQ != NOX::NLN::StatusTest::quantity_meshtying) return -1.0;

  // export the constraint solution
  Teuchos::RCP<Epetra_Vector> lagincr_ptr =
      gstate_ptr_->ExtractModelEntries(INPAR::STR::model_meshtying, xOld);
  Teuchos::RCP<const Epetra_Vector> lagnew_ptr =
      gstate_ptr_->ExtractModelEntries(INPAR::STR::model_meshtying, xNew);

  lagincr_ptr->Update(1.0, *lagnew_ptr, -1.0);
  Teuchos::RCP<const NOX::Epetra::Vector> lagincr_nox_ptr =
      Teuchos::rcp(new NOX::Epetra::Vector(lagincr_ptr, NOX::Epetra::Vector::CreateView));

  double updatenorm = -1.0;

  updatenorm = lagincr_nox_ptr->norm(type);
  // do scaling if desired
  if (isScaled) updatenorm /= static_cast<double>(lagincr_nox_ptr->length());

  return updatenorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double CONTACT::MtNoxInterface::GetPreviousLagrangeMultiplierNorms(const Epetra_Vector& xOld,
    NOX::NLN::StatusTest::QuantityType chQ, NOX::Abstract::Vector::NormType type,
    bool isScaled) const
{
  if (chQ != NOX::NLN::StatusTest::quantity_meshtying) return -1.0;

  // export the constraint solution
  Teuchos::RCP<Epetra_Vector> lagold_ptr =
      gstate_ptr_->ExtractModelEntries(INPAR::STR::model_meshtying, xOld);

  Teuchos::RCP<const NOX::Epetra::Vector> lagold_nox_ptr =
      Teuchos::rcp(new NOX::Epetra::Vector(lagold_ptr, NOX::Epetra::Vector::CreateView));

  double lagoldnorm = -1.0;

  lagoldnorm = lagold_nox_ptr->norm(type);
  // do scaling if desired
  if (isScaled) lagoldnorm /= static_cast<double>(lagold_nox_ptr->length());

  return lagoldnorm;
}
