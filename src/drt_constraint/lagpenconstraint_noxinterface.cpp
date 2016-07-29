/*---------------------------------------------------------------------*/
/*!
\file lagpenconstraint_noxinterface.cpp

\brief Concrete mplementation of all the %NOX::NLN::CONSTRAINT::Interface::Required
       (pure) virtual routines.

\level 3

\maintainer Marc Hirschvogel

\date July 29, 2016

*/
/*---------------------------------------------------------------------*/

#include "lagpenconstraint_noxinterface.H"

#include "../solver_nonlin_nox/nox_nln_aux.H"

#include "../linalg/linalg_utils.H"

#include <NOX_Epetra_Vector.H>
#include <Epetra_Vector.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LAGPENCONSTRAINT::NoxInterface::NoxInterface()
    : isinit_(false),
      issetup_(false),
      gstate_ptr_(Teuchos::null)
{
  // should stay empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LAGPENCONSTRAINT::NoxInterfacePrec::NoxInterfacePrec()
    : isinit_(false),
      issetup_(false)
{
  // should stay empty
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LAGPENCONSTRAINT::NoxInterface::Init(
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr)
{
  issetup_ = false;

  gstate_ptr_ = gstate_ptr;

  // set flag at the end
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LAGPENCONSTRAINT::NoxInterfacePrec::Init()
{
  issetup_ = false;

  // set flag at the end
  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LAGPENCONSTRAINT::NoxInterface::Setup()
{
  CheckInit();

  // set flag at the end
  issetup_ = true;
}
/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LAGPENCONSTRAINT::NoxInterfacePrec::Setup()
{
  CheckInit();

  // set flag at the end
  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double LAGPENCONSTRAINT::NoxInterface::GetConstraintRHSNorms(
    const Epetra_Vector& F,
    const NOX::NLN::StatusTest::QuantityType& chQ,
    const NOX::Abstract::Vector::NormType& type,
    const bool& isScaled) const
{

  if (chQ != NOX::NLN::StatusTest::quantity_lag_pen_constraint)
    return -1.0;

  Teuchos::RCP<Epetra_Vector> constrRhs =
      gstate_ptr_->ExtractModelEntries(INPAR::STR::model_lag_pen_constraint,F);

  // no constraint contributions present
  if (constrRhs.is_null())
    return 0.0;

  Teuchos::RCP<const NOX::Epetra::Vector> constrRhs_nox =
      Teuchos::rcp(new NOX::Epetra::Vector(constrRhs,NOX::Epetra::Vector::CreateView));

  double constrNorm = -1.0;

  constrNorm = constrRhs_nox->norm(type);
  if (isScaled)
    constrNorm /= static_cast<double>(constrRhs_nox->length());

  return constrNorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double LAGPENCONSTRAINT::NoxInterface::GetLagrangeMultiplierUpdateRMS(
    const Epetra_Vector& xNew,
    const Epetra_Vector& xOld,
    const double& aTol,
    const double& rTol,
    const NOX::NLN::StatusTest::QuantityType& chQ,
    const bool& disable_implicit_weighting) const
{

  if (chQ != NOX::NLN::StatusTest::quantity_lag_pen_constraint)
    return -1.0;

  double rms = -1.0;

  // export the constraint solution
  Teuchos::RCP<Epetra_Vector> lagincr_ptr =
      gstate_ptr_->ExtractModelEntries(INPAR::STR::model_lag_pen_constraint,xOld);
  Teuchos::RCP<const Epetra_Vector> lagnew_ptr =
      gstate_ptr_->ExtractModelEntries(INPAR::STR::model_lag_pen_constraint,xNew);

  lagincr_ptr->Update(1.0,*lagnew_ptr,-1.0);
  Teuchos::RCP<const NOX::Epetra::Vector> lagincr_nox_ptr =
      Teuchos::rcp(new NOX::Epetra::Vector(lagincr_ptr,NOX::Epetra::Vector::CreateView));

  rms = NOX::NLN::AUX::RootMeanSquareNorm(aTol,rTol,
      lagnew_ptr,lagincr_ptr,disable_implicit_weighting);

  return rms;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double LAGPENCONSTRAINT::NoxInterface::GetLagrangeMultiplierUpdateNorms(
    const Epetra_Vector& xNew,
    const Epetra_Vector& xOld,
    const NOX::NLN::StatusTest::QuantityType& chQ,
    const NOX::Abstract::Vector::NormType& type,
    const bool& isScaled) const
{

  if (chQ != NOX::NLN::StatusTest::quantity_lag_pen_constraint)
    return -1.0;

  // export the constraint solution
  Teuchos::RCP<Epetra_Vector> lagincr_ptr =
      gstate_ptr_->ExtractModelEntries(INPAR::STR::model_lag_pen_constraint,xOld);
  Teuchos::RCP<const Epetra_Vector> lagnew_ptr =
      gstate_ptr_->ExtractModelEntries(INPAR::STR::model_lag_pen_constraint,xNew);

  lagincr_ptr->Update(1.0,*lagnew_ptr,-1.0);
  Teuchos::RCP<const NOX::Epetra::Vector> lagincr_nox_ptr =
      Teuchos::rcp(new NOX::Epetra::Vector(lagincr_ptr,NOX::Epetra::Vector::CreateView));

  double updatenorm = -1.0;

  updatenorm = lagincr_nox_ptr->norm(type);
  // do scaling if desired
  if (isScaled)
    updatenorm /= static_cast<double>(lagincr_nox_ptr->length());

  return updatenorm;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double LAGPENCONSTRAINT::NoxInterface::GetPreviousLagrangeMultiplierNorms(
    const Epetra_Vector& xOld,
    const NOX::NLN::StatusTest::QuantityType& chQ,
    const NOX::Abstract::Vector::NormType& type,
    const bool& isScaled) const
{

  if (chQ != NOX::NLN::StatusTest::quantity_lag_pen_constraint)
    return -1.0;

  // export the constraint solution
  Teuchos::RCP<Epetra_Vector> lagold_ptr =
      gstate_ptr_->ExtractModelEntries(INPAR::STR::model_lag_pen_constraint,xOld);

  Teuchos::RCP<const NOX::Epetra::Vector> lagold_nox_ptr =
      Teuchos::rcp(new NOX::Epetra::Vector(lagold_ptr,NOX::Epetra::Vector::CreateView));

  double lagoldnorm = -1.0;

  lagoldnorm = lagold_nox_ptr->norm(type);
  // do scaling if desired
  if (isScaled)
    lagoldnorm /= static_cast<double>(lagold_nox_ptr->length());

  return lagoldnorm;
}




/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool LAGPENCONSTRAINT::NoxInterfacePrec::IsSaddlePointSystem() const
{
//  std::cout << "IsSaddlePointSystem" << std::endl;
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool LAGPENCONSTRAINT::NoxInterfacePrec::IsCondensedSystem() const
{
//  std::cout << "IsCondensedSystem" << std::endl;
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void LAGPENCONSTRAINT::NoxInterfacePrec::FillMapsForPreconditioner(std::vector<Teuchos::RCP<Epetra_Map> >& maps) const
{
//  std::cout << "FillMapsForPreconditioner" << std::endl;
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool LAGPENCONSTRAINT::NoxInterfacePrec::computePreconditioner(
    const Epetra_Vector& x,
    Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
//  std::cout << "computePreconditioner" << std::endl;
  CheckInitSetup();
  // currently not supported
  // ToDo add the scaled thickness conditioning (STC) approach here
  return false;
}
