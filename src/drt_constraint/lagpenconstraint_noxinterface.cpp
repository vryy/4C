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
      issetup_(false),
      gstate_ptr_(Teuchos::null)
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
void LAGPENCONSTRAINT::NoxInterfacePrec::Init(
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr)
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
    NOX::NLN::StatusTest::QuantityType chQ,
    NOX::Abstract::Vector::NormType type,
    bool isScaled) const
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
    double aTol,
    double rTol,
    NOX::NLN::StatusTest::QuantityType chQ,
    bool disable_implicit_weighting) const
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
    NOX::NLN::StatusTest::QuantityType chQ,
    NOX::Abstract::Vector::NormType type,
    bool isScaled) const
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
    NOX::NLN::StatusTest::QuantityType chQ,
    NOX::Abstract::Vector::NormType type,
    bool isScaled) const
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

  Teuchos::RCP<const DRT::DiscretizationInterface> dis = gstate_ptr_->GetDiscret();

  // ---------------------------------------------------------------------------
  // check type of constraint conditions (Lagrange multiplier vs. penalty)
  // ---------------------------------------------------------------------------
  bool have_lag_constraint = false;
  std::vector<DRT::Condition*> lagcond_volconstr3d(0);
  std::vector<DRT::Condition*> lagcond_areaconstr3d(0);
  std::vector<DRT::Condition*> lagcond_areaconstr2d(0);
  std::vector<DRT::Condition*> lagcond_mpconline2d(0);
  std::vector<DRT::Condition*> lagcond_mpconplane3d(0);
  std::vector<DRT::Condition*> lagcond_mpcnormcomp3d(0);
  dis->GetCondition("VolumeConstraint_3D",lagcond_volconstr3d);
  dis->GetCondition("AreaConstraint_3D",lagcond_areaconstr3d);
  dis->GetCondition("AreaConstraint_2D",lagcond_areaconstr2d);
  dis->GetCondition("MPC_NodeOnLine_2D",lagcond_mpconline2d);
  dis->GetCondition("MPC_NodeOnPlane_3D",lagcond_mpconplane3d);
  dis->GetCondition("MPC_NormalComponent_3D",lagcond_mpcnormcomp3d);
  if (
         lagcond_volconstr3d.size()  or
         lagcond_areaconstr3d.size() or
         lagcond_areaconstr2d.size() or
         lagcond_mpconline2d.size()  or
         lagcond_mpconplane3d.size() or
         lagcond_mpcnormcomp3d.size()
      )
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
void LAGPENCONSTRAINT::NoxInterfacePrec::FillMapsForPreconditioner(std::vector<Teuchos::RCP<Epetra_Map> >& maps) const
{
//  std::cout << "FillMapsForPreconditioner" << std::endl;
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool LAGPENCONSTRAINT::NoxInterfacePrec::computePreconditioner(
    const Epetra_Vector& x,
    Epetra_Operator& M,
    Teuchos::ParameterList* precParams)
{
//  std::cout << "computePreconditioner" << std::endl;
  CheckInitSetup();
  // currently not supported
  // ToDo add the scaled thickness conditioning (STC) approach here
  return false;
}
