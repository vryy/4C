/*-----------------------------------------------------------*/
/*!
\file str_impl_ost.cpp

\maintainer Philipp Farah

\date Dec 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_impl_ost.H"

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_sparseoperator.H"

#include <Epetra_Vector.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::OneStepTheta::OneStepTheta()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::Setup()
{
  CheckInit();
  // Call the Setup() of the abstract base class first.
  Generic::Setup();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::SetState(const Epetra_Vector& x)
{
  dserror("Not yet implemented! (see the Statics integration for an example)");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::ApplyForce(const Epetra_Vector& x,
        Epetra_Vector& f)
{
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::ApplyStiff(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::ApplyForceStiff(
    const Epetra_Vector& x,
    Epetra_Vector& f,
    LINALG::SparseOperator& jac)
{
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::OneStepTheta::CalcRefNormForce(
    const enum NOX::Abstract::Vector::NormType& type)
{
  dserror("Not yet implemented! (see the Statics integration for an example)");
  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::OneStepTheta::GetIntParam() const
{
  dserror("Set the time integration parameter as return value!");
  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::UpdateStepState()
{
  dserror("Not yet implemented! (see the Statics integration for an example)");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::UpdateStepElement()
{
  dserror("Not yet implemented! (see the Statics integration for an example)");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::OneStepTheta::PredictConstDisConsistVelAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  dserror("Not yet implemented! (see the Statics integration for an example)");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::PredictConstVelConsistAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  dserror("Not yet implemented! (see the Statics integration for an example)");
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::OneStepTheta::PredictConstAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  dserror("Not yet implemented! (see the Statics integration for an example)");
  return false;
}
