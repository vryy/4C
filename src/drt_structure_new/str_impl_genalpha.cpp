/*-----------------------------------------------------------*/
/*!
\file str_impl_genalpha.cpp

\maintainer Philipp Farah

\date Dec 14, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "str_impl_genalpha.H"

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_sparseoperator.H"

#include <Epetra_Vector.h>


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::IMPLICIT::GenAlpha::GenAlpha()
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::Setup()
{
  CheckInit();

  issetup_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::SetState(const Epetra_Vector& x)
{
  CheckInitSetup();
  dserror("Not yet implemented! (see the Statics integration for an example)");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlpha::ApplyForce(const Epetra_Vector& x,
        Epetra_Vector& f)
{
  CheckInitSetup();
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlpha::ApplyStiff(
    const Epetra_Vector& x,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlpha::ApplyForceStiff(
    const Epetra_Vector& x,
    Epetra_Vector& f,
    LINALG::SparseOperator& jac)
{
  CheckInitSetup();
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::GenAlpha::CalcRefNormForce(
    const enum NOX::Abstract::Vector::NormType& type)
{
  CheckInitSetup();
  dserror("Not yet implemented! (see the Statics integration for an example)");
  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::IMPLICIT::GenAlpha::GetIntParam() const
{
  CheckInitSetup();
  dserror("Set the time integration parameter as return value!");
  return -1.0;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::UpdateStepState()
{
  CheckInitSetup();
  dserror("Not yet implemented! (see the Statics integration for an example)");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::UpdateStepElement()
{
  CheckInitSetup();
  dserror("Not yet implemented! (see the Statics integration for an example)");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::IMPLICIT::GenAlpha::PredictConstDisConsistVelAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  CheckInitSetup();
  dserror("Not yet implemented! (see the Statics integration for an example)");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlpha::PredictConstVelConsistAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  CheckInitSetup();
  dserror("Not yet implemented! (see the Statics integration for an example)");
  return false;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::IMPLICIT::GenAlpha::PredictConstAcc(
    Epetra_Vector& disnp,
    Epetra_Vector& velnp,
    Epetra_Vector& accnp) const
{
  CheckInitSetup();
  dserror("Not yet implemented! (see the Statics integration for an example)");
  return false;
}
