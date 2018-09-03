/*-----------------------------------------------------------*/
/*!
\file str_timint_locainterface.cpp

\maintainer Michael Hiermeier

\date Nov 27, 2015

\level 3

*/
/*-----------------------------------------------------------*/


#include "str_timint_locainterface.H"

#include <LOCA_Parameter_Vector.H>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::TIMINT::LocaInterface::LocaInterface() : displ_fac_(0.0)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::LocaInterface::Setup() { STR::TIMINT::NoxInterface::Setup(); }

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::TIMINT::LocaInterface::setParameters(const LOCA::ParameterVector& p)
{
  CheckInit();
  displ_fac_ = p.getValue("displ_fac_");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::LocaInterface::computeF(
    const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  CheckInitSetup();
  dserror("You shall not pass!");

  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::LocaInterface::computeJacobian(const Epetra_Vector& x, Epetra_Operator& jac)
{
  CheckInitSetup();
  dserror("You shall not pass!");
  return true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
bool STR::TIMINT::LocaInterface::computeFandJacobian(
    const Epetra_Vector& x, Epetra_Vector& rhs, Epetra_Operator& jac)
{
  CheckInitSetup();
  dserror("You shall not pass!");
  return true;
}
