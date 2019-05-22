/*-----------------------------------------------------------*/
/*!
\file solver_aztecoo_conditionnumber.cpp

\maintainer Martin Kronbichler

\brief Condition number via AztecOO

\level 3

*/
/*-----------------------------------------------------------*/

#include "solver_aztecoo_conditionnumber.H"
#include "../drt_lib/drt_dserror.H"

#include <AztecOO.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double LINALG::AztecOOConditionNumber::getStatus(int az_key) const
{
  if (not solver_) dserror("The solver has not yet been initialized!");

  return solver_->GetAztecStatus()[az_key];
}
