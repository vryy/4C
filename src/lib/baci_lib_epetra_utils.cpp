/*---------------------------------------------------------------------*/
/*! \file

\brief General utility methods for all Epetra objects


\level 1

*/
/*---------------------------------------------------------------------*/

#include "baci_lib_epetra_utils.H"

#include <Epetra_Object.h>

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void set_trace_back_mode(int tracebackmode)
{
  static const int global_mode = Epetra_Object::GetTracebackMode();
  Epetra_Object::SetTracebackMode(std::max(global_mode, tracebackmode));
}
