/*----------------------------------------------------------------------*/
/**
\file epetra_utils.cpp

\brief General utility methods for all Epetra objects

\level 1

\maintainer Michael Hiermeier

\date Feb 27, 2018
*/
/*----------------------------------------------------------------------*/

#include <Epetra_Object.h>
#include "epetra_utils.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void set_trace_back_mode( int tracebackmode )
{
  static const int global_mode = Epetra_Object::GetTracebackMode();
  Epetra_Object::SetTracebackMode( std::max( global_mode, tracebackmode ) );
}
