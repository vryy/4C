/*---------------------------------------------------------------------*/
/*! \file

\brief General utility methods for all Epetra objects

\maintainer Martin Kronbichler

\level 1

*/
/*---------------------------------------------------------------------*/

#ifndef SRC_DRT_LIB_EPETRA_UTILS_H_
#define SRC_DRT_LIB_EPETRA_UTILS_H_

#include "drt_dserror.H"

/// set and restore the trace back mode during the macro execution
void set_trace_back_mode(int tracebackmode);


/** \brief Macro to wrap any Epetra method which returns an error code
 *
 *  This macro can be used to wrap a Epetra function. At the beginning the
 *  internal stack trace of the Epetra_Object is temporally activated. If an
 *  error or a warning is detected (return value is negative or positive), this
 *  stack trace will be printed to the screen. Additionally the program is
 *  interrupted in an error case.
 *
 *  \note Any previously set trace back mode (i.e. before this macro is
 *  executed the very first time) is in the end of this macro
 *  recovered.                                               hiermeier 02/17 */
#define CATCH_EPETRA_ERROR(a)                                                                     \
  {                                                                                               \
    set_trace_back_mode(2);                                                                       \
    int epetra_error = a;                                                                         \
    if (epetra_error < 0)                                                                         \
    {                                                                                             \
      dserror("EPETRA_ERROR %d caught!", epetra_error);                                           \
    }                                                                                             \
    else if (epetra_error > 0)                                                                    \
    {                                                                                             \
      std::cout << __FILE__ ", " << __LINE__ << ": EPETRA_WARNING " << epetra_error << " caught!" \
                << std::endl;                                                                     \
    }                                                                                             \
    set_trace_back_mode(0);                                                                       \
  }


#endif /* SRC_DRT_LIB_EPETRA_UTILS_H_ */
