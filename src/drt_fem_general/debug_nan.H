/*!----------------------------------------------------------------------
\file debug_nan.H

\brief A set of utility functions to identify NaNs in vectors and matrices

A set of utility functions to identify NaNs in vectors and matrices. Note,
that the performed operations might be expensive and are meant to be
used during debugging, not in optimized runs.

-------------------------------------------------------------------------

\level 2
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

*----------------------------------------------------------------------*/
#ifndef DEBUG_NAN_H
#define DEBUG_NAN_H

#include <vector>


#include <Epetra_SerialDenseVector.h>
#include <Epetra_SerialDenseMatrix.h>
#include "Epetra_Vector.h"

namespace DRT
{
  namespace DEBUGGING
  {
    /// check for any NaN within this array
    void NaNChecker(const Epetra_SerialDenseVector& vec);

    /// check for any NaN within this array
    void NaNChecker(const Epetra_SerialDenseMatrix& mat);

    /// check for any NaN within this array
    void NaNChecker(const std::vector<double>& vec);

    /// check for any NaN within this array
    void NaNChecker(const std::vector<int>& vec);

    /// check for any NaN within this array
    void NaNChecker(const Epetra_Vector& vec);

  }  // namespace DEBUGGING
}  // namespace DRT


#endif  // #ifndef DEBUG_NAN_H
