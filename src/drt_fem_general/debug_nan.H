/*!----------------------------------------------------------------------
\file debug_nan.H

A set of utility functions to identify NaNs in vectors and matrices

A set of utility functions to identify NaNs in vectors and matrices. Note,
that the performed operations might be expensive and are meant to be
used during debugging, not in optimized runs.

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef DEBUG_NAN_H
#define DEBUG_NAN_H

namespace DRT
{
  /// namespace for debugging utilities
  namespace DEBUGGING
  {

    static void NaNChecker(const Epetra_SerialDenseVector& vec)
    {
      if (std::isnan(vec.Norm2()))
      {
        dserror("NaNs detected! Quitting...");
      }
    }

    static void NaNChecker(const Epetra_SerialDenseMatrix& mat)
    {
      if (std::isnan(mat.InfNorm()))
      {
        dserror("NaNs detected! Quitting...");
      }
    }

    static void NaNChecker(const std::vector<double>& vec)
    {
      for (std::size_t i = 0; i < vec.size(); i++)
      {
        if (std::isnan(vec[i]))
        {
          dserror("NaNs detected! Quitting...");
        }
      }
    }

    static void NaNChecker(const std::vector<int>& vec)
    {
      for (std::size_t i = 0; i < vec.size(); i++)
      {
        if (std::isnan(vec[i]))
        {
          dserror("NaNs detected! Quitting...");
        }
      }
    }

    static void NaNChecker(const Epetra_Vector& vec)
    {
      double norm = 0;
      vec.Norm1(&norm);
      if (std::isnan(norm))
      {
        dserror("NaNs detected! Quitting...");
      }
    }

  }
}


#endif  // #ifndef DEBUG_NAN_H
#endif  // #ifdef CCADISCRET
