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
      for (int i = 0; i < vec.Length(); i++)
      {
        if (std::isnan(vec(i)))
        {
          dserror("NaNs detected! Quitting...");
        }
      }
    }

    static void NaNChecker(const Epetra_SerialDenseMatrix& mat)
    {
      for (int m = 0; m < mat.M(); m++)
      {
        for (int n = 0; n < mat.N(); n++)
        {
          if (std::isnan(mat(m,n)))
          {
            dserror("NaNs detected! Quitting...");
          }
        }
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
      for (int i = 0; i < vec.MyLength(); i++)
      {
        if (std::isnan(vec[i]))
        {
          dserror("NaNs detected! Quitting...");
        }
      }
    }

  }
}


#endif  // #ifndef DEBUG_NAN_H
#endif  // #ifdef CCADISCRET
