/*!----------------------------------------------------------------------
\file debug_nan.cpp

\brief A set of utility functions to identify NaNs in vectors and matrices

A set of utility functions to identify NaNs in vectors and matrices. Note,
that the performed operations might be expensive and are meant to be
used during debugging, not in optimized runs.

-------------------------------------------------------------------------
</pre>

<pre>
\level 2
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "debug_nan.H"
#include "../drt_lib/drt_dserror.H"


void DRT::DEBUGGING::NaNChecker(const Epetra_SerialDenseVector& vec)
{
  for (int i = 0; i < vec.Length(); i++)
  {
    if (std::isnan(vec(i)))
    {
      dserror("NaNs detected! Quitting...");
    }
  }
}

void DRT::DEBUGGING::NaNChecker(const Epetra_SerialDenseMatrix& mat)
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

void DRT::DEBUGGING::NaNChecker(const std::vector<double>& vec)
{
  for (std::size_t i = 0; i < vec.size(); i++)
  {
    if (std::isnan(vec[i]))
    {
      dserror("NaNs detected! Quitting...");
    }
  }
}

void DRT::DEBUGGING::NaNChecker(const std::vector<int>& vec)
{
  for (std::size_t i = 0; i < vec.size(); i++)
  {
    if (std::isnan(vec[i]))
    {
      dserror("NaNs detected! Quitting...");
    }
  }
}

void DRT::DEBUGGING::NaNChecker(const Epetra_Vector& vec)
{
  for (int i = 0; i < vec.MyLength(); i++)
  {
    if (std::isnan(vec[i]))
    {
      dserror("NaNs detected! Quitting...");
    }
  }
}
