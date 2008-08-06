/*!
\file xfem.H

\brief definitions used in geometry contexts

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
#ifdef CCADISCRET
#ifndef VECTOR_DEFINITIONS_H
#define VECTOR_DEFINITIONS_H

#include <blitz/array.h>


//! generally usefull typedefs, when working with Blitz arrays
typedef blitz::Array<double,1> BlitzVec; ///< a short-cut for a one-dimensional BLITZ vector
typedef blitz::Array<double,2> BlitzMat; ///< a short-cut for a two-dimensional BLITZ matrix

typedef blitz::TinyVector<double,3> BlitzVec3; ///< a short-cut for a one-dimensional BLITZ vector of length 3
typedef blitz::TinyVector<double,2> BlitzVec2; ///< a short-cut for a one-dimensional BLITZ vector of length 2
typedef blitz::TinyVector<double,1> BlitzVec1; ///< a short-cut for a one-dimensional BLITZ vector of length 1

typedef blitz::TinyMatrix<double,2,2> BlitzMat2x2; ///< a short-cut for a two-dimensional BLITZ matrix of size 2x2
typedef blitz::TinyMatrix<double,3,2> BlitzMat3x2; ///< a short-cut for a two-dimensional BLITZ matrix of size 3x2
typedef blitz::TinyMatrix<double,3,3> BlitzMat3x3; ///< a short-cut for a two-dimensional BLITZ matrix of size 3x3

//! take a double C array of (hopefully) length 3 and convert it to a blitz TinyVector of length 3
inline BlitzVec3 toBlitzArray(const double* x)
{
  BlitzVec3 blitz_x;
  blitz_x(0) = x[0];
  blitz_x(1) = x[1];
  blitz_x(2) = x[2];
  return blitz_x;
}

#endif  // #ifndef VECTOR_DEFINITIONS_H
#endif  // #ifdef CCADISCRET
