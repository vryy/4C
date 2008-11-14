/*!
\file vector_definitions.H

\brief definitions used in geometry contexts

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions and of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de) 
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de                   

-------------------------------------------------------------------------
</pre>

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
#include "../drt_lib/linalg_fixedsizematrix.H"


typedef blitz::Array<double,1> BlitzVec; ///< a short-cut for a one-dimensional BLITZ vector
typedef blitz::Array<double,2> BlitzMat; ///< a short-cut for a two-dimensional BLITZ matrix

typedef blitz::TinyVector<double,3> BlitzVec3; ///< a short-cut for a one-dimensional BLITZ vector of length 3
typedef blitz::TinyVector<double,2> BlitzVec2; ///< a short-cut for a one-dimensional BLITZ vector of length 2
typedef blitz::TinyVector<double,1> BlitzVec1; ///< a short-cut for a one-dimensional BLITZ vector of length 1

typedef blitz::TinyMatrix<double,2,2> BlitzMat2x2; ///< a short-cut for a two-dimensional BLITZ matrix of size 2x2
typedef blitz::TinyMatrix<double,3,2> BlitzMat3x2; ///< a short-cut for a two-dimensional BLITZ matrix of size 3x2
typedef blitz::TinyMatrix<double,3,3> BlitzMat3x3; ///< a short-cut for a two-dimensional BLITZ matrix of size 3x3

namespace LINALG
{
  typedef LINALG::FixedSizeSerialDenseMatrix<3,1> Vec3; ///< a short-cut for a one-dimensional fixed-size vector of length 3
  typedef LINALG::FixedSizeSerialDenseMatrix<2,1> Vec2; ///< a short-cut for a one-dimensional fixed-size vector of length 2
  typedef LINALG::FixedSizeSerialDenseMatrix<1,1> Vec1; ///< a short-cut for a one-dimensional fixed-size vector of length 1
  
  typedef LINALG::FixedSizeSerialDenseMatrix<2,2> Mat2x2; ///< a short-cut for a two-dimensional fixed-size matrix of size 2x2
  typedef LINALG::FixedSizeSerialDenseMatrix<3,2> Mat3x2; ///< a short-cut for a two-dimensional fixed-size matrix of size 3x2
  typedef LINALG::FixedSizeSerialDenseMatrix<3,3> Mat3x3; ///< a short-cut for a two-dimensional fixed-size matrix of size 3x3
}
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
