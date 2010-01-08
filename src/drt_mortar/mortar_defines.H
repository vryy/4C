/*!----------------------------------------------------------------------
\file mortar_defines.H

\brief A set of preprocessor defines and utility functions for mortar methods

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

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
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef MORTAR_DEFINES_H
#define MORTAR_DEFINES_H

using namespace std;

/************************************************************************/
/* Mortar algorithm options                                             */
/************************************************************************/

#define MORTARONELOOP            /* flag for computing D and M in one loop */
//#define MORTARBOUNDMOD         /* flag for modification of boundary nodes */
//#define MORTARPETROVGALERKIN   /* flag for PetrovGalerkin approach for LM */
//#define MORTARCONTACT          /* stuff that is contact-specific (this has to go away soon) */

/************************************************************************/
/* Mortar algorithm parameters                                          */
/************************************************************************/

// MORTAR PROJECTION (2D/3D)
#define MORTARMAXITER 10       /* max. no. iterations for local Newton */
#define MORTARCONVTOL 1.0e-12  /* convergence tolerance for local Newton */

// MORTAR PROJECTION AND INTEGRATION (2D)
#define MORTARPROJTOL 0.05     /* projection tolerance for overlap */
#define MORTARPROJLIM 1.0e-8   /* exact projection limit (no tolerance!) */

// MORTAR PROJECTION AND INTEGRATION (3D)
#define MORTARCLIPTOL 1.0e-12  /* tolerance for polygon clipping */
#define MORTARINTLIM 1.0e-12   /* min(area-%) cell/slave for integration */

/************************************************************************/
/* Mortar utility functions                                             */
/************************************************************************/
/*!
\brief MORTAR: namespace of the BACI mortar coupling algorithms
*/
namespace MORTAR
{
  /*!
  \brief Template to swap 2 instances of type kind
  */
  template<typename kind>
  void Swap(kind& a, kind& b)
  {
    kind tmp = a;
    a = b;
    b = tmp;
    return;
  }
  
  /*!
  \brief Sort vector in ascending order
  
  This routine is taken from Trilinos MOERTEL package.
  
  \param dlist (in): vector to be sorted (unsorted on input, sorted on output)
  \param N (in):     length of vector to be sorted
  \param list2 (in): another vector which is sorted accordingly  
  */  
  void Sort(double* dlist, int N, int* list2);

} // namespace MORTAR

#endif  // #ifndef MORTAR_DEFINES_H
#endif  // #ifdef CCADISCRET
