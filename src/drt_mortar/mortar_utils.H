/*!----------------------------------------------------------------------
\file mortar_utils.H

\brief A set of utility functions for mortar methods

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
#ifndef MORTAR_UTILS_H
#define MORTAR_UTILS_H

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_CrsMatrix.h>
#include "../linalg/linalg_sparsematrix.H"

using namespace std;
using namespace Teuchos;


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

  /*!
  \brief Transform the row map of a matrix

  */
  RCP<LINALG::SparseMatrix> MatrixRowTransformGIDs(RCP<LINALG::SparseMatrix> inmat,
  		                                             RCP<Epetra_Map> newrowmap);

  /*!
  \brief Transform the column map of a matrix

  */
  RCP<LINALG::SparseMatrix> MatrixColTransformGIDs(RCP<LINALG::SparseMatrix> inmat,
  		                                             RCP<Epetra_Map> newdomainmap);

  /*!
  \brief Transform the row and column maps of a matrix

  */
  RCP<LINALG::SparseMatrix> MatrixRowColTransformGIDs(RCP<LINALG::SparseMatrix> inmat,
  		                                                RCP<Epetra_Map> newrowmap,
  		                                                RCP<Epetra_Map> newdomainmap);

  /*!
  \brief Transform the row map of a matrix

  */
  RCP<LINALG::SparseMatrix> MatrixRowTransform(RCP<LINALG::SparseMatrix> inmat,
  		                                         RCP<Epetra_Map> permrowmap,
  		                                         RCP<Epetra_Map> newrowmap);

  /*!
  \brief Transform the column map of a matrix

  */
  RCP<LINALG::SparseMatrix> MatrixColTransform(RCP<LINALG::SparseMatrix> inmat,
  		                                         RCP<Epetra_Map> permdomainmap,
  		                                         RCP<Epetra_Map> newdomainmap);

  /*!
  \brief Transform the row and column maps of a matrix

  */
  RCP<LINALG::SparseMatrix> MatrixRowColTransform(RCP<LINALG::SparseMatrix> inmat,
  		                                            RCP<Epetra_Map> permrowmap,
  		                                            RCP<Epetra_Map> newrowmap,
  		                                            RCP<Epetra_Map> permdomainmap,
  		                                            RCP<Epetra_Map> newdomainmap);

  /*!
  \brief Parallel redistribution of a matrix according to permutation map

  */
  RCP<Epetra_CrsMatrix> Redistribute(const LINALG::SparseMatrix& src,
                                     const Epetra_Map& permsrcmap);

} // namespace MORTAR

#endif  // #ifndef MORTAR_UTILS_H
#endif  // #ifdef CCADISCRET
