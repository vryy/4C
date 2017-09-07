/*!----------------------------------------------------------------------
\file mortar_utils.H

\brief A set of utility functions for mortar methods

\level 1

<pre>
\maintainer Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
#ifndef MORTAR_UTILS_H
#define MORTAR_UTILS_H

#include "mortar_coupling3d_classes.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#include <Epetra_Export.h>
#include <Epetra_CrsMatrix.h>


// forward declarations
namespace LINALG
{
  class SparseMatrix;
}

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
  \brief Transform the row map of a matrix (only GIDs)

  This method changes the row map of an input matrix to a new row map
  with different GID numbering. However, the parallel distribution of
  this new row map is exactly the same as in the old row map. Thus, this
  is simply a processor-local 1:1 matching of old and new GIDs.

  */
  Teuchos::RCP<LINALG::SparseMatrix> MatrixRowTransformGIDs(Teuchos::RCP<const LINALG::SparseMatrix> inmat,
                                                            Teuchos::RCP<const Epetra_Map> newrowmap);

  /*!
  \brief Transform the column map of a matrix (only GIDs)

  This method changes the column map of an input matrix to a new column
  map with different GID numbering (and the domain map, accordingly).
  However, the parallel distribution of the new domain map is exactly
  the same as in the old domain map. Thus, this is simply a processor-local
  1:1 matching of old and new GIDs.

  */
  Teuchos::RCP<LINALG::SparseMatrix> MatrixColTransformGIDs(Teuchos::RCP<const LINALG::SparseMatrix> inmat,
                                                            Teuchos::RCP<const Epetra_Map> newdomainmap);

  /*!
  \brief Transform the row and column maps of a matrix (only GIDs)

  This method changes the row and column maps of an input matrix to new
  row and column maps with different GID numbering (and the domain map,
  accordingly). However, the parallel distribution of the new row and
  domain maps is exactly the same as in the old ones. Thus, this is simply
  a processor-local 1:1 matching of old and new GIDs.

  */
  Teuchos::RCP<LINALG::SparseMatrix> MatrixRowColTransformGIDs(Teuchos::RCP<const LINALG::SparseMatrix> inmat,
                                                               Teuchos::RCP<const Epetra_Map> newrowmap,
                                                               Teuchos::RCP<const Epetra_Map> newdomainmap);

  /*!
  \brief Transform the row map of a matrix (parallel distribution)

  This method changes the row map of an input matrix to a new
  row map with identical GIDs but different parallel distribution.

  */
  Teuchos::RCP<LINALG::SparseMatrix> MatrixRowTransform(Teuchos::RCP<const LINALG::SparseMatrix> inmat,
                                                        Teuchos::RCP<const Epetra_Map> newrowmap);

  /*!
  \brief Transform the column map of a matrix (parallel distribution)

  This method changes the column map of an input matrix to a new
  column map with identical GIDs but different parallel distribution
  (and the domain map, accordingly).

  */
  Teuchos::RCP<LINALG::SparseMatrix> MatrixColTransform(Teuchos::RCP<const LINALG::SparseMatrix> inmat,
                                                        Teuchos::RCP<const Epetra_Map> newdomainmap);

  /*!
  \brief Transform the row and column maps of a matrix (parallel distribution)

  This method changes the row and column maps of an input matrix
  to new row and column maps with identical GIDs but different
  parallel distribution (and the domain map, accordingly).

  */
  Teuchos::RCP<LINALG::SparseMatrix> MatrixRowColTransform(Teuchos::RCP<const LINALG::SparseMatrix> inmat,
                                                           Teuchos::RCP<const Epetra_Map> newrowmap,
                                                           Teuchos::RCP<const Epetra_Map> newdomainmap);

  /*!
  \brief Parallel redistribution of a sparse matrix

  Helper method for the MatrixTransform() methods above.

  */
  Teuchos::RCP<Epetra_CrsMatrix> Redistribute(const LINALG::SparseMatrix& src,
                                              const Epetra_Map& permrowmap,
                                              const Epetra_Map& permdomainmap);

  /*!
  \brief Convex hull points are sorted in order to obtain final clip polygon

  \param out (in): bool to switch output on/off
  \param transformed (in): coordinates of vertex objects transformed into auxiliary plane
  \param collconvexhull (in): vector of vertex objects to be sorted
  \param respoly (out): vector of vertex objects for result polygon
  \param tol (in): clipping tolerance for close vertices detection
  \return number of removed points from collconvexhull

  */
  int SortConvexHullPoints(bool out,
                           Epetra_SerialDenseMatrix& transformed,
                           std::vector<Vertex>& collconvexhull,
                           std::vector<Vertex>& respoly,
                           double& tol);

namespace UTILS
{
  /*!
  \brief copy the ghosting of dis_src to all discretizations with names in
         vector dis_tar. They are fetched from the global problem, which is
         not very safe, so be sure, what you do. Material pointers can be
         added according to link_materials
  */
  void CreateVolumeGhosting(
      const DRT::Discretization& dis_src,
      const std::vector<std::string> dis_tar,
      std::vector<std::pair<int,int> > material_links,
      bool check_on_in=true,
      bool check_on_exit=true
      );


  /*!
  \brief Prepare mortar element for nurbs case


  store knot vector, zerosized information and normal factor
  */
  void PrepareNURBSElement(
      DRT::Discretization& discret,
      Teuchos::RCP<DRT::Element> ele,
      Teuchos::RCP<MORTAR::MortarElement> cele,
      int dim);

  /*!
  \brief Prepare mortar node for nurbs case

  store control point weight

  */
  void PrepareNURBSNode(
      DRT::Node* node,
      Teuchos::RCP<MORTAR::MortarNode> mnode);
}
} // namespace MORTAR

#endif  // #ifndef MORTAR_UTILS_H
