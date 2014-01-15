/*----------------------------------------------------------------------*/
/*!
\file drt_meshfree_utils.H

\brief service methods for a given meshfree discretisations

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/
#include <vector>

#ifndef DRT_MESHFREE_UTILS_H
#define DRT_MESHFREE_UTILS_H

namespace LINALG
{
  class SerialDenseMatrix;
}

namespace DRT
{

namespace MESHFREE
{

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Reduces node position to dimensions of the face in which the nodes
   *        lie. Herefore, node position are collected in a SerialDenseMatrix
   *        and given back in another SerialDenseMatrix.
   *
   * \param xyz         (in) : original node positions; size: dim x number of nodes
   * \param xyz_reduced (out): reduced node positions; size: face dim x number of nodes
   * \param dims        (out): dimensions in which face lies
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  std::vector<int> ReduceDimensionOfFaceNodes(
    const LINALG::SerialDenseMatrix & xyz,
    LINALG::SerialDenseMatrix       & xyz_reduced
    );

} // namespace MESHFREE

} // namespace DRT

#endif // ifndef DRT_MESHFREE_UTILS_H
