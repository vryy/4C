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

#ifndef DRT_MESHFREE_UTILS_H
#define DRT_MESHFREE_UTILS_H

#include <Teuchos_RCP.hpp>

namespace LINALG
{
  class SerialDenseMatrix;
}

namespace DRT
{

namespace MESHFREE
{

  /*--------------------------------------------------------------------------*
   *                                                                          *
   * enums to specify error code                                              *
   *                                                                nis Feb14 *
   *--------------------------------------------------------------------------*/

  //! type of errors - to be extended ad lib
  enum errorcode {
    fine,
    det_zero
  };

  void OutputMeshfreeError(
    int err
    );

  /*------------------------------------------------------------------------*/
  /*!
   * \brief Reduces node position to dimensions of the face in which the nodes
   *        lie. Herefore, node position are collected in a SerialDenseMatrix
   *        and given back in another SerialDenseMatrix.
   *
   * \param xyz         (in/out)
   *///                                                   (public) nis Jan14
  /*------------------------------------------------------------------------*/
  void MovePointOfReferenceToFace(
    Teuchos::RCP<LINALG::SerialDenseMatrix>& xyz
    );

  void MovePointOfReferenceToFace(
    Teuchos::RCP<LINALG::SerialDenseMatrix>& xyz1,
    Teuchos::RCP<LINALG::SerialDenseMatrix>& xyz2
    );

  int ReduceCloudDimensionOnFaces(
    Teuchos::RCP<LINALG::SerialDenseMatrix>& xyz
    );

  int ReduceCloudDimensionOnFaces(
    Teuchos::RCP<LINALG::SerialDenseMatrix>& xyz1,
    Teuchos::RCP<LINALG::SerialDenseMatrix>& xyz2
    );

} // namespace MESHFREE

} // namespace DRT

#endif // ifndef DRT_MESHFREE_UTILS_H
