/*----------------------------------------------------------------------*/
/*! \file

\brief collection of service methods for intersection computations


\level 2

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/


#include "../drt_geometry/position_array.H"

/*!
 * \brief create an often used array with 3D nodal positions
 */
LINALG::SerialDenseMatrix GEO::InitialPositionArray(const DRT::Element* const ele)
{
  const int numnode = ele->NumNode();
  LINALG::SerialDenseMatrix xyze(3, numnode);
  const DRT::Node* const* nodes = ele->Nodes();
  if (nodes == NULL)
  {
    dserror("element has no nodal pointers, so getting a position array doesn't make sense!");
  }
  for (int inode = 0; inode < numnode; inode++)
  {
    const double* x = nodes[inode]->X();
    std::copy(x, x + 3, &xyze(0, inode));
  }
  return xyze;
}


/*!
 * \brief create an often used array with 3D nodal positions
 */
void GEO::InitialPositionArray(Epetra_SerialDenseMatrix& xyze, const DRT::Element* const ele)
{
  const int numnode = ele->NumNode();
  xyze.Shape(3, numnode);
  const DRT::Node* const* nodes = ele->Nodes();
  if (nodes == NULL)
  {
    dserror("element has no nodal pointers, so getting a position array doesn't make sense!");
  }
  for (int inode = 0; inode < numnode; inode++)
  {
    const double* x = nodes[inode]->X();
    std::copy(x, x + 3, &xyze(0, inode));
  }
}



/*!
\brief  fill array with current nodal positions

\return array with element nodal positions (3,numnode)
 */
LINALG::SerialDenseMatrix GEO::getCurrentNodalPositions(
    const DRT::Element* const ele,  ///< element with nodal pointers
    const std::map<int, LINALG::Matrix<3, 1>>&
        currentcutterpositions  ///< current positions of all cutter nodes
)
{
  const int numnode = ele->NumNode();
  LINALG::SerialDenseMatrix xyze(3, numnode);
  const int* nodeids = ele->NodeIds();
  for (int inode = 0; inode < numnode; ++inode)
  {
    const LINALG::Matrix<3, 1>& x = currentcutterpositions.find(nodeids[inode])->second;
    xyze(0, inode) = x(0);
    xyze(1, inode) = x(1);
    xyze(2, inode) = x(2);
  }
  return xyze;
}



/*!
\brief  fill array with current nodal positions

\return array with element nodal positions (3,numnode)
 */
LINALG::SerialDenseMatrix GEO::getCurrentNodalPositions(
    const Teuchos::RCP<const DRT::Element> ele,  ///< pointer on element
    const std::map<int, LINALG::Matrix<3, 1>>&
        currentpositions  ///< current positions of all cutter nodes
)
{
  const int numnode = ele->NumNode();
  LINALG::SerialDenseMatrix xyze(3, numnode);
  const int* nodeids = ele->NodeIds();
  for (int inode = 0; inode < numnode; ++inode)
  {
    const LINALG::Matrix<3, 1>& x = currentpositions.find(nodeids[inode])->second;
    xyze(0, inode) = x(0);
    xyze(1, inode) = x(1);
    xyze(2, inode) = x(2);
  }
  return xyze;
}
