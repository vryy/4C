/*!----------------------------------------------------------------------
\file position_array.H

\brief collection of service methods for intersection computations

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
</pre>
*----------------------------------------------------------------------*/
#ifndef POSITION_ARRAY_H
#define POSITION_ARRAY_H


#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "Epetra_SerialDenseMatrix.h"


namespace GEO
{
  /*!
   * a common task is to get an array of the positions of all nodes of this element
   *
   * template version
   *
   * \note xyze is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   *
   * \author a.ger
   * \date 02/08
   */
  template <class M>
  void fillInitialPositionArray(const DRT::Element* const ele, M& xyze)
  {
    const int numnode = ele->NumNode();

    const DRT::Node* const* nodes = ele->Nodes();
    dsassert(nodes != NULL,
        "element has no nodal pointers, so getting a position array doesn't make sense!");

    for (int inode = 0; inode < numnode; inode++)
    {
      const double* x = nodes[inode]->X();
      xyze(0, inode) = x[0];
      xyze(1, inode) = x[1];
      xyze(2, inode) = x[2];
    }
    return;
  }


  /*!
   * a common task is to get an array of the positions of all nodes of this element
   *
   * template version
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   *
   * \author a.ger
   * \date 02/08
   */
  template <DRT::Element::DiscretizationType distype, class M>
  void fillInitialPositionArray(const DRT::Element* const ele, M& xyze)
  {
    dsassert(distype == ele->Shape(), "mismatch in distype");
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    const DRT::Node* const* nodes = ele->Nodes();
    dsassert(nodes != NULL,
        "element has no nodal pointers, so getting a position array doesn't make sense!");

    for (int inode = 0; inode < numnode; inode++)
    {
      const double* x = nodes[inode]->X();
      xyze(0, inode) = x[0];
      xyze(1, inode) = x[1];
      xyze(2, inode) = x[2];
    }
    return;
  }


  /*!
   * a common task is to get an array of the positions of all nodes of this element
   *
   * template version
   *
   * \note array is defined as (dim,numnode) with user-specified number of space dimensions
   *       that are of interest for the element application calling this method
   *
   * \return Array with 1, 2 or 3 dimensional position of all element nodes in the coordinate system
   * of the nodes
   *
   * \author gjb
   * \date 12/08
   */
  template <DRT::Element::DiscretizationType distype, int dim, class M>
  void fillInitialPositionArray(const DRT::Element* const ele, M& xyze)
  {
    dsassert(distype == ele->Shape(), "mismatch in distype");
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    const DRT::Node* const* nodes = ele->Nodes();
    dsassert(nodes != NULL,
        "element has no nodal pointers, so getting a position array doesn't make sense!");

    dsassert((dim > 0) && (dim < 4), "Illegal number of space dimensions");

    for (int inode = 0; inode < numnode; inode++)
    {
      const double* x = nodes[inode]->X();
      // copy the values in the current column
      std::copy(x, x + dim, &xyze(0, inode));
      // fill the remaining entries of the column with zeros, if the given matrix has
      // the wrong row dimension (this is primarily for safety reasons)
      std::fill(&xyze(0, inode) + dim, &xyze(0, inode) + xyze.M(), 0.0);
    }
    return;
  }

  /*!
   * a common task is to get an array of the positions of all nodes of this element
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   *
   * \author schott
   * \date 01/15
   */
  void InitialPositionArray(Epetra_SerialDenseMatrix& xyze, const DRT::Element* const ele);

  /*!
   * a common task is to get an array of the positions of all nodes of this element
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   *
   * \author a.ger
   * \date 02/08
   */
  LINALG::SerialDenseMatrix InitialPositionArray(const DRT::Element* const
          ele  ///< pointer to element, whose nodes we evaluate for their position
  );


  /*!
   * fill array with current nodal positions
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   *
   * \author a.ger
   * \date 02/08
   */
  LINALG::SerialDenseMatrix getCurrentNodalPositions(
      const DRT::Element* const ele,  ///< element with nodal pointers
      const std::map<int, LINALG::Matrix<3, 1>>&
          currentcutterpositions  ///< current positions of all cutter nodes
  );


  /*!
   * fill array with current nodal positions
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the
   * nodes
   *
   * \author u.may
   * \date 09/09
   */
  LINALG::SerialDenseMatrix getCurrentNodalPositions(
      const Teuchos::RCP<const DRT::Element> ele,  ///< pointer on element
      const std::map<int, LINALG::Matrix<3, 1>>&
          currentpositions  ///< current positions of all cutter nodes
  );


  /*!
   * fill array with current nodal positions
   *
   * \note array is defined as (3,numnode)
   *
   * \author a.ger
   * \date 11/08
   */
  template <class M>
  void fillCurrentNodalPositions(const DRT::Element* const ele,  ///< element with nodal pointers
      const std::map<int, LINALG::Matrix<3, 1>>&
          currentpositions,  ///< current positions of all cutter nodes
      M& xyze                ///< position array (3,numnode) to be filled
  )
  {
    const int numnode = ele->NumNode();
    const int* nodeids = ele->NodeIds();
    for (int inode = 0; inode < numnode; ++inode)
    {
      const LINALG::Matrix<3, 1>& x = currentpositions.find(nodeids[inode])->second;
      xyze(0, inode) = x(0);
      xyze(1, inode) = x(1);
      xyze(2, inode) = x(2);
    }
  }

}  // namespace GEO

#endif
