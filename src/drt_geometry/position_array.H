/*!
\file position_array.H

\brief collection of service methods for intersection computations

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

*/


#ifndef POSITION_ARRAY_H
#define POSITION_ARRAY_H


#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"


namespace GEO
{
  /*!
   * a common task is to get an array of the positions of all nodes of this element
   *
   * template version
   *
   * \note xyze is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the nodes
   *
   * \author a.ger
   * \date 02/08
   */
  template <class M>
  void fillInitialPositionArray(
      const DRT::Element* const ele,
      M&                        xyze
  )
  {
    const int numnode = ele->NumNode();

    const DRT::Node*const* nodes = ele->Nodes();
    dsassert(nodes != NULL, "element has no nodal pointers, so getting a position array doesn't make sense!");

    for (int inode=0; inode<numnode; inode++)
    {
      const double* x = nodes[inode]->X();
      xyze(0,inode) = x[0];
      xyze(1,inode) = x[1];
      xyze(2,inode) = x[2];
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
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the nodes
   *
   * \author a.ger
   * \date 02/08
   */
  template <DRT::Element::DiscretizationType distype, class M>
  void fillInitialPositionArray(
          const DRT::Element* const ele,
          M&                        xyze
          )
  {
      dsassert(distype == ele->Shape(), "mismatch in distype");
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      const DRT::Node*const* nodes = ele->Nodes();
      dsassert(nodes != NULL, "element has no nodal pointers, so getting a position array doesn't make sense!");

      for (int inode=0; inode<numnode; inode++)
      {
          const double* x = nodes[inode]->X();
          xyze(0,inode) = x[0];
          xyze(1,inode) = x[1];
          xyze(2,inode) = x[2];
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
   * \return Array with 1, 2 or 3 dimensional position of all element nodes in the coordinate system of the nodes
   *
   * \author gjb
   * \date 12/08
   */
  template <DRT::Element::DiscretizationType distype, int dim, class M>
  void fillInitialPositionArray(
          const DRT::Element* const ele,
          M&                        xyze
          )
  {
      dsassert(distype == ele->Shape(), "mismatch in distype");
      const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

      const DRT::Node*const* nodes = ele->Nodes();
      dsassert(nodes != NULL, "element has no nodal pointers, so getting a position array doesn't make sense!");

      dsassert((dim > 0) && (dim < 4),"Illegal number of space dimensions");

      for (int inode=0; inode<numnode; inode++)
      {
          const double* x = nodes[inode]->X();
          for(int idim=0; idim<dim; idim++)
          {
            xyze(idim,inode) = x[idim];
          }
      }
      return;
  }


  /*!
    * a common task is to get an array of the positions of all nodes of this element
    *
    * \note array is defined as (3,numnode)
    *
    * \return Array with 3 dimensional position of all element nodes in the coordinate system of the nodes
    *
    * \author a.ger
    * \date 02/08
    */
   LINALG::SerialDenseMatrix InitialPositionArray(
       const DRT::Element* const ele           ///< pointer to element, whose nodes we evaluate for their position
       );


  /*!
   * fill array with current nodal positions
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the nodes
   *
   * \author a.ger
   * \date 02/08
   */
  LINALG::SerialDenseMatrix getCurrentNodalPositions(
      const DRT::Element* const               ele,                   ///< element with nodal pointers
      const std::map<int,LINALG::Matrix<3,1> >&    currentcutterpositions ///< current positions of all cutter nodes
      );
  
  
  /*!
   * fill array with current nodal positions
   *
   * \note array is defined as (3,numnode)
   *
   * \return Array with 3 dimensional position of all element nodes in the coordinate system of the nodes
   *
   * \author u.may
   * \date 09/09
   */
  LINALG::SerialDenseMatrix getCurrentNodalPositions(
      const Teuchos::RCP<const DRT::Element>     ele,              ///< pointer on element
      const std::map<int,LINALG::Matrix<3,1> >&  currentpositions  ///< current positions of all cutter nodes
      );
  
  
  /*!
   * fill array with current nodal positions
   *
   * \note array is defined as (3,numnode)
   *
   * \author a.ger
   * \date 11/08
   */
  template<class M>
  void fillCurrentNodalPositions(
      const DRT::Element* const             ele,                ///< element with nodal pointers
      const std::map<int,LINALG::Matrix<3,1> >&  currentpositions,   ///< current positions of all cutter nodes
      M&                                    xyze                ///< position array (3,numnode) to be filled
      )
  {
      const int numnode = ele->NumNode();
      const int* nodeids = ele->NodeIds();
      for (int inode = 0; inode < numnode; ++inode)
      {
        const LINALG::Matrix<3,1>& x = currentpositions.find(nodeids[inode])->second;
        xyze(0,inode) = x(0);
        xyze(1,inode) = x(1);
        xyze(2,inode) = x(2);
      }
  }

}

#endif
