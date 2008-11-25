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

#ifdef CCADISCRET
#ifndef POSITION_ARRAY_H
#define POSITION_ARRAY_H

#include "../drt_lib/drt_node.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_fixedsizematrix.H"


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
      const DRT::Element*      ele,
      M&                       xyze
  )
  {
    const DRT::Element::DiscretizationType distype = ele->Shape();
    const int numnode = DRT::UTILS::getNumberOfElementNodes(distype);
    
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
}


#endif
#endif


