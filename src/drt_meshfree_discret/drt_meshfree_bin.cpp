/*!--------------------------------------------------------------------------
\file drt_meshfree_bin.cpp
\brief

<pre>
-----------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
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

-----------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>

*--------------------------------------------------------------------------*/

#include "drt_meshfree_bin.H"

/*--------------------------------------------------------------------------*
 |  ctor                                               (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeBin::MeshfreeBin(int id, int owner)
  : DRT::Element::Element(id,owner)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                          (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeBin::MeshfreeBin(const DRT::MESHFREE::MeshfreeBin& old)
  : DRT::Element::Element(old)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  dtor                                               (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeBin::~MeshfreeBin()
{
  return;
}

/*--------------------------------------------------------------------------*
 | Delete a single node from the element               (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeBin::DeleteNode(int gid)
{
  for (unsigned int i = 0; i<nodeid_.size(); i++){
    if (nodeid_[i]==gid){
      nodeid_.erase(nodeid_.begin()+i);
      node_.erase(node_.begin()+i);
      return;
    }
  }
  dserror("Connectivity issues: No node with specified gid to delete in element. ");
  return;
}
