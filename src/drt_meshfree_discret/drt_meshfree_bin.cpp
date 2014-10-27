/*!--------------------------------------------------------------------------
\file drt_meshfree_bin.cpp
\brief

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
