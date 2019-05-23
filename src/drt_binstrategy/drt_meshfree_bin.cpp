/*!--------------------------------------------------------------------------
\brief

\maintainer Jonas Eichinger

\level 2

*--------------------------------------------------------------------------*/

#include "drt_meshfree_bin.H"

#include "../drt_mortar/mortar_element.H"

/*--------------------------------------------------------------------------*
 |  ctor                                               (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
template <typename ELEMENT>
DRT::MESHFREE::MeshfreeBin<ELEMENT>::MeshfreeBin(int id, int owner) : ELEMENT(id, owner)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                          (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
template <typename ELEMENT>
DRT::MESHFREE::MeshfreeBin<ELEMENT>::MeshfreeBin(const DRT::MESHFREE::MeshfreeBin<ELEMENT>& old)
    : ELEMENT(old)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  dtor                                               (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
template <typename ELEMENT>
DRT::MESHFREE::MeshfreeBin<ELEMENT>::~MeshfreeBin()
{
  return;
}

/*--------------------------------------------------------------------------*
 | Delete a single node from the element               (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
template <typename ELEMENT>
void DRT::MESHFREE::MeshfreeBin<ELEMENT>::DeleteNode(int gid)
{
  for (unsigned int i = 0; i < ELEMENT::nodeid_.size(); i++)
  {
    if (ELEMENT::nodeid_[i] == gid)
    {
      ELEMENT::nodeid_.erase(ELEMENT::nodeid_.begin() + i);
      ELEMENT::node_.erase(ELEMENT::node_.begin() + i);
      return;
    }
  }
  dserror("Connectivity issues: No node with specified gid to delete in element. ");
  return;
}

/*--------------------------------------------------------------------------*
 | Explicit instantiations                                kronbichler 03/15 |
 *--------------------------------------------------------------------------*/
template class DRT::MESHFREE::MeshfreeBin<DRT::Element>;
template class DRT::MESHFREE::MeshfreeBin<DRT::FaceElement>;
template class DRT::MESHFREE::MeshfreeBin<MORTAR::MortarElement>;
