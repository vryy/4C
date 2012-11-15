/*!--------------------------------------------------------------------------
\file drt_meshfree_cell.cpp
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
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*--------------------------------------------------------------------------*/

#include "drt_meshfree_cell.H"
#include "drt_meshfree_node.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_node.H"

/*--------------------------------------------------------------------------*
 |  ctor                                                 (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::Cell::Cell(int id, int owner)
  : DRT::MESHFREE::MeshfreeBin(id,owner)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                            (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::Cell::Cell(const DRT::MESHFREE::Cell& old)
  : DRT::MESHFREE::MeshfreeBin(old),
    knotid_(old.knotid_),
    knot_(old.knot_)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  dtor                                                 (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::Cell::~Cell()
{
  return;
}


/*--------------------------------------------------------------------------*
 |  << operator                                                   nis Jan12 |
 *--------------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::MESHFREE::Cell& cell)
{
  cell.Print(os);
  return os;
}

/*--------------------------------------------------------------------------*
 |  print element                                        (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::Cell::Print(ostream& os) const
{
  DRT::Element::Print(os);
  const int nknot = NumKnot();
  const int* knotids = KnotIds();
  if (nknot > 0)
  {
    os << " Knots ";
    for (int i=0; i<nknot; ++i) os << std::setw(10) << knotids[i] << " ";
  }
  return;
}

/*--------------------------------------------------------------------------*
 |  set knot numbers to element                          (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::Cell::SetKnotIds(const int nknot, const int* knots)
{
  knotid_.resize(nknot);
  for (int i=0; i<nknot; ++i) knotid_[i] = knots[i];
  knot_.resize(0);
  return;
}

/*--------------------------------------------------------------------------*
 |  set knot numbers to element                          (public) nis Jan12 |
 *--------------------------------------------------------------------------*/

void DRT::MESHFREE::Cell::SetKnotIds(const std::string& distype, DRT::INPUT::LineDefinition* linedef)
{
  linedef->ExtractIntVector(distype,knotid_);
  for (unsigned i=0; i<knotid_.size(); ++i)
    knotid_[i] -= 1;
  knot_.resize(0);
}

/*----------------------------------------------------------------------*
 |  Pack data  (public)                                       nis Jan12 |
 *----------------------------------------------------------------------*/
void DRT::MESHFREE::Cell::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class DRT::MESHFREE::MeshfreeBin
  DRT::MESHFREE::MeshfreeBin::Pack(data);
  // add vector knotid_
  AddtoPack(data,knotid_);
  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data  (public)                                     nis Jan12 |
 *----------------------------------------------------------------------*/
void DRT::MESHFREE::Cell::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class DRT::MESHFREE::MeshfreeBin
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::MESHFREE::MeshfreeBin::Unpack(basedata);
  // extract knotid_
  ExtractfromPack(position,data,knotid_);
  // knot_ is NOT communicated
  knot_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  Build knot pointers (protected)                           nis Jan12 |
 *----------------------------------------------------------------------*/
bool DRT::MESHFREE::Cell::BuildKnotPointers(std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >& knots)
{
  int        nknot   = NumKnot();
  const int* knotids = KnotIds();
  knot_.resize(nknot);
  for (int i=0; i<nknot; ++i)
  {
    std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >::const_iterator curr = knots.find(knotids[i]);
    // this knot is not on this proc
    if (curr==knots.end()) dserror("Meshfree cell %d cannot find knot %d",Id(),knotids[i]);
    else
      knot_[i] = curr->second.get();
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  Build knot pointers (protected)                           nis Jan12 |
 *----------------------------------------------------------------------*/
bool DRT::MESHFREE::Cell::BuildKnotPointers(DRT::MESHFREE::MeshfreeNode** knots)
{
  knot_.resize(NumKnot());
  for (int i=0; i<NumKnot(); ++i) knot_[i] = knots[i];
  return true;
}

//DRT::MESHFREE::Cell::MeshfreeDiscretizationType DRT::MESHFREE::StringToDistype(std::string name)
//{
//  static std::map<std::string,DRT::MESHFREE::Cell::MeshfreeDiscretizationType> gid2distype;
//  if (gid2distype.size()==0)
//  {
//    gid2distype["HEX8"]     = DRT::MESHFREE::Cell::hex8;
//    gid2distype["TET4"]     = DRT::MESHFREE::Cell::tet4;
//    gid2distype["QUAD4"]    = DRT::MESHFREE::Cell::quad4;
//    gid2distype["TRI3"]     = DRT::MESHFREE::Cell::tri3;
//    gid2distype["LINE2"]    = DRT::MESHFREE::Cell::line2;
//    gid2distype["POINT1"]   = DRT::MESHFREE::Cell::point1;
//    gid2distype["DIS_NONE"] = DRT::MESHFREE::Cell::dis_none;
//    gid2distype["MAX_DISTYPE"] = DRT::MESHFREE::Cell::max_distype;
//  }
//
//  std::map<std::string,DRT::MESHFREE::Cell::MeshfreeDiscretizationType>::iterator i;
//  i = gid2distype.find(name);
//  if (i!=gid2distype.end())
//    return i->second;
//  dserror("unsupported distype '%s'",name.c_str());
//  return DRT::MESHFREE::Cell::dis_none;
//}
