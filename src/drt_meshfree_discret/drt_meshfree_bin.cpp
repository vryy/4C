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

/// class MeshfreeTransportType
DRT::MESHFREE::MeshfreeBinType DRT::MESHFREE::MeshfreeBinType::instance_;

DRT::ParObject* DRT::MESHFREE::MeshfreeBinType::Create( const std::vector<char> & data )
{
  DRT::MESHFREE::MeshfreeBin* object =
    new DRT::MESHFREE::MeshfreeBin(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::MESHFREE::MeshfreeBinType::Create( const std::string eletype,
                                                                         const std::string eledistype,
                                                                         const int id,
                                                                         const int owner )
{
  if (eletype=="MESHFREEBIN")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::MESHFREE::MeshfreeBin(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::MESHFREE::MeshfreeBinType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::MESHFREE::MeshfreeBin(id,owner));
  return ele;
}



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
 |  clone-ctor (public)                                          ghamm 11/12|
 *--------------------------------------------------------------------------*/
DRT::Element* DRT::MESHFREE::MeshfreeBin::Clone() const
{
  DRT::MESHFREE::MeshfreeBin* newele = new DRT::MESHFREE::MeshfreeBin(*this);
  return newele;
}


/*--------------------------------------------------------------------------*
 |  << operator                                                 ghamm 11/12 |
 *--------------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::MESHFREE::MeshfreeBin& bin)
{
  bin.Print(os);
  return os;
}

/*--------------------------------------------------------------------------*
 |  print element                                      (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeBin::Print(ostream& os) const
{
  os << "MeshfreeBin ";
  DRT::Element::Print(os);
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
  dserror("Connectivity issues: No node with secified gid to delete in element. ");
  return;
}

/*--------------------------------------------------------------------------*
 | Pack data                                           (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeBin::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class DRT::Element
  DRT::Element::Pack(data);
  // nothing to pack for meshfree bin
  return;
}

/*--------------------------------------------------------------------------*
 | Unpack data                                         (public) ghamm 11/12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeBin::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class DRT::Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::Element::Unpack(basedata);
  // nothing to extract for meshfree bin
  return;
}
