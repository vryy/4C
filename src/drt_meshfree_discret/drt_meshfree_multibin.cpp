/*!--------------------------------------------------------------------------
\file drt_meshfree_multibin.cpp
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

#include "drt_meshfree_multibin.H"

/// class MeshfreeMultiBinType
DRT::MESHFREE::MeshfreeMultiBinType DRT::MESHFREE::MeshfreeMultiBinType::instance_;

DRT::ParObject* DRT::MESHFREE::MeshfreeMultiBinType::Create( const std::vector<char> & data )
{
  DRT::MESHFREE::MeshfreeMultiBin* object =
    new DRT::MESHFREE::MeshfreeMultiBin(-1,-1);
  object->Unpack(data);
  return object;
}

Teuchos::RCP<DRT::Element> DRT::MESHFREE::MeshfreeMultiBinType::Create( const std::string eletype,
                                                                         const std::string eledistype,
                                                                         const int id,
                                                                         const int owner )
{
  if (eletype=="MESHFREEMULTIBIN")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::MESHFREE::MeshfreeMultiBin(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::MESHFREE::MeshfreeMultiBinType::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::MESHFREE::MeshfreeMultiBin(id,owner));
  return ele;
}



/*--------------------------------------------------------------------------*
 |  ctor                                               (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeMultiBin::MeshfreeMultiBin(int id, int owner)
  : DRT::MESHFREE::MeshfreeBin(id,owner)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  copy-ctor                                          (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeMultiBin::MeshfreeMultiBin(const DRT::MESHFREE::MeshfreeMultiBin& old)
  : DRT::MESHFREE::MeshfreeBin(old),
    associatedwalleleid_(old.associatedwalleleid_),
    associatedwallele_(old.associatedwallele_),
    associatedfluideleid_(old.associatedfluideleid_),
    associatedfluidele_(old.associatedfluidele_)
{
  return;
}

/*--------------------------------------------------------------------------*
 |  dtor                                               (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
DRT::MESHFREE::MeshfreeMultiBin::~MeshfreeMultiBin()
{
  return;
}

/*--------------------------------------------------------------------------*
 |  clone-ctor (public)                                          ghamm 04/13|
 *--------------------------------------------------------------------------*/
DRT::Element* DRT::MESHFREE::MeshfreeMultiBin::Clone() const
{
  DRT::MESHFREE::MeshfreeMultiBin* newele = new DRT::MESHFREE::MeshfreeMultiBin(*this);
  return newele;
}

/*--------------------------------------------------------------------------*
 |  << operator                                                 ghamm 04/13 |
 *--------------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::MESHFREE::MeshfreeMultiBin& bin)
{
  bin.Print(os);
  return os;
}

/*--------------------------------------------------------------------------*
 |  print element                                      (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::Print(ostream& os) const
{
  os << "MeshfreeMultiBin ";
  DRT::Element::Print(os);

  const int nwallele = NumAssociatedWallEle();
  const int* walleleids = AssociatedWallEleIds();
  if (nwallele > 0)
  {
    os << " Associated wall elements ";
    for (int j=0; j<nwallele; ++j) os << std::setw(10) << walleleids[j] << " ";
  }

  const int nfluidele = NumAssociatedFluidEle();
  const int* wfluideleids = AssociatedFluidEleIds();
  if (nfluidele > 0)
  {
    os << " Associated fluid elements ";
    for (int j=0; j<nfluidele; ++j) os << std::setw(10) << wfluideleids[j] << " ";
  }

  return;
}

/*--------------------------------------------------------------------------*
 | Delete a single wall element from the bin           (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::DeleteAssociatedWallEle(int gid)
{
  for (unsigned int i = 0; i<associatedwalleleid_.size(); i++){
    if (associatedwalleleid_[i]==gid){
      associatedwalleleid_.erase(associatedwalleleid_.begin()+i);
      associatedwallele_.erase(associatedwallele_.begin()+i);
      return;
    }
  }
  dserror("Connectivity issues: No wall element with specified gid to delete in bin. ");
  return;
}

/*--------------------------------------------------------------------------*
 | Delete a single fluid element from the bin          (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::DeleteAssociatedFluidEle(int gid)
{
  for (unsigned int i = 0; i<associatedfluideleid_.size(); i++){
    if (associatedfluideleid_[i]==gid){
      associatedfluideleid_.erase(associatedfluideleid_.begin()+i);
      associatedfluidele_.erase(associatedfluidele_.begin()+i);
      return;
    }
  }
  dserror("Connectivity issues: No fluid element with specified gid to delete in bin. ");
  return;
}

/*----------------------------------------------------------------------*
 |  Build wall element pointers                    (public) ghamm 04/13 |
 *----------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeMultiBin::BuildWallElePointers(DRT::Element** walleles)
{
  associatedwallele_.resize(NumAssociatedWallEle());
  for (int i=0; i<NumAssociatedWallEle(); ++i) associatedwallele_[i] = walleles[i];
  return true;
}

/*----------------------------------------------------------------------*
 |  Build fluid element pointers                   (public) ghamm 04/13 |
 *----------------------------------------------------------------------*/
bool DRT::MESHFREE::MeshfreeMultiBin::BuildFluidElePointers(DRT::Element** fluideles)
{
  associatedfluidele_.resize(NumAssociatedFluidEle());
  for (int i=0; i<NumAssociatedFluidEle(); ++i) associatedfluidele_[i] = fluideles[i];
  return true;
}

/*--------------------------------------------------------------------------*
 | Pack data                                           (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class DRT::Element
  DRT::Element::Pack(data);
  // add vector associatedwalleleid_
  AddtoPack(data,associatedwalleleid_);
  // add vector associatedfluideleid_
  AddtoPack(data,associatedfluideleid_);
  return;
}

/*--------------------------------------------------------------------------*
 | Unpack data                                         (public) ghamm 04/13 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeMultiBin::Unpack(const std::vector<char>& data)
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
  // extract associatedwalleleid_
  ExtractfromPack(position,data,associatedwalleleid_);
  // extract associatedfluideleid_
  ExtractfromPack(position,data,associatedfluideleid_);
  // associatedwallele_ is NOT communicated
  associatedwallele_.clear();
  // associatedwallele_ is NOT communicated
  associatedfluidele_.clear();
  return;
}
