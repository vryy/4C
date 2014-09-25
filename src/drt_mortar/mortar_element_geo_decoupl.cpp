/*!----------------------------------------------------------------------
\file mortar_element_geo_decoupl.cpp
\brief A mortar coupling element with decoupled geometric (=point) and dof (=node) information

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#include "mortar_element_geo_decoupl.H"

MORTAR::MortarElementGeoDecouplType MORTAR::MortarElementGeoDecouplType::instance_;


DRT::ParObject* MORTAR::MortarElementGeoDecouplType::Create( const std::vector<char> & data )
{
  MORTAR::MortarElementGeoDecoupl* ele = new MORTAR::MortarElementGeoDecoupl(0,
                                                         0,DRT::Element::dis_none,
                                                         0,0,0,NULL,false);
  ele->Unpack(data);
  return ele;
}


Teuchos::RCP<DRT::Element> MORTAR::MortarElementGeoDecouplType::Create( const int id, const int owner )
{
  //return Teuchos::rcp( new MortarElementGeoDecoupl( id, owner ) );
  return Teuchos::null;
}


void MORTAR::MortarElementGeoDecouplType::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
}

void MORTAR::MortarElementGeoDecouplType::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                            ghamm 09/14|
 *----------------------------------------------------------------------*/
MORTAR::MortarElementGeoDecoupl::MortarElementGeoDecoupl(int id, int owner,
                           const DRT::Element::DiscretizationType& shape,
                           const int numnode,
                           const int* nodeids,
                           const int numpoint,
                           const int* pointids,
                           const bool isslave,
                           bool isnurbs) :
DRT::Element(id,owner),  // necessary due to virtual inheritance from DRT::Element
MORTAR::MortarElement(id,owner,shape,numnode,nodeids,isslave,isnurbs),
DRT::MESHFREE::Cell(id,owner)
{
  SetPointIds(numpoint,pointids);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       ghamm 09/14|
 *----------------------------------------------------------------------*/
MORTAR::MortarElementGeoDecoupl::MortarElementGeoDecoupl(const MORTAR::MortarElementGeoDecoupl& old) :
DRT::Element(old),  // necessary due to virtual inheritance from DRT::Element
MORTAR::MortarElement(old),
DRT::MESHFREE::Cell(old)
{
  // not yet used and thus not necessarily consistent
  dserror("ERROR: MortarElement copy-ctor not yet implemented");

  return;
}

/*----------------------------------------------------------------------*
 |  clone-ctor (public)                                      ghamm 09/14|
 *----------------------------------------------------------------------*/
DRT::Element* MORTAR::MortarElementGeoDecoupl::Clone() const
{
  MORTAR::MortarElementGeoDecoupl* newele = new MORTAR::MortarElementGeoDecoupl(*this);
  return newele;
}

/*----------------------------------------------------------------------*
 |  << operator                                              ghamm 09/14|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const MORTAR::MortarElementGeoDecoupl& element)
{
  element.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print element (public)                                   ghamm 09/14|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElementGeoDecoupl::Print(std::ostream& os) const
{
  os << "Mortar Element Geo Decoupl ";
  DRT::Element::Print(os);

  const int npoint = DRT::MESHFREE::Cell::NumPoint();
  const int* pointids = DRT::MESHFREE::Cell::PointIds();
  if (npoint > 0)
  {
    os << " Points ";
    for (int i=0; i<npoint; ++i) os << std::setw(10) << pointids[i] << " ";
  }

  if (isslave_) os << " Slave  ";
  else          os << " Master ";

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           ghamm 09/14|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElementGeoDecoupl::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class MortarElement
  MORTAR::MortarElement::Pack(data);
  // add base class Cell
  DRT::MESHFREE::Cell::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           ghamm 09/14|
 *----------------------------------------------------------------------*/
void MORTAR::MortarElementGeoDecoupl::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class DRT::Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  // extract base class MortarElement
  MORTAR::MortarElement::Unpack(basedata);
  // extract base class Cell
  DRT::MESHFREE::Cell::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

