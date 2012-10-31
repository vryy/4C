/*!----------------------------------------------------------------------
\file so_hex8_thermo.cpp
\brief

<pre>
   Maintainer: Caroline Danowski
               danowski@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15253
</pre>

*----------------------------------------------------------------------*/

#include "so3_thermo.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 | ctor (public)                                              dano 08/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::So3_Thermo(
  int id,
  int owner
  ):
  so3_ele(id,owner),
  intpoints_(distype)
{
  numgpt_ = intpoints_.NumPoints();
  ishigherorder_ = DRT::UTILS::secondDerivativesZero<distype>();
  kintype_ = so3_thermo_nonlinear;  // TODO 2012-10-26 default is defined here!!
  return;
}


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                        dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::So3_Thermo(
  const DRT::ELEMENTS::So3_Thermo<so3_ele,distype>& old
  ):
  so3_ele(old),
  intpoints_(distype),
  ishigherorder_(old.ishigherorder_),
  kintype_(old.kintype_)
{
  numgpt_ = intpoints_.NumPoints();
  return;
}


/*----------------------------------------------------------------------*
 | deep copy this instance of Solid3 and return pointer to   dano 08/12 |
 | it (public)                                                          |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::Clone() const
{
  DRT::ELEMENTS::So3_Thermo< so3_ele, distype>* newelement
    = new DRT::ELEMENTS::So3_Thermo< so3_ele, distype>(*this);

  return newelement;
}


/*----------------------------------------------------------------------*
 | pack data (public)                                        dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::Pack(
  DRT::PackBuffer& data
  ) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data,type);
  // data_
  so3_ele::AddtoPack(data,data_);
  // kintype_
  so3_ele::AddtoPack(data,kintype_);
  // detJ_
  so3_ele::AddtoPack(data,detJ_);

  // invJ_
  const int size = (int)invJ_.size();
  so3_ele::AddtoPack(data,size);
  for (int i=0; i<size; ++i)
    so3_ele::AddtoPack(data,invJ_[i]);

  // add base class Element
  so3_ele::Pack(data);

  return;

}  // Pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                      dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::Unpack(
  const vector<char>& data
  )
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  so3_ele::ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class element data_
  vector<char> tmp(0);
  so3_ele::ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);
  // kintype_
//  so3_ele::ExtractfromPack(position,data,kintype_);
  kintype_ = static_cast<KinematicType>( so3_ele::ExtractInt(position,data) );
  // detJ_
  so3_ele::ExtractfromPack(position,data,detJ_);
  // invJ_
  int size = 0;
  so3_ele::ExtractfromPack(position,data,size);
  invJ_.resize(size, LINALG::Matrix<nsd_,nsd_>(true));
  for (int i=0; i<size; ++i)
    so3_ele::ExtractfromPack(position,data,invJ_[i]);

  // extract base class Element
  vector<char> basedata(0);
  so3_ele::ExtractfromPack(position,data,basedata);
  so3_ele::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;

}  // Unpack()


/*----------------------------------------------------------------------*
 | print this element (public)                               dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::Print(ostream& os) const
{
  os << "So3_Thermo ";
  return;
}


/*----------------------------------------------------------------------*
 | read this element, get the material (public)              dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::ReadElement(
  const std::string& eletype,
  const std::string& eledistype,
  DRT::INPUT::LineDefinition* linedef
  )
{
  so3_ele::ReadElement(eletype,eledistype,linedef);

  Teuchos::RCP<MAT::Material> mat = so3_ele::Material();

  // read thermo-materials
  switch (mat->MaterialType())
  {
  // TODO 2012-10-30 plastic materials are already set in L181: so_hex8_input ReadElement
  // materials without history need no separate Setup()
  default :
  break;
  }

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);

  // geometrically linear
  if(buffer=="linear")
    kintype_ = so3_thermo_linear;
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer=="nonlinear")
    kintype_ = so3_thermo_nonlinear;
  else
    dserror("Reading of SO3_THERMO element failed! KINEM unknown");

  return true;
}  // ReadElement()


/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                           dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
inline DRT::Node** DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::Nodes()
{
  return so3_ele::Nodes();
}


/*----------------------------------------------------------------------*
 | get the material from so3 (public)                        dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
inline RCP<MAT::Material>  DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::Material(
  ) const
{
  return so3_ele::Material();
}


/*----------------------------------------------------------------------*
 | get the node Ids from so3 (public)                        dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
inline int DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::Id() const
{
  return so3_ele::Id();
}


/*----------------------------------------------------------------------*/
// include the file at the end of so3_thermo.cpp
#include "so3_thermo_fwd.hpp"