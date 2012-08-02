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

//for ReadElement()
// TODO
//#include "../drt_mat/structporo.H"
//for secondDerivativesZero
#include "../drt_fem_general/drt_utils_shapefunctions_service.H"


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
  ishigherorder_(old.ishigherorder_)
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
  DRT::ELEMENTS::So3_Thermo< so3_ele, distype>* newelement =
      new DRT::ELEMENTS::So3_Thermo< so3_ele, distype>(*this);

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

  // TODO 2012-07-26 what do we have to pack???

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data,type);
  // data_
  so3_ele::AddtoPack(data,data_);

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
  // TODO 2012-07-26 what has to be packed/unpacked???

  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  so3_ele::ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // data_
  vector<char> tmp(0);
  so3_ele::ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  // detJ_
  so3_ele::ExtractfromPack(position,data,detJ_);
  // invJ_
  int size = 0;
  so3_ele::ExtractfromPack(position,data,size);
  invJ_.resize(size, LINALG::Matrix<numdim_,numdim_>(true));
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
  so3_ele::ReadElement(eletype,eledistype,linedef );

  Teuchos::RCP<MAT::Material> mat = so3_ele::Material();

  // TODO 2012-07-26 replace by thermo-materials
//  if(mat->MaterialType() == INPAR::MAT::m_structporo)
//  {
//    MAT::StructPoro* actmat = static_cast<MAT::StructPoro*>(mat.get());
//    if(actmat == NULL)
//      dserror("StructPoro Material Type expected for porous media!");
//    actmat->Setup(numgpt_);
//  }

  return true;
}  // ReadElement()


/*----------------------------------------------------------------------*
 | get the nodes (public)                                    dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
inline DRT::Node** DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::Nodes()
{
  return so3_ele::Nodes();
}


/*----------------------------------------------------------------------*
 | get the material (public)                                 dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
inline RCP<MAT::Material>  DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::Material(
  ) const
{
  return so3_ele::Material();
}


/*----------------------------------------------------------------------*
 | get the nodes (public)                                    dano 08/12 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
inline int DRT::ELEMENTS::So3_Thermo<so3_ele,distype>::Id() const
{
  return so3_ele::Id();
}


/*----------------------------------------------------------------------*/
// TODO check if include has to placed at the end of the file 2012-07-26
#include "so3_thermo_fwd.hpp"
