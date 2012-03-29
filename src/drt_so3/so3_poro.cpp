/*!----------------------------------------------------------------------
\file so_hex8_poro.cpp
\brief

<pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
</pre>

*----------------------------------------------------------------------*/

#include "so3_poro.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Poro<distype>::So3_Poro(int id, int owner):
DRT::Element(id,owner),
intpoints_(distype)
{
  numgpt_ = intpoints_.NumPoints();
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Poro<distype>::So3_Poro(const DRT::ELEMENTS::So3_Poro<distype>& old):
DRT::Element(old),
intpoints_(distype)
{
  numgpt_ = intpoints_.NumPoints();
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro<distype>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // data_
  AddtoPack(data,data_);

  // detJ_
  AddtoPack(data,detJ_);

  // invJ_
  const int size = (int)invJ_.size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
    AddtoPack(data,invJ_[i]);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro<distype>::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  // detJ_
  ExtractfromPack(position,data,detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position,data,size);
  invJ_.resize(size, LINALG::Matrix<numdim_,numdim_>(true));
  for (int i=0; i<size; ++i)
    ExtractfromPack(position,data,invJ_[i]);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              vuong 03/12|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro<distype>::Print(ostream& os) const
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
  template<DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::So3_Poro<distype>::ReadElement(const std::string& eletype,
                                             const std::string& eledistype,
                                             DRT::INPUT::LineDefinition* linedef)
{
  return true;
}

#include "so3_poro_fwd.hpp"

