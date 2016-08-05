/*!----------------------------------------------------------------------
\file so3_scatra.cpp

\brief Solid-scatra elements base class

\level 2

<pre>
   \maintainer Thon Moritz
               thon@mhpc.mw.tum.de
               089 - 289-10264
</pre>

*----------------------------------------------------------------------*/


#include "so3_scatra.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Scatra<so3_ele,distype>::So3_Scatra(int id, int owner):
so3_ele(id,owner),
intpoints_(distype)
{
  numgpt_ = intpoints_.NumPoints();
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Scatra<so3_ele,distype>::So3_Scatra(const DRT::ELEMENTS::So3_Scatra<so3_ele,distype>& old):
so3_ele(old),
intpoints_(distype)
{
  numgpt_ = intpoints_.NumPoints();
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            vuong 03/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::So3_Scatra<so3_ele,distype>::Clone() const
{
  DRT::ELEMENTS::So3_Scatra< so3_ele, distype>* newelement =
      new DRT::ELEMENTS::So3_Scatra< so3_ele, distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Scatra<so3_ele,distype>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data,type);
  // data_
  //so3_ele::AddtoPack(data,data_);

  // detJ_
  //so3_ele::AddtoPack(data,detJ_);

  // invJ_
  //const int size = (int)invJ_.size();
  //so3_ele::AddtoPack(data,size);
  //for (int i=0; i<size; ++i)
  //  so3_ele::AddtoPack(data,invJ_[i]);


  // add base class Element
  so3_ele::Pack(data);



  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Scatra<so3_ele,distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  so3_ele::ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // data_
  //vector<char> tmp(0);
  //so3_ele::ExtractfromPack(position,data,tmp);
  //data_.Unpack(tmp);

  // detJ_
  //so3_ele::ExtractfromPack(position,data,detJ_);
  // invJ_
  //int size = 0;
  //so3_ele::ExtractfromPack(position,data,size);
  //invJ_.resize(size, LINALG::Matrix<numdim_,numdim_>(true));
  //for (int i=0; i<size; ++i)
  //  so3_ele::ExtractfromPack(position,data,invJ_[i]);


  // extract base class Element
  std::vector<char> basedata(0);
  so3_ele::ExtractfromPack(position,data,basedata);

  so3_ele::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              vuong 03/12|
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Scatra<so3_ele,distype>::Print(std::ostream& os) const
{
  os << "So3_scatra ";
  os<<" Discretization type: "<<DRT::DistypeToString(distype).c_str();
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::So3_Scatra<so3_ele,distype>::ReadElement(const std::string& eletype,
                                             const std::string& eledistype,
                                             DRT::INPUT::LineDefinition* linedef)
{
  so3_ele::ReadElement(eletype,eledistype,linedef );

  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
inline DRT::Node** DRT::ELEMENTS::So3_Scatra<so3_ele,distype>::Nodes()
{
  return so3_ele::Nodes();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
inline Teuchos::RCP<MAT::Material>  DRT::ELEMENTS::So3_Scatra<so3_ele,distype>::Material() const
{
  return so3_ele::Material();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
inline int DRT::ELEMENTS::So3_Scatra<so3_ele,distype>::Id() const
{
  return so3_ele::Id();
}


template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8,DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex8fbar,DRT::Element::hex8>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_hex27,DRT::Element::hex27>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet4,DRT::Element::tet4>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_tet10,DRT::Element::tet10>;
template class DRT::ELEMENTS::So3_Scatra<DRT::ELEMENTS::So_weg6,DRT::Element::wedge6>;

