/*----------------------------------------------------------------------*/
/*!
 \file wall1_poro.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/

#include "wall1_poro.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"

//for secondDerivativesZero
#include "../drt_fem_general/drt_utils_shapefunctions_service.H"

#include "../drt_mat/structporo.H"

#include "../drt_fem_general/drt_utils_gausspoints.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Wall1_Poro<distype>::Wall1_Poro(int id, int owner) :
DRT::ELEMENTS::Wall1(id,owner),
data_(),
intpoints_(distype)
{
  numgpt_ = intpoints_.NumPoints();
  ishigherorder_ = DRT::UTILS::secondDerivativesZero<distype>();

  invJ_.resize(numgpt_, LINALG::Matrix<numdim_,numdim_>(true));
  detJ_.resize(numgpt_, 0.0);
  xsi_.resize(numgpt_, LINALG::Matrix<numdim_,1>(true));

  init_=false;

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        popp 07/10|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Wall1_Poro<distype>::Wall1_Poro(const DRT::ELEMENTS::Wall1_Poro<distype>& old) :
DRT::ELEMENTS::Wall1(old),
invJ_(old.invJ_),
detJ_(old.detJ_),
data_(old.data_),
xsi_(old.xsi_),
intpoints_(distype),
ishigherorder_(old.ishigherorder_),
init_(old.init_)
{
  numgpt_ = intpoints_.NumPoints();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::Wall1_Poro<distype>::Clone() const
{
  DRT::ELEMENTS::Wall1_Poro<distype>* newelement = new DRT::ELEMENTS::Wall1_Poro<distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::Pack(DRT::PackBuffer& data) const
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
  int size = (int)invJ_.size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
    AddtoPack(data,invJ_[i]);

  // xsi_
  size = (int)xsi_.size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
    AddtoPack(data,xsi_[i]);

  // add base class Element
  DRT::ELEMENTS::Wall1::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // data_
  std::vector<char> tmp(0);
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

  // xsi_
  size = 0;
  ExtractfromPack(position,data,size);
  xsi_.resize(size, LINALG::Matrix<numdim_,1>(true));
  for (int i=0; i<size; ++i)
    ExtractfromPack(position,data,xsi_[i]);

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  DRT::ELEMENTS::Wall1::Unpack(basedata);


  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  init_=true;

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::Print(ostream& os) const
{
  os << "Wall1_Poro ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::Wall1_Poro<distype>::ReadElement(const std::string& eletype,
                                             const std::string& eledistype,
                                             DRT::INPUT::LineDefinition* linedef)
{
  Wall1::ReadElement(eletype,eledistype,linedef );

  Teuchos::RCP<MAT::Material> mat = Wall1::Material();

  if(mat->MaterialType() == INPAR::MAT::m_structporo)
  {
    MAT::StructPoro* actmat = static_cast<MAT::StructPoro*>(mat.get());
    if(actmat == NULL)
      dserror("StructPoro Material Type expected for porous media!");
    actmat->Setup(numgpt_);
  }

  return true ;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::VisNames(std::map<string,int>& names)
{

  if (Material()->MaterialType() == INPAR::MAT::m_structporo)
  {
    string porosity = "porosity";
    names[porosity] = 1; // scalar
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::Wall1_Poro<distype>::VisData(const string& name, std::vector<double>& data)
{
  Wall1::VisData(name, data);

  if (Material()->MaterialType() == INPAR::MAT::m_structporo)
  {
    MAT::StructPoro* structporo = static_cast <MAT::StructPoro*>(Material().get());
      if (name=="porosity")
      {
        if ((int)data.size()!=1) dserror("size mismatch");
          data[0] = structporo->PorosityAv();
      }
  }
  return true;
}

template class DRT::ELEMENTS::Wall1_Poro<DRT::Element::quad4>;
template class DRT::ELEMENTS::Wall1_Poro<DRT::Element::quad9>;
