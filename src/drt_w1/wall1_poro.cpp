/*----------------------------------------------------------------------*/
/*!
 \file wall1_poro.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

#include "wall1_poro.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"

//for secondDerivativesZero
#include "../drt_fem_general/drt_utils_shapefunctions_service.H"

#include "../drt_mat/structporo.H"
#include "../drt_mat/structporo_reaction.H"
#include "../drt_mat/fluidporo.H"

#include "../drt_fem_general/drt_utils_gausspoints.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Wall1_Poro<distype>::Wall1_Poro(int id, int owner) :
DRT::ELEMENTS::Wall1(id,owner),
data_(),
intpoints_(distype),
weights_(true),
myknots_(numdim_)
{
  numgpt_ = intpoints_.NumPoints();
  ishigherorder_ = DRT::UTILS::secondDerivativesZero<distype>();

  invJ_.resize(numgpt_, LINALG::Matrix<numdim_,numdim_>(true));
  detJ_.resize(numgpt_, 0.0);
  xsi_.resize(numgpt_, LINALG::Matrix<numdim_,1>(true));

  init_=false;

  scatracoupling_=false;

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
init_(old.init_),
scatracoupling_(old.scatracoupling_),
weights_(old.weights_),
myknots_(old.myknots_)
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

  // scatracoupling_
  AddtoPack(data,scatracoupling_);

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

  // scatracoupling_
  scatracoupling_ = (bool)( ExtractInt(position,data) );

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
std::vector<Teuchos::RCP<DRT::Element> >  DRT::ELEMENTS::Wall1_Poro<distype>::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Wall1Line,Wall1_Poro>(DRT::UTILS::buildLines,this);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element> >  DRT::ELEMENTS::Wall1_Poro<distype>::Surfaces()
{
  std::vector<Teuchos::RCP<Element> > surfaces(1);
  surfaces[0]= Teuchos::rcp(this, false);
  return surfaces;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::Print(std::ostream& os) const
{
  os << "Wall1_Poro ";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << data_;
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::Wall1_Poro<distype>::ReadElement(const std::string& eletype,
                                             const std::string& eledistype,
                                             DRT::INPUT::LineDefinition* linedef)
{
  // read base element
  Wall1::ReadElement(eletype,eledistype,linedef );

  //setup poro material
  Teuchos::RCP<MAT::StructPoro> poromat = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(Material());
  if(poromat==Teuchos::null)
    dserror("material assigned to poro element is not a poro material!");
  poromat->PoroSetup(numgpt_, linedef);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::GetMaterials( )
{
  if(fluidmat_==Teuchos::null)
  {
    // get global id of the structure element
    int id = Id();
    //access fluid discretization
    Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis("porofluid");
    if(fluiddis==Teuchos::null)
      dserror("no fluid discretization");

    //get corresponding fluid element
    DRT::Element* fluidele = fluiddis->gElement(id);
    if (fluidele == NULL)
      dserror("Fluid element %i not on local processor", id);

    //get fluid material
    fluidmat_ = Teuchos::rcp_dynamic_cast<MAT::FluidPoro>(fluidele->Material());
    if(fluidmat_->MaterialType() != INPAR::MAT::m_fluidporo)
      dserror("invalid fluid material for poroelasticity");
  }

  if(structmat_==Teuchos::null)
  {
    //get structure material
    if(Material()==Teuchos::null) dserror("no structure material available!");
    structmat_ = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(Material());
    if(structmat_->MaterialType() != INPAR::MAT::m_structporo and
       structmat_->MaterialType() != INPAR::MAT::m_structpororeaction)
      dserror("invalid structure material for poroelasticity");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)               |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_Poro<distype>::VisNames(std::map<std::string,int>& names)
{
  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
  so3mat->VisNames(names);

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                         |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::Wall1_Poro<distype>::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name,data))
    return true;

  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
  return so3mat->VisData(name, data, numgpt_,this->Id());
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template class DRT::ELEMENTS::Wall1_Poro<DRT::Element::quad4>;
template class DRT::ELEMENTS::Wall1_Poro<DRT::Element::quad9>;
template class DRT::ELEMENTS::Wall1_Poro<DRT::Element::nurbs4>;
template class DRT::ELEMENTS::Wall1_Poro<DRT::Element::nurbs9>;
