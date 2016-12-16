/*----------------------------------------------------------------------*/
/*!
 \file wall1_poro_p1.cpp

 \brief a 2D wall element for solid-part of porous medium using p1 (mixed) approach

\level 2

 \maintainer Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "wall1_poro_p1.H"
#include "wall1_poro_p1_eletypes.H"

#include "../drt_lib/drt_utils_factory.H"

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Wall1_PoroP1<distype>::Wall1_PoroP1(int id, int owner) :
DRT::ELEMENTS::Wall1_Poro<distype>(id,owner)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         vuong 07/13 |
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::Wall1_PoroP1<distype>::Wall1_PoroP1(const DRT::ELEMENTS::Wall1_PoroP1<distype>& old) :
DRT::ELEMENTS::Wall1_Poro<distype>(old)
{

  return;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::Wall1_PoroP1<distype>::Clone() const
{
  DRT::ELEMENTS::Wall1_PoroP1<distype>* newelement = new DRT::ELEMENTS::Wall1_PoroP1<distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  my::AddtoPack(data,type);

  // add base class Element
  my::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  my::ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  my::ExtractfromPack(position,data,basedata);
  my::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;

  return;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element> >  DRT::ELEMENTS::Wall1_PoroP1<distype>::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Wall1Line,Wall1_PoroP1>(DRT::UTILS::buildLines,this);
}


/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element> >  DRT::ELEMENTS::Wall1_PoroP1<distype>::Surfaces()
{
  std::vector<Teuchos::RCP<Element> > surfaces(1);
  surfaces[0]= Teuchos::rcp(this, false);
  return surfaces;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Wall1_PoroP1<distype>::Print(std::ostream& os) const
{
  os << "Wall1_PoroP1 ";
  Element::Print(os);
  std::cout << std::endl;
  std::cout << my::data_;
  return;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::Wall1_PoroP1<distype>::UniqueParObjectId() const
{
  switch(distype)
  {
  case DRT::Element::tri3:
    return DRT::ELEMENTS::WallTri3PoroP1Type::Instance().UniqueParObjectId();
    break;
  case DRT::Element::quad4:
    return DRT::ELEMENTS::WallQuad4PoroP1Type::Instance().UniqueParObjectId();
    break;
  case DRT::Element::quad9:
    return DRT::ELEMENTS::WallQuad9PoroP1Type::Instance().UniqueParObjectId();
    break;
  default:
    dserror("unknown element type");
    break;
  }
  return -1;
}

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ElementType & DRT::ELEMENTS::Wall1_PoroP1<distype>::ElementType() const
{
  {
    switch(distype)
    {
    case DRT::Element::tri3:
      return DRT::ELEMENTS::WallTri3PoroP1Type::Instance();
      break;
    case DRT::Element::quad4:
      return DRT::ELEMENTS::WallQuad4PoroP1Type::Instance();
      break;
    case DRT::Element::quad9:
      return DRT::ELEMENTS::WallQuad9PoroP1Type::Instance();
      break;
    default:
      dserror("unknown element type");
      break;
    }
    return DRT::ELEMENTS::WallQuad4PoroP1Type::Instance();
  }
}

/*----------------------------------------------------------------------*
 *                                                            vuong 07/13|
 *----------------------------------------------------------------------*/
template class DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::tri3>;
template class DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad4>;
template class DRT::ELEMENTS::Wall1_PoroP1<DRT::Element::quad9>;
