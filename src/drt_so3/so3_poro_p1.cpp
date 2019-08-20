/*----------------------------------------------------------------------*/
/*! \file

 \brief implementation of the 3D solid-poro element (p1, mixed approach)

 \level 2

\maintainer  Christoph Ager
 *----------------------------------------------------------------------*/

#include "so3_poro_p1.H"
#include "so3_poro_p1_eletypes.H"

#include "so_surface.H"
#include "so_line.H"

#include "../drt_lib/drt_utils_factory.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::So3_Poro_P1(int id, int owner)
    : So3_Poro<so3_ele, distype>(id, owner), init_porosity_(Teuchos::null), is_init_porosity_(false)
{
  return;
}


/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       vuong 03/12|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::So3_Poro_P1(
    const DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>& old)
    : So3_Poro<so3_ele, distype>(old),
      init_porosity_(old.init_porosity_),
      is_init_porosity_(old.is_init_porosity_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Clone() const
{
  DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>* newelement =
      new DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data, type);

  data.AddtoPack<int>(is_init_porosity_);

  if (is_init_porosity_) DRT::ParObject::AddtoPack<my::numnod_, 1>(data, *init_porosity_);

  // add base class Element
  my::Pack(data);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  so3_ele::ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  is_init_porosity_ = DRT::ParObject::ExtractInt(position, data);

  if (is_init_porosity_)
  {
    init_porosity_ = Teuchos::rcp(new LINALG::Matrix<my::numnod_, 1>(true));
    DRT::ParObject::ExtractfromPack<my::numnod_, 1>(position, data, *init_porosity_);
  }


  // extract base class Element
  std::vector<char> basedata(0);
  my::ExtractfromPack(position, data, basedata);
  my::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                 vuong 11/13|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Volumes()
{
  std::vector<Teuchos::RCP<Element>> volumes(1);
  volumes[0] = Teuchos::rcp(this, false);
  return volumes;
}

/*----------------------------------------------------------------------*
|  get vector of surfaces (public)                           vuong 11/13|
|  surface normals always point outward                                 |
*----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface, DRT::Element>(
      DRT::UTILS::buildSurfaces, this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             vuong 11/13|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine, DRT::Element>(
      DRT::UTILS::buildLines, this);
}

/*----------------------------------------------------------------------*
 |  print this element (public)                              vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Print(std::ostream& os) const
{
  os << "So3_Poro_P1 ";
  os << DRT::DistypeToString(distype).c_str() << " ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::UniqueParObjectId() const
{
  switch (distype)
  {
    case DRT::Element::hex8:
      return So_hex8PoroP1Type::Instance().UniqueParObjectId();
      break;
    case DRT::Element::tet4:
      return So_tet4PoroP1Type::Instance().UniqueParObjectId();
      break;
    default:
      dserror("unknown element type!");
      break;
  }
  return -1;
}

/*----------------------------------------------------------------------*
 |                                                           vuong 03/12|
 *----------------------------------------------------------------------*/
template <class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ElementType& DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::ElementType() const
{
  switch (distype)
  {
    case DRT::Element::tet4:
      return So_tet4PoroP1Type::Instance();
    // case DRT::Element::tet10:
    //  return So_tet10PoroType::Instance();
    case DRT::Element::hex8:
      return So_hex8PoroP1Type::Instance();
    // case DRT::Element::hex27:
    //  return So_hex27PoroType::Instance();
    default:
      dserror("unknown element type!");
      break;
  }
  return So_hex8PoroP1Type::Instance();
}

#include "so3_poro_p1_fwd.hpp"
