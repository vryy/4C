/*----------------------------------------------------------------------*/
/*! \file

 \brief implementation of the 3D solid-poro element (p1, mixed approach)

 \level 2

 *----------------------------------------------------------------------*/

#include "baci_so3_poro_p1.hpp"

#include "baci_comm_utils_factory.hpp"
#include "baci_so3_line.hpp"
#include "baci_so3_poro_p1_eletypes.hpp"
#include "baci_so3_surface.hpp"

BACI_NAMESPACE_OPEN

template <class so3_ele, CORE::FE::CellType distype>
DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::So3_Poro_P1(int id, int owner)
    : So3_Poro<so3_ele, distype>(id, owner), init_porosity_(Teuchos::null), is_init_porosity_(false)
{
}

template <class so3_ele, CORE::FE::CellType distype>
DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::So3_Poro_P1(
    const DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>& old)
    : So3_Poro<so3_ele, distype>(old),
      init_porosity_(old.init_porosity_),
      is_init_porosity_(old.is_init_porosity_)
{
}

template <class so3_ele, CORE::FE::CellType distype>
DRT::Element* DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>(*this);
  return newelement;
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data, type);

  data.AddtoPack<int>(is_init_porosity_);

  if (is_init_porosity_) CORE::COMM::ParObject::AddtoPack<Base::numnod_, 1>(data, *init_porosity_);

  // add base class Element
  Base::Pack(data);
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  is_init_porosity_ = CORE::COMM::ParObject::ExtractInt(position, data);

  if (is_init_porosity_)
  {
    init_porosity_ = Teuchos::rcp(new CORE::LINALG::Matrix<Base::numnod_, 1>(true));
    CORE::COMM::ParObject::ExtractfromPack<Base::numnod_, 1>(position, data, *init_porosity_);
  }


  // extract base class Element
  std::vector<char> basedata(0);
  Base::ExtractfromPack(position, data, basedata);
  Base::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
}

template <class so3_ele, CORE::FE::CellType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Surfaces()
{
  return CORE::COMM::ElementBoundaryFactory<StructuralSurface, DRT::Element>(
      CORE::COMM::buildSurfaces, *this);
}

template <class so3_ele, CORE::FE::CellType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<StructuralLine, DRT::Element>(
      CORE::COMM::buildLines, *this);
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::Print(std::ostream& os) const
{
  os << "So3_Poro_P1 ";
  os << CORE::FE::CellTypeToString(distype).c_str() << " ";
  Element::Print(os);
}

template <class so3_ele, CORE::FE::CellType distype>
int DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::UniqueParObjectId() const
{
  switch (distype)
  {
    case CORE::FE::CellType::hex8:
      return So_hex8PoroP1Type::Instance().UniqueParObjectId();
      break;
    case CORE::FE::CellType::tet4:
      return So_tet4PoroP1Type::Instance().UniqueParObjectId();
      break;
    default:
      dserror("unknown element type!");
      break;
  }
  return -1;
}

template <class so3_ele, CORE::FE::CellType distype>
DRT::ElementType& DRT::ELEMENTS::So3_Poro_P1<so3_ele, distype>::ElementType() const
{
  switch (distype)
  {
    case CORE::FE::CellType::tet4:
      return So_tet4PoroP1Type::Instance();
    case CORE::FE::CellType::hex8:
      return So_hex8PoroP1Type::Instance();
    default:
      dserror("unknown element type!");
      break;
  }
  return So_hex8PoroP1Type::Instance();
}

BACI_NAMESPACE_CLOSE

#include "baci_so3_poro_p1_fwd.hpp"
