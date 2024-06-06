/*----------------------------------------------------------------------*/
/*! \file

 \brief implementation of the 3D solid-poro element (p1, mixed approach)

 \level 2

 *----------------------------------------------------------------------*/

#include "4C_so3_poro_p1.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_poro_p1_eletypes.hpp"
#include "4C_so3_surface.hpp"

FOUR_C_NAMESPACE_OPEN

template <class so3_ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::So3PoroP1(int id, int owner)
    : So3Poro<so3_ele, distype>(id, owner), init_porosity_(Teuchos::null), is_init_porosity_(false)
{
}

template <class so3_ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::So3PoroP1(
    const Discret::ELEMENTS::So3PoroP1<so3_ele, distype>& old)
    : So3Poro<so3_ele, distype>(old),
      init_porosity_(old.init_porosity_),
      is_init_porosity_(old.is_init_porosity_)
{
}

template <class so3_ele, Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::So3PoroP1<so3_ele, distype>(*this);
  return newelement;
}

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::Pack(
    Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data, type);

  data.AddtoPack<int>(is_init_porosity_);

  if (is_init_porosity_)
    Core::Communication::ParObject::AddtoPack<Base::numnod_, 1>(data, *init_porosity_);

  // add base class Element
  Base::Pack(data);
}

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  is_init_porosity_ = Core::Communication::ParObject::ExtractInt(position, data);

  if (is_init_porosity_)
  {
    init_porosity_ = Teuchos::rcp(new Core::LinAlg::Matrix<Base::numnod_, 1>(true));
    Core::Communication::ParObject::ExtractfromPack<Base::numnod_, 1>(
        position, data, *init_porosity_);
  }


  // extract base class Element
  std::vector<char> basedata(0);
  Base::ExtractfromPack(position, data, basedata);
  Base::Unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
}

template <class so3_ele, Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::Surfaces()
{
  return Core::Communication::ElementBoundaryFactory<StructuralSurface, Core::Elements::Element>(
      Core::Communication::buildSurfaces, *this);
}

template <class so3_ele, Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::Lines()
{
  return Core::Communication::ElementBoundaryFactory<StructuralLine, Core::Elements::Element>(
      Core::Communication::buildLines, *this);
}

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::Print(std::ostream& os) const
{
  os << "So3_Poro_P1 ";
  os << Core::FE::CellTypeToString(distype).c_str() << " ";
  Core::Elements::Element::Print(os);
}

template <class so3_ele, Core::FE::CellType distype>
int Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::UniqueParObjectId() const
{
  switch (distype)
  {
    case Core::FE::CellType::hex8:
      return SoHex8PoroP1Type::Instance().UniqueParObjectId();
      break;
    case Core::FE::CellType::tet4:
      return SoTet4PoroP1Type::Instance().UniqueParObjectId();
      break;
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  return -1;
}

template <class so3_ele, Core::FE::CellType distype>
Core::Elements::ElementType& Discret::ELEMENTS::So3PoroP1<so3_ele, distype>::ElementType() const
{
  switch (distype)
  {
    case Core::FE::CellType::tet4:
      return SoTet4PoroP1Type::Instance();
    case Core::FE::CellType::hex8:
      return SoHex8PoroP1Type::Instance();
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  return SoHex8PoroP1Type::Instance();
}

FOUR_C_NAMESPACE_CLOSE

#include "4C_so3_poro_p1_fwd.hpp"
