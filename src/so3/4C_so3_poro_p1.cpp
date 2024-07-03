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

template <class So3Ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::So3PoroP1(int id, int owner)
    : So3Poro<So3Ele, distype>(id, owner), init_porosity_(Teuchos::null), is_init_porosity_(false)
{
}

template <class So3Ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::So3PoroP1(
    const Discret::ELEMENTS::So3PoroP1<So3Ele, distype>& old)
    : So3Poro<So3Ele, distype>(old),
      init_porosity_(old.init_porosity_),
      is_init_porosity_(old.is_init_porosity_)
{
}

template <class So3Ele, Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::So3PoroP1<So3Ele, distype>(*this);
  return newelement;
}

template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::pack(
    Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  So3Ele::add_to_pack(data, type);

  data.add_to_pack<int>(is_init_porosity_);

  if (is_init_porosity_)
    Core::Communication::ParObject::add_to_pack<Base::numnod_, 1>(data, *init_porosity_);

  // add base class Element
  Base::pack(data);
}

template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  is_init_porosity_ = Core::Communication::ParObject::extract_int(position, data);

  if (is_init_porosity_)
  {
    init_porosity_ = Teuchos::rcp(new Core::LinAlg::Matrix<Base::numnod_, 1>(true));
    Core::Communication::ParObject::extract_from_pack<Base::numnod_, 1>(
        position, data, *init_porosity_);
  }


  // extract base class Element
  std::vector<char> basedata(0);
  Base::extract_from_pack(position, data, basedata);
  Base::unpack(basedata);

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
}

template <class So3Ele, Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::Surfaces()
{
  return Core::Communication::ElementBoundaryFactory<StructuralSurface, Core::Elements::Element>(
      Core::Communication::buildSurfaces, *this);
}

template <class So3Ele, Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::Lines()
{
  return Core::Communication::ElementBoundaryFactory<StructuralLine, Core::Elements::Element>(
      Core::Communication::buildLines, *this);
}

template <class So3Ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::print(std::ostream& os) const
{
  os << "So3_Poro_P1 ";
  os << Core::FE::CellTypeToString(distype).c_str() << " ";
  Core::Elements::Element::print(os);
}

template <class So3Ele, Core::FE::CellType distype>
int Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::UniqueParObjectId() const
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

template <class So3Ele, Core::FE::CellType distype>
Core::Elements::ElementType& Discret::ELEMENTS::So3PoroP1<So3Ele, distype>::ElementType() const
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
