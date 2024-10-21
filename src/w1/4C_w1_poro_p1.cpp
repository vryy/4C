#include "4C_w1_poro_p1.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_utils_factory.hpp"
#include "4C_w1_poro_p1_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

template <Core::FE::CellType distype>
Discret::ELEMENTS::Wall1PoroP1<distype>::Wall1PoroP1(int id, int owner)
    : Discret::ELEMENTS::Wall1Poro<distype>(id, owner)
{
}

template <Core::FE::CellType distype>
Discret::ELEMENTS::Wall1PoroP1<distype>::Wall1PoroP1(
    const Discret::ELEMENTS::Wall1PoroP1<distype>& old)
    : Discret::ELEMENTS::Wall1Poro<distype>(old)
{
}

template <Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::Wall1PoroP1<distype>::clone() const
{
  auto* newelement = new Discret::ELEMENTS::Wall1PoroP1<distype>(*this);
  return newelement;
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1PoroP1<distype>::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // add base class Element
  Base::pack(data);
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1PoroP1<distype>::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // extract base class Element
  std::vector<char> basedata(0);
  extract_from_pack(buffer, basedata);
  Core::Communication::UnpackBuffer basedata_buffer(basedata);
  Base::unpack(basedata_buffer);

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}

template <Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Wall1PoroP1<distype>::lines()
{
  return Core::Communication::element_boundary_factory<Wall1Line, Wall1PoroP1>(
      Core::Communication::buildLines, *this);
}

template <Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::Wall1PoroP1<distype>::surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1PoroP1<distype>::print(std::ostream& os) const
{
  os << "Wall1_PoroP1 ";
  Core::Elements::Element::print(os);
  std::cout << std::endl;
}

template <Core::FE::CellType distype>
int Discret::ELEMENTS::Wall1PoroP1<distype>::unique_par_object_id() const
{
  switch (distype)
  {
    case Core::FE::CellType::tri3:
      return Discret::ELEMENTS::WallTri3PoroP1Type::instance().unique_par_object_id();
    case Core::FE::CellType::quad4:
      return Discret::ELEMENTS::WallQuad4PoroP1Type::instance().unique_par_object_id();
    case Core::FE::CellType::quad9:
      return Discret::ELEMENTS::WallQuad9PoroP1Type::instance().unique_par_object_id();
    default:
      FOUR_C_THROW("unknown element type");
  }
  return -1;
}

template <Core::FE::CellType distype>
Core::Elements::ElementType& Discret::ELEMENTS::Wall1PoroP1<distype>::element_type() const
{
  {
    switch (distype)
    {
      case Core::FE::CellType::tri3:
        return Discret::ELEMENTS::WallTri3PoroP1Type::instance();
      case Core::FE::CellType::quad4:
        return Discret::ELEMENTS::WallQuad4PoroP1Type::instance();
      case Core::FE::CellType::quad9:
        return Discret::ELEMENTS::WallQuad9PoroP1Type::instance();
      default:
        FOUR_C_THROW("unknown element type");
    }
    return Discret::ELEMENTS::WallQuad4PoroP1Type::instance();
  }
}

template class Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::Wall1PoroP1<Core::FE::CellType::quad9>;

FOUR_C_NAMESPACE_CLOSE
