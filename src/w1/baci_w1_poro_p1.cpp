/*----------------------------------------------------------------------------*/
/*! \file
\brief A 2D wall element for solid-part of porous medium using p1 (mixed) approach.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "baci_w1_poro_p1.hpp"

#include "baci_comm_utils_factory.hpp"
#include "baci_w1_poro_p1_eletypes.hpp"

FOUR_C_NAMESPACE_OPEN

template <CORE::FE::CellType distype>
DRT::ELEMENTS::Wall1PoroP1<distype>::Wall1PoroP1(int id, int owner)
    : DRT::ELEMENTS::Wall1Poro<distype>(id, owner)
{
}

template <CORE::FE::CellType distype>
DRT::ELEMENTS::Wall1PoroP1<distype>::Wall1PoroP1(const DRT::ELEMENTS::Wall1PoroP1<distype>& old)
    : DRT::ELEMENTS::Wall1Poro<distype>(old)
{
}

template <CORE::FE::CellType distype>
DRT::Element* DRT::ELEMENTS::Wall1PoroP1<distype>::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::Wall1PoroP1<distype>(*this);
  return newelement;
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Wall1PoroP1<distype>::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  Base::AddtoPack(data, type);

  // add base class Element
  Base::Pack(data);
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Wall1PoroP1<distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  Base::ExtractfromPack(position, data, basedata);
  Base::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
}

template <CORE::FE::CellType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Wall1PoroP1<distype>::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<Wall1Line, Wall1PoroP1>(CORE::COMM::buildLines, *this);
}

template <CORE::FE::CellType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::Wall1PoroP1<distype>::Surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}

template <CORE::FE::CellType distype>
void DRT::ELEMENTS::Wall1PoroP1<distype>::Print(std::ostream& os) const
{
  os << "Wall1_PoroP1 ";
  Element::Print(os);
  std::cout << std::endl;
}

template <CORE::FE::CellType distype>
int DRT::ELEMENTS::Wall1PoroP1<distype>::UniqueParObjectId() const
{
  switch (distype)
  {
    case CORE::FE::CellType::tri3:
      return DRT::ELEMENTS::WallTri3PoroP1Type::Instance().UniqueParObjectId();
    case CORE::FE::CellType::quad4:
      return DRT::ELEMENTS::WallQuad4PoroP1Type::Instance().UniqueParObjectId();
    case CORE::FE::CellType::quad9:
      return DRT::ELEMENTS::WallQuad9PoroP1Type::Instance().UniqueParObjectId();
    default:
      dserror("unknown element type");
  }
  return -1;
}

template <CORE::FE::CellType distype>
DRT::ElementType& DRT::ELEMENTS::Wall1PoroP1<distype>::ElementType() const
{
  {
    switch (distype)
    {
      case CORE::FE::CellType::tri3:
        return DRT::ELEMENTS::WallTri3PoroP1Type::Instance();
      case CORE::FE::CellType::quad4:
        return DRT::ELEMENTS::WallQuad4PoroP1Type::Instance();
      case CORE::FE::CellType::quad9:
        return DRT::ELEMENTS::WallQuad9PoroP1Type::Instance();
      default:
        dserror("unknown element type");
    }
    return DRT::ELEMENTS::WallQuad4PoroP1Type::Instance();
  }
}

template class DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::quad4>;
template class DRT::ELEMENTS::Wall1PoroP1<CORE::FE::CellType::quad9>;

FOUR_C_NAMESPACE_CLOSE
