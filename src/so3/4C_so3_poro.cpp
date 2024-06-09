/*----------------------------------------------------------------------*/
/*! \file
\brief implementation of the 3D solid-poro element


\level 2

*----------------------------------------------------------------------*/

#include "4C_so3_poro.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_so3_line.hpp"
#include "4C_so3_poro_eletypes.hpp"
#include "4C_so3_surface.hpp"

FOUR_C_NAMESPACE_OPEN

template <class so3_ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3Poro<so3_ele, distype>::So3Poro(int id, int owner)
    : so3_ele(id, owner),
      intpoints_(distype),
      init_(false),
      scatra_coupling_(false),
      isNurbs_(false),
      weights_(true),
      myknots_(numdim_),
      fluid_mat_(Teuchos::null),
      fluidmulti_mat_(Teuchos::null),
      struct_mat_(Teuchos::null)
{
  numgpt_ = intpoints_.NumPoints();

  invJ_.resize(numgpt_, Core::LinAlg::Matrix<numdim_, numdim_>(true));
  detJ_.resize(numgpt_, 0.0);
  xsi_.resize(numgpt_, Core::LinAlg::Matrix<numdim_, 1>(true));
  anisotropic_permeability_directions_.resize(3, std::vector<double>(3, 0.0));
  anisotropic_permeability_nodal_coeffs_.resize(3, std::vector<double>(numnod_, 0.0));
}

template <class so3_ele, Core::FE::CellType distype>
Discret::ELEMENTS::So3Poro<so3_ele, distype>::So3Poro(
    const Discret::ELEMENTS::So3Poro<so3_ele, distype>& old)
    : so3_ele(old),
      invJ_(old.invJ_),
      detJ_(old.detJ_),
      xsi_(old.xsi_),
      intpoints_(distype),
      init_(old.init_),
      scatra_coupling_(old.scatra_coupling_),
      isNurbs_(old.isNurbs_),
      weights_(old.weights_),
      myknots_(old.myknots_),
      fluid_mat_(old.fluid_mat_),
      fluidmulti_mat_(old.fluidmulti_mat_),
      struct_mat_(old.struct_mat_),
      anisotropic_permeability_directions_(old.anisotropic_permeability_directions_),
      anisotropic_permeability_nodal_coeffs_(old.anisotropic_permeability_nodal_coeffs_)
{
  numgpt_ = intpoints_.NumPoints();
}

template <class so3_ele, Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::So3Poro<so3_ele, distype>::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::So3Poro<so3_ele, distype>(*this);
  return newelement;
}

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Poro<so3_ele, distype>::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::add_to_pack(data, type);

  // detJ_
  so3_ele::add_to_pack(data, detJ_);

  // invJ_
  auto size = static_cast<int>(invJ_.size());
  so3_ele::add_to_pack(data, size);
  for (int i = 0; i < size; ++i) so3_ele::add_to_pack(data, invJ_[i]);

  // xsi_
  size = static_cast<int>(xsi_.size());
  so3_ele::add_to_pack(data, size);
  for (int i = 0; i < size; ++i) so3_ele::add_to_pack(data, xsi_[i]);

  // scatra_coupling_
  so3_ele::add_to_pack(data, scatra_coupling_);

  // isNurbs_
  so3_ele::add_to_pack(data, isNurbs_);

  // anisotropic_permeability_directions_
  size = static_cast<int>(anisotropic_permeability_directions_.size());
  so3_ele::add_to_pack(data, size);
  for (int i = 0; i < size; ++i)
    so3_ele::add_to_pack(data, anisotropic_permeability_directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = static_cast<int>(anisotropic_permeability_nodal_coeffs_.size());
  so3_ele::add_to_pack(data, size);
  for (int i = 0; i < size; ++i)
    so3_ele::add_to_pack(data, anisotropic_permeability_nodal_coeffs_[i]);

  // add base class Element
  so3_ele::Pack(data);
}

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Poro<so3_ele, distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // detJ_
  so3_ele::extract_from_pack(position, data, detJ_);

  // invJ_
  int size = 0;
  so3_ele::extract_from_pack(position, data, size);
  invJ_.resize(size, Core::LinAlg::Matrix<numdim_, numdim_>(true));
  for (int i = 0; i < size; ++i) so3_ele::extract_from_pack(position, data, invJ_[i]);

  // xsi_
  size = 0;
  so3_ele::extract_from_pack(position, data, size);
  xsi_.resize(size, Core::LinAlg::Matrix<numdim_, 1>(true));
  for (int i = 0; i < size; ++i) so3_ele::extract_from_pack(position, data, xsi_[i]);

  // scatra_coupling_
  scatra_coupling_ = static_cast<bool>(so3_ele::ExtractInt(position, data));

  // isNurbs_
  isNurbs_ = static_cast<bool>(so3_ele::ExtractInt(position, data));

  // anisotropic_permeability_directions_
  size = 0;
  so3_ele::extract_from_pack(position, data, size);
  anisotropic_permeability_directions_.resize(size, std::vector<double>(3, 0.0));
  for (int i = 0; i < size; ++i)
    so3_ele::extract_from_pack(position, data, anisotropic_permeability_directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = 0;
  so3_ele::extract_from_pack(position, data, size);
  anisotropic_permeability_nodal_coeffs_.resize(size, std::vector<double>(numnod_, 0.0));
  for (int i = 0; i < size; ++i)
    so3_ele::extract_from_pack(position, data, anisotropic_permeability_nodal_coeffs_[i]);

  // extract base class Element
  std::vector<char> basedata(0);
  so3_ele::extract_from_pack(position, data, basedata);
  so3_ele::Unpack(basedata);

  init_ = true;

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
}

template <class so3_ele, Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::So3Poro<so3_ele, distype>::Surfaces()
{
  return Core::Communication::ElementBoundaryFactory<StructuralSurface, Core::Elements::Element>(
      Core::Communication::buildSurfaces, *this);
}

template <class so3_ele, Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>>
Discret::ELEMENTS::So3Poro<so3_ele, distype>::Lines()
{
  return Core::Communication::ElementBoundaryFactory<StructuralLine, Core::Elements::Element>(
      Core::Communication::buildLines, *this);
}

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Poro<so3_ele, distype>::Print(std::ostream& os) const
{
  os << "So3_poro ";
  os << Core::FE::CellTypeToString(distype).c_str() << " ";
  Core::Elements::Element::Print(os);
}

template <class so3_ele, Core::FE::CellType distype>
bool Discret::ELEMENTS::So3Poro<so3_ele, distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, Input::LineDefinition* linedef)
{
  // read base element
  so3_ele::ReadElement(eletype, eledistype, linedef);

  // setup poro material
  Teuchos::RCP<Mat::StructPoro> poromat = Teuchos::rcp_dynamic_cast<Mat::StructPoro>(material());
  if (poromat == Teuchos::null) FOUR_C_THROW("no poro material assigned to poro element!");
  poromat->PoroSetup(numgpt_, linedef);

  read_anisotropic_permeability_directions_from_element_line_definition(linedef);
  read_anisotropic_permeability_nodal_coeffs_from_element_line_definition(linedef);

  return true;
}

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Poro<so3_ele, distype>::
    read_anisotropic_permeability_directions_from_element_line_definition(
        Input::LineDefinition* linedef)
{
  for (int dim = 0; dim < 3; ++dim)
  {
    std::string definition_name = "POROANISODIR" + std::to_string(dim + 1);
    if (linedef->HaveNamed(definition_name))
      linedef->ExtractDoubleVector(definition_name, anisotropic_permeability_directions_[dim]);
  }
}

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Poro<so3_ele, distype>::
    read_anisotropic_permeability_nodal_coeffs_from_element_line_definition(
        Input::LineDefinition* linedef)
{
  for (int dim = 0; dim < 3; ++dim)
  {
    std::string definition_name = "POROANISONODALCOEFFS" + std::to_string(dim + 1);
    if (linedef->HaveNamed(definition_name))
      linedef->ExtractDoubleVector(definition_name, anisotropic_permeability_nodal_coeffs_[dim]);
  }
}

template <class so3_ele, Core::FE::CellType distype>
void Discret::ELEMENTS::So3Poro<so3_ele, distype>::VisNames(std::map<std::string, int>& names)
{
  so3_ele::VisNames(names);
}

template <class so3_ele, Core::FE::CellType distype>
bool Discret::ELEMENTS::So3Poro<so3_ele, distype>::VisData(
    const std::string& name, std::vector<double>& data)
{
  return so3_ele::VisData(name, data);
}

template <class so3_ele, Core::FE::CellType distype>
int Discret::ELEMENTS::So3Poro<so3_ele, distype>::UniqueParObjectId() const
{
  switch (distype)
  {
    case Core::FE::CellType::tet4:
      return SoTet4PoroType::Instance().UniqueParObjectId();
    case Core::FE::CellType::tet10:
      return SoTet10PoroType::Instance().UniqueParObjectId();
    case Core::FE::CellType::hex8:
      return SoHex8PoroType::Instance().UniqueParObjectId();
    case Core::FE::CellType::hex27:
      return SoHex27PoroType::Instance().UniqueParObjectId();
    case Core::FE::CellType::nurbs27:
      return SoNurbs27PoroType::Instance().UniqueParObjectId();
    default:
      FOUR_C_THROW("unknown element type!");
      break;
  }
  return -1;
}

template <class so3_ele, Core::FE::CellType distype>
Core::Elements::ElementType& Discret::ELEMENTS::So3Poro<so3_ele, distype>::ElementType() const
{
  switch (distype)
  {
    case Core::FE::CellType::tet4:
      return SoTet4PoroType::Instance();
    case Core::FE::CellType::tet10:
      return SoTet10PoroType::Instance();
    case Core::FE::CellType::hex8:
      return SoHex8PoroType::Instance();
    case Core::FE::CellType::hex27:
      return SoHex27PoroType::Instance();
    case Core::FE::CellType::nurbs27:
      return SoNurbs27PoroType::Instance();
    default:
      FOUR_C_THROW("unknown element type!");
  }
}

template <class so3_ele, Core::FE::CellType distype>
inline Core::Nodes::Node** Discret::ELEMENTS::So3Poro<so3_ele, distype>::Nodes()
{
  return so3_ele::Nodes();
}

template <class so3_ele, Core::FE::CellType distype>
inline Teuchos::RCP<Core::Mat::Material> Discret::ELEMENTS::So3Poro<so3_ele, distype>::material()
    const
{
  return so3_ele::Material();
}

template <class so3_ele, Core::FE::CellType distype>
inline int Discret::ELEMENTS::So3Poro<so3_ele, distype>::id() const
{
  return so3_ele::Id();
}

FOUR_C_NAMESPACE_CLOSE

#include "4C_so3_poro_fwd.hpp"
