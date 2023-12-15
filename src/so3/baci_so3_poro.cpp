/*----------------------------------------------------------------------*/
/*! \file
\brief implementation of the 3D solid-poro element


\level 2

*----------------------------------------------------------------------*/

#include "baci_so3_poro.H"

#include "baci_comm_utils_factory.H"
#include "baci_io_linedefinition.H"
#include "baci_lib_discret.H"
#include "baci_mat_structporo.H"
#include "baci_so3_line.H"
#include "baci_so3_poro_eletypes.H"
#include "baci_so3_surface.H"

BACI_NAMESPACE_OPEN

template <class so3_ele, CORE::FE::CellType distype>
DRT::ELEMENTS::So3_Poro<so3_ele, distype>::So3_Poro(int id, int owner)
    : so3_ele(id, owner),
      data_(),
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

  invJ_.resize(numgpt_, CORE::LINALG::Matrix<numdim_, numdim_>(true));
  detJ_.resize(numgpt_, 0.0);
  xsi_.resize(numgpt_, CORE::LINALG::Matrix<numdim_, 1>(true));
  anisotropic_permeability_directions_.resize(3, std::vector<double>(3, 0.0));
  anisotropic_permeability_nodal_coeffs_.resize(3, std::vector<double>(numnod_, 0.0));
}

template <class so3_ele, CORE::FE::CellType distype>
DRT::ELEMENTS::So3_Poro<so3_ele, distype>::So3_Poro(
    const DRT::ELEMENTS::So3_Poro<so3_ele, distype>& old)
    : so3_ele(old),
      data_(old.data_),
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

template <class so3_ele, CORE::FE::CellType distype>
DRT::Element* DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Clone() const
{
  auto* newelement = new DRT::ELEMENTS::So3_Poro<so3_ele, distype>(*this);
  return newelement;
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data, type);

  // data_
  so3_ele::AddtoPack(data, data_);

  // detJ_
  so3_ele::AddtoPack(data, detJ_);

  // invJ_
  auto size = static_cast<int>(invJ_.size());
  so3_ele::AddtoPack(data, size);
  for (int i = 0; i < size; ++i) so3_ele::AddtoPack(data, invJ_[i]);

  // xsi_
  size = static_cast<int>(xsi_.size());
  so3_ele::AddtoPack(data, size);
  for (int i = 0; i < size; ++i) so3_ele::AddtoPack(data, xsi_[i]);

  // scatra_coupling_
  so3_ele::AddtoPack(data, scatra_coupling_);

  // isNurbs_
  so3_ele::AddtoPack(data, isNurbs_);

  // anisotropic_permeability_directions_
  size = static_cast<int>(anisotropic_permeability_directions_.size());
  so3_ele::AddtoPack(data, size);
  for (int i = 0; i < size; ++i) so3_ele::AddtoPack(data, anisotropic_permeability_directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = static_cast<int>(anisotropic_permeability_nodal_coeffs_.size());
  so3_ele::AddtoPack(data, size);
  for (int i = 0; i < size; ++i)
    so3_ele::AddtoPack(data, anisotropic_permeability_nodal_coeffs_[i]);

  // add base class Element
  so3_ele::Pack(data);
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // data_
  std::vector<char> tmp(0);
  so3_ele::ExtractfromPack(position, data, tmp);
  data_.Unpack(tmp);

  // detJ_
  so3_ele::ExtractfromPack(position, data, detJ_);

  // invJ_
  int size = 0;
  so3_ele::ExtractfromPack(position, data, size);
  invJ_.resize(size, CORE::LINALG::Matrix<numdim_, numdim_>(true));
  for (int i = 0; i < size; ++i) so3_ele::ExtractfromPack(position, data, invJ_[i]);

  // xsi_
  size = 0;
  so3_ele::ExtractfromPack(position, data, size);
  xsi_.resize(size, CORE::LINALG::Matrix<numdim_, 1>(true));
  for (int i = 0; i < size; ++i) so3_ele::ExtractfromPack(position, data, xsi_[i]);

  // scatra_coupling_
  scatra_coupling_ = static_cast<bool>(so3_ele::ExtractInt(position, data));

  // isNurbs_
  isNurbs_ = static_cast<bool>(so3_ele::ExtractInt(position, data));

  // anisotropic_permeability_directions_
  size = 0;
  so3_ele::ExtractfromPack(position, data, size);
  anisotropic_permeability_directions_.resize(size, std::vector<double>(3, 0.0));
  for (int i = 0; i < size; ++i)
    so3_ele::ExtractfromPack(position, data, anisotropic_permeability_directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = 0;
  so3_ele::ExtractfromPack(position, data, size);
  anisotropic_permeability_nodal_coeffs_.resize(size, std::vector<double>(numnod_, 0.0));
  for (int i = 0; i < size; ++i)
    so3_ele::ExtractfromPack(position, data, anisotropic_permeability_nodal_coeffs_[i]);

  // extract base class Element
  std::vector<char> basedata(0);
  so3_ele::ExtractfromPack(position, data, basedata);
  so3_ele::Unpack(basedata);

  init_ = true;

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);
}

template <class so3_ele, CORE::FE::CellType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Surfaces()
{
  return CORE::COMM::ElementBoundaryFactory<StructuralSurface, DRT::Element>(
      CORE::COMM::buildSurfaces, *this);
}

template <class so3_ele, CORE::FE::CellType distype>
std::vector<Teuchos::RCP<DRT::Element>> DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Lines()
{
  return CORE::COMM::ElementBoundaryFactory<StructuralLine, DRT::Element>(
      CORE::COMM::buildLines, *this);
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Print(std::ostream& os) const
{
  os << "So3_poro ";
  os << CORE::FE::CellTypeToString(distype).c_str() << " ";
  Element::Print(os);
}

template <class so3_ele, CORE::FE::CellType distype>
bool DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, DRT::INPUT::LineDefinition* linedef)
{
  // read base element
  so3_ele::ReadElement(eletype, eledistype, linedef);

  // setup poro material
  Teuchos::RCP<MAT::StructPoro> poromat = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(Material());
  if (poromat == Teuchos::null) dserror("no poro material assigned to poro element!");
  poromat->PoroSetup(numgpt_, linedef);

  ReadAnisotropicPermeabilityDirectionsFromElementLineDefinition(linedef);
  ReadAnisotropicPermeabilityNodalCoeffsFromElementLineDefinition(linedef);

  return true;
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::
    ReadAnisotropicPermeabilityDirectionsFromElementLineDefinition(
        DRT::INPUT::LineDefinition* linedef)
{
  for (int dim = 0; dim < 3; ++dim)
  {
    std::string definition_name = "POROANISODIR" + std::to_string(dim + 1);
    if (linedef->HaveNamed(definition_name))
      linedef->ExtractDoubleVector(definition_name, anisotropic_permeability_directions_[dim]);
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::
    ReadAnisotropicPermeabilityNodalCoeffsFromElementLineDefinition(
        DRT::INPUT::LineDefinition* linedef)
{
  for (int dim = 0; dim < 3; ++dim)
  {
    std::string definition_name = "POROANISONODALCOEFFS" + std::to_string(dim + 1);
    if (linedef->HaveNamed(definition_name))
      linedef->ExtractDoubleVector(definition_name, anisotropic_permeability_nodal_coeffs_[dim]);
  }
}

template <class so3_ele, CORE::FE::CellType distype>
void DRT::ELEMENTS::So3_Poro<so3_ele, distype>::VisNames(std::map<std::string, int>& names)
{
  so3_ele::VisNames(names);
}

template <class so3_ele, CORE::FE::CellType distype>
bool DRT::ELEMENTS::So3_Poro<so3_ele, distype>::VisData(
    const std::string& name, std::vector<double>& data)
{
  return so3_ele::VisData(name, data);
}

template <class so3_ele, CORE::FE::CellType distype>
int DRT::ELEMENTS::So3_Poro<so3_ele, distype>::UniqueParObjectId() const
{
  switch (distype)
  {
    case CORE::FE::CellType::tet4:
      return So_tet4PoroType::Instance().UniqueParObjectId();
    case CORE::FE::CellType::tet10:
      return So_tet10PoroType::Instance().UniqueParObjectId();
    case CORE::FE::CellType::hex8:
      return So_hex8PoroType::Instance().UniqueParObjectId();
    case CORE::FE::CellType::hex27:
      return So_hex27PoroType::Instance().UniqueParObjectId();
    case CORE::FE::CellType::nurbs27:
      return So_nurbs27PoroType::Instance().UniqueParObjectId();
    default:
      dserror("unknown element type!");
      break;
  }
  return -1;
}

template <class so3_ele, CORE::FE::CellType distype>
DRT::ElementType& DRT::ELEMENTS::So3_Poro<so3_ele, distype>::ElementType() const
{
  switch (distype)
  {
    case CORE::FE::CellType::tet4:
      return So_tet4PoroType::Instance();
    case CORE::FE::CellType::tet10:
      return So_tet10PoroType::Instance();
    case CORE::FE::CellType::hex8:
      return So_hex8PoroType::Instance();
    case CORE::FE::CellType::hex27:
      return So_hex27PoroType::Instance();
    case CORE::FE::CellType::nurbs27:
      return So_nurbs27PoroType::Instance();
    default:
      dserror("unknown element type!");
      break;
  }
  return So_hex8PoroType::Instance();
}

template <class so3_ele, CORE::FE::CellType distype>
inline DRT::Node** DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Nodes()
{
  return so3_ele::Nodes();
}

template <class so3_ele, CORE::FE::CellType distype>
inline Teuchos::RCP<MAT::Material> DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Material() const
{
  return so3_ele::Material();
}

template <class so3_ele, CORE::FE::CellType distype>
inline int DRT::ELEMENTS::So3_Poro<so3_ele, distype>::Id() const
{
  return so3_ele::Id();
}

BACI_NAMESPACE_CLOSE

#include "baci_so3_poro_fwd.hpp"
