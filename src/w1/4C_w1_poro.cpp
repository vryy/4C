/*----------------------------------------------------------------------------*/
/*! \file
\brief 2D wall element for structure part of porous medium.

\level 2


*/
/*---------------------------------------------------------------------------*/

#include "4C_w1_poro.hpp"

#include "4C_comm_utils_factory.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_fluidporo.hpp"
#include "4C_mat_fluidporo_multiphase.hpp"
#include "4C_mat_structporo.hpp"
#include "4C_poroelast_utils.hpp"

FOUR_C_NAMESPACE_OPEN

template <Core::FE::CellType distype>
Discret::ELEMENTS::Wall1Poro<distype>::Wall1Poro(int id, int owner)
    : Discret::ELEMENTS::Wall1(id, owner), intpoints_(distype), weights_(true), myknots_(numdim_)
{
  numgpt_ = intpoints_.NumPoints();

  invJ_.resize(numgpt_, Core::LinAlg::Matrix<numdim_, numdim_>(true));
  detJ_.resize(numgpt_, 0.0);
  xsi_.resize(numgpt_, Core::LinAlg::Matrix<numdim_, 1>(true));
  anisotropic_permeability_directions_.resize(2, std::vector<double>(2, 0.0));
  anisotropic_permeability_nodal_coeffs_.resize(2, std::vector<double>(numnod_, 0.0));

  init_ = false;

  scatra_coupling_ = false;
}

template <Core::FE::CellType distype>
Discret::ELEMENTS::Wall1Poro<distype>::Wall1Poro(const Discret::ELEMENTS::Wall1Poro<distype>& old)
    : Discret::ELEMENTS::Wall1(old),
      invJ_(old.invJ_),
      detJ_(old.detJ_),
      xsi_(old.xsi_),
      intpoints_(distype),
      init_(old.init_),
      scatra_coupling_(old.scatra_coupling_),
      weights_(old.weights_),
      myknots_(old.myknots_),
      anisotropic_permeability_directions_(old.anisotropic_permeability_directions_),
      anisotropic_permeability_nodal_coeffs_(old.anisotropic_permeability_nodal_coeffs_)
{
  numgpt_ = intpoints_.NumPoints();
}

template <Core::FE::CellType distype>
Core::Elements::Element* Discret::ELEMENTS::Wall1Poro<distype>::Clone() const
{
  auto* newelement = new Discret::ELEMENTS::Wall1Poro<distype>(*this);
  return newelement;
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1Poro<distype>::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // detJ_
  AddtoPack(data, detJ_);

  // invJ_
  int size = static_cast<int>(invJ_.size());
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, invJ_[i]);

  // xsi_
  size = static_cast<int>(xsi_.size());
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, xsi_[i]);

  // scatra_coupling_
  AddtoPack(data, scatra_coupling_);

  // anisotropic_permeability_directions_
  size = static_cast<int>(anisotropic_permeability_directions_.size());
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, anisotropic_permeability_directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = static_cast<int>(anisotropic_permeability_nodal_coeffs_.size());
  AddtoPack(data, size);
  for (int i = 0; i < size; ++i) AddtoPack(data, anisotropic_permeability_nodal_coeffs_[i]);

  // add base class Element
  Discret::ELEMENTS::Wall1::Pack(data);
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1Poro<distype>::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // detJ_
  ExtractfromPack(position, data, detJ_);

  // invJ_
  int size = 0;
  ExtractfromPack(position, data, size);
  invJ_.resize(size, Core::LinAlg::Matrix<numdim_, numdim_>(true));
  for (int i = 0; i < size; ++i) ExtractfromPack(position, data, invJ_[i]);

  // xsi_
  size = 0;
  ExtractfromPack(position, data, size);
  xsi_.resize(size, Core::LinAlg::Matrix<numdim_, 1>(true));
  for (int i = 0; i < size; ++i) ExtractfromPack(position, data, xsi_[i]);

  // scatra_coupling_
  scatra_coupling_ = static_cast<bool>(ExtractInt(position, data));

  // anisotropic_permeability_directions_
  size = 0;
  ExtractfromPack(position, data, size);
  anisotropic_permeability_directions_.resize(size, std::vector<double>(3, 0.0));
  for (int i = 0; i < size; ++i)
    ExtractfromPack(position, data, anisotropic_permeability_directions_[i]);

  // anisotropic_permeability_nodal_coeffs_
  size = 0;
  ExtractfromPack(position, data, size);
  anisotropic_permeability_nodal_coeffs_.resize(size, std::vector<double>(numnod_, 0.0));
  for (int i = 0; i < size; ++i)
    ExtractfromPack(position, data, anisotropic_permeability_nodal_coeffs_[i]);

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Discret::ELEMENTS::Wall1::Unpack(basedata);


  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", static_cast<int>(data.size()), position);

  init_ = true;
}

template <Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Wall1Poro<distype>::Lines()
{
  return Core::Communication::ElementBoundaryFactory<Wall1Line, Wall1Poro>(
      Core::Communication::buildLines, *this);
}

template <Core::FE::CellType distype>
std::vector<Teuchos::RCP<Core::Elements::Element>> Discret::ELEMENTS::Wall1Poro<distype>::Surfaces()
{
  return {Teuchos::rcpFromRef(*this)};
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1Poro<distype>::Print(std::ostream& os) const
{
  os << "Wall1_Poro ";
  Element::Print(os);
  std::cout << std::endl;
}

template <Core::FE::CellType distype>
bool Discret::ELEMENTS::Wall1Poro<distype>::ReadElement(
    const std::string& eletype, const std::string& eledistype, Input::LineDefinition* linedef)
{
  // read base element
  Wall1::ReadElement(eletype, eledistype, linedef);

  // setup poro material
  Teuchos::RCP<Mat::StructPoro> poromat = Teuchos::rcp_dynamic_cast<Mat::StructPoro>(Material());
  if (poromat == Teuchos::null)
    FOUR_C_THROW("material assigned to poro element is not a poro material!");
  poromat->PoroSetup(numgpt_, linedef);

  read_anisotropic_permeability_directions_from_element_line_definition(linedef);
  read_anisotropic_permeability_nodal_coeffs_from_element_line_definition(linedef);

  return true;
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1Poro<distype>::
    read_anisotropic_permeability_directions_from_element_line_definition(
        Input::LineDefinition* linedef)
{
  for (int dim = 0; dim < 2; ++dim)
  {
    std::string definition_name = "POROANISODIR" + std::to_string(dim + 1);
    if (linedef->HaveNamed(definition_name))
      linedef->ExtractDoubleVector(definition_name, anisotropic_permeability_directions_[dim]);
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1Poro<distype>::
    read_anisotropic_permeability_nodal_coeffs_from_element_line_definition(
        Input::LineDefinition* linedef)
{
  for (int dim = 0; dim < 2; ++dim)
  {
    std::string definition_name = "POROANISONODALCOEFFS" + std::to_string(dim + 1);
    if (linedef->HaveNamed(definition_name))
      linedef->ExtractDoubleVector(definition_name, anisotropic_permeability_nodal_coeffs_[dim]);
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1Poro<distype>::get_materials()
{
  // get structure material
  if (struct_mat_ == Teuchos::null)
  {
    struct_mat_ = Teuchos::rcp_dynamic_cast<Mat::StructPoro>(Material());
    if (struct_mat_ == Teuchos::null) FOUR_C_THROW("cast to poro material failed");

    if (struct_mat_->MaterialType() != Core::Materials::m_structporo and
        struct_mat_->MaterialType() != Core::Materials::m_structpororeaction and
        struct_mat_->MaterialType() != Core::Materials::m_structpororeactionECM)
      FOUR_C_THROW("invalid structure material for poroelasticity");
  }

  // get fluid material
  if (fluid_mat_ == Teuchos::null)
  {
    // access second material in structure element
    if (NumMaterial() > 1)
    {
      fluid_mat_ = Teuchos::rcp_dynamic_cast<Mat::FluidPoro>(Material(1));
      if (fluid_mat_ == Teuchos::null) return;
      // FOUR_C_THROW("cast to fluid poro material failed");
      if (fluid_mat_->MaterialType() != Core::Materials::m_fluidporo)
        FOUR_C_THROW("invalid fluid material for poroelasticity");
    }
    else
      FOUR_C_THROW("no second material defined for element %i", Id());
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1Poro<distype>::get_materials_pressure_based()
{
  // get structure material
  if (struct_mat_ == Teuchos::null)
  {
    struct_mat_ = Teuchos::rcp_dynamic_cast<Mat::StructPoro>(Material());
    if (struct_mat_ == Teuchos::null) FOUR_C_THROW("cast to poro material failed");

    if (struct_mat_->MaterialType() != Core::Materials::m_structporo and
        struct_mat_->MaterialType() != Core::Materials::m_structpororeaction and
        struct_mat_->MaterialType() != Core::Materials::m_structpororeactionECM)
      FOUR_C_THROW("invalid structure material for poroelasticity");
  }

  // Get Fluid-multiphase-Material
  if (fluidmulti_mat_ == Teuchos::null)
  {
    // access second material in structure element
    if (NumMaterial() > 1)
    {
      fluidmulti_mat_ = Teuchos::rcp_dynamic_cast<Mat::FluidPoroMultiPhase>(Material(1));
      if (fluidmulti_mat_ == Teuchos::null)
        FOUR_C_THROW("cast to multiphase fluid poro material failed");
      if (fluidmulti_mat_->MaterialType() != Core::Materials::m_fluidporo_multiphase and
          fluidmulti_mat_->MaterialType() != Core::Materials::m_fluidporo_multiphase_reactions)
        FOUR_C_THROW("invalid fluid material for poro-multiphase-elasticity");
      if (fluidmulti_mat_->NumFluidPhases() == 0)
      {
        FOUR_C_THROW(
            "NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE = 0 currently not supported since this requires "
            "an adaption of the definition of the solid pressure");
      }
    }
    else
      FOUR_C_THROW("no second material defined for element %i", Id());
  }
}

template <Core::FE::CellType distype>
void Discret::ELEMENTS::Wall1Poro<distype>::VisNames(std::map<std::string, int>& names)
{
  SolidMaterial()->VisNames(names);
}

template <Core::FE::CellType distype>
bool Discret::ELEMENTS::Wall1Poro<distype>::VisData(
    const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (Wall1::VisData(name, data)) return true;

  return SolidMaterial()->VisData(name, data, numgpt_, this->Id());
}

template class Discret::ELEMENTS::Wall1Poro<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::Wall1Poro<Core::FE::CellType::quad4>;
template class Discret::ELEMENTS::Wall1Poro<Core::FE::CellType::quad9>;
template class Discret::ELEMENTS::Wall1Poro<Core::FE::CellType::nurbs4>;
template class Discret::ELEMENTS::Wall1Poro<Core::FE::CellType::nurbs9>;

FOUR_C_NAMESPACE_CLOSE
