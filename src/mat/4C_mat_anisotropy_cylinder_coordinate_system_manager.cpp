/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of a cylinder coordinate system manager

\level 3


*/
/*----------------------------------------------------------------------*/

#include "4C_mat_anisotropy_cylinder_coordinate_system_manager.hpp"

#include "4C_comm_parobject.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_anisotropy.hpp"
#include "4C_mat_anisotropy_utils.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

Mat::CylinderCoordinateSystemManager::CylinderCoordinateSystemManager() = default;

void Mat::CylinderCoordinateSystemManager::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::ParObject::AddtoPack(data, radial_);
  Core::Communication::ParObject::AddtoPack(data, axial_);
  Core::Communication::ParObject::AddtoPack(data, circumferential_);
  Core::Communication::ParObject::AddtoPack(data, static_cast<int>(is_defined_));
}

void Mat::CylinderCoordinateSystemManager::Unpack(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  Core::Communication::ParObject::ExtractfromPack(position, data, radial_);
  Core::Communication::ParObject::ExtractfromPack(position, data, axial_);
  Core::Communication::ParObject::ExtractfromPack(position, data, circumferential_);
  is_defined_ = static_cast<bool>(Core::Communication::ParObject::ExtractInt(position, data));
}

void Mat::CylinderCoordinateSystemManager::read_from_element_line_definition(
    Input::LineDefinition* linedef)
{
  if (linedef->HaveNamed("RAD") and linedef->HaveNamed("AXI") and linedef->HaveNamed("CIR"))
  {
    ReadAnisotropyFiber(linedef, "RAD", radial_);
    ReadAnisotropyFiber(linedef, "AXI", axial_);
    ReadAnisotropyFiber(linedef, "CIR", circumferential_);
    is_defined_ = true;
  }
}

void Mat::CylinderCoordinateSystemManager::evaluate_local_coordinate_system(
    Core::LinAlg::Matrix<3, 3>& cosy) const
{
  for (int i = 0; i < 3; ++i)
  {
    cosy(i, 0) = GetRad()(i);
    cosy(i, 1) = GetAxi()(i);
    cosy(i, 2) = GetCir()(i);
  }
}

const Mat::CylinderCoordinateSystemManager&
Mat::Anisotropy::get_element_cylinder_coordinate_system() const
{
  return element_cylinder_coordinate_system_manager_.value();
}

const Mat::CylinderCoordinateSystemManager& Mat::Anisotropy::get_gp_cylinder_coordinate_system(
    const int gp) const
{
  return gp_cylinder_coordinate_system_managers_[gp];
}
FOUR_C_NAMESPACE_CLOSE
