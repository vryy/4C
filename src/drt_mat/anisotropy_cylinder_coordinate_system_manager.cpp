/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of a cylinder coordinate system manager

\level 3

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#include "anisotropy_cylinder_coordinate_system_manager.H"
#include "../drt_lib/drt_parobject.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"
#include "anisotropy_utils.H"
#include "anisotropy.H"

MAT::CylinderCoordinateSystemManager::CylinderCoordinateSystemManager() = default;

void MAT::CylinderCoordinateSystemManager::Pack(DRT::PackBuffer& data) const
{
  DRT::ParObject::AddtoPack(data, radial_);
  DRT::ParObject::AddtoPack(data, axial_);
  DRT::ParObject::AddtoPack(data, circumferential_);
  DRT::ParObject::AddtoPack(data, static_cast<int>(isDefined_));
}

void MAT::CylinderCoordinateSystemManager::Unpack(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  DRT::ParObject::ExtractfromPack(position, data, radial_);
  DRT::ParObject::ExtractfromPack(position, data, axial_);
  DRT::ParObject::ExtractfromPack(position, data, circumferential_);
  isDefined_ = static_cast<bool>(DRT::ParObject::ExtractInt(position, data));
}

void MAT::CylinderCoordinateSystemManager::ReadFromElementLineDefinition(
    DRT::INPUT::LineDefinition* linedef)
{
  if (linedef->HaveNamed("RAD") and linedef->HaveNamed("AXI") and linedef->HaveNamed("CIR"))
  {
    ReadAnisotropyFiber(linedef, "RAD", radial_);
    ReadAnisotropyFiber(linedef, "AXI", axial_);
    ReadAnisotropyFiber(linedef, "CIR", circumferential_);
    isDefined_ = true;
  }
}

void MAT::CylinderCoordinateSystemManager::EvaluateLocalCoordinateSystem(
    LINALG::Matrix<3, 3>& cosy) const
{
  for (int i = 0; i < 3; ++i)
  {
    cosy(i, 0) = GetRad()(i);
    cosy(i, 1) = GetAxi()(i);
    cosy(i, 2) = GetCir()(i);
  }
}

const MAT::CylinderCoordinateSystemManager& MAT::Anisotropy::GetElementCylinderCoordinateSystem()
    const
{
  return elementCylinderCoordinateSystemManager_.get();
}

const MAT::CylinderCoordinateSystemManager& MAT::Anisotropy::GetGPCylinderCoordinateSystem(
    const int gp) const
{
  return gpCylinderCoordinateSystemManagers_[gp];
}