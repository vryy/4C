/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of a cylinder coordinate system anisotropy extension to be used by anisotropic
materials with @MAT::Anisotropy

\level 3


*/
/*----------------------------------------------------------------------*/

#include "baci_mat_anisotropy_extension_cylinder_cosy.hpp"

#include "baci_comm_parobject.hpp"
#include "baci_mat_anisotropy.hpp"
#include "baci_mat_anisotropy_coordinate_system_provider.hpp"

FOUR_C_NAMESPACE_OPEN

MAT::CylinderCoordinateSystemAnisotropyExtension::CylinderCoordinateSystemAnisotropyExtension()
    : cosyLocation_(CosyLocation::None)
{
}

void MAT::CylinderCoordinateSystemAnisotropyExtension::PackAnisotropy(
    CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::ParObject::AddtoPack(data, static_cast<int>(cosyLocation_));
}

void MAT::CylinderCoordinateSystemAnisotropyExtension::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  cosyLocation_ = static_cast<CosyLocation>(CORE::COMM::ParObject::ExtractInt(position, data));
}

void MAT::CylinderCoordinateSystemAnisotropyExtension::OnGlobalDataInitialized()
{
  if (GetAnisotropy()->HasGPCylinderCoordinateSystem())
  {
    cosyLocation_ = CosyLocation::GPCosy;
  }
  else if (GetAnisotropy()->HasElementCylinderCoordinateSystem())
  {
    cosyLocation_ = CosyLocation::ElementCosy;
  }
  else
  {
    cosyLocation_ = CosyLocation::None;
  }
}

void MAT::CylinderCoordinateSystemAnisotropyExtension::OnGlobalElementDataInitialized()
{
  // do nothing
}

void MAT::CylinderCoordinateSystemAnisotropyExtension::OnGlobalGPDataInitialized()
{
  // do nothing
}

const MAT::CylinderCoordinateSystemProvider&
MAT::CylinderCoordinateSystemAnisotropyExtension::GetCylinderCoordinateSystem(int gp) const
{
  if (cosyLocation_ == CosyLocation::None)
  {
    dserror("No cylinder coordinate system defined!");
  }

  if (cosyLocation_ == CosyLocation::ElementCosy)
  {
    return GetAnisotropy()->GetElementCylinderCoordinateSystem();
  }

  return GetAnisotropy()->GetGPCylinderCoordinateSystem(gp);
}

Teuchos::RCP<MAT::CoordinateSystemProvider>
MAT::CylinderCoordinateSystemAnisotropyExtension::GetCoordinateSystemProvider(int gp) const
{
  auto cosy = Teuchos::rcp(new CoordinateSystemHolder());

  if (cosyLocation_ != CosyLocation::None)
    cosy->SetCylinderCoordinateSystemProvider(Teuchos::rcpFromRef(GetCylinderCoordinateSystem(gp)));

  return cosy;
}
FOUR_C_NAMESPACE_CLOSE
