/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration of a cylinder coordinate system anisotropy extension to be used by anisotropic
materials with @MAT::Anisotropy

\level 3


*/
/*----------------------------------------------------------------------*/

#include "anisotropy.H"
#include "anisotropy_extension_cylinder_cosy.H"
#include "../drt_lib/drt_parobject.H"

MAT::CylinderCoordinateSystemAnisotropyExtension::CylinderCoordinateSystemAnisotropyExtension()
    : cosyLocation_(CosyLocation::None)
{
}

void MAT::CylinderCoordinateSystemAnisotropyExtension::PackAnisotropy(DRT::PackBuffer& data) const
{
  DRT::ParObject::AddtoPack(data, static_cast<int>(cosyLocation_));
}

void MAT::CylinderCoordinateSystemAnisotropyExtension::UnpackAnisotropy(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  cosyLocation_ = static_cast<CosyLocation>(DRT::ParObject::ExtractInt(position, data));
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