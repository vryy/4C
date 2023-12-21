/*! \file

\brief Helpers for solid elements

\level 1
*/

#ifndef BACI_SOLID_ELE_UTILS_H
#define BACI_SOLID_ELE_UTILS_H

#include "baci_config.H"

#include "baci_inpar_structure.H"
#include "baci_lib_element.H"
#include "baci_solid_ele_properties.H"

BACI_NAMESPACE_OPEN

namespace STR::UTILS
{
  void NodalBlockInformationSolid(DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

  void Pk2ToCauchy(const CORE::LINALG::Matrix<6, 1>& pk2, const CORE::LINALG::Matrix<3, 3>& defgrd,
      CORE::LINALG::Matrix<6, 1>& cauchy);

  CORE::LINALG::Matrix<6, 1> GreenLagrangeToEulerAlmansi(
      const CORE::LINALG::Matrix<6, 1>& gl, const CORE::LINALG::Matrix<3, 3>& defgrd);

  namespace READELEMENT
  {
    int ReadElementMaterial(INPUT::LineDefinition* linedef);

    INPAR::STR::KinemType ReadElementKinematicType(INPUT::LineDefinition* linedef);

    DRT::ELEMENTS::ElementTechnology ReadElementTechnology(INPUT::LineDefinition* linedef);

    DRT::ELEMENTS::PrestressTechnology ReadPrestressTechnology(INPUT::LineDefinition* linedef);
  }  // namespace READELEMENT

}  // namespace STR::UTILS

BACI_NAMESPACE_CLOSE

#endif  // SOLID_ELE_UTILS_H
