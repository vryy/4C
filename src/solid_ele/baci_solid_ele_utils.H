/*! \file

\brief Helpers for solid elements

\level 1
*/

#ifndef BACI_SOLID_ELE_UTILS_H
#define BACI_SOLID_ELE_UTILS_H

#include "baci_inpar_structure.H"
#include "baci_lib_element.H"
#include "baci_solid_ele_calc_eas.H"

namespace STR::ELEMENTS
{
  enum class EasType;
}

namespace STR::UTILS
{
  void NodalBlockInformationSolid(DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

  void Pk2ToCauchy(const CORE::LINALG::Matrix<6, 1>& pk2, const CORE::LINALG::Matrix<3, 3>& defgrd,
      CORE::LINALG::Matrix<6, 1>& cauchy);

  CORE::LINALG::Matrix<6, 1> GreenLagrangeToEulerAlmansi(
      const CORE::LINALG::Matrix<6, 1>& gl, const CORE::LINALG::Matrix<3, 3>& defgrd);

  namespace READELEMENT
  {
    int ReadElementMaterial(DRT::INPUT::LineDefinition* linedef);

    INPAR::STR::KinemType ReadElementKinematicType(DRT::INPUT::LineDefinition* linedef);

    void ReadAndSetEAS(DRT::INPUT::LineDefinition* linedef, ::STR::ELEMENTS::EasType& eastype,
        std::set<INPAR::STR::EleTech>& eletech);
  }  // namespace READELEMENT

}  // namespace STR::UTILS

#endif  // SOLID_ELE_UTILS_H
