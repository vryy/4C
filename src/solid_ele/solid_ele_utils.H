/*----------------------------------------------------------------------*/
/*! \file

\brief little helpers for solid elements

\level 1

*----------------------------------------------------------------------*/

#ifndef SOLID_ELE_UTILS_H
#define SOLID_ELE_UTILS_H

#include "discretization_fem_general_utils_integration.H"
#include "lib_element.H"
#include "inpar_structure.H"
#include "solid_ele_calc_eas.H"

namespace STR::ELEMENTS
{
  enum class EasType;
}

namespace STR::UTILS
{
  template <CORE::DRT::UTILS::GaussRule3D rule>
  struct IntRuleToNquad
  {
  };
  template <>
  struct IntRuleToNquad<CORE::DRT::UTILS::GaussRule3D::hex_8point>
  {
    static const int ngp = 8;
  };
  template <>
  struct IntRuleToNquad<CORE::DRT::UTILS::GaussRule3D::hex_27point>
  {
    static const int ngp = 27;
  };
  template <>
  struct IntRuleToNquad<CORE::DRT::UTILS::GaussRule3D::hex_18point>
  {
    static const int ngp = 18;
  };
  template <>
  struct IntRuleToNquad<CORE::DRT::UTILS::GaussRule3D::pyramid_1point>
  {
    static const int ngp = 1;
  };
  template <>
  struct IntRuleToNquad<CORE::DRT::UTILS::GaussRule3D::pyramid_8point>
  {
    static const int ngp = 8;
  };
  template <>
  struct IntRuleToNquad<CORE::DRT::UTILS::GaussRule3D::tet_1point>
  {
    static const int ngp = 1;
  };
  template <>
  struct IntRuleToNquad<CORE::DRT::UTILS::GaussRule3D::tet_4point>
  {
    static const int ngp = 4;
  };
  template <>
  struct IntRuleToNquad<CORE::DRT::UTILS::GaussRule3D::tet_11point>
  {
    static const int ngp = 11;
  };



  int DisTypeToNgpOptGaussRule(DRT::Element::DiscretizationType distype);

  void NodalBlockInformationSolid(DRT::Element* dwele, int& numdf, int& dimns, int& nv, int& np);

  void Pk2ToCauchy(const LINALG::Matrix<6, 1>& pk2, const LINALG::Matrix<3, 3>& defgrd,
      LINALG::Matrix<6, 1>& cauchy);

  LINALG::Matrix<6, 1> GreenLagrangeToEulerAlmansi(
      const LINALG::Matrix<6, 1>& gl, const LINALG::Matrix<3, 3>& defgrd);


  namespace ELEMENTDEFINITION
  {
    void AddElementDefinitionElementNodes(
        DRT::INPUT::LineDefinition& defs, std::basic_string<char> name, int numnodes);

    void AddElementDefinitionElementMaterial(DRT::INPUT::LineDefinition& defs);

    void AddElementDefinitionElementKinematics(DRT::INPUT::LineDefinition& defs);

    template <int nsd>
    void AddElementDefinitionCylinderCoordinates(DRT::INPUT::LineDefinition& defs);

    template <int nsd>
    void AddElementDefinitionFiber(DRT::INPUT::LineDefinition& defs);

    void AddElementDefinitionEAS(DRT::INPUT::LineDefinition& defs);

    void AddElementDefinitionFBAR(DRT::INPUT::LineDefinition& defs);
  }  // namespace ELEMENTDEFINITION

  namespace READELEMENT
  {
    int ReadElementMaterial(DRT::INPUT::LineDefinition* linedef);

    INPAR::STR::KinemType ReadElementKinematicType(DRT::INPUT::LineDefinition* linedef);

    void ReadAndSetEAS(DRT::INPUT::LineDefinition* linedef, ::STR::ELEMENTS::EasType& eastype,
        std::set<INPAR::STR::EleTech>& eletech);


  }  // namespace READELEMENT

}  // namespace STR::UTILS

#endif  // SOLID_ELE_UTILS_H
