/*!
\file unit_geometry_pair_line_to_volume_segmentation_geometry_functions.H

\brief Create the geometry for the unit tests.

\level 1
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef UNIT_GEOMETRY_PAIR_LINE_TO_VOLUME_SEGMENTATION_LINE_ALONG_ELEMENT_SURFACE_H_
#define UNIT_GEOMETRY_PAIR_LINE_TO_VOLUME_SEGMENTATION_LINE_ALONG_ELEMENT_SURFACE_H_

#include "unit_geometry_pair_line_to_volume_segmentation.H"


namespace GEOMETRYPAIR
{
  /**
   * \brief The following code part is generated with meshpy. The function defines element
   * coordinates for unit test examples.
   */
  void xtest_line_along_element_surface_geometry(
      std::vector<Teuchos::RCP<DRT::Element>>& line_elements,
      std::vector<Teuchos::RCP<DRT::Element>>& volume_elements,
      std::vector<LINALG::TMatrix<double, 12, 1>>& q_line_elements,
      std::vector<LINALG::TMatrix<double, 9, 1>>& q_rot_line_elements,
      std::vector<LINALG::TMatrix<double, 24, 1>>& q_volume_elements,
      std::vector<Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double,
          GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>>>& geometry_pairs)
  {
    // Create the elements.
    const int dummy_node_ids[2] = {0, 1};
    line_elements.push_back(Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0)));
    line_elements.back()->SetNodeIds(2, dummy_node_ids);
    volume_elements.push_back(Teuchos::rcp(new DRT::ELEMENTS::So_hex8(1, 0)));
    volume_elements.push_back(Teuchos::rcp(new DRT::ELEMENTS::So_hex8(2, 0)));

    // Positional and tangent DOFs of the line(s).
    q_line_elements.push_back(LINALG::TMatrix<double, 12, 1>());
    q_line_elements.back()(0) = -0.5;
    q_line_elements.back()(1) = 0.0;
    q_line_elements.back()(2) = 0.0;
    q_line_elements.back()(3) = 0.7071067811865477;
    q_line_elements.back()(4) = 0.7071067811865475;
    q_line_elements.back()(5) = 0.0;
    q_line_elements.back()(6) = 0.5;
    q_line_elements.back()(7) = 0.0;
    q_line_elements.back()(8) = 0.0;
    q_line_elements.back()(9) = 0.7071067811865477;
    q_line_elements.back()(10) = 0.7071067811865475;
    q_line_elements.back()(11) = 0.0;
    q_rot_line_elements.push_back(LINALG::TMatrix<double, 9, 1>());
    q_rot_line_elements.back()(0) = 1.4821898202742554;
    q_rot_line_elements.back()(1) = 0.6139431255689366;
    q_rot_line_elements.back()(2) = 0.6139431255689366;
    q_rot_line_elements.back()(3) = 1.4821898202742554;
    q_rot_line_elements.back()(4) = 0.6139431255689366;
    q_rot_line_elements.back()(5) = 0.6139431255689366;
    q_rot_line_elements.back()(6) = 0.0;
    q_rot_line_elements.back()(7) = 0.0;
    q_rot_line_elements.back()(8) = 0.0;

    // Positional DOFs of the solid(s).
    q_volume_elements.push_back(LINALG::TMatrix<double, 24, 1>());
    q_volume_elements.back()(0) = -0.5;
    q_volume_elements.back()(1) = -0.5;
    q_volume_elements.back()(2) = 1.0;
    q_volume_elements.back()(3) = -0.5;
    q_volume_elements.back()(4) = -0.5;
    q_volume_elements.back()(5) = 0.0;
    q_volume_elements.back()(6) = -0.5;
    q_volume_elements.back()(7) = 0.5;
    q_volume_elements.back()(8) = 0.0;
    q_volume_elements.back()(9) = -0.5;
    q_volume_elements.back()(10) = 0.5;
    q_volume_elements.back()(11) = 1.0;
    q_volume_elements.back()(12) = 0.5;
    q_volume_elements.back()(13) = -0.5;
    q_volume_elements.back()(14) = 1.0;
    q_volume_elements.back()(15) = 0.5;
    q_volume_elements.back()(16) = -0.5;
    q_volume_elements.back()(17) = 0.0;
    q_volume_elements.back()(18) = 0.5;
    q_volume_elements.back()(19) = 0.5;
    q_volume_elements.back()(20) = 0.0;
    q_volume_elements.back()(21) = 0.5;
    q_volume_elements.back()(22) = 0.5;
    q_volume_elements.back()(23) = 1.0;
    q_volume_elements.push_back(LINALG::TMatrix<double, 24, 1>());
    q_volume_elements.back()(0) = -0.5;
    q_volume_elements.back()(1) = -0.5;
    q_volume_elements.back()(2) = 0.0;
    q_volume_elements.back()(3) = -0.5;
    q_volume_elements.back()(4) = -0.5;
    q_volume_elements.back()(5) = -1.0;
    q_volume_elements.back()(6) = -0.5;
    q_volume_elements.back()(7) = 0.5;
    q_volume_elements.back()(8) = -1.0;
    q_volume_elements.back()(9) = -0.5;
    q_volume_elements.back()(10) = 0.5;
    q_volume_elements.back()(11) = 0.0;
    q_volume_elements.back()(12) = 0.5;
    q_volume_elements.back()(13) = -0.5;
    q_volume_elements.back()(14) = 0.0;
    q_volume_elements.back()(15) = 0.5;
    q_volume_elements.back()(16) = -0.5;
    q_volume_elements.back()(17) = -1.0;
    q_volume_elements.back()(18) = 0.5;
    q_volume_elements.back()(19) = 0.5;
    q_volume_elements.back()(20) = -1.0;
    q_volume_elements.back()(21) = 0.5;
    q_volume_elements.back()(22) = 0.5;
    q_volume_elements.back()(23) = 0.0;
  }

  /**
   * \brief The following code part is generated with meshpy. The function defines element
   * coordinates for unit test examples.
   */
  void xtest_multiple_intersections_hex27_geometry(
      std::vector<Teuchos::RCP<DRT::Element>>& line_elements,
      std::vector<Teuchos::RCP<DRT::Element>>& volume_elements,
      std::vector<LINALG::TMatrix<double, 12, 1>>& q_line_elements,
      std::vector<LINALG::TMatrix<double, 9, 1>>& q_rot_line_elements,
      std::vector<LINALG::TMatrix<double, 81, 1>>& q_volume_elements,
      std::vector<Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double,
          GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>>>& geometry_pairs)
  {
    // Create the elements.
    const int dummy_node_ids[2] = {0, 1};
    line_elements.push_back(Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0)));
    line_elements.back()->SetNodeIds(2, dummy_node_ids);
    volume_elements.push_back(Teuchos::rcp(new DRT::ELEMENTS::So_hex27(1, 0)));

    // Positional and tangent DOFs of the line(s).
    q_line_elements.push_back(LINALG::TMatrix<double, 12, 1>());
    q_line_elements.back()(0) = -0.4;
    q_line_elements.back()(1) = -0.5;
    q_line_elements.back()(2) = 0.515;
    q_line_elements.back()(3) = 0.13302222155948917;
    q_line_elements.back()(4) = 0.754406506735489;
    q_line_elements.back()(5) = 0.6427876096865395;
    q_line_elements.back()(6) = -0.3135063030374552;
    q_line_elements.back()(7) = -0.009469868904796452;
    q_line_elements.back()(8) = 0.47142212862617094;
    q_line_elements.back()(9) = 0.13302222155948917;
    q_line_elements.back()(10) = 0.754406506735489;
    q_line_elements.back()(11) = 0.6427876096865395;
    q_rot_line_elements.push_back(LINALG::TMatrix<double, 9, 1>());
    q_rot_line_elements.back()(0) = -1.9910582384581363;
    q_rot_line_elements.back()(1) = -1.6706962335427067;
    q_rot_line_elements.back()(2) = -0.7246859334892601;
    q_rot_line_elements.back()(3) = -1.9910582384581363;
    q_rot_line_elements.back()(4) = -1.6706962335427067;
    q_rot_line_elements.back()(5) = -0.7246859334892601;
    q_rot_line_elements.back()(6) = 1.228315266067868;
    q_rot_line_elements.back()(7) = 1.1255435698194929;
    q_rot_line_elements.back()(8) = 1.0142058073830034;

    // Positional DOFs of the solid(s).
    q_volume_elements.push_back(LINALG::TMatrix<double, 81, 1>());
    q_volume_elements.back()(0) = -0.45;
    q_volume_elements.back()(1) = -0.425;
    q_volume_elements.back()(2) = 0.6;
    q_volume_elements.back()(3) = -0.5;
    q_volume_elements.back()(4) = -0.5;
    q_volume_elements.back()(5) = -0.5;
    q_volume_elements.back()(6) = -0.5;
    q_volume_elements.back()(7) = 0.5;
    q_volume_elements.back()(8) = -0.5;
    q_volume_elements.back()(9) = -0.5;
    q_volume_elements.back()(10) = 0.5;
    q_volume_elements.back()(11) = 0.5;
    q_volume_elements.back()(12) = 0.5;
    q_volume_elements.back()(13) = -0.5;
    q_volume_elements.back()(14) = 0.5;
    q_volume_elements.back()(15) = 0.5;
    q_volume_elements.back()(16) = -0.5;
    q_volume_elements.back()(17) = -0.5;
    q_volume_elements.back()(18) = 0.5;
    q_volume_elements.back()(19) = 0.5;
    q_volume_elements.back()(20) = -0.5;
    q_volume_elements.back()(21) = 0.5;
    q_volume_elements.back()(22) = 0.5;
    q_volume_elements.back()(23) = 0.5;
    q_volume_elements.back()(24) = -0.5;
    q_volume_elements.back()(25) = -0.5;
    q_volume_elements.back()(26) = 0.0;
    q_volume_elements.back()(27) = -0.5;
    q_volume_elements.back()(28) = 0.0;
    q_volume_elements.back()(29) = -0.5;
    q_volume_elements.back()(30) = -0.5;
    q_volume_elements.back()(31) = 0.5;
    q_volume_elements.back()(32) = 0.0;
    q_volume_elements.back()(33) = -0.4;
    q_volume_elements.back()(34) = -0.1;
    q_volume_elements.back()(35) = 0.45;
    q_volume_elements.back()(36) = 0.0;
    q_volume_elements.back()(37) = -0.5;
    q_volume_elements.back()(38) = 0.5;
    q_volume_elements.back()(39) = 0.0;
    q_volume_elements.back()(40) = -0.5;
    q_volume_elements.back()(41) = -0.5;
    q_volume_elements.back()(42) = 0.0;
    q_volume_elements.back()(43) = 0.5;
    q_volume_elements.back()(44) = -0.5;
    q_volume_elements.back()(45) = 0.0;
    q_volume_elements.back()(46) = 0.5;
    q_volume_elements.back()(47) = 0.5;
    q_volume_elements.back()(48) = 0.5;
    q_volume_elements.back()(49) = -0.5;
    q_volume_elements.back()(50) = 0.0;
    q_volume_elements.back()(51) = 0.5;
    q_volume_elements.back()(52) = 0.0;
    q_volume_elements.back()(53) = -0.5;
    q_volume_elements.back()(54) = 0.5;
    q_volume_elements.back()(55) = 0.5;
    q_volume_elements.back()(56) = 0.0;
    q_volume_elements.back()(57) = 0.5;
    q_volume_elements.back()(58) = 0.0;
    q_volume_elements.back()(59) = 0.5;
    q_volume_elements.back()(60) = -0.5;
    q_volume_elements.back()(61) = 0.0;
    q_volume_elements.back()(62) = 0.0;
    q_volume_elements.back()(63) = 0.0;
    q_volume_elements.back()(64) = -0.5;
    q_volume_elements.back()(65) = 0.0;
    q_volume_elements.back()(66) = 0.0;
    q_volume_elements.back()(67) = 0.0;
    q_volume_elements.back()(68) = -0.5;
    q_volume_elements.back()(69) = 0.0;
    q_volume_elements.back()(70) = 0.5;
    q_volume_elements.back()(71) = 0.0;
    q_volume_elements.back()(72) = 0.0;
    q_volume_elements.back()(73) = 0.0;
    q_volume_elements.back()(74) = 0.5;
    q_volume_elements.back()(75) = 0.5;
    q_volume_elements.back()(76) = 0.0;
    q_volume_elements.back()(77) = 0.0;
    q_volume_elements.back()(78) = 0.0;
    q_volume_elements.back()(79) = 0.0;
    q_volume_elements.back()(80) = 0.0;
  }

  /**
   * \brief The following code part is generated with meshpy. The function defines element
   * coordinates for unit test examples.
   */
  void xtest_multiple_intersections_tet10_geometry(
      std::vector<Teuchos::RCP<DRT::Element>>& line_elements,
      std::vector<Teuchos::RCP<DRT::Element>>& volume_elements,
      std::vector<LINALG::TMatrix<double, 12, 1>>& q_line_elements,
      std::vector<LINALG::TMatrix<double, 9, 1>>& q_rot_line_elements,
      std::vector<LINALG::TMatrix<double, 30, 1>>& q_volume_elements,
      std::vector<Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double,
          GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>>>& geometry_pairs)
  {
    // Create the elements.
    const int dummy_node_ids[2] = {0, 1};
    line_elements.push_back(Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0)));
    line_elements.back()->SetNodeIds(2, dummy_node_ids);
    volume_elements.push_back(Teuchos::rcp(new DRT::ELEMENTS::So_tet10(1, 0)));

    // Positional and tangent DOFs of the line(s).
    q_line_elements.push_back(LINALG::TMatrix<double, 12, 1>());
    q_line_elements.back()(0) = 0.3984274010784635;
    q_line_elements.back()(1) = 0.112954622498459;
    q_line_elements.back()(2) = 0.4240495331815059;
    q_line_elements.back()(3) = -0.7377092435057686;
    q_line_elements.back()(4) = -0.6531007789277388;
    q_line_elements.back()(5) = -0.17101007166283447;
    q_line_elements.back()(6) = 0.10157259892153653;
    q_line_elements.back()(7) = 0.487045377501541;
    q_line_elements.back()(8) = 0.27595046681849406;
    q_line_elements.back()(9) = -0.9253422132715935;
    q_line_elements.back()(10) = -0.2495541519597498;
    q_line_elements.back()(11) = -0.28541988994686207;
    q_rot_line_elements.push_back(LINALG::TMatrix<double, 9, 1>());
    q_rot_line_elements.back()(0) = 0.9576420862321825;
    q_rot_line_elements.back()(1) = -2.1136886557085592;
    q_rot_line_elements.back()(2) = -1.6586847488203964;
    q_rot_line_elements.back()(3) = -0.5899119884260623;
    q_rot_line_elements.back()(4) = 2.304864621069406;
    q_rot_line_elements.back()(5) = 1.9641024934984173;
    q_rot_line_elements.back()(6) = -0.4343496589483586;
    q_rot_line_elements.back()(7) = -0.03800067116433282;
    q_rot_line_elements.back()(8) = 2.2410533701016377;

    // Positional DOFs of the solid(s).
    q_volume_elements.push_back(LINALG::TMatrix<double, 30, 1>());
    q_volume_elements.back()(0) = -1.0;
    q_volume_elements.back()(1) = 1.224646799147353e-16;
    q_volume_elements.back()(2) = -0.5;
    q_volume_elements.back()(3) = 0.2;
    q_volume_elements.back()(4) = 0.3;
    q_volume_elements.back()(5) = 0.6;
    q_volume_elements.back()(6) = 0.5;
    q_volume_elements.back()(7) = -0.8660254037844386;
    q_volume_elements.back()(8) = -0.5;
    q_volume_elements.back()(9) = 0.5;
    q_volume_elements.back()(10) = 0.8660254037844386;
    q_volume_elements.back()(11) = -0.5;
    q_volume_elements.back()(12) = -0.5;
    q_volume_elements.back()(13) = 6.123233995736765e-17;
    q_volume_elements.back()(14) = 0.0;
    q_volume_elements.back()(15) = 0.45;
    q_volume_elements.back()(16) = -0.1330127018922193;
    q_volume_elements.back()(17) = 0.1;
    q_volume_elements.back()(18) = -0.25;
    q_volume_elements.back()(19) = -0.43301270189221924;
    q_volume_elements.back()(20) = -0.5;
    q_volume_elements.back()(21) = -0.25;
    q_volume_elements.back()(22) = 0.43301270189221935;
    q_volume_elements.back()(23) = -0.5;
    q_volume_elements.back()(24) = 0.15;
    q_volume_elements.back()(25) = 0.6330127018922194;
    q_volume_elements.back()(26) = 0.05;
    q_volume_elements.back()(27) = 0.5;
    q_volume_elements.back()(28) = 0.0;
    q_volume_elements.back()(29) = -0.5;
  }

}  // namespace GEOMETRYPAIR


#endif
