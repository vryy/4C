/*!
\file unit_beaminteraction_beam_to_solid_volume_meshtying_mortar_geometry.H

\brief Create the geometry for the unit tests.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef UNIT_BEAMINTERACTION_BEAM_TO_SOLIV_VOLUME_MESHTYING_MORTAR_GEOMETRY_H_
#define UNIT_BEAMINTERACTION_BEAM_TO_SOLIV_VOLUME_MESHTYING_MORTAR_GEOMETRY_H_

#include "unit_beaminteraction_beam_to_solid_volume_meshtying_mortar_geometry.H"

namespace BEAMINTERACTION
{
  /**
   * \brief The following code part is generated with meshpy. The function defines element
   * coordinates for unit test examples.
   */
  void xtest_beam_in_hex8_geometry(std::vector<Teuchos::RCP<DRT::Element>>& line_elements,
      std::vector<Teuchos::RCP<DRT::Element>>& volume_elements,
      std::vector<LINALG::TMatrix<double, 12, 1>>& q_line_elements,
      std::vector<LINALG::TMatrix<double, 9, 1>>& q_rot_line_elements,
      std::vector<LINALG::TMatrix<double, 24, 1>>& q_volume_elements)
  {
    // Create the elements.
    const int dummy_node_ids[2] = {0, 1};
    line_elements.push_back(Teuchos::rcp(new DRT::ELEMENTS::Beam3r(0, 0)));
    line_elements.back()->SetNodeIds(2, dummy_node_ids);
    volume_elements.push_back(Teuchos::rcp(new DRT::ELEMENTS::So_hex8(1, 0)));

    // Positional and tangent DOFs of the line(s).
    q_line_elements.push_back(LINALG::TMatrix<double, 12, 1>());
    q_line_elements.back()(0) = -0.5;
    q_line_elements.back()(1) = 0.0;
    q_line_elements.back()(2) = 0.0;
    q_line_elements.back()(3) = 0.5773502691896255;
    q_line_elements.back()(4) = 0.5773502691896258;
    q_line_elements.back()(5) = 0.577350269189626;
    q_line_elements.back()(6) = 0.09999999999999998;
    q_line_elements.back()(7) = 0.0;
    q_line_elements.back()(8) = 0.0;
    q_line_elements.back()(9) = 0.8017837257372733;
    q_line_elements.back()(10) = -0.5345224838248488;
    q_line_elements.back()(11) = 0.2672612419124244;
    q_rot_line_elements.push_back(LINALG::TMatrix<double, 9, 1>());
    q_rot_line_elements.back()(0) = 1.674352746442651;
    q_rot_line_elements.back()(1) = 0.1425949677148126;
    q_rot_line_elements.back()(2) = 1.0831163124736984;
    q_rot_line_elements.back()(3) = 1.4331999091513161;
    q_rot_line_elements.back()(4) = -0.6560404572957742;
    q_rot_line_elements.back()(5) = -0.2491376152457331;
    q_rot_line_elements.back()(6) = 0.0;
    q_rot_line_elements.back()(7) = 0.0;
    q_rot_line_elements.back()(8) = 0.0;

    // Positional DOFs of the solid(s).
    q_volume_elements.push_back(LINALG::TMatrix<double, 24, 1>());
    q_volume_elements.back()(0) = -0.5;
    q_volume_elements.back()(1) = -0.5;
    q_volume_elements.back()(2) = 0.5;
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
  }
}  // namespace BEAMINTERACTION

#endif
