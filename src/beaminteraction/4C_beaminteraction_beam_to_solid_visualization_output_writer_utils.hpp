// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VISUALIZATION_OUTPUT_WRITER_UTILS_HPP
#define FOUR_C_BEAMINTERACTION_BEAM_TO_SOLID_VISUALIZATION_OUTPUT_WRITER_UTILS_HPP


#include "4C_config.hpp"

#include "4C_linalg_multi_vector.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

#include <memory>
#include <unordered_map>
#include <vector>


FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  template <unsigned int rows, unsigned int cols, class ValueType>
  class Matrix;
}
namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE
namespace BEAMINTERACTION
{
  class BeamToSolidOutputWriterVisualization;
}
namespace GEOMETRYPAIR
{
  class FaceElement;
}


namespace BEAMINTERACTION
{
  /**
   * \brief Add the nodal forces of beams and solid discretization to the output writer.
   * @param visualization (in/out) Output writer containting the visualization.
   * @param discret_ptr (in) Pointer to the discretization.
   * @param displacement (in) Global displacement vector.
   * @param force (in) Global force vector.
   * @param write_unique_ids (in) If unique IDs should be written.
   */
  void add_beam_interaction_nodal_forces(
      const std::shared_ptr<BEAMINTERACTION::BeamToSolidOutputWriterVisualization>& visualization,
      const std::shared_ptr<const Core::FE::Discretization>& discret_ptr,
      const std::shared_ptr<const Core::LinAlg::MultiVector<double>>& displacement,
      const Core::LinAlg::MultiVector<double>& force, const bool write_unique_ids = false);

  /**
   * \brief Add the averaged normal fields to the output writer.
   * @param output_writer_base_ptr (in/out) Output writer containting the visualization.
   * @param face_elements (in) Map with all face elements that should be written to the output
   * writer.
   * @param condition_coupling_id (in) Id of the coupling condition that this normal filed is part
   * of.
   * @param write_unique_ids (in) If unique IDs should be written.
   */
  void add_averaged_nodal_normals(
      BEAMINTERACTION::BeamToSolidOutputWriterVisualization& output_writer_base_ptr,
      const std::unordered_map<int, std::shared_ptr<GEOMETRYPAIR::FaceElement>>& face_elements,
      const int condition_coupling_id, const bool write_unique_ids = false);

  /**
   * \brief Split the discrete coupling forces into beam and solid ones and then calculate the
   * resultants.
   *
   * Also compute the moment resultants around the origin.
   *
   * @param discret (in) discretization.
   * @param force (in) Global coupling force vector.
   * @param displacement (in) Global displacement vector.
   * @param beam_resultant (out) Matrix with the force and moment resultants for the beam nodes.
   * @param solid_resultant (out) Matrix with the force and moment resultants for the solid nodes.
   */
  void get_global_coupling_force_resultants(const Core::FE::Discretization& discret,
      const Core::LinAlg::MultiVector<double>& force,
      const Core::LinAlg::MultiVector<double>& displacement,
      Core::LinAlg::Matrix<3, 2, double>& beam_resultant,
      Core::LinAlg::Matrix<3, 2, double>& solid_resultant);

  /**
   * \brief Add node contributions to the force and moment data arrays.
   * @param local_force (in) Force on the node.
   * @param local_position (in) Position of the node.
   * @param resultant (in/out) Array with the force and moment resultants.
   */
  void get_node_coupling_force_resultants(const std::vector<double>& local_force,
      const std::vector<double>& local_position, Core::LinAlg::Matrix<3, 2, double>& resultant);

}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
