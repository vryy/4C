/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the beam-to-solid vtu output writers.

\level 3

\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VTU_OUTPUT_WRITER_UTILS_H
#define BEAM_TO_SOLID_VTU_OUTPUT_WRITER_UTILS_H


#include <Teuchos_RCP.hpp>
#include <unordered_map>


// Forward declarations.
class Epetra_MultiVector;
namespace DRT
{
  class Discretization;
}
namespace BEAMINTERACTION
{
  class BeamToSolidVtuOutputWriterVisualization;
}
namespace GEOMETRYPAIR
{
  class FaceElement;
}


namespace BEAMINTERACTION
{
  /**
   * \brief Add the nodal forces of beams and solid discretization to the output writer.
   * @param output_writer_base_ptr (in/out) Output writer containting the visualization.
   * @param discret_ptr (in) Pointer to the discretization.
   * @param displacement (in) Global displacement vector.
   * @param force (in) Global force vector.
   */
  void AddBeamInteractionNodalForces(
      const Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>&
          output_writer_base_ptr,
      const Teuchos::RCP<const DRT::Discretization>& discret_ptr,
      const Teuchos::RCP<const Epetra_MultiVector>& displacement,
      const Teuchos::RCP<const Epetra_MultiVector>& force);

  /**
   * \brief Add the averaged normal fields to the output writer.
   * @param output_writer_base_ptr (in/out) Output writer containting the visualization.
   * @param face_elements (in) Map with all face elements that should be written to the output
   * writer.
   * @param condition_coupling_id (in) Id of the coupling condition that this normal filed is part
   * of.
   */
  void AddAveragedNodalNormals(
      const Teuchos::RCP<BEAMINTERACTION::BeamToSolidVtuOutputWriterVisualization>&
          output_writer_base_ptr,
      const std::unordered_map<int, Teuchos::RCP<GEOMETRYPAIR::FaceElement>>& face_elements,
      const int condition_coupling_id);

}  // namespace BEAMINTERACTION

#endif
