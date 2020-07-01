/*----------------------------------------------------------------------*/
/*! \file

\brief Utility functions for the beam-to-solid vtu output writers.

\level 3

*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VTU_OUTPUT_WRITER_UTILS_H
#define BEAM_TO_SOLID_VTU_OUTPUT_WRITER_UTILS_H


#include <Teuchos_RCP.hpp>
#include <vector>
#include <unordered_map>


// Forward declarations.
class Epetra_MultiVector;
namespace LINALG
{
  template <unsigned int rows, unsigned int cols, class value_type>
  class Matrix;
}
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

  /**
   * \brief Split the discrete coupling forces into beam and solid ones and then calculate the
   * resultants.
   *
   * Also compute the moment resultants around the origin.
   *
   * @param discret (in) Discretization.
   * @param force (in) Global coupling force vector.
   * @param displacement (in) Global displacement vector.
   * @param beam_resultant (out) Matrix with the force and moment resultants for the beam nodes.
   * @param solid_resultant (out) Matrix with the force and moment resultants for the solid nodes.
   */
  void GetGlobalCouplingForceResultants(const DRT::Discretization& discret,
      const Epetra_MultiVector& force, const Epetra_MultiVector& displacement,
      LINALG::Matrix<3, 2, double>& beam_resultant, LINALG::Matrix<3, 2, double>& solid_resultant);

  /**
   * \brief Add node contributions to the force and moment data arrays.
   * @param local_force (in) Force on the node.
   * @param local_position (in) Position of the node.
   * @param resultant (in/out) Array with the force and moment resultants.
   */
  void GetNodeCouplingForceResultants(const std::vector<double>& local_force,
      const std::vector<double>& local_position, LINALG::Matrix<3, 2, double>& resultant);

}  // namespace BEAMINTERACTION

#endif
