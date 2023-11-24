/*----------------------------------------------------------------------*/
/*! \file

\brief Object that stores the relevant data for a single output file.

\level 3

*/
// End doxygen header.


#ifndef BACI_BEAMINTERACTION_BEAM_TO_SOLID_VTU_OUTPUT_WRITER_VISUALIZATION_H
#define BACI_BEAMINTERACTION_BEAM_TO_SOLID_VTU_OUTPUT_WRITER_VISUALIZATION_H


#include "baci_io_visualization_manager.H"

#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace STR
{
  namespace TIMINT
  {
    class ParamsRuntimeOutput;
  }
}  // namespace STR
namespace DRT
{
  class Discretization;
}
class Epetra_Map;
class Epetra_MultiVector;


namespace BEAMINTERACTION
{
  /**
   * \brief Object to write visualization for beam-to-solid interactions.
   *
   * This object represents the data that will be written to a single collection of output
   * files to disc, i.e. one file that will be opened with ParaView. The data that can be stored, is
   * point geometry, cell connected to the points, point data and cell data.
   */
  class BeamToSolidVtuOutputWriterVisualization : public IO::VisualizationManager
  {
   public:
    /**
     * \brief Constructor. The underlying RuntimeVtuWriter will be set up in this method.
     * @param writer_full_name (in) Full name of this visualization on disc.
     * @param visualization_output_params (in) Global visualization parameter pointer.
     * @param restart_time (in) Initial time of the simulation.
     */
    BeamToSolidVtuOutputWriterVisualization(const std::string& writer_full_name,
        Teuchos::RCP<const STR::TIMINT::ParamsRuntimeOutput> visualization_output_params,
        double restart_time);

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVtuOutputWriterVisualization() = default;

    /**
     * \brief Add all nodes of a discretization to the output.
     *
     * Only the row nodes will be added on this rank. If this is called a map from the global node
     * DOFs on this rank will be created and stored in the object.
     *
     * @param discret (in) Pointer to the discretization.
     */
    void AddDiscretizationNodalReferencePosition(
        const Teuchos::RCP<const ::DRT::Discretization>& discret);

    /**
     * \brief Add global DOF based data to the writer.
     *
     * The data will automatically be extracted, and added on the correct rank.
     *
     * @param data_name (in) name of the added data.
     * @param vector (in) Global state vector. The size of this vector has to be 3 * n_nodes.
     */
    void AddDiscretizationNodalData(
        const std::string& data_name, const Teuchos::RCP<const Epetra_MultiVector>& vector);

    /**
     * \brief Write the object to disc.
     * @param timestep_number (in) Number of this time step.
     * @param time (in) Time of this time step.
     */
    void Write(const unsigned int timestep_number, const double time);

   private:
    //! Global parameters of visualization output.
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeOutput> visualization_output_params_;

    //! Full name of this visualization.
    const std::string writer_full_name_;

    //! Discretization based on which global dof data can be written.
    Teuchos::RCP<const ::DRT::Discretization> discret_;

    //! Map for nodal GID of discretization.
    Teuchos::RCP<Epetra_Map> node_gid_map_;
  };

}  // namespace BEAMINTERACTION

#endif
