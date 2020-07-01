/*----------------------------------------------------------------------*/
/*! \file

\brief Object that stores the relevant data for a single output file.

\level 3

*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VTU_OUTPUT_WRITER_VISUALIZATION_H
#define BEAM_TO_SOLID_VTU_OUTPUT_WRITER_VISUALIZATION_H


#include "../drt_io/runtime_vtu_writer.H"

#include <Teuchos_RCP.hpp>


// Forward declarations.
namespace STR
{
  namespace TIMINT
  {
    class ParamsRuntimeVtkOutput;
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
  class BeamToSolidVtuOutputWriterVisualization : public RuntimeVtuWriter
  {
   public:
    /**
     * \brief Constructor. The underlying RuntimeVtuWriter will be set up in this method.
     * @param writer_full_name (in) Full name of this visualization on disc.
     * @param vtk_params (in) Global VTK parameter pointer.
     * @param restart_time (in) Initial time of the simulation.
     */
    BeamToSolidVtuOutputWriterVisualization(const std::string& writer_full_name,
        Teuchos::RCP<const STR::TIMINT::ParamsRuntimeVtkOutput> vtk_params, double restart_time);

    /**
     * \brief Destructor.
     */
    virtual ~BeamToSolidVtuOutputWriterVisualization(){};

    /**
     * \brief Return a mutable reference to the point coordinate vector of the vtu writer.
     * @param additional_reserve (in) If this is larger than 0, the allocated length of the vector
     * will be increased by additional_reserve.
     * @return Point coordinate vector.
     */
    std::vector<double>& GetMutablePointCoordinateVector(int additional_reserve = -1);

    /**
     * \brief Return a mutable reference to the field data vector.
     * @param data_name (in) Name of the data vector to return. If the name does not exist, an error
     * will be thrown.
     * @param additional_reserve (in) If this is larger than 0, the allocated length of the vector
     * will be increased by additional_reserve.
     * @return Point coordinate vector.
     */
    std::vector<double>& GetMutableFieldDataVector(std::string data_name);

    /**
     * \brief Return a mutable reference to the point data vector.
     * @param data_name (in) Name of the data vector to return. If the name does not exist, an error
     * will be thrown.
     * @param additional_reserve (in) If this is larger than 0, the allocated length of the vector
     * will be increased by additional_reserve.
     * @return Point coordinate vector.
     */
    std::vector<double>& GetMutablePointDataVector(
        std::string data_name, int additional_reserve = -1);

    /**
     * \brief Return a mutable reference to the cell data vector.
     * @param data_name (in) Name of the data vector to return. If the name does not exist, an error
     * will be thrown.
     * @param additional_reserve (in) If this is larger than 0, the allocated length of the vector
     * will be increased by additional_reserve.
     * @return Point coordinate vector.
     */
    std::vector<double>& GetMutableCellDataVector(
        std::string data_name, int additional_reserve = -1);

    /**
     * \brief Add an new (empty) field data vector to the object.
     * @param data_name (in) Name of the new data. If this name already exists, an error will be
     * thrown.
     * @param reserve (in) If this is larger than 0, the created vector will allocate the number of
     * entries.
     */
    void AddFieldDataVector(std::string data_name, int reserve = -1);

    /**
     * \brief Add an new (empty) point data vector to the object.
     * @param data_name (in) Name of the new data. If this name already exists, an error will be
     * thrown.
     * @param n_dim (in) Number of spatial dimensions of the data. When writing the object, the
     * dimension of this vector should be n_dim * n_nodes.
     * @param reserve (in) If this is larger than 0, the created vector will allocate the number of
     * entries.
     */
    void AddPointDataVector(std::string data_name, unsigned int n_dim, int reserve = -1);

    /**
     * \brief Add an new (empty) cell data vector to the object.
     * @param data_name (in) Name of the new data. If this name already exists, an error will be
     * thrown.
     * @param n_dim (in) Number of spatial dimensions of the data. When writing the object, the
     * dimension of this vector should be n_dim * n_nodes.
     * @param reserve (in) If this is larger than 0, the created vector will allocate the number of
     * entries.
     */
    void AddCellDataVector(std::string data_name, unsigned int n_dim, int reserve = -1);

    /**
     * \brief Add all nodes of a discretization to the output.
     *
     * Only the row nodes will be added on this rank. If this is called a map from the global node
     * DOFs on this rank will be created and stored in the object.
     *
     * @param discret (in) Pointer to the discretization.
     */
    void AddDiscretizationNodalReferencePosition(
        const Teuchos::RCP<const DRT::Discretization>& discret);

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
    //! Global parameters of vtk output.
    Teuchos::RCP<const STR::TIMINT::ParamsRuntimeVtkOutput> vtk_params_;

    //! Full name of this visualization.
    const std::string writer_full_name_;

    //! Discretization based on which global dof data can be written.
    Teuchos::RCP<const DRT::Discretization> discret_;

    //! Map for nodal GID of discretization.
    Teuchos::RCP<Epetra_Map> node_gid_map_;
  };

}  // namespace BEAMINTERACTION

#endif
