/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Define a base class for objects that write visualization data to disk

\level 0

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_IO_VISUALIZATION_WRITER_BASE_HPP
#define FOUR_C_IO_VISUALIZATION_WRITER_BASE_HPP

#include "baci_config.hpp"

#include "baci_io_visualization_data.hpp"
#include "baci_io_visualization_parameters.hpp"
#include "baci_utils_exceptions.hpp"

#include <Epetra_Comm.h>

FOUR_C_NAMESPACE_OPEN

namespace IO
{
  /**
   * @brief This is a base class that defines the API for how visualization data is actually written
   * to disk
   *
   * The default implementation in this base class defines interfaces that allow the user to write
   * data to disk, i.e., create the files that can be opened with post-processing tools (e.g.,
   * ParaView). There are two different ways to do so:
   * 1. Call the following methods in that specific order (if some data structures are empty, simply
   *    pass an empty container or don't call the method):
   *      \sa InitializeTimeStep
   *      \sa WriteFieldDataToDisk
   *      \sa WriteGeometryToDisk
   *      \sa WritePointDataVectorToDisk
   *      \sa WriteCellDataVectorToDisk
   *      \sa FinalizeTimeStep
   *    This abstraction allows to treat the special case where the output vectors are so large that
   *    they create memory issues and one would like to create each individual point and cell data
   *    vector, fill it up, write it to disk and delete it (i.e., to not have to store the whole
   *    output data at the same time)
   * 2. Pass a fully populated visualization_data to \sa WriteVisualizationDataToDisk. The base
   *    implementation of this method calls the methods described in the previous point and passes
   *    the correct data.
   *
   * When deriving from this class one should aim for an implementation that only overwrites the
   * methods described in point 1, even if the writer is mainly called through \sa
   * WriteVisualizationDataToDisk. This ensures that the procedure described in point 1 is possible
   * with the derived writer. Only in the very special case that the abstraction described in point
   * 1 is not possible due to certain constraints of the writer/data format the \sa
   * WriteVisualizationDataToDisk should be overwritten in a derived class.
   */
  class VisualizationWriterBase
  {
   public:
    /**
     * @brief Default constructor
     */
    VisualizationWriterBase(IO::VisualizationParameters parameters, const Epetra_Comm& comm,
        std::string visualization_data_name);

    /**
     * @brief Default destructor
     */
    virtual ~VisualizationWriterBase() = default;

    /**
     * @brief Write a VisualizationData object to disk.
     *
     * @param visualization_data (in) Data container holding all visualization data
     * @param visualziation_time (in) Time value of current step, this is not necessarily the same
     * as the simulation time
     * @param visualization_step (in) Time step counter of current time step (does not have to be
     * continuous) this is not necessarily the same as the simulation time step counter
     */
    virtual void WriteVisualizationDataToDisk(const VisualizationData& visualization_data,
        double visualziation_time, int visualization_step);

    /**
     * @brief Initialize the current time step
     *
     * @param visualziation_time (in) Time value of current step, this is not necessarily the same
     * as the simulation time
     * @param visualization_step (in) Time step counter of current time step (does not have to be
     * continuous) this is not necessarily the same as the simulation time step counter
     */
    virtual void InitializeTimeStep(
        const double visualziation_time, const int visualization_step) = 0;

    /**
     * @brief Write all fields contained in the field data map to disk
     *
     * @param field_data_map (in) Map containing all field data vectors
     */
    virtual void WriteFieldDataToDisk(
        const std::map<std::string, visualization_vector_type_variant>& field_data_map) = 0;

    /**
     * @brief Write the full geometry, i.e., points, cells, faces and the respective connectivity to
     * disk
     *
     * @param point_coordinates (in) Vector containing the point coordinates
     * @param point_cell_connectivity (in) Vector containing the connectivity array
     * @param cell_offset (in) Vector containing the offsets in the connectivity array
     * @param cell_types (in) Vector containing the cell types
     * @param face_connectivity (in) Vector containing the face connectivity array
     * @param face_offset (in) Vector containing the offsets in the face connectivity array
     */
    virtual void WriteGeometryToDisk(const std::vector<double>& point_coordinates,
        const std::vector<IO::index_type>& point_cell_connectivity,
        const std::vector<IO::index_type>& cell_offset, const std::vector<uint8_t>& cell_types,
        const std::vector<IO::index_type>& face_connectivity,
        const std::vector<IO::index_type>& face_offset) = 0;

    /**
     * @brief Write a single point data vector to disk
     *
     * @param data (in) Vector containing the data
     * @param num_components_per_point (in) Number of components per point
     * @param name (in) Name of the data
     */
    virtual void WritePointDataVectorToDisk(const visualization_vector_type_variant& data,
        unsigned int num_components_per_point, const std::string& name) = 0;

    /**
     * @brief Write a single cell data vector to disk
     *
     * @param data (in) Vector containing the data
     * @param num_components_per_point (in) Number of components per cell
     * @param name (in) Name of the data
     */
    virtual void WriteCellDataVectorToDisk(const visualization_vector_type_variant& data,
        unsigned int num_components_per_point, const std::string& name) = 0;

    /**
     * @brief Finalize the write operations for the current time step
     */
    virtual void FinalizeTimeStep() = 0;

   protected:
    //! Visualization parameters
    IO::VisualizationParameters parameters_;

    // MPI communicator
    const Epetra_Comm& comm_;

    //! Specific name for the output data written by this object
    std::string visualization_data_name_;
  };
}  // namespace IO

FOUR_C_NAMESPACE_CLOSE

#endif