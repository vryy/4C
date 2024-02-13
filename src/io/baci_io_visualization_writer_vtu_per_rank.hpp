/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization data to disk in the pvtu format, where each MPI rank writes its own file

\level 0

*/
/*-----------------------------------------------------------------------------------------------*/

#ifndef BACI_IO_VISUALIZATION_WRITER_VTU_PER_RANK_HPP
#define BACI_IO_VISUALIZATION_WRITER_VTU_PER_RANK_HPP

#include "baci_config.hpp"

#include "baci_io_visualization_writer_base.hpp"
#include "baci_io_vtu_writer.hpp"

BACI_NAMESPACE_OPEN

namespace IO
{
  class VisualizationWriterVtuPerRank : public VisualizationWriterBase
  {
   public:
    /**
     * @brief Default constructor
     */
    VisualizationWriterVtuPerRank(const IO::VisualizationParameters& parameters,
        const Epetra_Comm& comm, std::string visualization_data_name);

    /**
     * @brief Default destructor
     */
    ~VisualizationWriterVtuPerRank() override = default;

    /**
     * @brief Initialize the current time step (derived)
     */
    void InitializeTimeStep(const double visualziation_time, const int visualization_step) override;

    /**
     * @brief Write all fields contained in the field data map to disk (derived)
     */
    void WriteFieldDataToDisk(
        const std::map<std::string, visualization_vector_type_variant>& field_data_map) override;

    /**
     * @brief Write the full geometry, i.e., points, cells, faces and the respective connectivity to
     * disk (derived)
     */
    void WriteGeometryToDisk(const std::vector<double>& point_coordinates,
        const std::vector<IO::index_type>& point_cell_connectivity,
        const std::vector<IO::index_type>& cell_offset, const std::vector<uint8_t>& cell_types,
        const std::vector<IO::index_type>& face_connectivity,
        const std::vector<IO::index_type>& face_offset) override;

    /**
     * @brief Write a single point data vector to disk (derived)
     */
    void WritePointDataVectorToDisk(const visualization_vector_type_variant& data,
        unsigned int num_components_per_point, const std::string& name) override;

    /**
     * @brief Write a single cell data vector to disk (derived)
     */
    void WriteCellDataVectorToDisk(const visualization_vector_type_variant& data,
        unsigned int num_components_per_point, const std::string& name) override;

    /**
     * @brief Finalize the write operations for the current time step (derived)
     */
    void FinalizeTimeStep() override;

    //! At the moment, simply point to the "old" vtu writer here, the functionality of this
    //! object will be implemented in this class in the future
    VtuWriter vtu_writer_;
  };
}  // namespace IO

BACI_NAMESPACE_CLOSE

#endif