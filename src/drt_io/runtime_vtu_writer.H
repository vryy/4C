/*-----------------------------------------------------------------------------------------------*/
/*!

\brief Write visualization output in vtk/vtu format at runtime

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/
#ifndef IO_RUNTIME_VTU_WRITER_H
#define IO_RUNTIME_VTU_WRITER_H

/*-----------------------------------------------------------------------------------------------*/
/* headers */

#include <Teuchos_RCP.hpp>

#include <string>
#include <vector>
#include <map>
#include <stdint.h>

/*-----------------------------------------------------------------------------------------------*/
/* forward declarations */
class VtuWriter;

namespace DRT
{
  class Discretization;
}

/*-----------------------------------------------------------------------------------------------*/
/* namespace */

/*!
 * \brief This object allows to write visualization output
 *        - in vtk/vtu format
 *        - at runtime
 *        - in parallel
 *        - binary-encoded
 *
 * \author grill
 * \date 03/17
 */
class RuntimeVtuWriter
{
 public:
  /// Constructor
  RuntimeVtuWriter();

  /// Destructor
  virtual ~RuntimeVtuWriter() {}

  /** \brief reset current simulation time and time step number
   *
   *  \author grill
   *  \date 03/17 */
  void Initialize(unsigned int myrank, unsigned int num_processors,
      unsigned int max_number_timesteps_to_be_written,
      const std::string& path_existing_output_directory, const std::string& simulation_name,
      const std::string& geometry_name, const std::string& restart_name, double restart_time,
      bool write_binary_output);

  /** \brief reset current simulation time and time step number and geometry name
   *
   *  \author grill
   *  \date 03/17 */
  void SetupForNewTimeStepAndGeometry(
      double time, unsigned int timestep, const std::string& geometryname);

  /** \brief reset geometry data for unstructured grid all at once
   *
   *  \author grill
   *  \date 03/17 */
  void ResetGeometry(const std::vector<double>& point_coordinates,
      const std::vector<uint8_t>& cell_types, const std::vector<int32_t>& cell_offset);

  /** \brief get mutable point coordinate vector for unstructured grid
   *
   *         this allows you to append geometry data in multiple steps
   *
   *  \author grill
   *  \date 03/17 */
  std::vector<double>& GetMutablePointCoordinateVector();

  /** \brief get mutable cell type vector for unstructured grid
   *
   *         this allows you to append geometry data in multiple steps
   *
   *  \author grill
   *  \date 03/17 */
  std::vector<uint8_t>& GetMutableCellTypeVector();

  /** \brief get mutable cell offset vector for unstructured grid
   *
   *         this allows you to append geometry data in multiple steps
   *
   *  \author grill
   *  \date 03/17 */
  std::vector<int32_t>& GetMutableCellOffsetVector();

  /**
   * \brief Return a mutable reference to the point data map.
   */
  std::map<std::string, std::pair<std::vector<double>, unsigned int>>& GetMutablePointDataMap();

  /**
   * \brief Return a mutable reference to the cell data map.
   */
  std::map<std::string, std::pair<std::vector<double>, unsigned int>>& GetMutableCellDataMap();

  /** \brief append visualization data vector with n components per point
   *
   *  \author grill
   *  \date 03/17 */
  // Todo template <typename T>, double, int,
  void AppendVisualizationPointDataVector(const std::vector<double>& datavalues,
      const unsigned int num_components_per_point, const std::string& dataname);

  /** \brief append visualization data vector with n components per cell
   *
   *  \author grill
   *  \date 03/17 */
  // Todo template <typename T>, double, int,
  void AppendVisualizationCellDataVector(const std::vector<double>& datavalues,
      const unsigned int num_components_per_cell, const std::string& dataname);


  /** \brief write all required VTK files to filesystem
   *
   *  \author grill
   *  \date 03/17 */
  void WriteFiles();

  /** \brief write a VTK collection file summarizing all written files of this object to filesystem
   *
   *  \author grill
   *  \date 03/17 */
  void WriteCollectionFileOfAllWrittenFiles(const std::string& collectionfilename);


 private:
  /** \brief write all stored data vectors to VTK files
   *
   *  \author grill
   *  \date 03/17 */
  void WriteAllPointAndCellDataVectorsToFiles();


 private:
  //! the actual vtu writer object
  Teuchos::RCP<VtuWriter> vtu_writer_;


  //! @name member variables required for the definition of an unstructured grid
  //! @{

  //! point coordinates with three components each
  std::vector<double> point_coordinates_;

  //! type of vtk cell (see list in ...)
  std::vector<uint8_t> cell_types_;

  //! index of first point of each cell (assumed to be continuous => no connectivity information
  //! required)
  std::vector<int32_t> cell_offset_;

  //! @}


  //! key: result name, entry: (solution data vector, num_components)
  std::map<std::string, std::pair<std::vector<double>, unsigned int>> point_data_vectors_;

  //! key: result name, entry: (solution data vector, num_components)
  std::map<std::string, std::pair<std::vector<double>, unsigned int>> cell_data_vectors_;
};

#endif
