/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief Write visualization output in vtk/vtu format at runtime

\level 3

\maintainer Maximilian Grill
*/
/*-----------------------------------------------------------------------------------------------*/

/* headers */
#include "runtime_vtu_writer.H"

#include "vtu_writer.H"

#include "../drt_lib/drt_dserror.H"


#include <iomanip>
#include <sstream>
#include <vector>

#include <boost/filesystem.hpp>

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
RuntimeVtuWriter::RuntimeVtuWriter() { vtu_writer_ = Teuchos::rcp(new VtuWriter()); }

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void RuntimeVtuWriter::Initialize(unsigned int myrank, unsigned int num_processors,
    unsigned int max_number_timesteps_to_be_written,
    const std::string& path_existing_output_directory, const std::string& simulation_name,
    const std::string& geometry_name, const std::string& restart_name, double restart_time,
    bool write_binary_output)
{
  vtu_writer_->Initialize(myrank, num_processors, max_number_timesteps_to_be_written,
      path_existing_output_directory, (simulation_name + "-vtk-files"), geometry_name, restart_name,
      restart_time, write_binary_output);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void RuntimeVtuWriter::SetupForNewTimeStepAndGeometry(
    double time, unsigned int timestep, const std::string& geometryname)
{
  vtu_writer_->ResetTimeAndTimeStep(time, timestep);

  vtu_writer_->ResetGeometryName(geometryname);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void RuntimeVtuWriter::ResetGeometry(const std::vector<double>& point_coordinates,
    const std::vector<uint8_t>& cell_types, const std::vector<int32_t>& cell_offset)
{
  point_coordinates_ = point_coordinates;
  cell_types_ = cell_types;
  cell_offset_ = cell_offset;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
// Todo template <typename T>, double, int,
void RuntimeVtuWriter::AppendVisualizationPointDataVector(const std::vector<double>& datavalues,
    const unsigned int num_components_per_point, const std::string& dataname)
{
  point_data_vectors_[dataname] = std::make_pair(datavalues, num_components_per_point);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
// Todo template <typename T>: double, (unsigned) int
void RuntimeVtuWriter::AppendVisualizationCellDataVector(const std::vector<double>& datavalues,
    const unsigned int num_components_per_point, const std::string& dataname)
{
  cell_data_vectors_[dataname] = std::make_pair(datavalues, num_components_per_point);
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
std::vector<double>& RuntimeVtuWriter::GetMutablePointCoordinateVector()
{
  return point_coordinates_;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
std::vector<uint8_t>& RuntimeVtuWriter::GetMutableCellTypeVector() { return cell_types_; }

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
std::vector<int32_t>& RuntimeVtuWriter::GetMutableCellOffsetVector() { return cell_offset_; }

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
std::map<std::string, std::vector<double>>& RuntimeVtuWriter::GetMutableFieldDataMap()
{
  return field_data_vectors_;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
std::map<std::string, std::pair<std::vector<double>, unsigned int>>&
RuntimeVtuWriter::GetMutablePointDataMap()
{
  return point_data_vectors_;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
std::map<std::string, std::pair<std::vector<double>, unsigned int>>&
RuntimeVtuWriter::GetMutableCellDataMap()
{
  return cell_data_vectors_;
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void RuntimeVtuWriter::WriteFiles()
{
  vtu_writer_->InitializeVtkFileStreamsForNewGeometryAndOrTimeStep();

  vtu_writer_->WriteVtkHeaders();

  vtu_writer_->WriteVtkFieldDataAndOrTimeAndOrCycle(field_data_vectors_);

  vtu_writer_->WriteGeometryUnstructuredGridContiguous(
      point_coordinates_, cell_offset_, cell_types_);

  WriteAllPointAndCellDataVectorsToFiles();

  vtu_writer_->WriteVtkFooters();
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void RuntimeVtuWriter::WriteAllPointAndCellDataVectorsToFiles()
{
  /* note: the implementation of the vtu writer expects to first write all point data and then
   *       write all cell data */

  // write all collected point data vectors to files
  typedef std::map<std::string, std::pair<std::vector<double>, unsigned int>> point_data_map;

  for (point_data_map::const_iterator point_data_iter = point_data_vectors_.begin();
       point_data_iter != point_data_vectors_.end(); ++point_data_iter)
  {
    vtu_writer_->WritePointDataVector(
        point_data_iter->second.first, point_data_iter->second.second, point_data_iter->first);
  }

  // write all collected cell data vectors to files
  typedef std::map<std::string, std::pair<std::vector<double>, unsigned int>> cell_data_map;

  for (cell_data_map::const_iterator cell_data_iter = cell_data_vectors_.begin();
       cell_data_iter != cell_data_vectors_.end(); ++cell_data_iter)
  {
    vtu_writer_->WriteCellDataVector(
        cell_data_iter->second.first, cell_data_iter->second.second, cell_data_iter->first);
  }
}

/*-----------------------------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------------------------*/
void RuntimeVtuWriter::WriteCollectionFileOfAllWrittenFiles(const std::string& collectionfilename)
{
  vtu_writer_->WriteVtkCollectionFileForAllWrittenMasterFiles(collectionfilename);
}
