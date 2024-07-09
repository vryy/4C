/*----------------------------------------------------------------------*/
/*! \file

\brief VTU writer

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_io_vtu_writer.hpp"

#include "4C_io_pstream.hpp"
#include "4C_io_vtk_writer_base.hpp"
#include "4C_utils_exceptions.hpp"

#include <iomanip>
#include <iostream>
#include <sstream>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
VtuWriter::VtuWriter(unsigned int myrank, unsigned int num_processors,
    unsigned int max_number_timesteps_to_be_written,
    const std::string& path_existing_working_directory,
    const std::string& name_new_vtk_subdirectory, const std::string& geometry_name,
    const std::string& restart_name, const double restart_time, bool write_binary_output)
    : VtkWriterBase(myrank, num_processors, max_number_timesteps_to_be_written,
          path_existing_working_directory, name_new_vtk_subdirectory, geometry_name, restart_name,
          restart_time, write_binary_output)
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtuWriter::writer_string() const
{
  static std::string name("UnstructuredGrid");
  return name;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtuWriter::writer_opening_tag() const
{
  static std::string tag("<UnstructuredGrid>");
  return tag;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtuWriter::writer_p_opening_tag() const
{
  static std::string tag("<PUnstructuredGrid GhostLevel=\"0\">");
  return tag;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::vector<std::string>& VtuWriter::writer_p_piece_tags() const
{
  static std::vector<std::string> tags;
  tags.clear();

  for (size_t iproc = 0; iproc < numproc_; ++iproc)
  {
    std::stringstream stream;
    stream << "<Piece Source=\"" << filename_base_
           << get_part_of_file_name_indicating_processor_id(iproc) << ".vtu\"/>";
    tags.push_back(std::string(stream.str()));
  }
  return tags;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtuWriter::writer_suffix() const
{
  static std::string name(".vtu");
  return name;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtuWriter::writer_p_suffix() const
{
  static std::string name(".pvtu");
  return name;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void VtuWriter::write_geometry_unstructured_grid(const std::vector<double>& point_coordinates,
    const std::vector<Core::IO::index_type>& point_cell_connectivity,
    const std::vector<Core::IO::index_type>& cell_offset, const std::vector<uint8_t>& cell_types,
    const std::vector<Core::IO::index_type>& face_connectivity,
    const std::vector<Core::IO::index_type>& face_offset)
{
  // always assume 3D for now Todo maybe use this as template to allow for 2D case
  const unsigned int num_spatial_dimensions = 3;

  const unsigned int num_points = point_coordinates.size() / num_spatial_dimensions;

  const unsigned int num_cells = cell_types.size();

  // some sanity checks
  if (point_coordinates.size() % num_spatial_dimensions != 0)
    FOUR_C_THROW("VtuWriter assumes 3D point coordinates here! Extend to 2D if needed");

  if (cell_offset.size() != cell_types.size())
    FOUR_C_THROW(
        "VtuWriter: number of specified cell types does not equal number of "
        "specified index offsets");


  // step 0: tell the master file that we specify point coordinates
  /*----------------------------------------------------------------------*/
  if (myrank_ == 0)
  {
    throw_error_if_invalid_file_stream(currentmasterout_);

    currentmasterout_ << "    <PPoints>\n";
    currentmasterout_ << "      <PDataArray type=\"Float64\" NumberOfComponents=\""
                      << num_spatial_dimensions << "\"/>\n";
    currentmasterout_ << "    </PPoints>\n";
  }


  // step 1: write point coordinates into file
  /*----------------------------------------------------------------------*/
  throw_error_if_invalid_file_stream(currentout_);

  currentout_ << "    <Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << num_cells
              << "\" >\n"
              << "      <Points>\n"
              << "        <DataArray type=\"Float64\" NumberOfComponents=\""
              << num_spatial_dimensions << "\"";

  if (write_binary_output_)
  {
    currentout_ << " format=\"binary\">\n";
    LibB64::writeCompressedBlock(point_coordinates, currentout_);
  }
  else
  {
    currentout_ << " format=\"ascii\">\n";

    int counter = 1;
    for (std::vector<double>::const_iterator it = point_coordinates.begin();
         it != point_coordinates.end(); ++it)
    {
      currentout_ << std::setprecision(15) << std::scientific << *it;

      // single space between dimensions, new line upon completion of a point
      if (counter % num_spatial_dimensions != 0)
        currentout_ << " ";
      else
        currentout_ << '\n';

      counter++;
    }

    currentout_ << std::resetiosflags(std::ios::scientific);
  }


  currentout_ << "        </DataArray>\n"
              << "      </Points>\n\n";



  // step 2: write mesh-point topology into file
  /*----------------------------------------------------------------------*/
  currentout_ << "      <Cells>\n"
              << "        <DataArray type=\"Int32\" Name=\"connectivity\"";

  if (write_binary_output_)
  {
    currentout_ << " format=\"binary\">\n";
    LibB64::writeCompressedBlock(point_cell_connectivity, currentout_);
  }
  else
  {
    currentout_ << " format=\"ascii\">\n";

    for (std::vector<int32_t>::const_iterator it = point_cell_connectivity.begin();
         it != point_cell_connectivity.end(); ++it)
      currentout_ << *it << " ";
  }

  currentout_ << "\n        </DataArray>\n";



  // step 3: write indices where individual cells begin
  /*----------------------------------------------------------------------*/
  currentout_ << "        <DataArray type=\"Int32\" Name=\"offsets\"";

  if (write_binary_output_)
  {
    currentout_ << " format=\"binary\">\n";
    LibB64::writeCompressedBlock(cell_offset, currentout_);
  }
  else
  {
    currentout_ << " format=\"ascii\">\n";
    for (std::vector<int32_t>::const_iterator it = cell_offset.begin(); it != cell_offset.end();
         ++it)
      currentout_ << *it << " ";
  }

  currentout_ << "\n        </DataArray>\n";



  // step 4: write cell types
  /*----------------------------------------------------------------------*/
  currentout_ << "        <DataArray type=\"UInt8\" Name=\"types\"";
  if (write_binary_output_)
  {
    currentout_ << " format=\"binary\">\n";
    LibB64::writeCompressedBlock(cell_types, currentout_);
  }
  else
  {
    currentout_ << " format=\"ascii\">\n";
    for (std::vector<uint8_t>::const_iterator it = cell_types.begin(); it != cell_types.end(); ++it)
      currentout_ << (unsigned int)*it << " ";
  }
  currentout_ << "\n        </DataArray>\n";

  // step 5: write face data if required
  if (face_offset.size() > 0)
  {
    // Face connectivity
    currentout_ << R"(        <DataArray type="Int32" Name="faces")";
    if (write_binary_output_)
    {
      currentout_ << " format=\"binary\">\n";
      LibB64::writeCompressedBlock(face_connectivity, currentout_);
    }
    else
    {
      currentout_ << " format=\"ascii\">\n";
      for (const int value : face_connectivity) currentout_ << value << " ";
    }
    currentout_ << "\n        </DataArray>\n";

    // Face offsets
    currentout_ << R"(        <DataArray type="Int32" Name="faceoffsets")";
    if (write_binary_output_)
    {
      currentout_ << " format=\"binary\">\n";
      LibB64::writeCompressedBlock(face_offset, currentout_);
    }
    else
    {
      currentout_ << " format=\"ascii\">\n";
      for (const int value : face_offset) currentout_ << value << " ";
    }
    currentout_ << "\n        </DataArray>\n";
  }

  currentout_ << "      </Cells>\n\n";
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void VtuWriter::write_point_data_vector(const Core::IO::visualization_vector_type_variant& data,
    unsigned int num_components_per_point, const std::string& name)
{
  // start the point data section that will be written subsequently
  if (currentPhase_ == INIT)
  {
    throw_error_if_invalid_file_stream(currentout_);
    currentout_ << "  <PointData>\n";

    if (myrank_ == 0)
    {
      throw_error_if_invalid_file_stream(currentmasterout_);
      currentmasterout_ << "    <PPointData>\n";
    }

    currentPhase_ = POINTS;
  }

  if (currentPhase_ != POINTS)
    FOUR_C_THROW(
        "VtuWriter cannot write point data at this stage. Most likely, cell and "
        "point data fields are mixed. First, all point data needs to be written, "
        "then all cell data!");

  this->write_data_array(data, num_components_per_point, name);

  if (myrank_ == 0)
    Core::IO::cout(Core::IO::debug)
        << "\nVtuWriter: point data " << name << " written." << Core::IO::endl;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void VtuWriter::write_cell_data_vector(const Core::IO::visualization_vector_type_variant& data,
    unsigned int num_components_per_cell, const std::string& name)
{
  // if required, end the point data section
  if (currentPhase_ == POINTS)
  {
    throw_error_if_invalid_file_stream(currentout_);
    currentout_ << "  </PointData>\n";

    if (myrank_ == 0)
    {
      throw_error_if_invalid_file_stream(currentmasterout_);
      currentmasterout_ << "    </PPointData>\n";
    }
  }

  // start the cell data section that will be written subsequently
  if (currentPhase_ == INIT || currentPhase_ == POINTS)
  {
    throw_error_if_invalid_file_stream(currentout_);
    currentout_ << "  <CellData>\n";

    if (myrank_ == 0)
    {
      throw_error_if_invalid_file_stream(currentmasterout_);
      currentmasterout_ << "    <PCellData>\n";
    }

    currentPhase_ = CELLS;
  }

  if (currentPhase_ != CELLS)
    FOUR_C_THROW(
        "VtuWriter cannot write cell data at this stage. Most likely, cell and "
        "point data fields are mixed. First, all point data needs to be written, "
        "then all cell data!");

  this->write_data_array(data, num_components_per_cell, name);

  if (myrank_ == 0)
    Core::IO::cout(Core::IO::debug)
        << "\nVtuWriter: cell data " << name << " written." << Core::IO::endl;
}

FOUR_C_NAMESPACE_CLOSE
