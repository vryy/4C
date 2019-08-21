/*----------------------------------------------------------------------*/
/*! \file

\brief VTU writer

\level 2

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#include "../drt_io/io_pstream.H"

#include "vtu_writer.H"

#include "vtk_writer_base.H"

#include "../drt_lib/drt_dserror.H"

#include <sstream>
#include <iostream>
#include <iomanip>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
VtuWriter::VtuWriter() : VtkWriterBase()
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtuWriter::WriterString() const
{
  static std::string name("UnstructuredGrid");
  return name;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtuWriter::WriterOpeningTag() const
{
  static std::string tag("<UnstructuredGrid>");
  return tag;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtuWriter::WriterPOpeningTag() const
{
  static std::string tag("<PUnstructuredGrid GhostLevel=\"0\">");
  return tag;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::vector<std::string>& VtuWriter::WriterPPieceTags() const
{
  static std::vector<std::string> tags;
  tags.clear();

  for (size_t iproc = 0; iproc < numproc_; ++iproc)
  {
    std::stringstream stream;
    stream << "<Piece Source=\"" << filename_base_ << GetPartOfFileNameIndicatingProcessorId(iproc)
           << ".vtu\"/>";
    tags.push_back(std::string(stream.str()));
  }
  return tags;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtuWriter::WriterSuffix() const
{
  static std::string name(".vtu");
  return name;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtuWriter::WriterPSuffix() const
{
  static std::string name(".pvtu");
  return name;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void VtuWriter::WriteGeometryUnstructuredGridContiguous(
    const std::vector<double>& point_coordinates, const std::vector<int32_t>& cell_offset,
    const std::vector<uint8_t>& cell_types)
{
  // we assumed contiguous order of coordinates in this format, so fill the
  // vector with ascending numbers from 0 to num_points here

  // always assume 3D for now Todo use as template parameter
  const unsigned int num_spatial_dimensions = 3;

  const unsigned int num_points = point_coordinates.size() / num_spatial_dimensions;

  std::vector<int32_t> point_cell_connectivity;
  point_cell_connectivity.reserve(num_points);

  for (unsigned int i = 0; i < num_points; ++i) point_cell_connectivity.push_back(i);

  WriteGeometryUnstructuredGrid(
      point_coordinates, point_cell_connectivity, cell_offset, cell_types);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void VtuWriter::WriteGeometryUnstructuredGrid(const std::vector<double>& point_coordinates,
    const std::vector<int32_t>& point_cell_connectivity, const std::vector<int32_t>& cell_offset,
    const std::vector<uint8_t>& cell_types)
{
  // always assume 3D for now Todo maybe use this as template to allow for 2D case
  const unsigned int num_spatial_dimensions = 3;

  const unsigned int num_points = point_coordinates.size() / num_spatial_dimensions;

  const unsigned int num_cells = cell_types.size();

  // some sanity checks
  if (point_coordinates.size() % num_spatial_dimensions != 0)
    dserror("VtuWriter assumes 3D point coordinates here! Extend to 2D if needed");

  if (point_cell_connectivity.size() != num_points)
    dserror("VtuWriter: length of connectivity vector must equal number of points");

  if (cell_offset.size() != cell_types.size())
    dserror(
        "VtuWriter: number of specified cell types does not equal number of "
        "specified index offsets");


  // step 0: tell the master file that we specify point coordinates
  /*----------------------------------------------------------------------*/
  if (myrank_ == 0)
  {
    ThrowErrorIfInvalidFileStream(currentmasterout_);

    currentmasterout_ << "    <PPoints>\n";
    currentmasterout_ << "      <PDataArray type=\"Float64\" NumberOfComponents=\""
                      << num_spatial_dimensions << "\"/>\n";
    currentmasterout_ << "    </PPoints>\n";
  }


  // step 1: write point coordinates into file
  /*----------------------------------------------------------------------*/
  ThrowErrorIfInvalidFileStream(currentout_);

  currentout_ << "    <Piece NumberOfPoints=\"" << num_points << "\" NumberOfCells=\"" << num_cells
              << "\" >\n"
              << "      <Points>\n"
              << "        <DataArray type=\"Float64\" NumberOfComponents=\""
              << num_spatial_dimensions << "\"";

  if (write_binary_output_)
  {
    currentout_ << " format=\"binary\">\n";
    LIBB64::writeCompressedBlock(point_coordinates, currentout_);
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
    LIBB64::writeCompressedBlock(point_cell_connectivity, currentout_);
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
    LIBB64::writeCompressedBlock(cell_offset, currentout_);
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
    LIBB64::writeCompressedBlock(cell_types, currentout_);
  }
  else
  {
    currentout_ << " format=\"ascii\">\n";
    for (std::vector<uint8_t>::const_iterator it = cell_types.begin(); it != cell_types.end(); ++it)
      currentout_ << (unsigned int)*it << " ";
  }
  currentout_ << "\n        </DataArray>\n";

  currentout_ << "      </Cells>\n\n";
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void VtuWriter::WritePointDataVector(
    const std::vector<double>& data, unsigned int num_components_per_point, const std::string& name)
{
  // start the point data section that will be written subsequently
  if (currentPhase_ == INIT)
  {
    ThrowErrorIfInvalidFileStream(currentout_);
    currentout_ << "  <PointData>\n";

    if (myrank_ == 0)
    {
      ThrowErrorIfInvalidFileStream(currentmasterout_);
      currentmasterout_ << "    <PPointData>\n";
    }

    currentPhase_ = POINTS;
  }

  if (currentPhase_ != POINTS)
    dserror(
        "VtuWriter cannot write point data at this stage. Most likely, cell and "
        "point data fields are mixed. First, all point data needs to be written, "
        "then all cell data!");

  this->WriteDataArray(data, num_components_per_point, name);

  if (myrank_ == 0)
    IO::cout(IO::debug) << "\nVtuWriter: point data " << name << " written." << IO::endl;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void VtuWriter::WriteCellDataVector(
    const std::vector<double>& data, unsigned int num_components_per_cell, const std::string& name)
{
  // if required, end the point data section
  if (currentPhase_ == POINTS)
  {
    ThrowErrorIfInvalidFileStream(currentout_);
    currentout_ << "  </PointData>\n";

    if (myrank_ == 0)
    {
      ThrowErrorIfInvalidFileStream(currentmasterout_);
      currentmasterout_ << "    </PPointData>\n";
    }
  }

  // start the cell data section that will be written subsequently
  if (currentPhase_ == INIT || currentPhase_ == POINTS)
  {
    ThrowErrorIfInvalidFileStream(currentout_);
    currentout_ << "  <CellData>\n";

    if (myrank_ == 0)
    {
      ThrowErrorIfInvalidFileStream(currentmasterout_);
      currentmasterout_ << "    <PCellData>\n";
    }

    currentPhase_ = CELLS;
  }

  if (currentPhase_ != CELLS)
    dserror(
        "VtuWriter cannot write cell data at this stage. Most likely, cell and "
        "point data fields are mixed. First, all point data needs to be written, "
        "then all cell data!");

  this->WriteDataArray(data, num_components_per_cell, name);

  if (myrank_ == 0)
    IO::cout(IO::debug) << "\nVtuWriter: cell data " << name << " written." << IO::endl;
}
