/*----------------------------------------------------------------------*/
/*!
\file vtp_writer.cpp

\brief VTP writer

\level 2

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/

#include "../drt_io/io_pstream.H"

#include "vtp_writer.H"

#include "vtk_writer_base.H"

#include "../drt_lib/drt_dserror.H"

#include <sstream>
#include <iostream>
#include <iomanip>


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
VtpWriter::VtpWriter() : VtkWriterBase()
{
  // empty constructor
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtpWriter::WriterString() const
{
  static std::string name("PolyData");
  return name;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtpWriter::WriterOpeningTag() const
{
  static std::string tag("<PolyData>");
  return tag;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtpWriter::WriterPOpeningTag() const
{
  static std::string tag("<PPolyData GhostLevel=\"0\">");
  return tag;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::vector<std::string>& VtpWriter::WriterPPieceTags() const
{
  static std::vector<std::string> tags;
  tags.clear();

  for (size_t iproc = 0; iproc < numproc_; ++iproc)
  {
    std::stringstream stream;
    stream << "<Piece Source=\"" << filename_base_ << GetPartOfFileNameIndicatingProcessorId(iproc)
           << ".vtp\"/>";
    tags.push_back(std::string(stream.str()));
  }
  return tags;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtpWriter::WriterSuffix() const
{
  static std::string name(".vtp");
  return name;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string& VtpWriter::WriterPSuffix() const
{
  static std::string name(".pvtp");
  return name;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void VtpWriter::WriteGeometryPolyData(const std::vector<double>& point_coordinates)
{
  // always assume 3D for now Todo maybe use this as template to allow for 2D case
  const unsigned int num_spatial_dimensions = 3;

  const unsigned int num_points = point_coordinates.size() / num_spatial_dimensions;

  // some sanity checks
  if (point_coordinates.size() % num_spatial_dimensions != 0)
    dserror("VtpWriter assumes 3D point coordinates here! Extend to 2D if needed");


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

  currentout_ << "    <Piece NumberOfPoints=\"" << num_points << "\" >\n"
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
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void VtpWriter::WritePointDataVector(
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
        "VtpWriter cannot write point data at this stage. Most likely, cell and "
        "point data fields are mixed. First, all point data needs to be written, "
        "then all cell data!");

  this->WriteDataArray(data, num_components_per_point, name);

  if (myrank_ == 0)
    IO::cout(IO::debug) << "\nVtpWriter: point data " << name << " written." << IO::endl;
}
