/*----------------------------------------------------------------------*/
/*!

\brief VTP filter specialized for particle output

\level 2

\maintainer Martin Kronbichler
*-----------------------------------------------------------------------*/


#include "post_drt_vtp_writer.H"

#include <sstream>

#include "../drt_lib/drt_dserror.H"
#include "../post_drt_common/post_drt_common.H"
extern "C"
{
#include "../pss_full/pss_table.h"
}


namespace
{
  template <typename T>
  class make_vector
  {
   public:
    make_vector<T>& operator<<(const T& val)
    {
      data_.push_back(val);
      return *this;
    }
    operator std::vector<T>() const { return data_; }

   private:
    std::vector<T> data_;
  };
}  // namespace



PostVtpWriter::PostVtpWriter(PostField* field, const std::string& filename)
    : PostVtkWriter(field, filename)
{
}


const std::string& PostVtpWriter::WriterString() const
{
  static std::string name("PolyData");
  return name;
}

const std::string& PostVtpWriter::WriterOpeningTag() const
{
  static std::string tag("<PolyData>");
  return tag;
}

const std::string& PostVtpWriter::WriterPOpeningTag() const
{
  static std::string tag("<PPolyData GhostLevel=\"0\">");
  return tag;
}

const std::vector<std::string>& PostVtpWriter::WriterPPieceTags() const
{
  static std::vector<std::string> tags;
  tags.clear();
  for (size_t i = 0; i < numproc_; ++i)
  {
    std::stringstream stream;
    stream << "<Piece Source=\"" << filenamebase_ << "-" << i << ".vtp\"/>";
    tags.push_back(std::string(stream.str()));
  }
  return tags;
}

const std::string& PostVtpWriter::WriterSuffix() const
{
  static std::string name(".vtp");
  return name;
}

const std::string& PostVtpWriter::WriterPSuffix() const
{
  static std::string name(".pvtp");
  return name;
}

void PostVtpWriter::WriteGeo()
{
  if (timestep_ == 0 && myrank_ == 0)
  {
    int maxnodeid = field_->problem()->get_max_nodeid("particle");
    std::cout << "at most " << maxnodeid + 1
              << " particle(s) will be written in vtk PolyData format..." << std::endl;
  }
  // find maximum number of possible particles during simulation


  // read displacement field which is used for exclusively describing the geometry
  Teuchos::RCP<PostResult> result = Teuchos::rcp(new PostResult(field_));
  if (not result->next_result("displacement")) dserror("first displacement result is missing");

  // jump to the correct location in the data vector
  for (int i = 0; i < timestep_; ++i)
  {
    if (not result->next_result("displacement"))
    {
      // get current output step for info in dserror
      std::vector<int> solstep;
      std::vector<double> dummy;
      {
        PostResult result = PostResult(field_);
        result.get_result_timesandsteps(field_->name(), dummy, solstep);
      }
      int outputstep = solstep[timestep_];
      dserror("displacement for step %i is missing", outputstep);
    }
  }

  const Teuchos::RCP<Epetra_Vector> disp = result->read_result("displacement");

  // the number of nodes can be directly deduced from the state vector as
  // particle problems are always 3D
  const int nnodes = disp->MyLength() / 3;

  // fill in the coordinates
  std::vector<double> coordinates;
  coordinates.reserve(3 * nnodes);

  for (int n = 0; n < disp->MyLength(); ++n) coordinates.push_back((*disp)[n]);

  // write node coordinates into file
  currentout_ << "<Piece NumberOfPoints=\"" << nnodes << "\" >\n"
              << "  <Points>\n"
              << "    <DataArray type=\"Float64\" NumberOfComponents=\"3\"";

  if (!write_binary_output_)
  {
    currentout_ << " format=\"binary\">\n";
    LIBB64::writeCompressedBlock(coordinates, currentout_);
  }
  else
  {
    currentout_ << " format=\"ascii\">\n";
    for (std::vector<double>::const_iterator it = coordinates.begin(); it != coordinates.end();
         ++it)
      currentout_ << std::setprecision(15) << std::scientific << *it << " ";
    currentout_ << std::resetiosflags(std::ios::scientific);
  }
  currentout_ << "    </DataArray>\n"
              << "  </Points>\n\n";

  // avoid too much memory consumption -> clear coordinates vector now that we're done
  {
    std::vector<double> empty;
    empty.swap(coordinates);
  }

  if (myrank_ == 0)
  {
    currentmasterout_ << "    <PPoints>\n";
    currentmasterout_ << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n";
    currentmasterout_ << "    </PPoints>\n";
  }
}

void PostVtpWriter::WriteDofResultStep(std::ofstream& file, const Teuchos::RCP<Epetra_Vector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
    const std::string& groupname, const std::string& name, const int numdf, const int from,
    const bool fillzeros)
{
  if (myrank_ == 0 && timestep_ == 0) std::cout << "writing dof-based field " << name << std::endl;

  std::vector<double> solution;
  solution.reserve(data->MyLength());

  for (int lid = 0; lid < data->MyLength(); ++lid) solution.push_back((*data)[lid]);

  // start the scalar fields that will later be written
  if (currentPhase_ == INIT)
  {
    currentout_ << "  <PointData>\n";  // Scalars=\"scalars\">\n";
    if (myrank_ == 0)
    {
      currentmasterout_ << "    <PPointData>\n";  // Scalars=\"scalars\">\n";
    }
    currentPhase_ = POINTS;
  }

  if (currentPhase_ != POINTS)
    dserror(
        "Cannot write point data at this stage. Most likely cell and point data fields are mixed.");

  this->WriteSolutionVector(solution, 3, name, file);
}

void PostVtpWriter::WriteNodalResultStep(std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
    const std::string& groupname, const std::string& name, const int numdf)
{
  if (myrank_ == 0 && timestep_ == 0) std::cout << "writing node-based field " << name << std::endl;

  const int ncomponents = numdf;

  // count number of nodes for each processor
  const int nnodes = data->Map().NumMyElements();

  std::vector<double> solution;
  solution.reserve(ncomponents * nnodes);

  for (int lid = 0; lid < nnodes; ++lid)
  {
    for (int idf = 0; idf < numdf; ++idf)
    {
      Epetra_Vector* column = (*data)(idf);
      solution.push_back((*column)[lid]);
    }
  }  // loop over all nodes

  // start the scalar fields that will later be written
  if (currentPhase_ == INIT)
  {
    currentout_ << "  <PointData>\n";  // Scalars=\"scalars\">\n";
    if (myrank_ == 0)
    {
      currentmasterout_ << "    <PPointData>\n";  // Scalars=\"scalars\">\n";
    }
    currentPhase_ = POINTS;
  }

  if (currentPhase_ != POINTS)
    dserror(
        "Cannot write point data at this stage. Most likely cell and point data fields are mixed.");

  this->WriteSolutionVector(solution, ncomponents, name, file);
}

void PostVtpWriter::WriteElementResultStep(std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
    const std::string& groupname, const std::string& name, const int numdf, const int from)
{
  dserror("no element results expected for discretization based on nodes -> use vtu instead");
}
