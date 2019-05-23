/*----------------------------------------------------------------------*/
/*!

\brief VTK filter

\level 2

\maintainer Martin Kronbichler
*-----------------------------------------------------------------------*/


#include "post_drt_vtk_writer.H"

#include <iomanip>
#include <boost/filesystem.hpp>

#include "../post_drt_common/post_drt_common.H"
#include "../pss_full/pss_cpp.h"
#include "../drt_lib/drt_discret.H"


// deactivate for ascii output. Only do this for debugging.
//#define BIN_VTK_OUT

PostVtkWriter::PostVtkWriter(PostField *field, const std::string &filename)
    : PostWriterBase(field, filename),
      currentPhase_(INIT),
      time_(std::numeric_limits<double>::min()),
      timestep_(0),
      cycle_(std::numeric_limits<int>::max())
{
  if (field->problem()->outputtype() == "bin")
    write_binary_output_ = true;  ///< write binary output
  else
    write_binary_output_ = false;
}



void PostVtkWriter::WriteVtkHeader()
{
  if (!currentout_) dserror("Invalid output stream");

  // TODO: might need BigEndian on some systems
  const std::string byteorder = "LittleEndian";

  currentout_ << "<?xml version=\"1.0\" ?> \n";
  currentout_ << "<!-- \n";
  currentout_ << "# vtk DataFile Version 3.0\n";
  currentout_ << "-->\n";
  currentout_ << "<VTKFile type=\"" << this->WriterString() << "\" version=\"0.1\"";
  currentout_ << " compressor=\"vtkZLibDataCompressor\"";
  currentout_ << " byte_order=\"" << byteorder << "\"";
  currentout_ << ">\n";
  currentout_ << this->WriterOpeningTag() << "\n";

  // Print output time and cycle
  if (time_ != std::numeric_limits<double>::min() || cycle_ != std::numeric_limits<int>::max())
  {
    currentout_ << "<FieldData>\n";

    if (time_ != std::numeric_limits<double>::min())
      currentout_
          << "<DataArray type=\"Float32\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\">"
          << time_ << "</DataArray>\n";
    if (cycle_ != std::numeric_limits<int>::max())
      currentout_
          << "<DataArray type=\"Float32\" Name=\"CYCLE\" NumberOfTuples=\"1\" format=\"ascii\">"
          << cycle_ << "</DataArray>\n";

    currentout_ << "</FieldData>\n";
  }

  // Also start master file on processor 0
  if (myrank_ == 0)
  {
    if (!currentmasterout_) dserror("Invalid output stream");

    currentmasterout_ << "<?xml version=\"1.0\" ?> \n";
    currentmasterout_ << "<!-- \n";
    currentmasterout_ << "# vtk DataFile Version 3.0\n";
    currentmasterout_ << "-->\n";
    currentmasterout_ << "<VTKFile type=\"P" << this->WriterString() << "\" version=\"0.1\"";
    currentmasterout_ << " byte_order=\"" << byteorder << "\"";
    currentmasterout_ << ">\n";
    currentmasterout_ << "  " << this->WriterPOpeningTag() << "\n";
  }

  currentPhase_ = INIT;
}



void PostVtkWriter::WriteVtkFooter()
{
  if (!currentout_) dserror("Invalid output stream");

  // end the scalar fields
  if (currentPhase_ == POINTS)
  {
    currentout_ << "  </PointData>\n";
    if (myrank_ == 0)
    {
      currentmasterout_ << "    </PPointData>\n";
    }
    currentPhase_ = FINAL;
  }
  else if (currentPhase_ == CELLS)
  {
    currentout_ << "  </CellData>\n";
    if (myrank_ == 0)
    {
      currentmasterout_ << "    </PCellData>\n";
    }
    currentPhase_ = FINAL;
  }
  else
  {
    dserror("No data was written or writer was already in final phase.");
  }

  currentout_ << "</Piece>\n";

  currentout_ << "</" << this->WriterString() << ">\n";
  currentout_ << "</VTKFile>\n";

  currentout_ << std::flush;

  // Also start master file on processor 0
  typedef std::vector<std::string> pptags_type;
  const pptags_type &ppiecetags = this->WriterPPieceTags();
  if (myrank_ == 0)
  {
    if (!currentmasterout_) dserror("Invalid output stream");
    if (numproc_ != ppiecetags.size()) dserror("Incorrect number of Pieces.");

    for (pptags_type::const_iterator it = ppiecetags.begin(); it != ppiecetags.end(); ++it)
      currentmasterout_ << "    " << *it << "\n";
    currentmasterout_ << "  </P" << this->WriterString() << ">\n";
    currentmasterout_ << "</VTKFile>\n";

    currentmasterout_ << std::flush;
  }
  return;
}



void PostVtkWriter::WriteSpecialField(SpecialFieldInterface &special,
    PostResult &result,  ///< result group in the control file
    const ResultType restype, const std::string &groupname,
    const std::vector<std::string> &fieldnames, const std::string &outinfo)
{
  // Vtk writes everything into the same file, so create to each output the
  // pointer to the same output writer
  std::vector<Teuchos::RCP<std::ofstream>> files(fieldnames.size());
  for (unsigned int i = 0; i < fieldnames.size(); ++i) files[i] = Teuchos::rcp(&currentout_, false);

  bool foundit = false;
  PostResult activeresult(result.field());
  while (activeresult.next_result(groupname))
  {
    if (map_has_map(activeresult.group(), groupname.c_str()))
    {
      foundit = true;
      break;
    }
  }
  // should always find the correct result
  if (!foundit) dserror("Internal error when trying to identify output type %s", groupname.c_str());

  // jump to the correct location in the data vector. Some fields might only
  // be stored once, so need to catch that case as well
  bool once = false;
  for (int i = 0; i < timestep_; ++i)
    if (not activeresult.next_result(groupname))
    {
      once = true;
      break;
    }
  if (once)
  {
    activeresult = PostResult(field_);
    activeresult.next_result(groupname);
  }

  std::map<std::string, std::vector<std::ofstream::pos_type>> resultfilepos;
  special(files, activeresult, resultfilepos, groupname, fieldnames);
}



void PostVtkWriter::WriteSolutionVector(const std::vector<double> &solution, const int ncomponents,
    const std::string &name, std::ofstream &file) const
{
  file << "    <DataArray type=\"Float64\" Name=\"" << name << "\"";
  if (ncomponents > 1) file << " NumberOfComponents=\"" << ncomponents << "\"";
  if (write_binary_output_)
  {
    file << " format=\"binary\">\n";
    LIBB64::writeCompressedBlock(solution, file);
  }
  else
  {
    file << " format=\"ascii\">\n";
    int counter = 1;
    for (std::vector<double>::const_iterator it = solution.begin(); it != solution.end(); ++it)
    {
      file << std::setprecision(15) << std::scientific << *it << " ";
      if (counter % ncomponents == 0) file << '\n';
      counter++;
    }
    file << std::resetiosflags(std::ios::scientific);
  }

  file << "    </DataArray>\n";

  std::ofstream &masterfile = const_cast<std::ofstream &>(currentmasterout_);
  if (myrank_ == 0)
  {
    masterfile << "      <PDataArray type=\"Float64\" Name=\"" << name << "\"";
    if (ncomponents > 1) masterfile << " NumberOfComponents=\"" << ncomponents << "\"";
    masterfile << " format=\"ascii\"/>\n";
  }
  return;
}


void PostVtkWriter::WriteResult(const std::string groupname, const std::string name,
    const ResultType restype, const int numdf, const int from, const bool fillzeros)
{
  Teuchos::RCP<PostResult> result = Teuchos::rcp(new PostResult(field_));
  // only write results which exist in the first result step
  bool foundit = false;
  if (result->next_result())
  {
    if (map_has_map(result->group(), groupname.c_str()))
    {
      foundit = true;
    }
  }
  if (!foundit) return;

  // jump to the correct location in the data vector. Some fields might only
  // be stored once, so need to catch that case as well
  bool once = false;
  for (int i = 0; i < timestep_; ++i)
    if (not result->next_result(groupname))
    {
      once = true;
      break;
    }

  if (once)
  {
    // recreate PostResult, go one step and throw away the old one
    result = Teuchos::rcp(new PostResult(field_));
    result->next_result(groupname);
  }

  if (not(field_->problem()->SpatialApproximation() == "Polynomial" or
          field_->problem()->SpatialApproximation() == "Meshfree" or
          field_->problem()->SpatialApproximation() == "HDG" or
          field_->problem()->SpatialApproximation() == "Nurbs"))
    dserror(
        "Only polynomial or meshfree or Nurbs approximations can be written with the VTK filter");

  // need dummy structure that is required for the generic writer interface
  // but not needed by the vtk writer.
  std::map<std::string, std::vector<std::ofstream::pos_type>> dummy;

  switch (restype)
  {
    case dofbased:
    {
      const Teuchos::RCP<Epetra_Vector> data = result->read_result(groupname);
      this->WriteDofResultStep(currentout_, data, dummy, groupname, name, numdf, from, fillzeros);
      break;
    }
    case nodebased:
    {
      const Teuchos::RCP<Epetra_MultiVector> data = result->read_multi_result(groupname);
      this->WriteNodalResultStep(currentout_, data, dummy, groupname, name, numdf);
      break;
    }
    case elementbased:
    {
      const Teuchos::RCP<Epetra_MultiVector> data = result->read_multi_result(groupname);
      this->WriteElementResultStep(currentout_, data, dummy, groupname, name, numdf, from);
      break;
    }
    default:
    {
      dserror("Result type not yet implemented");
      break;
    }
  }
}



void PostVtkWriter::WriteFiles(PostFilterBase &filter)
{
  PostResult result = PostResult(field_);

  // timesteps when the solution is written
  const std::vector<double> soltime = result.get_result_times(field_->name());
  ntdigits_ = LIBB64::ndigits(soltime.size());
  npdigits_ = LIBB64::ndigits(field_->discretization()->Comm().NumProc());
  std::vector<std::pair<double, std::string>> filenames;

  const std::string dirname = filename_ + "-files";
  boost::filesystem::create_directories(dirname);

  for (timestep_ = 0; timestep_ < (int)soltime.size(); ++timestep_)
  {
    this->WriterPrepTimestep();

    {
      std::ostringstream tmpstream;
      tmpstream << field_->name() << "-" << std::setfill('0') << std::setw(ntdigits_) << timestep_;
      filenamebase_ = tmpstream.str();
    }

    time_ = soltime[timestep_];
    filenames.push_back(
        std::pair<double, std::string>(time_, filenamebase_ + this->WriterPSuffix()));

    {
      std::ostringstream tmpstream;
      tmpstream << dirname << "/" << filenamebase_ << "-" << std::setfill('0')
                << std::setw(npdigits_) << myrank_ << this->WriterSuffix();
      currentout_.close();
      currentout_.open(tmpstream.str().c_str());
    }

    if (myrank_ == 0)
    {
      currentmasterout_.close();
      currentmasterout_.open((dirname + "/" + filenamebase_ + this->WriterPSuffix()).c_str());
    }

    WriteVtkHeader();

    WriteGeo();

    filter.WriteAllResults(field_);

    WriteVtkFooter();
  }

  WriteVtkMasterFile(filenames, dirname);
}



void PostVtkWriter::WriteFilesChangingGeom(PostFilterBase &filter)
{
  std::vector<int> solstep;
  std::vector<double> soltime;
  {
    PostResult result = PostResult(field_);
    result.get_result_timesandsteps(field_->name(), soltime, solstep);
  }

  unsigned int ntdigits = LIBB64::ndigits(soltime.size());
  unsigned int npdigits = LIBB64::ndigits(field_->discretization()->Comm().NumProc());
  std::vector<std::pair<double, std::string>> filenames;

  const std::string dirname = filename_ + "-files";
  boost::filesystem::create_directories(dirname);

  for (timestep_ = 0; timestep_ < (int)soltime.size(); ++timestep_)
  {
    filenamebase_ = field_->name() + "-" + LIBB64::int2string(timestep_, ntdigits);
    time_ = soltime[timestep_];
    filenames.push_back(
        std::pair<double, std::string>(time_, filenamebase_ + this->WriterPSuffix()));

    currentout_.close();
    currentout_.open((dirname + "/" + filenamebase_ + "-" + LIBB64::int2string(myrank_, npdigits) +
                      this->WriterSuffix())
                         .c_str());

    if (myrank_ == 0)
    {
      currentmasterout_.close();
      currentmasterout_.open((dirname + "/" + filenamebase_ + this->WriterPSuffix()).c_str());
    }

    int fieldpos = field_->field_pos();
    std::string fieldname = field_->name();
    field_->problem()->re_read_mesh(fieldpos, fieldname, solstep[timestep_]);

    WriteVtkHeader();

    WriteGeo();

    filter.WriteAllResults(field_);

    WriteVtkFooter();
  }

  WriteVtkMasterFile(filenames, dirname);
}



void PostVtkWriter::WriteVtkMasterFile(
    const std::vector<std::pair<double, std::string>> &filenames, const std::string &dirname) const
{
  // finally, write a single masterfile
  if (myrank_ == 0)
  {
    size_t pos = dirname.find_last_of("/");
    if (pos == dirname.npos)
      pos = 0ul;
    else
      pos++;
    std::string relative_dirname(dirname.substr(pos));

    std::ofstream masterfile((filename_ + "-" + field_->name() + ".pvd").c_str());

    masterfile << "<?xml version=\"1.0\"?>\n";

    masterfile << "<!--\n";
    masterfile << "-->\n";

    masterfile << "<VTKFile type=\"Collection\" version=\"0.1\" ByteOrder=\"LittleEndian\">\n";
    masterfile << "  <Collection>\n";

    for (unsigned int i = 0; i < filenames.size(); ++i)
      masterfile << "    <DataSet timestep=\"" << filenames[i].first
                 << "\" group=\"\" part=\"0\" file=\"" << relative_dirname << "/"
                 << filenames[i].second << "\"/>\n";

    masterfile << "  </Collection>\n";
    masterfile << "</VTKFile>\n";

    masterfile.flush();
  }
}
