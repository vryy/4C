/*----------------------------------------------------------------------*/
/*! \file

\brief VTK filter

\level 2

*-----------------------------------------------------------------------*/


#include "4C_post_vtk_writer.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_io_legacy_table.hpp"
#include "4C_post_common.hpp"

#include <filesystem>
#include <iomanip>

FOUR_C_NAMESPACE_OPEN

// deactivate for ascii output. Only do this for debugging.
// #define BIN_VTK_OUT

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



void PostVtkWriter::write_vtk_header()
{
  if (!currentout_) FOUR_C_THROW("Invalid output stream");

  // TODO: might need BigEndian on some systems
  const std::string byteorder = "LittleEndian";

  currentout_ << "<?xml version=\"1.0\" ?> \n";
  currentout_ << "<!-- \n";
  currentout_ << "# vtk DataFile Version 3.0\n";
  currentout_ << "-->\n";
  currentout_ << "<VTKFile type=\"" << this->writer_string() << "\" version=\"0.1\"";
  currentout_ << " compressor=\"vtkZLibDataCompressor\"";
  currentout_ << " byte_order=\"" << byteorder << "\"";
  currentout_ << ">\n";
  currentout_ << this->writer_opening_tag() << "\n";

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
    if (!currentmasterout_) FOUR_C_THROW("Invalid output stream");

    currentmasterout_ << "<?xml version=\"1.0\" ?> \n";
    currentmasterout_ << "<!-- \n";
    currentmasterout_ << "# vtk DataFile Version 3.0\n";
    currentmasterout_ << "-->\n";
    currentmasterout_ << "<VTKFile type=\"P" << this->writer_string() << "\" version=\"0.1\"";
    currentmasterout_ << " byte_order=\"" << byteorder << "\"";
    currentmasterout_ << ">\n";
    currentmasterout_ << "  " << this->writer_p_opening_tag() << "\n";
  }

  currentPhase_ = INIT;
}



void PostVtkWriter::write_vtk_footer()
{
  if (!currentout_) FOUR_C_THROW("Invalid output stream");

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
    FOUR_C_THROW("No data was written or writer was already in final phase.");
  }

  currentout_ << "</Piece>\n";

  currentout_ << "</" << this->writer_string() << ">\n";
  currentout_ << "</VTKFile>\n";

  currentout_ << std::flush;

  // Also start master file on processor 0
  typedef std::vector<std::string> pptags_type;
  const pptags_type &ppiecetags = this->writer_p_piece_tags();
  if (myrank_ == 0)
  {
    if (!currentmasterout_) FOUR_C_THROW("Invalid output stream");
    if (numproc_ != ppiecetags.size()) FOUR_C_THROW("Incorrect number of Pieces.");

    for (pptags_type::const_iterator it = ppiecetags.begin(); it != ppiecetags.end(); ++it)
      currentmasterout_ << "    " << *it << "\n";
    currentmasterout_ << "  </P" << this->writer_string() << ">\n";
    currentmasterout_ << "</VTKFile>\n";

    currentmasterout_ << std::flush;
  }
  return;
}



void PostVtkWriter::write_special_field(SpecialFieldInterface &special,
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
  if (!foundit)
    FOUR_C_THROW("Internal error when trying to identify output type %s", groupname.c_str());

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



void PostVtkWriter::write_solution_vector(const std::vector<double> &solution,
    const int ncomponents, const std::string &name, std::ofstream &file) const
{
  using namespace FourC;

  file << "    <DataArray type=\"Float64\" Name=\"" << name << "\"";
  if (ncomponents > 1) file << " NumberOfComponents=\"" << ncomponents << "\"";
  if (write_binary_output_)
  {
    file << " format=\"binary\">\n";
    LibB64::writeCompressedBlock(solution, file);
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


void PostVtkWriter::write_result(const std::string groupname, const std::string name,
    const ResultType restype, const int numdf, const int from, const bool fillzeros)
{
  using namespace FourC;

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
  if (not(field_->problem()->spatial_approximation_type() ==
              Core::FE::ShapeFunctionType::polynomial or
          field_->problem()->spatial_approximation_type() == Core::FE::ShapeFunctionType::hdg or
          field_->problem()->spatial_approximation_type() == Core::FE::ShapeFunctionType::nurbs))
    FOUR_C_THROW(
        "Undefined spatial approximation type or the VTK filter is not yet implemented for the "
        "given type.");
  // need dummy structure that is required for the generic writer interface
  // but not needed by the vtk writer.
  std::map<std::string, std::vector<std::ofstream::pos_type>> dummy;

  switch (restype)
  {
    case dofbased:
    {
      const Teuchos::RCP<Epetra_Vector> data = result->read_result(groupname);
      this->write_dof_result_step(
          currentout_, data, dummy, groupname, name, numdf, from, fillzeros);
      break;
    }
    case nodebased:
    {
      const Teuchos::RCP<Epetra_MultiVector> data = result->read_multi_result(groupname);
      this->write_nodal_result_step(currentout_, data, dummy, groupname, name, numdf);
      break;
    }
    case elementbased:
    {
      const Teuchos::RCP<Epetra_MultiVector> data = result->read_multi_result(groupname);
      this->write_element_result_step(currentout_, data, dummy, groupname, name, numdf, from);
      break;
    }
    default:
    {
      FOUR_C_THROW("Result type not yet implemented");
      break;
    }
  }
}



void PostVtkWriter::write_files(PostFilterBase &filter)
{
  using namespace FourC;

  PostResult result = PostResult(field_);

  // timesteps when the solution is written
  const std::vector<double> soltime = result.get_result_times(field_->name());
  ntdigits_ = LibB64::ndigits(soltime.size());
  npdigits_ = LibB64::ndigits(field_->discretization()->get_comm().NumProc());
  std::vector<std::pair<double, std::string>> filenames;

  const std::string dirname = filename_ + "-files";
  std::filesystem::create_directories(dirname);

  for (timestep_ = 0; timestep_ < (int)soltime.size(); ++timestep_)
  {
    this->writer_prep_timestep();

    {
      std::ostringstream tmpstream;
      tmpstream << field_->name() << "-" << std::setfill('0') << std::setw(ntdigits_) << timestep_;
      filenamebase_ = tmpstream.str();
    }

    time_ = soltime[timestep_];
    filenames.push_back(
        std::pair<double, std::string>(time_, filenamebase_ + this->writer_p_suffix()));

    {
      std::ostringstream tmpstream;
      tmpstream << dirname << "/" << filenamebase_ << "-" << std::setfill('0')
                << std::setw(npdigits_) << myrank_ << this->writer_suffix();
      currentout_.close();
      currentout_.open(tmpstream.str().c_str());
    }

    if (myrank_ == 0)
    {
      currentmasterout_.close();
      currentmasterout_.open((dirname + "/" + filenamebase_ + this->writer_p_suffix()).c_str());
    }

    write_vtk_header();

    write_geo();

    filter.write_all_results(field_);

    write_vtk_footer();
  }

  write_vtk_master_file(filenames, dirname);
}


void PostVtkWriter::write_vtk_master_file(
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

FOUR_C_NAMESPACE_CLOSE
