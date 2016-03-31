/*----------------------------------------------------------------------*/
/*!
\file post_drt_vtk_writer.cpp

\brief VTK filter

\maintainer Martin Kronbichler
*-----------------------------------------------------------------------*/


#include "post_drt_vtk_writer.H"

#include <iomanip>
#include <boost/filesystem.hpp>

#include "../post_drt_common/post_drt_common.H"
#include "../pss_full/pss_cpp.h"
#include "../drt_lib/drt_discret.H"


// deactivate for ascii output. Only do this for debugging.
#define BIN_VTK_OUT


namespace LIBB64
{
  // functions taken from the libb64 project, http://sourceforge.net/projects/libb64
  //
  // libb64 is in the public domain
  namespace base64
  {
    typedef enum
    {
      step_A, step_B, step_C
    } base64_encodestep;

    typedef struct
    {
      base64_encodestep step;
      char result;
    } base64_encodestate;

    void base64_init_encodestate(base64_encodestate *state_in)
    {
      state_in->step = step_A;
      state_in->result = 0;
    }

    inline
    char base64_encode_value(char value_in)
    {
      static const char *encoding
        = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
      if (value_in > 63) return '=';
      return encoding[(int)value_in];
    }

    int base64_encode_block(const char *plaintext_in,
                            int length_in,
                            char *code_out,
                            base64_encodestate *state_in)
    {
      const char *plainchar = plaintext_in;
      const char *const plaintextend = plaintext_in + length_in;
      char *codechar = code_out;
      char result;
      char fragment;

      result = state_in->result;

      switch (state_in->step)
        {
          while (1)
            {
            case step_A:
              if (plainchar == plaintextend)
                {
                  state_in->result = result;
                  state_in->step = step_A;
                  return codechar - code_out;
                }
              fragment = *plainchar++;
              result = (fragment & 0x0fc) >> 2;
              *codechar++ = base64_encode_value(result);
              result = (fragment & 0x003) << 4;
            case step_B:
              if (plainchar == plaintextend)
                {
                  state_in->result = result;
                  state_in->step = step_B;
                  return codechar - code_out;
                }
              fragment = *plainchar++;
              result |= (fragment & 0x0f0) >> 4;
              *codechar++ = base64_encode_value(result);
              result = (fragment & 0x00f) << 2;
            case step_C:
              if (plainchar == plaintextend)
                {
                  state_in->result = result;
                  state_in->step = step_C;
                  return codechar - code_out;
                }
              fragment = *plainchar++;
              result |= (fragment & 0x0c0) >> 6;
              *codechar++ = base64_encode_value(result);
              result  = (fragment & 0x03f) >> 0;
              *codechar++ = base64_encode_value(result);
            }
        }
      /* control should not reach here */
      return codechar - code_out;
    }

    int base64_encode_blockend(char *code_out, base64_encodestate *state_in)
    {
      char *codechar = code_out;

      switch (state_in->step)
        {
        case step_B:
          *codechar++ = base64_encode_value(state_in->result);
          *codechar++ = '=';
          *codechar++ = '=';
          break;
        case step_C:
          *codechar++ = base64_encode_value(state_in->result);
          *codechar++ = '=';
          break;
        case step_A:
          break;
        }
      *codechar++ = '\0';

      return codechar - code_out;
    }
  }



  /**
   * Do a base64 encoding of the given data.
   *
   * The function allocates memory as necessary and returns a pointer to
   * it. The calling function must release this memory again.
   */
  char *
  encode_block (const char *data,
                const int   data_size)
  {
    base64::base64_encodestate state;
    base64::base64_init_encodestate(&state);

    char *encoded_data = new char[2*data_size+1];

    const int encoded_length_data
      = base64::base64_encode_block (data, data_size,
                                     encoded_data, &state);
    base64::base64_encode_blockend (encoded_data + encoded_length_data,
                                    &state);

    return encoded_data;
  }


  /**
   \brief Helper function to determine output file string from time step number
   */
  std::string int2string(const unsigned int i,
                         const unsigned int digits)
  {
    dsassert(i<std::pow(10,digits), "Invalid digits information");
    if (digits == 0 || digits > 9)
      return "invalid_digit";

    std::string digitstring (digits, '0');
    unsigned int divisor = 1;
    for (unsigned int d=0; d<digits; ++d, divisor *= 10)
      digitstring[digits-1-d] = '0' + (i%(divisor*10))/divisor;
    return digitstring;
  }

  /**
   \brief Helper function to determine output file string from time step number
   */
  unsigned int ndigits (unsigned int number)
  {
    // start numbering from 0, so need count digits based on number one less
    if (number > 1)
      number -= 1;
    unsigned int digits = 0;
    while (number > 0) {
      digits++;
      number /= 10;
    }
    return digits;
  }


} // end of LIBB64 namespace



VtkWriter::VtkWriter(PostField* field,
                     const std::string &filename)
:
    PostWriterBase(field, filename),
    currentPhase_(INIT),
    time_ (std::numeric_limits<double>::min()),
    timestep_ (0),
    cycle_ (std::numeric_limits<int>::max())
{}



void
VtkWriter::WriteVtkHeader ()
{
  if (!currentout_)
    dserror("Invalid output stream");

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
  if (time_ != std::numeric_limits<double>::min() ||
      cycle_ != std::numeric_limits<int>::max()) {
    currentout_ << "<FieldData>\n";

    if (time_ != std::numeric_limits<double>::min())
      currentout_ << "<DataArray type=\"Float32\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\">"
      << time_
      << "</DataArray>\n";
    if (cycle_ != std::numeric_limits<int>::max())
      currentout_ << "<DataArray type=\"Float32\" Name=\"CYCLE\" NumberOfTuples=\"1\" format=\"ascii\">"
      << cycle_
      << "</DataArray>\n";

    currentout_ << "</FieldData>\n";
  }

  // Also start master file on processor 0
  if (myrank_ == 0) {
    if (!currentmasterout_)
      dserror("Invalid output stream");

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



void
VtkWriter::WriteVtkFooter()
{
  if (!currentout_)
    dserror("Invalid output stream");

  // end the scalar fields
  if (currentPhase_ == POINTS) {
    currentout_ << "  </PointData>\n";
    if (myrank_ == 0) {
      currentmasterout_ << "    </PPointData>\n";
    }
    currentPhase_ = FINAL;
  } else if (currentPhase_ == CELLS) {
    currentout_ << "  </CellData>\n";
    if (myrank_ == 0) {
      currentmasterout_ << "    </PCellData>\n";
    }
    currentPhase_ = FINAL;
  } else {
    dserror("No data was written or writer was already in final phase.");
  }

  currentout_ << "</Piece>\n";

  currentout_ << "</" << this->WriterString() << ">\n";
  currentout_ << "</VTKFile>\n";

  currentout_ << std::flush;

  // Also start master file on processor 0
  typedef std::vector<std::string> pptags_type;
  const pptags_type& ppiecetags = this->WriterPPieceTags();
  if (myrank_ == 0) {
    if (!currentmasterout_)
      dserror("Invalid output stream");
    if (numproc_ != ppiecetags.size())
      dserror("Incorrect number of Pieces.");

    for (pptags_type::const_iterator it = ppiecetags.begin(); it != ppiecetags.end(); ++it)
      currentmasterout_ << "    " << *it << "\n";
    currentmasterout_ << "  </P" << this->WriterString() << ">\n";
    currentmasterout_ << "</VTKFile>\n";

    currentmasterout_ << std::flush;
  }
  return;
}



void
VtkWriter::WriteSpecialField (
      SpecialFieldInterface &special,
      PostResult& result,   ///< result group in the control file
      const ResultType  restype,
      const std::string &groupname,
      const std::vector<std::string> &fieldnames,
      const std::string &outinfo)
{
  // Vtk writes everything into the same file, so create to each output the
  // pointer to the same output writer
  std::vector<Teuchos::RCP<std::ofstream> > files(fieldnames.size());
  for (unsigned int i=0; i<fieldnames.size(); ++i)
    files[i] = Teuchos::rcp(&currentout_, false);

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
    dserror("Internal error when trying to identify output type %s",
            groupname.c_str());

  // jump to the correct location in the data vector. Some fields might only
  // be stored once, so need to catch that case as well
  bool once = false;
  for (int i=0; i<timestep_; ++i)
    if ( not activeresult.next_result(groupname) )
    {
      once = true;
      break;
    }
  if (once)
  {
    activeresult = PostResult(field_);
    activeresult.next_result(groupname);
  }

  std::map<std::string, std::vector<std::ofstream::pos_type> > resultfilepos;
  special(files,activeresult,resultfilepos,groupname,fieldnames);
}



void
VtkWriter::WriteSolutionVector (const std::vector<double> &solution,
                                const int ncomponents,
                                const std::string &name,
                                std::ofstream &file) const
{
  file << "    <DataArray type=\"Float64\" Name=\"" << name << "\"";
  if (ncomponents > 1)
    file << " NumberOfComponents=\"" << ncomponents << "\"";
#ifdef BIN_VTK_OUT
  file << " format=\"binary\">\n";
  LIBB64::writeCompressedBlock(solution, file);
#else
  file << " format=\"ascii\">\n";
  for (std::vector<double>::const_iterator it = solution.begin(); it != solution.end(); ++it)
    file << std::setprecision(15) << std::scientific << *it << " ";
  file << std::resetiosflags(std::ios::scientific);
#endif

  file << "    </DataArray>\n";

  std::ofstream & masterfile = const_cast<std::ofstream&>(currentmasterout_);
  if (myrank_ == 0) {
    masterfile << "      <PDataArray type=\"Float64\" Name=\"" << name << "\"";
    if (ncomponents > 1)
      masterfile << " NumberOfComponents=\"" << ncomponents << "\"";
    masterfile << " format=\"ascii\"/>\n";
  }
  return;
}


void
VtkWriter::WriteResult(const std::string groupname,
                       const std::string name,
                       const ResultType restype,
                       const int numdf,
                       const int from,
                       const bool fillzeros)
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
  for (int i=0; i<timestep_; ++i)
    if ( not result->next_result(groupname) )
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

  if ( not (field_->problem()->SpatialApproximation()=="Polynomial" or
            field_->problem()->SpatialApproximation()=="Meshfree" or
            field_->problem()->SpatialApproximation()=="HDG" or
            field_->problem()->SpatialApproximation()=="Nurbs"))
    dserror("Only polynomial or meshfree or Nurbs approximations can be written with the VTK filter");

  // need dummy structure that is required for the generic writer interface
  // but not needed by the vtk writer.
  std::map<std::string, std::vector<std::ofstream::pos_type> > dummy;

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



void
VtkWriter::WriteFiles(PostFilterBase &filter)
{
  PostResult result = PostResult(field_);

  // timesteps when the solution is written
  const std::vector<double> soltime = result.get_result_times(field_->name());
  ntdigits_ = LIBB64::ndigits(soltime.size());
  npdigits_ = LIBB64::ndigits(field_->discretization()->Comm().NumProc());
  std::vector<std::pair<double, std::string> > filenames;

  const std::string dirname = filename_ + "-files";
  boost::filesystem::create_directories(dirname);

  for (timestep_=0; timestep_<(int)soltime.size(); ++timestep_) {

    this->WriterPrepTimestep();

    {
      std::ostringstream tmpstream;
      tmpstream << field_->name() << "-" << std::setfill('0') << std::setw(ntdigits_) << timestep_;
      filenamebase_ = tmpstream.str();
    }

    time_ = soltime[timestep_];
    filenames.push_back(std::pair<double,std::string>(time_, filenamebase_+this->WriterPSuffix()));

    {
      std::ostringstream tmpstream;
      tmpstream << dirname << "/" << filenamebase_ << "-" << std::setfill('0') << std::setw(npdigits_) << myrank_ << this->WriterSuffix();
      currentout_.close();
      currentout_.open(tmpstream.str().c_str());
    }

    if (myrank_ == 0) {
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



void
VtkWriter::WriteFilesChangingGeom(PostFilterBase &filter)
{
  std::vector<int> solstep;
  std::vector<double> soltime;
  {
    PostResult result = PostResult(field_);
    result.get_result_timesandsteps(field_->name(), soltime, solstep);
  }

  unsigned int ntdigits = LIBB64::ndigits(soltime.size());
  unsigned int npdigits = LIBB64::ndigits(field_->discretization()->Comm().NumProc());
  std::vector<std::pair<double, std::string> > filenames;

  const std::string dirname = filename_ + "-files";
  boost::filesystem::create_directories(dirname);

  for (timestep_=0; timestep_<(int)soltime.size(); ++timestep_) {
    filenamebase_ = field_->name() + "-" + LIBB64::int2string(timestep_,ntdigits);
    time_ = soltime[timestep_];
    filenames.push_back(std::pair<double,std::string>(time_, filenamebase_+this->WriterPSuffix()));

    currentout_.close();
    currentout_.open((dirname + "/" + filenamebase_ + "-" +
        LIBB64::int2string(myrank_,npdigits) + this->WriterSuffix()).c_str());

    if (myrank_ == 0) {
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



void
VtkWriter::WriteVtkMasterFile(const std::vector<std::pair<double, std::string> > &filenames,
                              const std::string &dirname) const
{
  // finally, write a single masterfile
  if (myrank_ == 0) {
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

    masterfile
        << "<VTKFile type=\"Collection\" version=\"0.1\" ByteOrder=\"LittleEndian\">\n";
    masterfile << "  <Collection>\n";

    for (unsigned int i = 0; i < filenames.size(); ++i)
      masterfile << "    <DataSet timestep=\"" << filenames[i].first
                 << "\" group=\"\" part=\"0\" file=\"" << relative_dirname << "/"
                 << filenames[i].second
                 << "\"/>\n";

    masterfile << "  </Collection>\n";
    masterfile << "</VTKFile>\n";

    masterfile.flush();
  }

}
