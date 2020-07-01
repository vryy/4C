/*----------------------------------------------------------------------------*/
/** \file

\brief Input/output routines for the multi discretization wrapper

\level 3

*/
/*----------------------------------------------------------------------------*/

#include "xstr_multi_io.H"
#include "xstr_multi_discretization_wrapper.H"

#include "../drt_io/io_control.H"

#include "../pss_full/pss_cpp.h"

#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XSTR::IO::MultiDiscretizationReader::MultiDiscretizationReader(
    const Teuchos::RCP<XSTR::MultiDiscretizationWrapper>& dis_wrapper,
    const Teuchos::RCP<::IO::InputControl>& input, int step)
    : ::IO::DiscretizationReader(), dis_wrapper_(Teuchos::rcp(dis_wrapper.get(), false))
{
  XSTR::MultiDiscretizationWrapper::XDisMap::const_iterator cit;
  for (cit = dis_wrapper_->DiscretMap().begin(); cit != dis_wrapper_->DiscretMap().end(); ++cit)
  {
    readers_[cit->first] =
        Teuchos::rcp(new XSTR::IO::DiscretizationReader(cit->second, input, step));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XSTR::IO::MultiDiscretizationReader::MultiDiscretizationReader(
    const Teuchos::RCP<XSTR::MultiDiscretizationWrapper>& dis_wrapper, int step)
    : ::IO::DiscretizationReader(), dis_wrapper_(dis_wrapper)
{
  XSTR::MultiDiscretizationWrapper::XDisMap::const_iterator cit;
  for (cit = dis_wrapper_->DiscretMap().begin(); cit != dis_wrapper_->DiscretMap().end(); ++cit)
  {
    readers_[cit->first] = Teuchos::rcp(new XSTR::IO::DiscretizationReader(cit->second, step));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationReader::ReadVector(
    Teuchos::RCP<Epetra_Vector> full_vec, std::string name)
{
  // get the vector type of the first discretization
  const enum ::IO::VectorType vt = readers_.begin()->second->VectorType(name);

  XSTR::MultiDiscretizationWrapper::XDisMap::const_iterator cit;
  for (cit = dis_wrapper_->DiscretMap().begin(); cit != dis_wrapper_->DiscretMap().end(); ++cit)
  {
    Teuchos::RCP<Epetra_MultiVector> partial_vec =
        Teuchos::rcp(CreatePartialMultiVector(vt, 1, *cit->second));

    // read the partial vector
    readers_.at(cit->first)->ReadVector(partial_vec, name);

    // insert it into the full vector
    dis_wrapper_->InsertVector(*partial_vec, cit->first, *full_vec);

    // clear the partial vector and go to the next discretization
    partial_vec = Teuchos::null;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Comm& XSTR::IO::MultiDiscretizationReader::Comm() const
{
  return dis_wrapper_->Comm();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationReader::ReadMultiVector(
    Teuchos::RCP<Epetra_MultiVector> full_vec, std::string name)
{
  // get the vector type of the first discretization
  const enum ::IO::VectorType vt = readers_.begin()->second->VectorType(name);
  const int num_vecs = full_vec->NumVectors();

  XSTR::MultiDiscretizationWrapper::XDisMap::const_iterator cit;
  for (cit = dis_wrapper_->DiscretMap().begin(); cit != dis_wrapper_->DiscretMap().end(); ++cit)
  {
    Teuchos::RCP<Epetra_MultiVector> partial_vec =
        Teuchos::rcp(CreatePartialMultiVector(vt, num_vecs, *cit->second));
    // read the partial vector
    readers_.at(cit->first)->ReadMultiVector(partial_vec, name);

    // insert it into the full vector
    dis_wrapper_->InsertVector(*partial_vec, cit->first, *full_vec, vt);

    // clear the partial vector and go to the next discretization
    partial_vec = Teuchos::null;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Epetra_MultiVector* XSTR::IO::CreatePartialMultiVector(
    const enum ::IO::VectorType& vt, const int& num_vecs, const DRT::DiscretizationInterface& dis)
{
  Epetra_MultiVector* partial_vec = NULL;
  switch (vt)
  {
    case ::IO::dofvector:
    {
      if (dis.NumDofSets() > 1) dserror("Currently we support only one dofset!");
      const Epetra_Map* dof_row_map = dis.DofRowMap(0);
      partial_vec = new Epetra_MultiVector(*dof_row_map, num_vecs);
      break;
    }
    case ::IO::nodevector:
    {
      const Epetra_Map* node_row_map = dis.NodeRowMap();
      partial_vec = new Epetra_MultiVector(*node_row_map, num_vecs);
      break;
    }
    case ::IO::elementvector:
    {
      const Epetra_Map* ele_row_map = dis.ElementRowMap();
      partial_vec = new Epetra_MultiVector(*ele_row_map, num_vecs);
      break;
    }
  }
  return partial_vec;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationReader::ReadSerialDenseMatrix(
    Teuchos::RCP<std::map<int, Teuchos::RCP<Epetra_SerialDenseMatrix>>> mapdata, std::string name)
{
  dserror("Currently unsupported!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XSTR::IO::MultiDiscretizationReader::ReadInt(std::string name)
{
  /* this is supposed to be a redundant information, thus we get it only from
   * the standard structure discretization */
  return readers_.at(XFEM::structure)->ReadInt(name);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double XSTR::IO::MultiDiscretizationReader::ReadDouble(std::string name)
{
  /* this is supposed to be a redundant information, thus we get it only from
   * the standard structure discretization */
  return readers_.at(XFEM::structure)->ReadDouble(name);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationReader::ReadMesh(int step)
{
  XReaderMap::iterator it;
  for (it = readers_.begin(); it != readers_.end(); ++it) it->second->ReadMesh(step);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationReader::ReadNodesOnly(int step)
{
  XReaderMap::iterator it;
  for (it = readers_.begin(); it != readers_.end(); ++it) it->second->ReadNodesOnly(step);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationReader::ReadHistoryData(int step)
{
  XReaderMap::iterator it;
  for (it = readers_.begin(); it != readers_.end(); ++it) it->second->ReadNodesOnly(step);

  // re-do the wrapper internal FillComplete routines
  dis_wrapper_->FillComplete(false, false, false, true, true);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationReader::ReadRedundantDoubleVector(
    Teuchos::RCP<std::vector<double>>& doublevec, const std::string name)
{
  /* this is supposed to be a redundant information, thus we get it only from
   * the standard structure discretization */
  return readers_.at(XFEM::structure)->ReadRedundantDoubleVector(doublevec, name);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationReader::ReadRedundantIntVector(
    Teuchos::RCP<std::vector<int>>& intvec, const std::string name)
{
  /* this is supposed to be a redundant information, thus we get it only from
   * the standard structure discretization */
  return readers_.at(XFEM::structure)->ReadRedundantIntVector(intvec, name);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int XSTR::IO::MultiDiscretizationReader::GetNumOutputProc(int step)
{
  /* this is supposed to be a redundant information, thus we get it only from
   * the standard structure discretization */
  return readers_.at(XFEM::structure)->GetNumOutputProc(step);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XSTR::IO::DiscretizationReader::DiscretizationReader(Teuchos::RCP<DRT::DiscretizationInterface> dis,
    Teuchos::RCP<::IO::InputControl> input, int step)
    : ::IO::DiscretizationReader(
          Teuchos::rcp_dynamic_cast<DRT::Discretization>(dis, true), input, step)
{
  // should stay empty, see base class
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XSTR::IO::DiscretizationReader::DiscretizationReader(
    Teuchos::RCP<DRT::DiscretizationInterface> dis, int step)
    : ::IO::DiscretizationReader(Teuchos::rcp_dynamic_cast<DRT::Discretization>(dis, true), step)
{
  // should stay empty, see base class
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
enum ::IO::VectorType XSTR::IO::DiscretizationReader::VectorType(const std::string& name)
{
  enum ::IO::VectorType vector_type = ::IO::dofvector;

  // This function will throw an error, if the discretization contains no vector
  // with the given name!
  MAP* result = map_read_map(RestartStepMap(), name.c_str());

  std::string str_vectortype = map_read_string(result, "type");
  if (str_vectortype == "dof")
    vector_type = ::IO::dofvector;
  else if (str_vectortype == "node")
    vector_type = ::IO::nodevector;
  else if (str_vectortype == "element")
    vector_type = ::IO::elementvector;
  else
    dserror("Unknown vector type: %s", str_vectortype.c_str());

  return vector_type;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XSTR::IO::MultiDiscretizationWriter::MultiDiscretizationWriter(
    const Teuchos::RCP<XSTR::MultiDiscretizationWrapper>& dis_wrapper)
    : ::IO::DiscretizationWriter(), dis_wrapper_(dis_wrapper), curr_step_(-1)
{
  XSTR::MultiDiscretizationWrapper::XDisMap::const_iterator cit;
  for (cit = dis_wrapper_->DiscretMap().begin(); cit != dis_wrapper_->DiscretMap().end(); ++cit)
  {
    writers_[cit->first] = Teuchos::rcp<XSTR::IO::DiscretizationWriter>(
        new XSTR::IO::DiscretizationWriter(cit->second));
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::NewStep(const int step, const double time)
{
  // store the current step counter
  curr_step_ = step;

  // initialize a new step section for each wrapped discretization
  XWriterMap::iterator it;
  for (it = writers_.begin(); it != writers_.end(); ++it) it->second->NewStep(step, time);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::WriteDouble(const std::string name, const double value)
{
  // we write redundant information only into the structure output file
  writers_.at(XFEM::structure)->WriteDouble(name, value);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::WriteInt(const std::string name, const int value)
{
  // we write redundant information only into the structure output file
  writers_.at(XFEM::structure)->WriteInt(name, value);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::WriteVector(const std::string name,
    Teuchos::RCP<const Epetra_MultiVector> full_vec, enum ::IO::VectorType vt)
{
  const int num_vecs = full_vec->NumVectors();
  XSTR::MultiDiscretizationWrapper::XDisMap::const_iterator cit;
  for (cit = dis_wrapper_->DiscretMap().begin(); cit != dis_wrapper_->DiscretMap().end(); ++cit)
  {
    Teuchos::RCP<Epetra_MultiVector> partial_vec =
        Teuchos::rcp(CreatePartialMultiVector(vt, num_vecs, *cit->second));

    // extract partial vector from the full vector
    dis_wrapper_->ExtractVector(*full_vec, cit->first, *partial_vec, vt);

    // read the partial vector
    writers_.at(cit->first)->WriteVector(name, partial_vec, vt);

    // clear the partial vector and go to the next discretization
    partial_vec = Teuchos::null;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::WriteVector(const std::string name,
    const std::vector<char>& vec, const Epetra_Map& elemap, enum ::IO::VectorType vt)
{
  dserror("Currently unsupported!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::CreateNewResultAndMeshFile()
{
  XWriterMap::const_iterator cit;
  for (cit = writers_.begin(); cit != writers_.end(); ++cit)
    cit->second->CreateNewResultAndMeshFile();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::WriteMesh(const int step, const double time)
{
  XWriterMap::const_iterator cit;
  for (cit = writers_.begin(); cit != writers_.end(); ++cit) cit->second->WriteMesh(step, time);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::WriteMesh(
    const int step, const double time, std::string name_base_file)
{
  dserror("Currently unsupported!");
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::WriteElementData(bool writeowner)
{
  XWriterMap::const_iterator cit;
  for (cit = writers_.begin(); cit != writers_.end(); ++cit)
  {
    std::vector<char> buffer(0);
    unsigned bufferlength = cit->second->AugmentControlFile(curr_step_, buffer);
    cit->second->WriteElementData(writeowner);
    cit->second->AugmentControlFile(&buffer[0], bufferlength);
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::WriteNodeData(bool writeowner)
{
  XWriterMap::const_iterator cit;
  for (cit = writers_.begin(); cit != writers_.end(); ++cit) cit->second->WriteNodeData(writeowner);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::WriteRedundantDoubleVector(
    const std::string name, Teuchos::RCP<std::vector<double>> doublevec)
{
  // we write redundant information only into the structure output file
  writers_.at(XFEM::structure)->WriteRedundantDoubleVector(name, doublevec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::WriteRedundantIntVector(
    const std::string name, Teuchos::RCP<std::vector<int>> intvec)
{
  // we write redundant information only into the structure output file
  writers_.at(XFEM::structure)->WriteRedundantIntVector(name, intvec);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::OverwriteResultFile()
{
  XWriterMap::const_iterator cit;
  for (cit = writers_.begin(); cit != writers_.end(); ++cit) cit->second->OverwriteResultFile();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::NewResultFile(int numb_run)
{
  XWriterMap::const_iterator cit;
  for (cit = writers_.begin(); cit != writers_.end(); ++cit) cit->second->NewResultFile(numb_run);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::NewResultFile(std::string name_appendix, int numb_run)
{
  XWriterMap::const_iterator cit;
  for (cit = writers_.begin(); cit != writers_.end(); ++cit)
    cit->second->NewResultFile(name_appendix, numb_run);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::ClearMapCache()
{
  XWriterMap::const_iterator cit;
  for (cit = writers_.begin(); cit != writers_.end(); ++cit) cit->second->ClearMapCache();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::IO::OutputControl> XSTR::IO::MultiDiscretizationWriter::Output() const
{
  /* all writers share the same output control file, thus we take the one of
   * the structure writer */
  return writers_.at(XFEM::structure)->Output();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::MultiDiscretizationWriter::SetOutput(Teuchos::RCP<::IO::OutputControl> output)
{
  // we set the same output control file pointer in all wrapped writers
  XWriterMap::const_iterator cit;
  for (cit = writers_.begin(); cit != writers_.end(); ++cit) cit->second->SetOutput(output);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
const Epetra_Comm& XSTR::IO::MultiDiscretizationWriter::Comm() const
{
  return dis_wrapper_->Comm();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
XSTR::IO::DiscretizationWriter::DiscretizationWriter(Teuchos::RCP<DRT::DiscretizationInterface> dis)
    : ::IO::DiscretizationWriter(Teuchos::rcp_dynamic_cast<DRT::Discretization>(dis, true))
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
unsigned XSTR::IO::DiscretizationWriter::AugmentControlFile(int step, std::vector<char>& ebuffer)
{
  if (Comm().MyPID() != 0) return 0;

  std::fstream& icontrol = output_->ControlFile();
  if (icontrol.is_open()) icontrol.close();
  std::stringstream full_file_name;
  full_file_name << output_->FileName() << ".control";
  icontrol.open(full_file_name.str().c_str(), std::fstream::in);

  // start at the beginning
  icontrol.seekg(0, icontrol.beg);
  unsigned cline = 0;
  std::stringstream search_stream;
  std::string line;

  while (std::getline(icontrol, line))
  {
    ++cline;
    // found result line
    if (line.find("result:", 0) != std::string::npos)
    {
      ++cline;
      // get next line and check it
      std::getline(icontrol, line);
      search_stream.clear();
      search_stream.str("field = \"" + dis_->Name() + "\"");

      if (line.find(search_stream.str(), 0) == std::string::npos) continue;

      cline += 2;
      // if successful, skip the next line
      icontrol.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      // get next line and check it
      std::getline(icontrol, line);

      // clear the search_stream
      search_stream.clear();
      search_stream.str(std::string());

      search_stream << "step = " << step;
      if (line.find(search_stream.str(), 0) != std::string::npos)
      {
        ++cline;
        // go one line further and stop
        icontrol.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        break;
      }
    }
  }
  // check if the search was successful
  if (icontrol.eof()) dserror("We couldn't find the right result section!");

  // -----------------------
  // read data of first block
  // -----------------------
  // get the length and initialize a buffer
  int slength = icontrol.tellg();
  char sbuffer[slength];

  // go back to start first
  icontrol.seekg(0, icontrol.beg);
  icontrol.read(sbuffer, slength);

  // -----------------------
  // read data of second block
  // -----------------------
  icontrol.seekg(0, icontrol.end);
  int elength = static_cast<unsigned>(icontrol.tellg()) - slength;
  ebuffer.resize(elength);

  // go to the end of the first block
  icontrol.seekg(slength, icontrol.beg);
  icontrol.read(&ebuffer[0], elength);

  //  std::cout << " :::: " << dis_->Name() << " :::: \n";
  //  std::cout << "=== sbuffer ===\n";
  //  std::cout << sbuffer.get();
  //
  //  std::cout << "=== ebuffer ===\n";
  //  std::cout << ebuffer.get();

  // close the file and reopen it to add the sbuffer
  icontrol.close();
  icontrol.open(full_file_name.str().c_str(), std::fstream::out | std::fstream::trunc);
  icontrol.write(sbuffer, slength);
  icontrol << std::flush;

  return static_cast<unsigned>(elength);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void XSTR::IO::DiscretizationWriter::AugmentControlFile(const char* buffer, unsigned bufferlength)
{
  if (Comm().MyPID() != 0) return;

  output_->ControlFile().write(buffer, bufferlength);
}
