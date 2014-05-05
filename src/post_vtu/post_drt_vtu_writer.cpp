/*----------------------------------------------------------------------*/
/*!
\file post_drt_vtu.cpp

\brief VTU filter

<pre>
Maintainer: Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*/
/*----------------------------------------------------------------------*/



#include <string>

#include <zlib.h>
#include <stdint.h>

#include "post_drt_vtu_writer.H"
#include "../post_drt_common/post_drt_common.H"
#include "../pss_full/pss_cpp.h"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../linalg/linalg_utils.H"



/*----------------------------------------------------------------------*/
/*
 \brief Converts between our element types and the VTK numbers
 */
/*----------------------------------------------------------------------*/
uint8_t vtk_element_types [DRT::Element::max_distype] =
{
    // the VTK element types are from the documentation of vtkCellType,
    // e.g. at http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
    // this list must be kept in sync with the element types since we use this
    // for index translation
    0,  // dis_none
    9,  // quad4
    23, // quad8
    28, // quad9
    5,  // tri3
    22, // tri6
    12, // hex8
    25, // hex20
    29, // hex27
    10, // tet4
    24, // tet10
    13, // wedge6
    26, // wedge15
    14, // pyramid5
    3,  // line2
    21, // line3
    21, // line4 -> mapped onto line3
    21, // line5 -> mapped onto line3
    21, // line6 -> mapped onto line3
    1,  // point1
    static_cast<uint8_t>(-1), // nurbs2, not yet implemented
    static_cast<uint8_t>(-1), // nurbs3, not yet implemented
    static_cast<uint8_t>(-1), // nurbs4, not yet implemented
    static_cast<uint8_t>(-1), // nurbs9, not yet implemented
    static_cast<uint8_t>(-1), // nurbs8, not yet implemented
    static_cast<uint8_t>(-1)  // nurbs27, not yet implemented
};



namespace
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



  template <typename T>
  void writeCompressedBlock (const std::vector<T> &data,
                             std::ostream         &out)
  {
    if (!data.empty()) {
      uLongf compressed_data_length = compressBound(data.size() * sizeof(T));
      char *compressed_data = new char[compressed_data_length];
      int err = compress2((Bytef *) compressed_data, &compressed_data_length,
                          (const Bytef *) &data[0], data.size() * sizeof(T),
                          Z_BEST_COMPRESSION);
      if (err != Z_OK)
        dserror("zlib compression failed");

      // now encode the compression header
      const uint32_t compression_header[4] =
          { 1, /* number of blocks */
            (uint32_t) (data.size() * sizeof(T)), /* size of block */
            (uint32_t) (data.size() * sizeof(T)), /* size of last block */
            (uint32_t) compressed_data_length };  /* list of compressed sizes of blocks */

      char *encoded_header = encode_block((char *) &compression_header[0],
                                          4 * sizeof(compression_header[0]));
      out << encoded_header;
      delete[] encoded_header;

      // next do the compressed
      // data encoding in base64
      char *encoded_data = encode_block(compressed_data,
                                        compressed_data_length);
      delete[] compressed_data;

      out << encoded_data;
      out << std::endl;
      delete[] encoded_data;
    }
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


} // end of anonymous namespace



VtuWriter::VtuWriter(PostField* field,
                     const std::string &filename)
:
    PostWriterBase(field, filename),
    numproc_ (field->problem()->comm()->NumProc()),
    time_ (std::numeric_limits<double>::min()),
    timestep_ (0),
    cycle_ (std::numeric_limits<int>::max())
{}



void
VtuWriter::WriteVtuHeader ()
{
  if (!currentout_)
    dserror("Invalid output stream");

  // TODO: might need BigEndian on some systems
  const std::string byteorder = "LittleEndian";

  currentout_ << "<?xml version=\"1.0\" ?> \n";
  currentout_ << "<!-- \n";
  currentout_ << "# vtk DataFile Version 3.0\n";
  currentout_ << "-->\n";
  currentout_ << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
  currentout_ << " compressor=\"vtkZLibDataCompressor\"";
  currentout_ << " byte_order=\"" << byteorder << "\"";
  currentout_ << ">\n";
  currentout_ << "<UnstructuredGrid>\n";

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
    currentmasterout_ << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"";
    currentmasterout_ << " byte_order=\"" << byteorder << "\"";
    currentmasterout_ << ">\n";
    currentmasterout_ << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
  }
}



void
VtuWriter::WriteVtuFooter(const std::string &filenamebase)
{
  if (!currentout_)
    dserror("Invalid output stream");

  currentout_ << "  </PointData>\n";
  currentout_ << "</Piece>\n";

  currentout_ << "</UnstructuredGrid>\n";
  currentout_ << "</VTKFile>\n";

  currentout_ << std::flush;

  // Also start master file on processor 0
  if (myrank_ == 0) {
    if (!currentmasterout_)
      dserror("Invalid output stream");

    currentmasterout_ << "    </PPointData>\n";
    for (int i=0; i<numproc_; ++i)
      currentmasterout_ << "    <Piece Source=\""
                        << filenamebase << "-" << i << ".vtu\"/>\n";

    currentmasterout_ << "  </PUnstructuredGrid>\n";
    currentmasterout_ << "</VTKFile>\n";

    currentmasterout_ << std::flush;
  }
}



void
VtuWriter::WriteGeo()
{
  Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  // count number of nodes and number for each processor; output is completely independent of
  // the number of processors involved
  int nelements = dis->NumMyRowElements();
  int nnodes = 0;
  for (int e=0; e<nelements; ++e)
    nnodes += dis->lRowElement(e)->NumNode();

  // do not need to store connectivity indices here because we create a
  // contiguous array by the order in which we fill the coordinates (otherwise
  // need to adjust filling in the coordinates).
  std::vector<double> coordinates;
  coordinates.reserve(3*nnodes);
  std::vector<uint8_t> celltypes;
  celltypes.reserve(nelements);
  std::vector<int32_t> celloffset;
  celloffset.reserve(nelements);

  // loop over my elements and write the data
  int outNodeId = 0;
  for (int e=0; e<dis->NumMyRowElements(); ++e) {
    const DRT::Element* ele = dis->lRowElement(e);
    celltypes.push_back(vtk_element_types[ele->Shape()]);
    const DRT::Node* const* nodes = ele->Nodes();
    for (int n=0; n<ele->NumNode(); ++n) {
      for (int d=0; d<3; ++d)
        coordinates.push_back(nodes[n]->X()[d]);
    }
    outNodeId += ele->NumNode();
    celloffset.push_back(outNodeId);
  }
  dsassert((int)coordinates.size() == 3*nnodes, "internal error");

  // step 1: write node coordinates into file
  currentout_ << "<Piece NumberOfPoints=\"" << nnodes
      <<"\" NumberOfCells=\"" << nelements << "\" >\n"
      << "  <Points>\n"
      << "    <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"binary\">\n";
  writeCompressedBlock(coordinates, currentout_);
  currentout_ << "    </DataArray>\n"
      << "  </Points>\n\n";

  // avoid too much memory consumption -> clear coordinates vector now that we're done
  {
    std::vector<double> empty;
    empty.swap(coordinates);
  }

  // step 2: write mesh-node topology into file. we assumed contiguous order of coordinates
  // in this format, so fill the vector only here
  currentout_ << "  <Cells>\n"
              << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n";
  {
    std::vector<int32_t> connectivity;
    connectivity.reserve(nnodes);
    for (int i=0; i<nnodes; ++i)
      connectivity.push_back(i);
    writeCompressedBlock(connectivity, currentout_);
  }
  currentout_ << "    </DataArray>\n";

  // step 3: write start indices for individual cells
  currentout_ << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n";
  writeCompressedBlock(celloffset, currentout_);
  currentout_ << "    </DataArray>\n";

  // step 4: write cell types
  currentout_ << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">\n";
  writeCompressedBlock(celltypes, currentout_);
  currentout_ << "    </DataArray>\n";

  currentout_ << "  </Cells>\n\n";

  // start the scalar fields that will later be written
  currentout_ << "  <PointData Scalars=\"scalars\">\n";


  if (myrank_ == 0) {
    currentmasterout_ << "    <PPoints>\n";
    currentmasterout_ << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n";
    currentmasterout_ << "    </PPoints>\n";
    currentmasterout_ << "    <PPointData Scalars=\"scalars\">\n";
  }
}


void
VtuWriter::WriteSpecialField (
      SpecialFieldInterface &special,
      PostResult& result,   ///< result group in the control file
      const ResultType  restype,
      const std::string &groupname,
      const std::vector<std::string> &fieldnames,
      const std::string &outinfo)
{
  // Vtu writes everything into the same file, so create to each output the
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
VtuWriter::WriteResult(const std::string groupname,
                       const std::string name,
                       const ResultType restype,
                       const int numdf,
                       const int from,
                       const bool fillzeros)
{
  PostResult result(field_);
  bool foundit = false;
  while (result.next_result(groupname))
  {
    if (map_has_map(result.group(), groupname.c_str()))
    {
      foundit = true;
      break;
    }
  }
  if (!foundit) return;

  // jump to the correct location in the data vector. Some fields might only
  // be stored once, so need to catch that case as well
  bool once = false;
  for (int i=0; i<timestep_; ++i)
    if ( not result.next_result(groupname) )
    {
      once = true;
      break;
    }
  if (once)
  {
    result = PostResult(field_);
    result.next_result(groupname);
  }

  if ( not (field_->problem()->SpatialApproximation()=="Polynomial" or
            field_->problem()->SpatialApproximation()=="Meshfree" or
            field_->problem()->SpatialApproximation()=="HDG") )
    dserror("Only polynomial or meshfree approximations can be written with the VTU filter");

  // need dummy structure that is required for the generic writer interface
  // but not needed by the vtu writer.
  std::map<std::string, std::vector<std::ofstream::pos_type> > dummy;

  switch (restype)
    {
    case dofbased:
      {
        const Teuchos::RCP<Epetra_Vector> data = result.read_result(groupname);
        WriteDofResultStep(currentout_, data, dummy, groupname, name, numdf, from, fillzeros);
        break;
      }
    case nodebased:
      {
        const Teuchos::RCP<Epetra_MultiVector> data = result.read_multi_result(groupname);
        WriteNodalResultStep(currentout_, data, dummy, groupname, name, numdf);
        break;
      }
    case elementbased:
      {
        const Teuchos::RCP<Epetra_MultiVector> data = result.read_multi_result(groupname);
        WriteElementResultStep(currentout_, data, dummy, groupname, name, numdf, from);
        break;
      }
    default:
      dserror("Result type not yet implemented");
    }
}



void
VtuWriter::WriteDofResultStep(
    std::ofstream& file,
    const Teuchos::RCP<Epetra_Vector> &data,
    std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
    const std::string& groupname,
    const std::string& name,
    const int numdf,
    const int from,
    const bool fillzeros) const
{
  if (myrank_==0 && timestep_ == 0)
    std::cout << "writing dof-based field " << name <<std::endl;

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  // Here is the only thing we need to do for parallel computations: We need read access to all dofs
  // on the row elements, so need to get the DofColMap to have this access
  const Epetra_Map* colmap = dis->DofColMap(0);
  const Epetra_BlockMap& vecmap = data->Map();
  int offset = vecmap.MinAllGID()-dis->DofRowMap(0)->MinAllGID();
  if (fillzeros)
    offset = 0;

  Teuchos::RCP<Epetra_Vector> ghostedData;
  if (colmap->SameAs(vecmap))
    ghostedData = data;
  else {
    ghostedData = LINALG::CreateVector(*colmap,false);
    LINALG::Export(*data,*ghostedData);
  }

  int ncomponents = numdf;
  if (numdf > 1 && numdf == field_->problem()->num_dim())
    ncomponents = 3;

  // count number of nodes and number of elements for each processor
  int nnodes = 0;
  for (int e=0; e<dis->NumMyRowElements(); ++e)
    nnodes += dis->lRowElement(e)->NumNode();

  std::vector<double> solution;
  solution.reserve(ncomponents*nnodes);

  std::vector<int> nodedofs;
  for (int e=0; e<dis->NumMyRowElements(); ++e) {
    const DRT::Element* ele = dis->lRowElement(e);
    for (int n=0; n<ele->NumNode(); ++n) {
      nodedofs.clear();

      // local storage position of desired dof gid
      dis->Dof(ele->Nodes()[n], nodedofs);

      for (int d=0; d<numdf; ++d) {
        const int lid = ghostedData->Map().LID(nodedofs[d+from]+offset);
        if (lid > -1)
          solution.push_back((*ghostedData)[lid]);
        else {
          if(fillzeros)
            solution.push_back(0.);
          else
            dserror("received illegal dof local id: %d", lid);
        }
      }
      for (int d=numdf; d<ncomponents; ++d)
        solution.push_back(0.);
    }
  }
  dsassert((int)solution.size() == ncomponents*nnodes, "internal error");

  WriteSolutionVector(solution, ncomponents, name, file);
}



void
VtuWriter::WriteNodalResultStep(
    std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
    const std::string& groupname,
    const std::string& name,
    const int numdf) const
{
  if (myrank_==0 && timestep_ == 0)
    std::cout << "writing node-based field " << name <<std::endl;

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  // Here is the only thing we need to do for parallel computations: We need read access to all dofs
  // on the row elements, so need to get the NodeColMap to have this access
  const Epetra_Map* colmap = dis->NodeColMap();
  const Epetra_BlockMap& vecmap = data->Map();

  dsassert(colmap->MaxAllGID() == vecmap.MaxAllGID() &&
           colmap->MinAllGID() == vecmap.MinAllGID(),
           "Given data vector does not seem to match discretization node map");

  Teuchos::RCP<Epetra_MultiVector> ghostedData;
  if (colmap->SameAs(vecmap))
    ghostedData = data;
  else {
    ghostedData = Teuchos::rcp(new Epetra_MultiVector(*colmap,data->NumVectors(),
                                                      false));
    LINALG::Export(*data,*ghostedData);
  }

  int ncomponents = numdf;
  if (numdf > 1 && numdf == field_->problem()->num_dim())
    ncomponents = 3;

  // count number of nodes and number of elements for each processor
  int nnodes = 0;
  for (int e=0; e<dis->NumMyRowElements(); ++e)
    nnodes += dis->lRowElement(e)->NumNode();

  std::vector<double> solution;
  solution.reserve(ncomponents*nnodes);

  for (int e=0; e<dis->NumMyRowElements(); ++e) {
    const DRT::Element* ele = dis->lRowElement(e);
    for (int n=0; n<ele->NumNode(); ++n) {

      for (int idf=0; idf<numdf; ++idf)
        {
          Epetra_Vector* column = (*ghostedData)(idf);

          int lid = ghostedData->Map().LID(ele->Nodes()[n]->Id());

          if (lid > -1)
            solution.push_back((*column)[lid]);
          else {
            dserror("received illegal node local id: %d", lid);
          }
        }
      for (int d=numdf; d<ncomponents; ++d)
        solution.push_back(0.);
    }
  }
  dsassert((int)solution.size() == ncomponents*nnodes, "internal error");

  WriteSolutionVector(solution, ncomponents, name, file);
}



void
VtuWriter::WriteElementResultStep(
    std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
    const std::string& groupname,
    const std::string& name,
    const int numdf,
    const int from) const
{
  if (myrank_==0 && timestep_ == 0)
    std::cout << "writing element-based field " << name << std::endl;

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  int ncomponents = numdf;
  if (numdf > 1 && numdf == field_->problem()->num_dim())
    ncomponents = 3;

  // count number of nodes and number of elements for each processor
  int nnodes = 0;
  for (int e=0; e<dis->NumMyRowElements(); ++e)
    nnodes += dis->lRowElement(e)->NumNode();

  std::vector<double> solution;
  solution.reserve(ncomponents*nnodes);

  const int numcol = data->NumVectors();
  if (numdf+from > numcol)
    dserror("violated column range of Epetra_MultiVector: %d",numcol);

  Teuchos::RCP<Epetra_MultiVector> importedData;
  if (dis->ElementRowMap()->SameAs(data->Map()))
    importedData = data;
  else {
    importedData = Teuchos::rcp(new Epetra_MultiVector(*dis->ElementRowMap(),
                                                       data->NumVectors(),
                                                       false));
    LINALG::Export(*data,*importedData);
  }

  for (int e=0; e<dis->NumMyRowElements(); ++e) {
    const DRT::Element* ele = dis->lRowElement(e);
    for (int n=0; n<ele->NumNode(); ++n) {
      for (int d=0; d<numdf; ++d) {
        Epetra_Vector* column = (*importedData)(d+from);
        solution.push_back((*column)[e]);
      }
      for (int d=numdf; d<ncomponents; ++d)
        solution.push_back(0.);
    }
  }
  dsassert((int)solution.size() == ncomponents*nnodes, "internal error");

  WriteSolutionVector(solution, ncomponents, name, file);
}



void
VtuWriter::WriteSolutionVector (const std::vector<double> &solution,
                                const int ncomponents,
                                const std::string &name,
                                std::ofstream &file) const
{
  file << "    <DataArray type=\"Float64\" Name=\"" << name << "\"";
  if (ncomponents > 1)
    file << " NumberOfComponents=\"" << ncomponents << "\"";
  file << " format=\"binary\">\n";

  writeCompressedBlock(solution, file);

  file << "    </DataArray>\n";

  std::ofstream & masterfile = const_cast<std::ofstream&>(currentmasterout_);
  if (myrank_ == 0) {
    masterfile << "      <PDataArray type=\"Float64\" Name=\"" << name << "\"";
    if (ncomponents > 1)
      masterfile << " NumberOfComponents=\"" << ncomponents << "\"";
    masterfile << " format=\"ascii\"/>\n";
  }
}



void
VtuWriter::WriteFiles(PostFilterBase &filter)
{
  PostResult result = PostResult(field_);

  // timesteps when the solution is written
  const std::vector<double> soltime = result.get_result_times(field_->name());
  unsigned int ntdigits = ndigits(soltime.size()),
      npdigits = ndigits(field_->discretization()->Comm().NumProc());
  std::vector<std::pair<double, std::string> > filenames;
  for (timestep_=0; timestep_<(int)soltime.size(); ++timestep_) {
    const std::string filename_base = filename_ + "-" + field_->name() + "-" + int2string(timestep_,ntdigits);
    time_ = soltime[timestep_];
    filenames.push_back(std::pair<double,std::string>(time_, filename_base+".pvtu"));

    currentout_.close();
    currentout_.open((filename_base+"-" + int2string(myrank_,npdigits) + ".vtu").c_str());

    if (myrank_ == 0) {
      currentmasterout_.close();
      currentmasterout_.open((filename_base+".pvtu").c_str());
    }

    WriteVtuHeader();

    WriteGeo();

    filter.WriteAllResults(field_);

    WriteVtuFooter(filename_base);
  }


  // finally, write a single masterfile
  if (myrank_ == 0) {
    std::ofstream masterfile((filename_ + "-" + field_->name() + ".pvd").c_str());

    masterfile << "<?xml version=\"1.0\"?>\n";

    masterfile << "<!--\n";
    masterfile << "-->\n";

    masterfile
        << "<VTKFile type=\"Collection\" version=\"0.1\" ByteOrder=\"LittleEndian\">\n";
    masterfile << "  <Collection>\n";

    for (unsigned int i = 0; i < filenames.size(); ++i)
      masterfile << "    <DataSet timestep=\"" << filenames[i].first
                 << "\" group=\"\" part=\"0\" file=\"" << filenames[i].second
                 << "\"/>\n";

    masterfile << "  </Collection>\n";
    masterfile << "</VTKFile>\n";

    masterfile.flush();
  }

}
