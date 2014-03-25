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
      dsassert(err == Z_OK, "zlib compression failed");

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

} // end of anonymous namespace



VtuWriter::VtuWriter(PostField* field,
                     const std::string filename)
  : field_(field),
    filename_(filename),
    myrank_(((field->problem())->comm())->MyPID()),
    time_ (std::numeric_limits<double>::min()),
    timestep_ (0),
    cycle_ (std::numeric_limits<int>::max())
{}



void
VtuWriter::WriteVtuHeader ()
{
  if (myrank_ != 0)
    return;

  if (!currentout_)
    dserror("Invalid output stream");
  currentout_ << "<?xml version=\"1.0\" ?> \n";
  currentout_ << "<!-- \n";
  currentout_ << "# vtk DataFile Version 3.0\n";
  currentout_ << "-->\n";
  currentout_ << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
  currentout_ << " compressor=\"vtkZLibDataCompressor\"";
  currentout_ << " byte_order=\"LittleEndian\""; // TODO: might need BigEndian on some systems
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
}



void
VtuWriter::WriteVtuFooter()
{
  if (myrank_ != 0)
    return;

  if (!currentout_)
    dserror("Invalid output stream");

  currentout_ << "  </PointData>\n";
  currentout_ << "</Piece>\n";

  currentout_ << "</UnstructuredGrid>\n";
  currentout_ << "</VTKFile>\n";

  currentout_ << std::flush;
}



void
VtuWriter::WriteGeo()
{
  Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  // count number of nodes and number of elements globally
  int nelements = dis->NumGlobalElements();
  int myNnodes = 0, nnodes = 0;
  for (int e=0; e<dis->NumMyRowElements(); ++e)
    myNnodes += dis->lRowElement(e)->NumNode();
  dis->Comm().SumAll(&myNnodes, &nnodes, 1);

  // do not need to store connectivity indices here because we create a
  // contiguous array by the order in which we fill the coordinates (otherwise
  // need to adjust filling in the coordinates).
  std::vector<double> coordinates;
  coordinates.reserve(3*myNnodes);
  std::vector<int> mycelltypes;
  mycelltypes.reserve(dis->NumMyRowElements());
  std::vector<int> mycelloffset;
  mycelloffset.reserve(dis->NumMyRowElements());

  // loop over my elements and write the data
  int outNodeId = 0;
  for (int e=0; e<dis->NumMyRowElements(); ++e) {
    const DRT::Element* ele = dis->lRowElement(e);
    mycelltypes.push_back(vtk_element_types[ele->Shape()]);
    const DRT::Node* const* nodes = ele->Nodes();
    for (int n=0; n<ele->NumNode(); ++n) {
      for (int d=0; d<3; ++d)
        coordinates.push_back(nodes[n]->X()[d]);
    }
    outNodeId += ele->NumNode();
    mycelloffset.push_back(outNodeId);
  }
  dsassert((int)coordinates.size() == 3*myNnodes, "internal error");

  // in parallel, need to send data from other processors to proc 0 and possibly
  // shift their node ids

  std::vector<int32_t> celloffset;
  std::vector<uint8_t> celltypes;
  std::vector<float> owner;
  if (myrank_ == 0) {
    coordinates.reserve(nnodes*3);

    // set owner for data of first processor
    owner.reserve(nnodes*3);
    owner.resize(myNnodes, 0.);

    celltypes.reserve(nelements);
    for (unsigned int i=0; i<mycelltypes.size(); ++i)
      celltypes.push_back(static_cast<uint8_t>(mycelltypes[i]));

    celloffset.reserve(nelements);
    celloffset.insert(celloffset.begin(), mycelloffset.begin(), mycelloffset.end());
  }

#ifdef PARALLEL
  DRT::Exporter exporter(dis->Comm());
  if (myrank_ > 0) {
    MPI_Request request[3];
    exporter.ISend(myrank_, 0, &coordinates[0],  coordinates.size(),  3*myrank_+0, request[0]);
    exporter.ISend(myrank_, 0, &mycelltypes[0],  mycelltypes.size(),  3*myrank_+1, request[1]);
    exporter.ISend(myrank_, 0, &mycelloffset[0], mycelloffset.size(), 3*myrank_+2, request[2]);
    for (int i=0; i<3; ++i)
      exporter.Wait(request[i]);
  }
  else {
    std::vector<double> othercoordinates;
    std::vector<int> otherintegers;
    for (int i=1; i<dis->Comm().NumProc(); ++i) {
      int length = 0;
      exporter.Receive(i, 3*i+0, othercoordinates, length);
      coordinates.insert(coordinates.end(), othercoordinates.begin(), othercoordinates.end());

      // add ownership for ith processor
      owner.resize(owner.size()+(othercoordinates.size()/3), i);

      exporter.Receive(i, 3*i+1, otherintegers, length);
      for (int k=0; k<length; ++k)
        celltypes.push_back(static_cast<uint8_t>(otherintegers[k]));
      exporter.Receive(i, 3*i+2, otherintegers, length);
      // need to shift indices by the number of entries from previous processors
      const int shift = celloffset.empty() ? 0 : celloffset.back();
      for (int k=0; k<length; ++k)
        celloffset.push_back(otherintegers[k]+shift);
    }
  }
#endif

  // now do the actual write on processor 0
  if (myrank_ != 0)
    return;

  // check dimensions
  dsassert((int)coordinates.size() == nnodes*3, "Incorrect communication");
  dsassert((int)celltypes.size() == nelements, "Incorrect communication");
  dsassert((int)celloffset.size() == nelements, "Incorrect communication");


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

  // step 5: start data fields by writing out the owner of the cells
  currentout_ << "  <PointData Scalars=\"scalars\">\n";
  currentout_ << "    <DataArray type=\"Float32\" Name=\"owner\" format=\"binary\">\n";
  writeCompressedBlock(owner, currentout_);
  currentout_ << "    </DataArray>\n";
}


// what do I want to re-use from other writers (post_ensight, post_generic)
// generic fluid, structure, acou, etc writer
// maybe even eigenstress writer
// the 'driver' routine in the main file that sets up which fields to write
// how to do: Introduce a base interface in post_drt_common that holds the actual writer object and distinguishes steps
// then finally create the actual writer interface

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

  // jump to the correct location in the data vector
  for (int i=0; i<timestep_; ++i)
    if ( not result.next_result(groupname) )
      dserror("Could not locate time step %d for field %s in result file", timestep_, groupname.c_str());

  if (restype != dofbased)
    dserror("Given result type %i currently not implemented", restype);

  const Teuchos::RCP<Epetra_Vector> data = result.read_result(groupname);

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

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

  // count number of nodes and number of elements globally
  const int myrank = dis->Comm().MyPID();
  int myNnodes = 0, nnodes = 0;
  for (int e=0; e<dis->NumMyRowElements(); ++e)
    myNnodes += dis->lRowElement(e)->NumNode();
  dis->Comm().SumAll(&myNnodes, &nnodes, 1);

  // get the data on the respective owner in complete analogy to the geometry
  std::vector<double> solution;
  solution.reserve(ncomponents*myNnodes);

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
  dsassert((int)solution.size() == ncomponents*myNnodes, "internal error");

#ifdef PARALLEL
  DRT::Exporter exporter(dis->Comm());
  if (myrank > 0) {
    MPI_Request request;
    exporter.ISend(myrank, 0, &solution[0],  solution.size(),  myrank, request);
    exporter.Wait(request);
  }
  else {
    solution.reserve(numdf*nnodes);
    std::vector<double> othersolution;
    for (int i=1; i<dis->Comm().NumProc(); ++i) {
      int length = 0;
      exporter.Receive(i, i, othersolution, length);
      solution.insert(solution.end(), othersolution.begin(), othersolution.end());
    }
  }
#endif

  if (myrank == 0) {
    dsassert((int)solution.size() == ncomponents*nnodes, "Incorrect communication");
    currentout_ << "    <DataArray type=\"Float64\" Name=\"" << name << "\"";
    if (numdf > 1)
      currentout_ << " NumberOfComponents=\"3\"";
    currentout_ << " format=\"binary\">\n";

    writeCompressedBlock(solution, currentout_);

    currentout_<< "    </DataArray>\n";
  }
}



namespace
{
  /**
   \brief Helper function to determine output file string from time step number
   */
  std::string int2string(const unsigned int i,
                         const unsigned int digits)
  {
    dsassert(i<std::pow(10,digits), "Invalid digits information");
    std::string string;
    switch (digits)
    {
    // note: no break at the end of case, so fall through
    case 10:
      string += '0' + i/1000000000;
    case 9:
      string += '0' + (i%1000000000)/100000000;
    case 8:
      string += '0' + (i%100000000)/10000000;
    case 7:
      string += '0' + (i%10000000)/1000000;
    case 6:
      string += '0' + (i%1000000)/100000;
    case 5:
      string += '0' + (i%100000)/10000;
    case 4:
      string += '0' + (i%10000)/1000;
    case 3:
      string += '0' + (i%1000)/100;
    case 2:
      string += '0' + (i%100)/10;
    case 1:
      string += '0' + i%10;
      break;
    default:
      string += "invalid_digit";
    }
    return string;
  }
}


void
VtuWriter::WriteFiles()
{
  PostResult result = PostResult(field_);

  // timesteps when the solution is written
  const std::vector<double> soltime = result.get_result_times(field_->name());
  unsigned int ndigits = 0, solsize = soltime.size();
  while (solsize > 0) {
    ndigits++;
    solsize /= 10;
  }
  std::vector<std::pair<double, std::string> > filenames;
  for (timestep_=0; timestep_<(int)soltime.size(); ++timestep_) {
    const std::string filename = filename_ + "-" + field_->name() + "-" + int2string(timestep_,ndigits) + ".vtu";
    time_ = soltime[timestep_];
    filenames.push_back(std::pair<double,std::string>(time_, filename));
    currentout_.close();
    currentout_.open(filename.c_str());

    WriteVtuHeader();

    WriteGeo();

    WriteAllResults(field_);

    WriteVtuFooter();
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

