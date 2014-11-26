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


#include "post_drt_vtu_writer.H"

#include <sstream>
#include <boost/static_assert.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../post_drt_common/post_drt_common.H"

// deactivate for ascii output. Only do this for debugging.
#define BIN_VTK_OUT

namespace
{
  template <typename T>
  class make_vector {
  public:
    make_vector<T>& operator<< (const T& val) {
      data_.push_back(val);
      return *this;
    }
    operator std::vector<T>() const {
      return data_;
    }
  private:
    std::vector<T> data_;
  };
}



/*----------------------------------------------------------------------*/
/*
 \brief Converts between our element types and the VTK numbers and the respective numbering of degrees of freedom
 */
/*----------------------------------------------------------------------*/
std::pair<uint8_t,std::vector<int> > vtk_element_types [] =
{
    // the VTK element types are from the documentation of vtkCellType,
    // e.g. at http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
    // this list must be kept in sync with the element types since we use this
    // for index translation
    std::pair<uint8_t,std::vector<int> > (0, std::vector<int>()),                                                    // dis_none
    std::pair<uint8_t,std::vector<int> > (9, make_vector<int>() << 0 << 1 << 2 << 3),                                // quad4
    std::pair<uint8_t,std::vector<int> > (30, make_vector<int>() << 0 << 1 << 4 << 3 << 2 << 5),                     // quad6
    std::pair<uint8_t,std::vector<int> > (23, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7),           // quad8
    std::pair<uint8_t,std::vector<int> > (28, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8),      // quad9
    std::pair<uint8_t,std::vector<int> > (5, make_vector<int>() << 0 << 1 << 2),                                     // tri3
    std::pair<uint8_t,std::vector<int> > (22, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5),                     // tri6
    std::pair<uint8_t,std::vector<int> > (12, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7),           // hex8
    std::pair<uint8_t,std::vector<int> > (12, make_vector<int>() << 0 << 1 << 2 << 3 << 9 << 10 << 11 << 12),        // hex18
    std::pair<uint8_t,std::vector<int> > (25, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8        // hex20
                                          << 9 << 10 << 11 << 16 << 17 << 18 << 19 << 12 << 13 << 14 << 15),
    std::pair<uint8_t,std::vector<int> > (29, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8        // hex27
                                          << 9 << 10 << 11 << 16 << 17 << 18 << 19 << 12 << 13 << 14 << 15
                                          << 24 << 22 << 21 << 23 << 20 << 25 << 26),
    std::pair<uint8_t,std::vector<int> > (10, make_vector<int>() << 0 << 1 << 2 << 3),                               // tet4
    std::pair<uint8_t,std::vector<int> > (24, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9), // tet10
    std::pair<uint8_t,std::vector<int> > (13, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5),                     // wedge6
    std::pair<uint8_t,std::vector<int> > (26, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8        // wedge15
                                          << 12 << 13 << 14 << 9 << 10 << 11),
    std::pair<uint8_t,std::vector<int> > (14, make_vector<int>() << 0 << 1 << 2 << 3 << 4), // pyramid5
    std::pair<uint8_t,std::vector<int> > (3, make_vector<int>() << 0 << 1),                 // line2
    std::pair<uint8_t,std::vector<int> > (21, make_vector<int>() << 0 << 1 << 2),           // line3
    std::pair<uint8_t,std::vector<int> > (35, make_vector<int>() << 0 << 1 << 2 << 3),      // line4
    std::pair<uint8_t,std::vector<int> > (35, make_vector<int>() << 0 << 1 << 2 << 3),      // line5 -> mapped onto line4
    std::pair<uint8_t,std::vector<int> > (35, make_vector<int>() << 0 << 1 << 2 << 3),      // line6 -> mapped onto line4
    std::pair<uint8_t,std::vector<int> > (1, make_vector<int>() << 0),                      // point1
    std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>()),    // nurbs2, not yet implemented
    std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>()),    // nurbs3, not yet implemented
    std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>()),    // nurbs4, not yet implemented
    std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>()),    // nurbs9, not yet implemented
    std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>()),    // nurbs8, not yet implemented
    std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>())     // nurbs27, not yet implemented
};

VtuWriter::VtuWriter(PostField* field, const std::string &filename) :
    VtkWriter(field, filename)
{
  BOOST_STATIC_ASSERT_MSG((ssize_t)(sizeof(vtk_element_types)/sizeof(std::pair<uint8_t,std::vector<int> >) == (ssize_t)DRT::Element::max_distype), "The number of element types defined by DRT::Element::DiscretizationType does not match the number of element types supported by the post vtu filter.") __attribute__((unused));
}


const std::string& VtuWriter::WriterString() const
{
  static std::string name("UnstructuredGrid");
  return name;
}

const std::string& VtuWriter::WriterOpeningTag() const
{
  static std::string tag("<UnstructuredGrid>");
  return tag;
}

const std::string& VtuWriter::WriterPOpeningTag() const
{
  static std::string tag("<PUnstructuredGrid GhostLevel=\"0\">");
  return tag;
}

const std::vector<std::string>& VtuWriter::WriterPPieceTags() const
{
  static std::vector<std::string> tags;
  tags.clear();
  for (size_t i=0; i<numproc_; ++i)
  {
    std::stringstream stream;
    stream << "<Piece Source=\"" << filenamebase_ << "-" << i << ".vtu\"/>";
    tags.push_back(std::string(stream.str()));
  }
  return tags;
}

const std::string& VtuWriter::WriterSuffix() const
{
  static std::string name(".vtu");
  return name;
}

const std::string& VtuWriter::WriterPSuffix() const
{
  static std::string name(".pvtu");
  return name;
}

void
VtuWriter::WriteGeo()
{
  Teuchos::RCP<DRT::Discretization> dis = this->GetField()->discretization();

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
    celltypes.push_back(vtk_element_types[ele->Shape()].first);
    const std::vector<int> &numbering = vtk_element_types[ele->Shape()].second;
    const DRT::Node* const* nodes = ele->Nodes();
    for (int n=0; n<ele->NumNode(); ++n) {
      for (int d=0; d<3; ++d)
        coordinates.push_back(nodes[numbering[n]]->X()[d]);
    }
    outNodeId += ele->NumNode();
    celloffset.push_back(outNodeId);
  }
  dsassert((int)coordinates.size() == 3*nnodes, "internal error");

  // step 1: write node coordinates into file
  currentout_ << "<Piece NumberOfPoints=\"" << nnodes
      <<"\" NumberOfCells=\"" << nelements << "\" >\n"
      << "  <Points>\n"
      << "    <DataArray type=\"Float64\" NumberOfComponents=\"3\"";

#ifdef BIN_VTK_OUT
  currentout_ << " format=\"binary\">\n";
  LIBB64:: writeCompressedBlock(coordinates, currentout_);
#else
  currentout_ << " format=\"ascii\">\n";
  for (std::vector<double>::const_iterator it = coordinates.begin(); it != coordinates.end(); ++it)
    currentout_ << std::setprecision(15) << std::scientific << *it << " ";
  currentout_ << std::resetiosflags(std::ios::scientific);
#endif


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
              << "    <DataArray type=\"Int32\" Name=\"connectivity\"";
#ifdef BIN_VTK_OUT
  currentout_ << " format=\"binary\">\n";
#else
  currentout_ << " format=\"ascii\">\n";
#endif
  {
    std::vector<int32_t> connectivity;
    connectivity.reserve(nnodes);
    for (int i=0; i<nnodes; ++i)
      connectivity.push_back(i);
#ifdef BIN_VTK_OUT
    LIBB64::writeCompressedBlock(connectivity, currentout_);
#else
  for (std::vector<int32_t>::const_iterator it = connectivity.begin(); it != connectivity.end(); ++it)
    currentout_ << *it << " ";
#endif

  }
  currentout_ << "    </DataArray>\n";

  // step 3: write start indices for individual cells
  currentout_ << "    <DataArray type=\"Int32\" Name=\"offsets\"";
#ifdef BIN_VTK_OUT
  currentout_ << " format=\"binary\">\n";
  LIBB64::writeCompressedBlock(celloffset, currentout_);
#else
  currentout_ << " format=\"ascii\">\n";
  for (std::vector<int32_t>::const_iterator it = celloffset.begin(); it != celloffset.end(); ++it)
    currentout_ << *it << " ";
#endif
  currentout_ << "    </DataArray>\n";

  // step 4: write cell types
  currentout_ << "    <DataArray type=\"UInt8\" Name=\"types\"";
#ifdef BIN_VTK_OUT
  currentout_ << " format=\"binary\">\n";
  LIBB64::writeCompressedBlock(celltypes, currentout_);
#else
  currentout_ << " format=\"ascii\">\n";
  for (std::vector<uint8_t>::const_iterator it = celltypes.begin(); it != celltypes.end(); ++it)
    currentout_ << (unsigned int)*it << " ";
#endif
  currentout_ << "    </DataArray>\n";

  currentout_ << "  </Cells>\n\n";

  if (myrank_ == 0) {
    currentmasterout_ << "    <PPoints>\n";
    currentmasterout_ << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n";
    currentmasterout_ << "    </PPoints>\n";
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
    const bool fillzeros)
{
  if (myrank_==0 && timestep_ == 0)
    std::cout << "writing dof-based field " << name <<std::endl;

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  // Here is the only thing we need to do for parallel computations: We need read access to all dofs
  // on the row elements, so need to get the DofColMap to have this access
  const Epetra_Map* colmap = dis->DofColMap(0);
  const Epetra_BlockMap& vecmap = data->Map();

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
    const std::vector<int> &numbering = vtk_element_types[ele->Shape()].second;
    for (int n=0; n<ele->NumNode(); ++n) {
      nodedofs.clear();

      // local storage position of desired dof gid
      dis->Dof(ele->Nodes()[numbering[n]], nodedofs);

      for (int d=0; d<numdf; ++d) {
        const int lid = ghostedData->Map().LID(nodedofs[d+from]);
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

  // start the scalar fields that will later be written
  if (currentPhase_ == INIT) {
    currentout_ << "  <PointData>\n"; // Scalars=\"scalars\">\n";
    if (myrank_ == 0) {
      currentmasterout_ << "    <PPointData>\n"; // Scalars=\"scalars\">\n";
    }
    currentPhase_ = POINTS;
  }

  if (currentPhase_ != POINTS)
    dserror("Cannot write point data at this stage. Most likely cell and point data fields are mixed.");

  this->WriteSolutionVector(solution, ncomponents, name, file);
}



void
VtuWriter::WriteNodalResultStep(
    std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
    const std::string& groupname,
    const std::string& name,
    const int numdf)
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
    const std::vector<int> &numbering = vtk_element_types[ele->Shape()].second;
    for (int n=0; n<ele->NumNode(); ++n) {

      for (int idf=0; idf<numdf; ++idf)
        {
          Epetra_Vector* column = (*ghostedData)(idf);

          int lid = ghostedData->Map().LID(ele->Nodes()[numbering[n]]->Id());

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

  // start the scalar fields that will later be written
  if (currentPhase_ == INIT) {
    currentout_ << "  <PointData>\n"; // Scalars=\"scalars\">\n";
    if (myrank_ == 0) {
      currentmasterout_ << "    <PPointData>\n"; // Scalars=\"scalars\">\n";
    }
    currentPhase_ = POINTS;
  }

  if (currentPhase_ != POINTS)
    dserror("Cannot write point data at this stage. Most likely cell and point data fields are mixed.");

  this->WriteSolutionVector(solution, ncomponents, name, file);
}



void
VtuWriter::WriteElementResultStep(
    std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type> >& resultfilepos,
    const std::string& groupname,
    const std::string& name,
    const int numdf,
    const int from)
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

  // start the scalar fields that will later be written
  if (currentPhase_ == POINTS) {
    // end the scalar fields
    currentout_ << "  </PointData>\n";
    if (myrank_ == 0) {
      currentmasterout_ << "    </PPointData>\n";
    }
  }

  if (currentPhase_ == INIT || currentPhase_ == POINTS) {
    currentout_ << "  <CellData>\n"; // Scalars=\"scalars\">\n";
    if (myrank_ == 0) {
      currentmasterout_ << "    <PCellData>\n"; // Scalars=\"scalars\">\n";
    }
    currentPhase_ = CELLS;
  }

  if (currentPhase_ != CELLS)
    dserror("Cannot write cell data at this stage. Most likely cell and point data fields are mixed.");

  this->WriteSolutionVector(solution, ncomponents, name, file);
}
