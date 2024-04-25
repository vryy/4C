/*----------------------------------------------------------------------*/
/*! \file

\brief node based VTU filter

\level 3

*/
/*----------------------------------------------------------------------*/


#include "4C_post_vtk_vtu_writer_node_based.hpp"

#include "4C_beam3_base.hpp"
#include "4C_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_discretization_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_element.hpp"
#include "4C_lib_element_vtk_cell_type_register.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_nurbs_discret.hpp"
#include "4C_post_common.hpp"
#include "4C_post_vtk_vtu_writer.hpp"
#include "4C_utils_exceptions.hpp"

#include <sstream>

FOUR_C_NAMESPACE_OPEN

PostVtuWriterNode::PostVtuWriterNode(PostField* field, const std::string& filename)
    : PostVtuWriter(field, filename)
{
  static_assert(29 == static_cast<int>(CORE::FE::CellType::max_distype),
      "The number of element types defined by CORE::FE::CellType does not match the "
      "number of element types supported by the post vtu filter.");
  if (myrank_ != 0) FOUR_C_THROW("Node based filtering only works in serial mode");
}


const std::string& PostVtuWriterNode::WriterString() const
{
  static std::string name("UnstructuredGrid");
  return name;
}

const std::string& PostVtuWriterNode::WriterOpeningTag() const
{
  static std::string tag("<UnstructuredGrid>");
  return tag;
}

const std::string& PostVtuWriterNode::WriterPOpeningTag() const
{
  static std::string tag("<PUnstructuredGrid GhostLevel=\"0\">");
  return tag;
}

const std::vector<std::string>& PostVtuWriterNode::WriterPPieceTags() const
{
  static std::vector<std::string> tags;
  tags.clear();
  for (size_t i = 0; i < numproc_; ++i)
  {
    std::stringstream stream;
    stream << "<Piece Source=\"" << filenamebase_ << "-" << i << ".vtu\"/>";
    tags.push_back(std::string(stream.str()));
  }
  return tags;
}

const std::string& PostVtuWriterNode::WriterSuffix() const
{
  static std::string name(".vtu");
  return name;
}

const std::string& PostVtuWriterNode::WriterPSuffix() const
{
  static std::string name(".pvtu");
  return name;
}


void PostVtuWriterNode::WriteGeo()
{
  using namespace FourC;

  Teuchos::RCP<DRT::Discretization> dis = this->GetField()->discretization();

  // count number of nodes and number for each processor; output is completely independent of
  // the number of processors involved
  int nelements = dis->NumMyRowElements();
  int nnodes = dis->NumMyRowNodes();
  std::vector<int32_t> connectivity;
  connectivity.reserve(nnodes);
  std::vector<double> coordinates;
  coordinates.reserve(3 * nnodes);
  std::vector<uint8_t> celltypes;
  celltypes.reserve(nelements);
  std::vector<int32_t> celloffset;
  celloffset.reserve(nelements);

  // loop over my elements and write the data
  int outNodeId = 0;
  for (int e = 0; e < nelements; ++e)
  {
    const DRT::Element* ele = dis->lRowElement(e);
    // check for beam element that potentially needs special treatment due to Hermite interpolation
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    if (ele->IsNurbsElement())
    {
      WriteGeoNurbsEle(ele, celltypes, outNodeId, celloffset, coordinates);
    }
    else if (beamele != nullptr)
    {
      WriteGeoBeamEle(beamele, celltypes, outNodeId, celloffset, coordinates);
    }
    else
    {
      celltypes.push_back(
          DRT::ELEMENTS::GetVtkCellTypeFromFourCElementShapeType(ele->Shape()).first);
      const std::vector<int>& numbering =
          DRT::ELEMENTS::GetVtkCellTypeFromFourCElementShapeType(ele->Shape()).second;
      const DRT::Node* const* nodes = ele->Nodes();
      for (int n = 0; n < ele->NumNode(); ++n)
      {
        connectivity.push_back(nodes[numbering[n]]->LID());
      }
      outNodeId += ele->NumNode();
      celloffset.push_back(outNodeId);
    }
  }

  outNodeId = 0;
  for (int n = 0; n < nnodes; n++)
  {
    for (int d = 0; d < 3; ++d)
    {
      coordinates.push_back(dis->lRowNode(n)->X()[d]);
      outNodeId++;
    }
  }

  FOUR_C_ASSERT((int)coordinates.size() == 3 * outNodeId, "internal error");

  // step 1: write node coordinates into file
  currentout_ << "<Piece NumberOfPoints=\"" << nnodes << "\" NumberOfCells=\"" << nelements
              << "\" >\n"
              << "  <Points>\n"
              << "    <DataArray type=\"Float64\" NumberOfComponents=\"3\"";

  if (write_binary_output_)
  {
    currentout_ << " format=\"binary\">\n";
    LIBB64::writeCompressedBlock(coordinates, currentout_);
  }
  else
  {
    currentout_ << " format=\"ascii\">\n";
    int counter = 1;
    for (std::vector<double>::const_iterator it = coordinates.begin(); it != coordinates.end();
         ++it)
    {
      currentout_ << std::setprecision(15) << std::scientific << *it << " ";
      // dimension is hard coded to three, thus
      if (counter % 3 == 0) currentout_ << '\n';
      counter++;
    }
    currentout_ << std::resetiosflags(std::ios::scientific);
  }



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
  if (write_binary_output_)
    currentout_ << " format=\"binary\">\n";
  else
    currentout_ << " format=\"ascii\">\n";


  if (write_binary_output_)
    LIBB64::writeCompressedBlock(connectivity, currentout_);
  else
  {
    for (std::vector<int32_t>::const_iterator it = connectivity.begin(); it != connectivity.end();
         ++it)
      currentout_ << *it << " ";
  }



  currentout_ << "    </DataArray>\n";

  // step 3: write start indices for individual cells
  currentout_ << "    <DataArray type=\"Int32\" Name=\"offsets\"";
  if (write_binary_output_)
  {
    currentout_ << " format=\"binary\">\n";
    LIBB64::writeCompressedBlock(celloffset, currentout_);
  }
  else
  {
    currentout_ << " format=\"ascii\">\n";
    for (std::vector<int32_t>::const_iterator it = celloffset.begin(); it != celloffset.end(); ++it)
      currentout_ << *it << " ";
  }
  currentout_ << "    </DataArray>\n";

  // step 4: write cell types
  currentout_ << "    <DataArray type=\"UInt8\" Name=\"types\"";
  if (write_binary_output_)
  {
    currentout_ << " format=\"binary\">\n";
    LIBB64::writeCompressedBlock(celltypes, currentout_);
  }
  else
  {
    currentout_ << " format=\"ascii\">\n";
    for (std::vector<uint8_t>::const_iterator it = celltypes.begin(); it != celltypes.end(); ++it)
      currentout_ << (unsigned int)*it << " ";
  }

  currentout_ << "    </DataArray>\n";

  currentout_ << "  </Cells>\n\n";

  if (myrank_ == 0)
  {
    currentmasterout_ << "    <PPoints>\n";
    currentmasterout_ << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n";
    currentmasterout_ << "    </PPoints>\n";
  }
}



void PostVtuWriterNode::WriteDofResultStep(std::ofstream& file,
    const Teuchos::RCP<Epetra_Vector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
    const std::string& groupname, const std::string& name, const int numdf, const int from,
    const bool fillzeros)
{
  using namespace FourC;

  if (myrank_ == 0 && timestep_ == 0) std::cout << "writing dof-based field " << name << std::endl;

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  // For parallel computations, we need to access all dofs on the elements, including the
  // nodes owned by other processors. Therefore, we need to import that data here.
  const Epetra_BlockMap& vecmap = data->Map();
  const Epetra_Map* colmap = dis->DofColMap(0);

  int offset = vecmap.MinAllGID() - dis->DofRowMap()->MinAllGID();
  if (fillzeros) offset = 0;

  Teuchos::RCP<Epetra_Vector> ghostedData;
  if (colmap->SameAs(vecmap))
    ghostedData = data;
  else
  {
    // There is one more complication: The map of the vector and the map governed by the
    // degrees of freedom in the discretization might be offset (e.g. pressure for fluid).
    // Therefore, we need to adjust the numbering in the vector to the numbering in the
    // discretization by the variable 'offset'.
    std::vector<int> gids(vecmap.NumMyElements());
    for (int i = 0; i < vecmap.NumMyElements(); ++i)
      gids[i] = vecmap.MyGlobalElements()[i] - offset;
    Teuchos::RCP<Epetra_Map> rowmap = Teuchos::rcp(new Epetra_Map(
        vecmap.NumGlobalElements(), vecmap.NumMyElements(), gids.data(), 0, vecmap.Comm()));
    Teuchos::RCP<Epetra_Vector> dofvec = CORE::LINALG::CreateVector(*rowmap, false);
    for (int i = 0; i < vecmap.NumMyElements(); ++i) (*dofvec)[i] = (*data)[i];

    ghostedData = CORE::LINALG::CreateVector(*colmap, true);
    CORE::LINALG::Export(*dofvec, *ghostedData);
  }

  int ncomponents = numdf;
  if (numdf > 1 && numdf == field_->problem()->num_dim()) ncomponents = 3;

  // count number of nodes for each processor
  int nnodes = dis->NumMyRowNodes();
  std::vector<double> solution;

  solution.reserve(ncomponents * nnodes);

  std::vector<int> nodedofs;

  for (int i = 0; i < nnodes; i++)
  {
    nodedofs.clear();

    // local storage position of desired dof gid
    dis->Dof(dis->lRowNode(i), nodedofs);
    for (int d = 0; d < numdf; ++d)
    {
      const int lid = ghostedData->Map().LID(nodedofs[d + from]);
      if (lid > -1)
        solution.push_back((*ghostedData)[lid]);
      else
      {
        if (fillzeros)
          solution.push_back(0.);
        else
          FOUR_C_THROW("received illegal dof local id: %d", lid);
      }
    }

    for (int d = numdf; d < ncomponents; ++d) solution.push_back(0.);

  }  // loop over all nodes

  FOUR_C_ASSERT((int)solution.size() == ncomponents * nnodes, "internal error");


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
    FOUR_C_THROW(
        "Cannot write point data at this stage. Most likely cell and point data fields are mixed.");

  this->WriteSolutionVector(solution, ncomponents, name, file);
}


void PostVtuWriterNode::WriteNodalResultStep(std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
    const std::string& groupname, const std::string& name, const int numdf)
{
  using namespace FourC;

  if (myrank_ == 0 && timestep_ == 0) std::cout << "writing node-based field " << name << std::endl;

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  // Here is the only thing we need to do for parallel computations: We need read access to all dofs
  // on the row elements, so need to get the NodeColMap to have this access
  const Epetra_Map* colmap = dis->NodeColMap();
  const Epetra_BlockMap& vecmap = data->Map();

  FOUR_C_ASSERT(
      colmap->MaxAllGID() == vecmap.MaxAllGID() && colmap->MinAllGID() == vecmap.MinAllGID(),
      "Given data vector does not seem to match discretization node map");

  Teuchos::RCP<Epetra_MultiVector> ghostedData;
  if (colmap->SameAs(vecmap))
    ghostedData = data;
  else
  {
    ghostedData = Teuchos::rcp(new Epetra_MultiVector(*colmap, data->NumVectors(), false));
    CORE::LINALG::Export(*data, *ghostedData);
  }

  int ncomponents = numdf;
  if (numdf > 1 && numdf == field_->problem()->num_dim()) ncomponents = 3;



  // count number of nodes for each processor
  int nnodes = dis->NumMyRowNodes();

  std::vector<double> solution;
  solution.reserve(ncomponents * nnodes);

  for (int i = 0; i < nnodes; i++)
  {
    for (int idf = 0; idf < numdf; ++idf)
    {
      Epetra_Vector* column = (*ghostedData)(idf);
      solution.push_back((*column)[i]);
    }
    for (int d = numdf; d < ncomponents; ++d) solution.push_back(0.);
  }  // loop over all nodes

  FOUR_C_ASSERT((int)solution.size() == ncomponents * nnodes, "internal error");


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
    FOUR_C_THROW(
        "Cannot write point data at this stage. Most likely cell and point data fields are mixed.");

  this->WriteSolutionVector(solution, ncomponents, name, file);
}



void PostVtuWriterNode::WriteElementResultStep(std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
    const std::string& groupname, const std::string& name, const int numdf, const int from)
{
  if (myrank_ == 0 && timestep_ == 0)
  {
    std::cout << "WARNING: Cannot write element-based quantity in node based vtu filter "
              << std::endl;
    std::cout << "Skipping field " << name.c_str() << std::endl;
  }
}

void PostVtuWriterNode::WriteGeoNurbsEle(const DRT::Element* ele, std::vector<uint8_t>& celltypes,
    int& outNodeId, std::vector<int32_t>& celloffset, std::vector<double>& coordinates)
{
  FOUR_C_THROW("VTU node based filter cannot handle NURBS elements");
}

void PostVtuWriterNode::WriteGeoBeamEle(const DRT::ELEMENTS::Beam3Base* beamele,
    std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
    std::vector<double>& coordinates)
{
  FOUR_C_THROW("VTU node based filter cannot handle beam elements");
}

void PostVtuWriterNode::WirteDofResultStepNurbsEle(const DRT::Element* ele, int ncomponents,
    const int numdf, std::vector<double>& solution, Teuchos::RCP<Epetra_Vector> ghostedData,
    const int from, const bool fillzeros)
{
  FOUR_C_THROW("VTU node based filter cannot handle NURBS elements");
}

void PostVtuWriterNode::WriteDofResultStepBeamEle(const DRT::ELEMENTS::Beam3Base* beamele,
    const int& ncomponents, const int& numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_Vector>& ghostedData, const int& from, const bool fillzeros)
{
  FOUR_C_THROW("VTU node based filter cannot handle beam elements");
}

void PostVtuWriterNode::WriteNodalResultStepNurbsEle(const DRT::Element* ele, int ncomponents,
    const int numdf, std::vector<double>& solution, Teuchos::RCP<Epetra_MultiVector> ghostedData)
{
  FOUR_C_THROW("VTU node based filter cannot handle NURBS elements");
}

FOUR_C_NAMESPACE_CLOSE
