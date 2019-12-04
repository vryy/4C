/*----------------------------------------------------------------------*/
/*! \file

\brief VTU filter

\level 2

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/


#include "post_drt_vtu_writer.H"

#include <sstream>

#include "../drt_lib/drt_element_vtk_cell_type_register.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../post_drt_common/post_drt_common.H"

#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_lib/drt_element.H"

#include "../drt_beam3/beam3_base.H"


PostVtuWriter::PostVtuWriter(PostField* field, const std::string& filename)
    : PostVtkWriter(field, filename)
{
}

const std::string& PostVtuWriter::WriterString() const
{
  static std::string name("UnstructuredGrid");
  return name;
}

const std::string& PostVtuWriter::WriterOpeningTag() const
{
  static std::string tag("<UnstructuredGrid>");
  return tag;
}

const std::string& PostVtuWriter::WriterPOpeningTag() const
{
  static std::string tag("<PUnstructuredGrid GhostLevel=\"0\">");
  return tag;
}

const std::vector<std::string>& PostVtuWriter::WriterPPieceTags() const
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

const std::string& PostVtuWriter::WriterSuffix() const
{
  static std::string name(".vtu");
  return name;
}

const std::string& PostVtuWriter::WriterPSuffix() const
{
  static std::string name(".pvtu");
  return name;
}

void PostVtuWriter::WriteGeo()
{
  Teuchos::RCP<DRT::Discretization> dis = this->GetField()->discretization();

  // count number of nodes and number for each processor; output is completely independent of
  // the number of processors involved
  int nelements = dis->NumMyRowElements();
  int nnodes = 0;
  for (int e = 0; e < nelements; ++e) nnodes += dis->lRowElement(e)->NumNode();

  // do not need to store connectivity indices here because we create a
  // contiguous array by the order in which we fill the coordinates (otherwise
  // need to adjust filling in the coordinates).
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

    if (DRT::UTILS::IsNurbsDisType(ele->Shape()))
    {
      WriteGeoNurbsEle(ele, celltypes, outNodeId, celloffset, coordinates);
    }
    else if (beamele != NULL)
    {
      WriteGeoBeamEle(beamele, celltypes, outNodeId, celloffset, coordinates);
    }
    else
    {
      celltypes.push_back(
          DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(ele->Shape()).first);
      const std::vector<int>& numbering =
          DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(ele->Shape()).second;
      const DRT::Node* const* nodes = ele->Nodes();
      for (int n = 0; n < ele->NumNode(); ++n)
        for (int d = 0; d < 3; ++d) coordinates.push_back(nodes[numbering[n]]->X()[d]);
      outNodeId += ele->NumNode();
      celloffset.push_back(outNodeId);
    }
  }
  dsassert((int)coordinates.size() == 3 * outNodeId, "internal error");

  // step 1: write node coordinates into file
  currentout_ << "<Piece NumberOfPoints=\"" << outNodeId << "\" NumberOfCells=\"" << nelements
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

  {
    std::vector<int32_t> connectivity;
    connectivity.reserve(outNodeId);
    for (int i = 0; i < outNodeId; ++i) connectivity.push_back(i);
    if (write_binary_output_)
      LIBB64::writeCompressedBlock(connectivity, currentout_);
    else
    {
      for (std::vector<int32_t>::const_iterator it = connectivity.begin(); it != connectivity.end();
           ++it)
        currentout_ << *it << " ";
    }
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



void PostVtuWriter::WriteDofResultStep(std::ofstream& file, const Teuchos::RCP<Epetra_Vector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
    const std::string& groupname, const std::string& name, const int numdf, const int from,
    const bool fillzeros)
{
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
        vecmap.NumGlobalElements(), vecmap.NumMyElements(), &gids[0], 0, vecmap.Comm()));
    Teuchos::RCP<Epetra_Vector> dofvec = LINALG::CreateVector(*rowmap, false);
    for (int i = 0; i < vecmap.NumMyElements(); ++i) (*dofvec)[i] = (*data)[i];

    ghostedData = LINALG::CreateVector(*colmap, true);
    LINALG::Export(*dofvec, *ghostedData);
  }

  int ncomponents = numdf;
  if (numdf > 1 && numdf == field_->problem()->num_dim()) ncomponents = 3;

  // count number of nodes and number of elements for each processor
  int nnodes = 0;
  for (int e = 0; e < dis->NumMyRowElements(); ++e) nnodes += dis->lRowElement(e)->NumNode();

  std::vector<double> solution;
  solution.reserve(ncomponents * nnodes);

  std::vector<int> nodedofs;
  for (int e = 0; e < dis->NumMyRowElements(); ++e)
  {
    const DRT::Element* ele = dis->lRowElement(e);
    // check for beam element that potentially needs special treatment due to Hermite interpolation
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    if (DRT::UTILS::IsNurbsDisType(ele->Shape()))
    {
      WirteDofResultStepNurbsEle(ele, ncomponents, numdf, solution, ghostedData, from, fillzeros);
    }
    else if (beamele != NULL)
    {
      WriteDofResultStepBeamEle(
          beamele, ncomponents, numdf, solution, ghostedData, from, fillzeros);
    }
    else
    {
      const std::vector<int>& numbering =
          DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(ele->Shape()).second;

      for (int n = 0; n < ele->NumNode(); ++n)
      {
        nodedofs.clear();

        // local storage position of desired dof gid
        dis->Dof(ele->Nodes()[numbering[n]], nodedofs);

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
              dserror("received illegal dof local id: %d", lid);
          }
        }
        for (int d = numdf; d < ncomponents; ++d) solution.push_back(0.);
      }
    }

  }  // loop over all elements
     //  dsassert((int)solution.size() == ncomponents*nnodes, "internal error");

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



void PostVtuWriter::WriteNodalResultStep(std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
    const std::string& groupname, const std::string& name, const int numdf)
{
  if (myrank_ == 0 && timestep_ == 0) std::cout << "writing node-based field " << name << std::endl;

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  // Here is the only thing we need to do for parallel computations: We need read access to all dofs
  // on the row elements, so need to get the NodeColMap to have this access
  const Epetra_Map* colmap = dis->NodeColMap();
  const Epetra_BlockMap& vecmap = data->Map();

  dsassert(colmap->MaxAllGID() == vecmap.MaxAllGID() && colmap->MinAllGID() == vecmap.MinAllGID(),
      "Given data vector does not seem to match discretization node map");

  Teuchos::RCP<Epetra_MultiVector> ghostedData;
  if (colmap->SameAs(vecmap))
    ghostedData = data;
  else
  {
    ghostedData = Teuchos::rcp(new Epetra_MultiVector(*colmap, data->NumVectors(), false));
    LINALG::Export(*data, *ghostedData);
  }

  int ncomponents = numdf;
  if (numdf > 1 && numdf == field_->problem()->num_dim()) ncomponents = 3;


  // count number of nodes and number of elements for each processor
  int nnodes = 0;
  for (int e = 0; e < dis->NumMyRowElements(); ++e) nnodes += dis->lRowElement(e)->NumNode();

  std::vector<double> solution;
  solution.reserve(ncomponents * nnodes);

  for (int e = 0; e < dis->NumMyRowElements(); ++e)
  {
    const DRT::Element* ele = dis->lRowElement(e);

    if (DRT::UTILS::IsNurbsDisType(ele->Shape()))
    {
      WriteNodalResultStepNurbsEle(ele, ncomponents, numdf, solution, ghostedData);
    }
    else
    {
      const std::vector<int>& numbering =
          DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(ele->Shape()).second;
      for (int n = 0; n < ele->NumNode(); ++n)
      {
        for (int idf = 0; idf < numdf; ++idf)
        {
          Epetra_Vector* column = (*ghostedData)(idf);

          int lid = ghostedData->Map().LID(ele->Nodes()[numbering[n]]->Id());

          if (lid > -1)
            solution.push_back((*column)[lid]);
          else
          {
            dserror("received illegal node local id: %d", lid);
          }
        }
        for (int d = numdf; d < ncomponents; ++d) solution.push_back(0.);
      }
    }
  }  // loop over all elements
     //  dsassert((int)solution.size() == ncomponents*nnodes, "internal error");


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



void PostVtuWriter::WriteElementResultStep(std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
    const std::string& groupname, const std::string& name, const int numdf, const int from)
{
  if (myrank_ == 0 && timestep_ == 0)
    std::cout << "writing element-based field " << name << std::endl;

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  int ncomponents = numdf;
  if (numdf > 1 && numdf == field_->problem()->num_dim()) ncomponents = 3;

  // count number of nodes and number of elements for each processor
  int nnodes = 0;
  for (int e = 0; e < dis->NumMyRowElements(); ++e) nnodes += dis->lRowElement(e)->NumNode();

  std::vector<double> solution;
  solution.reserve(ncomponents * nnodes);

  const int numcol = data->NumVectors();
  if (numdf + from > numcol) dserror("violated column range of Epetra_MultiVector: %d", numcol);

  Teuchos::RCP<Epetra_MultiVector> importedData;
  if (dis->ElementRowMap()->SameAs(data->Map()))
    importedData = data;
  else
  {
    importedData =
        Teuchos::rcp(new Epetra_MultiVector(*dis->ElementRowMap(), data->NumVectors(), false));
    LINALG::Export(*data, *importedData);
  }

  for (int e = 0; e < dis->NumMyRowElements(); ++e)
  {
    const DRT::Element* ele = dis->lRowElement(e);
    for (int n = 0; n < ele->NumNode(); ++n)
    {
      for (int d = 0; d < numdf; ++d)
      {
        Epetra_Vector* column = (*importedData)(d + from);
        solution.push_back((*column)[e]);
      }
      for (int d = numdf; d < ncomponents; ++d) solution.push_back(0.);
    }
  }
  dsassert((int)solution.size() == ncomponents * nnodes, "internal error");

  // start the scalar fields that will later be written
  if (currentPhase_ == POINTS)
  {
    // end the scalar fields
    currentout_ << "  </PointData>\n";
    if (myrank_ == 0)
    {
      currentmasterout_ << "    </PPointData>\n";
    }
  }

  if (currentPhase_ == INIT || currentPhase_ == POINTS)
  {
    currentout_ << "  <CellData>\n";  // Scalars=\"scalars\">\n";
    if (myrank_ == 0)
    {
      currentmasterout_ << "    <PCellData>\n";  // Scalars=\"scalars\">\n";
    }
    currentPhase_ = CELLS;
  }

  if (currentPhase_ != CELLS)
    dserror(
        "Cannot write cell data at this stage. Most likely cell and point data fields are mixed.");

  this->WriteSolutionVector(solution, ncomponents, name, file);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void PostVtuWriter::WriteGeoNurbsEle(const DRT::Element* ele, std::vector<uint8_t>& celltypes,
    int& outNodeId, std::vector<int32_t>& celloffset, std::vector<double>& coordinates) const
{
  switch (ele->Shape())
  {
    case DRT::Element::nurbs2:
    {
      WriteGeoNurbsEle<DRT::Element::nurbs2>(ele, celltypes, outNodeId, celloffset, coordinates);
      break;
    }
    case DRT::Element::nurbs3:
    {
      WriteGeoNurbsEle<DRT::Element::nurbs3>(ele, celltypes, outNodeId, celloffset, coordinates);
      break;
    }
    case DRT::Element::nurbs4:
    {
      WriteGeoNurbsEle<DRT::Element::nurbs4>(ele, celltypes, outNodeId, celloffset, coordinates);
      break;
    }
    case DRT::Element::nurbs9:
    {
      WriteGeoNurbsEle<DRT::Element::nurbs9>(ele, celltypes, outNodeId, celloffset, coordinates);
      break;
    }
    case DRT::Element::nurbs8:
    {
      WriteGeoNurbsEle<DRT::Element::nurbs8>(ele, celltypes, outNodeId, celloffset, coordinates);
      break;
    }
    case DRT::Element::nurbs27:
    {
      WriteGeoNurbsEle<DRT::Element::nurbs27>(ele, celltypes, outNodeId, celloffset, coordinates);
      break;
    }
    default:
      dserror("VTK output not yet implemented for given NURBS element");
      exit(EXIT_FAILURE);
  }

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType nurbs_type>
void PostVtuWriter::WriteGeoNurbsEle(const DRT::Element* ele, std::vector<uint8_t>& celltypes,
    int& outNodeId, std::vector<int32_t>& celloffset, std::vector<double>& coordinates) const
{
  const unsigned NUMNODES = DRT::UTILS::DisTypeToNumNodePerEle<nurbs_type>::numNodePerElement;
  const unsigned DIM = DRT::UTILS::DisTypeToDim<nurbs_type>::dim;

  const DRT::Element::DiscretizationType mapped_dis_type =
      MapNurbsDisTypeToLagrangeDisType(nurbs_type);

  const std::vector<int>& numbering =
      DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(mapped_dis_type).second;

  Teuchos::RCP<const DRT::Discretization> dis = this->GetField()->discretization();

  celltypes.push_back(DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(mapped_dis_type).first);

  LINALG::Matrix<NUMNODES, 1> weights;
  const DRT::Node* const* nodes = ele->Nodes();
  for (unsigned inode = 0; inode < NUMNODES; inode++)
  {
    const DRT::NURBS::ControlPoint* cp =
        dynamic_cast<const DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }
  const DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(dis.get());

  if (not nurbsdis) dserror("Cast to NURBS discretization failed.\n");

  std::vector<Epetra_SerialDenseVector> myknots(3);

  const bool zero_ele = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, ele->Id());

  if (zero_ele) return;

  for (int n = 0; n < ele->NumNode(); ++n)
  {
    LINALG::Matrix<NUMNODES, 1> funct;
    LINALG::Matrix<DIM, NUMNODES> deriv;  // dummy
    LINALG::Matrix<3, 1> gpa;

    gpa = DRT::UTILS::getNodeCoordinates(numbering[n], mapped_dis_type);

    DRT::NURBS::UTILS::nurbs_get_funct_deriv(funct, deriv, gpa, myknots, weights, nurbs_type);

    double X[3] = {0., 0., 0.};
    for (unsigned i = 0; i < 3; ++i)
      for (unsigned m = 0; m < NUMNODES; ++m) X[i] += funct(m) * (nodes[m]->X()[i]);

    for (unsigned d = 0; d < 3; ++d) coordinates.push_back(X[d]);
  }
  outNodeId += ele->NumNode();
  celloffset.push_back(outNodeId);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
DRT::Element::DiscretizationType PostVtuWriter::MapNurbsDisTypeToLagrangeDisType(
    const DRT::Element::DiscretizationType nurbs_dis_type) const
{
  switch (nurbs_dis_type)
  {
    case DRT::Element::nurbs2:
      return DRT::Element::line2;
    case DRT::Element::nurbs3:
      return DRT::Element::line3;
    case DRT::Element::nurbs4:
      return DRT::Element::quad4;
    case DRT::Element::nurbs9:
      return DRT::Element::quad9;
    case DRT::Element::nurbs8:
      return DRT::Element::hex8;
    case DRT::Element::nurbs27:
      return DRT::Element::hex27;
    default:
      dserror("No known mapping from NURBS to Lagrange.");
      exit(EXIT_FAILURE);
  }
}

void PostVtuWriter::WriteGeoBeamEle(const DRT::ELEMENTS::Beam3Base* beamele,
    std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
    std::vector<double>& coordinates)
{
  /* visualize the beam centerline as a sequence of straight line segments (POLY_LINE)
   * which is supported as vtkCellType number 4 (see also list in GetVtkElementType) */
  celltypes.push_back(4);

  /* loop over the chosen visualization points (equidistant distribution in the element
   * parameter space xi \in [-1,1] ) and determine their interpolated initial positions r */
  LINALG::Matrix<3, 1> r;
  double xi = 0.0;

  for (unsigned int i = 0; i < BEAMSVTUVISUALSUBSEGMENTS + 1; ++i)
  {
    r.Clear();
    xi = -1.0 + i * 2.0 / BEAMSVTUVISUALSUBSEGMENTS;

    beamele->GetRefPosAtXi(r, xi);

    for (int d = 0; d < 3; ++d) coordinates.push_back(r(d));
  }

  outNodeId += BEAMSVTUVISUALSUBSEGMENTS + 1;
  celloffset.push_back(outNodeId);
}

void PostVtuWriter::WirteDofResultStepNurbsEle(const DRT::Element* ele, int ncomponents,
    const int numdf, std::vector<double>& solution, Teuchos::RCP<Epetra_Vector> ghostedData,
    const int from, const bool fillzeros) const
{
  switch (ele->Shape())
  {
    case DRT::Element::nurbs2:
    {
      WirteDofResultStepNurbsEle<DRT::Element::nurbs2>(
          ele, ncomponents, numdf, solution, ghostedData, from, fillzeros);
      break;
    }
    case DRT::Element::nurbs3:
    {
      WirteDofResultStepNurbsEle<DRT::Element::nurbs3>(
          ele, ncomponents, numdf, solution, ghostedData, from, fillzeros);
      break;
    }
    case DRT::Element::nurbs4:
    {
      WirteDofResultStepNurbsEle<DRT::Element::nurbs4>(
          ele, ncomponents, numdf, solution, ghostedData, from, fillzeros);
      break;
    }
    case DRT::Element::nurbs9:
    {
      WirteDofResultStepNurbsEle<DRT::Element::nurbs9>(
          ele, ncomponents, numdf, solution, ghostedData, from, fillzeros);
      break;
    }
    case DRT::Element::nurbs8:
    {
      WirteDofResultStepNurbsEle<DRT::Element::nurbs8>(
          ele, ncomponents, numdf, solution, ghostedData, from, fillzeros);
      break;
    }
    case DRT::Element::nurbs27:
    {
      WirteDofResultStepNurbsEle<DRT::Element::nurbs27>(
          ele, ncomponents, numdf, solution, ghostedData, from, fillzeros);
      break;
    }
    default:
      dserror("VTK output not yet implemented for given NURBS element");
      break;
  }  // end switch shape

  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType nurbs_type>
void PostVtuWriter::WirteDofResultStepNurbsEle(const DRT::Element* ele, int ncomponents,
    const int numdf, std::vector<double>& solution, Teuchos::RCP<Epetra_Vector> ghostedData,
    const int from, const bool fillzeros) const
{
  const unsigned NUMNODES = DRT::UTILS::DisTypeToNumNodePerEle<nurbs_type>::numNodePerElement;
  const unsigned DIM = DRT::UTILS::DisTypeToDim<nurbs_type>::dim;

  const Teuchos::RCP<const DRT::Discretization> dis = field_->discretization();
  std::vector<int> nodedofs;

  const DRT::Element::DiscretizationType mapped_dis_type =
      MapNurbsDisTypeToLagrangeDisType(nurbs_type);

  const std::vector<int>& numbering =
      DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(mapped_dis_type).second;

  LINALG::Matrix<NUMNODES, 1> weights;
  const DRT::Node* const* nodes = ele->Nodes();
  for (unsigned inode = 0; inode < NUMNODES; inode++)
  {
    const DRT::NURBS::ControlPoint* cp =
        dynamic_cast<const DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }
  const DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(dis.get());

  if (not nurbsdis) dserror("Cast to NURBS discretization failed.\n");

  std::vector<Epetra_SerialDenseVector> myknots(3);
  const bool zero_ele = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, ele->Id());

  if (zero_ele) return;

  for (unsigned n = 0; n < NUMNODES; ++n)
  {
    LINALG::Matrix<NUMNODES, 1> funct;
    LINALG::Matrix<DIM, NUMNODES> deriv;
    LINALG::Matrix<3, 1> gpa;

    gpa = DRT::UTILS::getNodeCoordinates(numbering[n], mapped_dis_type);

    DRT::NURBS::UTILS::nurbs_get_funct_deriv(funct, deriv, gpa, myknots, weights, nurbs_type);

    double val[numdf];
    std::fill(val, val + numdf, 0.0);

    for (unsigned m = 0; m < NUMNODES; ++m)
    {
      nodedofs.clear();
      dis->Dof(ele->Nodes()[m], nodedofs);
      for (int d = 0; d < numdf; ++d)
      {
        const int lid = ghostedData->Map().LID(nodedofs[d + from]);
        if (lid > -1)
          val[d] += funct(m) * ((*ghostedData)[lid]);
        else
        {
          if (fillzeros)
            val[d] += 0.;
          else
            dserror("received illegal dof local id: %d", lid);
        }
      }
    }

    for (int d = 0; d < numdf; ++d) solution.push_back(val[d]);

    for (int d = numdf; d < ncomponents; ++d) solution.push_back(0.);
  }
}


void PostVtuWriter::WriteDofResultStepBeamEle(const DRT::ELEMENTS::Beam3Base* beamele,
    const int& ncomponents, const int& numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_Vector>& ghostedData, const int& from, const bool fillzeros)
{
  if (numdf != ncomponents or numdf != 3)
  {
    dserror(
        "writing of dof-based result for beams with Hermite centerline interpolation is "
        "restricted to centerline displacements where numdof = ncomponents = 3");
  }

  const Teuchos::RCP<DRT::Discretization> dis = this->GetField()->discretization();
  std::vector<int> nodedofs;
  std::vector<double> elementdofvals;

  for (int n = 0; n < beamele->NumNode(); ++n)
  {
    nodedofs.clear();

    // local storage position of desired dof gid
    dis->Dof(beamele->Nodes()[n], nodedofs);

    for (std::vector<int>::const_iterator it = nodedofs.begin(); it != nodedofs.end(); ++it)
    {
      const int lid = ghostedData->Map().LID(*it);
      if (lid > -1)
        elementdofvals.push_back((*ghostedData)[lid]);
      else
      {
        if (fillzeros)
          elementdofvals.push_back(0.);
        else
          dserror("received illegal dof local id: %d", lid);
      }
    }
  }

  /* visualize the beam centerline as a sequence of straight line segments (POLY_LINE)
   * which is supported as vtkCellType number 4 (see also list in GetVtkElementType)
   * loop over the chosen visualization points (equidistant distribution in the element
   * parameter space xi \in [-1,1] ) and determine their interpolated initial positions r */
  LINALG::Matrix<3, 1> pos, refpos;
  double xi = 0.0;

  for (unsigned int i = 0; i < BEAMSVTUVISUALSUBSEGMENTS + 1; ++i)
  {
    pos.Clear();
    refpos.Clear();
    xi = -1.0 + i * 2.0 / BEAMSVTUVISUALSUBSEGMENTS;

    // let the beam element do the interpolation
    beamele->GetRefPosAtXi(refpos, xi);
    beamele->GetPosAtXi(pos, xi, elementdofvals);

    for (int d = 0; d < 3; ++d) solution.push_back(pos(d) - refpos(d));
  }
}

void PostVtuWriter::WriteNodalResultStepNurbsEle(const DRT::Element* ele, int ncomponents,
    const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_MultiVector> ghostedData) const
{
  switch (ele->Shape())
  {
    case DRT::Element::nurbs2:
    {
      WriteNodalResultStepNurbsEle<DRT::Element::nurbs2>(
          ele, ncomponents, numdf, solution, ghostedData);
      break;
    }
    case DRT::Element::nurbs3:
    {
      WriteNodalResultStepNurbsEle<DRT::Element::nurbs3>(
          ele, ncomponents, numdf, solution, ghostedData);
      break;
    }
    case DRT::Element::nurbs4:
    {
      WriteNodalResultStepNurbsEle<DRT::Element::nurbs4>(
          ele, ncomponents, numdf, solution, ghostedData);
      break;
    }
    case DRT::Element::nurbs9:
    {
      WriteNodalResultStepNurbsEle<DRT::Element::nurbs9>(
          ele, ncomponents, numdf, solution, ghostedData);
      break;
    }
    case DRT::Element::nurbs8:
    {
      WriteNodalResultStepNurbsEle<DRT::Element::nurbs8>(
          ele, ncomponents, numdf, solution, ghostedData);
      break;
    }
    case DRT::Element::nurbs27:
    {
      WriteNodalResultStepNurbsEle<DRT::Element::nurbs27>(
          ele, ncomponents, numdf, solution, ghostedData);
      break;
    }
    default:
      dserror("VTK output not yet implemented for given NURBS element");
      break;
  }  // end switch ele shape
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType nurbs_type>
void PostVtuWriter::WriteNodalResultStepNurbsEle(const DRT::Element* ele, int ncomponents,
    const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_MultiVector> ghostedData) const
{
  const unsigned NUMNODES = DRT::UTILS::DisTypeToNumNodePerEle<nurbs_type>::numNodePerElement;
  const unsigned DIM = DRT::UTILS::DisTypeToDim<nurbs_type>::dim;

  const Teuchos::RCP<const DRT::Discretization> dis = field_->discretization();
  std::vector<int> nodedofs;

  const DRT::Element::DiscretizationType mapped_dis_type =
      MapNurbsDisTypeToLagrangeDisType(nurbs_type);

  const std::vector<int>& numbering =
      DRT::ELEMENTS::GetVtkCellTypeFromBaciElementShapeType(mapped_dis_type).second;

  LINALG::Matrix<NUMNODES, 1> weights;
  const DRT::Node* const* nodes = ele->Nodes();
  for (unsigned inode = 0; inode < NUMNODES; inode++)
  {
    const DRT::NURBS::ControlPoint* cp =
        dynamic_cast<const DRT::NURBS::ControlPoint*>(nodes[inode]);

    weights(inode) = cp->W();
  }

  const DRT::NURBS::NurbsDiscretization* nurbsdis =
      dynamic_cast<const DRT::NURBS::NurbsDiscretization*>(dis.get());

  if (not nurbsdis) dserror("Cast to NURBS discretization failed.\n");

  std::vector<Epetra_SerialDenseVector> myknots(3);
  bool zero_ele = (*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots, ele->Id());
  if (zero_ele) return;
  double val[numdf];
  for (unsigned n = 0; n < NUMNODES; ++n)
  {
    LINALG::Matrix<NUMNODES, 1> funct;
    LINALG::Matrix<DIM, NUMNODES> deriv;
    LINALG::Matrix<3, 1> gpa;

    gpa = DRT::UTILS::getNodeCoordinates(numbering[n], mapped_dis_type);

    DRT::NURBS::UTILS::nurbs_get_funct_deriv(funct, deriv, gpa, myknots, weights, nurbs_type);

    for (int i = 0; i < numdf; ++i) val[i] = 0.;

    for (int idf = 0; idf < numdf; ++idf)
    {
      Epetra_Vector* column = (*ghostedData)(idf);

      for (unsigned m = 0; m < NUMNODES; ++m)
      {
        int lid = ghostedData->Map().LID(ele->Nodes()[m]->Id());
        if (lid > -1)
          val[idf] += funct(m) * (*column)[lid];
        else
          dserror("received illegal node local id: %d", lid);
      }
    }

    for (int d = 0; d < numdf; ++d) solution.push_back(val[d]);

    for (int d = numdf; d < ncomponents; ++d) solution.push_back(0.);
  }
}

/*----------------------------------------------------------------------------*/
template void PostVtuWriter::WriteGeoNurbsEle<DRT::Element::nurbs2>(const DRT::Element* ele,
    std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
    std::vector<double>& coordinates) const;
template void PostVtuWriter::WriteGeoNurbsEle<DRT::Element::nurbs3>(const DRT::Element* ele,
    std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
    std::vector<double>& coordinates) const;
template void PostVtuWriter::WriteGeoNurbsEle<DRT::Element::nurbs4>(const DRT::Element* ele,
    std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
    std::vector<double>& coordinates) const;
template void PostVtuWriter::WriteGeoNurbsEle<DRT::Element::nurbs9>(const DRT::Element* ele,
    std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
    std::vector<double>& coordinates) const;
template void PostVtuWriter::WriteGeoNurbsEle<DRT::Element::nurbs8>(const DRT::Element* ele,
    std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
    std::vector<double>& coordinates) const;
template void PostVtuWriter::WriteGeoNurbsEle<DRT::Element::nurbs27>(const DRT::Element* ele,
    std::vector<uint8_t>& celltypes, int& outNodeId, std::vector<int32_t>& celloffset,
    std::vector<double>& coordinates) const;

/*----------------------------------------------------------------------------*/
template void PostVtuWriter::WirteDofResultStepNurbsEle<DRT::Element::nurbs2>(
    const DRT::Element* ele, int ncomponents, const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_Vector> ghostedData, const int from, const bool fillzeros) const;
template void PostVtuWriter::WirteDofResultStepNurbsEle<DRT::Element::nurbs3>(
    const DRT::Element* ele, int ncomponents, const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_Vector> ghostedData, const int from, const bool fillzeros) const;
template void PostVtuWriter::WirteDofResultStepNurbsEle<DRT::Element::nurbs4>(
    const DRT::Element* ele, int ncomponents, const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_Vector> ghostedData, const int from, const bool fillzeros) const;
template void PostVtuWriter::WirteDofResultStepNurbsEle<DRT::Element::nurbs9>(
    const DRT::Element* ele, int ncomponents, const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_Vector> ghostedData, const int from, const bool fillzeros) const;
template void PostVtuWriter::WirteDofResultStepNurbsEle<DRT::Element::nurbs8>(
    const DRT::Element* ele, int ncomponents, const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_Vector> ghostedData, const int from, const bool fillzeros) const;
template void PostVtuWriter::WirteDofResultStepNurbsEle<DRT::Element::nurbs27>(
    const DRT::Element* ele, int ncomponents, const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_Vector> ghostedData, const int from, const bool fillzeros) const;

/*----------------------------------------------------------------------------*/
template void PostVtuWriter::WriteNodalResultStepNurbsEle<DRT::Element::nurbs2>(
    const DRT::Element* ele, int ncomponents, const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_MultiVector> ghostedData) const;
template void PostVtuWriter::WriteNodalResultStepNurbsEle<DRT::Element::nurbs3>(
    const DRT::Element* ele, int ncomponents, const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_MultiVector> ghostedData) const;
template void PostVtuWriter::WriteNodalResultStepNurbsEle<DRT::Element::nurbs4>(
    const DRT::Element* ele, int ncomponents, const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_MultiVector> ghostedData) const;
template void PostVtuWriter::WriteNodalResultStepNurbsEle<DRT::Element::nurbs9>(
    const DRT::Element* ele, int ncomponents, const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_MultiVector> ghostedData) const;
template void PostVtuWriter::WriteNodalResultStepNurbsEle<DRT::Element::nurbs8>(
    const DRT::Element* ele, int ncomponents, const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_MultiVector> ghostedData) const;
template void PostVtuWriter::WriteNodalResultStepNurbsEle<DRT::Element::nurbs27>(
    const DRT::Element* ele, int ncomponents, const int numdf, std::vector<double>& solution,
    Teuchos::RCP<Epetra_MultiVector> ghostedData) const;
