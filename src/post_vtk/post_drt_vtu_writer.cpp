/*----------------------------------------------------------------------*/
/*!
\file post_drt_vtu_writer.cpp

\brief VTU filter

\level 2

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/


#include "post_drt_vtu_writer.H"

#include <sstream>
#include <boost/static_assert.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../post_drt_common/post_drt_common.H"

#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_lib/drt_element.H"

#include "../drt_beam3/beam3_base.H"


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


VtuWriter::VtuWriter(PostField* field, const std::string &filename) :
    VtkWriter(field, filename)
{
  BOOST_STATIC_ASSERT_MSG( 29 == DRT::Element::max_distype, "The number of element types defined by DRT::Element::DiscretizationType does not match the number of element types supported by the post vtu filter.") BACI_ATTRIBUTE_UNUSED;

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
  for (int e=0; e<nelements; ++e)
  {
    const DRT::Element* ele = dis->lRowElement(e);
    // check for beam element that potentially needs special treatment due to Hermite interpolation
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    if (ele->IsNurbsElement())
    {
      WriteGeoNurbsEle(ele,celltypes,outNodeId,celloffset,coordinates);
    }
    else if (beamele!=NULL)
    {
      WriteGeoBeamEle(beamele,celltypes,outNodeId,celloffset,coordinates);
    }
    else
    {
      celltypes.push_back(GetVtkElementType(ele->Shape()).first);
      const std::vector<int> &numbering = GetVtkElementType(ele->Shape()).second;
      const DRT::Node* const* nodes = ele->Nodes();
      for (int n=0; n<ele->NumNode(); ++n)
        for (int d=0; d<3; ++d)
          coordinates.push_back(nodes[numbering[n]]->X()[d]);
      outNodeId += ele->NumNode();
      celloffset.push_back(outNodeId);
    }
  }
  dsassert((int)coordinates.size() == 3*outNodeId, "internal error");

  // step 1: write node coordinates into file
  currentout_ << "<Piece NumberOfPoints=\"" << outNodeId
      <<"\" NumberOfCells=\"" << nelements << "\" >\n"
      << "  <Points>\n"
      << "    <DataArray type=\"Float64\" NumberOfComponents=\"3\"";

  if(write_binary_output_)
  {
    currentout_ << " format=\"binary\">\n";
    LIBB64:: writeCompressedBlock(coordinates, currentout_);
  }
  else
  {
    currentout_ << " format=\"ascii\">\n";
    int counter = 1;
    for (std::vector<double>::const_iterator it = coordinates.begin(); it != coordinates.end(); ++it)
    {
      currentout_ << std::setprecision(15) << std::scientific << *it << " ";
      // dimension is hard coded to three, thus
      if(counter%3==0)
        currentout_ << '\n';
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
  if(write_binary_output_)
    currentout_ << " format=\"binary\">\n";
  else
  currentout_ << " format=\"ascii\">\n";

  {
    std::vector<int32_t> connectivity;
    connectivity.reserve(outNodeId);
    for (int i=0; i<outNodeId; ++i)
      connectivity.push_back(i);
    if(write_binary_output_)
      LIBB64::writeCompressedBlock(connectivity, currentout_);
    else
    {
      for (std::vector<int32_t>::const_iterator it = connectivity.begin(); it != connectivity.end(); ++it)
        currentout_ << *it << " ";
    }

  }
  currentout_ << "    </DataArray>\n";

  // step 3: write start indices for individual cells
  currentout_ << "    <DataArray type=\"Int32\" Name=\"offsets\"";
  if(write_binary_output_)
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
  if(write_binary_output_)
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

  // For parallel computations, we need to access all dofs on the elements, including the
  // nodes owned by other processors. Therefore, we need to import that data here.
  const Epetra_BlockMap& vecmap = data->Map();
  const Epetra_Map* colmap = dis->DofColMap(0);

  int offset = vecmap.MinAllGID() - dis->DofRowMap()->MinAllGID();
  if(fillzeros)
    offset=0;

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
    for (int i=0; i<vecmap.NumMyElements(); ++i)
      gids[i] = vecmap.MyGlobalElements()[i] - offset;
    Teuchos::RCP<Epetra_Map> rowmap = Teuchos::rcp(new Epetra_Map(vecmap.NumGlobalElements(),
                                       vecmap.NumMyElements(),
                                       &gids[0],
                                       0,
                                       vecmap.Comm()));
    Teuchos::RCP<Epetra_Vector> dofvec = LINALG::CreateVector(*rowmap, false);
    for (int i=0; i<vecmap.NumMyElements(); ++i)
      (*dofvec)[i] = (*data)[i];

    ghostedData = LINALG::CreateVector(*colmap,true);
    LINALG::Export(*dofvec,*ghostedData);
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
  for (int e=0; e<dis->NumMyRowElements(); ++e)
  {
    const DRT::Element* ele = dis->lRowElement(e);
    // check for beam element that potentially needs special treatment due to Hermite interpolation
    const DRT::ELEMENTS::Beam3Base* beamele = dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele);

    if (ele->IsNurbsElement())
    {
      WirteDofResultStepNurbsEle(ele,ncomponents,numdf,solution,ghostedData,from,fillzeros);
    }
    else if (beamele!=NULL)
    {
      WriteDofResultStepBeamEle(beamele,ncomponents,numdf,solution,ghostedData,from,fillzeros);
    }
    else
    {
      const std::vector<int> &numbering = GetVtkElementType(ele->Shape()).second;
      for (int n=0; n<ele->NumNode(); ++n)
      {
        nodedofs.clear();

        // local storage position of desired dof gid
        dis->Dof(ele->Nodes()[numbering[n]], nodedofs);

        for (int d=0; d<numdf; ++d)
        {
          const int lid = ghostedData->Map().LID(nodedofs[d+from]);
          if (lid > -1)
            solution.push_back((*ghostedData)[lid]);
          else
          {
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

  } // loop over all elements
  dsassert((int)solution.size() == ncomponents*nnodes, "internal error");

  // start the scalar fields that will later be written
  if (currentPhase_ == INIT)
  {
    currentout_ << "  <PointData>\n"; // Scalars=\"scalars\">\n";
    if (myrank_ == 0)
    {
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
  else
  {
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

  for (int e=0; e<dis->NumMyRowElements(); ++e)
  {
  const DRT::Element* ele = dis->lRowElement(e);

  if (ele->IsNurbsElement()==false)
  {
    const std::vector<int> &numbering = GetVtkElementType(ele->Shape()).second;
    for (int n=0; n<ele->NumNode(); ++n)
    {
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
  else
    WriteNodalResultStepNurbsEle(ele,ncomponents,numdf,solution,ghostedData);
  } // loop over all elements
  dsassert((int)solution.size() == ncomponents*nnodes, "internal error");


  // start the scalar fields that will later be written
  if (currentPhase_ == INIT)
  {
    currentout_ << "  <PointData>\n"; // Scalars=\"scalars\">\n";
    if (myrank_ == 0)
    {
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

  for (int e=0; e<dis->NumMyRowElements(); ++e)
  {
    const DRT::Element* ele = dis->lRowElement(e);
    for (int n=0; n<ele->NumNode(); ++n)
    {
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
    currentout_ << "  <CellData>\n"; // Scalars=\"scalars\">\n";
    if (myrank_ == 0)
    {
      currentmasterout_ << "    <PCellData>\n"; // Scalars=\"scalars\">\n";
    }
    currentPhase_ = CELLS;
  }

  if (currentPhase_ != CELLS)
    dserror("Cannot write cell data at this stage. Most likely cell and point data fields are mixed.");

  this->WriteSolutionVector(solution, ncomponents, name, file);
}

void
VtuWriter::WriteGeoNurbsEle(const DRT::Element* ele,std::vector<uint8_t>& celltypes,
      int& outNodeId,std::vector<int32_t>& celloffset,std::vector<double>& coordinates)
{
  Teuchos::RCP<DRT::Discretization> dis = this->GetField()->discretization();

  switch (ele->Shape())
  {
  case DRT::Element::nurbs27:
  {
    const std::vector<int> &numbering = GetVtkElementType(DRT::Element::hex27).second;

    celltypes.push_back(GetVtkElementType(DRT::Element::hex27).first);

    LINALG::Matrix<27,1> weights;
    const DRT::Node* const* nodes = ele->Nodes();
    for (int inode=0; inode<27; inode++)
    {
      const DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<const DRT::NURBS::ControlPoint* > (nodes[inode]);

      weights(inode) = cp->W();
    }
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*dis));

    if(nurbsdis==NULL)
    {
      dserror("So_nurbs27 appeared in non-nurbs discretisation\n");
    }
    std::vector<Epetra_SerialDenseVector> myknots(3);
    bool zero_ele=(*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,ele->Id());
    if (zero_ele)
      return;
    for (int n=0; n<ele->NumNode(); ++n)
    {
      LINALG::Matrix<27,1> funct;
      LINALG::Matrix<3,27> deriv;
      LINALG::Matrix<3,1> gpa;
      gpa = DRT::UTILS::getNodeCoordinates(numbering[n],DRT::Element::hex27);
      DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv
        (funct                 ,
         deriv                 ,
         gpa                   ,
         myknots               ,
         weights               ,
         DRT::Element::nurbs27);
      double X[3]={0.,0.,0.};
      for (int i=0;i<3;++i)
        for (int m=0;m<27;++m)
          X[i] += funct(m)*(nodes[m]->X()[i]);

      for (int d=0; d<3; ++d)
        coordinates.push_back(X[d]);
    }
    outNodeId += ele->NumNode();
    celloffset.push_back(outNodeId);
  }
  break;

  default:
    dserror("VTK output not yet implemented for given NURBS element");
    break;
  }

  return;
}

void
VtuWriter::WriteGeoBeamEle(const DRT::ELEMENTS::Beam3Base* beamele,
                           std::vector<uint8_t>& celltypes,
                           int& outNodeId,
                           std::vector<int32_t>& celloffset,
                           std::vector<double>& coordinates)
{
  /* visualize the beam centerline as a sequence of straight line segments (POLY_LINE)
   * which is supported as vtkCellType number 4 (see also list in GetVtkElementType) */
  celltypes.push_back(4);

  /* loop over the chosen visualization points (equidistant distribution in the element
   * parameter space xi \in [-1,1] ) and determine their interpolated initial positions r */
  LINALG::Matrix<3,1> r;
  double xi=0.0;

  for (unsigned int i=0; i<BEAMSVTUVISUALSUBSEGMENTS+1; ++i)
  {
    r.Clear();
    xi= -1.0 + i*2.0/BEAMSVTUVISUALSUBSEGMENTS;

    beamele->GetRefPosAtXi(r,xi);

    for (int d=0; d<3; ++d)
      coordinates.push_back(r(d));
  }

  outNodeId += BEAMSVTUVISUALSUBSEGMENTS+1;
  celloffset.push_back(outNodeId);
}

void
VtuWriter::WirteDofResultStepNurbsEle(const DRT::Element* ele, int ncomponents,const int numdf,
      std::vector<double>& solution,Teuchos::RCP<Epetra_Vector> ghostedData,
      const int from, const bool fillzeros)
{
  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();
  std::vector<int> nodedofs;

  switch (ele->Shape())
  {
  case DRT::Element::nurbs27:
  {
    const std::vector<int> &numbering = GetVtkElementType(DRT::Element::hex27).second;

    LINALG::Matrix<27,1> weights;
    const DRT::Node* const* nodes = ele->Nodes();
    for (int inode=0; inode<27; inode++)
    {
      const DRT::NURBS::ControlPoint* cp
      =
          dynamic_cast<const DRT::NURBS::ControlPoint* > (nodes[inode]);

      weights(inode) = cp->W();
    }
    DRT::NURBS::NurbsDiscretization* nurbsdis
    =
        dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*dis));

    if(nurbsdis==NULL)
    {
      dserror("So_nurbs27 appeared in non-nurbs discretisation\n");
    }
    std::vector<Epetra_SerialDenseVector> myknots(3);
    bool zero_ele=(*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,ele->Id());
    if (zero_ele)
      return;
    for (int n=0; n<ele->NumNode(); ++n)
    {
      LINALG::Matrix<27,1> funct;
      LINALG::Matrix<3,27> deriv;
      LINALG::Matrix<3,1> gpa;
      gpa = DRT::UTILS::getNodeCoordinates(numbering[n],DRT::Element::hex27);
      DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv
      (funct                 ,
          deriv                 ,
          gpa                   ,
          myknots               ,
          weights               ,
          DRT::Element::nurbs27);
      double val[numdf] ; for (int d=0;d<numdf;++d) val[d]=0.;
      for (int m=0;m<27;++m)
      {
        nodedofs.clear();
        dis->Dof(ele->Nodes()[m], nodedofs);
        for (int d=0; d<numdf; ++d)
        {
          const int lid = ghostedData->Map().LID(nodedofs[d+from]);
          if (lid > -1)
            val[d] += funct(m)*((*ghostedData)[lid]);
          else
          {
            if(fillzeros)
              val[d] += 0.;
            else
              dserror("received illegal dof local id: %d", lid);
          }
        }
      }
      for (int d=0;d<numdf;++d)
        solution.push_back(val[d]);
    }
    break;
  } // end case nurbs27
  default:
    dserror("VTK output not yet implemented for given NURBS element");
    break;
  }// end switch shape
  return;
}

void
VtuWriter::WriteDofResultStepBeamEle(const DRT::ELEMENTS::Beam3Base* beamele,
                                    const int& ncomponents,
                                    const int& numdf,
                                    std::vector<double>& solution,
                                    Teuchos::RCP<Epetra_Vector>& ghostedData,
                                    const int& from,
                                    const bool fillzeros)
{
  if (numdf != ncomponents or numdf != 3)
  {
    dserror("writing of dof-based result for beams with Hermite centerline interpolation is "
            "restricted to centerline displacements where numdof = ncomponents = 3");
  }

  const Teuchos::RCP<DRT::Discretization> dis = this->GetField()->discretization();
  std::vector<int> nodedofs;
  std::vector<double> elementdofvals;

  for (int n=0; n<beamele->NumNode(); ++n)
  {
    nodedofs.clear();

    // local storage position of desired dof gid
    dis->Dof(beamele->Nodes()[n], nodedofs);

    for (std::vector<int>::const_iterator it=nodedofs.begin(); it!=nodedofs.end(); ++it)
    {
      const int lid = ghostedData->Map().LID(*it);
      if (lid > -1)
        elementdofvals.push_back((*ghostedData)[lid]);
      else
      {
        if(fillzeros)
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
  LINALG::Matrix<3,1> pos, refpos;
  double xi=0.0;

  for (unsigned int i=0; i<BEAMSVTUVISUALSUBSEGMENTS+1; ++i)
  {
    pos.Clear();
    refpos.Clear();
    xi= -1.0 + i*2.0/BEAMSVTUVISUALSUBSEGMENTS;

    // let the beam element do the interpolation
    beamele->GetRefPosAtXi(refpos,xi);
    beamele->GetPosAtXi(pos,xi,elementdofvals);

    for (int d=0; d<3; ++d)
      solution.push_back(pos(d) - refpos(d));
  }
}

void
VtuWriter::WriteNodalResultStepNurbsEle(const DRT::Element* ele, int ncomponents,const int numdf,
    std::vector<double>& solution,Teuchos::RCP<Epetra_MultiVector> ghostedData)
{
  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  switch (ele->Shape())
  {
  case DRT::Element::nurbs27:
  {
  const std::vector<int> &numbering = GetVtkElementType(DRT::Element::hex27).second;

  LINALG::Matrix<27,1> weights;
  const DRT::Node* const* nodes = ele->Nodes();
  for (int inode=0; inode<27; inode++)
  {
    const DRT::NURBS::ControlPoint* cp
    =
        dynamic_cast<const DRT::NURBS::ControlPoint* > (nodes[inode]);

    weights(inode) = cp->W();
  }
  DRT::NURBS::NurbsDiscretization* nurbsdis
  =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*dis));

  if(nurbsdis==NULL)
  {
    dserror("So_nurbs27 appeared in non-nurbs discretisation\n");
  }
  std::vector<Epetra_SerialDenseVector> myknots(3);
  bool zero_ele=(*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,ele->Id());
  if (zero_ele)
    return;
  double val[numdf];
  for (int n=0; n<ele->NumNode(); ++n)
  {
    LINALG::Matrix<27,1> funct;
    LINALG::Matrix<3,27> deriv;
    LINALG::Matrix<3,1> gpa;
    gpa = DRT::UTILS::getNodeCoordinates(numbering[n],DRT::Element::hex27);
    DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv
    (funct                 ,
        deriv                 ,
        gpa                   ,
        myknots               ,
        weights               ,
        DRT::Element::nurbs27);

    for (int i=0;i<numdf;++i)val[i]=0.;
    for (int idf=0; idf<numdf; ++idf)
    {
      Epetra_Vector* column = (*ghostedData)(idf);

      for (int m=0;m<27;++m)
      {
        int lid = ghostedData->Map().LID(ele->Nodes()[m]->Id());
        if (lid>-1)
          val[idf] += funct(m)* (*column)[lid];
        else
          dserror("received illegal node local id: %d", lid);
      }
    }
    for (int d=0;d<numdf;++d)
      solution.push_back(val[d]);

    for (int d=numdf; d<ncomponents; ++d)
      solution.push_back(0.);
  }
  break;
  }// end nurbs27 case
  default:
    dserror("VTK output not yet implemented for given NURBS element");
    break;
  } // end switch ele shape
  return;
}

std::pair<uint8_t,std::vector<int> >
VtuWriter::GetVtkElementType(int bacieletype)
{
  // the VTK element types are from the documentation of vtkCellType,
  // e.g. at http://www.vtk.org/doc/nightly/html/vtkCellType_8h.html
  // this list must be kept in sync with the element types since we use this
  // for index translation
  switch (bacieletype){
  case 0: // dis_none
    return  std::pair<uint8_t,std::vector<int> > (0, std::vector<int>());
  case 1: // quad4
    return  std::pair<uint8_t,std::vector<int> > (9, make_vector<int>() << 0 << 1 << 2 << 3);
  case 2: // quad6
    return  std::pair<uint8_t,std::vector<int> > (30, make_vector<int>() << 0 << 1 << 4 << 3 << 2 << 5);
  case 3: // quad8
    return  std::pair<uint8_t,std::vector<int> > (23, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7);
  case 4: // quad9
    return  std::pair<uint8_t,std::vector<int> > (28, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8);
  case 5: // tri3
    return  std::pair<uint8_t,std::vector<int> > (5, make_vector<int>() << 0 << 1 << 2);
  case 6: // tri6
    return  std::pair<uint8_t,std::vector<int> > (22, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5);
  case 7: // hex8
    return  std::pair<uint8_t,std::vector<int> > (12, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7);
  case 8: // hex16
    return  std::pair<uint8_t,std::vector<int> > (12, make_vector<int>() << 0 << 1 << 2 << 3 << 8 << 9 << 10 << 11);
  case 9: // hex18
    return  std::pair<uint8_t,std::vector<int> > (12, make_vector<int>() << 0 << 1 << 2 << 3 << 9 << 10 << 11 << 12);
  case 10: // hex20
    return  std::pair<uint8_t,std::vector<int> > (25, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8
        << 9 << 10 << 11 << 16 << 17 << 18 << 19 << 12 << 13 << 14 << 15);
  case 11: // hex27
    return  std::pair<uint8_t,std::vector<int> > (29, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8
        << 9 << 10 << 11 << 16 << 17 << 18 << 19 << 12 << 13 << 14 << 15
        << 24 << 22 << 21 << 23 << 20 << 25 << 26);
  case 12: // tet4
    return  std::pair<uint8_t,std::vector<int> > (10, make_vector<int>() << 0 << 1 << 2 << 3);
  case 13: // tet10
    return  std::pair<uint8_t,std::vector<int> > (24, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8 << 9);
  case 14: // wedge6
    return  std::pair<uint8_t,std::vector<int> > (13, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5);
  case 15: // wedge15
    return  std::pair<uint8_t,std::vector<int> > (26, make_vector<int>() << 0 << 1 << 2 << 3 << 4 << 5 << 6 << 7 << 8
        << 12 << 13 << 14 << 9 << 10 << 11);
  case 16: // pyramid5
    return  std::pair<uint8_t,std::vector<int> > (14, make_vector<int>() << 0 << 1 << 2 << 3 << 4);
  case 17: // line2
    return  std::pair<uint8_t,std::vector<int> > (3, make_vector<int>() << 0 << 1);
  case 18: // line3
    return  std::pair<uint8_t,std::vector<int> > (21, make_vector<int>() << 0 << 1 << 2);
  case 19: // line4
    return  std::pair<uint8_t,std::vector<int> > (35, make_vector<int>() << 0 << 1 << 2 << 3);
  case 20: // line5 -> mapped onto line4
    return  std::pair<uint8_t,std::vector<int> > (35, make_vector<int>() << 0 << 1 << 2 << 3);
  case 21: // line6 -> mapped onto line4
    return  std::pair<uint8_t,std::vector<int> > (35, make_vector<int>() << 0 << 1 << 2 << 3);
  case 22: // point1
    return  std::pair<uint8_t,std::vector<int> > (1, make_vector<int>() << 0);
  case 23: // nurbs2, not yet implemented
    return  std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>());
  case 24: // nurbs3, not yet implemented
    return  std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>());
  case 25: // nurbs4, not yet implemented
    return  std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>());
  case 26: // nurbs9, not yet implemented
    return  std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>());
  case 27: // nurbs8, not yet implemented
    return  std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>());
  case 28: // nurbs27, not yet implemented
    return  std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>());
  case -1: // return the number of available element types +1
    return  std::pair<uint8_t,std::vector<int> > (29, make_vector<int>() << 0);
  default:
    dserror("Unknown element type");
    return std::pair<uint8_t,std::vector<int> > (static_cast<uint8_t>(-1), std::vector<int>());
}

}
