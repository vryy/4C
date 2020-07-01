/*----------------------------------------------------------------------*/
/*! \file

\brief VTI filter


\level 2
*/
/*----------------------------------------------------------------------*/

#include "post_drt_vti_writer.H"

#include <sstream>

#include <Epetra_Vector.h>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../post_drt_common/post_drt_common.H"

#ifdef TOL_N
#undef TOL_N
#endif
#define TOL_N 1.0e-14

template <class T>
struct less_tol
{
  bool operator()(const T& x, const T& y) const { return x < y - TOL_N; }
};

PostVtiWriter::PostVtiWriter(PostField* field, const std::string& filename)
    : PostVtkWriter(field, filename)
{
  for (size_t i = 0; i < sizeof(origin_) / sizeof(origin_[0]); ++i) origin_[i] = 0.0;
  for (size_t i = 0; i < sizeof(spacing_) / sizeof(spacing_[0]); ++i) spacing_[i] = 0.0;
  for (size_t i = 0; i < sizeof(globalextent_) / sizeof(globalextent_[0]); ++i)
    globalextent_[i] = 0;
  for (size_t i = 0; i < sizeof(localextent_) / sizeof(localextent_[0]); ++i) localextent_[i] = 0;
}


const std::string& PostVtiWriter::WriterString() const
{
  static std::string name("ImageData");
  return name;
}

const std::string& PostVtiWriter::WriterOpeningTag() const
{
  static std::string tag;
  std::stringstream stream;
  stream << "<ImageData WholeExtent=\"" << localextent_[0] << " " << localextent_[1] << " "
         << localextent_[2] << " " << localextent_[3] << " " << localextent_[4] << " "
         << localextent_[5] << "\" Origin=\"" << origin_[0] << " " << origin_[1] << " "
         << origin_[2] << "\" Spacing=\"" << spacing_[0] << " " << spacing_[1] << " " << spacing_[2]
         << "\">";
  tag.assign(stream.str());

  return tag;
}

const std::string& PostVtiWriter::WriterPOpeningTag() const
{
  static std::string tag;
  std::stringstream stream;
  stream << "<PImageData WholeExtent=\"" << globalextent_[0] << " " << globalextent_[1] << " "
         << globalextent_[2] << " " << globalextent_[3] << " " << globalextent_[4] << " "
         << globalextent_[5] << "\" GhostLevel=\"0\" Origin=\"" << origin_[0] << " " << origin_[1]
         << " " << origin_[2] << "\" Spacing=\"" << spacing_[0] << " " << spacing_[1] << " "
         << spacing_[2] << "\">";
  tag.assign(stream.str());

  return tag;
}

const std::vector<std::string>& PostVtiWriter::WriterPPieceTags() const
{
  static std::vector<std::string> tags;
  tags.clear();

  int allextents[numproc_][6];
  field_->problem()->comm()->GatherAll((int*)localextent_, (int*)allextents, 6);

  if (myrank_ == 0)
  {
    for (size_t i = 0; i < numproc_; ++i)
    {
      std::stringstream stream;
      stream << "<Piece Extent=\"" << allextents[i][0] << " " << allextents[i][1] << " "
             << allextents[i][2] << " " << allextents[i][3] << " " << allextents[i][4] << " "
             << allextents[i][5] << "\" Source=\"" << filenamebase_ << "-" << std::setfill('0')
             << std::setw(npdigits_) << i << ".vti\"/>";
      tags.push_back(stream.str());
    }
  }

  return tags;
}

const std::string& PostVtiWriter::WriterSuffix() const
{
  static std::string name(".vti");
  return name;
}

const std::string& PostVtiWriter::WriterPSuffix() const
{
  static std::string name(".pvti");
  return name;
}



void PostVtiWriter::WriteGeo()
{
  // There is no such thing as geometry for ImageData, however we need to prepare some things

  // start the piece
  currentout_ << "<Piece Extent=\"" << localextent_[0] << " " << localextent_[1] << " "
              << localextent_[2] << " " << localextent_[3] << " " << localextent_[4] << " "
              << localextent_[5] << "\">\n";
  return;
}



void PostVtiWriter::WriteDofResultStep(std::ofstream& file, const Teuchos::RCP<Epetra_Vector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
    const std::string& groupname, const std::string& name, const int numdf, const int from,
    const bool fillzeros)
{
  if (myrank_ == 0 && timestep_ == 0) std::cout << "writing dof-based field " << name << std::endl;

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  // Here is the only thing we need to do for parallel computations: We need read access to all dofs
  // on the row elements, so need to get the DofColMap to have this access
  const Epetra_Map* colmap = dis->DofColMap(0);
  const Epetra_BlockMap& vecmap = data->Map();

  // TODO: wic, once the vtu pressure writer is fixed, apply the same solution here.
  const int offset = (fillzeros) ? 0 : vecmap.MinAllGID() - dis->DofRowMap(0)->MinAllGID();

  Teuchos::RCP<Epetra_Vector> ghostedData;
  if (colmap->SameAs(vecmap))
    ghostedData = data;
  else
  {
    ghostedData = LINALG::CreateVector(*colmap, false);
    LINALG::Export(*data, *ghostedData);
  }

  const int ncomponents = (numdf > 1 && numdf == field_->problem()->num_dim()) ? 3 : numdf;

  std::vector<double> solution;
  if (fillzeros)
    solution.resize(ncomponents * ((localextent_[1] - localextent_[0] + 1) *
                                      (localextent_[3] - localextent_[2] + 1) *
                                      (localextent_[5] - localextent_[4] + 1)),
        0.0);
  else
    solution.resize(ncomponents * ((localextent_[1] - localextent_[0] + 1) *
                                      (localextent_[3] - localextent_[2] + 1) *
                                      (localextent_[5] - localextent_[4] + 1)));

  std::vector<int> nodedofs;
  for (int e = 0; e < dis->NumMyColElements(); ++e)
  {
    const DRT::Element* ele = dis->lColElement(e);

    for (int n = 0; n < ele->NumNode(); ++n)
    {
      nodedofs.clear();

      // node gid for insertion in solution
      const int ngid = ele->Nodes()[n]->Id();
      const int inpos = ncomponents * idmapping_.find(ngid)->second;

      // local storage position of desired dof gid
      dis->Dof(ele->Nodes()[n], nodedofs);

      for (int d = 0; d < numdf; ++d)
      {
        const int lid = ghostedData->Map().LID(nodedofs[d + from] + offset);
        if (lid > -1)
        {
          solution[inpos + d] = (*ghostedData)[lid];
        }
        else
        {
          if (!fillzeros) dserror("received illegal dof local id: %d", lid);
        }
      }
      for (int d = numdf; d < ncomponents; ++d) solution[inpos + d] = 0.0;
    }
  }

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

  return;
}



void PostVtiWriter::WriteNodalResultStep(std::ofstream& file,
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

  const int ncomponents = (numdf > 1 && numdf == field_->problem()->num_dim()) ? 3 : numdf;

  std::vector<double> solution(ncomponents * ((localextent_[1] - localextent_[0] + 1) *
                                                 (localextent_[3] - localextent_[2] + 1) *
                                                 (localextent_[5] - localextent_[4] + 1)));

  for (int e = 0; e < dis->NumMyColElements(); ++e)
  {
    const DRT::Element* ele = dis->lColElement(e);

    for (int n = 0; n < ele->NumNode(); ++n)
    {
      const int gid = ele->Nodes()[n]->Id();
      const int inpos = ncomponents * idmapping_.find(gid)->second;

      for (int idf = 0; idf < numdf; ++idf)
      {
        Epetra_Vector* column = (*ghostedData)(idf);

        const int lid = ghostedData->Map().LID(gid);

        if (lid > -1)
        {
          solution[inpos + idf] = (*column)[lid];
        }
        else
        {
          dserror("received illegal node local id: %d", lid);
        }
      }
      for (int d = numdf; d < ncomponents; ++d) solution[inpos + d] = 0.;
    }
  }

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

  return;
}



void PostVtiWriter::WriteElementResultStep(std::ofstream& file,
    const Teuchos::RCP<Epetra_MultiVector>& data,
    std::map<std::string, std::vector<std::ofstream::pos_type>>& resultfilepos,
    const std::string& groupname, const std::string& name, const int numdf, const int from)
{
  if (myrank_ == 0 && timestep_ == 0)
    std::cout << "writing element-based field " << name << std::endl;

  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();

  const int ncomponents = (numdf > 1 && numdf == field_->problem()->num_dim()) ? 3 : numdf;

  std::vector<double> solution(
      ncomponents * ((localextent_[1] - localextent_[0]) * (localextent_[3] - localextent_[2]) *
                        (localextent_[5] - localextent_[4])),
      0.0);

  const int numcol = data->NumVectors();
  if (numdf + from > numcol) dserror("violated column range of Epetra_MultiVector: %d", numcol);

  Teuchos::RCP<Epetra_MultiVector> importedData;
  if (dis->ElementColMap()->SameAs(data->Map()))
    importedData = data;
  else
  {
    importedData =
        Teuchos::rcp(new Epetra_MultiVector(*dis->ElementColMap(), data->NumVectors(), false));
    LINALG::Export(*data, *importedData);
  }

  for (int e = 0; e < dis->NumMyColElements(); ++e)
  {
    const DRT::Element* ele = dis->lColElement(e);
    const int egid = ele->Id();
    const int inpos = ncomponents * (eidmapping_.find(egid)->second);
    for (int d = 0; d < numdf; ++d)
    {
      Epetra_Vector* column = (*importedData)(d + from);
      solution[inpos + d] = (*column)[e];
    }
  }

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

  return;
}



void PostVtiWriter::WriterPrepTimestep()
{
  const Teuchos::RCP<DRT::Discretization> dis = field_->discretization();
  // collect all possible values of the x-, y- and z-coordinate
  typedef std::set<double, less_tol<double>> set_tol;
  set_tol collected_coords[3];
  for (int n = 0; n < dis->NumMyColNodes(); ++n)
  {
    const double* coord = dis->lColNode(n)->X();
    for (int i = 0; i < 3; ++i) collected_coords[i].insert(coord[i]);
  }

  // determine local and global domain size
  double lorigin[3], lextent[3], gorigin[3], gextent[3];
  for (int i = 0; i < 3; ++i)
  {
    lorigin[i] = *collected_coords[i].begin();
    lextent[i] = *collected_coords[i].rbegin();
  }
  field_->discretization()->Comm().MinAll(lorigin, gorigin, sizeof(lorigin) / sizeof(lorigin[0]));
  field_->discretization()->Comm().MaxAll(lextent, gextent, sizeof(lextent) / sizeof(lextent[0]));

  // determine spacing and check whether it is consistent for ImageData
  for (int i = 0; i < 3; ++i)
    spacing_[i] = (*(collected_coords[i].rbegin()) - *(collected_coords[i].begin())) /
                  (collected_coords[i].size() - 1);
  for (int i = 0; i < 3; ++i)
  {
    int k = 0;
    for (set_tol::const_iterator it = collected_coords[i].begin(); it != collected_coords[i].end();
         ++it, ++k)
      if (*it - lorigin[i] - k * spacing_[i] > TOL_N || *it - lorigin[i] - k * spacing_[i] < -TOL_N)
        dserror(
            "The mesh is not a uniform rectangular grid: grid coordinate[%d]: %lf, "
            "node coordinate[%d]: %lf, difference: %e",
            i, lorigin[i] + k * spacing_[i], i, *it, *it - lorigin[i] - k * spacing_[i]);
  }

  // determine extents
  for (int i = 0; i < 3; ++i)
  {
    localextent_[0 + 2 * i] = round((lorigin[i] - gorigin[i]) / spacing_[i]);
    localextent_[1 + 2 * i] = round((lextent[i] - gorigin[i]) / spacing_[i]);
    globalextent_[0 + 2 * i] = 0;
    globalextent_[1 + 2 * i] = round((gextent[i] - gorigin[i]) / spacing_[i]);
  }

  // finally create the mapping
  int nx = localextent_[1] - localextent_[0] + 1;
  int ny = localextent_[3] - localextent_[2] + 1;
  idmapping_.clear();
  for (int n = 0; n < dis->NumMyColNodes(); ++n)
  {
    const double* coord = dis->lColNode(n)->X();
    int i = round((coord[0] - lorigin[0]) / spacing_[0]);
    int j = round((coord[1] - lorigin[1]) / spacing_[1]);
    int k = round((coord[2] - lorigin[2]) / spacing_[2]);
    idmapping_[dis->NodeColMap()->GID(n)] = (k * ny + j) * nx + i;
  }

  // create element id mapping by the coordinates of the lowest node coordinate
  nx = localextent_[1] - localextent_[0];
  ny = localextent_[3] - localextent_[2];
  eidmapping_.clear();
  for (int e = 0; e < dis->NumMyColElements(); ++e)
  {
    const DRT::Element* ele = dis->lColElement(e);
    double mincoord[] = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()};
    for (int n = 0; n < ele->NumNode(); ++n)
    {
      const double* coord = ele->Nodes()[n]->X();
      mincoord[0] = std::min(mincoord[0], coord[0]);
      mincoord[1] = std::min(mincoord[1], coord[1]);
      mincoord[2] = std::min(mincoord[2], coord[2]);
    }
    int i = round((mincoord[0] - lorigin[0]) / spacing_[0]);
    int j = round((mincoord[1] - lorigin[1]) / spacing_[1]);
    int k = round((mincoord[2] - lorigin[2]) / spacing_[2]);
    eidmapping_[dis->ElementColMap()->GID(e)] = (k * ny + j) * nx + i;
  }
  return;
}
