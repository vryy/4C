/*----------------------------------------------------------------------*/
/*! \file

\brief Class containing geometric operations usually needed for the coupling of an embedded
body using a binning strategy to pre-sort the beam elements, for which an octree search needs to be
performed afterwards

\level 3

*----------------------------------------------------------------------*/
#include "baci_fbi_immersed_geometry_coupler_binning.hpp"

#include "baci_binstrategy.hpp"
#include "baci_binstrategy_utils.hpp"
#include "baci_discretization_fem_general_extract_values.hpp"
#include "baci_lib_discret_faces.hpp"
#include "baci_lib_element.hpp"
#include "baci_lib_node.hpp"
#include "baci_linalg_fixedsizematrix.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN
/*----------------------------------------------------------------------*/

FBI::FBIBinningGeometryCoupler::FBIBinningGeometryCoupler()
    : FBIGeometryCoupler::FBIGeometryCoupler(),
      binstrategy_(),
      bintoelemap_(),
      binrowmap_(Teuchos::null)
{
}
/*----------------------------------------------------------------------*/
void FBI::FBIBinningGeometryCoupler::SetupBinning(
    std::vector<Teuchos::RCP<DRT::Discretization>>& discretizations,
    Teuchos::RCP<const Epetra_Vector> structure_displacement)
{
  Teuchos::RCP<const Epetra_Vector> disp2 =
      Teuchos::rcp(new const Epetra_Vector(*(discretizations[1]->DofColMap())));

  std::vector<Teuchos::RCP<const Epetra_Vector>> disp_vec = {structure_displacement, disp2};

  PartitionGeometry(discretizations, structure_displacement);
}
/*----------------------------------------------------------------------*/
void FBI::FBIBinningGeometryCoupler::PartitionGeometry(
    std::vector<Teuchos::RCP<DRT::Discretization>>& discretizations,
    Teuchos::RCP<const Epetra_Vector> structure_displacement)
{
  Teuchos::RCP<const Epetra_Vector> disp2 =
      Teuchos::rcp(new const Epetra_Vector(*(discretizations[1]->DofColMap())));

  std::vector<Teuchos::RCP<const Epetra_Vector>> disp_vec = {structure_displacement, disp2};

  // nodes, that are owned by a proc, are distributed to the bins of this proc
  std::vector<std::map<int, std::vector<int>>> nodesinbin(2);

  std::map<int, std::set<int>> bintorowelemap_fluid;

  binstrategy_->DistributeRowElementsToBinsUsingEleAABB(
      discretizations[0], bintoelemap_, structure_displacement);

  binstrategy_->BinDiscret()->FillComplete(false, false, false);

  std::set<int> colbins;

  // first, add default one layer ghosting

  std::vector<int> binvec(27);
  for (auto i = 0; i < binstrategy_->BinDiscret()->ElementRowMap()->NumMyElements(); ++i)
  {
    auto currbin = binstrategy_->BinDiscret()->lRowElement(i);
    int it = currbin->Id();
    {
      binstrategy_->GetNeighborAndOwnBinIds(it, binvec);
      colbins.insert(binvec.begin(), binvec.end());
      binvec.clear();
    }
  }


  // extend ghosting of bin discretization
  binstrategy_->ExtendGhostingOfBinningDiscretization(binrowmap_, colbins, true);

  // assign Elements to bins
  binstrategy_->RemoveAllElesFromBins();
  binstrategy_->AssignElesToBins(discretizations[0], bintoelemap_);
}
/*----------------------------------------------------------------------*/
void FBI::FBIBinningGeometryCoupler::UpdateBinning(
    Teuchos::RCP<DRT::Discretization>& structure_discretization,
    Teuchos::RCP<const Epetra_Vector> structure_column_displacement)
{
  binstrategy_->DistributeColElementsToBinsUsingEleAABB(
      structure_discretization, bintoelemap_, structure_column_displacement);


  // assign Elements to bins
  binstrategy_->RemoveAllElesFromBins();
  binstrategy_->AssignElesToBins(structure_discretization, bintoelemap_);
}
/*----------------------------------------------------------------------*/
void FBI::FBIBinningGeometryCoupler::Setup(
    std::vector<Teuchos::RCP<DRT::Discretization>>& discretizations,
    Teuchos::RCP<const Epetra_Vector> structure_displacement)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("FBI::FBICoupler::Setup");
  Teuchos::TimeMonitor monitor(*t);

  SetupBinning(discretizations, structure_displacement);

  FBI::FBIGeometryCoupler::Setup(discretizations, structure_displacement);
}
/*----------------------------------------------------------------------*/

Teuchos::RCP<std::map<int, std::vector<int>>> FBI::FBIBinningGeometryCoupler::Search(
    std::vector<Teuchos::RCP<DRT::Discretization>>& discretizations,
    Teuchos::RCP<const Epetra_Vector>& column_structure_displacement)
{
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("FBI::FBIBinningCoupler::Search");
  Teuchos::TimeMonitor monitor(*t);

  UpdateBinning(discretizations[0], column_structure_displacement);
  // Vector to hand elements pointers to the bridge object
  Teuchos::RCP<std::map<int, std::vector<int>>> pairids =
      Teuchos::rcp(new std::map<int, std::vector<int>>);

  pairids = FBI::FBIGeometryCoupler::Search(discretizations, column_structure_displacement);

  return pairids;
}

/*----------------------------------------------------------------------*/

void FBI::FBIBinningGeometryCoupler::ComputeCurrentPositions(DRT::Discretization& dis,
    Teuchos::RCP<std::map<int, CORE::LINALG::Matrix<3, 1>>> positions,
    Teuchos::RCP<const Epetra_Vector> disp) const
{
  positions->clear();
  std::vector<int> src_dofs(
      9);  // todo this does not work for all possible elements, does it? Variable size?
  std::vector<double> mydisp(3, 0.0);

  const Epetra_Map* bincolmap = binstrategy_->BinDiscret()->ElementColMap();
  std::vector<int> colbinvec;
  colbinvec.reserve(bincolmap->NumMyElements());

  for (int lid = 0; lid < bincolmap->NumMyElements(); ++lid)
  {
    DRT::Element* currbin = binstrategy_->BinDiscret()->lColElement(lid);
    colbinvec.push_back(currbin->Id());
  }

  std::set<DRT::Element*> beam_element_list;

  binstrategy_->GetBinContent(beam_element_list, {BINSTRATEGY::UTILS::Beam}, colbinvec, false);

  for (std::set<DRT::Element*>::iterator element = beam_element_list.begin();
       element != beam_element_list.end(); element++)
  {
    DRT::Node** node_list = (*element)->Nodes();
    unsigned int numnode = (*element)->NumNode();
    for (unsigned int i = 0; i < numnode; i++)
    {
      const DRT::Node* node = node_list[i];
      if (disp != Teuchos::null)
      {
        // get the DOF numbers of the current node
        dis.Dof(node, 0, src_dofs);
        // get the current displacements
        CORE::FE::ExtractMyValues(*disp, mydisp, src_dofs);

        for (int d = 0; d < 3; ++d) (*positions)[node->Id()](d) = node->X()[d] + mydisp.at(d);
      }
    }
  }
}

/*----------------------------------------------------------------------*/

void FBI::FBIBinningGeometryCoupler::SetBinning(Teuchos::RCP<BINSTRATEGY::BinningStrategy> binning)
{
  binstrategy_ = binning;
  binstrategy_->BinDiscret()->FillComplete(false, false, false);
  binrowmap_ = Teuchos::rcp(new Epetra_Map(*(binstrategy_->BinDiscret()->ElementRowMap())));
};

FOUR_C_NAMESPACE_CLOSE
