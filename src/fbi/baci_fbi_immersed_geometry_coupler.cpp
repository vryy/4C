/*----------------------------------------------------------------------*/
/*! \file

\brief Class containing geometric operations usually needed for the coupling of an embedded body.
The current implementation does not scale at all!

\level 3

*----------------------------------------------------------------------*/
#include "baci_fbi_immersed_geometry_coupler.H"

#include "baci_binstrategy.H"
#include "baci_binstrategy_utils.H"
#include "baci_discretization_geometry_searchtree.H"
#include "baci_discretization_geometry_searchtree_service.H"
#include "baci_global_data.H"
#include "baci_inpar_fluid.H"
#include "baci_lib_discret_faces.H"
#include "baci_lib_element.H"
#include "baci_lib_node.H"
#include "baci_lib_utils.H"
#include "baci_lib_utils_parallel.H"
#include "baci_linalg_fixedsizematrix.H"
#include "baci_linalg_utils_densematrix_communication.H"

#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/

FBI::FBIGeometryCoupler::FBIGeometryCoupler()
    : fluidpositions_(new std::map<int, CORE::LINALG::Matrix<3, 1>>),
      beampositions_(new std::map<int, CORE::LINALG::Matrix<3, 1>>),
      searchtree_(new CORE::GEO::SearchTree(5)),
      searchradius_(DRT::Problem::Instance()
                        ->FBIParams()
                        .sublist("BEAM TO FLUID MESHTYING")
                        .get<double>("SEARCH_RADIUS")),
      edgebased_fluidstabilization_(false)
{
  edgebased_fluidstabilization_ =
      (INPUT::IntegralValue<INPAR::FLUID::StabType>(
           DRT::Problem::Instance()->FluidDynamicParams().sublist("RESIDUAL-BASED STABILIZATION"),
           "STABTYPE") == INPAR::FLUID::stabtype_edgebased);
}
/*----------------------------------------------------------------------*/
void FBI::FBIGeometryCoupler::Setup(std::vector<Teuchos::RCP<DRT::Discretization>>& discretizations,
    Teuchos::RCP<const Epetra_Vector> structure_displacement)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("FBI::FBICoupler::Setup");
  Teuchos::TimeMonitor monitor(*t);

  fluidpositions_ = Teuchos::rcp(new std::map<int, CORE::LINALG::Matrix<3, 1>>);
  beampositions_ = Teuchos::rcp(new std::map<int, CORE::LINALG::Matrix<3, 1>>);

  // todo Specific for fixed grids.. we will have to do something here for ALE (overload?)
  ComputeFixedPositions(*discretizations[1], fluidpositions_);

  // Computes a bounding box for the fluid elements, within which the search will be done
  CORE::LINALG::Matrix<3, 2> fluidBox = CORE::GEO::getXAABBofPositions(*fluidpositions_);

  // Sets-up the searchtree (octtree) for the fluid elements in the given bounding box
  searchtree_->initializeTree(
      fluidBox, *discretizations[1], CORE::GEO::TreeType(CORE::GEO::OCTTREE));
}
/*----------------------------------------------------------------------*/

Teuchos::RCP<std::map<int, std::vector<int>>> FBI::FBIGeometryCoupler::Search(
    std::vector<Teuchos::RCP<DRT::Discretization>>& discretizations,
    Teuchos::RCP<const Epetra_Vector>& column_structure_displacement)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("FBI::FBICoupler::Search");
  Teuchos::TimeMonitor monitor(*t);

  // Vector to hand elements pointers to the bridge object
  Teuchos::RCP<std::map<int, std::vector<int>>> pairids =
      Teuchos::rcp(new std::map<int, std::vector<int>>);

  // todo Specific to 'linearized penalty'. Maybe have to do something for structure+beam in
  // discretization.
  ComputeCurrentPositions(*discretizations[0], beampositions_, column_structure_displacement);

  // todo Maybe have to do something for structure+beam in discretization.
  std::map<int, CORE::LINALG::Matrix<3, 1>>::const_iterator beamnodeiterator;

  // loop over beam nodes
  for (beamnodeiterator = beampositions_->begin(); beamnodeiterator != beampositions_->end();
       beamnodeiterator++)
  {
    const CORE::LINALG::Matrix<3, 1>& curbeamnodeposition = beamnodeiterator->second;

    // search for all fluid elements in the given radius
    std::map<int, std::set<int>> closeeles = searchtree_->searchElementsInRadius(
        *discretizations[1], *fluidpositions_, curbeamnodeposition, searchradius_, 0);

    // loop over the map of beam node-IDs and fluid elements within the search radius
    for (std::map<int, std::set<int>>::const_iterator closefluideles = closeeles.begin();
         closefluideles != closeeles.end(); closefluideles++)
    {
      const DRT::Node* const beamnode = discretizations[0]->gNode(beamnodeiterator->first);
      const DRT::Element* const* beamelements = beamnode->Elements();

      // loop over the set of beam elements adjacent to the current beam node (this leads to
      // duplicate pairs)
      for (int beamelementsnumber = 0; beamelementsnumber < beamnode->NumElement();
           beamelementsnumber++)
      {
        // loop over the gids of the fluid elements
        for (std::set<int>::const_iterator fluideleIter = (closefluideles->second).begin();
             fluideleIter != (closefluideles->second).end(); fluideleIter++)
        {
          if (discretizations[1]->gElement(*fluideleIter)->Owner() ==
              discretizations[0]->Comm().MyPID())
          {
            // store pairs because we have to create them on the beam element owner and we are
            // currently on the fluid element owner
            (*pairids)[beamelements[beamelementsnumber]->Id()].push_back(*fluideleIter);
          }
        }
      }
    }
  }
  return pairids;
}

/*----------------------------------------------------------------------*/

// todo Maybe we can use DRT::UTILS::GhostDiscretizationOnAllProcs instead
// todo Needs to be adapted as soon as problems can contain beam and general structure nodes
void FBI::FBIGeometryCoupler::ExtendBeamGhosting(DRT::Discretization& discretization)
{
  // DRT::UTILS::GhostDiscretizationOnAllProcs(structure_->Discretization());
  std::vector<int> allproc(discretization.Comm().NumProc());
  for (int i = 0; i < discretization.Comm().NumProc(); ++i) allproc[i] = i;

  // fill my own row node ids
  const Epetra_Map* noderowmap = discretization.NodeRowMap();
  std::vector<int> sdata(noderowmap->NumMyElements());
  for (int i = 0; i < noderowmap->NumMyElements(); ++i) sdata[i] = noderowmap->GID(i);

  // gather all gids of nodes redundantly
  std::vector<int> rdata;
  CORE::LINALG::Gather<int>(
      sdata, rdata, (int)allproc.size(), allproc.data(), discretization.Comm());

  // build completely overlapping map of nodes (on ALL processors)
  Teuchos::RCP<Epetra_Map> newnodecolmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, discretization.Comm()));
  sdata.clear();
  rdata.clear();

  // fill my own row element ids
  const Epetra_Map* elerowmap = discretization.ElementRowMap();
  sdata.resize(elerowmap->NumMyElements());
  for (int i = 0; i < elerowmap->NumMyElements(); ++i) sdata[i] = elerowmap->GID(i);

  // gather all gids of elements redundantly
  rdata.resize(0);
  CORE::LINALG::Gather<int>(
      sdata, rdata, (int)allproc.size(), allproc.data(), discretization.Comm());

  // build complete overlapping map of elements (on ALL processors)
  Teuchos::RCP<Epetra_Map> newelecolmap =
      Teuchos::rcp(new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, discretization.Comm()));
  sdata.clear();
  rdata.clear();
  allproc.clear();

  // redistribute the discretization of the interface according to the
  // new column layout
  discretization.ExportColumnNodes(*newnodecolmap);
  discretization.ExportColumnElements(*newelecolmap);

  discretization.FillComplete(true, false, false);
}

/*----------------------------------------------------------------------*/

void FBI::FBIGeometryCoupler::PreparePairCreation(
    std::vector<Teuchos::RCP<DRT::Discretization>>& discretizations,
    Teuchos::RCP<std::map<int, std::vector<int>>> pairids)
{
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("FBI::FBICoupler::PreparePairCreation");
  Teuchos::TimeMonitor monitor(*t);

  std::vector<std::vector<int>> element_senddata(discretizations[0]->Comm().NumProc());
  std::vector<std::vector<int>> node_senddata(discretizations[0]->Comm().NumProc());
  std::vector<std::vector<int>> pairids_to_send(element_senddata.size());
  std::vector<std::vector<int>> pairids_to_recv;
  std::vector<int> element_recvdata;
  std::vector<int> node_recvdata;
  std::vector<int> nodegids;

  for (int i = 0; i < discretizations[0]->Comm().NumProc(); ++i)
  {
    element_senddata[i] = std::vector<int>();
    pairids_to_send[i] = std::vector<int>();
  }

  // Create data containing the ids of beam-fluid element pairs as well as the corresponding fluid
  // element gids which have to be sent to the beam element owner in order to create a pair
  int owner;
  std::map<int, std::vector<int>>::const_iterator beamelementiterator;
  for (beamelementiterator = pairids->begin(); beamelementiterator != pairids->end();
       beamelementiterator++)
  {
    DRT::Element* beamele = discretizations[0]->gElement(beamelementiterator->first);
    if (!beamele) dserror("There is no element with gid %i", beamelementiterator->first);
    owner = beamele->Owner();
    for (std::vector<int>::const_iterator fluideleIter = beamelementiterator->second.begin();
         fluideleIter != (beamelementiterator->second).end(); fluideleIter++)
    {
      element_senddata[owner].push_back(*fluideleIter);

      pairids_to_send[owner].push_back(beamelementiterator->first);

      pairids_to_send[owner].push_back(*fluideleIter);
    }
  }

  // Communicate pair ids
  CORE::LINALG::AllToAllCommunication(discretizations[0]->Comm(), pairids_to_send, pairids_to_recv);

  // bring pair_ids in correct format
  pairids->clear();
  for (int proc = 0; proc < (int)pairids_to_recv.size(); proc++)
  {
    for (int pair = 0; pair < (int)pairids_to_recv[proc].size() - 1; pair = pair + 2)
    {
      (*pairids)[pairids_to_recv[proc][pair]].push_back(pairids_to_recv[proc][pair + 1]);
    }
  }

  // Communicate element gids
  CORE::LINALG::AllToAllCommunication(
      discretizations[0]->Comm(), element_senddata, element_recvdata);


  // Add my current column elements to the set for the map
  const Epetra_Map* elecolmap = discretizations[1]->ElementColMap();
  for (int i = 0; i < elecolmap->NumMyElements(); ++i)
  {
    int gid = elecolmap->GID(i);
    DRT::Element* ele = discretizations[1]->gElement(gid);
    if (!ele) dserror("Cannot find element with gid %", gid);
    element_recvdata.push_back(gid);
  }

  // build overlapping column map of the elements
  Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(new Epetra_Map(
      -1, (int)element_recvdata.size(), element_recvdata.data(), 0, discretizations[1]->Comm()));


  if (!newelecolmap->SameAs(*elecolmap))
  {
    // get node gids of all nodes within the elements to be communicated
    for (int proc = 0; proc < (int)element_senddata.size(); proc++)
    {
      for (unsigned int ele = 0; ele < element_senddata[proc].size(); ele++)
      {
        DRT::Element* element = discretizations[1]->gElement(element_senddata[proc][ele]);
        if (!element) dserror("Cannot find node with gid %", element_senddata[proc][ele]);
        for (int node = 0; node < element->NumNode(); node++)
        {
          node_senddata[proc].push_back((element->NodeIds())[node]);
        }
      }
    }

    // communicate node gids
    CORE::LINALG::AllToAllCommunication(discretizations[0]->Comm(), node_senddata, node_recvdata);

    // add new node gids to overlapping column map
    const Epetra_Map* nodecolmap = discretizations[1]->NodeColMap();
    for (int i = 0; i < nodecolmap->NumMyElements(); ++i)
    {
      int gid = nodecolmap->GID(i);
      DRT::Node* node = discretizations[1]->gNode(gid);
      if (!node) dserror("Cannot find node with gid %", gid);
      node_recvdata.push_back(gid);
    }

    // build complete overlapping map of elements (on ALL processors)
    Teuchos::RCP<Epetra_Map> newnodecolmap = Teuchos::rcp(new Epetra_Map(
        -1, (int)node_recvdata.size(), node_recvdata.data(), 0, discretizations[1]->Comm()));

    // export nodes and elements
    discretizations[1]->ExportColumnNodes(*newnodecolmap);
    discretizations[1]->ExportColumnElements(*newelecolmap);

    Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discretizations[1], true)
        ->FillCompleteFaces(true, true, true, edgebased_fluidstabilization_);
  }

  element_senddata.clear();
  element_recvdata.clear();
  node_senddata.clear();
  node_recvdata.clear();
}

/*----------------------------------------------------------------------*/
void FBI::FBIGeometryCoupler::ComputeFixedPositions(DRT::Discretization& dis,
    Teuchos::RCP<std::map<int, CORE::LINALG::Matrix<3, 1>>> positions) const
{
  positions->clear();
  for (int lid = 0; lid < dis.NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = dis.lColNode(lid);

    for (int d = 0; d < 3; ++d) (*positions)[node->Id()](d) = node->X()[d];
  }
}
/*----------------------------------------------------------------------*/

void FBI::FBIGeometryCoupler::ComputeCurrentPositions(DRT::Discretization& dis,
    Teuchos::RCP<std::map<int, CORE::LINALG::Matrix<3, 1>>> positions,
    Teuchos::RCP<const Epetra_Vector> disp) const
{
  positions->clear();
  std::vector<int> src_dofs(
      9);  // todo this does not work for all possible elements, does it? Variable size?
  std::vector<double> mydisp(3, 0.0);

  for (int lid = 0; lid < dis.NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = dis.lColNode(lid);
    if (disp != Teuchos::null)
    {
      // get the DOF numbers of the current node
      dis.Dof(node, 0, src_dofs);
      // get the current displacements
      DRT::UTILS::ExtractMyValues(*disp, mydisp, src_dofs);

      for (int d = 0; d < 3; ++d) (*positions)[node->Id()](d) = node->X()[d] + mydisp.at(d);
    }
  }
}

/*----------------------------------------------------------------------*/

void FBI::FBIGeometryCoupler::SetBinning(Teuchos::RCP<BINSTRATEGY::BinningStrategy> binning){};

BACI_NAMESPACE_CLOSE
