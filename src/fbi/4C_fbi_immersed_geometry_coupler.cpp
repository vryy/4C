/*----------------------------------------------------------------------*/
/*! \file

\brief Class containing geometric operations usually needed for the coupling of an embedded body.
The current implementation does not scale at all!

\level 3

*----------------------------------------------------------------------*/
#include "4C_fbi_immersed_geometry_coupler.hpp"

#include "4C_binstrategy.hpp"
#include "4C_binstrategy_utils.hpp"
#include "4C_fem_discretization_faces.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_fem_general_node.hpp"
#include "4C_fem_geometry_searchtree.hpp"
#include "4C_fem_geometry_searchtree_service.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_utils_densematrix_communication.hpp"
#include "4C_rebalance_binning_based.hpp"

#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/

FBI::FBIGeometryCoupler::FBIGeometryCoupler()
    : fluidpositions_(new std::map<int, Core::LinAlg::Matrix<3, 1>>),
      beampositions_(new std::map<int, Core::LinAlg::Matrix<3, 1>>),
      searchtree_(new Core::Geo::SearchTree(5)),
      searchradius_(Global::Problem::instance()
                        ->fbi_params()
                        .sublist("BEAM TO FLUID MESHTYING")
                        .get<double>("SEARCH_RADIUS")),
      edgebased_fluidstabilization_(false)
{
  edgebased_fluidstabilization_ = (Core::UTILS::IntegralValue<Inpar::FLUID::StabType>(
                                       Global::Problem::instance()->fluid_dynamic_params().sublist(
                                           "RESIDUAL-BASED STABILIZATION"),
                                       "STABTYPE") == Inpar::FLUID::stabtype_edgebased);
}
/*----------------------------------------------------------------------*/
void FBI::FBIGeometryCoupler::setup(
    std::vector<Teuchos::RCP<Core::FE::Discretization>>& discretizations,
    Teuchos::RCP<const Epetra_Vector> structure_displacement)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("FBI::FBICoupler::Setup");
  Teuchos::TimeMonitor monitor(*t);

  fluidpositions_ = Teuchos::rcp(new std::map<int, Core::LinAlg::Matrix<3, 1>>);
  beampositions_ = Teuchos::rcp(new std::map<int, Core::LinAlg::Matrix<3, 1>>);

  // todo Specific for fixed grids.. we will have to do something here for ALE (overload?)
  compute_fixed_positions(*discretizations[1], fluidpositions_);

  // Computes a bounding box for the fluid elements, within which the search will be done
  Core::LinAlg::Matrix<3, 2> fluidBox = Core::Geo::getXAABBofPositions(*fluidpositions_);

  // Sets-up the searchtree (octtree) for the fluid elements in the given bounding box
  searchtree_->initialize_tree(
      fluidBox, *discretizations[1], Core::Geo::TreeType(Core::Geo::OCTTREE));
}
/*----------------------------------------------------------------------*/

Teuchos::RCP<std::map<int, std::vector<int>>> FBI::FBIGeometryCoupler::search(
    std::vector<Teuchos::RCP<Core::FE::Discretization>>& discretizations,
    Teuchos::RCP<const Epetra_Vector>& column_structure_displacement)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("FBI::FBICoupler::Search");
  Teuchos::TimeMonitor monitor(*t);

  // Vector to hand elements pointers to the bridge object
  Teuchos::RCP<std::map<int, std::vector<int>>> pairids =
      Teuchos::rcp(new std::map<int, std::vector<int>>);

  // todo Specific to 'linearized penalty'. Maybe have to do something for structure+beam in
  // discretization.
  compute_current_positions(*discretizations[0], beampositions_, column_structure_displacement);

  // todo Maybe have to do something for structure+beam in discretization.
  std::map<int, Core::LinAlg::Matrix<3, 1>>::const_iterator beamnodeiterator;

  // loop over beam nodes
  for (beamnodeiterator = beampositions_->begin(); beamnodeiterator != beampositions_->end();
       beamnodeiterator++)
  {
    const Core::LinAlg::Matrix<3, 1>& curbeamnodeposition = beamnodeiterator->second;

    // search for all fluid elements in the given radius
    std::map<int, std::set<int>> closeeles = searchtree_->search_elements_in_radius(
        *discretizations[1], *fluidpositions_, curbeamnodeposition, searchradius_, 0);

    // loop over the map of beam node-IDs and fluid elements within the search radius
    for (std::map<int, std::set<int>>::const_iterator closefluideles = closeeles.begin();
         closefluideles != closeeles.end(); closefluideles++)
    {
      const Core::Nodes::Node* const beamnode = discretizations[0]->g_node(beamnodeiterator->first);
      const Core::Elements::Element* const* beamelements = beamnode->elements();

      // loop over the set of beam elements adjacent to the current beam node (this leads to
      // duplicate pairs)
      for (int beamelementsnumber = 0; beamelementsnumber < beamnode->num_element();
           beamelementsnumber++)
      {
        // loop over the gids of the fluid elements
        for (std::set<int>::const_iterator fluideleIter = (closefluideles->second).begin();
             fluideleIter != (closefluideles->second).end(); fluideleIter++)
        {
          if (discretizations[1]->g_element(*fluideleIter)->owner() ==
              discretizations[0]->get_comm().MyPID())
          {
            // store pairs because we have to create them on the beam element owner and we are
            // currently on the fluid element owner
            (*pairids)[beamelements[beamelementsnumber]->id()].push_back(*fluideleIter);
          }
        }
      }
    }
  }
  return pairids;
}

/*----------------------------------------------------------------------*/

// todo Maybe we can use Core::Rebalance::GhostDiscretizationOnAllProcs instead
// todo Needs to be adapted as soon as problems can contain beam and general structure nodes
void FBI::FBIGeometryCoupler::extend_beam_ghosting(Core::FE::Discretization& discretization)
{
  // Core::Rebalance::GhostDiscretizationOnAllProcs(structure_->discretization());
  std::vector<int> allproc(discretization.get_comm().NumProc());
  for (int i = 0; i < discretization.get_comm().NumProc(); ++i) allproc[i] = i;

  // fill my own row node ids
  const Epetra_Map* noderowmap = discretization.node_row_map();
  std::vector<int> sdata(noderowmap->NumMyElements());
  for (int i = 0; i < noderowmap->NumMyElements(); ++i) sdata[i] = noderowmap->GID(i);

  // gather all gids of nodes redundantly
  std::vector<int> rdata;
  Core::LinAlg::Gather<int>(
      sdata, rdata, (int)allproc.size(), allproc.data(), discretization.get_comm());

  // build completely overlapping map of nodes (on ALL processors)
  Teuchos::RCP<Epetra_Map> newnodecolmap = Teuchos::rcp(
      new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, discretization.get_comm()));
  sdata.clear();
  rdata.clear();

  // fill my own row element ids
  const Epetra_Map* elerowmap = discretization.element_row_map();
  sdata.resize(elerowmap->NumMyElements());
  for (int i = 0; i < elerowmap->NumMyElements(); ++i) sdata[i] = elerowmap->GID(i);

  // gather all gids of elements redundantly
  rdata.resize(0);
  Core::LinAlg::Gather<int>(
      sdata, rdata, (int)allproc.size(), allproc.data(), discretization.get_comm());

  // build complete overlapping map of elements (on ALL processors)
  Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(
      new Epetra_Map(-1, (int)rdata.size(), rdata.data(), 0, discretization.get_comm()));
  sdata.clear();
  rdata.clear();
  allproc.clear();

  // redistribute the discretization of the interface according to the
  // new column layout
  discretization.export_column_nodes(*newnodecolmap);
  discretization.export_column_elements(*newelecolmap);

  discretization.fill_complete(true, false, false);
}

/*----------------------------------------------------------------------*/

void FBI::FBIGeometryCoupler::prepare_pair_creation(
    std::vector<Teuchos::RCP<Core::FE::Discretization>>& discretizations,
    Teuchos::RCP<std::map<int, std::vector<int>>> pairids)
{
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("FBI::FBICoupler::PreparePairCreation");
  Teuchos::TimeMonitor monitor(*t);

  std::vector<std::vector<int>> element_senddata(discretizations[0]->get_comm().NumProc());
  std::vector<std::vector<int>> node_senddata(discretizations[0]->get_comm().NumProc());
  std::vector<std::vector<int>> pairids_to_send(element_senddata.size());
  std::vector<std::vector<int>> pairids_to_recv;
  std::vector<int> element_recvdata;
  std::vector<int> node_recvdata;
  std::vector<int> nodegids;

  for (int i = 0; i < discretizations[0]->get_comm().NumProc(); ++i)
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
    Core::Elements::Element* beamele = discretizations[0]->g_element(beamelementiterator->first);
    if (!beamele) FOUR_C_THROW("There is no element with gid %i", beamelementiterator->first);
    owner = beamele->owner();
    for (std::vector<int>::const_iterator fluideleIter = beamelementiterator->second.begin();
         fluideleIter != (beamelementiterator->second).end(); fluideleIter++)
    {
      element_senddata[owner].push_back(*fluideleIter);

      pairids_to_send[owner].push_back(beamelementiterator->first);

      pairids_to_send[owner].push_back(*fluideleIter);
    }
  }

  // Communicate pair ids
  Core::LinAlg::AllToAllCommunication(
      discretizations[0]->get_comm(), pairids_to_send, pairids_to_recv);

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
  Core::LinAlg::AllToAllCommunication(
      discretizations[0]->get_comm(), element_senddata, element_recvdata);


  // Add my current column elements to the set for the map
  const Epetra_Map* elecolmap = discretizations[1]->element_col_map();
  for (int i = 0; i < elecolmap->NumMyElements(); ++i)
  {
    int gid = elecolmap->GID(i);
    Core::Elements::Element* ele = discretizations[1]->g_element(gid);
    if (!ele) FOUR_C_THROW("Cannot find element with gid %", gid);
    element_recvdata.push_back(gid);
  }

  // build overlapping column map of the elements
  Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(new Epetra_Map(-1,
      (int)element_recvdata.size(), element_recvdata.data(), 0, discretizations[1]->get_comm()));


  if (!newelecolmap->SameAs(*elecolmap))
  {
    // get node gids of all nodes within the elements to be communicated
    for (int proc = 0; proc < (int)element_senddata.size(); proc++)
    {
      for (unsigned int ele = 0; ele < element_senddata[proc].size(); ele++)
      {
        Core::Elements::Element* element =
            discretizations[1]->g_element(element_senddata[proc][ele]);
        if (!element) FOUR_C_THROW("Cannot find node with gid %", element_senddata[proc][ele]);
        for (int node = 0; node < element->num_node(); node++)
        {
          node_senddata[proc].push_back((element->node_ids())[node]);
        }
      }
    }

    // communicate node gids
    Core::LinAlg::AllToAllCommunication(
        discretizations[0]->get_comm(), node_senddata, node_recvdata);

    // add new node gids to overlapping column map
    const Epetra_Map* nodecolmap = discretizations[1]->node_col_map();
    for (int i = 0; i < nodecolmap->NumMyElements(); ++i)
    {
      int gid = nodecolmap->GID(i);
      Core::Nodes::Node* node = discretizations[1]->g_node(gid);
      if (!node) FOUR_C_THROW("Cannot find node with gid %", gid);
      node_recvdata.push_back(gid);
    }

    // build complete overlapping map of elements (on ALL processors)
    Teuchos::RCP<Epetra_Map> newnodecolmap = Teuchos::rcp(new Epetra_Map(
        -1, (int)node_recvdata.size(), node_recvdata.data(), 0, discretizations[1]->get_comm()));

    // export nodes and elements
    discretizations[1]->export_column_nodes(*newnodecolmap);
    discretizations[1]->export_column_elements(*newelecolmap);

    Teuchos::rcp_dynamic_cast<Core::FE::DiscretizationFaces>(discretizations[1], true)
        ->fill_complete_faces(true, true, true, edgebased_fluidstabilization_);
  }

  element_senddata.clear();
  element_recvdata.clear();
  node_senddata.clear();
  node_recvdata.clear();
}

/*----------------------------------------------------------------------*/
void FBI::FBIGeometryCoupler::compute_fixed_positions(Core::FE::Discretization& dis,
    Teuchos::RCP<std::map<int, Core::LinAlg::Matrix<3, 1>>> positions) const
{
  positions->clear();
  for (int lid = 0; lid < dis.num_my_col_nodes(); ++lid)
  {
    const Core::Nodes::Node* node = dis.l_col_node(lid);

    for (int d = 0; d < 3; ++d) (*positions)[node->id()](d) = node->x()[d];
  }
}
/*----------------------------------------------------------------------*/

void FBI::FBIGeometryCoupler::compute_current_positions(Core::FE::Discretization& dis,
    Teuchos::RCP<std::map<int, Core::LinAlg::Matrix<3, 1>>> positions,
    Teuchos::RCP<const Epetra_Vector> disp) const
{
  positions->clear();
  std::vector<int> src_dofs(
      9);  // todo this does not work for all possible elements, does it? Variable size?
  std::vector<double> mydisp(3, 0.0);

  for (int lid = 0; lid < dis.num_my_col_nodes(); ++lid)
  {
    const Core::Nodes::Node* node = dis.l_col_node(lid);
    if (disp != Teuchos::null)
    {
      // get the DOF numbers of the current node
      dis.dof(node, 0, src_dofs);
      // get the current displacements
      Core::FE::ExtractMyValues(*disp, mydisp, src_dofs);

      for (int d = 0; d < 3; ++d) (*positions)[node->id()](d) = node->x()[d] + mydisp.at(d);
    }
  }
}

/*----------------------------------------------------------------------*/

void FBI::FBIGeometryCoupler::set_binning(
    Teuchos::RCP<Core::Binstrategy::BinningStrategy> binning){};

FOUR_C_NAMESPACE_CLOSE
