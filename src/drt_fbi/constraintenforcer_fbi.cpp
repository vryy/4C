/*----------------------------------------------------------------------*/
/*! \file
\file constraintenforcer_fbi.cpp

\brief Adapter used to transfer beam and fluid data between meshes.

\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/

#include "constraintenforcer_fbi.H"
#include "ad_fbi_constraintbridge.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_geometry_pair/geometry_pair.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_parallel.cpp"
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/searchtree.H"
#include "../drt_adapter/ad_str_fbiwrapper.H"
#include "../drt_adapter/ad_fld_fluidbeam_immersed.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_beaminteraction/beaminteraction_calc_utils.H"  // todo put this into bridge to keep everything beam specific in there
#include "../drt_inpar/inpar_fbi.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_lib/drt_discret_faces.H"

#include <mpi.h>

ADAPTER::FBIConstraintenforcer::FBIConstraintenforcer(
    Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge)
    : fluid_(Teuchos::null),
      structure_(Teuchos::null),
      searchtree_(new GEO::SearchTree(5)),
      fluidpositions_(new std::map<int, LINALG::Matrix<3, 1>>),
      beampositions_(new std::map<int, LINALG::Matrix<3, 1>>),
      discretizations_(new std::vector<Teuchos::RCP<DRT::Discretization>>),
      bridge_(bridge),
      edgebased_fluidstabilization(false),
      column_structure_displacement_(Teuchos::null),
      column_structure_velocity_(Teuchos::null),
      column_fluid_velocity_(Teuchos::null)
{
  edgebased_fluidstabilization =
      (DRT::INPUT::IntegralValue<INPAR::FLUID::StabType>(
           DRT::Problem::Instance()->FluidDynamicParams().sublist("RESIDUAL-BASED STABILIZATION"),
           "STABTYPE") == INPAR::FLUID::stabtype_edgebased);
}

void ADAPTER::FBIConstraintenforcer::Setup(Teuchos::RCP<ADAPTER::FSIStructureWrapper> structure,
    Teuchos::RCP<ADAPTER::FluidMovingBoundary> fluid)
{
  fluid_ = fluid;
  structure_ = structure;
  printf("After structure\n");
  discretizations_->push_back(structure_->Discretization());
  discretizations_->push_back(fluid_->Discretization());
  printf("After discretization\n");
  bridge_->Setup(structure_->Discretization()->DofRowMap(), fluid_->Discretization()->DofRowMap());
  fluidpositions_ = Teuchos::rcp(new std::map<int, LINALG::Matrix<3, 1>>);
  beampositions_ = Teuchos::rcp(new std::map<int, LINALG::Matrix<3, 1>>);
  // todo Specific for fixed grids.. we will have to do something here for ALE (overload?)
  ComputeFixedPositions(fluid_->Discretization(), fluidpositions_);
  // Computes a bounding box for the fluid elements, within which the search will be done
  LINALG::Matrix<3, 2> fluidBox = GEO::getXAABBofPositions(*fluidpositions_);
  // Sets-up the searchtree (octtree) for the fluid elements in the given bounding box
  searchtree_->initializeTree(fluidBox, *(fluid_->Discretization()), GEO::TreeType(GEO::OCTTREE));
  if (structure_->Discretization()->Comm().NumProc() > 1) ExtendBeamGhosting();
}

void ADAPTER::FBIConstraintenforcer::Evaluate()
{
  std::cout << "Fluid velocity is "
            << *(Teuchos::rcp_dynamic_cast<ADAPTER::FluidBeamImmersed>(fluid_, true)->Velnp())
            << std::endl;

  column_structure_displacement_ =
      DRT::UTILS::GetColVersionOfRowVector(structure_->Discretization(), structure_->Dispnp());
  column_structure_velocity_ =
      DRT::UTILS::GetColVersionOfRowVector(structure_->Discretization(), structure_->Velnp());
  column_fluid_velocity_ = DRT::UTILS::GetColVersionOfRowVector(fluid_->Discretization(),
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidBeamImmersed>(fluid_, true)->Velnp());

  std::cout << "Structure velocity is " << *(structure_->Velnp()) << std::endl;
  // Vector to hand elements pointers to the bridge object
  std::map<int, std::vector<int>> pairids;
  // todo search radius is hard-coded for the first few tries
  const double searchradius = 0.4;
  bridge_->ResetBridge();
  // todo Specific to 'linearized penalty'. Maybe have to do something for structure+beam in
  // discretization.
  ComputeCurrentPositions(
      structure_->Discretization(), beampositions_, column_structure_displacement_);
  // todo Maybe have to do something for structure+beam in discretization.
  std::map<int, LINALG::Matrix<3, 1>>::const_iterator beamnodeiterator;
  for (beamnodeiterator = beampositions_->begin(); beamnodeiterator != beampositions_->end();
       beamnodeiterator++)
  {
    const LINALG::Matrix<3, 1>& curbeamnodeposition = beamnodeiterator->second;
    std::map<int, std::set<int>> closeeles =
        searchtree_->searchElementsInRadius(*(fluid_->Discretization()), *fluidpositions_,
            curbeamnodeposition, searchradius, 0);  // todo check label

    // loop over the map of beam node-IDs and fluid elements within the search radius
    for (std::map<int, std::set<int>>::const_iterator closefluideles = closeeles.begin();
         closefluideles != closeeles.end(); closefluideles++)
    {
      const DRT::Node* const beamnode =
          (structure_->Discretization())->gNode(beamnodeiterator->first);
      const DRT::Element* const* beamelements = beamnode->Elements();

      // loop over the set of beam elements adjacent to the current beam node
      for (int beamelementsnumber = 0; beamelementsnumber < beamnode->NumElement();
           beamelementsnumber++)
      {
        // loop over the gids of the fluid elements
        for (std::set<int>::const_iterator fluideleIter = (closefluideles->second).begin();
             fluideleIter != (closefluideles->second).end(); fluideleIter++)
        {
          if (fluid_->Discretization()->gElement(*fluideleIter)->Owner() ==
              structure_->Discretization()->Comm().MyPID())
          {
            (pairids)[beamelements[beamelementsnumber]->Id()].push_back(*fluideleIter);
          }
        }
      }
    }
  }
  CreatePairs(pairids);
  bridge_->Evaluate(discretizations_);
  return;
}

void ADAPTER::FBIConstraintenforcer::CreatePairs(std::map<int, std::vector<int>> pairids)
{
  if ((structure_->Discretization())->Comm().NumProc() > 1)
  {
    PreparePairCreation(pairids);
    column_structure_displacement_ =
        DRT::UTILS::GetColVersionOfRowVector(structure_->Discretization(), structure_->Dispnp());
    column_structure_velocity_ =
        DRT::UTILS::GetColVersionOfRowVector(structure_->Discretization(), structure_->Velnp());
    column_fluid_velocity_ = DRT::UTILS::GetColVersionOfRowVector(fluid_->Discretization(),
        Teuchos::rcp_dynamic_cast<ADAPTER::FluidBeamImmersed>(fluid_, true)->Velnp());
  }

  std::vector<DRT::Element const*> ele_ptrs(2);
  Teuchos::RCP<std::vector<double>> beam_dofvec = Teuchos::rcp(new std::vector<double>);
  Teuchos::RCP<std::vector<double>> fluid_dofvec = Teuchos::rcp(new std::vector<double>);

  std::map<int, std::vector<int>>::const_iterator beamelementiterator;
  for (beamelementiterator = pairids.begin(); beamelementiterator != pairids.end();
       beamelementiterator++)
  {
    ele_ptrs[0] = (structure_->Discretization())->gElement(beamelementiterator->first);
    if (ele_ptrs[0]->Owner() != structure_->Discretization()->Comm().MyPID())
      dserror(
          "For now we can only create the pair on the beam owner, but beam element owner is %i and "
          "we are on proc %i \n",
          ele_ptrs[0]->Owner(), structure_->Discretization()->Comm().MyPID());

    for (std::vector<int>::const_iterator fluideleIter = beamelementiterator->second.begin();
         fluideleIter != (beamelementiterator->second).end(); fluideleIter++)
    {
      DRT::Element* fluidele = (fluid_->Discretization())->gElement(*fluideleIter);
      ele_ptrs[1] = fluidele;
      ExtractCurrentElementDofs(ele_ptrs, beam_dofvec, fluid_dofvec);
      bridge_->CreatePair(ele_ptrs, beam_dofvec, fluid_dofvec);
    }
  }
}

void ADAPTER::FBIConstraintenforcer::PreparePairCreation(std::map<int, std::vector<int>>& pairids)
{
  //  std::vector<int> allproc(structure_->Discretization()->Comm().NumProc());
  std::vector<std::vector<int>> element_senddata(structure_->Discretization()->Comm().NumProc());
  std::vector<std::vector<int>> node_senddata(structure_->Discretization()->Comm().NumProc());
  std::vector<std::vector<int>> pairids_to_send(element_senddata.size());
  std::vector<std::vector<int>> pairids_to_recv;
  std::vector<int> element_recvdata;
  std::vector<int> node_recvdata;
  std::vector<int> nodegids;

  for (int i = 0; i < structure_->Discretization()->Comm().NumProc(); ++i)
  {
    element_senddata[i] = std::vector<int>();
    pairids_to_send[i] = std::vector<int>();
  }

  // Create data containing the ids of beam-fluid element pairs as well as the corresponding fluid
  // element gids which have to be sent to the beam element owner in order to create a pair
  int owner;
  std::map<int, std::vector<int>>::const_iterator beamelementiterator;
  for (beamelementiterator = pairids.begin(); beamelementiterator != pairids.end();
       beamelementiterator++)
  {
    DRT::Element* beamele = structure_->Discretization()->gElement(beamelementiterator->first);
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
  LINALG::AllToAllCommunication(
      structure_->Discretization()->Comm(), pairids_to_send, pairids_to_recv);

  pairids.clear();
  for (int proc = 0; proc < (int)pairids_to_recv.size(); proc++)
  {
    for (int pair = 0; pair < (int)pairids_to_recv[proc].size() - 1; pair = pair + 2)
    {
      pairids[pairids_to_recv[proc][pair]].push_back(pairids_to_recv[proc][pair + 1]);
    }
  }

  // Communicate element gids
  LINALG::AllToAllCommunication(
      structure_->Discretization()->Comm(), element_senddata, element_recvdata);


  // Add my current column elements to the set for the map
  const Epetra_Map* elecolmap = fluid_->Discretization()->ElementColMap();
  for (int i = 0; i < elecolmap->NumMyElements(); ++i)
  {
    int gid = elecolmap->GID(i);
    DRT::Element* ele = fluid_->Discretization()->gElement(gid);
    if (!ele) dserror("ERROR: Cannot find element with gid %", gid);
    element_recvdata.push_back(gid);
  }

  // build overlapping column map of the elements
  Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(new Epetra_Map(
      -1, (int)element_recvdata.size(), &element_recvdata[0], 0, fluid_->Discretization()->Comm()));


  for (int proc = 0; proc < (int)element_senddata.size(); proc++)
  {
    for (unsigned int ele = 0; ele < element_senddata[proc].size(); ele++)
    {
      DRT::Element* element = fluid_->Discretization()->gElement(element_senddata[proc][ele]);
      if (!element) dserror("ERROR: Cannot find node with gid %", element_senddata[proc][ele]);
      for (int node = 0; node < element->NumNode(); node++)
      {
        node_senddata[proc].push_back((element->NodeIds())[node]);
      }
    }
  }

  LINALG::AllToAllCommunication(structure_->Discretization()->Comm(), node_senddata, node_recvdata);

  const Epetra_Map* nodecolmap = fluid_->Discretization()->NodeColMap();
  for (int i = 0; i < nodecolmap->NumMyElements(); ++i)
  {
    int gid = nodecolmap->GID(i);
    DRT::Node* node = fluid_->Discretization()->gNode(gid);
    if (!node) dserror("ERROR: Cannot find node with gid %", gid);
    node_recvdata.push_back(gid);
  }

  // build complete overlapping map of elements (on ALL processors)
  Teuchos::RCP<Epetra_Map> newnodecolmap = Teuchos::rcp(new Epetra_Map(
      -1, (int)node_recvdata.size(), &node_recvdata[0], 0, fluid_->Discretization()->Comm()));

  fluid_->Discretization()->ExportColumnNodes(*newnodecolmap);
  fluid_->Discretization()->ExportColumnElements(*newelecolmap);

  // We need this cast in case we want to use edge-based stabilization.. Maybe only to it in that
  // particular case?
  Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(fluid_->Discretization(), true)
      ->FillCompleteFaces(true, true, true, true);
  element_senddata.clear();
  element_recvdata.clear();
  node_senddata.clear();
  node_recvdata.clear();
}

// todo Maybe we can use DRT::UTILS::GhostDiscretizationOnAllProcs instead
// todo Needs to be adapted as soon as problems can contain beam and general structure nodes
void ADAPTER::FBIConstraintenforcer::ExtendBeamGhosting()
{
  // DRT::UTILS::GhostDiscretizationOnAllProcs(structure_->Discretization());
  std::vector<int> allproc(structure_->Discretization()->Comm().NumProc());
  for (int i = 0; i < structure_->Discretization()->Comm().NumProc(); ++i) allproc[i] = i;

  // fill my own row node ids
  const Epetra_Map* noderowmap = structure_->Discretization()->NodeRowMap();
  std::vector<int> sdata(noderowmap->NumMyElements());
  for (int i = 0; i < noderowmap->NumMyElements(); ++i) sdata[i] = noderowmap->GID(i);

  // gather all gids of nodes redundantly
  std::vector<int> rdata;
  LINALG::Gather<int>(
      sdata, rdata, (int)allproc.size(), &allproc[0], structure_->Discretization()->Comm());

  // build completely overlapping map of nodes (on ALL processors)
  Teuchos::RCP<Epetra_Map> newnodecolmap = Teuchos::rcp(
      new Epetra_Map(-1, (int)rdata.size(), &rdata[0], 0, structure_->Discretization()->Comm()));
  sdata.clear();
  rdata.clear();

  // fill my own row element ids
  const Epetra_Map* elerowmap = structure_->Discretization()->ElementRowMap();
  sdata.resize(elerowmap->NumMyElements());
  for (int i = 0; i < elerowmap->NumMyElements(); ++i) sdata[i] = elerowmap->GID(i);

  // gather all gids of elements redundantly
  rdata.resize(0);
  LINALG::Gather<int>(
      sdata, rdata, (int)allproc.size(), &allproc[0], structure_->Discretization()->Comm());

  // build complete overlapping map of elements (on ALL processors)
  Teuchos::RCP<Epetra_Map> newelecolmap = Teuchos::rcp(
      new Epetra_Map(-1, (int)rdata.size(), &rdata[0], 0, structure_->Discretization()->Comm()));
  sdata.clear();
  rdata.clear();
  allproc.clear();

  // redistribute the discretization of the interface according to the
  // new column layout
  structure_->Discretization()->ExportColumnNodes(*newnodecolmap);
  structure_->Discretization()->ExportColumnElements(*newelecolmap);

  structure_->Discretization()->FillComplete(true, false, false);

  Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_, true)
      ->SetupMultiMapExtractor();
}


Teuchos::RCP<Epetra_Vector> ADAPTER::FBIConstraintenforcer::MasterToSlave()
{
  const Teuchos::ParameterList& fbi = DRT::Problem::Instance()->FBIParams();
  if (DRT::INPUT::IntegralValue<int>(fbi, "COUPLING") != INPAR::FBI::BeamToFluidCoupling::solid)
  {
    Teuchos::rcp_dynamic_cast<ADAPTER::FluidBeamImmersed>(fluid_, true)
        ->SetCouplingContributions(AssembleMasterStiffness());
    fluid_->ApplyInterfaceValues(AssembleMasterForce());
  }

  return Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_, true)
      ->ExtractInterfaceVelnp();
};


// This is the interface to the outside world regarding the force.
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIConstraintenforcer::SlaveToMaster()
{
  return AssembleSlaveForce();
};

void ADAPTER::FBIConstraintenforcer::ComputeFixedPositions(Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<std::map<int, LINALG::Matrix<3, 1>>> positions) const
{
  positions->clear();
  for (int lid = 0; lid < dis->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = dis->lColNode(lid);

    for (int d = 0; d < 3; ++d) (*positions)[node->Id()](d) = node->X()[d];
  }
}

void ADAPTER::FBIConstraintenforcer::ComputeCurrentPositions(Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<std::map<int, LINALG::Matrix<3, 1>>> positions,
    Teuchos::RCP<const Epetra_Vector> disp) const
{
  positions->clear();
  std::vector<int> src_dofs(
      6);  // todo this does not work for all possible elements, does it? Variable size?
  std::vector<double> mydisp(3, 0.0);
  for (int lid = 0; lid < dis->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = dis->lColNode(lid);
    if (disp != Teuchos::null)
    {
      // get the DOF numbers of the current node
      dis->Dof(node, 0, src_dofs);
      // get the current displacements
      DRT::UTILS::ExtractMyValues(*disp, mydisp, src_dofs);

      for (int d = 0; d < 3; ++d) (*positions)[node->Id()](d) = node->X()[d] + mydisp.at(d);
    }
  }
}

// todo take a good look at dispn vs. dispnp and veln vs. velnp!
void ADAPTER::FBIConstraintenforcer::ExtractCurrentElementDofs(
    std::vector<DRT::Element const*> elements, Teuchos::RCP<std::vector<double>>& beam_dofvec,
    Teuchos::RCP<std::vector<double>>& fluid_dofvec) const
{
  std::vector<int> tmp_dofs(6);
  std::vector<int> src_dofs(3);
  std::vector<double> myvel(3, 0.0);
  std::vector<double> vel_tmp;

  // extract the current position of the beam element from the displacement vector
  BEAMINTERACTION::UTILS::ExtractPosDofVecAbsoluteValues(*(structure_->Discretization()),
      elements[0], column_structure_displacement_,
      *beam_dofvec);  // todo check time step todo get "interface" displacements only for beam
                      // elements
  ;
  // extract veclocity of the beam element
  BEAMINTERACTION::UTILS::GetCurrentElementDis(
      *(structure_->Discretization()), elements[0], column_structure_velocity_, vel_tmp);

  for (double val : vel_tmp) beam_dofvec->push_back(val);

  vel_tmp.clear();
  // extract the current positions and velocities of the fluid element todo only valid for fixed
  // grid, not for ALE
  fluid_dofvec->clear();
  const DRT::Node* const* fluidnodes = elements[1]->Nodes();
  for (int lid = 0; lid < elements[1]->NumNode(); ++lid)
  {
    for (int dim = 0; dim < 3; dim++)
    {
      fluid_dofvec->push_back(fluidnodes[lid]->X()[dim]);
      //      printf("position of fluid element %i is %f\n", lid, fluidnodes[lid]->X()[d]);
    }
  }

  BEAMINTERACTION::UTILS::GetCurrentElementDis(
      *(fluid_->Discretization()), elements[1], column_fluid_velocity_, vel_tmp);

  for (unsigned int i = 0; i < vel_tmp.size(); i++)
  {
    if ((i + 1) % 4) fluid_dofvec->push_back(vel_tmp[i]);
    //    if (val > 05) printf("fluid velocity is greater than 0.5\n");
  }
}
