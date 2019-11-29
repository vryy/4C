/*----------------------------------------------------------------------*/
/*! \file
\file constraintenforcer_fbi.cpp

\brief Adapter used to transfer beam and fluid data between meshes.

\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/

#include "constraintenforcer_fbi.H"
#include "ad_fbi_constraintbridge.H"

#include "../drt_geometry_pair/geometry_pair.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_geometry/searchtree.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/searchtree.H"
#include "../drt_adapter/ad_str_fbiwrapper.H"
#include "../drt_adapter/ad_fld_fluidbeam_immersed.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_beaminteraction/beaminteraction_calc_utils.H"  // todo put this into bridge to keep everything beam specific in there

ADAPTER::FBIConstraintenforcer::FBIConstraintenforcer(
    Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge)
    : fluid_(Teuchos::null),
      structure_(Teuchos::null),
      searchtree_(new GEO::SearchTree(5)),
      fluidpositions_(new std::map<int, LINALG::Matrix<3, 1>>),
      beampositions_(new std::map<int, LINALG::Matrix<3, 1>>),
      discretizations_(new std::vector<Teuchos::RCP<DRT::Discretization>>),
      bridge_(bridge)
{
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
}

void ADAPTER::FBIConstraintenforcer::Evaluate()
{
  std::cout << "Fluid velocity is "
            << *(Teuchos::rcp_dynamic_cast<ADAPTER::FluidBeamImmersed>(fluid_, true)->Velnp())
            << std::endl;
  std::cout << "Structure velocity is " << *(structure_->Velnp()) << std::endl;
  // Vector to hand elements pointers to the bridge object
  std::vector<DRT::Element const*> ele_ptrs(2);
  Teuchos::RCP<std::vector<double>> beam_dofvec = Teuchos::rcp(new std::vector<double>);
  Teuchos::RCP<std::vector<double>> fluid_dofvec = Teuchos::rcp(new std::vector<double>);
  // todo search radius is hard-coded for the first few tries
  const double searchradius = 0.5;
  bridge_->ResetBridge();
  // todo Specific to 'linearized penalty'. Maybe have to do something for structure+beam in
  // discretization.
  ComputeCurrentPositions(structure_->Discretization(), beampositions_, structure_->Dispn());
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
        ele_ptrs[0] = beamelements[beamelementsnumber];
        // loop over the gids of the fluid elements
        for (std::set<int>::const_iterator fluideleIter = (closefluideles->second).begin();
             fluideleIter != (closefluideles->second).end(); fluideleIter++)
        {
          DRT::Element* fluidele = (fluid_->Discretization())->gElement(*fluideleIter);
          ele_ptrs[1] = fluidele;
          ExtractCurrentElementDofs(ele_ptrs, beam_dofvec, fluid_dofvec);
          bridge_->CreatePair(ele_ptrs, beam_dofvec, fluid_dofvec);
        }
      }
    }
  }
  bridge_->Evaluate(discretizations_);  // Do projection & Evaluate local matrices (needs pairs) &
                                        // Assemble into global matrices (needs local matrices)
  return;
}

// todo check if this is even possible for Epetra_vectors/sparse matrices

Teuchos::RCP<Epetra_Vector> ADAPTER::FBIConstraintenforcer::MasterToSlave()
{
  //  Teuchos::rcp_dynamic_cast<ADAPTER::FluidBeamImmersed>(fluid_, true)
  //      ->SetCouplingContributions(AssembleMasterStiffness());
  //  fluid_->ApplyInterfaceValues(AssembleMasterForce());
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
      9);  // todo this does not work for all possible elements, does it? Variable size?
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
  std::vector<int> tmp_dofs(9);
  std::vector<int> src_dofs(3);
  std::vector<double> myvel(3, 0.0);
  std::vector<double> vel_tmp;

  // extract the current position of the beam element from the displacement vector
  BEAMINTERACTION::UTILS::ExtractPosDofVecAbsoluteValues(*(structure_->Discretization()),
      elements[0], structure_->Dispnp(),
      *beam_dofvec);  // todo check time step todo get "interface" displacements only for beam
                      // elements

  // extract veclocity of the beam element
  BEAMINTERACTION::UTILS::GetCurrentElementDis(
      *(structure_->Discretization()), elements[0], structure_->Velnp(), vel_tmp);

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

  BEAMINTERACTION::UTILS::GetCurrentElementDis(*(fluid_->Discretization()), elements[1],
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidBeamImmersed>(fluid_, true)->Velnp(), vel_tmp);

  for (unsigned int i = 0; i < vel_tmp.size(); i++)
  {
    if ((i + 1) % 4) fluid_dofvec->push_back(vel_tmp[i]);
    //    if (val > 05) printf("fluid velocity is greater than 0.5\n");
  }
}
