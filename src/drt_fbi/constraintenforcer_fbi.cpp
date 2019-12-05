/*----------------------------------------------------------------------*/
/*! \file

\brief Abstract class to be overloaded by different constraint enforcement techniques for fluid-beam
interaction.

\level 2

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/

#include "constraintenforcer_fbi.H"

#include "ad_fbi_constraintbridge.H"
#include "ad_fbi_constraintbridge_penalty.H"
#include "immersed_geometry_coupler_fbi.H"

#include "../drt_adapter/ad_fld_fbi_movingboundary.H"
#include "../drt_adapter/ad_str_fbiwrapper.H"
#include "../drt_beaminteraction/beaminteraction_calc_utils.H"  // todo put this into bridge to keep everything beam specific in there
#include "../drt_geometry_pair/geometry_pair.H"
#include "../drt_inpar/inpar_fbi.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_node.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_utils_parallel.H"
#include "../drt_lib/drt_discret_faces.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"

#include <iostream>

ADAPTER::FBIConstraintenforcer::FBIConstraintenforcer(
    Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge,
    Teuchos::RCP<FBI::FBIGeometryCoupler> geometrycoupler)
    : fluid_(Teuchos::null),
      structure_(Teuchos::null),
      discretizations_(new std::vector<Teuchos::RCP<DRT::Discretization>>),
      bridge_(bridge),
      geometrycoupler_(geometrycoupler),
      column_structure_displacement_(Teuchos::null),
      column_structure_velocity_(Teuchos::null),
      column_fluid_velocity_(Teuchos::null)
{
}

/*----------------------------------------------------------------------*/

void ADAPTER::FBIConstraintenforcer::Setup(Teuchos::RCP<ADAPTER::FSIStructureWrapper> structure,
    Teuchos::RCP<ADAPTER::FluidMovingBoundary> fluid)
{
  fluid_ = fluid;
  structure_ = structure;
  discretizations_->push_back(structure_->Discretization());
  discretizations_->push_back(fluid_->Discretization());

  bridge_->Setup(structure_->Discretization()->DofRowMap(), fluid_->Discretization()->DofRowMap());
  geometrycoupler_->Setup(discretizations_);

  if (structure_->Discretization()->Comm().NumProc() > 1)
  {
    geometrycoupler_->ExtendBeamGhosting(structure->Discretization());

    // After ghosting we need to explicitly set up the MultiMapExtractor again
    Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_, true)
        ->SetupMultiMapExtractor();
  }
  std::ofstream log;
  if ((*discretizations_)[1]->Comm().MyPID() == 0)
  {
    std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
    s.append(".penalty");
    log.open(s.c_str(), std::ofstream::out);
    log << "Time Step ViolationNorm FluidViolationNorm StructureViolationNorm" << std::endl;
    log.close();
  }
}

/*----------------------------------------------------------------------*/

void ADAPTER::FBIConstraintenforcer::Evaluate()
{
  // We use the column vectors here, because currently the search is based on neighboring nodes,
  // but the element pairs are created using the elements needing all information on all their
  // DOFs
  column_structure_displacement_ =
      DRT::UTILS::GetColVersionOfRowVector(structure_->Discretization(), structure_->Dispnp());
  column_structure_velocity_ =
      DRT::UTILS::GetColVersionOfRowVector(structure_->Discretization(), structure_->Velnp());
  column_fluid_velocity_ = DRT::UTILS::GetColVersionOfRowVector(fluid_->Discretization(),
      Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)->Velnp());

  // Before each search we delete all pair and segment information
  bridge_->ResetBridge();

  // Do the search in the geometrycoupler_ and return the possible pair ids
  Teuchos::RCP<std::map<int, std::vector<int>>> pairids = geometrycoupler_->Search(discretizations_,
      column_structure_displacement_, column_structure_velocity_,
      column_fluid_velocity_);  // todo make this a vector? At some point we probably need the ale
                                // displacements as well

  // For now we need to separate the pair creation from the search, since the search takes place
  // on the fluid elements owner, while (for now) the pair has to be created on the beam element
  // owner
  CreatePairs(pairids);

  // Create all needed matrix and vector contributions based on the current state
  bridge_->Evaluate(discretizations_);
}

/*----------------------------------------------------------------------*/

Teuchos::RCP<Epetra_Vector> ADAPTER::FBIConstraintenforcer::StructureToFluid()
{
  // todo only access the parameter list once

  // Check if we want to couple the fluid
  const Teuchos::ParameterList& fbi = DRT::Problem::Instance()->FBIParams();
  if (DRT::INPUT::IntegralValue<int>(fbi, "COUPLING") != INPAR::FBI::BeamToFluidCoupling::solid)
  {
    // Assemble the fluid stiffness matrix and hand it to the fluid solver
    Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)
        ->SetCouplingContributions(AssembleFluidStiffness());

    // Assemble the fluid force vector and hand it to the fluid solver
    fluid_->ApplyInterfaceValues(AssembleFluidForce());
  }

  // return the current struture velocity
  return Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_, true)
      ->ExtractInterfaceVelnp();
};

/*----------------------------------------------------------------------*/

// return the structure force
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIConstraintenforcer::FluidToStructure()
{
  return AssembleStructureForce();
};

/*----------------------------------------------------------------------*/

// For now we need to separate the pair creation from the search, since the search takes place on
// the fluid elements owner, while (for now) the pair has to be created on the beam element owner
void ADAPTER::FBIConstraintenforcer::CreatePairs(
    Teuchos::RCP<std::map<int, std::vector<int>>> pairids)
{
  if ((structure_->Discretization())->Comm().NumProc() > 1)
  {
    // The geometrycoupler takes care of all MPI communication that needs to be done before the
    // pairs can finally be created
    geometrycoupler_->PreparePairCreation(discretizations_, pairids);

    column_structure_displacement_ =
        DRT::UTILS::GetColVersionOfRowVector(structure_->Discretization(), structure_->Dispnp());
    column_structure_velocity_ =
        DRT::UTILS::GetColVersionOfRowVector(structure_->Discretization(), structure_->Velnp());
    column_fluid_velocity_ = DRT::UTILS::GetColVersionOfRowVector(fluid_->Discretization(),
        Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)->Velnp());
  }


  std::vector<DRT::Element const*> ele_ptrs(2);
  Teuchos::RCP<std::vector<double>> beam_dofvec = Teuchos::rcp(new std::vector<double>);
  Teuchos::RCP<std::vector<double>> fluid_dofvec = Teuchos::rcp(new std::vector<double>);

  // loop over all (embedded) beam elements
  std::map<int, std::vector<int>>::const_iterator beamelementiterator;
  for (beamelementiterator = pairids->begin(); beamelementiterator != pairids->end();
       beamelementiterator++)
  {
    // add beam elements to the element pair pointer
    ele_ptrs[0] = (structure_->Discretization())->gElement(beamelementiterator->first);


    if (ele_ptrs[0]->Owner() != structure_->Discretization()->Comm().MyPID())
      dserror(
          "For now we can only create the pair on the beam owner, but beam element owner is %i "
          "and "
          "we are on proc %i \n",
          ele_ptrs[0]->Owner(), structure_->Discretization()->Comm().MyPID());

    // loop over all fluid elements, in which the beam element might lie
    for (std::vector<int>::const_iterator fluideleIter = beamelementiterator->second.begin();
         fluideleIter != (beamelementiterator->second).end(); fluideleIter++)
    {
      DRT::Element* fluidele = (fluid_->Discretization())->gElement(*fluideleIter);

      // add fluid element to the element pair pointer
      ele_ptrs[1] = fluidele;

      // Extract current element dofs, i.e. positions and velocities
      ExtractCurrentElementDofs(ele_ptrs, beam_dofvec, fluid_dofvec);

      // Finally tell the bridge to create the pair
      bridge_->CreatePair(ele_ptrs, beam_dofvec, fluid_dofvec);
    }
  }
}

/*----------------------------------------------------------------------*/

// todo take a good look at dispn vs. dispnp and veln vs. velnp!
void ADAPTER::FBIConstraintenforcer::ExtractCurrentElementDofs(
    std::vector<DRT::Element const*> elements, Teuchos::RCP<std::vector<double>>& beam_dofvec,
    Teuchos::RCP<std::vector<double>>& fluid_dofvec) const
{
  std::vector<double> vel_tmp;

  // extract the current position of the beam element from the displacement vector
  BEAMINTERACTION::UTILS::ExtractPosDofVecAbsoluteValues(*(structure_->Discretization()),
      elements[0], column_structure_displacement_,
      *beam_dofvec);  // todo check time step todo get "interface" displacements only for beam
                      // elements todo Only get centerline velocities?!
  ;
  // extract veclocity of the beam element
  BEAMINTERACTION::UTILS::ExtractPosDofVecValues(
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
    }
  }

  // extract current fluid velocities
  BEAMINTERACTION::UTILS::GetCurrentElementDis(
      *(fluid_->Discretization()), elements[1], column_fluid_velocity_, vel_tmp);

  // todo This is a very crude way to separate the pressure from the velocity dofs.. maybe just
  // use a extractor?
  for (unsigned int i = 0; i < vel_tmp.size(); i++)
  {
    if ((i + 1) % 4) fluid_dofvec->push_back(vel_tmp[i]);
  }
}

/*----------------------------------------------------------------------*/

void ADAPTER::FBIConstraintenforcer::PrintViolation()
{
  Teuchos::RCP<Epetra_Vector> violation = LINALG::CreateVector(
      Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)->Velnp()->Map());

  int err =
      Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)
          ->GetCff()
          ->Multiply(false,
              *(Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)->Velnp()), *violation);

  if (err != 0) dserror(" Matrix vector product threw error code %i ", err);

  err = violation->Update(1.0, *AssembleFluidForce(), 0.0);
  if (err != 0) dserror(" Epetra_Vector update threw error code %i ", err);

  double norm, normf, norms;
  double norm_vel;

  Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)->Velnp()->MaxValue(&norm_vel);

  violation->Norm2(&norm);
  if (norm_vel > 1e-15) normf = norm / norm_vel;

  Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_, true)
      ->Velnp()
      ->MaxValue(&norm_vel);
  if (norm_vel > 1e-15) norms = norm / norm_vel;

  std::ofstream log;
  if ((*discretizations_)[1]->Comm().MyPID() == 0)
  {
    std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
    s.append(".penalty");
    log.open(s.c_str(), std::ofstream::app);
    log << " " << norm << " " << normf << " " << norms << std::endl;
  }
}
