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
#include "immersed_geometry_coupler_fbi.H"

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
    Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_, true)
        ->SetupMultiMapExtractor();
  }
}

void ADAPTER::FBIConstraintenforcer::Evaluate()
{
  column_structure_displacement_ =
      DRT::UTILS::GetColVersionOfRowVector(structure_->Discretization(), structure_->Dispnp());
  column_structure_velocity_ =
      DRT::UTILS::GetColVersionOfRowVector(structure_->Discretization(), structure_->Velnp());
  column_fluid_velocity_ = DRT::UTILS::GetColVersionOfRowVector(fluid_->Discretization(),
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidBeamImmersed>(fluid_, true)->Velnp());

  bridge_->ResetBridge();

  Teuchos::RCP<std::map<int, std::vector<int>>> pairids = geometrycoupler_->Search(discretizations_,
      column_structure_displacement_, column_structure_velocity_,
      column_fluid_velocity_);  // todo make this a vector?

  CreatePairs(pairids);

  bridge_->Evaluate(discretizations_);
}

Teuchos::RCP<Epetra_Vector> ADAPTER::FBIConstraintenforcer::StructureToFluid()
{
  const Teuchos::ParameterList& fbi = DRT::Problem::Instance()->FBIParams();
  if (DRT::INPUT::IntegralValue<int>(fbi, "COUPLING") != INPAR::FBI::BeamToFluidCoupling::solid)
  {
    Teuchos::rcp_dynamic_cast<ADAPTER::FluidBeamImmersed>(fluid_, true)
        ->SetCouplingContributions(AssembleFluidStiffness());
    fluid_->ApplyInterfaceValues(AssembleFluidForce());
  }

  return Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_, true)
      ->ExtractInterfaceVelnp();
};


// This is the interface to the outside world regarding the force.
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIConstraintenforcer::FluidToStructure()
{
  return AssembleStructureForce();
};

// For now we need to separate the pair creation from the search, since the search takes place on
// the fluid elements owner, while (for now) the pair has to be created on the beam element owner
void ADAPTER::FBIConstraintenforcer::CreatePairs(
    Teuchos::RCP<std::map<int, std::vector<int>>> pairids)
{
  if ((structure_->Discretization())->Comm().NumProc() > 1)
  {
    geometrycoupler_->PreparePairCreation(discretizations_, pairids);
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
  for (beamelementiterator = pairids->begin(); beamelementiterator != pairids->end();
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

  BEAMINTERACTION::UTILS::GetCurrentElementDis(
      *(fluid_->Discretization()), elements[1], column_fluid_velocity_, vel_tmp);

  // todo This is a very crude way to seperate the pressure from the velocity dofs.. maybe just use
  // a extractor?
  for (unsigned int i = 0; i < vel_tmp.size(); i++)
  {
    if ((i + 1) % 4) fluid_dofvec->push_back(vel_tmp[i]);
  }
}
