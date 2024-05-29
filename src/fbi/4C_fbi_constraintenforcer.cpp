/*----------------------------------------------------------------------*/
/*! \file

\brief Abstract class to be overloaded by different constraint enforcement techniques for fluid-beam
interaction.

\level 3

*----------------------------------------------------------------------*/

#include "4C_fbi_constraintenforcer.hpp"

#include "4C_adapter_fld_fbi_movingboundary.hpp"
#include "4C_adapter_str_fbiwrapper.hpp"
#include "4C_beaminteraction_calc_utils.hpp"  // todo put this into bridge to keep everything beam specific in there
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_binstrategy.hpp"
#include "4C_discretization_fem_general_element.hpp"
#include "4C_fbi_adapter_constraintbridge.hpp"
#include "4C_fbi_adapter_constraintbridge_penalty.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_output_params.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_params.hpp"
#include "4C_fbi_immersed_geometry_coupler.hpp"
#include "4C_fluid_utils.hpp"
#include "4C_geometry_pair.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fbi.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_discret_faces.hpp"
#include "4C_lib_node.hpp"
#include "4C_linalg_blocksparsematrix.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_mapextractor.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_rebalance_binning_based.hpp"

#include <iostream>

FOUR_C_NAMESPACE_OPEN

ADAPTER::FBIConstraintenforcer::FBIConstraintenforcer(
    Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge,
    Teuchos::RCP<FBI::FBIGeometryCoupler> geometrycoupler)
    : fluid_(Teuchos::null),
      structure_(Teuchos::null),
      discretizations_(),
      bridge_(bridge),
      geometrycoupler_(geometrycoupler),
      column_structure_displacement_(Teuchos::null),
      column_structure_velocity_(Teuchos::null),
      column_fluid_velocity_(Teuchos::null),
      velocity_pressure_splitter_(Teuchos::rcp(new CORE::LINALG::MapExtractor()))
{
}

/*----------------------------------------------------------------------*/

void ADAPTER::FBIConstraintenforcer::Setup(Teuchos::RCP<ADAPTER::FSIStructureWrapper> structure,
    Teuchos::RCP<ADAPTER::FluidMovingBoundary> fluid)
{
  fluid_ = fluid;
  structure_ = structure;
  discretizations_.push_back(structure_->discretization());
  discretizations_.push_back(fluid_->discretization());

  CORE::LINALG::CreateMapExtractorFromDiscretization(
      *(fluid_->discretization()), 3, *velocity_pressure_splitter_);

  bool meshtying =
      (GLOBAL::Problem::Instance()->FluidDynamicParams().get<std::string>("MESHTYING") != "no");

  Teuchos::RCP<CORE::LINALG::SparseOperator> fluidmatrix(Teuchos::null);

  if (meshtying)
  {
    if (structure_->discretization()->Comm().NumProc() > 1)
      FOUR_C_THROW(
          "Currently fluid mesh tying can only be used for serial computations, since offproc "
          "assembly is not supported. Once the coupling matrices are computed by the fluid element "
          "owner, this will change.");

    fluidmatrix = (Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)->GetMeshtying())
                      ->init_system_matrix();
  }
  else
  {
    fluidmatrix = Teuchos::rcp(
        new CORE::LINALG::SparseMatrix(*(fluid_->discretization()->dof_row_map()), 30, true, true,
            CORE::LINALG::SparseMatrix::FE_MATRIX));  // todo Is there a better estimator?
  }

  bridge_->Setup(structure_->discretization()->dof_row_map(),
      fluid_->discretization()->dof_row_map(), fluidmatrix, meshtying);
  if (structure_->discretization()->Comm().NumProc() > 1)
  {
    geometrycoupler_->ExtendBeamGhosting(*(structure->discretization()));

    // After ghosting we need to explicitly set up the MultiMapExtractor again
    Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_, true)
        ->setup_multi_map_extractor();
  }

  geometrycoupler_->Setup(
      discretizations_, CORE::REBALANCE::GetColVersionOfRowVector(
                            structure_->discretization(), structure_->Dispnp()));
}

/*----------------------------------------------------------------------*/

void ADAPTER::FBIConstraintenforcer::Evaluate()
{
  // We use the column vectors here, because currently the search is based on neighboring nodes,
  // but the element pairs are created using the elements needing all information on all their
  // DOFs
  column_structure_displacement_ =
      CORE::REBALANCE::GetColVersionOfRowVector(structure_->discretization(), structure_->Dispnp());
  column_structure_velocity_ =
      CORE::REBALANCE::GetColVersionOfRowVector(structure_->discretization(), structure_->Velnp());
  column_fluid_velocity_ = CORE::REBALANCE::GetColVersionOfRowVector(fluid_->discretization(),
      Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)->Velnp());

  geometrycoupler_->UpdateBinning(discretizations_[0], column_structure_displacement_);

  // Before each search we delete all pair and segment information
  bridge_->Clear();
  bridge_->ResetBridge();

  // Do the search in the geometrycoupler_ and return the possible pair ids
  Teuchos::RCP<std::map<int, std::vector<int>>> pairids = geometrycoupler_->Search(discretizations_,
      column_structure_displacement_);  // todo make this a vector? At some point we probably
                                        // need the ale displacements as well

  // For now we need to separate the pair creation from the search, since the search takes place
  // on the fluid elements owner, while (for now) the pair has to be created on the beam element
  // owner
  create_pairs(pairids);

  // Create all needed matrix and vector contributions based on the current state
  bridge_->Evaluate(discretizations_[0], discretizations_[1],
      Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)->Velnp(), structure_->Velnp());
}

/*----------------------------------------------------------------------*/

Teuchos::RCP<Epetra_Vector> ADAPTER::FBIConstraintenforcer::StructureToFluid(int step)
{
  // todo only access the parameter list once

  // Check if we want to couple the fluid
  const Teuchos::ParameterList& fbi = GLOBAL::Problem::Instance()->FBIParams();
  if (Teuchos::getIntegralValue<INPAR::FBI::BeamToFluidCoupling>(fbi, "COUPLING") !=
          INPAR::FBI::BeamToFluidCoupling::solid &&
      fbi.get<int>("STARTSTEP") < step)
  {
    // Assemble the fluid stiffness matrix and hand it to the fluid solver
    Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)
        ->set_coupling_contributions(assemble_fluid_coupling_matrix());

    // Assemble the fluid force vector and hand it to the fluid solver
    fluid_->apply_interface_values(assemble_fluid_coupling_residual());
  }

  // return the current struture velocity
  return Teuchos::rcp_dynamic_cast<ADAPTER::FBIStructureWrapper>(structure_, true)
      ->extract_interface_velnp();
};

/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintenforcer::recompute_coupling_without_pair_creation()
{
  // Before each search we delete all pair and segment information
  bridge_->ResetBridge();

  reset_all_pair_states();

  // Create all needed matrix and vector contributions based on the current state
  bridge_->Evaluate(discretizations_[0], discretizations_[1],
      Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)->Velnp(), structure_->Velnp());
};

/*----------------------------------------------------------------------*/
// return the structure force
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIConstraintenforcer::FluidToStructure()
{
  return assemble_structure_coupling_residual();
};

/*----------------------------------------------------------------------*/

// For now we need to separate the pair creation from the search, since the search takes place on
// the fluid elements owner, while (for now) the pair has to be created on the beam element owner
void ADAPTER::FBIConstraintenforcer::create_pairs(
    Teuchos::RCP<std::map<int, std::vector<int>>> pairids)
{
  if ((structure_->discretization())->Comm().NumProc() > 1)
  {
    // The geometrycoupler takes care of all MPI communication that needs to be done before the
    // pairs can finally be created
    geometrycoupler_->PreparePairCreation(discretizations_, pairids);

    column_structure_displacement_ = CORE::REBALANCE::GetColVersionOfRowVector(
        structure_->discretization(), structure_->Dispnp());
    column_structure_velocity_ = CORE::REBALANCE::GetColVersionOfRowVector(
        structure_->discretization(), structure_->Velnp());
    column_fluid_velocity_ = CORE::REBALANCE::GetColVersionOfRowVector(fluid_->discretization(),
        Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)->Velnp());
  }


  std::vector<CORE::Elements::Element const*> ele_ptrs(2);
  std::vector<double> beam_dofvec = std::vector<double>();
  std::vector<double> fluid_dofvec = std::vector<double>();

  // loop over all (embedded) beam elements
  std::map<int, std::vector<int>>::const_iterator beamelementiterator;
  for (beamelementiterator = pairids->begin(); beamelementiterator != pairids->end();
       beamelementiterator++)
  {
    // add beam elements to the element pair pointer
    ele_ptrs[0] = (structure_->discretization())->gElement(beamelementiterator->first);


    if (ele_ptrs[0]->Owner() != structure_->discretization()->Comm().MyPID())
      FOUR_C_THROW(
          "For now we can only create the pair on the beam owner, but beam element owner is %i "
          "and "
          "we are on proc %i \n",
          ele_ptrs[0]->Owner(), structure_->discretization()->Comm().MyPID());

    // loop over all fluid elements, in which the beam element might lie
    for (std::vector<int>::const_iterator fluideleIter = beamelementiterator->second.begin();
         fluideleIter != (beamelementiterator->second).end(); fluideleIter++)
    {
      CORE::Elements::Element* fluidele = (fluid_->discretization())->gElement(*fluideleIter);

      // add fluid element to the element pair pointer
      ele_ptrs[1] = fluidele;

      // Extract current element dofs, i.e. positions and velocities
      extract_current_element_dofs(ele_ptrs, beam_dofvec, fluid_dofvec);

      // Finally tell the bridge to create the pair
      bridge_->CreatePair(ele_ptrs, beam_dofvec, fluid_dofvec);
    }
  }
}
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintenforcer::reset_all_pair_states()
{
  // Get current state
  column_structure_displacement_ =
      CORE::REBALANCE::GetColVersionOfRowVector(structure_->discretization(), structure_->Dispnp());
  column_structure_velocity_ =
      CORE::REBALANCE::GetColVersionOfRowVector(structure_->discretization(), structure_->Velnp());
  column_fluid_velocity_ = CORE::REBALANCE::GetColVersionOfRowVector(fluid_->discretization(),
      Teuchos::rcp_dynamic_cast<ADAPTER::FBIFluidMB>(fluid_, true)->Velnp());

  std::vector<CORE::Elements::Element const*> ele_ptrs(2);
  std::vector<double> beam_dofvec = std::vector<double>();
  std::vector<double> fluid_dofvec = std::vector<double>();

  for (auto pairiterator = bridge_->GetPairs()->begin(); pairiterator != bridge_->GetPairs()->end();
       pairiterator++)
  {
    ele_ptrs[0] = (*pairiterator)->Element1();
    ele_ptrs[1] = (*pairiterator)->Element2();

    // Extract current element dofs, i.e. positions and velocities
    extract_current_element_dofs(ele_ptrs, beam_dofvec, fluid_dofvec);

    // Finally tell the bridge to create the pair
    bridge_->ResetPair(beam_dofvec, fluid_dofvec, *pairiterator);
  }
}
/*----------------------------------------------------------------------*/

void ADAPTER::FBIConstraintenforcer::extract_current_element_dofs(
    std::vector<CORE::Elements::Element const*> elements, std::vector<double>& beam_dofvec,
    std::vector<double>& fluid_dofvec) const
{
  std::vector<double> vel_tmp;

  // extract the current position of the beam element from the displacement vector
  BEAMINTERACTION::UTILS::ExtractPosDofVecAbsoluteValues(*(structure_->discretization()),
      elements[0], column_structure_displacement_,
      beam_dofvec);  // todo get "interface" displacements only for beam
                     // elements
  // extract velocity of the beam element
  BEAMINTERACTION::UTILS::ExtractPosDofVecValues(
      *(structure_->discretization()), elements[0], column_structure_velocity_, vel_tmp);

  for (double val : vel_tmp) beam_dofvec.push_back(val);

  vel_tmp.clear();
  // extract the current positions and velocities of the fluid element todo only valid for fixed
  // grid, not for ALE
  fluid_dofvec.clear();
  const DRT::Node* const* fluidnodes = elements[1]->Nodes();
  for (int lid = 0; lid < elements[1]->num_node(); ++lid)
  {
    for (int dim = 0; dim < 3; dim++)
    {
      fluid_dofvec.push_back(fluidnodes[lid]->X()[dim]);
    }
  }

  // extract current fluid velocities
  BEAMINTERACTION::UTILS::GetCurrentElementDis(
      *(fluid_->discretization()), elements[1], column_fluid_velocity_, vel_tmp);

  // todo This is a very crude way to separate the pressure from the velocity dofs.. maybe just
  // use an extractor?
  for (unsigned int i = 0; i < vel_tmp.size(); i++)
  {
    if ((i + 1) % 4) fluid_dofvec.push_back(vel_tmp[i]);
  }
}

/*----------------------------------------------------------------------*/

void ADAPTER::FBIConstraintenforcer::SetBinning(Teuchos::RCP<BINSTRATEGY::BinningStrategy> binning)
{
  geometrycoupler_->SetBinning(binning);
};

FOUR_C_NAMESPACE_CLOSE
