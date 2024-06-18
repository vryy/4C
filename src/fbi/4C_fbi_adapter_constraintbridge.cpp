/*----------------------------------------------------------------------*/
/*! \file

\brief Abstract class to be overloaded by different adapter implementations connecting the
constraint enforcement technique with a discretization approach for Fluid-beam interaction.

\level 2

*----------------------------------------------------------------------*/
#include "4C_fbi_adapter_constraintbridge.hpp"

#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_pair_factory.hpp"
#include "4C_fbi_beam_to_fluid_meshtying_params.hpp"
#include "4C_fbi_fluid_assembly_strategy.hpp"
#include "4C_fbi_fluidblockmatrix_assembly_strategy.hpp"
#include "4C_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fbi.hpp"
#include "4C_linalg_sparseoperator.hpp"

FOUR_C_NAMESPACE_OPEN

Adapter::FBIConstraintBridge::FBIConstraintBridge()
    : beam_interaction_params_(Teuchos::null),
      assemblystrategy_(Teuchos::null),
      meshtying_pairs_(
          Teuchos::rcp(new std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>)),
      geometry_evaluation_data_(Teuchos::null){};
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIConstraintBridge::setup(const Epetra_Map* beam_map, const Epetra_Map* fluid_map,
    Teuchos::RCP<Core::LinAlg::SparseOperator> fluidmatrix, bool fluidmeshtying)
{
  // Create the beaminteraction data container and set the parameters
  beam_interaction_params_ = Teuchos::rcp(new FBI::BeamToFluidMeshtyingParams());
  beam_interaction_params_->Init();
  beam_interaction_params_->setup();

  const Teuchos::ParameterList& geometry_parameter_list =
      Global::Problem::Instance()->FBIParams().sublist("BEAM TO FLUID MESHTYING");

  // Create the beaminteraction data container and set the parameters
  geometry_evaluation_data_ = Teuchos::rcp<GEOMETRYPAIR::LineTo3DEvaluationData>(
      new GEOMETRYPAIR::LineTo3DEvaluationData(geometry_parameter_list));

  if (fluidmeshtying)
  {
    // For the option condensed smat this can be changed by creating a FEMatrix instead of a
    // CRSMatrix!
    if (beam_interaction_params_->get_contact_discretization() ==
        Inpar::FBI::BeamToFluidDiscretization::mortar)
      FOUR_C_THROW("Fluid Meshtying is not supported when using a mortar discretization!");

    assemblystrategy_ = Teuchos::rcp<FBI::UTILS::FBIBlockAssemblyStrategy>(
        new FBI::UTILS::FBIBlockAssemblyStrategy());
  }
  else
    assemblystrategy_ =
        Teuchos::rcp<FBI::UTILS::FBIAssemblyStrategy>(new FBI::UTILS::FBIAssemblyStrategy());
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIConstraintBridge::CreatePair(
    const std::vector<Core::Elements::Element const*> elements,
    const std::vector<double> beam_centerline_dofvec, const std::vector<double> fluid_nodal_dofvec)
{
  // create a new beaminteratcion pair
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> newinteractionpair =
      FBI::PairFactory::CreatePair(elements, GetParams());

  // create the underlying geometrypair doing the integration (segment or gauss point projection
  // based)
  newinteractionpair->CreateGeometryPair(elements[0], elements[1], GetGeometryData());
  newinteractionpair->Init(GetParams(), elements);
  newinteractionpair->setup();

  // hand in the current position and velocities of the participating elements
  ResetPair(beam_centerline_dofvec, fluid_nodal_dofvec, newinteractionpair);

  // add to list of current contact pairs
  meshtying_pairs_->push_back(newinteractionpair);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIConstraintBridge::ResetPair(const std::vector<double> beam_centerline_dofvec,
    const std::vector<double> fluid_nodal_dofvec,
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> interactionpair)
{
  // hand in the current position and velocities of the participating elements
  interactionpair->ResetState(beam_centerline_dofvec, fluid_nodal_dofvec);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::FBIConstraintBridge::Clear()
{
  // Delete all pairs and segmentation information
  meshtying_pairs_->clear();
  geometry_evaluation_data_->Clear();
}

FOUR_C_NAMESPACE_CLOSE
