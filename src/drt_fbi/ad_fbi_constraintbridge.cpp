/*----------------------------------------------------------------------*/
/*! \file

\brief Abstract class to be overloaded by different adapter implementations connecting the
constraint enforcement technique with a discretization approach for Fluid-beam interaction.

\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/
#include "ad_fbi_constraintbridge.H"
#include "../drt_geometry_pair/geometry_pair_line_to_3D_evaluation_data.H"
#include "../drt_beaminteraction/beam_contact_pair.H"
#include "beam_to_fluid_meshtying_pair_factory.H"
#include "beam_to_fluid_meshtying_params.H"
#include "../drt_lib/drt_globalproblem.H"

ADAPTER::FBIConstraintBridge::FBIConstraintBridge()
    : meshtying_pairs_(
          Teuchos::rcp(new std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>)),
      beam_interaction_params_(Teuchos::null),
      geometry_evaluation_data_(Teuchos::null){};
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridge::Setup(const Epetra_Map* beam_map, const Epetra_Map* fluid_map)
{
  beam_interaction_params_ = Teuchos::rcp(new FBI::BeamToFluidMeshtyingParams());
  beam_interaction_params_->Init();
  beam_interaction_params_->Setup();

  const Teuchos::ParameterList& geometry_parameter_list =
      DRT::Problem::Instance()->FBIParams().sublist("BEAM TO FLUID MESHTYING");

  geometry_evaluation_data_ = Teuchos::rcp<GEOMETRYPAIR::LineTo3DEvaluationData>(
      new GEOMETRYPAIR::LineTo3DEvaluationData(geometry_parameter_list));
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/// Creates a volume_beam_meshtying pair
void ADAPTER::FBIConstraintBridge::CreatePair(const std::vector<DRT::Element const*> elements,
    const Teuchos::RCP<std::vector<double>> beam_centerline_dofvec,
    const Teuchos::RCP<std::vector<double>> fluid_nodal_dofvec)
{
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> newinteractionpair =
      FBI::PairFactory::CreatePair(elements, GetParams());
  newinteractionpair->CreateGeometryPair(GetGeometryData());
  newinteractionpair->Init(GetParams(), elements);
  newinteractionpair->Setup();
  ResetPair(beam_centerline_dofvec, fluid_nodal_dofvec, newinteractionpair);
  // add to list of current contact pairs
  meshtying_pairs_->push_back(newinteractionpair);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridge::ResetPair(
    const Teuchos::RCP<const std::vector<double>> beam_centerline_dofvec,
    const Teuchos::RCP<const std::vector<double>> fluid_nodal_dofvec,
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> interactionpair)
{
  interactionpair->ResetState(*beam_centerline_dofvec, *fluid_nodal_dofvec);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridge::ResetBridge()
{
  meshtying_pairs_->clear();
  geometry_evaluation_data_->Reset();
}
