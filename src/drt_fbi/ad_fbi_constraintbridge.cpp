/*----------------------------------------------------------------------*/
/*! \file

\brief Abstract class to be overloaded by different adapter implementations connecting the
constraint enforcement technique with a discretization approach for Fluid-beam interaction.

\level 2

*----------------------------------------------------------------------*/
#include "ad_fbi_constraintbridge.H"
#include "fluid_assembly_strategy.H"
#include "fluidblockmatrix_assembly_strategy.H"
#include "../drt_geometry_pair/geometry_pair_line_to_3D_evaluation_data.H"
#include "../drt_beaminteraction/beam_contact_pair.H"
#include "beam_to_fluid_meshtying_pair_factory.H"
#include "beam_to_fluid_meshtying_params.H"
#include "../drt_inpar/inpar_fbi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_sparseoperator.H"

ADAPTER::FBIConstraintBridge::FBIConstraintBridge()
    : beam_interaction_params_(Teuchos::null),
      assemblystrategy_(Teuchos::null),
      meshtying_pairs_(
          Teuchos::rcp(new std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>>)),
      geometry_evaluation_data_(Teuchos::null){};
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridge::Setup(const Epetra_Map* beam_map, const Epetra_Map* fluid_map,
    Teuchos::RCP<LINALG::SparseOperator> fluidmatrix, bool fluidmeshtying)
{
  // Create the beaminteraction data container and set the parameters
  beam_interaction_params_ = Teuchos::rcp(new FBI::BeamToFluidMeshtyingParams());
  beam_interaction_params_->Init();
  beam_interaction_params_->Setup();

  const Teuchos::ParameterList& geometry_parameter_list =
      DRT::Problem::Instance()->FBIParams().sublist("BEAM TO FLUID MESHTYING");

  // Create the beaminteraction data container and set the parameters
  geometry_evaluation_data_ = Teuchos::rcp<GEOMETRYPAIR::LineTo3DEvaluationData>(
      new GEOMETRYPAIR::LineTo3DEvaluationData(geometry_parameter_list));

  if (fluidmeshtying)
  {
    // For the option condensed smat this can be changed by creating a FEMatrix instead of a
    // CRSMatrix!
    if (beam_interaction_params_->GetContactDiscretization() ==
        INPAR::FBI::BeamToFluidDiscretization::mortar)
      dserror("Fluid Meshtying is not supported when using a mortar discretization!");

    assemblystrategy_ = Teuchos::rcp<FBI::UTILS::FBIBlockAssemblyStrategy>(
        new FBI::UTILS::FBIBlockAssemblyStrategy());
  }
  else
    assemblystrategy_ =
        Teuchos::rcp<FBI::UTILS::FBIAssemblyStrategy>(new FBI::UTILS::FBIAssemblyStrategy());
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridge::CreatePair(const std::vector<DRT::Element const*> elements,
    const Teuchos::RCP<std::vector<double>> beam_centerline_dofvec,
    const Teuchos::RCP<std::vector<double>> fluid_nodal_dofvec)
{
  // create a new beaminteratcion pair
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> newinteractionpair =
      FBI::PairFactory::CreatePair(elements, GetParams());

  // create the underlying geometrypair doing the integration (segment or gauss point projection
  // based)
  newinteractionpair->CreateGeometryPair(GetGeometryData());
  newinteractionpair->Init(GetParams(), elements);
  newinteractionpair->Setup();

  // hand in the current position and velocities of the participating elements
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
  // hand in the current position and velocities of the participating elements
  interactionpair->ResetState(*beam_centerline_dofvec, *fluid_nodal_dofvec);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIConstraintBridge::Clear()
{
  // Delete all pairs and segmentation information
  meshtying_pairs_->clear();
  geometry_evaluation_data_->Clear();
}
