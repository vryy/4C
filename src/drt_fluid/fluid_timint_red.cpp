/*-----------------------------------------------------------*/
/*! \file

\brief Basic time integration driver for reduced models


\level 2

*/
/*-----------------------------------------------------------*/

#include "fluid_timint_red.H"
#include "fluid_coupling_red_models.H"
#include "fluid_volumetric_surfaceFlow_condition.H"
#include "../drt_adapter/ad_art_net.H"
#include "fluid_meshtying.H"
#include "../drt_lib/drt_locsys.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntRedModels::TimIntRedModels(const Teuchos::RCP<DRT::Discretization>& actdis,
    const Teuchos::RCP<LINALG::Solver>& solver, const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis, solver, params, output, alefluid),
      traction_vel_comp_adder_bc_(Teuchos::null),
      coupled3D_redDbc_art_(Teuchos::null),
      ART_timeInt_(Teuchos::null),
      coupled3D_redDbc_airways_(Teuchos::null),
      airway_imp_timeInt_(Teuchos::null),
      vol_surf_flow_bc_(Teuchos::null),
      vol_surf_flow_bc_maps_(Teuchos::null),
      vol_flow_rates_bc_extractor_(Teuchos::null),
      strong_redD_3d_coupling_(false)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize algorithm                                rasthofer 04/14 |
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::Init()
{
  // Vectors associated to boundary conditions
  // -----------------------------------------

  // create the volumetric-surface-flow condition
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispn_);
  }

  vol_surf_flow_bc_ = Teuchos::rcp(new UTILS::FluidVolumetricSurfaceFlowWrapper(discret_, dta_));

  // evaluate the map of the womersley bcs
  vol_flow_rates_bc_extractor_ = Teuchos::rcp(new FLD::UTILS::VolumetricFlowMapExtractor());
  vol_flow_rates_bc_extractor_->Setup(*discret_);
  vol_surf_flow_bc_maps_ =
      Teuchos::rcp(new Epetra_Map(*(vol_flow_rates_bc_extractor_->VolumetricSurfaceFlowCondMap())));

  // -------------------------------------------------------------------
  // Initialize the reduced models
  // -------------------------------------------------------------------

  strong_redD_3d_coupling_ = false;
  if (params_->get<std::string>("Strong 3D_redD coupling", "no") == "yes")
    strong_redD_3d_coupling_ = true;

  {
    ART_timeInt_ = dyn_art_net_drt(true);
    // Check if one-dimensional artery network problem exist
    if (ART_timeInt_ != Teuchos::null)
    {
      IO::DiscretizationWriter output_redD(ART_timeInt_->Discretization());
      discret_->ClearState();
      discret_->SetState("velaf", zeros_);
      if (alefluid_)
      {
        discret_->SetState("dispnp", dispnp_);
      }
      coupled3D_redDbc_art_ =
          Teuchos::rcp(new UTILS::Fluid_couplingWrapper<ADAPTER::ArtNet>(discret_,
              ART_timeInt_->Discretization(), ART_timeInt_, output_redD, dta_, ART_timeInt_->Dt()));
    }


    airway_imp_timeInt_ = dyn_red_airways_drt(true);
    // Check if one-dimensional artery network problem exist
    if (airway_imp_timeInt_ != Teuchos::null)
    {
      IO::DiscretizationWriter output_redD(airway_imp_timeInt_->Discretization());
      discret_->ClearState();
      discret_->SetState("velaf", zeros_);
      if (alefluid_)
      {
        discret_->SetState("dispnp", dispnp_);
      }
      coupled3D_redDbc_airways_ =
          Teuchos::rcp(new UTILS::Fluid_couplingWrapper<AIRWAY::RedAirwayImplicitTimeInt>(discret_,
              airway_imp_timeInt_->Discretization(), airway_imp_timeInt_, output_redD, dta_,
              airway_imp_timeInt_->Dt()));
    }


    zeros_->PutScalar(0.0);  // just in case of change
  }

  traction_vel_comp_adder_bc_ = Teuchos::rcp(new UTILS::TotalTractionCorrector(discret_, dta_));


  // ------------------------------------------------------------------------------
  // Check, if features are used with the locsys manager that are not supported,
  // or better, not implemented yet.
  // ------------------------------------------------------------------------------
  if (locsysman_ != Teuchos::null)
  {
    // Models
    if ((ART_timeInt_ != Teuchos::null) or (airway_imp_timeInt_ != Teuchos::null))
    {
      dserror("No problem types involving airways are supported for use with locsys conditions!");
    }
  }

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                    bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntRedModels::~TimIntRedModels() { return; }

/*----------------------------------------------------------------------*
 | evaluate special boundary conditions                        bk 12/13 |
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::DoProblemSpecificBoundaryConditions()
{
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
  }

  // Check if one-dimensional artery network problem exist
  if (ART_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_art_->EvaluateDirichlet(velnp_, *(dbcmaps_->CondMap()), time_);
  }
  // update the 3D-to-reduced_D coupling data
  // Check if one-dimensional artery network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_airways_->EvaluateDirichlet(velnp_, *(dbcmaps_->CondMap()), time_);
  }

  // Evaluate the womersley velocities
  vol_surf_flow_bc_->EvaluateVelocities(velnp_, time_);

  return;
}

/*----------------------------------------------------------------------*
| Update3DToReduced in AssembleMatAndRHS                       bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntRedModels::Update3DToReducedMatAndRHS()
{
  discret_->ClearState();

  discret_->SetState("velaf", velnp_);
  discret_->SetState("hist", hist_);

  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
  }

  // Check if one-dimensional artery network problem exist
  if (ART_timeInt_ != Teuchos::null)
  {
    if (strong_redD_3d_coupling_)
    {
      coupled3D_redDbc_art_->LoadState();
      coupled3D_redDbc_art_->FlowRateCalculation(time_, dta_);
      coupled3D_redDbc_art_->ApplyBoundaryConditions(time_, dta_, theta_);
    }
    coupled3D_redDbc_art_->UpdateResidual(residual_);
  }
  // update the 3D-to-reduced_D coupling data
  // Check if one-dimensional artery network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    if (strong_redD_3d_coupling_)
    {
      coupled3D_redDbc_airways_->LoadState();
      coupled3D_redDbc_airways_->FlowRateCalculation(time_, dta_);
      coupled3D_redDbc_airways_->ApplyBoundaryConditions(time_, dta_, theta_);
    }
    coupled3D_redDbc_airways_->UpdateResidual(residual_);
  }

  //----------------------------------------------------------------------
  // add the traction velocity component
  //----------------------------------------------------------------------

  traction_vel_comp_adder_bc_->EvaluateVelocities(velnp_, time_, theta_, dta_);
  traction_vel_comp_adder_bc_->UpdateResidual(residual_);

  discret_->ClearState();
  return;
}

/*----------------------------------------------------------------------*
| call Update3DToReducedMatAndRHS                              bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntRedModels::SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams)
{
  // these are the only routines that have to be called in AssembleMatAndRHS
  // before Evaluate in the RedModels case
  Update3DToReducedMatAndRHS();

  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector of ReducedD problem to binio   ismail 01/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::OutputReducedD()
{
  // output of solution
  if (step_ % upres_ == 0)
  {
    // write reduced model problem
    // Check if one-dimensional artery network problem exist
    if (ART_timeInt_ != Teuchos::null)
    {
      Teuchos::RCP<Teuchos::ParameterList> redD_export_params;
      redD_export_params = Teuchos::rcp(new Teuchos::ParameterList());

      redD_export_params->set<int>("step", step_);
      redD_export_params->set<int>("upres", upres_);
      redD_export_params->set<int>("uprestart", uprestart_);
      redD_export_params->set<double>("time", time_);

      ART_timeInt_->Output(true, redD_export_params);
    }

    // Check if one-dimensional artery network problem exist
    if (airway_imp_timeInt_ != Teuchos::null)
    {
      Teuchos::RCP<Teuchos::ParameterList> redD_export_params;
      redD_export_params = Teuchos::rcp(new Teuchos::ParameterList());

      redD_export_params->set<int>("step", step_);
      redD_export_params->set<int>("upres", upres_);
      redD_export_params->set<int>("uprestart", uprestart_);
      redD_export_params->set<double>("time", time_);

      airway_imp_timeInt_->Output(true, redD_export_params);
    }
  }
  return;
}  // FLD::TimIntRedModels::OutputReducedD

/*----------------------------------------------------------------------*
 | read some additional data in restart                         bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_, step);

  vol_surf_flow_bc_->ReadRestart(reader);

  traction_vel_comp_adder_bc_->ReadRestart(reader);

  // Read restart of one-dimensional arterial network
  if (ART_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_art_->ReadRestart(reader);
  }
  // Check if zero-dimensional airway network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_airways_->ReadRestart(reader);
  }

  ReadRestartReducedD(step);

  return;
}

/*----------------------------------------------------------------------*
 |                                                          ismail 01/13|
 -----------------------------------------------------------------------*/
void FLD::TimIntRedModels::ReadRestartReducedD(int step)
{
  // Check if one-dimensional artery network problem exist
  if (ART_timeInt_ != Teuchos::null)
  {
    ART_timeInt_->ReadRestart(step, true);
  }

  // Check if one-dimensional artery network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    airway_imp_timeInt_->ReadRestart(step, true);
  }
}  // FLD::TimIntRedModels::ReadRestartReadRestart(int step)

/*----------------------------------------------------------------------*
 | do some additional steps in SetupMeshtying                   bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::SetupMeshtying()
{
  FluidImplicitTimeInt::SetupMeshtying();
  // Volume surface flow conditions are treated in the same way as Dirichlet condition.
  // Therefore, a volume surface flow condition cannot be defined on the same nodes as the
  // slave side of an internal interface
  // Solution:  Exclude those nodes of your surface
  // but:       The resulting inflow rate (based on the area)
  //            as well as the profile will be different
  //            since it is based on a different surface discretization!!

  if (vol_surf_flow_bc_maps_->NumGlobalElements() != 0)
  {
    meshtying_->CheckOverlappingBC(vol_surf_flow_bc_maps_);
    meshtying_->DirichletOnMaster(vol_surf_flow_bc_maps_);
  }
  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector to binio                        gammi 04/07|
 | overloading function                                         bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::Output()
{
  FluidImplicitTimeInt::Output();
  // output of solution
  if (step_ % upres_ == 0)
  {
    vol_surf_flow_bc_->Output(*output_);
    traction_vel_comp_adder_bc_->Output(*output_);

    if (uprestart_ != 0 && step_ % uprestart_ == 0)  // add restart data
    {
      // Check if one-dimensional artery network problem exist
      if (ART_timeInt_ != Teuchos::null)
      {
        coupled3D_redDbc_art_->WriteRestart(*output_);
      }
      // Check if zero-dimensional airway network problem exist
      if (airway_imp_timeInt_ != Teuchos::null)
      {
        coupled3D_redDbc_airways_->WriteRestart(*output_);
      }
    }
  }
  // write restart also when uprestart_ is not a integer multiple of upres_
  else if (uprestart_ > 0 && step_ % uprestart_ == 0)
  {
    // write reduced model problem
    // Check if one-dimensional artery network problem exist
    if (ART_timeInt_ != Teuchos::null)
    {
      coupled3D_redDbc_art_->WriteRestart(*output_);
    }
    // Check if zero-dimensional airway network problem exist
    if (airway_imp_timeInt_ != Teuchos::null)
    {
      coupled3D_redDbc_airways_->WriteRestart(*output_);
    }
  }

  OutputReducedD();

  return;
}  // TimIntRedModels::Output

/*----------------------------------------------------------------------*
 | read some additional data in restart                         bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::InsertVolumetricSurfaceFlowCondVector(
    Teuchos::RCP<Epetra_Vector> vel, Teuchos::RCP<Epetra_Vector> res)
{
  // -------------------------------------------------------------------
  // take surface volumetric flow rate into account
  //    Teuchos::RCP<Epetra_Vector> temp_vec = Teuchos::rcp(new
  //    Epetra_Vector(*vol_surf_flow_bc_maps_,true)); vol_surf_flow_bc_->InsertCondVector( *temp_vec
  //    , *residual_);
  // -------------------------------------------------------------------
  vol_flow_rates_bc_extractor_->InsertVolumetricSurfaceFlowCondVector(
      vol_flow_rates_bc_extractor_->ExtractVolumetricSurfaceFlowCondVector(vel), res);

  return;
}

/*----------------------------------------------------------------------*
 | prepare AVM3-based scale separation                         vg 10/08 |
 | overloaded in TimIntRedModelsModels and TimIntLoma        bk 12/13 |
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::AVM3Preparation()
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("           + avm3");

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // necessary here, because some application time integrations add something to the residual
  // before the Neumann loads are added
  residual_->PutScalar(0.0);

  // Maybe this needs to be inserted in case of impedanceBC + AVM3
  //  if (nonlinearbc_ && isimpedancebc_)
  //  {
  //    // add impedance Neumann loads
  //    impedancebc_->UpdateResidual(residual_);
  //  }

  AVM3AssembleMatAndRHS(eleparams);

  // apply Womersley as a Dirichlet BC
  LINALG::ApplyDirichlettoSystem(sysmat_, incvel_, residual_, zeros_, *(vol_surf_flow_bc_maps_));

  // get scale-separation matrix
  AVM3GetScaleSeparationMatrix();

  return;
}  // TimIntRedModels::AVM3Preparation

/*----------------------------------------------------------------------*
 | RedModels - specific BC in LinearRelaxationSolve            bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::CustomSolve(Teuchos::RCP<Epetra_Vector> relax)
{
  // apply Womersley as a Dirichlet BC
  LINALG::ApplyDirichlettoSystem(incvel_, residual_, relax, *(vol_surf_flow_bc_maps_));

  // apply Womersley as a Dirichlet BC
  sysmat_->ApplyDirichlet(*(vol_surf_flow_bc_maps_));

  return;
}

/*----------------------------------------------------------------------*
 | RedModels - prepare time step                            ismail 06/14|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::PrepareTimeStep()
{
  FluidImplicitTimeInt::PrepareTimeStep();

  discret_->ClearState();
  discret_->SetState("velaf", velnp_);
  discret_->SetState("hist", hist_);

  if (alefluid_) discret_->SetState("dispnp", dispnp_);

  // Check if one-dimensional artery network problem exist
  if (ART_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_art_->SaveState();
    coupled3D_redDbc_art_->FlowRateCalculation(time_, dta_);
    coupled3D_redDbc_art_->ApplyBoundaryConditions(time_, dta_, theta_);
  }


  // Check if one-dimensional artery network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_airways_->SaveState();
    coupled3D_redDbc_airways_->FlowRateCalculation(time_, dta_);
    coupled3D_redDbc_airways_->ApplyBoundaryConditions(time_, dta_, theta_);
  }

  discret_->ClearState();

  return;
}


/*----------------------------------------------------------------------*
 | Apply Womersley bc to shapederivatives                       bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::AssembleMatAndRHS()
{
  FluidImplicitTimeInt::AssembleMatAndRHS();

  if (shapederivatives_ != Teuchos::null)
  {
    // apply the womersley bc as a dirichlet bc
    shapederivatives_->ApplyDirichlet(*(vol_surf_flow_bc_maps_), false);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Apply Womersley bc to system                                 bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::ApplyDirichletToSystem()
{
  FluidImplicitTimeInt::ApplyDirichletToSystem();

  if (LocsysManager() != Teuchos::null)
  {
    // apply Womersley as a Dirichlet BC
    LINALG::ApplyDirichlettoSystem(
        sysmat_, incvel_, residual_, locsysman_->Trafo(), zeros_, *(vol_surf_flow_bc_maps_));
  }
  else
  {
    // apply Womersley as a Dirichlet BC
    LINALG::ApplyDirichlettoSystem(sysmat_, incvel_, residual_, zeros_, *(vol_surf_flow_bc_maps_));
  }
  return;
}
