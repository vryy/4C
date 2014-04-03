/*----------------------------------------------------------------------*/
/*!
\file fluid_timint_red.cpp
\brief TimIntRedModels

<pre>
Maintainer: Mahmoud Ismail
            ismail@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_timint_red.H"
#include "fluid_coupling_red_models.H"
#include "fluid_volumetric_surfaceFlow_condition.H"
#include "fluid_windkessel_optimization.H"
#include "../drt_art_net/artnetexplicitintegration.H"
#include "fluid_meshtying.H"
#include "../drt_lib/drt_locsys.H"
#include "../drt_io/io.H"


/*----------------------------------------------------------------------*
 |  Constructor (public)                                       bk 11/13 |
 *----------------------------------------------------------------------*/
FLD::TimIntRedModels::TimIntRedModels(
        const Teuchos::RCP<DRT::Discretization>&      actdis,
        const Teuchos::RCP<LINALG::Solver>&           solver,
        const Teuchos::RCP<Teuchos::ParameterList>&   params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output,
        bool                                          alefluid /*= false*/)
    : FluidImplicitTimeInt(actdis,solver,params,output,alefluid),
      impedancebc_(Teuchos::null),
      Wk_optimization_(Teuchos::null),
      vol_surf_flow_bc_(Teuchos::null),
      traction_vel_comp_adder_bc_(Teuchos::null),
      coupled3D_redDbc_art_(Teuchos::null),
      ART_exp_timeInt_(Teuchos::null),
      coupled3D_redDbc_airways_(Teuchos::null),
      airway_imp_timeInt_(Teuchos::null),
      vol_surf_flow_bcmaps_(Teuchos::null),
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
#ifdef D_ALE_BFLOW
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispn_);
  }
#endif

  vol_surf_flow_bc_     = Teuchos::rcp(new UTILS::FluidVolumetricSurfaceFlowWrapper(discret_, *output_, dta_) );

  // evaluate the map of the womersley bcs
  vol_surf_flow_bc_ -> EvaluateMapExtractor(vol_flow_rates_bc_extractor_);
  vol_surf_flow_bc_ -> EvaluateCondMap(vol_surf_flow_bcmaps_);
//  std::cout<< "velnp_first: " << *velnp_ << std::endl;
//  std::cout << "time_: " << time_ << std::endl;
  // Evaluate the womersley velocities
  vol_surf_flow_bc_->EvaluateVelocities(velnp_,time_);
//  std::cout<< "velnp_last: " << *velnp_ << std::endl;
  // -------------------------------------------------------------------
  // Initialize the reduced models
  // -------------------------------------------------------------------

  strong_redD_3d_coupling_ = false;
  if (params_->get<std::string>("Strong 3D_redD coupling","no") == "yes")   strong_redD_3d_coupling_ = true;

  {
    ART_exp_timeInt_ = dyn_art_net_drt(true);
    // Check if one-dimensional artery network problem exist
    if (ART_exp_timeInt_ != Teuchos::null)
    {
      IO::DiscretizationWriter output_redD(ART_exp_timeInt_->Discretization());
      discret_->ClearState();
      discret_->SetState("velaf", zeros_);
      if (alefluid_)
      {
        discret_->SetState("dispnp", dispnp_);
      }
      coupled3D_redDbc_art_=   Teuchos::rcp(new  UTILS::Fluid_couplingWrapper<ART::ArtNetExplicitTimeInt>
                                   ( discret_,
                                     ART_exp_timeInt_->Discretization(),
                                     ART_exp_timeInt_,
                                     output_redD,
                                     dta_,
                                     ART_exp_timeInt_->Dt()));
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
      coupled3D_redDbc_airways_ =   Teuchos::rcp(new  UTILS::Fluid_couplingWrapper<AIRWAY::RedAirwayImplicitTimeInt>
                                   ( discret_,
                                     airway_imp_timeInt_->Discretization(),
                                     airway_imp_timeInt_,
                                     output_redD,
                                     dta_,
                                     airway_imp_timeInt_->Dt()));

    }


    zeros_->PutScalar(0.0); // just in case of change
  }

  traction_vel_comp_adder_bc_ = Teuchos::rcp(new UTILS::TotalTractionCorrector(discret_, *output_, dta_) );

#ifdef D_ALE_BFLOW
  if (alefluid_)
  {
    discret_->ClearState();
    discret_->SetState("dispnp", dispnp_);
  }
#endif // D_ALE_BFLOW
  // construct impedance bc wrapper
  impedancebc_      = Teuchos::rcp(new UTILS::FluidImpedanceWrapper(discret_, *output_, dta_) );

  Wk_optimization_  = Teuchos::rcp(new UTILS::FluidWkOptimizationWrapper(discret_,
                                                                *output_,
                                                                impedancebc_,
                                                                dta_) );

  // ------------------------------------------------------------------------------
  // Check, if features are used with the locsys manager that are not supported,
  // or better, not implemented yet.
  // ------------------------------------------------------------------------------
  if (locsysman_ != Teuchos::null) {

    // Models
    if ((ART_exp_timeInt_ != Teuchos::null) or (airway_imp_timeInt_ != Teuchos::null)) {
      dserror("No problem types involving airways are supported for use with locsys conditions!");
    }
  }

  return;
}


/*----------------------------------------------------------------------*
| Destructor dtor (public)                                    bk 11/13 |
*----------------------------------------------------------------------*/
FLD::TimIntRedModels::~TimIntRedModels()
{
  return;
}

/*----------------------------------------------------------------------*
 | evaluate special boundary conditions                        bk 12/13 |
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::DoProblemSpecificBoundaryConditions()
{

#ifdef D_ALE_BFLOW
    if (alefluid_)
    {
      discret_->SetState("dispnp", dispnp_);
    }
#endif // D_ALE_BFLOW

  // Check if one-dimensional artery network problem exist
  if (ART_exp_timeInt_ != Teuchos::null)
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
  vol_surf_flow_bc_->EvaluateVelocities(velnp_,time_);

  return;
}

/*----------------------------------------------------------------------*
| Update3DToReduced in AssembleMatAndRHS                       bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntRedModels::Update3DToReducedMatAndRHS()
{

  // update impedance boundary condition
  impedancebc_->UpdateResidual(residual_);

  discret_->ClearState();

  discret_->SetState("velaf",velnp_);
  discret_->SetState("hist",hist_);

#ifdef D_ALE_BFLOW
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
  }
#endif // D_ALE_BFLOW

  // Check if one-dimensional artery network problem exist
  if (ART_exp_timeInt_ != Teuchos::null)
  {
    if (strong_redD_3d_coupling_)
    {
      coupled3D_redDbc_art_->FlowRateCalculation(time_,dta_);
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
      coupled3D_redDbc_airways_->FlowRateCalculation(time_,dta_);
      coupled3D_redDbc_airways_->ApplyBoundaryConditions(time_, dta_, theta_);
    }
    coupled3D_redDbc_airways_->UpdateResidual(residual_);
  }

  //----------------------------------------------------------------------
  // add the traction velocity component
  //----------------------------------------------------------------------

  traction_vel_comp_adder_bc_->EvaluateVelocities(velnp_,time_,theta_,dta_);
  traction_vel_comp_adder_bc_->UpdateResidual(residual_);

  discret_->ClearState();
  return;
}

/*----------------------------------------------------------------------*
| call Update3DToReducedMatAndRHS                              bk 11/13 |
*----------------------------------------------------------------------*/
void FLD::TimIntRedModels::SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams)
{
  //these are the only routines that have to be called in AssembleMatAndRHS
  //before Evaluate in the RedModels case
  Update3DToReducedMatAndRHS();

  return;
}

/*----------------------------------------------------------------------*
 | output of solution vector of ReducedD problem to binio   ismail 01/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::OutputReducedD()
{
  // output of solution
  if (step_%upres_ == 0)
  {
    // write reduced model problem
    // Check if one-dimensional artery network problem exist
    if (ART_exp_timeInt_ != Teuchos::null)
    {
      Teuchos::RCP<Teuchos::ParameterList> redD_export_params;
      redD_export_params = Teuchos::rcp(new Teuchos::ParameterList());

      redD_export_params->set<int>("step",step_);
      redD_export_params->set<int>("upres",upres_);
      redD_export_params->set<int>("uprestart",uprestart_);
      redD_export_params->set<double>("time",time_);

      ART_exp_timeInt_->Output(true, redD_export_params);
    }

    // Check if one-dimensional artery network problem exist
    if (airway_imp_timeInt_ != Teuchos::null)
    {
      Teuchos::RCP<Teuchos::ParameterList> redD_export_params;
      redD_export_params = Teuchos::rcp(new Teuchos::ParameterList());

      redD_export_params->set<int>("step",step_);
      redD_export_params->set<int>("upres",upres_);
      redD_export_params->set<int>("uprestart",uprestart_);
      redD_export_params->set<double>("time",time_);

      airway_imp_timeInt_->Output(true, redD_export_params);
    }
  }
  return;
}//FLD::TimIntRedModels::OutputReducedD

/*----------------------------------------------------------------------*
 | read some additional data in restart                         bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::ReadRestart(int step)
{
  FluidImplicitTimeInt::ReadRestart(step);

  IO::DiscretizationReader reader(discret_,step);
  // also read impedance bc information if required
  // Note: this method acts only if there is an impedance BC
  impedancebc_->ReadRestart(reader);

  Wk_optimization_->ReadRestart(reader);

  vol_surf_flow_bc_->ReadRestart(reader);

  traction_vel_comp_adder_bc_->ReadRestart(reader);

  // Read restart of one-dimensional arterial network
  if (ART_exp_timeInt_ != Teuchos::null)
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
  if (ART_exp_timeInt_ != Teuchos::null)
  {
    ART_exp_timeInt_->ReadRestart(step,true);
  }

  // Check if one-dimensional artery network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    airway_imp_timeInt_->ReadRestart(step,true);
  }
}//FLD::TimIntRedModels::ReadRestartReadRestart(int step)

/*----------------------------------------------------------------------*
 | do some additional steps in TimeUpdate                       bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::TimeUpdate()
{
  FluidImplicitTimeInt::TimeUpdate();


  // -------------------------------------------------------------------
  // treat impedance BC
  // note: these methods return without action, if the problem does not
  //       have impedance boundary conditions
  // -------------------------------------------------------------------
  discret_->ClearState();
  discret_->SetState("velaf",velnp_);
  discret_->SetState("hist",hist_);

#ifdef D_ALE_BFLOW
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispn_);
  }
#endif //D_ALE_BFLOW

  impedancebc_->FlowRateCalculation(time_,dta_);
  impedancebc_->OutflowBoundary(time_,dta_,theta_);

  // get the parameters needed to be optimized
  Teuchos::ParameterList WkOpt_params;
  WkOpt_params.set<double> ("total time", time_);
  WkOpt_params.set<double> ("time step size", dta_);
  impedancebc_->getResultsOfAPeriod(WkOpt_params);

  // update wind kessel optimization condition
  Wk_optimization_->Solve(WkOpt_params);

  // -------------------------------------------------------------------
  // treat the 3D-to-reduced_D coupling condition
  // note: these methods return without action, if the problem does not
  //       have any coupling boundary conditions
  // -------------------------------------------------------------------


#ifdef D_ALE_BFLOW
  if (alefluid_)
  {
    discret_->SetState("dispnp", dispnp_);
  }
#endif // D_ALE_BFLOW

  // Check if one-dimensional artery network problem exist
  if (ART_exp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_art_->FlowRateCalculation(time_,dta_);
    coupled3D_redDbc_art_->ApplyBoundaryConditions(time_, dta_, theta_);
    //    coupled3D_redDbc_art_->TimeUpdate();
  }


  // Check if one-dimensional artery network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_airways_->FlowRateCalculation(time_,dta_);
    coupled3D_redDbc_airways_->ApplyBoundaryConditions(time_, dta_, theta_);
    //    coupled3D_redDbc_airways_->TimeUpdate();
  }

  discret_->ClearState();

  // update the 3D-to-reduce_D coupling condition
  // update the 3D-to-reduced_D coupling data
  // Check if one-dimensional artery network problem exist
  if (ART_exp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_art_->TimeUpdate();
  }
  // update the 3D-to-reduced_D coupling data
  // Check if one-dimensional artery network problem exist
  if (airway_imp_timeInt_ != Teuchos::null)
  {
    coupled3D_redDbc_airways_->TimeUpdate();
  }
  return;
}

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
  if(vol_surf_flow_bcmaps_->NumGlobalElements() != 0)
  {
    meshtying_->CheckOverlappingBC(vol_surf_flow_bcmaps_);
    meshtying_->DirichletOnMaster(vol_surf_flow_bcmaps_);

    if(myrank_==0)
      std::cout << "Think about your flow rate defined on the interal interface!!" << std::endl << std::endl;
      dserror("Think about your flow rate defined on the interal interface!!\n"
              "It is working qualitatively!!");
  }
  return;
}

/*----------------------------------------------------------------------*
 | do some additional steps in UpdateIterIncrementally          bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::UpdateIterIncrementally(
  Teuchos::RCP<const Epetra_Vector> vel)  //!< input residual velocities

{
  FluidImplicitTimeInt::UpdateIterIncrementally(vel);
  // set the new solution we just got
  if (vel != Teuchos::null)
  {

    Teuchos::RCP<Epetra_Vector> aux = LINALG::CreateVector(
        *(discret_->DofRowMap(0)), true);

    *aux=*velnp_;
    //
    vol_flow_rates_bc_extractor_->InsertVolumetricSurfaceFlowCondVector(
        vol_flow_rates_bc_extractor_->ExtractVolumetricSurfaceFlowCondVector(
            velnp_), aux);
    *velnp_=*aux;

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
  if (step_%upres_ == 0)
  {
    vol_surf_flow_bc_->Output(*output_);
    traction_vel_comp_adder_bc_->Output(*output_);

    if (uprestart_ != 0 && step_%uprestart_ == 0) //add restart data
    {
      // also write impedance bc information if required
      // Note: this method acts only if there is an impedance BC
      impedancebc_->WriteRestart(*output_);

      Wk_optimization_->WriteRestart(*output_);
      // write reduced model problem
      // Check if one-dimensional artery network problem exist
      if (ART_exp_timeInt_ != Teuchos::null)
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
  else if (uprestart_ > 0 && step_%uprestart_ == 0)
  {
    // also write impedance bc information if required
    // Note: this method acts only if there is an impedance BC
    impedancebc_->WriteRestart(*output_);

    Wk_optimization_->WriteRestart(*output_);
    // write reduced model problem
    // Check if one-dimensional artery network problem exist
    if (ART_exp_timeInt_ != Teuchos::null)
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
} // TimIntRedModels::Output

/*----------------------------------------------------------------------*
 | read some additional data in restart                         bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::InsertVolumetricSurfaceFlowCondVector(Teuchos::RCP<Epetra_Vector> vel ,Teuchos::RCP<Epetra_Vector> res)
{

  // -------------------------------------------------------------------
  // take surface volumetric flow rate into account
  //    Teuchos::RCP<Epetra_Vector> temp_vec = Teuchos::rcp(new Epetra_Vector(*vol_surf_flow_bcmaps_,true));
  //    vol_surf_flow_bc_->InsertCondVector( *temp_vec , *residual_);
  // -------------------------------------------------------------------
  vol_flow_rates_bc_extractor_->InsertVolumetricSurfaceFlowCondVector(
    vol_flow_rates_bc_extractor_->ExtractVolumetricSurfaceFlowCondVector(vel),
    res);

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

  //necessary here, because some application time integrations add something to the residual
  //before the Neumann loads are added
  residual_->PutScalar(0.0);

  // add impedance Neumann loads
  impedancebc_->UpdateResidual(residual_);

  AVM3AssembleMatAndRHS(eleparams);

  // apply Womersley as a Dirichlet BC
  LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(vol_surf_flow_bcmaps_));

  // get scale-separation matrix
  AVM3GetScaleSeparationMatrix();

  return;
}// TimIntRedModels::AVM3Preparation

/*----------------------------------------------------------------------*
 | RedModels - specific BC in LinearRelaxationSolve            bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::CustomSolve(Teuchos::RCP<Epetra_Vector> relax)
{
  // apply Womersley as a Dirichlet BC
  LINALG::ApplyDirichlettoSystem(incvel_,residual_,relax,*(vol_surf_flow_bcmaps_));

  // apply Womersley as a Dirichlet BC
  sysmat_->ApplyDirichlet(*(vol_surf_flow_bcmaps_));

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
    shapederivatives_->ApplyDirichlet(*(vol_surf_flow_bcmaps_),false);
  }

  return;
}

/*----------------------------------------------------------------------*
 | Apply Womersley bc to system                                 bk 12/13|
 *----------------------------------------------------------------------*/
void FLD::TimIntRedModels::ApplyDirichletToSystem()
{

  FluidImplicitTimeInt::ApplyDirichletToSystem();

  if (LocsysManager() != Teuchos::null) {
    // apply Womersley as a Dirichlet BC
    LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,locsysman_->Trafo(),zeros_,*(vol_surf_flow_bcmaps_));

  }
  else
  {

    // apply Womersley as a Dirichlet BC
    LINALG::ApplyDirichlettoSystem(sysmat_,incvel_,residual_,zeros_,*(vol_surf_flow_bcmaps_));

  }
  return;
}
