/*!-----------------------------------------------------------------------------------------------*
\file scatra_timint_reinitialization.cpp

\brief Time integration scheme for PDE-based reinitialization of level sets
       for further information see
       "Level set transport and PDE-based reinitialization techniques" by Benedikt Schott

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "scatra_timint_implicit.H"

#include "../drt_io/io_gmsh.H"

#include <Teuchos_TimeMonitor.hpp>



/*----------------------------------------------------------------------*
 | contains the time loop for level set reinitialization    schott 05/11|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::TimeLoop_Reinit()
{
  //====================================================================================================
  // REINITIALIZATION                                                           Benedikt Schott 12/2010
  // pde-based reinitialization according to Sussman 1994
  //====================================================================================================


  // set new values:
  double STOP_TOL = 1e-008; // relative tolerance for stop criterion with relative L2-gradient-norm


  // -------------------------------------------------------------------
  // evaluate gradient error before reinitialization
  // -------------------------------------------------------------------
  bool         STOP               = false;
  double       Gradient_Error_old = EvaluateGradientNormError();

  if(Gradient_Error_old > STOP_TOL) ReinitializeInfo();
  else
  {
    cout << "... no reinitialization necessary.\n\n";
    STOP=true;
  }

  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:  + reinitialization time loop");

  while ((step_<stepmax_) and STOP==false)
  {
    PrepareTimeStep();

    // -------------------------------------------------------------------
    //                  solve nonlinear equation
    // -------------------------------------------------------------------
    NonlinearSolve();

    // -------------------------------------------------------------------
    //              check for steady state and update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    CheckSteadyState(Gradient_Error_old,STOP);

  }

} // ScaTraTimIntImpl::TimeLoop_Reinit



/*----------------------------------------------------------------------*
 | add parameters specific for time-integration scheme     schott 05/11 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AddReinitializationParameters(
  ParameterList& params)
{

  // action for elements
  params.set("action","reinitialize_levelset");

  // set type of scalar transport problem
  params.set<int>("scatratype",scatratype_);

  // provide displacement field in case of ALE
  params.set("isale",isale_);

  // factor for scaling the pseudo time step factor
  double fac = 1.0;


  // set params for reinitialization
  params.set("reinit_pseudotimestepfactor", fac*extraparams_->sublist("COMBUSTION PDE REINITIALIZATION").get<double>("PSEUDOTIMESTEP_FACTOR"));
  params.set("reinit_strategy", DRT::INPUT::IntegralValue<INPAR::SCATRA::ReinitializationStrategy>(extraparams_->sublist("COMBUSTION PDE REINITIALIZATION"),"REINIT_METHOD"));
  params.set("reinit_smoothed_sign_type", DRT::INPUT::IntegralValue<INPAR::SCATRA::SmoothedSignType>(extraparams_->sublist("COMBUSTION PDE REINITIALIZATION"),"SMOOTHED_SIGN_TYPE"));
  params.set("reinitswitch",reinitswitch_);
  params.set("reinit_epsilon_bandwidth", extraparams_->sublist("COMBUSTION PDE REINITIALIZATION").get<double>("EPSILON_BANDWIDTH"));
  params.set("reinit_penalty_method", DRT::INPUT::IntegralValue<INPAR::SCATRA::PenaltyMethod>(extraparams_->sublist("COMBUSTION PDE REINITIALIZATION"),"PENALTY_METHOD"));
  params.set("reinit_penalty_interface", extraparams_->sublist("COMBUSTION PDE REINITIALIZATION").get<double>("PENALTY_INTERFACE"));


  params.set("reinit_shock_capturing", DRT::INPUT::IntegralValue<int>(extraparams_->sublist("COMBUSTION PDE REINITIALIZATION"),"SHOCK_CAPTURING"));
  params.set("reinit_shock_capturing_diffusivity", extraparams_->sublist("COMBUSTION PDE REINITIALIZATION").get<double>("SHOCK_CAPTURING_DIFFUSIVITY"));


  // provide velocity field and potentially acceleration/pressure field
  // (export to column map necessary for parallel evaluation)
  AddMultiVectorToParameterList(params,"reinit convective velocity field",convel_);
  AddMultiVectorToParameterList(params,"reinit velocity field",vel_);

  if(timealgo_  == INPAR::SCATRA::timeint_one_step_theta)
  {
    params.set("theta_reinit", 1.0);
    if (discret_->Comm().MyPID() == 0) cout << "WARNING: THETA_REINIT is set to 1.0!!!" << endl;
  }


  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);
  discret_->SetState("phin",phin_);
  discret_->SetState("phistart", phistart_);

  return;
}


/*----------------------------------------------------------------------*
 | contains the assembly process for matrix and rhs        schott 05/11 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AssembleMatAndRHS_Boundary()
{
  //----------------------------------------------------------------------
  // apply Taylor Galerkin Outflow boundary conditions
  //----------------------------------------------------------------------
  if((reinitswitch_ == false) and ( timealgo_ == INPAR::SCATRA::timeint_tg2 or
                                    timealgo_ == INPAR::SCATRA::timeint_tg3    ))
  {
    // evaluate boundary conditions for characteristic galerkin level set transport
    // new parameter list
    ParameterList TaylorGalerkinBoundaryParams;

    // set action for elements
    TaylorGalerkinBoundaryParams.set("action", "levelset_TaylorGalerkin_boundary");
    TaylorGalerkinBoundaryParams.set("incremental solver",incremental_);
    TaylorGalerkinBoundaryParams.set("isale",isale_);

    TaylorGalerkinBoundaryParams.set<int>("scatratype",scatratype_);
    TaylorGalerkinBoundaryParams.set<double>("time_step_size", dta_);
    AddMultiVectorToParameterList(TaylorGalerkinBoundaryParams,"convective velocity field",convel_);
    AddMultiVectorToParameterList(TaylorGalerkinBoundaryParams,"velocity field",vel_);

    discret_->ClearState();

    discret_->SetState("phinp",phinp_);
    discret_->SetState("phin",phin_);

    // add element parameters according to time-integration scheme
    AddSpecificTimeIntegrationParameters(TaylorGalerkinBoundaryParams);


    discret_->EvaluateConditionUsingParentData(TaylorGalerkinBoundaryParams,
                                               sysmat_,
                                               Teuchos::null,
                                               residual_,
                                               Teuchos::null,
                                               Teuchos::null,
                                               "TaylorGalerkinOutflow");

    discret_->ClearState();
  }

  if(reinitswitch_ == true)
  {
    INPAR::SCATRA::TimeIntegrationScheme timealgo_reinit = DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(extraparams_->sublist("COMBUSTION PDE REINITIALIZATION"),"REINIT_TIMEINTEGR");

    if (   timealgo_reinit == INPAR::SCATRA::timeint_tg2
        or timealgo_reinit == INPAR::SCATRA::timeint_tg3
        )
    {
      // evaluate boundary conditions for characteristic galerkin reinitialization method
      // new parameter list
      ParameterList reinitCharacteristicParams;

      // set action for elements
      reinitCharacteristicParams.set("action", "reinitialize_levelset_boundary");
      reinitCharacteristicParams.set("incremental solver",incremental_);
      reinitCharacteristicParams.set("isale",isale_);

      reinitCharacteristicParams.set<int>("scatratype",INPAR::SCATRA::scatratype_condif);
      reinitCharacteristicParams.set<double>("pseudotimestepsize_factor", dta_);

      discret_->ClearState();
      discret_->SetState("phinp",phinp_);
      discret_->SetState("phin",phin_);


      discret_->EvaluateConditionUsingParentData(reinitCharacteristicParams,
                                                 sysmat_,
                                                 Teuchos::null,
                                                 residual_,
                                                 Teuchos::null,
                                                 Teuchos::null,
                                                 "ReinitializationTaylorGalerkin");

      discret_->ClearState();
    }
  }

  return;
} // ScaTraTimIntImpl::AssembleMatAndRHS_Boundary


/*----------------------------------------------------------------------*
 | contains the check for steady state                      schott 05/11|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CheckSteadyState(double&    Gradient_Error_old,
		                                            bool&      STOP)
{
  double current_Gradient_Error = EvaluateGradientNormError();

  INPAR::SCATRA::ReInitialStationaryCheck reinit_stationary_check = DRT::INPUT::IntegralValue<INPAR::SCATRA::ReInitialStationaryCheck>(extraparams_->sublist("COMBUSTION PDE REINITIALIZATION"),"STATIONARY_CHECK");

  if(reinit_stationary_check == INPAR::SCATRA::reinit_stationarycheck_L1normintegrated)
  {
    if(current_Gradient_Error < Gradient_Error_old)
    {
      Update(); // do a further reinitialization step
      Gradient_Error_old = current_Gradient_Error;
      reinitialization_accepted_ = true;

      // -------------------------------------------------------------------
      //                         output of solution
      // -------------------------------------------------------------------
      OutputReinitializationSteps();
    }
    else // the reinitialization step was not successful
    {
      // stop reinitialization
      STOP= true;
      if (myrank_ == 0)
      {
        cout << "reinitialization step not accepted" << endl;
      }
      // set last phi-field;
      phinp_->Update(1.0,*phin_,0.0);
    }
  }
  else
  {
    // check via numsteps
    STOP=false;
    Update();
    reinitialization_accepted_ = true;
    OutputReinitializationSteps();
  }

  return;
}



/*--------------------------------------------------------------------------------------------*
 |  call different time discretizations for the sussman reinitialization equatio schott 04/11 |
 *-------------------------------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntImpl::CallReinitialization()
{

  reinitialization_accepted_= false;

  INPAR::SCATRA::ReinitializationStrategy reinitstrategy = DRT::INPUT::IntegralValue<INPAR::SCATRA::ReinitializationStrategy>(extraparams_->sublist("COMBUSTION PDE REINITIALIZATION"),"REINIT_METHOD");
  INPAR::SCATRA::TimeIntegrationScheme timeintscheme = DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(extraparams_->sublist("COMBUSTION PDE REINITIALIZATION"),"REINIT_TIMEINTEGR");

  // reinitialize Gfunc
  switch(reinitstrategy)
  {
  case INPAR::SCATRA::reinitstrategy_pdebased_characteristic_galerkin:
    if(myrank_==0) std::cout << "COMBUST::Algorithm: reinitialization via PDE-based method" << std::endl;
    if(timeintscheme != INPAR::SCATRA::timeint_tg2) dserror("characteristic galerkin strategy should be used only with tg2 time integration scheme");
    break;
  case INPAR::SCATRA::reinitstrategy_pdebased_linear_convection:
    if(myrank_==0) std::cout << "COMBUST::Algorithm: reinitialization via PDE-based linear transport equation" << std::endl;
    if(timeintscheme != INPAR::SCATRA::timeint_one_step_theta) dserror("linear convection strategy should be used only with OST time integration scheme");
    break;
  case INPAR::SCATRA::reinitstrategy_none:
    if(myrank_==0) std::cout << "No reinitialization chosen" << std::endl;
    break;
  default: dserror("unknown type of reinitialization technique");
  }

  TimeLoop_Reinit();

  return reinitialization_accepted_;
}


/*-----------------------------------------------------------------------------------*
 |  print reinitialization info                                         schott 04/11 |
 *----------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ReinitializeInfo(
    bool printtoscreen,
    bool printtofile)
{
  if (myrank_ == 0) {
    cout << "\n---------------------------------------  REINITIALIZATION SOLVER  ----------------------------\n";
  }


  return;
} // ScaTraImplicitTimeInt::ReinitializeInfo




/*----------------------------------------------------------------------*
 |  calculate error in relative gradient norm               schott 12/10|
 *----------------------------------------------------------------------*/
double SCATRA::ScaTraTimIntImpl::EvaluateGradientNormError()
{


  // create the parameters for the discretization
  ParameterList p;

  // parameters for the elements
  p.set("action","calc_error_reinit");

  //    p.set("scatratype",scatratype_);
  // set type of scalar transport problem
  p.set<int>("scatratype",scatratype_);

  //provide displacement field in case of ALE
  p.set("isale",isale_);

  p.set<double>("L1 integrated gradient error", 0.0);
  p.set<double>("volume", 0.0);

  if (isale_)
    AddMultiVectorToParameterList(p,"dispnp",dispnp_);

  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);

  // get (squared) error values
  discret_->Evaluate(p,null,null,null,null,null);
  discret_->ClearState();


  // get local errors
  double locL1gradienterr   = p.get<double>("L1 integrated gradient error");
  double locvolume          = p.get<double>("volume");

  // initialize global errors
  double L1gradient_err     = 0.0;
  double volume             = 0.0;
  double rel_gradient_err   = 0.0;

  // sum over all processors
  discret_->Comm().SumAll(&locL1gradienterr,&L1gradient_err,1);			// sum over processors, each list (list of a processor) has length 1
  discret_->Comm().SumAll(&locvolume,&volume,1);


  // relative gradient error
  if(fabs(volume) > 1e-014) rel_gradient_err = L1gradient_err / volume;
  else dserror("volume is smaller than 1e-14");


  if (myrank_ == 0)
  {
    printf("\nConvergence check for reinitialization:\n");
    printf("absolute gradient error || (||grad(phi)||-1.0) ||_L1(Omega) %15.8e\n"
        "volume                                                      %15.8e\n"
        "relative gradient error                                     %15.8e\n\n",
        L1gradient_err,volume,rel_gradient_err);
  }


  return rel_gradient_err;
}

/*----------------------------------------------------------------------*
 | set phi vector due to reinitialization                  schott 05/11 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetPhiReinit(Teuchos::RCP<Epetra_Vector> phi)
{
  if (phi != Teuchos::null)
  {
    phinp_->Update(1.0,*phi,0.0);
    phin_->Update(1.0,*phi,0.0);
    phistart_->Update(1.0,*phi,0.0);
  }
  else
    dserror("vector phi does not exist");
  return;
}

/*----------------------------------------------------------------------*
 | set phi vector due to reinitialization                  schott 05/11 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetPhinp(Teuchos::RCP<Epetra_Vector> phinp)
{
  if (phinp != Teuchos::null)
  {
    phinp_->Update(1.0,*phinp,0.0);
  }
  else
    dserror("vector phi does not exist");
  return;
} // ScaTraTimIntImpl::SetPhinp

/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                        schott 01/11|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputReinit(const int globalstep, const double globaltime)
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + output of reinitialized solution");

  OutputToGmshReinit(globalstep, globaltime, step_, time_);

  return;
} // ScaTraTimIntImpl::OutputReinit



/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                        schott 01/11|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputReinitializationSteps()
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + output of reinitialized solution");

  OutputToGmshReinitializationSteps(step_, time_);

  return;
} // ScaTraTimIntImpl::Output


/*----------------------------------------------------------------------*
 | write state vectors to Gmsh postprocessing files        schott  12/12|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputToGmshReinitializationSteps(
    const int reinit_step,
    const double reinit_time
    ) const
{
  // turn on/off screen output for writing process of Gmsh postprocessing file
  const bool screen_out = true;

  // create Gmsh postprocessing file
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_scalar_reinit_step", reinit_step, 1000, screen_out, discret_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinp_Reinit_Step \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(discret_,phinp_,gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }

  gmshfilecontent.close();
  if (screen_out) std::cout << " done" << endl;
}




/*----------------------------------------------------------------------*
 | write state vectors to Gmsh postprocessing files        schott  12/12|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputToGmshReinit(
    const int globalstep,
    const double globaltime,
    const int step,
    const double time
    ) const
{
  // turn on/off screen output for writing process of Gmsh postprocessing file
  const bool screen_out = true;

  // create Gmsh postprocessing file
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_scalar_reinit", globalstep, 1000, screen_out, discret_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
      // add 'View' to Gmsh postprocessing file
      gmshfilecontent << "View \" " << "Phinp_Reinit_start\" {" << endl;
      // draw scalar field 'phin0_Reinit_Reference' for every element
      IO::GMSH::ScalarFieldToGmsh(discret_,phistart_,gmshfilecontent);
      gmshfilecontent << "};" << endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinp_Reinit_final \" {" << endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(discret_,phinp_,gmshfilecontent);
    gmshfilecontent << "};" << endl;
  }

  gmshfilecontent.close();
  if (screen_out) std::cout << " done" << endl;
}



#endif /* CCADISCRET       */
