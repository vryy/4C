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

#include "scatra_timint_implicit.H"
#include "scatra_ele_action.H"

#include "../drt_io/io_gmsh.H"

#include <Teuchos_TimeMonitor.hpp>


/*==========================================================================*
 |                                                                          |
 | protected:                                                               |
 |                                                                          |
 *==========================================================================*/

/*----------------------------------------------------------------------*
 | add parameters specific for time-integration scheme     schott 05/11 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AddReinitializationParameters(
    Teuchos::ParameterList& params)
{

  // action for elements
  params.set<int>("action",SCATRA::reinitialize_levelset);

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
//    if (discret_->Comm().MyPID() == 0) cout << "WARNING: THETA_REINIT is set to 1.0!!!" << endl;
  }


  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp",phinp_);
  discret_->SetState("phin",phin_);
  discret_->SetState("phistart", phistart_);

  return;
} // ScaTraTimIntImpl::AddReinitializationParameters


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
    Teuchos::ParameterList TaylorGalerkinBoundaryParams;

    // set action for elements
    TaylorGalerkinBoundaryParams.set<int>("action", SCATRA::bd_calc_TG_outflow);
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
      Teuchos::ParameterList reinitCharacteristicParams;

      // set action for elements
      reinitCharacteristicParams.set<int>("action",SCATRA::bd_reinitialize_levelset);
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

