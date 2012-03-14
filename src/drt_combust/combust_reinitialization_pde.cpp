/*!-----------------------------------------------------------------------------------------------*
\file combust_reinitialization_pde.cpp

\brief
<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

/*
#include "../drt_io/io_gmsh.H"

#include <Teuchos_TimeMonitor.hpp>
*/
#include "combust_reinitialization_pde.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"



#include "../drt_io/io_gmsh.H"
#include <Teuchos_TimeMonitor.hpp>


#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_sparseoperator.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Epetra_SerialDenseVector.h>


/*----------------------------------------------------------------------*
 | constructor                                              schott 03/12|
 *----------------------------------------------------------------------*/
COMBUST::ReinitializationPDE::ReinitializationPDE()
{
  Teuchos::ParameterList combustdyn = DRT::Problem::Instance()->CombustionDynamicParams();

  // make a copy (inside an rcp) containing also all sublists
  RCP<ParameterList> combustdynreinit_ = rcp(new ParameterList(combustdyn));

  // -------------------------------------------------------------------
  // overrule certain parameters
  // -------------------------------------------------------------------
  combustdynreinit_->set<double>   ("TIMESTEP"    ,combustdyn.sublist("COMBUSTION PDE REINITIALIZATION").get<double>("PSEUDOTIMESTEP_FACTOR"));
  // maximum simulation time
  combustdynreinit_->set<double>   ("MAXTIME"     ,combustdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  combustdynreinit_->set<int>      ("NUMSTEP"     ,combustdyn.sublist("COMBUSTION PDE REINITIALIZATION").get<int>("NUMPSEUDOSTEPS"));
  // restart
  combustdynreinit_->set           ("RESTARTEVRY" ,combustdyn.get<int>("RESTARTEVRY"));
  // solution output
  combustdynreinit_->set           ("UPRES"       ,combustdyn.get<int>("UPRES"));

  reinit_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(
      *combustdynreinit_,
      false, // is_ale
      0,     // use first scatra discretization
      DRT::Problem::Instance()->ScalarTransportFluidSolverParams(),
      true   //
      ));}

/*----------------------------------------------------------------------*
 | destructor                                               schott 03/12|
 *----------------------------------------------------------------------*/
 COMBUST::ReinitializationPDE::~ReinitializationPDE(){}


/*----------------------------------------------------------------------*
 | contains the time loop for level set reinitialization    schott 03/12|
 *----------------------------------------------------------------------*/
void COMBUST::ReinitializationPDE::TimeLoop_Reinit()
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
  double       Gradient_Error_old = reinit_->ScaTraField().EvaluateGradientNormError();

  if(Gradient_Error_old > STOP_TOL) ReinitializeInfo();
  else
  {
    cout << "... no reinitialization necessary.\n\n";
    STOP=true;
  }

  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:  + reinitialization time loop");

  while ((reinit_->ScaTraField().Step()<stepmax_) and STOP==false)
  {
    reinit_->ScaTraField().PrepareTimeStep();

    // -------------------------------------------------------------------
    //                  solve nonlinear equation
    // -------------------------------------------------------------------
    reinit_->ScaTraField().Solve();

    // -------------------------------------------------------------------
    //              check for steady state and update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    reinit_->ScaTraField().CheckSteadyState(Gradient_Error_old,STOP);

  }

} // ScaTraTimIntImpl::TimeLoop_Reinit


/*----------------------------------------------------------------------*
 | contains the check for steady state                      schott 05/11|
 *----------------------------------------------------------------------*/
/*void COMBUST::ReinitializationPDE::CheckSteadyState(double&    Gradient_Error_old,
		                                            bool&      STOP)
{
  double current_Gradient_Error = reinit_->EvaluateGradientNormError();

  INPAR::SCATRA::ReInitialStationaryCheck reinit_stationary_check = DRT::INPUT::IntegralValue<INPAR::SCATRA::ReInitialStationaryCheck>(extraparams_->sublist("COMBUSTION PDE REINITIALIZATION"),"STATIONARY_CHECK");

  if(reinit_stationary_check == INPAR::SCATRA::reinit_stationarycheck_L1normintegrated)
  {
    if(current_Gradient_Error < Gradient_Error_old)
    {
      reinit_->Update(); // do a further reinitialization step
      Gradient_Error_old = current_Gradient_Error;
      reinitialization_accepted_ = true;

      // -------------------------------------------------------------------
      //                         output of solution
      // -------------------------------------------------------------------
      reinit_->OutputReinitializationSteps();
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
    reinit_->Update();
    reinitialization_accepted_ = true;
    reinit_->OutputReinitializationSteps();
  }

  return;
}
*/


/*--------------------------------------------------------------------------------------------*
 |  call different time discretizations for the sussman reinitialization equatio schott 04/11 |
 *-------------------------------------------------------------------------------------------*/
bool COMBUST::ReinitializationPDE::CallReinitialization()
{

  bool reinitialization_accepted_= false;

  INPAR::SCATRA::ReinitializationStrategy reinitstrategy = DRT::INPUT::IntegralValue<INPAR::SCATRA::ReinitializationStrategy>(combustdynreinit_->sublist("COMBUSTION PDE REINITIALIZATION"),"REINIT_METHOD");
  INPAR::SCATRA::TimeIntegrationScheme timeintscheme = DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(combustdynreinit_->sublist("COMBUSTION PDE REINITIALIZATION"),"REINIT_TIMEINTEGR");

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
void COMBUST::ReinitializationPDE::ReinitializeInfo()
{
  if (myrank_ == 0) {
    cout << "\n---------------------------------------  REINITIALIZATION SOLVER  ----------------------------\n";
  }


  return;
} // ScaTraImplicitTimeInt::ReinitializeInfo




/*----------------------------------------------------------------------*
 |  calculate error in relative gradient norm               schott 12/10|
 *----------------------------------------------------------------------*/
/*double COMBUST::ReinitializationPDE::EvaluateGradientNormError()
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
}*/


