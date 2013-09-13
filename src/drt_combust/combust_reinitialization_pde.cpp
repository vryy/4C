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


#include <Teuchos_TimeMonitor.hpp>
#include "combust_reinitialization_pde.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_pstream.H"
#include "../drt_scatra/scatra_ele_action.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 | constructor                                              schott 03/12|
 *----------------------------------------------------------------------*/
COMBUST::ReinitializationPDE::ReinitializationPDE(const Epetra_Comm& comm)
{
  myrank_ = comm.MyPID();

  Teuchos::ParameterList combustdyn = DRT::Problem::Instance()->CombustionDynamicParams();

  // make a copy (inside an Teuchos::rcp) containing also all sublists
  combustdynreinit_ = Teuchos::rcp(new Teuchos::ParameterList(combustdyn));

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

  // get linear solver id from SCALAR TRANSPORT DYNAMIC
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for COMBUST problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

  reinit_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(
      *combustdynreinit_,
      false, // is_ale
      "scatra",
      DRT::Problem::Instance()->SolverParams(linsolvernumber),
      true   // reinitswitch
      )
  );

  stepmax_ = combustdyn.sublist("COMBUSTION PDE REINITIALIZATION").get<int>("NUMPSEUDOSTEPS");

  reinitialization_accepted_=false;

} // end constructor


/*----------------------------------------------------------------------*
 | get Phinp vector                                         schott 03/12|
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> COMBUST::ReinitializationPDE::Phinp()
{
  return ScaTraReinitField().Phinp();
}


/*----------------------------------------------------------------------*
 | destructor                                               schott 03/12|
 *----------------------------------------------------------------------*/
 COMBUST::ReinitializationPDE::~ReinitializationPDE(){}




/*----------------------------------------------------------------------*
 | contains the check for steady state                      schott 05/11|
 *----------------------------------------------------------------------*/
void COMBUST::ReinitializationPDE::CheckSteadyState(double&    Gradient_Error_old,
                                                    bool&      STOP)
{
  double current_Gradient_Error = EvaluateGradientNormError();

  INPAR::SCATRA::ReInitialStationaryCheck reinit_stationary_check
                = DRT::INPUT::IntegralValue<INPAR::SCATRA::ReInitialStationaryCheck>(combustdynreinit_->sublist("COMBUSTION PDE REINITIALIZATION"),"STATIONARY_CHECK");

  if(reinit_stationary_check == INPAR::SCATRA::reinit_stationarycheck_L1normintegrated)
  {
    if(current_Gradient_Error < Gradient_Error_old)
    {
      ScaTraReinitField().Update(); // do a further reinitialization step
      Gradient_Error_old = current_Gradient_Error;
      reinitialization_accepted_ = true;

      // -------------------------------------------------------------------
      //                         output of solution
      // -------------------------------------------------------------------
//      OutputReinitializationSteps();
    }
    else // the reinitialization step was not successful
    {
      // stop reinitialization
      STOP= true;
      if (myrank_ == 0)
      {
        std::cout << "reinitialization step not accepted" << std::endl;
      }
      // set last phi-field;
      ScaTraReinitField().Phinp()->Update(1.0,*ScaTraReinitField().Phin(),0.0);
    }
  }
  else if(reinit_stationary_check == INPAR::SCATRA::reinit_stationarycheck_numsteps)
  {
    // check via numsteps
    STOP=false;
    ScaTraReinitField().Update();
    reinitialization_accepted_ = true;
//    OutputReinitializationSteps();
  }
  else dserror("no valid stationary check");

  return;

} //COMBUST::ReinitializationPDE::CheckSteadyState



/*--------------------------------------------------------------------------------------------*
 |  call different time discretizations for the sussman reinitialization equatio schott 04/11 |
 *-------------------------------------------------------------------------------------------*/
bool COMBUST::ReinitializationPDE::CallReinitialization(Teuchos::RCP<Epetra_Vector> phinp)
{

  reinitialization_accepted_= false;


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
    break;
  }

  // reset the ScatraFieldReinit
  ScaTraReinitField().SetTimeStep(0.0, 0);

  // set the start phinp and phin field and phistart field
  SetPhiReinit(phinp);

  TimeLoop_Reinit();

  // Update the origial phinp vector in the level set scatra field
  if(reinitialization_accepted_)
  {
      phinp->Update(1.0,*(ScaTraReinitField().Phinp()),0.0);
  }

  return reinitialization_accepted_;

} //COMBUST::ReinitializationPDE::CallReinitialization


SCATRA::ScaTraTimIntImpl& COMBUST::ReinitializationPDE::ScaTraReinitField()
{
  return reinit_->ScaTraField();
}


/*----------------------------------------------------------------------*
 |  calculate error in relative gradient norm               schott 12/10|
 *----------------------------------------------------------------------*/
double COMBUST::ReinitializationPDE::EvaluateGradientNormError()
{

  Teuchos::RCP<DRT::Discretization> discret = ScaTraReinitField().Discretization();

  // create the parameters for the discretization
  Teuchos::ParameterList p;

  // parameters for the elements
  p.set<int>("action",SCATRA::calc_error_reinit);

  //    p.set("scatratype",scatratype_);
  // set type of scalar transport problem
  p.set<int>("scatratype",INPAR::SCATRA::scatratype_levelset);

  //provide displacement field in case of ALE
  p.set("isale",false);

  p.set<double>("L1 integrated gradient error", 0.0);
  p.set<double>("volume", 0.0);

  // set vector values needed by elements
  discret->ClearState();
  discret->SetState("phinp",ScaTraReinitField().Phinp());

  // get (squared) error values
  discret->Evaluate(p,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  discret->ClearState();


  // get local errors
  double locL1gradienterr   = p.get<double>("L1 integrated gradient error");
  double locvolume          = p.get<double>("volume");

  // initialize global errors
  double L1gradient_err     = 0.0;
  double volume             = 0.0;
  double rel_gradient_err   = 0.0;

  // sum over all processors
  discret->Comm().SumAll(&locL1gradienterr,&L1gradient_err,1);     // sum over processors, each list (list of a processor) has length 1
  discret->Comm().SumAll(&locvolume,&volume,1);


  // relative gradient error
  if(fabs(volume) > 1e-014) rel_gradient_err = L1gradient_err / volume;
  else dserror("volume is smaller than 1e-14");


  if (discret->Comm().MyPID() == 0)
  {
    printf("\nConvergence check for reinitialization:\n");
    printf("absolute gradient error || (||grad(phi)||-1.0) ||_L1(Omega) %15.8e\n"
        "volume                                                      %15.8e\n"
        "relative gradient error                                     %15.8e\n\n",
        L1gradient_err,volume,rel_gradient_err);
  }


  return rel_gradient_err;

} //COMBUST::ReinitializationPDE::EvaluateGradientNormError


/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                        schott 01/11|
 *----------------------------------------------------------------------*/
void COMBUST::ReinitializationPDE::OutputReinit(const int globalstep, const double globaltime)
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + output of reinitialized solution");

  OutputToGmshReinit(globalstep, globaltime, ScaTraReinitField().Step(), ScaTraReinitField().Time());

  return;
} // COMBUST::ReinitializationPDE::OutputReinit



/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                        schott 01/11|
 *----------------------------------------------------------------------*/
void COMBUST::ReinitializationPDE::OutputReinitializationSteps()
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + output of reinitialized solution");

  OutputToGmshReinitializationSteps(ScaTraReinitField().Step(), ScaTraReinitField().Time());

  return;
} //COMBUST::ReinitializationPDE::OutputReinitializationSteps


/*----------------------------------------------------------------------*
 | write state vectors to Gmsh postprocessing files        schott  12/12|
 *----------------------------------------------------------------------*/
void COMBUST::ReinitializationPDE::OutputToGmshReinitializationSteps(
    const int reinit_step,
    const double reinit_time
    )
{
  // turn on/off screen output for writing process of Gmsh postprocessing file
  const bool screen_out = true;

  Teuchos::RCP<DRT::Discretization> dis = ScaTraReinitField().Discretization();

  // create Gmsh postprocessing file
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_scalar_reinit_step", reinit_step, 1000, screen_out, dis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinp_Reinit_Step \" {" << std::endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(dis,ScaTraReinitField().Phinp(),gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  gmshfilecontent.close();
  if (screen_out) std::cout << " done" << std::endl;
} //COMBUST::ReinitializationPDE::OutputToGmshReinitializationSteps




/*----------------------------------------------------------------------*
 | write state vectors to Gmsh postprocessing files        schott  12/12|
 *----------------------------------------------------------------------*/
void COMBUST::ReinitializationPDE::OutputToGmshReinit(
    const int globalstep,
    const double globaltime,
    const int step,
    const double time
    )
{
  // turn on/off screen output for writing process of Gmsh postprocessing file
  const bool screen_out = true;

  Teuchos::RCP<DRT::Discretization> dis = ScaTraReinitField().Discretization();

  // create Gmsh postprocessing file
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("solution_field_scalar_reinit", globalstep, 1000, screen_out, dis->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  {
      // add 'View' to Gmsh postprocessing file
      gmshfilecontent << "View \" " << "Phinp_Reinit_start\" {" << std::endl;
      // draw scalar field 'phin0_Reinit_Reference' for every element
      IO::GMSH::ScalarFieldToGmsh(dis,ScaTraReinitField().PhiReinitStart(),gmshfilecontent);
      gmshfilecontent << "};" << std::endl;
  }
  {
    // add 'View' to Gmsh postprocessing file
    gmshfilecontent << "View \" " << "Phinp_Reinit_final \" {" << std::endl;
    // draw scalar field 'Phinp' for every element
    IO::GMSH::ScalarFieldToGmsh(dis,ScaTraReinitField().Phinp(),gmshfilecontent);
    gmshfilecontent << "};" << std::endl;
  }

  gmshfilecontent.close();
  if (screen_out) std::cout << " done" << std::endl;

} //COMBUST::ReinitializationPDE::OutputToGmshReinit





/*-----------------------------------------------------------------------------------*
 |  print reinitialization info                                         schott 04/11 |
 *----------------------------------------------------------------------------------*/
void COMBUST::ReinitializationPDE::ReinitializeInfo()
{
  if (myrank_ == 0) {
    std::cout << "\n---------------------------------------  REINITIALIZATION SOLVER  ----------------------------\n";
  }


  return;
} // COMBUST::ReinitializationPDE::ReinitializeInfo()



/*----------------------------------------------------------------------*
  | set phi vector due to reinitialization                  schott 05/11 |
 *----------------------------------------------------------------------*/
void COMBUST::ReinitializationPDE::SetPhiReinit(Teuchos::RCP<Epetra_Vector> phi)
{
  if (phi != Teuchos::null)
  {
    ScaTraReinitField().Phinp()->Update(1.0,*phi,0.0);
    ScaTraReinitField().Phin()->Update(1.0,*phi,0.0);
    ScaTraReinitField().PhiReinitStart()->Update(1.0,*phi,0.0);
  }
  else
    dserror("vector phi does not exist");
  return;
} //COMBUST::ReinitializationPDE::SetPhiReinit


/*----------------------------------------------------------------------------------*
  | redistribute reinitializer after redistribution of scatra field rasthofer 09/11 |
 *----------------------------------------------------------------------------------*/
void COMBUST::ReinitializationPDE::Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph)
{
  if(myrank_==0)
    IO::cout << "Redistributing PDE Reinitializer                      ... " << IO::endl;

  ScaTraReinitField().Redistribute(nodegraph);

  if(myrank_==0)
    IO::cout << "done" << IO::endl;
}


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
  double       Gradient_Error_old = EvaluateGradientNormError();

  if(Gradient_Error_old > STOP_TOL) ReinitializeInfo();
  else
  {
    std::cout << "... no reinitialization necessary.\n\n";
    STOP=true;
  }

  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:  + reinitialization time loop");

  while ((ScaTraReinitField().Step()<stepmax_) and STOP==false)
  {
    ScaTraReinitField().PrepareTimeStep();

    // -------------------------------------------------------------------
    //                  solve nonlinear equation
    // -------------------------------------------------------------------
    ScaTraReinitField().Solve();

    // -------------------------------------------------------------------
    //              check for steady state and update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    CheckSteadyState(Gradient_Error_old,STOP);

  }

} //COMBUST::ReinitializationPDE::TimeLoop_Reinit
