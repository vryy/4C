/*----------------------------------------------------------------------*/
/*!
\file xfluid_levelset_coupling_algorithm.cpp

\brief Basis of xfluid-levelset coupling.

\level 3

\maintainer Magnus Winter
            winter@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245

*/
/*----------------------------------------------------------------------*/

#include "../drt_levelset/levelset_algorithm.H"
#include "../drt_levelset/levelset_timint_ost.H"
#include "../drt_scatra/scatra_timint_ost.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fluid/fluid_utils_mapextractor.H" //needed?
#include "../drt_fluid_xfluid/xfluid.H"

#include "../drt_lib/drt_discret_xfem.H"


#include "xfluid_levelset_coupling_algorithm.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
XFLUIDLEVELSET::Algorithm::Algorithm(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn,
    const Teuchos::ParameterList& solverparams
    )
:  ScaTraFluidCouplingAlgorithm(comm,prbdyn,false,"scatra",solverparams),
   dt_(0.0),
   maxtime_(0.0),
   stepmax_(0),
   itmax_(0),
   ittol_(1.0),
   upres_(-1),
   write_center_of_mass_(false),
   surftensapprox_(DRT::INPUT::IntegralValue<INPAR::TWOPHASE::SurfaceTensionApprox>(prbdyn.sublist("SURFACE TENSION"),"SURFTENSAPPROX")),
   laplacebeltrami_(DRT::INPUT::IntegralValue<INPAR::TWOPHASE::LaplaceBeltramiCalc>(prbdyn.sublist("SURFACE TENSION"),"LAPLACE_BELTRAMI")),
   velnpi_(Teuchos::null),
   phinpi_(Teuchos::null),
   prbdyn_(prbdyn)
{
  // Needs to stay emtpy
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
XFLUIDLEVELSET::Algorithm::~Algorithm()
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::Init(
    const Teuchos::ParameterList&   prbdyn,         ///< parameter list for global problem
    const Teuchos::ParameterList&   scatradyn,      ///< parameter list for scalar transport subproblem
    const Teuchos::ParameterList&   solverparams,   ///< parameter list for scalar transport solver
    const std::string&              disname,        ///< name of scalar transport discretization
    const bool                      isale           ///< ALE flag
)
{
  // call Setup() in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Init(
      prbdyn,
      scatradyn,
      solverparams,
      disname,
      isale
      );


  //TODO: Combine TWOPHASE and XFLUIDLEVELSET. Create a Parent class, TWOFLUIDCOUPLING or use the existing ScaTraFluidCouplingAlgorithm.
  //        Derived classes will be FLUIDLEVELSET and XFLUIDLEVELSET.
  //        This could also incorporate ELCH and other FLUID-ScaTra coupling schemes.

  // ScaTra Field is not given an initial velocity field. Thus it has to be instantiated before the ScaTra field
  // is being solved for.

  // time-step length, maximum time and maximum number of steps
  dt_      = prbdyn_.get<double>("TIMESTEP");
  maxtime_ = prbdyn_.get<double>("MAXTIME");
  stepmax_ = prbdyn_.get<int>("NUMSTEP");

  //Output specific criterions
  write_center_of_mass_ = DRT::INPUT::IntegralValue<bool>(prbdyn_,"WRITE_CENTER_OF_MASS");

  // (preliminary) maximum number of iterations and tolerance for outer iteration
  ittol_ = prbdyn_.get<double>("CONVTOL");
  itmax_ = prbdyn_.get<int>("ITEMAX");

  upres_ = prbdyn_.get<int>("RESULTSEVRY");

  //Instantiate vectors contatining outer loop increment data
  fsvelincnorm_.reserve(itmax_);
  fspressincnorm_.reserve(itmax_);
  fsphiincnorm_.reserve(itmax_);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::Setup()
{
  // call Setup() in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Setup();

  // Fluid-Scatra Iteration vectors are initialized
  velnpi_ = Teuchos::rcp(new Epetra_Vector(FluidField()->StdVelnp()->Map()),true);//*fluiddis->DofRowMap()),true);
  velnpi_->Update(1.0,*FluidField()->StdVelnp(),0.0);
  phinpi_ = Teuchos::rcp(new Epetra_Vector(ScaTraField()->Phinp()->Map()),true);
  phinpi_->Update(1.0,*ScaTraField()->Phinp(),0.0);

  //Values of velocity field are transferred to ScaTra field. This function overwrites the previous initialization in the constructor
  //of ScaTraFluidCouplingAlgorithm(). This is necessary to correctly initialize a particle algorithm.
  ScaTraField()->Discretization()->ReplaceDofSet(1,Teuchos::rcp(new DRT::DofSet(Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(FluidField()->Discretization())->InitialDofSet())),false);
  SetFluidValuesInScaTra(true);

  // this cannot be done in Init() for some reason
  SetProblemSpecificParameters(prbdyn_);

  return;
}


/*---------------------------------------------------------------------------------------*
| public: algorithm for a instationary XTPF problem                         winter 10/14 |
*----------------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::TimeLoop()
{

  OutputInitialField();

  // time loop
  while (NotFinished())
  {
    IncrementTimeAndStep();

    // prepare time step
    PrepareTimeStep();

    // do outer iteration loop for particular type of algorithm
    OuterLoop();

    // update for next time step
    TimeUpdate();

    // write output to files
    Output();

  } // time loop

return;
}

/*------------------------------------------------------------------------------------------------*
 | public: algorithm for a stationary XTPF problem                                   winter 10/14 |
 *------------------------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::SolveStationaryProblem()
{
  if (Comm().MyPID()==0)
  {
    printf("------Stationary-Xfluid-LevelSet-XFEM------  time step ----------------------------------------\n");
  }

  // check time integration schemes of single fields
  // remark: this was already done in ScaTraFluidCouplingAlgorithm() before
  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary)
    dserror("Scatra time integration scheme is not stationary");

  // write Scatra output (fluid output has been already called in FluidField()->Integrate();
  ScaTraField()->Output();

  //Give Scatra Values to fluid.
  //Needed for curvature etc..
  SetScaTraValuesInFluid();

  // run the simulation, calls the xfluid-"integrate()" routine
  FluidField()->Integrate();

//  // solve level set equation
  if (Comm().MyPID()==0)
    std::cout << "/!\\ warning === Level-set field not solved for Fluid_XFEM_LevelSet problems" << std::endl;


  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::OuterLoop()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";

    printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d\n",
           Time(),maxtime_,dt_,ScaTraField()->MethodTitle().c_str(),Step(),stepmax_);
  }

  // initially solve scalar transport equation
  // (values for intermediate time steps were calculated at the end of PerpareTimeStep)
  DoScaTraField();

  //Prepare variables for convergence check.
  PrepareOuterIteration();

  while (stopnonliniter==false)
  {
    itnum++;

    // solve fluid flow equations
    DoFluidField();

    // solve scalar transport equation
    DoScaTraField();

    // check convergence and stop iteration loop if convergence is achieved
    stopnonliniter = ConvergenceCheck(itnum);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::DoFluidField()
{

  if (Comm().MyPID()==0)
    std::cout<<"\n****************************************\n              FLUID SOLVER\n****************************************\n";


  //Set relevant ScaTra values in Fluid field.
  SetScaTraValuesInFluid();

  //Solve the Fluid field.
  FluidField()->Solve();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::DoScaTraField()
{

  if (Comm().MyPID()==0)
    std::cout<<"\n****************************************\n        SCALAR TRANSPORT SOLVER\n****************************************\n";


  //Set relevant Fluid values in ScaTra field.
  SetFluidValuesInScaTra(false);

  //Solve the ScaTra field.
  ScaTraField()->Solve();

  return;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::TimeUpdate()
{
  // update scalar
  ScaTraField()->Update();

  // update fluid
  FluidField()->Update();

  return;
}

/*---------------------------------------------------------------------------------------*
| Prepares values and variables needed in the outer iteration                            |
*----------------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::PrepareOuterIteration()
{
//  //Update phi for outer loop convergence check
  phinpi_->Update(1.0,*ScaTraField()->Phinp(),0.0);
  velnpi_->Update(1.0,*FluidField()->StdVelnp(),0.0);

  //Clear the vectors containing the data for the partitioned increments
  fsvelincnorm_.clear();
  fspressincnorm_.clear();
  fsphiincnorm_.clear();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::PrepareTimeStep()
{
  // prepare scalar transport time step
  // (+ computation of initial scalar time derivative in first time step)
  ScaTraField()->PrepareTimeStep();

  // prepare fluid time step, among other things, predict velocity field
  FluidField()->PrepareTimeStep();

  // synchronicity check between fluid algorithm and level set algorithms
  if (FluidField()->Time() != Time())
    dserror("Time in Fluid time integration differs from time in two phase flow algorithm");
  if (ScaTraField()->Time() != Time())
    dserror("Time in ScaTra time integration differs from time in two phase flow algorithm");

  return;
}

/*----------------------------------------------------------------------*
 | Set relevant values from ScaTra field in the Fluid field.            |
 *----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::SetScaTraValuesInFluid()
{
  // set level set in fluid field

  switch(FluidField()->TimIntScheme())
  {
  case INPAR::FLUID::timeint_stationary:
  case INPAR::FLUID::timeint_one_step_theta:
  {
    // set the new level set value to the fluid
    Teuchos::RCP<FLD::XFluid> xfluid = Teuchos::rcp_dynamic_cast<FLD::XFluid>(FluidField(), true);

    Teuchos::RCP<Epetra_MultiVector>  smoothedgradphi;
    Teuchos::RCP<Epetra_Vector> nodalcurvature;

    if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_divgrad_normal)
    {
      smoothedgradphi = ScaTraField()->ReconstructGradientAtNodes(ScaTraField()->Phinp());
    }
    else if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_nodal_curvature)
    {
      nodalcurvature = Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->GetNodalCurvature(ScaTraField()->Phinp());
    }
    else if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_laplacebeltrami)
    {
      if(not (laplacebeltrami_==INPAR::TWOPHASE::matrix_non_smoothed))
        smoothedgradphi = ScaTraField()->ReconstructGradientAtNodes(ScaTraField()->Phinp());
    }

    xfluid->SetLevelSetField(ScaTraField()->Phinp(),
                             nodalcurvature,
                             smoothedgradphi,
                             ScaTraField()->Discretization());

  }
  break;
  default:
    dserror("Time integration scheme not supported");
    break;
  }

  return;
}
/*----------------------------------------------------------------------*
 | Set relevant values from Fluid field in the ScaTra field.            |
 *----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::SetFluidValuesInScaTra(bool init)
{

  switch(FluidField()->TimIntScheme())
  {
  case INPAR::FLUID::timeint_stationary:
  case INPAR::FLUID::timeint_one_step_theta:
  {

    const Teuchos::RCP<const Epetra_Vector> convel = FluidField()->StdVelnp();

    Teuchos::RCP<SCATRA::LevelSetAlgorithm> levelsetalgo = Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField());

    if(levelsetalgo==Teuchos::null)
      dserror("Casting ScaTraTimInt onto LevelSetAlgorithm failed. (Null pointer returned)");

    //Transform output:
    DRT::DofSet stddofset_noptr = Teuchos::rcp_dynamic_cast<FLD::XFluid>(FluidField())->DiscretisationXFEM()->InitialDofSet();
    Teuchos::RCP<const DRT::DofSet> stddofs = Teuchos::rcp(new DRT::DofSet(stddofset_noptr));

    levelsetalgo->SetVelocityField(convel,
        Teuchos::null,
        Teuchos::null,
        Teuchos::null, //FluidField()->FsVel(),
        false,
        init);

  }
  break;
  default:
    dserror("Time integration scheme not supported");
    break;
  }


  return;
}

/*------------------------------------------------------------------------------------------------*
 | Fluid - ScaTra partitioned convergence check                                                   |
 |                                                                                                |
 |    The increment between outer iterations is checked, and if lower than given tolerance        |
 |    the loop is exited. However, at least 2 solutions of each field is required to perform a    |
 |    convergence check.                                                                          |
 |                                                                                   winter 02/15 |
 *------------------------------------------------------------------------------------------------*/
bool XFLUIDLEVELSET::Algorithm::ConvergenceCheck(int itnum)
{
  // define flags for Fluid and ScaTra convergence check
   bool fluidstopnonliniter  = false;
   bool scatrastopnonliniter = false;

   bool notconverged = false;

  if (itmax_ <= 0)
    dserror("Set iterations to something reasonable!!!");

  double fsvelincnorm   = 1.0;
  double fspressincnorm = 1.0;
  double fsphiincnorm   = 1.0;
  //Get increment for outer loop of Fluid and ScaTra
  GetOuterLoopIncFluid(fsvelincnorm,fspressincnorm,itnum);
  GetOuterLoopIncScaTra(fsphiincnorm,itnum);

  fsvelincnorm_[itnum-1]=fsvelincnorm;
  fspressincnorm_[itnum-1]=fspressincnorm;
  fsphiincnorm_[itnum-1]=fsphiincnorm;

  if (Comm().MyPID()==0)
  {
    printf("\n|+------ TWO PHASE FLOW CONVERGENCE CHECK:  time step %2d, outer iteration %2d ------+|", Step(), itnum);
    printf("\n|- iter/itermax -|----tol-[Norm]---|-- fluid-inc --|-- press inc --|-- levset inc --|");
  }

  for(int k_itnum=0;k_itnum < itnum; k_itnum++)
  {
    if(k_itnum==0)
    {
      if (Comm().MyPID()==0)
      {
        printf("\n|     %2d/%2d      | %10.3E [L2] |       -       |       -       |   %10.3E   |",
            (k_itnum+1),itmax_,ittol_,fsphiincnorm_[k_itnum]);
      } // end if processor 0 for output
    }
    else
    {
      if (Comm().MyPID()==0)
      {
        printf("\n|     %2d/%2d      | %10.3E [L2] |  %10.3E   |  %10.3E   |   %10.3E   |",
            (k_itnum+1),itmax_,ittol_,fsvelincnorm_[k_itnum],fspressincnorm_[k_itnum],fsphiincnorm_[k_itnum]);
      } // end if processor 0 for output
    }
  }
  if (Comm().MyPID()==0)
    printf("\n|+---------------------------------------------------------------------------------+|\n");


//  // In combust the velocity and pressure component are not separated. To compare values the following output is made.
//  bool compare_with_combust = true;
//  if(compare_with_combust)
//  {
//    double velnormL2 = 1.0;
//    Teuchos::RCP<Epetra_Vector> velnpip = FluidField()->StdVelnp(); //Contains Fluid and Pressure
//    velnpip->Norm2(&velnormL2);
//    if (velnormL2 < 1e-5) velnormL2 = 1.0;
//
//    double fgvelnormL2 = 1.0;
//
//    // compute increment and L2-norm of increment
//    Teuchos::RCP<Epetra_Vector> incvel = Teuchos::rcp(new Epetra_Vector(velnpip->Map()),true);
//    incvel->Update(1.0,*velnpip,-1.0,*velnpi_,0.0);
//    incvel->Norm2(&fgvelnormL2);
//
//    if (Comm().MyPID()==0)
//    {
//      printf("\n|+------------------------ FGI ------------------------+|");
//      printf("\n|iter/itermax|----tol-[Norm]--|-fluid inc--|-g-func inc (i+1)-|");
//      printf("\n|   %2d/%2d    | %10.3E[L2] | %10.3E | %10.3E |",
//          itnum,itmax_,ittol_,fgvelnormL2/velnormL2,fsphinormL2/phinormL2);
//      printf("\n|+-----------------------------------------------------+|\n");
//    } // end if processor 0 for output
//  }

  if ((fsvelincnorm <= ittol_) and (fspressincnorm <= ittol_) and itnum > 1)
    fluidstopnonliniter = true;

  if ((fsphiincnorm <= ittol_))
    scatrastopnonliniter = true;


  //If tolerance or number of maximum iterations are reached
  if((fluidstopnonliniter and scatrastopnonliniter) or (itnum>=itmax_))
  {
    notconverged=true;
  }

  if (Comm().MyPID()==0)
  {
    if ((itnum == stepmax_) and (notconverged == true))
    {
      printf("|+---------------- not converged ----------------------+|");
      printf("\n|+-----------------------------------------------------+|\n");
    }
  }

  return notconverged;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::Output()
{

  if(Step()%upres_ == 0) //Only perform output for given RESULTSEVRY in Control Algo section of input.
  {
    FluidField()->Output();
    ScaTraField()->Output();
  }

  if(write_center_of_mass_)
  {
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->MassCenterUsingSmoothing();
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: output of initial field                                                winter 07/14 |
 *------------------------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::OutputInitialField()
{
  if (Step() == 0)
  {
    // output fluid initial state
    if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_stationary)
      FluidField()->Output();

    // output Levelset function initial state
    if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary)
      ScaTraField()->Output();
  }

  if(write_center_of_mass_)
  {
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->MassCenterUsingSmoothing();
  }

  return;
}

/*----------------------------------------------------------------------*
 | perform result test                                     winter 06/14 |
 *----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::TestResults()
{

  // perform result tests if required
  DRT::Problem::Instance()->AddFieldTest(FluidField()->CreateFieldTest());
  //DRT::Problem::Instance()->TestAll(Comm());

  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_gen_alpha)
  {
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->TestResults();
  }
  else
  {
    dserror("Unknown time integration for Level Set field in Two Phase Flow problems.");
  }

  return;
}

/* -------------------------------------------------------------------------------*
 | Restart a X-two phase problem                                              winter|
 * -------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::Restart(int step)
{
  if (Comm().MyPID()==0)
  {
    std::cout << "##################################################" << std::endl;
    std::cout << "#                                                #" << std::endl;
    std::cout << "#     Restart of T(wo) P(hase) F(low) problem    #" << std::endl;
//    if (restartscatrainput)
//    {
//      std::cout << "#                                                #" << std::endl;
//      std::cout << "#   -Restart with scalar field from input file   #" << std::endl;
//    }
    std::cout << "#                                                #" << std::endl;
    std::cout << "##################################################" << std::endl;

    std::cout << "##########################################################################" << std::endl;
    std::cout << "#                                                                        #" << std::endl;
    std::cout << "#     WARNING: ONLY RESTART FROM XFLUID AND SCATRA ALLOWED FOR NOW!!!    #" << std::endl;
    std::cout << "#                                                                        #" << std::endl;
    std::cout << "##########################################################################" << std::endl;
  }



  FluidField()->ReadRestart(step);
  ScaTraField()->ReadRestart(step);
  SetTimeStep(FluidField()->Time(),step);

  //Needed for particle restart
  SetFluidValuesInScaTra(true);

  return;
}

/* -------------------------------------------------------------------------------*
 | Set Problem Specific Parameters                                          winter|
 * -------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::SetProblemSpecificParameters(const Teuchos::ParameterList& prbdyn)
{
    surftensapprox_ = DRT::INPUT::IntegralValue<INPAR::TWOPHASE::SurfaceTensionApprox>(prbdyn.sublist("SURFACE TENSION"),"SURFTENSAPPROX");

    if(surftensapprox_==INPAR::TWOPHASE::surface_tension_approx_laplacebeltrami)
      laplacebeltrami_ = DRT::INPUT::IntegralValue<INPAR::TWOPHASE::LaplaceBeltramiCalc>(prbdyn.sublist("SURFACE TENSION"),"LAPLACE_BELTRAMI");

    //SAFETY-CHECKS
    if(DRT::INPUT::IntegralValue<bool>(prbdyn.sublist("SURFACE TENSION"),"L2_PROJECTION_SECOND_DERIVATIVES"))
      dserror("Second L2-projected derivatives can not be calculated as of now for the Level Set.");

    if(DRT::INPUT::IntegralValue<INPAR::TWOPHASE::SmoothGradPhi>(prbdyn.sublist("SURFACE TENSION"),"SMOOTHGRADPHI")!=INPAR::TWOPHASE::smooth_grad_phi_l2_projection)
      dserror("No other smoothing for the gradient of the level set other than L2 is allowed for now.");

    if(DRT::INPUT::IntegralValue<INPAR::TWOPHASE::NodalCurvatureCalc>(prbdyn.sublist("SURFACE TENSION"),"NODAL_CURVATURE")!=INPAR::TWOPHASE::l2_projected)
      dserror("No other way to calculate the nodal curvature than L2.");

    if(prbdyn.sublist("SURFACE TENSION").get<double>("SMOOTHING_PARAMETER")!=0.0)
     dserror("No smoothing available for now.");

    Teuchos::RCP<FLD::XFluid> xfluid = Teuchos::rcp_dynamic_cast<FLD::XFluid>(FluidField(), true);
    xfluid->InitTwoPhaseSurftensParameters(surftensapprox_,laplacebeltrami_);

  return;
}

void XFLUIDLEVELSET::Algorithm::GetOuterLoopIncFluid(double& fsvelincnorm, double& fspressincnorm, int itnum)
{
  Teuchos::RCP<const Epetra_Vector> velnpip = FluidField()->StdVelnp(); //Contains Fluid and Pressure

  //Extract velocity and pressure components.
  Teuchos::RCP<const LINALG::MapExtractor> velpresspliter = Teuchos::rcp_dynamic_cast<FLD::XFluid>(FluidField())->VelPresSplitterStd();

  Teuchos::RCP<const Epetra_Vector> onlyvel   = velpresspliter->ExtractOtherVector(velnpip);
  Teuchos::RCP<const Epetra_Vector> onlypress = velpresspliter->ExtractCondVector(velnpip);

  Teuchos::RCP<Epetra_Vector> onlyveli = Teuchos::rcp(new Epetra_Vector(onlyvel->Map()),true);
  Teuchos::RCP<Epetra_Vector> onlypressi = Teuchos::rcp(new Epetra_Vector(onlypress->Map()),true);

  if(itnum>1)
  {
    onlyveli   = velpresspliter->ExtractOtherVector(velnpi_);
    onlypressi = velpresspliter->ExtractCondVector(velnpi_);
  }

  double velnormL2   = 1.0;
  double pressnormL2 = 1.0;

  onlyvel->Norm2(&velnormL2);
  onlypress->Norm2(&pressnormL2);

  if (velnormL2 < 1e-5) velnormL2 = 1.0;
  if (pressnormL2 < 1e-5) pressnormL2 = 1.0;

  double fsvelnormL2   = 1.0;
  double fspressnormL2 = 1.0;

  // compute increment and L2-norm of increment
  //-----------------------------------------------------
  Teuchos::RCP<Epetra_Vector> incvel = Teuchos::rcp(new Epetra_Vector(onlyvel->Map()),true);
  incvel->Update(1.0,*onlyvel,-1.0,*onlyveli,0.0);
  incvel->Norm2(&fsvelnormL2);

  Teuchos::RCP<Epetra_Vector> incpress = Teuchos::rcp(new Epetra_Vector(onlypress->Map()),true);
  incpress->Update(1.0,*onlypress,-1.0,*onlypressi,0.0);
  incpress->Norm2(&fspressnormL2);
  //-----------------------------------------------------

  fsvelincnorm  =fsvelnormL2/velnormL2;
  fspressincnorm=fspressnormL2/pressnormL2;

#if DEBUG
//-------------------------
  std::cout << "fsvelnormL2: " << fsvelnormL2 << std::endl;
  std::cout << "velnormL2: " << velnormL2 << std::endl << std::endl;

  std::cout << "fspressnormL2: " << fspressnormL2 << std::endl;
  std::cout << "pressnormL2: " << pressnormL2 << std::endl<< std::endl;
//-------------------------
#endif

  velnpi_->Update(1.0,*velnpip,0.0);

}

void XFLUIDLEVELSET::Algorithm::GetOuterLoopIncScaTra(double& fsphiincnorm, int itnum)
{
  Teuchos::RCP<const Epetra_Vector> phinpip = ScaTraField()->Phinp();

  double phinormL2   = 1.0;

  phinpip->Norm2(&phinormL2);
  if (phinormL2 < 1e-5) phinormL2 = 1.0;
  double fsphinormL2   = 1.0;

  // compute increment and L2-norm of increment
  //-----------------------------------------------------
  Teuchos::RCP<Epetra_Vector> incphi = Teuchos::rcp(new Epetra_Vector(phinpip->Map(),true));
  incphi->Update(1.0,*phinpip,-1.0,*phinpi_,0.0);
  incphi->Norm2(&fsphinormL2);
  //-----------------------------------------------------

  fsphiincnorm = fsphinormL2/phinormL2;

#if DEBUG
//-------------------------
  std::cout << "fsphinormL2: " << fsphinormL2 << std::endl;
  std::cout << "phinormL2: " << phinormL2 << std::endl<< std::endl;
//-------------------------
#endif


  phinpi_->Update(1.0,*phinpip,0.0);

}
