/*!----------------------------------------------------------------------
\file fs3i_partitioned_2wc.cpp
\brief Algorithmic routines for partitioned solution approaches to
       fluid-structure-scalar-scalar interaction (FS3I) specifically
       related to two-way-coupled problem configurations

\level 2

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10364


*----------------------------------------------------------------------*/


#include "fs3i_partitioned_2wc.H"

#include "../drt_fsi/fsi_monolithic.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"
#include "../drt_scatra/scatra_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_fluid/fluid_timint_loma.H"

#include "../drt_scatra/scatra_timint_loma.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::PartFS3I_2WC::PartFS3I_2WC(const Epetra_Comm& comm)
  : PartFS3I(comm),
    itmax_(-1),
    ittol_(-1.0)
{
  // constructor is supposed to stay empty
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::Init()
{
  // call Init() in base class
  FS3I::PartFS3I::Init();

  //---------------------------------------------------------------------
  // get input parameters for two-way-coupled problems, which is
  // thermo-fluid-structure interaction, for the time being
  //---------------------------------------------------------------------
  const Teuchos::ParameterList& fs3idyn = DRT::Problem::Instance()->FS3IDynamicParams();
  ittol_ = fs3idyn.get<double>("CONVTOL");
  itmax_ = fs3idyn.get<int>("ITEMAX");

  // flag for constant thermodynamic pressure
  consthermpress_ = fs3idyn.get<std::string>("CONSTHERMPRESS");

  // define fluid- and structure-based scalar transport problem
  fluidscatra_     = scatravec_[0];
  structurescatra_ = scatravec_[1];

  // add proxy of fluid degrees of freedom to scatra discretization
  if(fluidscatra_->ScaTraField()->Discretization()->AddDofSet(fsi_->FluidField()->Discretization()->GetDofSetProxy()) != 1)
    dserror("Scatra discretization has illegal number of dofsets!");

  // add proxy of structure degrees of freedom to scatra discretization
  if(structurescatra_->ScaTraField()->Discretization()->AddDofSet(fsi_->StructureField()->Discretization()->GetDofSetProxy()) != 1)
    dserror("Scatra discretization has illegal number of dofsets!");

  // generate proxy of dof set for structure-based scalar transport
  // problem to be used by structure field
  Teuchos::RCP<DRT::DofSet> structurescatradofset = structurescatra_->ScaTraField()->Discretization()->GetDofSetProxy();

  // check number of dof sets in structure field
  if (fsi_->StructureField()->Discretization()->AddDofSet(structurescatradofset)!=1)
    dserror("Incorrect number of dof sets in structure field!");

  if (volume_fieldcouplings_[0]==INPAR::FS3I::coupling_nonmatch or volume_fieldcouplings_[1]==INPAR::FS3I::coupling_nonmatch )
    dserror("Mortar volume coupling is not tested for thermo-fs3i.");

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::Setup()
{
  // call Setup() in base class
  FS3I::PartFS3I::Setup();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::Timeloop()
{
  CheckIsInit();
  CheckIsSetup();

  InitialCalculations();

  while (NotFinished())
  {
    IncrementTimeAndStep();

    PrepareTimeStep();

    OuterLoop();

    TimeUpdateAndOutput();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::InitialCalculations()
{
  // set initial values for mesh displacement field
  SetMeshDisp();

  // set initial fluid velocity field for evaluation of initial scalar
  // time derivative in fluid-based scalar transport
  fluidscatra_->ScaTraField()->SetVelocityField(fsi_->FluidField()->Velnp(),
                                               Teuchos::null,
                                               Teuchos::null,
                                               Teuchos::null,
                                               1);

  // set initial value of thermodynamic pressure in fluid-based scalar
  // transport
  Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->SetInitialThermPressure();

  // energy conservation: compute initial time derivative of therm. pressure
  // mass conservation: compute initial mass (initial time deriv. assumed zero)
  if (consthermpress_=="No_energy")
    Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->ComputeInitialThermPressureDeriv();
  else if (consthermpress_=="No_mass")
    Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->ComputeInitialMass();

  // set initial scalar field and thermodynamic pressure for evaluation of
  // Neumann boundary conditions in fluid at beginning of first time step
  fsi_->FluidField()->SetScalarFields(fluidscatra_->ScaTraField()->Phinp(),
                                       Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->ThermPressNp(),
                                       Teuchos::null,
                                       fluidscatra_->ScaTraField()->Discretization());
  // prepare time loop for FSI
  fsi_->PrepareTimeloop();

  // output of initial state for scalar values
  //ScatraOutput();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::PrepareTimeStep()
{
  CheckIsInit();
  CheckIsSetup();

  // set mesh displacement and velocity fields for present time step
  SetFSISolution();

  // prepare time step for both fluid- and structure-based scatra field
  // (+ computation of initial scalar time derivative in first time step)
  fluidscatra_->ScaTraField()->PrepareTimeStep();
  structurescatra_->ScaTraField()->PrepareTimeStep();

  // predict thermodynamic pressure and time derivative
  // (if not constant or based on mass conservation)
  if (consthermpress_=="No_energy")
    Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->PredictThermPressure();

  // prepare time step for fluid, structure and ALE fields
  fsi_->PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::OuterLoop()
{
  CheckIsInit();
  CheckIsSetup();

#ifdef PARALLEL
  const Epetra_Comm& comm = scatravec_[0]->ScaTraField()->Discretization()->Comm();
#else
  Epetra_SerialComm comm;
#endif

  int  itnum = 0;
  bool stopnonliniter = false;

  if (comm.MyPID()==0)
  {
    std::cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";

    printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d\n",
           fluidscatra_->ScaTraField()->Time(),timemax_,dt_,fluidscatra_->ScaTraField()->MethodTitle().c_str(),fluidscatra_->ScaTraField()->Step(),numstep_);
  }

  // the following already done in PrepareTimeStep:
  // set FSI values required in scatra
  //SetFSIValuesInScaTra();

  // initially solve coupled scalar transport equation system
  // (values for intermediate time steps were calculated at the end of PrepareTimeStep)
  if (comm.MyPID()==0) std::cout<<"\n****************************************\n        SCALAR TRANSPORT SOLVER\n****************************************\n";
  ScatraEvaluateSolveIterUpdate();

  while (stopnonliniter==false)
  {
    itnum++;

    // in case of non-constant thermodynamic pressure: compute
    // (either based on energy conservation or based on mass conservation)
    if (consthermpress_=="No_energy")
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->ComputeThermPressure();
    else if (consthermpress_=="No_mass")
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->ComputeThermPressureFromMassCons();

    // set fluid- and structure-based scalar transport values required in FSI
    SetScaTraValuesInFSI();

    // solve FSI system
    if (comm.MyPID()==0) std::cout<<"\n****************************************\n               FSI SOLVER\n****************************************\n";
    fsi_->TimeStep(fsi_);

    // set FSI values required in scatra (will be done in the following
    // routine, for the time being)
    //SetFSIValuesInScaTra();
    // set mesh displacement and velocity fields
    SetFSISolution();

    // solve scalar transport equation
    if (comm.MyPID()==0) std::cout<<"\n****************************************\n        SCALAR TRANSPORT SOLVER\n****************************************\n";
    ScatraEvaluateSolveIterUpdate();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(itnum);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*void FS3I::PartFS3I_2WC::SetFSIValuesInScaTra()
{
  // set respective field vectors for velocity/pressure, acceleration
  // and discretization based on time-integration scheme
  if (fsi_->FluidField()->TimIntScheme() == INPAR::FLUID::timeint_afgenalpha)
    fluidscatra_->ScaTraField()->SetVelocityField(fsi_->FluidField()->Velaf(),
                                                fsi_->FluidField()->Accam(),
                                                Teuchos::null,
                                                fsi_->FluidField()->FsVel(),
                                                Teuchos::null,
                                                fsi_->FluidField()->Discretization());
  else
    fluidscatra_->ScaTraField()->SetVelocityField(fsi_->FluidField()->Velnp(),
                                                fsi_->FluidField()->Hist(),
                                                Teuchos::null,
                                                fsi_->FluidField()->FsVel(),
                                                Teuchos::null,
                                                fsi_->FluidField()->Discretization());
}*/


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::SetScaTraValuesInFSI()
{
    // set scalar and thermodynamic pressure values as well as time
    // derivatives and discretization based on time-integration scheme
    /*if (fsi_->FluidField()->TimIntScheme() == INPAR::FLUID::timeint_afgenalpha)
    {
      dynamic_cast<FLD::TimIntLoma&>(fsi_->FluidField()).SetIterScalarFields(fluidscatra_->ScaTraField()->Phiaf(),
                                           fluidscatra_->ScaTraField()->Phiam(),
                                           fluidscatra_->ScaTraField()->Phidtam(),
                                           Teuchos::null,
                                           fluidscatra_->ScaTraField()->ThermPressAf(),
                                           fluidscatra_->ScaTraField()->ThermPressAm(),
                                           fluidscatra_->ScaTraField()->ThermPressDtAf(),
                                           fluidscatra_->ScaTraField()->ThermPressDtAm(),
                                           fluidscatra_->ScaTraField()->Discretization());

      fsi_->StructureField()->ApplyCouplingState(structurescatra_->ScaTraField()->Phiaf(),"temperature");
    }
    else
    {*/
     Teuchos::rcp_dynamic_cast<FLD::TimIntLoma>(fsi_->FluidField())->SetIterScalarFields(
                                           FluidScalarToFluid(fluidscatra_->ScaTraField()->Phinp()),
                                           FluidScalarToFluid(fluidscatra_->ScaTraField()->Phin()),
                                           FluidScalarToFluid(fluidscatra_->ScaTraField()->Phidtnp()),
                                           Teuchos::null,
                                           Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->ThermPressNp(),
                                           Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->ThermPressN(),
                                           Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->ThermPressDtNp(),
                                           Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->ThermPressDtNp(),
                                           Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->Discretization());
      fsi_->StructureField()->Discretization()->SetState(1,"temperature", StructureToStructureScalar(structurescatra_->ScaTraField()->Phinp()) );
    //}
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FS3I::PartFS3I_2WC::ConvergenceCheck(int itnum)
{
#ifdef PARALLEL
  const Epetra_Comm& comm = scatravec_[0]->ScaTraField()->Discretization()->Comm();
#else
  Epetra_SerialComm comm;
#endif

  // define flags for fluid and scatra convergence check
  bool fluidstopnonliniter  = false;
  bool scatrastopnonliniter = false;

  // dump on screen
  if (comm.MyPID()==0) std::cout<<"\n****************************************\n  CONVERGENCE CHECK FOR ITERATION STEP\n****************************************\n";

  // fsi convergence check
  if (fsi_->NoxStatus() == NOX::StatusTest::Converged) fluidstopnonliniter = true;

  // scatra convergence check
  scatrastopnonliniter = ScatraConvergenceCheck(itnum);

  // warn if itemax is reached without convergence of FSI solver,
  // but proceed to next timestep
  if ((itnum == itmax_) and (fluidstopnonliniter == false))
  {
    fluidstopnonliniter=true;
    if (comm.MyPID() == 0)
    {
      printf("\n");
      printf(">>>>>> FSI solver not converged in itemax steps!\n");
      printf("\n");
    }
  }

  if (fluidstopnonliniter == true and scatrastopnonliniter == true) return true;
  else                                                              return false;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FS3I::PartFS3I_2WC::ScatraConvergenceCheck(int itnum)
{
#ifdef PARALLEL
  const Epetra_Comm& comm = scatravec_[0]->ScaTraField()->Discretization()->Comm();
#else
  Epetra_SerialComm comm;
#endif

  // define flags for convergence check for scatra fields
  bool scatra1stopnonliniter = false;
  bool scatra2stopnonliniter = false;

  // convergence check of scatra fields
  if (comm.MyPID() == 0)
  {
    std::cout<<"\n****************************************\n         SCALAR TRANSPORT CHECK\n****************************************\n";
    std::cout<<"\n****************************************\n   FLUID-BASED SCALAR TRANSPORT CHECK\n****************************************\n";
  }
  scatra1stopnonliniter = Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->ConvergenceCheck(itnum,itmax_,ittol_);

  if (comm.MyPID() == 0) std::cout<<"\n****************************************\n STRUCTURE-BASED SCALAR TRANSPORT CHECK\n****************************************\n";
    dserror("ConvergenceCheck in scatra currently only for loma scatra!Fix this!");
    //scatra2stopnonliniter = scatravec_[1]->ScaTraField()->ConvergenceCheck(itnum,itmax_,ittol_);

  if (scatra1stopnonliniter == true and scatra2stopnonliniter == true) return true;
  else                                                                 return false;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::TimeUpdateAndOutput()
{
  // prepare output for FSI
  fsi_->PrepareOutput();

  // update fluid- and structure-based scalar transport
  UpdateScatraFields();

  // in case of non-constant thermodynamic pressure: update
  if (consthermpress_=="No_energy" or consthermpress_=="No_mass")
    Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->UpdateThermPressure();

  // update structure, fluid and ALE
  fsi_->Update();

  // set scalar and thermodynamic pressure at n+1 and SCATRA trueresidual
  // for statistical evaluation and evaluation of Neumann boundary
  // conditions at the beginning of the subsequent time step
  fsi_->FluidField()->SetScalarFields(FluidScalarToFluid(fluidscatra_->ScaTraField()->Phinp()),
                                      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntLoma>(fluidscatra_->ScaTraField())->ThermPressNp(),
                                      FluidScalarToFluid(fluidscatra_->ScaTraField()->TrueResidual()),
                                      fluidscatra_->ScaTraField()->Discretization());

  // Note: The order is important here! Herein, control file entries are
  // written, defining the order in which the filters handle the
  // discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  //fsi_->FluidField()->StatisticsAndOutput();
  fsi_->Output();

  // output of fluid- and structure-based scalar transport
  ScatraOutput();

  return;
}


