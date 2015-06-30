/*!----------------------------------------------------------------------
\file fpsi_partitioned.cpp

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de

</pre>

*----------------------------------------------------------------------*/
#include "../drt_inpar/inpar_fpsi.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_poroelast/poroelast_utils.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_io/io_control.H"

// FPSI includes
#include "fpsi_utils.H"
#include "fpsi_partitioned.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FPSI::Partitioned::Partitioned(const Epetra_Comm& comm,
                               const Teuchos::ParameterList& timeparams)
    : FPSI_Base(comm,timeparams)
{
  //DRT::Problem* problem = DRT::Problem::Instance();

  SetDefaultParameters(timeparams);

  dt_      = timeparams.get<double> ("TIMESTEP");
  numstep_ = timeparams.get<int>    ("NUMSTEP");
  timemax_ = timeparams.get<double> ("MAXTIME");

  step_ = 0;
  time_ = 0.;

   ///////////////////////////////////////////////////////////////////////////////
  // Call the poroelast_subproblem
  poroelast_subproblem_ = POROELAST::UTILS::CreatePoroAlgorithm(timeparams, comm);
  if(poroelast_subproblem_ == Teuchos::null)
  {
    dserror("Red Alert!!! Call of the poroelast subproblem failed! ... Roundhouse-Kick!!! ");
  }
   /////////////////////////////////////////////////////////////////////////////
  // Call the fluid subproblem
  Teuchos::RCP< ::ADAPTER::FluidMovingBoundaryBaseAlgorithm> FluidMBBase =
  Teuchos::rcp(new ADAPTER::FluidMovingBoundaryBaseAlgorithm(timeparams,"FPSICoupling"));
  fluid_subproblem_ = FluidMBBase->MBFluidField();
  if(fluid_subproblem_ == Teuchos::null)
  {
    dserror("Red Alert!!! Call of the fluid subproblem failed! ... Roundhouse-Kick!!! ");
  }

   /////////////////////////////////////////////////////////////////////////////
  // Interface coupling
  coupsf_ = Teuchos::rcp(new ADAPTER::Coupling());
  if (DRT::INPUT::IntegralValue<int>(timeparams,"COUPMETHOD"))
    {
      matchingnodes_ = true;
      const int ndim = DRT::Problem::Instance()->NDim();
      coupsf_->SetupConditionCoupling(*poroelast_subproblem_->StructureField()->Discretization(),            //masterdis (const)
                                       poroelast_subproblem_->StructureField()->Interface()->FPSICondMap(),   //mastercondmap
                                      *fluid_subproblem_->Discretization(),                                  //slavedis (const)
                                       fluid_subproblem_->Interface()->FSICondMap(),                          //slavecondmap
                                      "FPSICoupling",                                                        //condname
                                       ndim);
      //numdof

      if (coupsf_->MasterDofMap()->NumGlobalElements()==0)
        dserror("No nodes in matching FPSI interface. Empty FPSI coupling condition?");
    }
    else
    {
      matchingnodes_ = false;
      dserror("COUPMETHOD not set in Input !!! Set COUPMETHOD in FPSI Dynamic !!!"
              "\n Partitioned standard FSI uses Mortar-Coupling instead.\n"
              "However, that doesn't work here  :-P ... so far ...");
    }

}

/*----------------------------------------------------------------------*/
/*                            TIMELOOP                                  */
/*----------------------------------------------------------------------*/
void FPSI::Partitioned::Timeloop()
{
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
/*                      IMPLEMENT TIMELOOP                              */
/*----------------------------------------------------------------------*/
void FPSI::Partitioned::InitialCalculations()
{
  Teuchos::RCP<std::ofstream> log;
  if (Comm().MyPID()==0)
  {
    std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
    s.append(".iteration");
    log = Teuchos::rcp(new std::ofstream(s.c_str()));
    (*log) << "# num procs      = " << Comm().NumProc() << "\n"
           << "#\n"
           << "# step | time | time/step | #nliter  |R|  #liter  User\n"
      ;
  }

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Partitioned::IncrementTimeAndStep()
{
  step_ += 1;
  time_ += dt_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Partitioned::PrepareTimeStep()
{
  PrintHeader();

  poroelast_subproblem_->PrepareTimeStep();
  fluid_subproblem_    ->PrepareTimeStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Partitioned::OuterLoop()
{
#ifdef PARALLEL
  const Epetra_Comm& comm = fluid_subproblem_->Discretization()->Comm();
#else
  Epetra_SerialComm comm;
#endif

  int  itnum = 0;
  bool stopnonliniter = false;

  if (comm.MyPID()==0)
  {
    std::cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";

    printf("TIME: %11.4E/%11.4E  DT = %11.4E  STEP = %4d/%4d\n",
           poroelast_subproblem_->FluidField()->Time(),timemax_,dt_,poroelast_subproblem_->FluidField()->Step(),numstep_);
  }


  while (stopnonliniter==false)
  {
    itnum++;

    // solve FSI system
    if (comm.MyPID()==0) std::cout<<"\n****************************************\n     FLUID SOLVER     \n****************************************\n";
    {
      fluid_subproblem_ -> NonlinearSolve();
    }

    if (comm.MyPID()==0) std::cout<<"\n****************************************\n   POROELAST SOLVER   \n****************************************\n";
    {
      poroelast_subproblem_ -> Solve();
    }

    stopnonliniter = ConvergenceCheck(itnum);

  } // nonliniter
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Partitioned::TimeUpdateAndOutput()
{
  poroelast_subproblem_   ->PrepareOutput();
  poroelast_subproblem_   ->Update();
  fluid_subproblem_       ->Update();

  idispn_ = poroelast_subproblem_->StructureField()->ExtractInterfaceDispn(true);
  iveln_  = FluidToStruct(fluid_subproblem_->ExtractInterfaceVeln());

  poroelast_subproblem_   ->Output();
  fluid_subproblem_       ->Output();
}


/*----------------------------------------------------------------------*/
/*                     REQUIRED FPSI_DYN METHODS                        */
/*----------------------------------------------------------------------*/
void FPSI::Partitioned::ReadRestart(int restartstep)
{
  poroelast_subproblem_     -> ReadRestart(restartstep);

  time_ = fluid_subproblem_->  ReadRestart(restartstep);
  SetTimeStep(time_,restartstep);
}

/*----------------------------------------------------------------------*
 | redistribute the FPSI interface                           thon 11/14 |
 *----------------------------------------------------------------------*/
void FPSI::Partitioned::RedistributeInterface()
{
  DRT::Problem* problem = DRT::Problem::Instance();
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();
  Teuchos::RCP<FPSI::UTILS> FPSI_UTILS = FPSI::UTILS::Instance();

  if(comm.NumProc() > 1) //if we have more than one processor, we need to redistribute at the FPSI interface
  {
    Teuchos::RCP<std::map<int,int> > Fluid_PoroFluid_InterfaceMap = FPSI_UTILS->Get_Fluid_PoroFluid_InterfaceMap();
    Teuchos::RCP<std::map<int,int> > PoroFluid_Fluid_InterfaceMap = FPSI_UTILS->Get_PoroFluid_Fluid_InterfaceMap();

    FPSI_UTILS->RedistributeInterface(problem->GetDis("fluid")    ,problem->GetDis("porofluid"),"FPSICoupling",*PoroFluid_Fluid_InterfaceMap);
    FPSI_UTILS->RedistributeInterface(problem->GetDis("ale")      ,problem->GetDis("porofluid"),"FPSICoupling",*PoroFluid_Fluid_InterfaceMap);
    FPSI_UTILS->RedistributeInterface(problem->GetDis("porofluid"),problem->GetDis("fluid")    ,"FPSICoupling",*Fluid_PoroFluid_InterfaceMap);
    FPSI_UTILS->RedistributeInterface(problem->GetDis("structure"),problem->GetDis("fluid")    ,"FPSICoupling",*Fluid_PoroFluid_InterfaceMap);

    // Material pointers need to be reset after redistribution.
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(problem->GetDis("structure"), problem->GetDis("porofluid"));
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Partitioned::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(poroelast_subproblem_->StructureField()->CreateFieldTest());
  problem->AddFieldTest(poroelast_subproblem_->FluidField()->CreateFieldTest());

  problem->TestAll(comm);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Partitioned::SetupSystem()
{
 poroelast_subproblem_ -> SetupSystem();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FPSI::Partitioned::StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FPSI::Partitioned::FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FPSI::Partitioned::StructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FPSI::Partitioned::FluidToStruct(Teuchos::RCP<const Epetra_Vector> iv) const
{
  if (matchingnodes_)
  {
    return coupsf_->SlaveToMaster(iv);
  }
  else
  {
    dserror("No matching nodes at the interface !! This is not supported in FPSI, yet.\n"
            "Standard FSI supports Mortar-Coupling in this case.");
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Partitioned::SetDefaultParameters(const Teuchos::ParameterList& fpsidynparams)
{
  // Get the top level parameter list
  //Teuchos::ParameterList& nlParams = list;

  // to do

  switch (DRT::INPUT::IntegralValue<int>(fpsidynparams,"COUPALGO"))
    {
    case partitioned:
    {
      // to do
      break;
    }
    }

}

/*----------------------------------------------------------------------*/
/*                      IMPLEMENT OUTER LOOP                            */
/*----------------------------------------------------------------------*/
bool FPSI::Partitioned::ConvergenceCheck(int itnum)
{
  // to do
  return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FPSI::Partitioned::InterfaceDisp()
{
  // extract displacements
  return poroelast_subproblem_->StructureField()->ExtractInterfaceDispnp(true);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FPSI::Partitioned::InterfaceForce()
{
  // extract forces
  return FluidToStruct(fluid_subproblem_->ExtractInterfaceForces());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FPSI::Partitioned::InterfaceVelocity(Teuchos::RCP<const Epetra_Vector> idisp_npone)
{
  // should calculate the tangential velocity component (Beavers-Joseph) and
  // the normal velocity to be applied to the fluid field

  return;
}

