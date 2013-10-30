/*----------------------------------------------------------------------*/
/*!
\file fpsi_monolithic.cpp

\brief General framework for monolithic fpsi solution schemes

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de

</pre>
*/

/*----------------------------------------------------------------------*/
// GENERAL includes
#include <sstream>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

// FPSI includes
#include "fpsi_defines.H"
#include "fpsi_monolithic.H"
#include "fpsi_utils.H"

// POROELAST includes
#include "../drt_poroelast/poroelast_utils.H"
#include "../drt_poroelast/poro_base.H"
#include "../drt_poroelast/poroelast_monolithic.H"

// FSI includes
#include "../drt_fsi/fsi_debugwriter.H"
#include "../drt_fsi/fsi_debugwriter.H"
#include "../drt_fsi/fsi_statustest.H"

// LINALG includes
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

// INPAR includes
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_solver.H"

// drt_lib includes
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_condition_utils.H"

// drt_adapter includes
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"

// STRUCTURE includes
#include "../drt_structure/stru_aux.H"

// FLUID includes
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

// ALE includes
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_ale/ale_lin.H"

// OTHER includes
#include "../drt_constraint/constraint_manager.H"
#include "../drt_io/io_control.H"
#include "../drt_lib/drt_dofset.H"

#include <iostream>
#include <fstream>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

FPSI::MonolithicBase::MonolithicBase(const Epetra_Comm& comm,
                                     const Teuchos::ParameterList& fpsidynparams,
                                     const Teuchos::ParameterList& poroelastdynparams)
  : FPSI_Base(comm,fpsidynparams)
{
  // create instance of poroelast subproblem
  poroelast_subproblem_ = Teuchos::rcp(new POROELAST::Monolithic(comm, poroelastdynparams));
  // ask base algorithm for the fluid time integrator
  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fluiddynparams = problem->FluidDynamicParams();
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fpsidynparams,fluiddynparams,"fluid",true));
  fluid_subproblem_ = fluid->FluidFieldrcp();
  // ask base algorithm for the ale time integrator
  Teuchos::RCP<ALE::AleBaseAlgorithm> ale = Teuchos::rcp(new ALE::AleBaseAlgorithm(fpsidynparams));
  ale_ = ale->AleFieldrcp();

  couple_porofluid_fluid_matching_ = Teuchos::rcp(new ADAPTER::Coupling());
  couple_porofluid_fluid_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupsf_  = Teuchos::rcp(new ADAPTER::Coupling());
  coupsa_  = Teuchos::rcp(new ADAPTER::Coupling());
  coupfa_  = Teuchos::rcp(new ADAPTER::Coupling());
  innercoupsf_ = Teuchos::rcp(new ADAPTER::Coupling());
  smallcoupsf_ = Teuchos::rcp(new ADAPTER::Coupling());

  Teuchos::RCP<FPSI::UTILS> FPSI_UTILS = FPSI::UTILS::Instance();

  Fluid_PoroFluid_InterfaceMap = FPSI_UTILS->Get_Fluid_PoroFluid_InterfaceMap();
  PoroFluid_Fluid_InterfaceMap = FPSI_UTILS->Get_PoroFluid_Fluid_InterfaceMap();

  // build a proxy of the fluid discretization for the structure field
  aledofset        = Teuchos::null;
  aledofset =                      AleField()->Discretization()->GetDofSetProxy();

  if (FluidField()->Discretization()->AddDofSet(aledofset) != 1)
  {
    dserror("Yippie-ei-yeah ... Schweinebacke ... ");
  }


} // MonolithicBase

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FPSI::MonolithicBase::~MonolithicBase()
{
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::ReadRestart(int step)
{
  PoroField()     ->ReadRestart(step);
  FluidField()    ->ReadRestart(step);
  AleField()      ->ReadRestart(step);

  SetTimeStep(FluidField()->Time(), FluidField()->Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::PrepareTimeStep()
{
  IncrementTimeAndStep();
  PrintHeader();

  PoroField()   ->PrepareTimeStep();
  FluidField()  ->PrepareTimeStep();
  AleField()    ->PrepareTimeStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::Update()
{
  PoroField()     ->Update();
  FluidField()    ->Update();
  AleField()      ->Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::PrepareOutput()
{
  PoroField()->PrepareOutput();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::Output()
{
  PoroField()     ->Output();
  FluidField()    ->Output();
  AleField()      ->Output();
  FluidField()    ->LiftDrag();
}

/*----------------------------------------------------------------------*/
/*                          Coupling Methods                            */
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::SmallFluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return smallcoupsf_->MasterToSlave(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::SmallStructToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return smallcoupsf_->SlaveToMaster(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::StructToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_->MasterToSlave(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::AleToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_->SlaveToMaster(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_->MasterToSlave(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::FluidToPorofluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return couple_porofluid_fluid_->SlaveToMaster(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::PorofluidToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return couple_porofluid_fluid_->MasterToSlave(iv);
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<  MonolithicBase -> Monolithic  >>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

FPSI::Monolithic::Monolithic(const Epetra_Comm& comm,
                             const Teuchos::ParameterList& fpsidynparams,
                             const Teuchos::ParameterList& poroelastdynparams)
    : MonolithicBase(comm,fpsidynparams,poroelastdynparams),
      printscreen_(true),
      printiter_(true),
      printerrfile_(true),
      errfile_(DRT::Problem::Instance()->ErrorFile()->Handle()),
      timer_(comm),
      isfirsttimestep_(true),
      islinesearch_(false),
      firstcall_(true)
{
 // Setup linear solver
 bool builtsolver = false;
 builtsolver = FPSI::Monolithic::SetupSolver();
 if(builtsolver == false)
   dserror("Setup of linear solver failed ... Roundhouse-Kick!");

 const Teuchos::ParameterList& sdynparams = DRT::Problem::Instance()->StructuralDynamicParams();
 solveradapttol_= (DRT::INPUT::IntegralValue<int>(sdynparams, "ADAPTCONV") == 1);
 solveradaptolbetter_ = (sdynparams.get<double> ("ADAPTCONV_BETTER"));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::SetupSystem()
{
  const int ndim = DRT::Problem::Instance()->NDim();

  ADAPTER::Coupling& couppff = PoroFluidFluidCoupling();
  ADAPTER::Coupling& coupsf  = StructureFluidCoupling();
  ADAPTER::Coupling& smallcoupsf  = SmallStructureFluidCoupling();
  ADAPTER::Coupling& coupsa  = StructureAleCoupling();
  ADAPTER::Coupling& coupfa  = FluidAleCoupling();

  const Epetra_Map* fluidnodemap     = FluidField() ->Discretization() ->NodeRowMap();
  const Epetra_Map* alenodemap       = AleField()   ->Discretization() ->NodeRowMap();

  // porous fluid to fluid
  couppff.SetupConditionCoupling(*PoroField()->FluidField()->Discretization(),
                                 PoroField() ->FluidField()->FPSIInterface()->FSICondMap(),
                                *FluidField()->Discretization(),
                                 FluidField()->FPSIInterface()->FSICondMap(),
                                "FSICoupling",
                                 ndim+1,
                                 false);

  // porous structure to fluid
  coupsf.SetupConditionCoupling(*PoroField()->StructureField()->Discretization(),
                                 PoroField() ->StructureField()->Interface()->FSICondMap(),
                                *FluidField()->Discretization(),
                                 FluidField()->Interface()->FSICondMap(),
                                "FSICoupling",
                                 ndim,
                                 false);

  smallcoupsf.SetupConditionCoupling(
                                *FluidField()->Discretization(),
                                 FluidField()->Interface()->FSICondMap(),
                                *PoroField()->StructureField()->Discretization(),
                                 PoroField() ->StructureField()->Interface()->FSICondMap(),
                                "FSICoupling",
                                 ndim,
                                 false);

   // porous structure to ale
  coupsa.SetupConditionCoupling(*PoroField()->StructureField()->Discretization(),
                                 PoroField() ->StructureField()->Interface()->FSICondMap(),
                                *AleField()  ->Discretization(),
                                 AleField()  ->Interface()->FSICondMap(),
                                 "FSICoupling",
                                 ndim,
                                 false);

  // fluid to ale
  coupfa.SetupCoupling(*FluidField()-> Discretization(),
                       *AleField()  -> Discretization(),
                       *fluidnodemap,
                       *alenodemap,
                        ndim,
                        false);

  FluidField()->SetMeshMap(coupfa.MasterDofMap());

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::Timeloop()
{
  PrepareTimeloop();

  while (NotFinished()) // while step < maxsteps and time < maxtime
  {
    PrepareTimeStep();
    SetupNewton();
    TimeStep();
    PrepareOutput();
    Update();
    Output();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::PrepareTimeloop()
{
  // check if maps were destroyed before entring the timeloop
  Extractor().CheckForValidMapExtractor();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FPSI::Monolithic::TimeStep()
{
  //////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////                                                                ///////////////
  ///////////////                                 LOOP                           ///////////////
  ///////////////                                                                ///////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  while ((((not Converged()) and (iter_ <= maximumiterations_)) or (iter_ <= minimumiterations_)) or islinesearch_ == true)
  {
    // start time measurement
    timer_.ResetStartTime();
    Epetra_Time timer(Comm());

    Evaluate(iterinc_);

    // create full monolithic FPSI right-hand-side vector
    // moved to evaluate()

    // create full monolithic FPSI tangent stiffness matrix and check if it is filled
    SetupSystemMatrix();
    if (not systemmatrix_->Filled())
    {
      dserror("Effective tangent matrix must be filled here !");
    }

    // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
    // is done in PrepareSystemForNewtonSolve() within Evaluate(iterinc_)
    LinearSolve();

    //build norms
    BuildConvergenceNorms();

    // print stuff
    if (islinesearch_ == false)
      PrintNewtonIter();

    // reset solver tolerance
    solver_->ResetTolerance();

    // increment equilibrium loop index
    if (islinesearch_ == false)
    {
      iter_ += 1;
      PoroField()->IncrementPoroIter();
    }

  }// end loop

  // correct iteration counter
    iter_ -= 1;

    // test whether max iterations was hit
    if ( (Converged()) and (Comm().MyPID()==0) )
    {
      if (linesearch_counter > 0.5)
        std::cout<<"            Evaluation of residual with scaled increment yields: "<<normofrhs_<<std::endl;
      islinesearch_ = false;
      linesearch_counter=0.;
    }
    else if (iter_ >= maximumiterations_)
    {
      dserror("Newton found no convergence in %d iterations", iter_);
    }

    PoroField()->RecoverLagrangeMultiplier();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FPSI::Monolithic::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::Evaluate");

  if (linesearch_ and islinesearch_ == false)
  {
    FluidField()->Discretization()->ClearState();
    linesearch_counter=0.;
    FluidField()->Discretization()->SetState(0,"dispnp",FluidField()->Dispnp());
    meshdispold_ = AleToFluid(AleField()->WriteAccessDispnp());
    porointerfacedisplacementsold_ = StructToAle(PoroField() -> StructureField() -> ExtractInterfaceDispnp());
  }


  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> ax;

  if (x!=Teuchos::null)
  {
    ExtractFieldVectors(x,sx,fx,ax);
  }
  else
  {
    dserror("No existing increment vector !");
  }

  PoroField() -> Evaluate(sx);
  PoroField() -> UpdatePoroIterinc(sx);

  Teuchos::RCP<Epetra_Vector> porointerfacedisplacements = StructToAle(PoroField() -> StructureField() -> ExtractInterfaceDispnp());
  AleField()  -> ApplyInterfaceDisplacements(porointerfacedisplacements);
  Teuchos::RCP<Epetra_Vector> fullax = Teuchos::rcp(new Epetra_Vector(*AleField()->DofRowMap(),true));
  AleField()->Interface()->InsertOtherVector(ax,fullax);
  AleField()  -> Evaluate(fullax,"iter");

  Teuchos::RCP<Epetra_Vector> aledisplacements = AleToFluid(AleField()->WriteAccessDispnp());
  FluidField()->ApplyMeshDisplacement(aledisplacements);
  FluidField()->UpdateNewton(fx);
  FluidField()->Evaluate(Teuchos::null);

  SetupRHS(iter_ == 1);

}// Evaluate

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem::Instance()->AddFieldTest(PoroField()->StructureField()->CreateFieldTest());
  DRT::Problem::Instance()->AddFieldTest(PoroField()->FluidField()->CreateFieldTest());
  DRT::Problem::Instance()->AddFieldTest(FluidField()->CreateFieldTest());
  DRT::Problem::Instance()->TestAll(comm);
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<              Methods concerning          >>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<                    solver                >>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

bool FPSI::Monolithic::SetupSolver()
{

#ifdef FPSIDIRECTSOLVE
    const Teuchos::ParameterList& fpsidynamicparams = DRT::Problem::Instance()->FPSIDynamicParams();

    const int linsolvernumber = fpsidynamicparams.get<int>("LINEAR_SOLVER");
    if (linsolvernumber == (-1))
      dserror("No linear solver defined for FPSI problem. Please set LINEAR_SOLVER in FPSI DYNAMIC to a valid number !");

    const Teuchos::ParameterList& solverparams =
        DRT::Problem::Instance()->SolverParams(linsolvernumber);
    const int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(
        solverparams, "SOLVER");
    if (solvertype != INPAR::SOLVER::umfpack)
      dserror("umfpack solver expected but received any other type !");

    solver_ = Teuchos::rcp(new LINALG::Solver(solverparams,
                                              Comm(),
                                              DRT::Problem::Instance()->ErrorFile()->Handle())
                                             );
#else
    dserror("Only direct solver implemented so far !");
#endif //#ifdef FPSIDIRECTSOLVE

    // Get the parameters for the Newton iteration
    maximumiterations_ = fpsidynamicparams.get<int> ("ITEMAX");
    minimumiterations_ = fpsidynamicparams.get<int> ("ITEMIN");
    normtypeinc_ = DRT::INPUT::IntegralValue<INPAR::FPSI::ConvergenceNorm>(
        fpsidynamicparams, "NORM_INC");
    normtypefres_ = DRT::INPUT::IntegralValue<INPAR::FPSI::ConvergenceNorm>(
        fpsidynamicparams, "NORM_RESF");
    combinedconvergence_ = DRT::INPUT::IntegralValue<INPAR::FPSI::BinaryOp>(
        fpsidynamicparams, "NORMCOMBI_RESFINC");

    toleranceiterinc_        = fpsidynamicparams.get<double> ("INCTOL");
    toleranceresidualforces_ = fpsidynamicparams.get<double> ("RESTOL");

    DRT::Problem* problem = DRT::Problem::Instance();
    const Teuchos::ParameterList& fpsidynparams = problem->FPSIDynamicParams();
    linesearch_ = DRT::INPUT::IntegralValue<int>(fpsidynparams,"LineSearch");
    if(linesearch_==1)
      dserror("Parameter 'LineSearch' is set to 'Yes' in the FPSI Dynamic section in your input-file.  \n"
              "Though the framework for a line search algorithm is implemented in fpsi_monolithic.cpp, \n"
              "a proper routine to reset the participating single fields is still required. In Chuck's \n"
              "experimental baci this was solved by performing an evaluate with the negative increment.\n"
              "However this has not yet been committed.\n");
    linesearch_counter = 0.;

    return true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FPSI::Monolithic::LinearSolve()
{

  if (solveradapttol_ and (iter_ > 1))
  {
    double worst = normofrhs_;
    double wanted = toleranceresidualforces_;
    solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
  }
  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fpsidynparams = problem->FPSIDynamicParams();
  if(Teuchos::getIntegralValue<int>(fpsidynparams,"FDCheck"))
  {
    FPSIFDCheck();
  }

  Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more
  PoroField()->ClearPoroIterinc();

  LINALG::ApplyDirichlettoSystem(
      sparse,
      iterinc_,
      rhs_,
      Teuchos::null,
      zeros_,
      *CombinedDBCMap()
      );

    ///////////////////////////////////////////////////////
    // begin linesearch_      ///////////////////////////////
    ///////////////////////////////////////////////////////
    if (linesearch_ and iter_ > 1)
    {
      rhs_ -> Norm2(&normofrhs_);
      if (normofrhs_ - normofrhsold_ > 1e-13)
      {
        if (linesearch_counter > 0.5)
          std::cout<<"            Evaluation of residual with bisected increment yields: "<<normofrhs_<<std::endl;

        islinesearch_ = true;
        iterinc_ -> Update(pow(0.5,(linesearch_counter)),*iterincold_,0.0);
        linesearch_counter = linesearch_counter + 1.0;
        std::cout<<"linesearch_ : "<<std::setprecision(1)<<static_cast<int>(linesearch_counter+0.5)<<" iterinc_ multiplied by "<<std::setprecision(4)<<pow(0.5,linesearch_counter)<<"   residual = "<<normofrhs_<<" > "<<normofrhsold_<<std::endl;

        // substract the old interinc_ from all fields (undo the update)
        Teuchos::RCP<Epetra_Vector> sx;
        Teuchos::RCP<Epetra_Vector> fx;
        Teuchos::RCP<const Epetra_Vector> constsx;
        Teuchos::RCP<const Epetra_Vector> constfx;
        Teuchos::RCP<const Epetra_Vector> ax;

        sx = Teuchos::rcp(new Epetra_Vector(*PoroField()->DofRowMap(), true));
        fx = Teuchos::rcp(new Epetra_Vector(*FluidField()->DofRowMap(), true));

        ExtractFieldVectors(iterinc_,constsx,constfx,ax);
        iterinc_ -> Norm2(&normofiterinc_);
        std::cout<<"            Norm of step back: "<<normofiterinc_<<std::endl;
        //PoroField()  ->ResetNewton(sx);
        //FluidField() ->ResetNewton(fx);
        sx -> Update(1.0,*constsx,0.0);
        fx -> Update(1.0,*constfx,0.0);
        sx -> Scale(-1.0);
        fx -> Scale(-1.0);
        PoroField()  ->Evaluate(sx);
        FluidField() ->UpdateNewton(Teuchos::rcp_dynamic_cast<const Epetra_Vector>(fx));
        //AleField()   ->ResetNewton(ax);

        FluidField()->ApplyMeshDisplacement(meshdispold_);
        AleField()  ->ApplyInterfaceDisplacements(porointerfacedisplacementsold_);


        // set iterinc_ to a fraction of the old iterinc_
        iterinc_ -> Update(pow(0.5,linesearch_counter),*iterincold_,0.0);
        iterinc_ -> Norm2(&normofiterinc_);
        std::cout<<"            Norm of old increment: "<<normofiterincold_<<"  Norm of bisected increment: "<<normofiterinc_<<std::endl;

      }
      else
      {
        islinesearch_ = false;
        if (linesearch_counter > 0.5)
          std::cout<<"            Evaluation of residual with bisected increment yields: "<<normofrhs_<<std::endl;
        linesearch_counter = 0.0;
      }
    }
    // end linesearch_

    // prepare linesearch_
    // copy the old iterinc_ before new solve
    if (linesearch_ and islinesearch_ == false)
    {
      rhsold_ =  LINALG::CreateVector(*DofRowMap(), true);
      rhsold_ -> Update(1.0,*rhs_,0.0);
      rhsold_ -> Norm2(&normofrhsold_);
      if (abs(normofrhs_ - normofrhsold_) > 1.e-12 and iter_>1)
        dserror(" wrong copy of rhs_ ");
    }
    // end prepare linesearch_


//  std::cout<<"iterinc_ before solve(): \n"<<*iterinc_<<std::endl;
  // standard solver call
  if (islinesearch_ == false)
  {
    solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, true, iter_ == 1);
  }
//  std::cout<<"iterinc_ after solve(): \n"<<*iterinc_<<std::endl;
//  std::cout<<"Full Monolithic RHS of FPSI Problem : \n"<<*rhs_<<std::endl;
//  rhs_->Norm2(&normofrhs_);
//  std::cout<<"Norm of RHS: "<<normofrhs_<<std::endl;

  if (islinesearch_ ==  false)
  {
    // check whether iterinc_ points in right direction
    Teuchos::RCP<Epetra_Vector> tempvec =  LINALG::CreateVector(*DofRowMap(), true);
    sparse -> Multiply(true,*rhs_,*tempvec);
    double climb = 0.0;
    tempvec -> Dot(*iterinc_,&climb);
    climb = -climb;

    if(climb > 0.0)
    {
      std::cout<<"########################################################################"<<std::endl;
      std::cout<<"##                                                                    ##"<<std::endl;
      std::cout<<"## WARNING: A*x-b=0 ; A^T*b*x > 0 ; increment vector multiplied by -1 ##"<<std::endl;
      std::cout<<"##                                                                    ##"<<std::endl;
      std::cout<<"##                       Value = "<<std::setprecision(9)<<climb<<"    ##"<<std::endl;
      std::cout<<"##                                                                    ##"<<std::endl;
      std::cout<<"########################################################################"<<std::endl;
      iterinc_ -> Update(-1.0,*iterinc_,0.0);
    }
  }

  if (linesearch_ and islinesearch_ == false)
      {
        iterincold_ =  LINALG::CreateVector(*DofRowMap(), true);
        iterincold_ -> Update(1.0,*iterinc_,0.0);
        iterincold_ -> Norm2(&normofiterincold_);
      }

}

Teuchos::RCP<Epetra_Map> FPSI::Monolithic::CombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map> scondmap =
      PoroField()->StructureField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> pfcondmap =
       PoroField()->FluidField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> fcondmap =
                       FluidField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> acondmap =
                         AleField()->GetDBCMapExtractor()->CondMap();
  Teuchos::RCP<Epetra_Map> tempmap =
      LINALG::MergeMap(scondmap, pfcondmap, false);
  Teuchos::RCP<Epetra_Map> condmap_0 =
        LINALG::MergeMap(tempmap, fcondmap, false);
  Teuchos::RCP<Epetra_Map> condmap =
       LINALG::MergeMap(condmap_0, acondmap, false);

  return condmap;
} // CombinedDBCMap()

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<              Matrix operations           >>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<                                          >>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Teuchos::RCP<LINALG::SparseMatrix> FPSI::Monolithic::PoroFluidCouplingMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_);
  //return k_pf_->Merge();
}

void FPSI::Monolithic::ApplyCouplingTerms(Teuchos::RCP<LINALG::SparseOperator>   k_pf,
                                            Teuchos::RCP<LINALG::SparseOperator> k_fp,
                                            Teuchos::RCP<LINALG::SparseMatrix>   p,
                                            Teuchos::RCP<LINALG::SparseMatrix>   f,
                                            Teuchos::RCP<LINALG::SparseOperator> a,
                                            Teuchos::RCP<LINALG::SparseOperator> k_fa
                                           )
{
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::ApplyCouplingTerms");

  Teuchos::RCP<LINALG::SparseMatrix> k_pf_porofluid = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid_);
  Teuchos::RCP<LINALG::SparseMatrix> k_fp_porofluid = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_fp_porofluid_);

  a              ->Zero();
  k_pf           ->Zero();
  k_fp           ->Zero();
  //k_pf_struct    ->Zero();
  k_pf_porofluid ->Zero();

  const ADAPTER::Coupling& couppff      = PoroFluidFluidCoupling();
  const ADAPTER::Coupling& coupsf       = StructureFluidCoupling();
//  const ADAPTER::Coupling& smallcoupsf  = SmallStructureFluidCoupling();

  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fpsidynparams       = problem->FPSIDynamicParams();
  INPAR::FPSI::PartitionedCouplingMethod method = DRT::INPUT::IntegralValue<INPAR::FPSI::PartitionedCouplingMethod>(fpsidynparams,"PARTITIONED");

  if(method != INPAR::FPSI::nocoupling)
  {

  // set general vector values needed by elements

   PoroField()->FluidField()->Discretization()->ClearState();

   PoroField()->FluidField()->Discretization()->SetState(0,"dispnp" ,PoroField()->FluidField()->Dispnp());
   PoroField()->FluidField()->Discretization()->SetState(0,"gridv"  ,PoroField()->FluidField()->GridVel());
   PoroField()->FluidField()->Discretization()->SetState(0,"dispn"  ,PoroField()->FluidField()->Dispn());
   PoroField()->FluidField()->Discretization()->SetState(0,"veln"   ,PoroField()->FluidField()->Veln());
   PoroField()->FluidField()->Discretization()->SetState(0,"velaf"  ,PoroField()->FluidField()->Velnp());
   PoroField()->FluidField()->Discretization()->SetState(0,"velnp"  ,PoroField()->FluidField()->Velnp());

   FluidField()->Discretization()->ClearState();

   FluidField()->Discretization()->SetState(0,"dispnp",FluidField()->Dispnp());
   FluidField()->Discretization()->SetState(0,"gridv" ,FluidField()->GridVel());
   FluidField()->Discretization()->SetState(0,"dispn" ,FluidField()->Dispn());
   FluidField()->Discretization()->SetState(0,"veln"  ,FluidField()->Veln());
   FluidField()->Discretization()->SetState(0,"velaf" ,FluidField()->Velnp());
   FluidField()->Discretization()->SetState(0,"velnp" ,FluidField()->Velnp());


  // create the parameters for the discretization
  Teuchos::ParameterList fparams;

  // action for elements
  fparams.set<int>("action", FLD::fpsi_coupling);
  fparams.set("timescale",PoroField()->FluidField()->ResidualScaling());
  fparams.set("timestep",Step());
  fparams.set("iter",iter_);
  fparams.set("dt",fpsidynparams.get<double>("TIMESTEP"));
  //std::cout<<"ResidualScaling() = "<<FluidField()->ResidualScaling()<<std::endl;
  //fparams.set("timescale",PoroField()->PoroFluidField()->TimeScaling());

  if (method == INPAR::FPSI::monolithic)
  {

  fparams.set<string>("fillblock","Porofluid_Freefluid");
  fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
  DRT::AssembleStrategy fluidstrategy(
      0,              // fluiddofset for row
      0,              // fluiddofset for column
      k_pf_porofluid, // coupling matrix with fluid rowmap
      Teuchos::null,  // no other matrix or vectors
      Teuchos::null ,
      Teuchos::null,
      Teuchos::null
  );

  FluidField()->Discretization()->EvaluateCondition( fparams, fluidstrategy,"FSICoupling" );
  k_pf_porofluid -> Complete(f->DomainMap(),f->RangeMap());
  Teuchos::RCP<LINALG::SparseMatrix> temp2 = Teuchos::rcp(new LINALG::SparseMatrix(*(PoroField()->FluidField()->DofRowMap()),81,false)); // temp matrix with porofluid rowmap
  {
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
  (*couplingrowtransform_)(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid),
                            1.0,
                            ADAPTER::CouplingSlaveConverter(couppff),
                           *temp2,
                            false);
  }

  k_pf_porofluid -> Zero();
  k_pf          -> UnComplete();
  temp2          -> Complete(f->DomainMap(),*PoroField()->FluidRangeMap());
  k_pf->Add(*temp2,         false,1.0,1.0);
  k_pf->Complete(f->DomainMap(),p->RangeMap());

  fparams.set<string>("fillblock","Porofluid_Structure");
  fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
  k_pf_porofluid_ -> Zero();
  DRT::AssembleStrategy fluidstrategy21(
      0,                   // porofluiddofset for row
      0,                   // structuredofset for column
      k_pf_porofluid_,     // coupling matrix with fluid rowmap
      Teuchos::null,       // no other matrix or vectors
      Teuchos::null ,
      Teuchos::null,
      Teuchos::null
  );

  FluidField()    -> Discretization()->EvaluateCondition( fparams, fluidstrategy21, "FSICoupling" );
  k_pf_porofluid_ -> Complete(f->DomainMap(),f->RangeMap());

  Teuchos::RCP<LINALG::SparseOperator> temp51 = Teuchos::rcp(new LINALG::SparseMatrix((p->RowMap()),81,false));
  {
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
  (*couplingrowcoltransform2_)(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid_),
      1.0,
      ADAPTER::CouplingSlaveConverter(couppff), // row converter: important to use slave converter
      ADAPTER::CouplingSlaveConverter(coupsf), //  col converter: important to use slave converter
      *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(temp51),
      false); // bool exactmatch = true (default)
  }

  temp51  -> Complete(p->DomainMap(),p->RangeMap());
  p->Add(*temp51,false,1.0,1.0);

  fparams.set<string>("fillblock","Fluid_Porofluid");
  fparams.set("InterfaceFacingElementMap", PoroFluid_Fluid_InterfaceMap);
  k_fp_porofluid->Zero();
  DRT::AssembleStrategy porofluidstrategy(
      0,                   // porofluiddofset for row
      0,                   // porofluiddofset for column
      k_fp_porofluid,      // porofluid-structure coupling matrix
      Teuchos::null ,      // no other matrix or vectors
      Teuchos::null ,
      Teuchos::null ,
      Teuchos::null
  );

  PoroField()     -> FluidField()->Discretization()->EvaluateCondition( fparams, porofluidstrategy, "FSICoupling" );
  k_fp_porofluid  -> Complete(*PoroField()->FluidDomainMap(),*PoroField()->FluidRangeMap());
  //std::cout<<"\n k_fp_porofluid matrix: \n"<<*k_fp_porofluid<<std::endl;
  Teuchos::RCP<LINALG::SparseMatrix> temp = Teuchos::rcp(new LINALG::SparseMatrix((f->RowMap()),81,false));
  {
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
  (*couplingrowtransform2_)(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_fp_porofluid),
                             1.0,
                             ADAPTER::CouplingMasterConverter(couppff),
                            *temp,
                             true);
  }

  temp  -> Complete(*PoroField()->FluidRangeMap(),f->RangeMap());

  k_fp -> UnComplete();
  k_fp -> Add(*temp,false,1.0,1.0);
  k_fp -> Complete(p->DomainMap(),f->RangeMap());

  fparams.set<string>("fillblock","Fluid_Structure");
  fparams.set("InterfaceFacingElementMap", PoroFluid_Fluid_InterfaceMap);
  k_pfs_ -> Zero();
  k_pfs_ -> UnComplete();

  DRT::AssembleStrategy structurestrategy(
      0,                   // porofluiddofset for row
      1,                   // structuredofset for column
      k_pfs_,              // coupling matrix with porofluid rowmap
      Teuchos::null,       // no other matrix or vectors
      Teuchos::null,
      Teuchos::null,
      Teuchos::null
  );

  PoroField() -> FluidField()->Discretization()->EvaluateCondition( fparams, structurestrategy, "FSICoupling" );
  k_pfs_      -> Complete(*PoroField()->StructureRangeMap(),*PoroField()->FluidRangeMap());
  Teuchos::RCP<LINALG::SparseMatrix> temp3 = Teuchos::rcp(new LINALG::SparseMatrix((f->RowMap()),81,false));
  {
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
  (*couplingrowtransform3_)(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pfs_),
                             1.0,
                             ADAPTER::CouplingMasterConverter(couppff),
                            *temp3,
                             true);
  }

  temp3  -> Complete(*PoroField()->StructureRangeMap(),f->RangeMap());
  //std::cout<<"temp3 matrix: \n"<<*temp3<<std::endl;

  k_fp -> UnComplete();
  k_fp -> Add(*temp3,false,1.0,1.0);
  k_fp -> Complete(p->RangeMap(),f->RangeMap());

  ///// Fluid_Structure (fluid part / linearization of tangentials with respect to displacements)
  k_pf_porofluid_ -> Zero();
  k_pf_porofluid_ -> UnComplete();

  DRT::AssembleStrategy structurestrategy2(
      0,                   // porofluiddofset for row
      0,                   // fluiddofset for column
      k_pf_porofluid_,     // coupling matrix with porofluid rowmap
      Teuchos::null,       // no other matrix or vectors
      Teuchos::null,
      Teuchos::null,
      Teuchos::null
   );

  fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
  FluidField()->Discretization()->EvaluateCondition( fparams, structurestrategy2, "FSICoupling" );

  k_pf_porofluid_ -> Complete(f->RangeMap(),f->RangeMap());
  Teuchos::RCP<LINALG::SparseMatrix> temp32 = Teuchos::rcp(new LINALG::SparseMatrix((f->RowMap()),81,false));
  {
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
  (*couplingcoltransform_)( f->RowMap(),
                            f->ColMap(),
                            *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid_),
                             1.0,
                             ADAPTER::CouplingSlaveConverter(coupsf), // row converter: important to use slave converter
                            *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(temp32),
                             false); // bool exactmatch = true (default)
  }

  temp32  -> Complete(*PoroField()->StructureRangeMap(),f->RangeMap());
  k_fp -> UnComplete();
  k_fp -> Add(*temp32,false,1.0,1.0);
  k_fp -> Complete(p->RangeMap(),f->RangeMap());

  fparams.set<string>("fillblock","Fluid_Fluid");
  fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
    k_pf_porofluid_ -> Zero();
    k_pf_porofluid_ -> UnComplete();

    DRT::AssembleStrategy fluidfluidstrategy(
        0,                   // fluiddofset for row
        0,                   // fluiddofset for column
        k_pf_porofluid_,     // porofluid-structure coupling matrix
        Teuchos::null,       // no other matrix or vectors
        Teuchos::null,
        Teuchos::null,
        Teuchos::null
    );

    FluidField()   ->Discretization()->EvaluateCondition( fparams, fluidfluidstrategy, "FSICoupling" );
    k_pf_porofluid -> Complete(f->RangeMap(),f->RangeMap());

    f -> UnComplete();
    f -> Add(*k_pf_porofluid_,false,1.0,1.0);
    f -> Complete(f->RangeMap(),f->RangeMap());

    fparams.set<string>("fillblock","Structure_Fluid");
    fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
    k_pf_porofluid_ -> Zero();
    k_pf_porofluid_ -> UnComplete();

    DRT::AssembleStrategy structurefluidstrategy(
        0,                   // fluid dofset for row
        0,                   // fluid dofset for column
        k_pf_porofluid_ ,    // coupling matrix with fluid rowmap
        Teuchos::null,       // no other matrix or vectors
        Teuchos::null,
        Teuchos::null,
        Teuchos::null
    );

    FluidField()   ->Discretization()->EvaluateCondition( fparams, structurefluidstrategy, "FSICoupling" );

    k_pf_porofluid_->Complete(f->DomainMap(), f->RangeMap());
    Teuchos::RCP<LINALG::SparseMatrix> temp4 = Teuchos::rcp(new LINALG::SparseMatrix((*PoroField()->StructureField()->DofRowMap()),81,false));
    {
    TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
    (*couplingrowtransform4_)(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid_),
                               1.0,
                               ADAPTER::CouplingSlaveConverter(coupsf), // important to use slave converter
                              *temp4,
                               false);
    }

    temp4  -> Complete(f->DomainMap(),*PoroField()->StructureRangeMap());

    k_pf -> UnComplete();
    k_pf -> Add(*temp4,false,1.0,1.0);
    k_pf -> Complete(f->DomainMap(), p->RangeMap());

    fparams.set<string>("fillblock","Structure_Structure");
    fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
        k_pf_porofluid_ -> Zero();
        k_pf_porofluid_ -> UnComplete();

        DRT::AssembleStrategy structurestructurestrategy(
            0,                   // fluid dofset for row
            0,                   // fluid dofset for column
            k_pf_porofluid_,     // coupling matrix with fluid rowmap
            Teuchos::null,       // no other matrix or vectors
            Teuchos::null,
            Teuchos::null,
            Teuchos::null
        );

        FluidField()   ->Discretization()->EvaluateCondition( fparams, structurestructurestrategy, "FSICoupling" );

    // condense linearization with respect to the ale mesh motion (interface structural displacements = interface ale displacements)
        k_pf_porofluid_->Complete(f->DomainMap(), f->RangeMap());
        //std::cout<<*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid_)<<std::endl;
        Teuchos::RCP<LINALG::SparseOperator> temp5 = Teuchos::rcp(new LINALG::SparseMatrix((p->RowMap()),81,false));
        {
        TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
        (*couplingrowcoltransform_)(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid_),
                                   1.0,
                                   ADAPTER::CouplingSlaveConverter(coupsf), // row converter: important to use slave converter
                                   ADAPTER::CouplingSlaveConverter(coupsf), // col converter: important to use slave converter
                                  *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(temp5),
                                   false); // bool exactmatch = true (default)
        }

        temp5  -> Complete(p->DomainMap(),p->RangeMap());
        p      -> Add(*temp5,false,1.0,1.0);

    // Process inner ale dofs
        fparams.set<string>("fillblock","Structure_Ale");
        fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
        a -> Zero();
        a -> UnComplete();

        DRT::AssembleStrategy structurealestrategy(
            0,                   // fluid dofset for row
            1,                   // ale dofset for column
            a,                   // coupling matrix with fluid rowmap
            Teuchos::null,       // no other matrix or vectors
            Teuchos::null,
            Teuchos::null,
            Teuchos::null
        );

        FluidField()   ->Discretization()->EvaluateCondition( fparams, structurealestrategy, "FSICoupling" );

  } // if monolithic
  else
  {
    if (isfirsttimestep_)
    std::cout<<"LINEARIZATION OF FPSI INTERFACE IS TURNED OFF \n"
               "IF YOU THINK THAT SUCKS. SET 'PARITIONED' to 'monolithic' \n"
               "IN THE FPSI SECTION OF YOUR INPUT FILE !!! "<<std::endl;
    isfirsttimestep_ = false;
  }


        //////////////////////////////////////
        ///////      __       ___       //////
        ///////      |_|  |_| |__       //////
        ///////      | \  | |  __|      //////
        //////                          //////
        //////////////////////////////////////

        Teuchos::RCP<Epetra_Vector> temprhs  = Teuchos::null;
        Teuchos::RCP<Epetra_Vector> temprhs2 = Teuchos::null;;

        fparams.set<string>("fillblock","conti");
        fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
        temprhs  = Teuchos::rcp(new Epetra_Vector(*FluidField()->DofRowMap(),true));
        temprhs2 = Teuchos::rcp(new Epetra_Vector(*PoroField() ->Interface().FullMap(),true));
        temprhs ->PutScalar(0.0);
        temprhs2->PutScalar(0.0);

        DRT::AssembleStrategy rhscontistrategy(
            0,                   // fluid dofset for row
            0,                   // fluid dofset for column
            Teuchos::null,
            Teuchos::null,
            temprhs,             // rhs vector
            Teuchos::null,
            Teuchos::null
        );

        FluidField()   ->Discretization()->EvaluateCondition( fparams, rhscontistrategy, "FSICoupling" );

        // extract FSI part of the fluid field
        temprhs = FluidField()->FPSIInterface()->ExtractFSICondVector(temprhs);
        // replace global fluid interface dofs through porofluid interface dofs
        temprhs = FluidToPorofluid(temprhs);
        // insert porofluid interface entries into vector with full porofield length (0: inner dofs of structure, 1: interface dofs of structure, 2: inner dofs of porofluid, 3: interface dofs of porofluid )
        PoroField()->Interface().InsertVector(temprhs,3,temprhs2);
        // add vector with full porofield length to global rhs
        Extractor().AddVector(temprhs2,0,rhs_);

        temprhs ->PutScalar(0.0);
        temprhs2->PutScalar(0.0);

        fparams.set<string>("fillblock","structure");
        fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
        temprhs  = Teuchos::rcp(new Epetra_Vector(*FluidField()->Interface()->FullMap(),true)); //
        temprhs2 = Teuchos::rcp(new Epetra_Vector(*PoroField() ->Interface().FullMap(),true));

        DRT::AssembleStrategy rhsstructurestrategy(
            0,                   // fluid dofset for row
            0,                   // fluid dofset for column
            Teuchos::null,
            Teuchos::null,
            temprhs,             // rhs vector
            Teuchos::null,
            Teuchos::null
        );

        FluidField()   ->Discretization()->EvaluateCondition( fparams, rhsstructurestrategy, "FSICoupling" );

        // extract FPSI part of the fluid field
        temprhs = FluidField()->Interface()->ExtractFSICondVector(temprhs); //
        // replace global fluid interface dofs through porofluid interface dofs
        temprhs = SmallFluidToStruct(temprhs); //
        // insert porofluid interface entries into vector with full porofield length (0: inner dofs of structure, 1: interface dofs of structure, 2: inner dofs of porofluid, 3: interface dofs of porofluid )
        PoroField()->Interface().InsertVector(temprhs,1,temprhs2);
        // add vector with full porofield length to global rhs
        Extractor().AddVector(temprhs2,0,rhs_);

        temprhs ->PutScalar(0.0);
        temprhs2->PutScalar(0.0);

        fparams.set<string>("fillblock","fluid");
        fparams.set("InterfaceFacingElementMap", PoroFluid_Fluid_InterfaceMap);
        temprhs  = Teuchos::rcp(new Epetra_Vector(*PoroField()->FluidField()->FPSIInterface()->FullMap(),true)); //
        temprhs2 = Teuchos::rcp(new Epetra_Vector(*FluidField()->FPSIInterface()->FullMap(),true));

        DRT::AssembleStrategy rhsfluidstrategy(
            0,                   // fluid dofset for row
            0,                   // fluid dofset for column
            Teuchos::null,
            Teuchos::null,
            temprhs,             // rhs vector
            Teuchos::null,
            Teuchos::null
        );

        PoroField()->FluidField()->Discretization()->EvaluateCondition( fparams, rhsfluidstrategy, "FSICoupling" );

        // extract FSI part of the fluid field
        temprhs = PoroField()->FluidField()->FPSIInterface()->ExtractFSICondVector(temprhs); //
        //std::cout<<*temprhs<<std::endl;
        // replace global fluid interface dofs through porofluid interface dofs
        temprhs = PorofluidToFluid(temprhs); //
        // insert porofluid interface entries into vector with full fluidfield length (0: inner dofs of fluid, 1: interface dofs of fluid )
        FluidField()->FPSIInterface()->InsertVector(temprhs,1,temprhs2);
        // add vector with full porofield length to global rhs
        Extractor().AddVector(temprhs2,1,rhs_);

        temprhs ->PutScalar(0.0);
        temprhs2->PutScalar(0.0);


        fparams.set<string>("fillblock","fluidfluid"); // (wot,tangentialfac*uot) part
        fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
        temprhs  = Teuchos::rcp(new Epetra_Vector(*FluidField()->FPSIInterface()->FullMap(),true)); //


        DRT::AssembleStrategy rhsfluidfluidstrategy(
            0,                   // fluid dofset for row
            0,                   // fluid dofset for column
            Teuchos::null,
            Teuchos::null,
            temprhs,             // rhs vector
            Teuchos::null,
            Teuchos::null
        );

        FluidField()->Discretization()->EvaluateCondition( fparams, rhsfluidfluidstrategy, "FSICoupling" );

        Extractor().AddVector(temprhs,1,rhs_);

        temprhs ->PutScalar(0.0);

        //////////////////////////////////////
        ///////                         //////
        //////    NEUMANN INTEGRATION   //////
        //////                          //////
        //////////////////////////////////////
        if (method == INPAR::FPSI::monolithic or method == INPAR::FPSI::RobinNeumann)
        {
          fparams.set<string>("fillblock","NeumannIntegration");
          fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
          k_pf_porofluid_ -> Zero();
          k_pf_porofluid_ -> UnComplete();
          DRT::AssembleStrategy rhsfluidfluidstrategy2(
              0,                   // fluid dofset for row
              0,                   // fluid dofset for column
              k_pf_porofluid_ ,    // coupling matrix with fluid rowmap
              Teuchos::null,
              temprhs,             // rhs vector
              Teuchos::null,
              Teuchos::null
          );
          temprhs->PutScalar(0.0);
          FluidField()->Discretization()->EvaluateCondition( fparams, rhsfluidfluidstrategy2, "NeumannIntegration" );
          //std::cout<<"temprhs : \n"<<*temprhs<<std::endl;
          Extractor().AddVector(temprhs,1,rhs_);
          k_pf_porofluid_ -> Complete(f->RangeMap(),f->RangeMap());
          f        -> Add(*k_pf_porofluid_ ,false,1.0,1.0);


          fparams.set<string>("fillblock","NeumannIntegration_Ale");
          fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
          k_fa -> Zero();
          k_fa -> UnComplete();
          DRT::AssembleStrategy rhsfluidfluidstrategy3(
              0,                   // fluid dofset for row
              1,                   // ale dofset for column
              k_fa,                // coupling matrix with fluid rowmap
              Teuchos::null,
              Teuchos::null,
              Teuchos::null,
              Teuchos::null
          );
          FluidField()->Discretization()->EvaluateCondition( fparams, rhsfluidfluidstrategy3, "NeumannIntegration" );
          //std::cout<<"Fluid_ALE: \n"<<Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(k_fa_)->Matrix(1,1)<<std::endl;


          fparams.set<string>("fillblock","NeumannIntegration_Struct");
          fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
          k_pf_porofluid_ -> Zero();
          k_pf_porofluid_ -> UnComplete();
          DRT::AssembleStrategy rhsfluidfluidstrategy4(
              0,                   // fluid dofset for row
              0,                   // ale dofset for column
              k_pf_porofluid_,    // coupling matrix with fluid rowmap
              Teuchos::null,
              Teuchos::null,
              Teuchos::null,
              Teuchos::null
          );
           FluidField()->Discretization()->EvaluateCondition(fparams, rhsfluidfluidstrategy4, "NeumannIntegration" );
           k_pf_porofluid_ -> Complete(f->RangeMap(),f->RangeMap());

           Teuchos::RCP<LINALG::SparseMatrix> temp33 = Teuchos::rcp(new LINALG::SparseMatrix((f->RowMap()),81,false));
           (*couplingcoltransform_)( f->RowMap(),
                                     f->ColMap(),
                                     *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid_),
                                      1.0,
                                      ADAPTER::CouplingSlaveConverter(coupsf), // row converter: important to use slave converter
                                     *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(temp33),
                                      false); // bool exactmatch = true (default)

           temp33  -> Complete(*PoroField()->StructureRangeMap(),f->RangeMap());
           //std::cout<<"Fluid_Structure: \n"<<*temp33<<std::endl;
           k_fp -> UnComplete();
           k_fp -> Add(*temp33,false,1.0,1.0);
           k_fp -> Complete(p->RangeMap(),f->RangeMap());

        }
        else // only fill rhs
        {
          fparams.set<string>("fillblock","NeumannIntegration");
          fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap);
          temprhs->PutScalar(0.0);
          FluidField()->Discretization()->EvaluateCondition( fparams, rhsfluidfluidstrategy, "NeumannIntegration" );
          Extractor().AddVector(temprhs,1,rhs_);
        }

        ////////////////////////////
        // DONE -> CLEAR STATES
        PoroField() ->FluidField()->Discretization()->ClearState();
        FluidField()->Discretization()->ClearState();

  } // if not nocoupling
}
//--------------------------------------------------------------------------------//
//--------------------------------------------------------------------------------//

Teuchos::RCP<LINALG::SparseMatrix> FPSI::Monolithic::FluidPoroCouplingMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_fp_);
}

void FPSI::Monolithic::ApplyFluidCouplMatrix(Teuchos::RCP<LINALG::SparseOperator> k_fp)
{
  k_fp->Zero();
}
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<                 Newton Loop              >>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<                   Methods                >>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool FPSI::Monolithic::Converged()
{
  // check for single norms
  bool convinc  = false;
  bool convfres = false;

  // residual increments
  switch (normtypeinc_)
  {
    case INPAR::FPSI::absoluteconvergencenorm:
      convinc = normofiterinc_ < toleranceiterinc_;
      break;
    default:
      dserror("Cannot check for convergence of residual values for any reason :-p !");
      break;
  }

  // residual forces
  switch (normtypefres_)
  {
    case INPAR::FPSI::absoluteconvergencenorm:
    convfres = normofrhs_ < toleranceresidualforces_;
    break;
    default:
    dserror("Cannot check for convergence of residual forces for any reason :-P !");
    break;
  }

  // combine increments and forces
  bool converged = false;
  if (combinedconvergence_ == INPAR::FPSI::bop_and)
    converged = convinc and convfres;
  else if (combinedconvergence_== INPAR::FPSI::bop_or)
    converged = convinc or convfres;
  else
    dserror("Something went terribly wrong with binary operator!");

  // return things
  return converged;
}  // Converged()



void FPSI::Monolithic::BuildConvergenceNorms()
{
  rhs_->Norm2(&normofrhs_);
  Teuchos::RCP<const Epetra_Vector> rhs_poro;
  Teuchos::RCP<const Epetra_Vector> rhs_fluid;
  Teuchos::RCP<const Epetra_Vector> rhs_fluidvelocity;
  Teuchos::RCP<const Epetra_Vector> rhs_fluidpressure;
  Teuchos::RCP<const Epetra_Vector> rhs_porofluidvelocity;
  Teuchos::RCP<const Epetra_Vector> rhs_porofluidpressure;
  Teuchos::RCP<const Epetra_Vector> rhs_porointerface;
  Teuchos::RCP<const Epetra_Vector> rhs_fluidinterface;
  Teuchos::RCP<const Epetra_Vector> rhs_porofluid;
  Teuchos::RCP<const Epetra_Vector> rhs_porostruct;
  Teuchos::RCP<const Epetra_Vector> rhs_ale;

  rhs_poro              = Extractor().ExtractVector(rhs_, 0);
  rhs_porostruct = PoroField()->Extractor().ExtractVector(rhs_poro, 0);
  rhs_porofluid  = PoroField()->Extractor().ExtractVector(rhs_poro, 1);
  rhs_porofluidvelocity = PoroField()->FluidField()->ExtractVelocityPart(rhs_porofluid);
  rhs_porofluidpressure = PoroField()->FluidField()->ExtractPressurePart(rhs_porofluid);
  rhs_porointerface     = PoroField()->FluidField()->FPSIInterface()->ExtractFSICondVector(rhs_porofluid);

  rhs_fluid             = Extractor().ExtractVector(rhs_, 1);
  rhs_fluidvelocity     = FluidField()->ExtractVelocityPart(rhs_fluid);
  rhs_fluidpressure     = FluidField()->ExtractPressurePart(rhs_fluid);
  rhs_fluidinterface    = FluidField()->FPSIInterface()->ExtractFSICondVector(rhs_fluid);


  rhs_ale               = Extractor().ExtractVector(rhs_, 2);//Extractor().ExtractVector(rhs_, 2);

  rhs_porostruct        ->Norm2(&normrhsporostruct_);
  rhs_fluid             ->Norm2(&normrhsfluid_);
  rhs_fluidvelocity     ->Norm2(&normrhsfluidvelocity_);
  rhs_fluidpressure     ->Norm2(&normrhsfluidpressure_);
  rhs_porofluidvelocity ->Norm2(&normrhsporofluidvelocity_);
  rhs_porofluidpressure ->Norm2(&normrhsporofluidpressure_);
  rhs_porointerface     ->Norm2(&normrhsporointerface_);
  rhs_fluidinterface    ->Norm2(&normrhsfluidinterface_);
  rhs_ale               ->Norm2(&normrhsale_);

if (islinesearch_ == false)
  iterinc_->Norm2(&normofiterinc_);

  Teuchos::RCP<const Epetra_Vector> iterincporo;
  Teuchos::RCP<const Epetra_Vector> iterincporostruct;
  Teuchos::RCP<const Epetra_Vector> iterincporofluid;
  Teuchos::RCP<const Epetra_Vector> iterincfluid;
  Teuchos::RCP<const Epetra_Vector> iterincfluidvelocity;
  Teuchos::RCP<const Epetra_Vector> iterincfluidpressure;
  Teuchos::RCP<const Epetra_Vector> iterincporofluidvelocity;
  Teuchos::RCP<const Epetra_Vector> iterincporofluidpressure;
  Teuchos::RCP<const Epetra_Vector> iterincale;
  Teuchos::RCP<const Epetra_Vector> iterincporointerface;
  Teuchos::RCP<const Epetra_Vector> iterincfluidinterface;

  iterincporo              = Extractor().ExtractVector(iterinc_, 0);
  iterincporostruct        = PoroField()->Extractor().ExtractVector(iterincporo,0);
  iterincporofluid         = PoroField()->Extractor().ExtractVector(iterincporo,1);
  iterincporofluidvelocity = PoroField() ->FluidField()->ExtractVelocityPart(iterincporofluid);
  iterincporofluidpressure = PoroField() ->FluidField()->ExtractPressurePart(iterincporofluid);
  iterincporointerface     = PoroField() ->Interface().ExtractVector(iterincporo,3);

  iterincfluid             = Extractor().ExtractVector(iterinc_, 1);
  iterincfluidvelocity     = FluidField()->ExtractVelocityPart(iterincfluid);
  iterincfluidpressure     = FluidField()->ExtractPressurePart(iterincfluid);
  iterincfluidinterface    = FluidField()->FPSIInterface()->ExtractFSICondVector(iterincfluid);

  iterincale               = Extractor().ExtractVector(iterinc_, 2);

  iterincporo              ->Norm2(&normofiterincporo_);
  iterincporostruct        ->Norm2(&normofiterincporostruct_);
  iterincporofluid         ->Norm2(&normofiterincporofluid_);
  iterincporofluidvelocity ->Norm2(&normofiterincporofluidvelocity_);
  iterincporofluidpressure ->Norm2(&normofiterincporofluidpressure_);
  iterincporointerface     ->Norm2(&normofiterincporointerface_);

  iterincfluid             ->Norm2(&normofiterincfluid_);
  iterincfluidvelocity     ->Norm2(&normofiterincfluidvelocity_);
  iterincfluidpressure     ->Norm2(&normofiterincfluidpressure_);
  iterincfluidinterface    ->Norm2(&normofiterincfluidinterface_);

  iterincale               ->Norm2(&normofiterincale_);


  return;
}//BuildConvergenceNorms

void FPSI::Monolithic::PrintNewtonIter()
{
  // print to standard out
  // replace myrank_ here general by Comm().MyPID()
  if ((Comm().MyPID() == 0) and printscreen_ and (Step()%printscreen_==0) and printiter_)
  {
    if (iter_ == 1)
      PrintNewtonIterHeader(stdout);
      PrintNewtonIterText(stdout);
  }

  // print to error file
  if (printerrfile_ and printiter_)
  {
    if (iter_ == 1)
      PrintNewtonIterHeader(errfile_);
    PrintNewtonIterText(errfile_);
  }

  return;
} // PrintNewtonIter()

void FPSI::Monolithic::PrintNewtonIterHeader(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(8) << "numiter";

  // different style due relative or absolute error checking
  // displacement
  switch (normtypefres_)
  {
    case INPAR::FPSI::absoluteconvergencenorm:
      oss << std::setw(11) << "abs-res";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::FPSI::absoluteconvergencenorm:
    oss <<std::setw(11)<< "abs-inc";
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  switch ( normtypefres_ )
  {
    case INPAR::FPSI::absoluteconvergencenorm :
    oss <<std::setw(16)<< "poro-s-res";
   // oss <<std::setw(15)<< "abs-f-res";
    oss <<std::setw(15)<< "poro-fvel-res";
    oss <<std::setw(15)<< "poro-fpres-res";
    oss <<std::setw(15)<< "fld-fvel-res";
    oss <<std::setw(15)<< "fld-fpres-res";
    //oss <<std::setw(15)<< "pinterface-res";
    //oss <<std::setw(15)<< "finterface-res";
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::FPSI::absoluteconvergencenorm :
    oss <<std::setw(15)<< "poro-s-inc";
   // oss <<std::setw(15)<< "abs-f-inc";
    oss <<std::setw(16)<< "poro-fvel-inc";
    oss <<std::setw(16)<< "poro-fpres-inc";
    oss <<std::setw(15)<< "fld-fvel-inc";
    oss <<std::setw(15)<< "fld-fpres-inc";
    oss <<std::setw(11)<< "ale-inc";
    oss <<std::setw(14)<< "poro-int-inc";
    oss <<std::setw(14)<< "fld-int-inc";
    break;
    default:
    dserror("Begin at the beginning and go on till you come to the end. Then stop.");
    break;
  }

//  // add solution time
//  oss << std::setw(14)<< "wct";

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile==NULL)
  dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;
}  // PrintNewtonIterHeader()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen                       |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void FPSI::Monolithic::PrintNewtonIterText(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter  state etc
  oss << std::setw(4) << iter_;

  // different style due relative or absolute error checking
  // displacement
  switch (normtypefres_)
  {
    case INPAR::FPSI::absoluteconvergencenorm:
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofrhs_;
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::FPSI::absoluteconvergencenorm :
    oss << std::setw(12) << std::setprecision(5) << std::scientific << normofiterinc_;
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  switch ( normtypefres_ )
  {
    case INPAR::FPSI::absoluteconvergencenorm :
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporostruct_;
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporofluidvelocity_;
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporofluidpressure_;
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsfluidvelocity_;
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsfluidpressure_;
    //oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporointerface_;
    //oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsfluidinterface_;
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  switch ( normtypeinc_ )
  {
    case INPAR::FPSI::absoluteconvergencenorm:
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normofiterincporostruct_;
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normofiterincporofluidvelocity_;
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normofiterincporofluidpressure_;
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normofiterincfluidvelocity_;
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normofiterincfluidpressure_;
    oss << std::setw(15) << std::setprecision(5) << std::scientific << normofiterincale_;
    oss << std::setw(13) << std::setprecision(5) << std::scientific << normofiterincfluidinterface_;
    oss << std::setw(14) << std::setprecision(5) << std::scientific << normofiterincporointerface_;
    break;
    default:
    dserror("You should not turn up here.");
    break;
  }

  // add solution time
  //oss << std::setw(14) << std::setprecision(2) << std::scientific << timer_.ElapsedTime();

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile==NULL)
  dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;

}  // PrintN

//
void FPSI::Monolithic::SetupNewton()
{
  // initialise equilibrium loop and norms
    iter_ = 1;
    normofrhs_                      = 0.0;
    normofiterinc_                  = 0.0;
    normrhsfluid_                   = 0.0;
    normofiterincfluid_             = 0.0;
    normrhsfluidvelocity_           = 0.0;
    normofiterincfluidvelocity_     = 0.0;
    normrhsporostruct_              = 0.0;
    normofiterincporostruct_        = 0.0;
    normofiterincporofluid_         = 0.0;
    normrhsfluidpressure_           = 0.0;
    normofiterincfluidpressure_     = 0.0;
    normrhsporofluidvelocity_       = 0.0;
    normofiterincporofluidvelocity_ = 0.0;
    normrhsporofluidpressure_       = 0.0;
    normofiterincporofluidpressure_ = 0.0;
    normrhsporointerface_           = 0.0;
    normofiterincporointerface_     = 0.0;
    normrhsfluidinterface_          = 0.0;
    normofiterincfluidinterface_    = 0.0;


    // incremental solution vector with length of all dofs
    iterinc_ =  LINALG::CreateVector(*DofRowMap(), true);
    iterinc_ -> PutScalar(0.0);

    // a zero vector of full length
    zeros_ = LINALG::CreateVector(*DofRowMap(), true);
    zeros_->PutScalar(0.0);

    PoroField() -> SetupNewton();
}


void FPSI::Monolithic::FPSIFDCheck()
{
  PoroField()->FluidField()->Discretization()->ClearState();

  PoroField()->FluidField()->Discretization()->SetState(0,"dispnp" ,PoroField()->FluidField()->Dispnp());
  PoroField()->FluidField()->Discretization()->SetState(0,"gridv"  ,PoroField()->FluidField()->GridVel());
  PoroField()->FluidField()->Discretization()->SetState(0,"dispn"  ,PoroField()->FluidField()->Dispn());
  PoroField()->FluidField()->Discretization()->SetState(0,"veln"   ,PoroField()->FluidField()->Veln());
  PoroField()->FluidField()->Discretization()->SetState(0,"velaf"  ,PoroField()->FluidField()->Velnp());
  PoroField()->FluidField()->Discretization()->SetState(0,"velnp"  ,PoroField()->FluidField()->Velnp());

  FluidField()->Discretization()->ClearState();

  FluidField()->Discretization()->SetState(0,"dispnp",FluidField()->Dispnp());
  FluidField()->Discretization()->SetState(0,"gridv" ,FluidField()->GridVel());
  FluidField()->Discretization()->SetState(0,"dispn" ,FluidField()->Dispn());
  FluidField()->Discretization()->SetState(0,"veln"  ,FluidField()->Veln());
  FluidField()->Discretization()->SetState(0,"velaf" ,FluidField()->Velnp());
  FluidField()->Discretization()->SetState(0,"velnp" ,FluidField()->Velnp());


  const double delta = 1e-7;
  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fpsidynparams = problem->FPSIDynamicParams();
  int columntocheck = fpsidynparams.get<int>("FDCheck_column");
  int rowtocheck    = fpsidynparams.get<int>("FDCheck_row");

  // matrices and vectors
  Teuchos::RCP<Epetra_Vector>  iterinc = LINALG::CreateVector(*DofRowMap(), true);
  const int dofs = iterinc->GlobalLength();
  iterinc->PutScalar(0.0);
  iterinc->ReplaceGlobalValue(0, 0, delta);

  Teuchos::RCP<Epetra_CrsMatrix>  stiff_approx = LINALG::CreateMatrix(*DofRowMap(), 81);

  Teuchos::RCP<Epetra_Vector> rhs_old = Teuchos::rcp(new Epetra_Vector(*DofRowMap(),true));
  rhs_old->Update(1.0, *rhs_, 0.0);

  Teuchos::RCP<Epetra_Vector> rhs_copy = Teuchos::rcp(new Epetra_Vector(*DofRowMap(),true));

  Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

  Teuchos::RCP<LINALG::SparseMatrix> sparse_copy = Teuchos::rcp(new LINALG::SparseMatrix(*(sparse->EpetraMatrix())));


  std::cout << "\n****************** FPSI finite difference check ******************" << std::endl;

  int dof_poro_struct = (PoroField()->StructureField()->Discretization()->NumGlobalNodes()) * 3;
  int dof_poro_fluid  = (PoroField()->FluidField()-> Discretization()->NumGlobalNodes()) * 4;
  int dof_fluid       = (FluidField()->                 Discretization()->NumGlobalNodes()) * 4;
  int dof_ale         = (AleField()->                   Discretization()->NumGlobalNodes()) * 3;

  std::cout << "poro structure field has " << dof_poro_struct << " DOFs" << std::endl;
  std::cout << "poro fluid field has     " << dof_poro_fluid  << " DOFs" << std::endl;
  std::cout << "fluid field has          " << dof_fluid       << " DOFs" << std::endl;
  std::cout << "ale field has            " << dof_ale         << " DOFs" << std::endl;
  std::cout << "in total                 " << dofs            << " DOFs" << std::endl;


  for (int i = 0; i < dofs; ++i) // loop over columns
  {
    if (CombinedDBCMap()->MyGID(i))
    {
      iterinc->ReplaceGlobalValue(i, 0, 0.0);
    }

    Evaluate(iterinc); // initial iterinc is varied at first dof (0-th element)
    SetupSystemMatrix();

    rhs_copy->Update(1.0, *rhs_, 0.0);
    iterinc_->PutScalar(0.0); // Useful? depends on solver and more
    PoroField()->ClearPoroIterinc();
    LINALG::ApplyDirichlettoSystem(sparse_copy, iterinc_, rhs_copy,
        Teuchos::null, zeros_, *CombinedDBCMap());

    rhs_copy->Update(-1.0, *rhs_old, 1.0); // finite difference approximation of partial derivative
    rhs_copy->Scale(-1.0 / delta);
    if(i == columntocheck)
    {
      std::cout<<"iterinc:  "<<*iterinc  <<std::endl;
      std::cout<<"rhs_old:  "<<*rhs_old  <<std::endl;
      std::cout<<"rhs_copy: "<<*rhs_copy <<std::endl;
      dserror("");
    }

    int* index = &i;
    for (int j = 0; j < dofs; ++j) // loop over rows
    {
      double value = (*rhs_copy)[j];
      stiff_approx->InsertGlobalValues(j, 1, &value, index); //int InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);

    } //j-loop (rows)

    if (not CombinedDBCMap()->MyGID(i))
      iterinc->ReplaceGlobalValue(i, 0, -delta);

    iterinc->ReplaceGlobalValue(i - 1, 0, 0.0);

    if (i != dofs - 1)
      iterinc->ReplaceGlobalValue(i + 1, 0, delta);

  } // i-loop (columns)

  Evaluate(iterinc);
  SetupSystemMatrix();


  stiff_approx->FillComplete();

  Teuchos::RCP<LINALG::SparseMatrix> temp = Teuchos::rcp(new LINALG::SparseMatrix(*stiff_approx));

  Teuchos::RCP<Epetra_CrsMatrix> stiff_approx_sparse = temp->EpetraMatrix();

  Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = sparse_copy->EpetraMatrix();

  // calc error (subtraction of sparse_crs and stiff_approx_sparse)
  for(int i = 0; i < dofs; i++)
  {
    int length;
    int numentries = sparse_crs->NumGlobalEntries(i);
    std::vector<double> values (numentries);
    std::vector<int>    indices(numentries);
    sparse_crs->ExtractGlobalRowCopy(i,numentries,length,&values[0],&indices[0]);

    for(int k = 0; k<numentries; k++)
    {
      values[k] = -values[k];
    }

    stiff_approx_sparse ->SumIntoGlobalValues(i,numentries,&values[0],&indices[0]);
  }


  stiff_approx_sparse ->FillComplete();
  sparse_crs ->FillComplete();



  //std::cout<<*CombinedDBCMap()<<std::endl;
  bool success = true;
  double error_max = 0.0;
  for (int i = 0; i < dofs; ++i)
  {
    if (not CombinedDBCMap()->MyGID(i)) // only if there is no dirichlet value on dof
    {
      for (int j = 0; j < dofs; ++j)
      {
        if (not CombinedDBCMap()->MyGID(j))
        {
          double stiff_approx_ij;
          double sparse_ij;
          double error_ij;

          // get error_crs entry ij
          int errornumentries;
          int errorlength = stiff_approx_sparse->NumGlobalEntries(i);
          std::vector<double> errorvalues(errorlength);
          std::vector<int> errorindices(errorlength);
          stiff_approx_sparse->ExtractGlobalRowCopy(i,errorlength,errornumentries,&errorvalues[0],&errorindices[0]);
          for(int k = 0; k < errorlength; ++k)
          {
            if(errorindices[k] == j)
            {
              error_ij=errorvalues[k];
              break;
            }
            else
              error_ij = 0.0;
          }

          // get sparse_ij entry ij
          int sparsenumentries;
          int sparselength = sparse_crs->NumGlobalEntries(i);
          std::vector<double> sparsevalues(sparselength);
          std::vector<int> sparseindices(sparselength);
          sparse_crs->ExtractGlobalRowCopy(i,sparselength,sparsenumentries,&sparsevalues[0],&sparseindices[0]);
          for(int k = 0; k < sparselength; ++k)
          {
            if(sparseindices[k] == j)
            {
              sparse_ij=sparsevalues[k];
              break;
            }
            else
              sparse_ij=0.0;
          }

          // get stiff_approx entry ij
          int approxnumentries;
          int approxlength = stiff_approx->NumGlobalEntries(i);
          std::vector<double> approxvalues(approxlength);
          std::vector<int> approxindices(approxlength);
          stiff_approx->ExtractGlobalRowCopy(i,approxlength,approxnumentries,&approxvalues[0],&approxindices[0]);
          for(int k = 0; k < approxlength; ++k)
          {
            if(approxindices[k] == j)
            {
              stiff_approx_ij=approxvalues[k];
              break;
            }
            else
              stiff_approx_ij=0.0;
          }

          double error = 0.0;
          if (abs(stiff_approx_ij) > 1e-6)
            error = error_ij / stiff_approx_ij;
          else if (abs(sparse_ij) > 1e-6)
            error = error_ij / sparse_ij;

          if(i == rowtocheck and j == columntocheck)
          {
            std::cout<<"K_approx value: "<<stiff_approx_ij<<std::endl;
            std::cout<<"K value : "<<sparse_ij<<std::endl;
            std::cout<<"error : "<<error<<std::endl;
            std::cout<<"error_ij : "<<error_ij<<std::endl;
            std::cout<<"i : "<<i<<std::endl;
            std::cout<<"j : "<<j<<std::endl;
          }

          if (abs(error) > abs(error_max))
            error_max = abs(error);

          if ((abs(error) > 1e-6))
          {
            if (abs(error_ij) > 1e-6)
            {
              std::cout << "finite difference check failed entry (" << i << "," << j
                  << ")! stiff: " << sparse_ij << ", approx: "
                  << stiff_approx_ij << " ,abs. error: " << error_ij
                  << " , rel. error: " << error << std::endl;

              success = false;
            }
          }
        }
      }
    }
  } // loop over dofs of succes check

  if(success)
  {
    std::cout << "finite difference check successful, max. rel. error: "
        << error_max << std::endl;
    std::cout << "****************** finite difference check done ***************\n\n"
        << std::endl;
  }
  else
    dserror("FPSIFDCheck failed");

  PoroField()->FluidField()->Discretization()->ClearState();
  FluidField()->Discretization()->ClearState();

  return;
}

