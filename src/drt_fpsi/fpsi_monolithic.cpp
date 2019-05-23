/*----------------------------------------------------------------------*/
/*!

\brief General framework for monolithic fpsi solution schemes

\level 3

\maintainer  Christoph Ager
*/

/*----------------------------------------------------------------------*/
// GENERAL includes
#include <Teuchos_TimeMonitor.hpp>

// FPSI includes
#include "fpsi_defines.H"
#include "fpsi_monolithic.H"
#include "fpsi_utils.H"

// POROELAST includes
#include "../drt_poroelast/poroelast_monolithic.H"

// LINALG includes
#include "../linalg/linalg_solver.H"

// drt_lib includes
#include "../drt_lib/drt_globalproblem.H"

// drt_adapter includes
#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"

// STRUCTURE includes
#include "../drt_structure/stru_aux.H"

// OTHER includes
#include "../drt_io/io_control.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

FPSI::MonolithicBase::MonolithicBase(const Epetra_Comm& comm,
    const Teuchos::ParameterList& fpsidynparams, const Teuchos::ParameterList& poroelastdynparams)
    : FPSI_Base(comm, fpsidynparams)
{
  // Creation of the subproblems
  // As for general FPSI problems overlapping FSI/FPSI - interfaces can occur, the maps in the
  // different MapExtractors of the Fields are built overlapping! This is necessary to get the right
  // coupling objects, and be able to extract and insert whole FPSI vectors! For building the block
  // matrices of the different fields the following procedure happens to get non overlapping blocks
  // (as wanted):
  // 1. the maps of the FSI and FPSI block are overlapping
  // 2. during the assembling procedure all dofs which belong to the overlapping interface will be
  // assembled into the FSI block (is before the FPSI
  //    condition in the MapExtractors)
  // 3. as there are no entries in of the overlapping interface in the FPSI block, FillComplete()
  // removes them and the FPSI block has no overlapping
  //    entries!

  // --> It's clear that this is not really easily comprehensible, but as the alternatives to this
  // approach requires quite some modifications in BACI
  //     (e.g. extra interface for FSI and FPSI for every Field ...), therefore this version should
  //     be used at the moment.


  // create instance of poroelast subproblem
  poroelast_subproblem_ = Teuchos::rcp(new POROELAST::Monolithic(comm, fpsidynparams));
  // ask base algorithm for the fluid time integrator
  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fluiddynparams = problem->FluidDynamicParams();
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fpsidynparams, fluiddynparams, "fluid", true));
  fluid_subproblem_ = Teuchos::rcp_dynamic_cast<ADAPTER::FluidFPSI>(fluid->FluidField());
  // ask base algorithm for the ale time integrator
  Teuchos::RCP<ADAPTER::AleBaseAlgorithm> ale = Teuchos::rcp(
      new ADAPTER::AleBaseAlgorithm(fpsidynparams, DRT::Problem::Instance()->GetDis("ale")));
  ale_ = Teuchos::rcp_dynamic_cast<ADAPTER::AleFpsiWrapper>(ale->AleField());
  if (ale_ == Teuchos::null) dserror("cast from ADAPTER::Ale to ADAPTER::AleFpsiWrapper failed");

  coupfa_ = Teuchos::rcp(new ADAPTER::Coupling());

  coupsf_fsi_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupsa_fsi_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupfa_fsi_ = Teuchos::rcp(new ADAPTER::Coupling());
  icoupfa_fsi_ = Teuchos::rcp(new ADAPTER::Coupling());

  Teuchos::RCP<FPSI::Utils> FPSI_UTILS = FPSI::Utils::Instance();

  Fluid_PoroFluid_InterfaceMap = FPSI_UTILS->Get_Fluid_PoroFluid_InterfaceMap();
  PoroFluid_Fluid_InterfaceMap = FPSI_UTILS->Get_PoroFluid_Fluid_InterfaceMap();

  // build a proxy of the fluid discretization for the ale field
  if (FluidField()->Discretization()->AddDofSet(
          AleField()->WriteAccessDiscretization()->GetDofSetProxy()) != 1)
  {
    dserror("unexpected numbers of dofsets in fluid field");
  }
  FluidField()->Discretization()->FillComplete(true, false, false);
}  // MonolithicBase

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FPSI::MonolithicBase::~MonolithicBase() {}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::ReadRestart(int step)
{
  PoroField()->ReadRestart(step);
  FluidField()->ReadRestart(step);
  AleField()->ReadRestart(step);

  SetTimeStep(FluidField()->Time(), FluidField()->Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::PrepareTimeStep()
{
  IncrementTimeAndStep();
  PrintHeader();

  PoroField()->PrepareTimeStep();
  PoroField()->SetupNewton();
  FluidField()->PrepareTimeStep();
  AleField()->PrepareTimeStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::Update()
{
  PoroField()->Update();
  FluidField()->Update();
  AleField()->Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::PrepareOutput() { PoroField()->PrepareOutput(); }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::Output()
{
  PoroField()->Output();
  FluidField()->Output();
  AleField()->Output();
}

/*----------------------------------------------------------------------*/
/*                          Coupling Methods                            */
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::FluidToAle(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->MasterToSlave(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::AleToFluid(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}
/// Just in use for problems with FSI-interface ///
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::StructToFluid_FSI(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_fsi_->MasterToSlave(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::FluidToStruct_FSI(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_fsi_->SlaveToMaster(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::StructToAle_FSI(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_fsi_->MasterToSlave(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::AleToStruct_FSI(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_fsi_->SlaveToMaster(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::FluidToAle_FSI(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_fsi_->MasterToSlave(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::AleToFluid_FSI(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_fsi_->SlaveToMaster(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::AleToFluidInterface_FSI(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_fsi_->SlaveToMaster(iv);
}
/// ---------------------------------------------- ///

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<  MonolithicBase -> Monolithic  >>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

FPSI::Monolithic::Monolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& fpsidynparams,
    const Teuchos::ParameterList& poroelastdynparams)
    : MonolithicBase(comm, fpsidynparams, poroelastdynparams),
      directsolve_(true),
      printscreen_(true),
      printiter_(true),
      printerrfile_(true),
      errfile_(DRT::Problem::Instance()->ErrorFile()->Handle()),
      timer_(comm),
      isfirsttimestep_(true),
      islinesearch_(false),
      firstcall_(true)
{
  const Teuchos::ParameterList& sdynparams = DRT::Problem::Instance()->StructuralDynamicParams();
  solveradapttol_ = (DRT::INPUT::IntegralValue<int>(sdynparams, "ADAPTCONV") == 1);
  solveradaptolbetter_ = (sdynparams.get<double>("ADAPTCONV_BETTER"));

  // hydraulic conductivity (needed for coupling in case of probtype fps3i)
  // is overwritten in class fs3i
  conductivity_ = 0.0;

  // Check if FSI-Interface exists and set flag
  // Will be used to jump over all sections, which are just for FSI condensation procedure required!

  if (FluidField()->Interface()->Map(FLD::UTILS::MapExtractor::cond_fsi)->NumGlobalElements())
  {
    FSI_Interface_exists_ = true;
    if (comm.MyPID() == 0)
      std::cout << "FPSI Calculation will be performed with FSI - Interface!" << std::endl;
  }
  else
  {
    FSI_Interface_exists_ = false;
    if (comm.MyPID() == 0)
      std::cout << "FPSI Calculation will skip all FSI parts as there is no FSI - Interface!"
                << std::endl;
  }

  // Check for valid predictors
  if (sdynparams.get<std::string>("PREDICT") != "ConstDis")
    dserror(
        "No Structural Predictor for FPSI implemented at the moment, choose <PREDICT = ConstDis> "
        "in you .dat file! \n --> Or feel free to add the missing terms coming from the predictors "
        "to BACI!");

  const Teuchos::ParameterList& fdynparams = DRT::Problem::Instance()->FluidDynamicParams();
  if (fdynparams.get<std::string>("PREDICTOR") != "steady_state")
    dserror(
        "No Fluid Predictor for FPSI implemented at the moment, choose <PREDICTOR = steady_state> "
        "in you .dat file! \n --> Or feel free to add the missing terms coming from the predictors "
        "to BACI!");

  active_FD_check_ = false;  // to avoid adding RHS of firstiter moreoften!
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::SetupSystem()
{
  const int ndim = DRT::Problem::Instance()->NDim();

  ADAPTER::Coupling& coupfa = FluidAleCoupling();

  const Epetra_Map* fluidnodemap = FluidField()->Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap = AleField()->Discretization()->NodeRowMap();

  coupfa.SetupCoupling(*FluidField()->Discretization(), *AleField()->Discretization(),
      *fluidnodemap, *alenodemap, ndim, false);
  FluidField()->SetMeshMap(coupfa.MasterDofMap());

  if (FSI_Interface_exists_) SetupSystem_FSI();

  // Setup the FPSI Coupling Adapter
  FPSICoupl() = Teuchos::rcp(new FPSI::FPSICoupling(poroelast_subproblem_, fluid_subproblem_, ale_,
      Fluid_PoroFluid_InterfaceMap, PoroFluid_Fluid_InterfaceMap));
  FPSICoupl()->SetConductivity(conductivity_);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::SetupSystem_FSI()
{
  // create local FPSI MapExtractors as the MapExtractors of the single Fields are nonoverlapping
  // and therefore not the whole FPSI Interface Map is availlable for FSI, FPSI Interface overlap!!!

  // right now we use matching meshes at the interface

  const int ndim = DRT::Problem::Instance()->NDim();

  ADAPTER::Coupling& coupsf_fsi = StructureFluidCoupling_FSI();
  ADAPTER::Coupling& coupsa_fsi = StructureAleCoupling_FSI();
  ADAPTER::Coupling& icoupfa_fsi = InterfaceFluidAleCoupling_FSI();

  // structure to fluid
  coupsf_fsi.SetupConditionCoupling(*PoroField()->StructureField()->Discretization(),
      PoroField()->StructureField()->Interface()->FSICondMap(), *FluidField()->Discretization(),
      FluidField()->Interface()->FSICondMap(), "FSICoupling", ndim);
  // structure to ale
  coupsa_fsi.SetupConditionCoupling(*PoroField()->StructureField()->Discretization(),
      PoroField()->StructureField()->Interface()->FSICondMap(), *AleField()->Discretization(),
      AleField()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // fluid to ale at the interface

  icoupfa_fsi.SetupConditionCoupling(*FluidField()->Discretization(),
      FluidField()->Interface()->FSICondMap(), *AleField()->Discretization(),
      AleField()->Interface()->FSICondMap(), "FSICoupling", ndim);


  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf_fsi.MasterDofMap()->SameAs(*coupsa_fsi.MasterDofMap()))
    dserror("fsi structure interface dof maps do not match");

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::Timeloop()
{
  PrepareTimeloop();

  while (NotFinished())  // while step < maxsteps and time < maxtime
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
  while ((((not Converged()) and (iter_ <= maximumiterations_)) or (iter_ <= minimumiterations_)) or
         islinesearch_ == true)
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
    // build norms
    BuildConvergenceNorms();

    // print stuff
    if (islinesearch_ == false) PrintNewtonIter();

    // reset solver tolerance
    solver_->ResetTolerance();

    // increment equilibrium loop index
    if (islinesearch_ == false)
    {
      iter_ += 1;
      PoroField()->IncrementPoroIter();
    }
  }  // end loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ((Converged()) and (Comm().MyPID() == 0))
  {
    if (linesearch_counter > 0.5)
      std::cout << "            Evaluation of residual with scaled increment yields: " << normofrhs_
                << std::endl;
    islinesearch_ = false;
    linesearch_counter = 0.;
  }
  else if (iter_ >= maximumiterations_)
  {
    dserror("Newton found no convergence in %d iterations", iter_);
  }

  PoroField()->RecoverLagrangeMultiplierAfterTimeStep();

  // recover Lagrange multiplier \lambda_{\Gamma} at the interface at the end of each time step
  // (i.e. condensed traction/forces onto the structure) needed for rhs in next time step
  if (FSI_Interface_exists_)
    RecoverLagrangeMultiplier();  // LagrangeMultiplier of the FSI interface!
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FPSI::Monolithic::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::Evaluate");

  if (linesearch_ and islinesearch_ == false)
  {
    FluidField()->Discretization()->ClearState();
    linesearch_counter = 0.;
    FluidField()->Discretization()->SetState(0, "dispnp", FluidField()->Dispnp());
    meshdispold_ = AleToFluid(AleField()->Dispnp());
    porointerfacedisplacementsold_ =
        FPSICoupl()->iPorostructToAle(PoroField()->StructureField()->ExtractInterfaceDispnp(true));
  }


  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> pfx;
  Teuchos::RCP<const Epetra_Vector> ax;

  if (x != Teuchos::null)
  {
    ExtractFieldVectors(x, sx, pfx, fx, ax, (iter_ == 1 and !active_FD_check_));
  }
  else
  {
    dserror("No existing increment vector !");
  }

  PoroField()->Evaluate(sx, pfx);

  Teuchos::RCP<Epetra_Vector> porointerfacedisplacements_FPSI =
      FPSICoupl()->iPorostructToAle(PoroField()->StructureField()->ExtractInterfaceDispnp(true));
  AleField()->ApplyInterfaceDisplacements(porointerfacedisplacements_FPSI);

  if (FSI_Interface_exists_)
  {
    Teuchos::RCP<Epetra_Vector> porointerfacedisplacements_FSI =
        StructToAle_FSI(PoroField()->StructureField()->ExtractInterfaceDispnp(false));
    AleField()->ApplyFSIInterfaceDisplacements(porointerfacedisplacements_FSI);
  }

  AleField()->WriteAccessDispnp()->Update(
      1.0, *ax, 1.0);  // displacement increments on the interfaces are zero!!!
  AleField()->Evaluate(Teuchos::null);

  Teuchos::RCP<const Epetra_Vector> aledisplacements = AleToFluid(AleField()->Dispnp());
  FluidField()->ApplyMeshDisplacement(aledisplacements);

  FluidField()->UpdateNewton(fx);

  FluidField()->Evaluate(Teuchos::null);

  // Evaluate FPSI Coupling Matrixes and RHS
  FPSICoupl()->EvaluateCouplingMatrixesRHS();

  SetupRHS(iter_ == 1);

}  // Evaluate

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

void FPSI::Monolithic::SetupSolver()
{
  const Teuchos::ParameterList& fpsidynamicparams = DRT::Problem::Instance()->FPSIDynamicParams();

  const int linsolvernumber = fpsidynamicparams.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror(
        "No linear solver defined for FPSI problem. Please set LINEAR_SOLVER in FPSI DYNAMIC to a "
        "valid number !");

  const Teuchos::ParameterList& solverparams =
      DRT::Problem::Instance()->SolverParams(linsolvernumber);
  const int solvertype =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams, "SOLVER");

  directsolve_ = (solvertype == INPAR::SOLVER::umfpack or solvertype == INPAR::SOLVER::superlu or
                  solvertype == INPAR::SOLVER::amesos_klu_nonsym);

  if (directsolve_)
    solver_ = Teuchos::rcp(
        new LINALG::Solver(solverparams, Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
  else
    // create a linear solver
    CreateLinearSolver();

  // Get the parameters for the Newton iteration
  maximumiterations_ = fpsidynamicparams.get<int>("ITEMAX");
  minimumiterations_ = fpsidynamicparams.get<int>("ITEMIN");
  normtypeinc_ =
      DRT::INPUT::IntegralValue<INPAR::FPSI::ConvergenceNorm>(fpsidynamicparams, "NORM_INC");
  normtypefres_ =
      DRT::INPUT::IntegralValue<INPAR::FPSI::ConvergenceNorm>(fpsidynamicparams, "NORM_RESF");
  combinedconvergence_ =
      DRT::INPUT::IntegralValue<INPAR::FPSI::BinaryOp>(fpsidynamicparams, "NORMCOMBI_RESFINC");

  // toleranceiterinc_        = fpsidynamicparams.get<double> ("INCTOL");
  // toleranceresidualforces_ = fpsidynamicparams.get<double> ("RESTOL");

  {
    std::istringstream tolresstream(
        Teuchos::getNumericStringParameter(fpsidynamicparams, "RESTOL"));
    double word;
    while (tolresstream >> word)
    {
      toleranceresidualforceslist_.push_back(word);
    }
    toleranceresidualforces_ = toleranceresidualforceslist_[0];

    std::istringstream tolincstream(
        Teuchos::getNumericStringParameter(fpsidynamicparams, "INCTOL"));
    while (tolincstream >> word)
    {
      toleranceiterinclist_.push_back(word);
    }
    toleranceiterinc_ = toleranceiterinclist_[0];
  }

  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fpsidynparams = problem->FPSIDynamicParams();
  linesearch_ = DRT::INPUT::IntegralValue<int>(fpsidynparams, "LineSearch");
  if (linesearch_ == 1)
    dserror(
        "Parameter 'LineSearch' is set to 'Yes' in the FPSI Dynamic section in your input-file.  \n"
        "Though the framework for a line search algorithm is implemented in fpsi_monolithic.cpp, \n"
        "a proper routine to reset the participating single fields is still required. In Chuck's \n"
        "experimental baci this was solved by performing an evaluate with the negative increment.\n"
        "However this has not yet been committed.\n");
  linesearch_counter = 0.;

  return;
}

/*----------------------------------------------------------------------*
 | create linear solver                                   vuong 08/15 |
 *----------------------------------------------------------------------*/
void FPSI::Monolithic::CreateLinearSolver()
{
  // get dynamic section
  const Teuchos::ParameterList& fpsidyn = DRT::Problem::Instance()->FPSIDynamicParams();

  // get the linear solver number
  const int linsolvernumber = fpsidyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror(
        "No linear solver defined for FPSI problem. Please set LINEAR_SOLVER in FPSI DYNAMIC to a "
        "valid number !");

  // get parameter list of structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // use solver blocks for structure
  // get the solver number used for structural solver
  const int slinsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (slinsolvernumber == (-1))
    dserror(
        "no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL "
        "DYNAMIC to a valid number!");

  // get parameter list of fluid dynamics
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  // use solver blocks for fluid
  // get the solver number used for fluid solver
  const int flinsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
  // check if the fluid solver has a valid solver number
  if (flinsolvernumber == (-1))
    dserror(
        "no linear solver defined for fluid field. Please set LINEAR_SOLVER in FLUID DYNAMIC to a "
        "valid number!");

  // get parameter list of structural dynamics
  const Teuchos::ParameterList& aledyn = DRT::Problem::Instance()->AleDynamicParams();
  // use solver blocks for structure
  // get the solver number used for structural solver
  const int alinsolvernumber = aledyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (alinsolvernumber == (-1))
    dserror(
        "no linear solver defined for ALE field. Please set LINEAR_SOLVER in ALE DYNAMIC to a "
        "valid number!");

  // get solver parameter list of linear Poroelasticity solver
  const Teuchos::ParameterList& fpsisolverparams =
      DRT::Problem::Instance()->SolverParams(linsolvernumber);

  const int solvertype =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(fpsisolverparams, "SOLVER");

  if (solvertype != INPAR::SOLVER::aztec_msr && solvertype != INPAR::SOLVER::belos)
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << " Note: the BGS2x2 preconditioner now " << std::endl;
    std::cout << " uses the structural solver and fluid solver blocks" << std::endl;
    std::cout << " for building the internal inverses" << std::endl;
    std::cout << " Remove the old BGS PRECONDITIONER BLOCK entries " << std::endl;
    std::cout << " in the dat files!" << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    dserror("aztec solver expected");
  }
  const int azprectype =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(fpsisolverparams, "AZPREC");

  // plausibility check
  switch (azprectype)
  {
    case INPAR::SOLVER::azprec_AMGnxn:
    {
      // no plausibility checks here
      // if you forget to declare an xml file you will get an error message anyway
    }
    break;
    default:
      dserror("AMGnxn preconditioner expected");
      break;
  }

  solver_ = Teuchos::rcp(new LINALG::Solver(
      fpsisolverparams, Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

  // use solver blocks for structure and fluid
  const Teuchos::ParameterList& ssolverparams =
      DRT::Problem::Instance()->SolverParams(slinsolvernumber);
  const Teuchos::ParameterList& fsolverparams =
      DRT::Problem::Instance()->SolverParams(flinsolvernumber);
  const Teuchos::ParameterList& asolverparams =
      DRT::Problem::Instance()->SolverParams(alinsolvernumber);

  // for now, use same solver parameters for poro fluid and free fluid

  // poro/structure
  solver_->PutSolverParamsToSubParams("Inverse1", ssolverparams);
  // poro fluid
  solver_->PutSolverParamsToSubParams("Inverse2", fsolverparams);
  // fluid
  solver_->PutSolverParamsToSubParams("Inverse3", fsolverparams);
  // ale
  solver_->PutSolverParamsToSubParams("Inverse4", asolverparams);

  // prescribe rigid body modes
  PoroField()->StructureField()->Discretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse1"));
  PoroField()->FluidField()->Discretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse2"));
  FluidField()->Discretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse3"));
  AleField()->WriteAccessDiscretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse4"));

  // fixing length of Inverse1 nullspace (solver/preconditioner ML)
  {
    std::string inv = "Inverse1";
    const Epetra_Map& oldmap = *(PoroField()->StructureField()->DofRowMap());
    const Epetra_Map& newmap =
        systemmatrix_->Matrix(structure_block_, structure_block_).EpetraMatrix()->RowMap();
    solver_->FixMLNullspace(&inv[0], oldmap, newmap, solver_->Params().sublist("Inverse1"));
  }
  // fixing length of Inverse2 nullspace (solver/preconditioner ML)
  {
    std::string inv = "Inverse2";
    const Epetra_Map& oldmap = *(PoroField()->FluidField()->DofRowMap());
    ;
    const Epetra_Map& newmap =
        systemmatrix_->Matrix(porofluid_block_, porofluid_block_).EpetraMatrix()->RowMap();
    solver_->FixMLNullspace(&inv[0], oldmap, newmap, solver_->Params().sublist("Inverse2"));
  }
  // fixing length of Inverse3 nullspace (solver/preconditioner ML)
  {
    std::string inv = "Inverse3";
    const Epetra_Map& oldmap = *(FluidField()->DofRowMap());
    const Epetra_Map& newmap =
        systemmatrix_->Matrix(fluid_block_, fluid_block_).EpetraMatrix()->RowMap();
    solver_->FixMLNullspace(&inv[0], oldmap, newmap, solver_->Params().sublist("Inverse3"));
  }
  // fixing length of Inverse4 nullspace (solver/preconditioner ML)
  {
    std::string inv = "Inverse4";
    const Epetra_Map& oldmap = *(AleField()->DofRowMap());
    ;
    const Epetra_Map& newmap =
        systemmatrix_->Matrix(ale_i_block_, ale_i_block_).EpetraMatrix()->RowMap();
    solver_->FixMLNullspace(&inv[0], oldmap, newmap, solver_->Params().sublist("Inverse4"));
  }
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
  if (Teuchos::getIntegralValue<int>(fpsidynparams, "FDCheck"))
  {
    FPSIFDCheck();
  }

  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more
  PoroField()->ClearPoroIterinc();

  if (directsolve_)
  {
    Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

    if (FSI_Interface_exists_)
    {
      // remove entries in condensed dofs from matrix and rhs...
      LINALG::ApplyDirichlettoSystem(
          sparse, iterinc_, rhs_, Teuchos::null, zeros_, *FluidField()->Interface()->FSICondMap());
    }

    LINALG::ApplyDirichlettoSystem(
        sparse, iterinc_, rhs_, Teuchos::null, zeros_, *CombinedDBCMap());

    // line search
    if (linesearch_)
      LineSearch(sparse);
    else
    {
      // standard solver call
      solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, true, iter_ == 1);
    }
  }
  else
  {
    if (FSI_Interface_exists_)
    {
      // remove entries in condensed dofs from matrix and rhs...
      LINALG::ApplyDirichlettoSystem(systemmatrix_, iterinc_, rhs_, Teuchos::null, zeros_,
          *FluidField()->Interface()->FSICondMap());
    }

    LINALG::ApplyDirichlettoSystem(
        systemmatrix_, iterinc_, rhs_, Teuchos::null, zeros_, *CombinedDBCMap());

    // standard solver call
    solver_->Solve(systemmatrix_->EpetraOperator(), iterinc_, rhs_, true, iter_ == 1);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::LineSearch(Teuchos::RCP<LINALG::SparseMatrix>& sparse)
{
  // Note: the line search code seems to be experimental and is not
  // working properly (perhaps just a sign is wrong somewhere ...)
  // I would not recommend to use it without thorough testing
  // vuong 08/15

  if (iter_ > 1)
  {
    rhs_->Norm2(&normofrhs_);
    if (normofrhs_ - normofrhsold_ > 1e-13)
    {
      if (linesearch_counter > 0.5)
        std::cout << "            Evaluation of residual with bisected increment yields: "
                  << normofrhs_ << std::endl;

      islinesearch_ = true;
      iterinc_->Update(pow(0.5, (linesearch_counter)), *iterincold_, 0.0);
      linesearch_counter = linesearch_counter + 1.0;
      std::cout << "linesearch_ : " << std::setprecision(1)
                << static_cast<int>(linesearch_counter + 0.5) << " iterinc_ multiplied by "
                << std::setprecision(4) << pow(0.5, linesearch_counter)
                << "   residual = " << normofrhs_ << " > " << normofrhsold_ << std::endl;

      // substract the old interinc_ from all fields (undo the update)
      Teuchos::RCP<Epetra_Vector> sx;
      Teuchos::RCP<Epetra_Vector> pfx;
      Teuchos::RCP<Epetra_Vector> fx;
      Teuchos::RCP<const Epetra_Vector> constsx;
      Teuchos::RCP<const Epetra_Vector> constfpx;
      Teuchos::RCP<const Epetra_Vector> constfx;
      Teuchos::RCP<const Epetra_Vector> ax;

      sx = Teuchos::rcp(new Epetra_Vector(*PoroField()->StructureField()->DofRowMap(), true));
      pfx = Teuchos::rcp(new Epetra_Vector(*PoroField()->FluidField()->DofRowMap(), true));
      fx = Teuchos::rcp(new Epetra_Vector(*FluidField()->DofRowMap(), true));

      ExtractFieldVectors(iterinc_, constsx, constfpx, constfx, ax, iter_ == 1);
      iterinc_->Norm2(&normofiterinc_);
      std::cout << "            Norm of step back: " << normofiterinc_ << std::endl;
      // PoroField()  ->ResetNewton(sx);
      // FluidField() ->ResetNewton(fx);
      sx->Update(1.0, *constsx, 0.0);
      pfx->Update(1.0, *constfpx, 0.0);
      fx->Update(1.0, *constfx, 0.0);
      sx->Scale(-1.0);
      pfx->Scale(-1.0);
      fx->Scale(-1.0);
      PoroField()->Evaluate(sx, pfx);
      FluidField()->UpdateNewton(Teuchos::rcp_dynamic_cast<const Epetra_Vector>(fx));
      // AleField()   ->ResetNewton(ax);

      FluidField()->ApplyMeshDisplacement(meshdispold_);
      AleField()->ApplyInterfaceDisplacements(porointerfacedisplacementsold_);


      // set iterinc_ to a fraction of the old iterinc_
      iterinc_->Update(pow(0.5, linesearch_counter), *iterincold_, 0.0);
      iterinc_->Norm2(&normofiterinc_);
      std::cout << "            Norm of old increment: " << normofiterincold_
                << "  Norm of bisected increment: " << normofiterinc_ << std::endl;
    }
    else
    {
      islinesearch_ = false;
      if (linesearch_counter > 0.5)
        std::cout << "            Evaluation of residual with bisected increment yields: "
                  << normofrhs_ << std::endl;
      linesearch_counter = 0.0;
    }
  }

  // prepare linesearch_
  // copy the old iterinc_ before new solve
  if (linesearch_ and islinesearch_ == false)
  {
    rhsold_ = LINALG::CreateVector(*DofRowMap(), true);
    rhsold_->Update(1.0, *rhs_, 0.0);
    rhsold_->Norm2(&normofrhsold_);
    if (abs(normofrhs_ - normofrhsold_) > 1.e-12 and iter_ > 1) dserror(" wrong copy of rhs_ ");
  }
  // end prepare linesearch_

  // standard solver call
  if (islinesearch_ == false)
  {
    solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, true, iter_ == 1);
  }

  if (islinesearch_ == false)
  {
    // check whether iterinc_ points in right direction
    Teuchos::RCP<Epetra_Vector> tempvec = LINALG::CreateVector(*DofRowMap(), true);
    sparse->Multiply(true, *rhs_, *tempvec);
    double climb = 0.0;
    tempvec->Dot(*iterinc_, &climb);
    climb = -climb;

    if (climb > 0.0)
    {
      std::cout << "########################################################################"
                << std::endl;
      std::cout << "##                                                                    ##"
                << std::endl;
      std::cout << "## WARNING: A*x-b=0 ; A^T*b*x > 0 ; increment vector multiplied by -1 ##"
                << std::endl;
      std::cout << "##                                                                    ##"
                << std::endl;
      std::cout << "##                       Value = " << std::setprecision(9) << climb << "    ##"
                << std::endl;
      std::cout << "##                                                                    ##"
                << std::endl;
      std::cout << "########################################################################"
                << std::endl;
      iterinc_->Update(-1.0, *iterinc_, 0.0);
    }
  }

  if (linesearch_ and islinesearch_ == false)
  {
    iterincold_ = LINALG::CreateVector(*DofRowMap(), true);
    iterincold_->Update(1.0, *iterinc_, 0.0);
    iterincold_->Norm2(&normofiterincold_);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> FPSI::Monolithic::CombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map> scondmap =
      PoroField()->StructureField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> pfcondmap =
      PoroField()->FluidField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> fcondmap = FluidField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> acondmap = AleField()->GetDBCMapExtractor()->CondMap();
  Teuchos::RCP<Epetra_Map> tempmap = LINALG::MergeMap(scondmap, pfcondmap, false);
  Teuchos::RCP<Epetra_Map> condmap_0 = LINALG::MergeMap(tempmap, fcondmap, false);
  Teuchos::RCP<Epetra_Map> condmap = LINALG::MergeMap(condmap_0, acondmap, false);

  return condmap;
}  // CombinedDBCMap()

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<                 Newton Loop              >>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<                   Methods                >>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool FPSI::Monolithic::Converged()
{
  // check for single norms
  bool convinc = false;
  bool convfres = false;

  // residual increments
  switch (normtypeinc_)
  {
    case INPAR::FPSI::absoluteconvergencenorm:
      convinc = normofiterinc_ < toleranceiterinc_;
      break;
    case INPAR::FPSI::absoluteconvergencenorm_sys_split:
      dserror(
          "Check for convergence of primary variables with type <absolut sys split> not "
          "implemented yet!");
      break;
    case INPAR::FPSI::relativconvergencenorm_sys:
      // increment convergence is checked relative to the average dofs of the different fields
      std::cout << " |ps| " << norm1_ps_ << " |pfv| " << norm1_pfv_ << " |pfp| " << norm1_pfp_
                << " |fv| " << norm1_fv_ << " |fp| " << norm1_fp_ << " |a| " << norm1_a_
                << std::endl;
      convinc =
          ((normofiterincporostruct_ / norm1_ps_ * sqrtnps_ < toleranceiterinc_) and  // poro struct
              (normofiterincporofluidvelocity_ / norm1_pfv_ * sqrtnpfv_ <
                  toleranceiterinc_) and  // poro fluid velocity
              (normofiterincporofluidpressure_ / norm1_pfp_ * sqrtnpfp_ <
                  toleranceiterinc_) and  // poro fluid pressure
              (normofiterincale_ / norm1_a_ * sqrtna_ < toleranceiterinc_) and  // ale
              (normofiterincfluidvelocity_ / norm1_fv_ * sqrtnfv_ <
                  toleranceiterinc_) and  // fluid velocity
              (normofiterincfluidpressure_ / norm1_fp_ * sqrtnfp_ <
                  toleranceiterinc_)  // fluid pressure
          );
      break;
    default:
      dserror("Cannot check for convergence of primary variables for any reason :-p !");
      break;
  }
  // residual forces
  switch (normtypefres_)
  {
    case INPAR::FPSI::absoluteconvergencenorm:
      convfres = normofrhs_ < toleranceresidualforces_;
      break;
    case INPAR::FPSI::absoluteconvergencenorm_sys_split:
      convfres = ((normrhsporostruct_ / sqrtnps_ < toleranceresidualforceslist_[2]) and
                  (normrhsfluidvelocity_ / sqrtnpfv_ < toleranceresidualforceslist_[0]) and
                  (normrhsporofluidpressure_ / sqrtnpfp_ < toleranceresidualforceslist_[1]) and
                  (normrhsale_ / sqrtna_ < toleranceresidualforceslist_[5]) and
                  (normrhsfluidvelocity_ / sqrtnfv_ < toleranceresidualforceslist_[3]) and
                  (normrhsfluidpressure_ / sqrtnfp_ < toleranceresidualforceslist_[4]));
      break;
    case INPAR::FPSI::relativconvergencenorm_sys:
      dserror(
          "Check for convergence of residual forces with type <relativ_sys> not implemented yet!");
      break;
    default:
      dserror("Cannot check for convergence of residual forces for any reason :-P !");
      break;
  }

  // combine increments and forces
  bool converged = false;
  if (combinedconvergence_ == INPAR::FPSI::bop_and)
    converged = convinc and convfres;
  else if (combinedconvergence_ == INPAR::FPSI::bop_or)
    converged = convinc or convfres;
  else
    dserror("Something went terribly wrong with binary operator!");

  // return things
  return converged;
}  // Converged()


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::BuildConvergenceNorms()
{
  rhs_->Norm2(&normofrhs_);
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


  rhs_porostruct = Extractor().ExtractVector(rhs_, structure_block_);
  rhs_porofluid = Extractor().ExtractVector(rhs_, porofluid_block_);
  rhs_porofluidvelocity = PoroField()->FluidField()->ExtractVelocityPart(rhs_porofluid);
  rhs_porofluidpressure = PoroField()->FluidField()->ExtractPressurePart(rhs_porofluid);
  rhs_porointerface =
      FPSICoupl()->PoroFluidFpsiVelPresExtractor()->ExtractCondVector(rhs_porofluid);

  rhs_fluid = Extractor().ExtractVector(rhs_, fluid_block_);
  //  Teuchos::RCP<const Epetra_Vector> rhs_fullfluid = Teuchos::rcp(new
  //  Epetra_Vector(*FluidField()->DofRowMap())); Teuchos::RCP<const Epetra_Vector> rhs_fsi =
  //  Teuchos::rcp(new
  //  Epetra_Vector(*FluidField()->Interface()->Map(FLD::UTILS::MapExtractor::cond_fsi),true));
  //  rhs_fullfluid = LINALG::MergeVector(rhs_fluid,rhs_fsi,false);

  rhs_fluidvelocity = FluidField()->ExtractVelocityPart(rhs_fluid);
  rhs_fluidpressure = FluidField()->ExtractPressurePart(rhs_fluid);
  rhs_fluidinterface = FPSICoupl()->FluidFpsiVelPresExtractor()->ExtractCondVector(rhs_fluid);

  rhs_ale = Extractor().ExtractVector(rhs_, ale_i_block_);  // Extractor().ExtractVector(rhs_, 2);

  rhs_porostruct->Norm2(&normrhsporostruct_);
  rhs_fluid->Norm2(&normrhsfluid_);
  rhs_fluidvelocity->Norm2(&normrhsfluidvelocity_);
  rhs_fluidpressure->Norm2(&normrhsfluidpressure_);
  rhs_porofluidvelocity->Norm2(&normrhsporofluidvelocity_);
  rhs_porofluidpressure->Norm2(&normrhsporofluidpressure_);
  rhs_porointerface->Norm2(&normrhsporointerface_);
  rhs_fluidinterface->Norm2(&normrhsfluidinterface_);
  rhs_ale->Norm2(&normrhsale_);

  // get length of the porostructural, porofluid, fluid and ale vector
  sqrtnfv_ = rhs_fluidvelocity->GlobalLength();  // correct length here
  sqrtnfp_ = rhs_fluidpressure->GlobalLength();
  sqrtnpfv_ = rhs_porofluidvelocity->GlobalLength();
  sqrtnpfp_ = rhs_porofluidpressure->GlobalLength();
  sqrtnps_ = rhs_porostruct->GlobalLength();
  sqrtna_ = rhs_ale->GlobalLength();
  sqrtnall_ = sqrtnfv_ + sqrtnfp_ + sqrtnpfv_ + sqrtnpfp_ + sqrtnps_ + sqrtna_;

  sqrtnfv_ = sqrt(sqrtnfv_);
  sqrtnfp_ = sqrt(sqrtnfp_);
  sqrtnpfv_ = sqrt(sqrtnpfv_);
  sqrtnpfp_ = sqrt(sqrtnpfp_);
  sqrtnps_ = sqrt(sqrtnps_);
  sqrtna_ = sqrt(sqrtna_);
  sqrtnall_ = sqrt(sqrtnall_);

  if (islinesearch_ == false) iterinc_->Norm2(&normofiterinc_);

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

  iterincporostruct = Extractor().ExtractVector(iterinc_, structure_block_);
  iterincporofluid = Extractor().ExtractVector(iterinc_, porofluid_block_);
  iterincporofluidvelocity = PoroField()->FluidField()->ExtractVelocityPart(iterincporofluid);
  iterincporofluidpressure = PoroField()->FluidField()->ExtractPressurePart(iterincporofluid);
  iterincporointerface =
      FPSICoupl()->PoroFluidFpsiVelPresExtractor()->ExtractCondVector(iterincporofluid);

  iterincfluid = Extractor().ExtractVector(iterinc_, fluid_block_);

  //  Teuchos::RCP<const Epetra_Vector> iterincfullfluid = Teuchos::rcp(new
  //  Epetra_Vector(*FluidField()->DofRowMap())); Teuchos::RCP<const Epetra_Vector> iterincfsi =
  //  Teuchos::rcp(new
  //  Epetra_Vector(*FluidField()->Interface()->Map(FLD::UTILS::MapExtractor::cond_fsi),true));
  //  iterincfullfluid = LINALG::MergeVector(iterincfluid,iterincfsi,false);

  iterincfluidvelocity = FluidField()->ExtractVelocityPart(iterincfluid);
  iterincfluidpressure = FluidField()->ExtractPressurePart(iterincfluid);
  iterincfluidinterface = FPSICoupl()->FluidFpsiVelPresExtractor()->ExtractCondVector(iterincfluid);

  iterincale = Extractor().ExtractVector(iterinc_, ale_i_block_);

  iterincporostruct->Norm2(&normofiterincporostruct_);
  iterincporofluid->Norm2(&normofiterincporofluid_);
  iterincporofluidvelocity->Norm2(&normofiterincporofluidvelocity_);
  iterincporofluidpressure->Norm2(&normofiterincporofluidpressure_);
  iterincporointerface->Norm2(&normofiterincporointerface_);

  iterincfluid->Norm2(&normofiterincfluid_);
  iterincfluidvelocity->Norm2(&normofiterincfluidvelocity_);
  iterincfluidpressure->Norm2(&normofiterincfluidpressure_);
  iterincfluidinterface->Norm2(&normofiterincfluidinterface_);

  iterincale->Norm2(&normofiterincale_);

  // Get Norm1 of dof values for each field

  Teuchos::RCP<const Epetra_Vector> porofluiddof;
  Teuchos::RCP<const Epetra_Vector> porofluidvelocity;
  Teuchos::RCP<const Epetra_Vector> porofluidpressure;
  Teuchos::RCP<const Epetra_Vector> porostructdisplacements;
  Teuchos::RCP<const Epetra_Vector> fluiddof;
  Teuchos::RCP<const Epetra_Vector> fluidvelocity;
  Teuchos::RCP<const Epetra_Vector> fluidpressure;
  Teuchos::RCP<const Epetra_Vector> aledisplacements;
  Teuchos::RCP<const Epetra_Vector> porointerface;
  Teuchos::RCP<const Epetra_Vector> fluidinterface;

  porofluiddof = PoroField()->FluidField()->Velnp();
  porofluidvelocity = PoroField()->FluidField()->ExtractVelocityPart(porofluiddof);
  porofluidpressure = PoroField()->FluidField()->ExtractPressurePart(porofluiddof);
  porostructdisplacements = PoroField()->StructureField()->Dispnp();
  fluiddof = FluidField()->Velnp();
  fluidvelocity = FluidField()->ExtractVelocityPart(fluiddof);
  fluidpressure = FluidField()->ExtractPressurePart(fluiddof);
  aledisplacements = AleField()->Dispnp();


  porofluidvelocity->Norm1(&norm1_pfv_);
  porofluidpressure->Norm1(&norm1_pfp_);
  porostructdisplacements->Norm1(&norm1_ps_);
  fluidvelocity->Norm1(&norm1_fv_);
  fluidpressure->Norm1(&norm1_fp_);
  aledisplacements->Norm1(&norm1_a_);
  norm1_alldof_ = norm1_pfv_ + norm1_pfp_ + norm1_ps_ + norm1_fv_ + norm1_fp_ + norm1_a_;

  // add small number to avoid division by 0 --> division by 10^-10 results anyway in 'not
  // converged'
  norm1_alldof_ += 1e-10;
  norm1_pfv_ += 1e-10;
  norm1_pfp_ += 1e-10;
  norm1_ps_ += 1e-10;
  norm1_fv_ += 1e-10;
  norm1_fp_ += 1e-10;
  norm1_a_ += 1e-10;

  return;
}  // BuildConvergenceNorms

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::PrintNewtonIter()
{
  // print to standard out
  // replace myrank_ here general by Comm().MyPID()
  if ((Comm().MyPID() == 0) and printscreen_ and (Step() % printscreen_ == 0) and printiter_)
  {
    if (iter_ == 1) PrintNewtonIterHeader(stdout);
    PrintNewtonIterText(stdout);
  }

  // print to error file
  if (printerrfile_ and printiter_)
  {
    if (iter_ == 1) PrintNewtonIterHeader(errfile_);
    PrintNewtonIterText(errfile_);
  }

  return;
}  // PrintNewtonIter()

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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
    case INPAR::FPSI::absoluteconvergencenorm_sys_split:
      oss << std::setw(11) << "abs-res_s";
      break;
    case INPAR::FPSI::relativconvergencenorm_sys:
      oss << std::setw(11) << "rel-res_s";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case INPAR::FPSI::absoluteconvergencenorm:
      oss << std::setw(11) << "abs-inc";
      break;
    case INPAR::FPSI::absoluteconvergencenorm_sys_split:
      oss << std::setw(11) << "abs-inc_s";
      break;
    case INPAR::FPSI::relativconvergencenorm_sys:
      oss << std::setw(11) << "rel-inc_s";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch (normtypefres_)
  {
    case INPAR::FPSI::absoluteconvergencenorm:
    case INPAR::FPSI::absoluteconvergencenorm_sys_split:
    case INPAR::FPSI::relativconvergencenorm_sys:
      oss << std::setw(16) << "poro-s-res";
      // oss <<std::setw(15)<< "abs-f-res";
      oss << std::setw(15) << "poro-fvel-res";
      oss << std::setw(15) << "poro-fpres-res";
      oss << std::setw(15) << "fld-fvel-res";
      oss << std::setw(15) << "fld-fpres-res";
      // oss <<std::setw(15)<< "pinterface-res";
      // oss <<std::setw(15)<< "finterface-res";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case INPAR::FPSI::absoluteconvergencenorm:
    case INPAR::FPSI::absoluteconvergencenorm_sys_split:
      oss << std::setw(15) << "poro-s-inc";
      // oss <<std::setw(15)<< "abs-f-inc";
      oss << std::setw(16) << "poro-fvel-inc";
      oss << std::setw(16) << "poro-fpres-inc";
      oss << std::setw(15) << "fld-fvel-inc";
      oss << std::setw(15) << "fld-fpres-inc";
      oss << std::setw(11) << "ale-inc";
      oss << std::setw(14) << "poro-int-inc";
      oss << std::setw(14) << "fld-int-inc";
      break;
    case INPAR::FPSI::relativconvergencenorm_sys:
      oss << std::setw(15) << "poro-s-inc";
      // oss <<std::setw(15)<< "abs-f-inc";
      oss << std::setw(16) << "poro-fvel-inc";
      oss << std::setw(16) << "poro-fpres-inc";
      oss << std::setw(15) << "fld-fvel-inc";
      oss << std::setw(15) << "fld-fpres-inc";
      oss << std::setw(11) << "ale-inc";
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
  if (ofile == NULL) dserror("no ofile available");
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
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normofrhs_;
      break;
    case INPAR::FPSI::absoluteconvergencenorm_sys_split:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normofrhs_ / sqrtnall_;
      break;
    case INPAR::FPSI::relativconvergencenorm_sys:
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case INPAR::FPSI::absoluteconvergencenorm:
      oss << std::setw(12) << std::setprecision(5) << std::scientific << normofiterinc_;
      break;
    case INPAR::FPSI::relativconvergencenorm_sys:
      oss << std::setw(12) << std::setprecision(5) << std::scientific
          << normofiterinc_ / norm1_alldof_ * sqrtnall_;
      break;
    case INPAR::FPSI::absoluteconvergencenorm_sys_split:
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch (normtypefres_)
  {
    case INPAR::FPSI::absoluteconvergencenorm:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporostruct_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporofluidvelocity_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporofluidpressure_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsfluidvelocity_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsfluidpressure_;
      // oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporointerface_;
      // oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsfluidinterface_;
      break;
    case INPAR::FPSI::absoluteconvergencenorm_sys_split:
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normrhsporostruct_ / sqrtnps_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normrhsporofluidvelocity_ / sqrtnpfv_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normrhsporofluidpressure_ / sqrtnpfp_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normrhsfluidvelocity_ / sqrtnfv_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normrhsfluidpressure_ / sqrtnfp_;
      // oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporointerface_;
      // oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsfluidinterface_;
      break;
    case INPAR::FPSI::relativconvergencenorm_sys:
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case INPAR::FPSI::absoluteconvergencenorm:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normofiterincporostruct_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincporofluidvelocity_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincporofluidpressure_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincfluidvelocity_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincfluidpressure_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normofiterincale_;
      oss << std::setw(13) << std::setprecision(5) << std::scientific
          << normofiterincfluidinterface_;
      oss << std::setw(14) << std::setprecision(5) << std::scientific
          << normofiterincporointerface_;
      break;
    case INPAR::FPSI::relativconvergencenorm_sys:
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincporostruct_ / norm1_ps_ * sqrtnps_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincporofluidvelocity_ / norm1_pfv_ * sqrtnpfv_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincporofluidpressure_ / norm1_pfp_ * sqrtnpfp_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincfluidvelocity_ / norm1_fv_ * sqrtnfv_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincfluidpressure_ / norm1_fp_ * sqrtnfp_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincale_ / norm1_a_ * sqrtna_;
      break;
    case INPAR::FPSI::absoluteconvergencenorm_sys_split:
    default:
      dserror("You should not turn up here.");
      break;
  }

  // add solution time
  // oss << std::setw(14) << std::setprecision(2) << std::scientific << timer_.ElapsedTime();

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile == NULL) dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;

}  // PrintN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::SetupNewton()
{
  // initialise equilibrium loop and norms
  iter_ = 1;
  normofrhs_ = 0.0;
  normofiterinc_ = 0.0;
  normrhsfluid_ = 0.0;
  normofiterincfluid_ = 0.0;
  normrhsfluidvelocity_ = 0.0;
  normofiterincfluidvelocity_ = 0.0;
  normrhsporostruct_ = 0.0;
  normofiterincporostruct_ = 0.0;
  normofiterincporofluid_ = 0.0;
  normrhsfluidpressure_ = 0.0;
  normofiterincfluidpressure_ = 0.0;
  normrhsporofluidvelocity_ = 0.0;
  normofiterincporofluidvelocity_ = 0.0;
  normrhsporofluidpressure_ = 0.0;
  normofiterincporofluidpressure_ = 0.0;
  normrhsporointerface_ = 0.0;
  normofiterincporointerface_ = 0.0;
  normrhsfluidinterface_ = 0.0;
  normofiterincfluidinterface_ = 0.0;
  sqrtnall_ = 1;
  sqrtnfv_ = 1;
  sqrtnfp_ = 1;
  sqrtnpfv_ = 1;
  sqrtnpfp_ = 1;
  sqrtnps_ = 1;
  sqrtna_ = 1;
  norm1_alldof_ = 1.0;
  norm1_fv_ = 1.0;
  norm1_fp_ = 1.0;
  norm1_pfv_ = 1.0;
  norm1_pfp_ = 1.0;
  norm1_ps_ = 1.0;
  norm1_a_ = 1.0;


  // incremental solution vector with length of all dofs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  iterinc_->PutScalar(0.0);

  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  zeros_->PutScalar(0.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::FPSIFDCheck()
{
  // FD check is nice to check your linearisations, but be aware that we did not linearize following
  // terms:
  // gridvelocity in convective term, dispnp for the stabilization
  active_FD_check_ = true;
  // Set states
  PoroField()->FluidField()->Discretization()->ClearState();

  PoroField()->FluidField()->Discretization()->SetState(
      0, "dispnp", PoroField()->FluidField()->Dispnp());
  PoroField()->FluidField()->Discretization()->SetState(
      0, "gridv", PoroField()->FluidField()->GridVel());
  PoroField()->FluidField()->Discretization()->SetState(
      0, "dispn", PoroField()->FluidField()->Dispn());
  PoroField()->FluidField()->Discretization()->SetState(
      0, "veln", PoroField()->FluidField()->Veln());
  PoroField()->FluidField()->Discretization()->SetState(
      0, "velaf", PoroField()->FluidField()->Velnp());
  PoroField()->FluidField()->Discretization()->SetState(
      0, "velnp", PoroField()->FluidField()->Velnp());

  FluidField()->Discretization()->ClearState();

  FluidField()->Discretization()->SetState(0, "dispnp", FluidField()->Dispnp());
  FluidField()->Discretization()->SetState(0, "gridv", FluidField()->GridVel());
  FluidField()->Discretization()->SetState(0, "dispn", FluidField()->Dispn());
  FluidField()->Discretization()->SetState(0, "veln", FluidField()->Veln());
  FluidField()->Discretization()->SetState(0, "velaf", FluidField()->Velnp());
  FluidField()->Discretization()->SetState(0, "velnp", FluidField()->Velnp());

  // define and set toggle parameter delta
  const double delta = 1e-8;

  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fpsidynparams = problem->FPSIDynamicParams();
  int columntocheck = fpsidynparams.get<int>("FDCheck_column");
  int rowtocheck = fpsidynparams.get<int>("FDCheck_row");
  //////////////////////////////////////////////////////////////7
  // matrices and vectors
  //////////////////////////////////////////////////////////////
  // build artificial iteration increment
  Teuchos::RCP<Epetra_Vector> iterinc = LINALG::CreateVector(*DofRowMap(), true);
  const int dofs = iterinc->GlobalLength();
  iterinc->PutScalar(0.0);
  iterinc->ReplaceGlobalValue(0, 0, delta);

  // build approximated FD stiffness matrix
  Teuchos::RCP<Epetra_CrsMatrix> stiff_approx = LINALG::CreateMatrix(*DofRowMap(), 81);

  // store old rhs
  Teuchos::RCP<Epetra_Vector> rhs_old = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));
  rhs_old->Update(1.0, *rhs_, 0.0);

  Teuchos::RCP<Epetra_Vector> rhs_copy = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

  Teuchos::RCP<LINALG::SparseMatrix> sparse_copy =
      Teuchos::rcp(new LINALG::SparseMatrix(sparse->EpetraMatrix(), LINALG::Copy));


  std::cout << "\n****************** FPSI finite difference check ******************" << std::endl;

  int dof_poro_struct = (PoroField()->StructureField()->Discretization()->NumGlobalNodes()) * 3;
  int dof_poro_fluid = (PoroField()->FluidField()->Discretization()->NumGlobalNodes()) * 4;
  int dof_fluid = (FluidField()->Discretization()->NumGlobalNodes()) * 4;
  int dof_ale = (AleField()->Discretization()->NumGlobalNodes()) * 3;

  std::cout << "poro structure field has " << dof_poro_struct << " DOFs" << std::endl;
  std::cout << "poro fluid field has     " << dof_poro_fluid << " DOFs" << std::endl;
  std::cout << "fluid field has          " << dof_fluid << " DOFs" << std::endl;
  std::cout << "ale field has            " << dof_ale << " DOFs" << std::endl;
  std::cout << "in total                 " << dofs << " DOFs" << std::endl;


  for (int i_loc = 0; i_loc < dofs; ++i_loc)  // loop over columns
  {
    int i = DofRowMap()->GID(i_loc);
    int im1 = DofRowMap()->GID(i_loc - 1);
    int ip1 = DofRowMap()->GID(i_loc + 1);
    std::cout << i << "... " << std::flush;
    if (CombinedDBCMap()->MyGID(i) or FluidField()->Interface()->FSICondMap()->MyGID(i))
    {
      iterinc->ReplaceGlobalValue(i, 0, 0.0);
    }
    Evaluate(iterinc);  // initial iterinc is varied at first dof (0-th element)
    SetupSystemMatrix();

    rhs_copy->Update(1.0, *rhs_, 0.0);

    iterinc_->PutScalar(0.0);  // Useful? depends on solver and more
    PoroField()->ClearPoroIterinc();
    LINALG::ApplyDirichlettoSystem(sparse_copy, iterinc_, rhs_copy, Teuchos::null, zeros_,
        *FluidField()->Interface()->FSICondMap());
    LINALG::ApplyDirichlettoSystem(
        sparse_copy, iterinc_, rhs_copy, Teuchos::null, zeros_, *CombinedDBCMap());

    rhs_copy->Update(-1.0, *rhs_old, 1.0);  // finite difference approximation of partial derivative
    rhs_copy->Scale(-1.0 / delta);

    if (i == columntocheck)
    {
      std::cout << "iterinc:  " << *iterinc << std::endl;
      std::cout << "rhs_old:  " << *rhs_old << std::endl;
      std::cout << "rhs_copy: " << *rhs_copy << std::endl;
      dserror("Stopped by FPSI - FDCheck!");
    }

    int* index = &i;
    for (int j_loc = 0; j_loc < dofs; ++j_loc)  // loop over rows
    {
      int j = DofRowMap()->GID(j_loc);
      double value = (*rhs_copy)[j_loc];
      stiff_approx->InsertGlobalValues(
          j, 1, &value, index);  // int InsertGlobalValues(int GlobalRow, int NumEntries, double*
                                 // Values, int* Indices);

    }  // j-loop (rows)

    if (not CombinedDBCMap()->MyGID(i) and not FluidField()->Interface()->FSICondMap()->MyGID(i))
      iterinc->ReplaceGlobalValue(i, 0, -delta);

    iterinc->ReplaceGlobalValue(im1, 0, 0.0);

    if (i_loc != dofs - 1) iterinc->ReplaceGlobalValue(ip1, 0, delta);

  }  // i-loop (columns)
  Evaluate(iterinc);
  SetupSystemMatrix();

  int err = stiff_approx->FillComplete();
  if (err) dserror("FD_Check: FillComplete failed with err-code: %d", err);

  Teuchos::RCP<LINALG::SparseMatrix> temp =
      Teuchos::rcp(new LINALG::SparseMatrix(stiff_approx, LINALG::Copy));

  Teuchos::RCP<Epetra_CrsMatrix> stiff_approx_sparse = temp->EpetraMatrix();

  Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = sparse_copy->EpetraMatrix();

  // calc error (subtraction of sparse_crs and stiff_approx_sparse)
  for (int i_loc = 0; i_loc < dofs; i_loc++)
  {
    int i = DofRowMap()->GID(i_loc);
    int length;
    int numentries = sparse_crs->NumGlobalEntries(i);
    std::vector<double> values(numentries);
    std::vector<int> indices(numentries);
    sparse_crs->ExtractGlobalRowCopy(i, numentries, length, &values[0], &indices[0]);

    for (int k = 0; k < numentries; k++)
    {
      values[k] = -values[k];
    }

    stiff_approx_sparse->SumIntoGlobalValues(i, numentries, &values[0], &indices[0]);
  }
  stiff_approx_sparse->FillComplete();
  sparse_crs->FillComplete();

  bool success = true;
  double error_max = 0.0;
  for (int i_loc = 0; i_loc < dofs; ++i_loc)
  {
    int i = DofRowMap()->GID(i_loc);
    if (not CombinedDBCMap()->MyGID(i) and
        not FluidField()->Interface()->FSICondMap()->MyGID(
            i))  // only if there is no dirichlet value on dof and dof is not condensed
    {
      for (int j_loc = 0; j_loc < dofs; ++j_loc)
      {
        int j = DofRowMap()->GID(j_loc);
        if (not CombinedDBCMap()->MyGID(j) and
            not FluidField()->Interface()->FSICondMap()->MyGID(j))
        {
          double stiff_approx_ij = 0.0;
          double sparse_ij = 0.0;
          double error_ij = 0.0;

          // get error_crs entry ij
          int errornumentries;
          int errorlength = stiff_approx_sparse->NumGlobalEntries(i);
          std::vector<double> errorvalues(errorlength);
          std::vector<int> errorindices(errorlength);
          stiff_approx_sparse->ExtractGlobalRowCopy(
              i, errorlength, errornumentries, &errorvalues[0], &errorindices[0]);
          for (int k = 0; k < errorlength; ++k)
          {
            if (errorindices[k] == j)
            {
              error_ij = errorvalues[k];
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
          sparse_crs->ExtractGlobalRowCopy(
              i, sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);
          for (int k = 0; k < sparselength; ++k)
          {
            if (sparseindices[k] == j)
            {
              sparse_ij = sparsevalues[k];
              break;
            }
            else
              sparse_ij = 0.0;
          }
          // get stiff_approx entry ijs
          int approxnumentries;
          int approxlength = stiff_approx->NumGlobalEntries(i);
          std::vector<double> approxvalues(approxlength);
          std::vector<int> approxindices(approxlength);
          stiff_approx->ExtractGlobalRowCopy(
              i, approxlength, approxnumentries, &approxvalues[0], &approxindices[0]);
          for (int k = 0; k < approxlength; ++k)
          {
            if (approxindices[k] == j)
            {
              stiff_approx_ij = approxvalues[k];
              break;
            }
            else
              stiff_approx_ij = 0.0;
          }

          double error = 0.0;
          if (abs(stiff_approx_ij) > 1e-3)
            error = error_ij / stiff_approx_ij;
          else if (abs(sparse_ij) > 1e-3)
            error = error_ij / sparse_ij;

          if (i == rowtocheck and j == columntocheck)
          {
            std::cout << "K_approx value: " << stiff_approx_ij << std::endl;
            std::cout << "K value : " << sparse_ij << std::endl;
            std::cout << "error : " << error << std::endl;
            std::cout << "error_ij : " << error_ij << std::endl;
            std::cout << "i : " << i << std::endl;
            std::cout << "j : " << j << std::endl;
          }

          if (abs(error) > abs(error_max)) error_max = abs(error);

          if ((abs(error) > 1e-4))
          {
            if (abs(error_ij) > 1e-4)
            {
              std::cout << "finite difference check failed entry (" << i << "," << j
                        << ")! stiff: " << sparse_ij << ", approx: " << stiff_approx_ij
                        << " ,abs. error: " << error_ij << " , rel. error: " << error << std::endl;

              success = false;
            }
          }
        }
      }
    }
  }  // loop over dofs of succes check

  if (success)
  {
    std::cout << "finite difference check successful, max. rel. error: " << error_max << std::endl;
    std::cout << "****************** finite difference check done ***************\n\n" << std::endl;
  }
  else
    std::cout << "FPSIFDCheck failed" << std::endl;
  // dserror("FPSIFDCheck failed");
  PoroField()->FluidField()->Discretization()->ClearState();
  FluidField()->Discretization()->ClearState();

  active_FD_check_ = false;
  return;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::ExtractColumnsfromSparse(Teuchos::RCP<Epetra_CrsMatrix> src,
    const Teuchos::RCP<const Epetra_Map>& colmap, Teuchos::RCP<Epetra_CrsMatrix> dst)
{
  dst->PutScalar(0.0);  // clear matrix
  int rows = src->NumGlobalRows();
  for (int row = 0; row < rows; ++row)
  {
    int g_row = src->RangeMap().GID(row);
    int numentries;
    int length = src->NumGlobalEntries(g_row);
    std::vector<double> values(length);
    std::vector<int> indices(length);
    src->ExtractGlobalRowCopy(g_row, length, numentries, &values[0], &indices[0]);
    for (int col = 0; col < length; ++col)  // loop over non-zero columns in active row
    {
      if (colmap->LID(indices[col]) != -1)
      {
        dst->InsertGlobalValues(
            g_row, 1, &values[col], &indices[col]);  // add column value of active row!
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::SetConductivity(double conduct)
{
  if (FPSICoupl() != Teuchos::null) FPSICoupl()->SetConductivity(conduct);
  conductivity_ = conduct;  // remove me...
}
