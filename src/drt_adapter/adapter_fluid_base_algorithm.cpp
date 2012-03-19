/*----------------------------------------------------------------------*/
/*!
\file adapter_fluid_base_algorithm.cpp

\brief Fluid Base Algorithm

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
 */
/*----------------------------------------------------------------------*/


#include "adapter_fluid_base_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_inpar/inpar_combust.H"
#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_poroelast.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "adapter_fluid_impl.H"
#include "adapter_fluid_projection.H"
#include "adapter_xfluid2_impl.H"
#include "adapter_fluid_genalpha.H"
#include "adapter_fluid_combust.H"
#include "adapter_fluid_lung.H"
#include "adapter_fluid_fluid_impl.H"
#include "adapter_fluid_poro.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidBaseAlgorithm::FluidBaseAlgorithm(const Teuchos::ParameterList& prbdyn, bool isale)
{
  SetupFluid(prbdyn, isale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidBaseAlgorithm::FluidBaseAlgorithm(
  const Teuchos::ParameterList& prbdyn,
  const Teuchos::RCP<DRT::Discretization> discret,
  const Teuchos::RCP<std::map<int,std::vector<int> > > pbcmapmastertoslave)
{
  SetupInflowFluid(prbdyn, discret,pbcmapmastertoslave);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::FluidBaseAlgorithm::~FluidBaseAlgorithm()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::SetupFluid(const Teuchos::ParameterList& prbdyn, bool& isale)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("ADAPTER::FluidBaseAlgorithm::SetupFluid");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  if (genprob.probtyp == prb_fluid_fluid_fsi)
    actdis = DRT::Problem::Instance()->Dis(genprob.numff,1);
  else
    actdis = DRT::Problem::Instance()->Dis(genprob.numff,0);


  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  RCP<map<int,vector<int> > > pbcmapmastertoslave
    =
    Teuchos::rcp(new map<int,vector<int> > ());

  if((genprob.probtyp != prb_fsi) and
     (genprob.probtyp != prb_combust))
  {
    PeriodicBoundaryConditions pbc(actdis);
    pbc.UpdateDofsForPeriodicBoundaryConditions();

    pbcmapmastertoslave = pbc.ReturnAllCoupledColNodes();
  }

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->HaveDofs())
  {
    if (genprob.probtyp == prb_fsi_xfem or
        genprob.probtyp == prb_fluid_xfem or
        genprob.probtyp == prb_fluid_xfem2 or
        genprob.probtyp == prb_combust)
    {
      actdis->FillComplete(false,false,false);
    }
    else
    {
      actdis->FillComplete();
    }
  }

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter> output =
    rcp(new IO::DiscretizationWriter(actdis));
  if (genprob.probtyp != prb_fsi_xfem and
      genprob.probtyp != prb_fluid_xfem and
      //genprob.probtyp != prb_fluid_xfem2 and
      genprob.probtyp != prb_combust and
      genprob.probtyp != prb_fluid_fluid and
      genprob.probtyp != prb_fluid_fluid_ale and
      genprob.probtyp != prb_fluid_fluid_fsi)
  {
    output->WriteMesh(0,0.0);
  }

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  //const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::null;

  switch(DRT::INPUT::IntegralValue<int>(fdyn,"MESHTYING"))
  {
    case INPAR::FLUID::condensed_bmat:
    case INPAR::FLUID::sps_pc:
    {
      // meshtying fluid (formulation as saddle point problem)

      const Teuchos::ParameterList& mshparams = DRT::Problem::Instance()->MeshtyingAndContactParams();
      const int mshsolver = mshparams.get<int>("LINEAR_SOLVER");             // meshtying solver (with block preconditioner, e.g. BGS 2x2)
      const int fluidsolver = fdyn.get<int>("LINEAR_SOLVER");           // fluid solver
      const int fluidpressuresolver = fdyn.get<int>("SIMPLER_SOLVER");  // fluid pressure solver
      if (mshsolver == (-1))
        dserror("no linear solver defined for fluid meshtying problem. Please set LINEAR_SOLVER in MESHTYING AND CONTACT to a valid number!");
      if (fluidsolver == (-1))
        dserror("no linear solver defined for fluid meshtying problem. Please set LINEAR_SOLVER in FLUID DYNAMIC to a valid number! This solver is used within block preconditioner (e.g. BGS2x2) as \"Inverse 1\".");
      if (fluidpressuresolver == (-1))
        dserror("no linear solver defined for fluid meshtying problem. Please set SIMPLER_SOLVER in FLUID DYNAMIC to a valid number! This solver is used within block preconditioner (e.g. BGS2x2) as \"Inverse 2\".");

      // check, if meshtying solver is used with a valid block preconditioner
      const int azprectype
        = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(
            DRT::Problem::Instance()->SolverParams(mshsolver),
            "AZPREC"
            );

      // plausibility check
      switch (azprectype)
      {
        case INPAR::SOLVER::azprec_BGS2x2:      // block preconditioners, that are implemented in BACI
        case INPAR::SOLVER::azprec_CheapSIMPLE:
          break;
        case INPAR::SOLVER::azprec_BGSnxn:      // block preconditioners from Teko
        case INPAR::SOLVER::azprec_TekoSIMPLE:
        {
    #ifdef HAVE_TEKO
          // check if structural solver and thermal solver are Stratimikos based (Teko expects stratimikos)
          int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(fluidsolver), "SOLVER");
          if (solvertype != INPAR::SOLVER::stratimikos_amesos &&
              solvertype != INPAR::SOLVER::stratimikos_aztec  &&
              solvertype != INPAR::SOLVER::stratimikos_belos)
          dserror("Teko expects a STRATIMIKOS solver object in SOLVER %i", fluidsolver);

          solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(DRT::Problem::Instance()->SolverParams(fluidpressuresolver), "SOLVER");
          if (solvertype != INPAR::SOLVER::stratimikos_amesos &&
              solvertype != INPAR::SOLVER::stratimikos_aztec  &&
              solvertype != INPAR::SOLVER::stratimikos_belos)
            dserror("Teko expects a STRATIMIKOS solver object in SOLVER %i",fluidpressuresolver);
    #else
          dserror("Teko preconditioners only available with HAVE_TEKO flag for TRILINOS_DEV (>Q1/2011)");
    #endif
        }
        break;
        default:
              dserror("Block Gauss-Seidel BGS2x2 preconditioner expected for fluid meshtying problem. Please set AZPREC to BGS2x2 in solver block %i",mshsolver);
              break;
      }

      // create solver objects
      solver =
        rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(mshsolver),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));

      solver->Params().set<bool>("MESHTYING",true);   // mark it as meshtying problem

      // set Inverse blocks for block preconditioner
      solver->PutSolverParamsToSubParams("Inverse1",
          DRT::Problem::Instance()->SolverParams(fluidsolver));

      solver->PutSolverParamsToSubParams("Inverse2",
          DRT::Problem::Instance()->SolverParams(fluidpressuresolver));
    }
    break;
    case INPAR::FLUID::condensed_smat:
    case INPAR::FLUID::condensed_bmat_merged:
    case INPAR::FLUID::sps_coupled:
    case INPAR::FLUID::coupling_iontransport_laplace:
    { // meshtying (no saddle point problem)
      const Teuchos::ParameterList& mshparams = DRT::Problem::Instance()->MeshtyingAndContactParams();
      const int mshsolver = mshparams.get<int>("LINEAR_SOLVER");             // meshtying solver (with block preconditioner, e.g. BGS 2x2)
      if (mshsolver == (-1))
        dserror("no linear solver defined for fluid meshtying problem. Please set LINEAR_SOLVER in MESHTYING AND CONTACT to a valid number!");

      solver =
        rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(mshsolver),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));
    }
    break;
    case INPAR::FLUID::no_meshtying: // no meshtying -> use FLUID SOLVER
    default:
    {
      // default: create solver using the fluid solver params from FLUID SOLVER block

      // get the solver number used for linear fluid solver
      const int linsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for fluid problem. Please set LINEAR_SOLVER in FLUID DYNAMIC to a valid number!");
      solver =
        rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));
      break;
    }
  }

  // compute null space information
  if (genprob.probtyp != prb_fsi_xfem and
      genprob.probtyp != prb_fluid_xfem and
      genprob.probtyp != prb_fluid_xfem2 and
      genprob.probtyp != prb_combust)
  {
    switch(DRT::INPUT::IntegralValue<int>(fdyn,"MESHTYING"))
    {
      case INPAR::FLUID::condensed_bmat:
      {
        // meshtying fluid
        // block Gauss Seidel or SIMPLE preconditioners
        actdis->ComputeNullSpaceIfNecessary(solver->Params().sublist("Inverse1"),true);
        actdis->ComputeNullSpaceIfNecessary(solver->Params().sublist("Inverse2"),true);
      }
      break;
      case INPAR::FLUID::sps_pc:
      {
        // meshtying fluid
        // pure saddle point problem. only SIMPLE type preconditioners available
        // the standard nullspace is computed for the constraint block within the SIMPLE block preconditioner class
        actdis->ComputeNullSpaceIfNecessary(solver->Params().sublist("Inverse1"),true);
      }
      break;
      default:
        // no block matrix
        actdis->ComputeNullSpaceIfNecessary(solver->Params(),true);
    }
  }

  // -------------------------------------------------------------------
  // create a second solver for SIMPLER preconditioner if chosen from input
  // -------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<int>(fdyn,"SIMPLER"))
  {
    // add Inverse1 block for velocity dofs
    Teuchos::ParameterList& inv1 = solver->Params().sublist("Inverse1");
    inv1 = solver->Params();
    inv1.remove("SIMPLER",false); // not necessary
    inv1.remove("Inverse1",false);

    // get the solver number used for SIMPLER SOLVER
    const int linsolvernumber_simpler = fdyn.get<int>("SIMPLER_SOLVER");
    if (linsolvernumber_simpler == (-1))
      dserror("no SIMPLER_SOLVER number set for fluid problem solved with SIMPLER. Please set SIMPLER_SOLVER in FLUID DYNAMIC to a valid number!");
    // add Inverse2 block for pressure dofs
    solver->PutSolverParamsToSubParams("Inverse2", DRT::Problem::Instance()->SolverParams(linsolvernumber_simpler));
    // use CheapSIMPLE preconditioner (hardwired, change me for others)
    solver->Params().sublist("CheapSIMPLE Parameters").set("Prec Type","CheapSIMPLE");
    solver->Params().set("FLUID",true);
  }

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  RCP<ParameterList> fluidtimeparams = rcp(new ParameterList());

  // --------------------provide info about periodic boundary conditions
  fluidtimeparams->set<RCP<map<int,vector<int> > > >("periodic bc",pbcmapmastertoslave);

  fluidtimeparams->set<int>("Simple Preconditioner",DRT::INPUT::IntegralValue<int>(fdyn,"SIMPLER"));
  // fluidtimeparams->set<int>("AMG(BS) Preconditioner",DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(DRT::Problem::Instance()->FluidSolverParams(),"AZPREC")); // probably can be removed

  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  fluidtimeparams->set<int>("number of velocity degrees of freedom" ,probsize.get<int>("DIM"));

  // physical type of fluid flow (incompressible, Boussinesq Approximation, varying density, loma)
  fluidtimeparams->set<int>("Physical Type",
      DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn,"PHYSICAL_TYPE"));

  // check correct setting
  if ((genprob.probtyp == prb_thermo_fsi or genprob.probtyp == prb_loma) and
      DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn,"PHYSICAL_TYPE") 
      != INPAR::FLUID::loma)
    dserror("Input parameter PHYSICAL_TYPE in section FLUID DYNAMIC needs to be 'Loma' for low-Mach-number flow and Thermo-fluid-structure interaction!");


  // -------------------------------------------------- time integration
  // note: here, the values are taken out of the problem-dependent ParameterList prbdyn
  // (which also can be fluiddyn itself!)

  // the default time step size
  fluidtimeparams->set<double> ("time step size"      ,prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  fluidtimeparams->set<double> ("total time"          ,prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  fluidtimeparams->set<int>    ("max number timesteps",prbdyn.get<int>("NUMSTEP"));

  // -------- additional parameters in list for generalized-alpha scheme
#if 1
  // parameter alpha_M
  fluidtimeparams->set<double> ("alpha_M", fdyn.get<double>("ALPHA_M"));
  // parameter alpha_F
  fluidtimeparams->set<double> ("alpha_F", fdyn.get<double>("ALPHA_F"));
  // parameter gamma
  fluidtimeparams->set<double> ("gamma",   fdyn.get<double>("GAMMA"));
#else
  // parameter alpha_M
  fluidtimeparams->set<double> ("alpha_M", 1.-prbdyn.get<double>("ALPHA_M"));
  // parameter alpha_F
  fluidtimeparams->set<double> ("alpha_F", 1.-prbdyn.get<double>("ALPHA_F"));
  // parameter gamma
  fluidtimeparams->set<double> ("gamma",   prbdyn.get<double>("GAMMA"));
#endif

  // ---------------------------------------------- nonlinear iteration
  // type of predictor
  fluidtimeparams->set<string>          ("predictor"                 ,fdyn.get<string>("PREDICTOR"));
  // set linearisation scheme
  fluidtimeparams->set<int>("Linearisation", DRT::INPUT::IntegralValue<INPAR::FLUID::LinearisationAction>(fdyn,"NONLINITER"));
  // set bool flag "Newton true or false" for combustion formulation and XFEM
  //fluidtimeparams->set<bool>("Use reaction terms for linearisation",
  //                           DRT::INPUT::IntegralValue<INPAR::FLUID::LinearisationAction>(fdyn,"NONLINITER")== INPAR::FLUID::Newton);
  // maximum number of nonlinear iteration steps
  fluidtimeparams->set<int>             ("max nonlin iter steps"     ,fdyn.get<int>("ITEMAX"));
  // maximum number of nonlinear iteration steps for initial stationary solution
  fluidtimeparams->set<int>             ("max nonlin iter steps init stat sol",fdyn.get<int>("INITSTATITEMAX"));
  // stop nonlinear iteration when both incr-norms are below this bound
  fluidtimeparams->set<double>          ("tolerance for nonlin iter" ,fdyn.get<double>("CONVTOL"));
  // set convergence check
  fluidtimeparams->set<string>          ("CONVCHECK"  ,fdyn.get<string>("CONVCHECK"));
  // set adaptive linear solver tolerance
  fluidtimeparams->set<bool>            ("ADAPTCONV",DRT::INPUT::IntegralValue<int>(fdyn,"ADAPTCONV")==1);
  fluidtimeparams->set<double>          ("ADAPTCONV_BETTER",fdyn.get<double>("ADAPTCONV_BETTER"));
  fluidtimeparams->set<bool>            ("INFNORMSCALING", (DRT::INPUT::IntegralValue<int>(fdyn,"INFNORMSCALING")==1));

  // ----------------------------------------------- restart and output
  // restart
  fluidtimeparams->set ("write restart every", prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  fluidtimeparams->set ("write solution every", prbdyn.get<int>("UPRES"));
  // flag for writing stresses
  fluidtimeparams->set ("write stresses"  ,DRT::INPUT::IntegralValue<int>(ioflags,"FLUID_STRESS"));
  // flag for writing wall shear stress
  fluidtimeparams->set ("write wall shear stresses"  ,DRT::INPUT::IntegralValue<int>(ioflags,"FLUID_WALL_SHEAR_STRESS"));

  // ---------------------------------------------------- lift and
  // drag
  fluidtimeparams->set<int>("liftdrag",DRT::INPUT::IntegralValue<int>(fdyn,"LIFTDRAG"));

  // -----------evaluate error for test flows with analytical solutions
  INPAR::FLUID::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn,"INITIALFIELD");
  fluidtimeparams->set<int>("eval err for analyt sol", initfield);

  // ------------------------------------------ form of convective term
  fluidtimeparams->set<string> ("form of convective term", fdyn.get<string>("CONVFORM"));

  // ------------------------------------ potential Neumann inflow terms
  fluidtimeparams->set<string> ("Neumann inflow",fdyn.get<string>("NEUMANNINFLOW"));


  // ------------------------------------ potential reduced_D 3D coupling method
  fluidtimeparams->set<string> ("Strong 3D_redD coupling",fdyn.get<string>("STRONG_REDD_3D_COUPLING_TYPE"));

  //--------------------------------------mesh tying for fluid
  fluidtimeparams->set<int>("MESHTYING",
      DRT::INPUT::IntegralValue<int>(fdyn,"MESHTYING"));

  //--------------------------------------analytical error evaluation
  fluidtimeparams->set<int>("calculate error",
      Teuchos::getIntegralValue<int>(fdyn,"CALCERROR"));

  // -----------------------sublist containing stabilization parameters
  fluidtimeparams->sublist("STABILIZATION")=fdyn.sublist("STABILIZATION");

  // -----------------------------get also scatra stabilization sublist
  const Teuchos::ParameterList& scatradyn =
    DRT::Problem::Instance()->ScalarTransportDynamicParams();
  fluidtimeparams->sublist("SCATRA STABILIZATION")=scatradyn.sublist("STABILIZATION");

  // ------------------------------------------- Robin scheme parameters
  if (genprob.probtyp == prb_fsi)
  {
    INPAR::FSI::PartitionedCouplingMethod method =
      DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(prbdyn,"PARTITIONED");
  }

  // --------------------------sublist containing turbulence parameters
  {
    fluidtimeparams->sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");
    fluidtimeparams->sublist("SUBGRID VISCOSITY")=fdyn.sublist("SUBGRID VISCOSITY");
    fluidtimeparams->sublist("MULTIFRACTAL SUBGRID SCALES")=fdyn.sublist("MULTIFRACTAL SUBGRID SCALES");
    fluidtimeparams->sublist("TURBULENT INFLOW")=fdyn.sublist("TURBULENT INFLOW");

    fluidtimeparams->sublist("TURBULENCE MODEL").set<string>("statistics outfile",DRT::Problem::Instance()->OutputControlFile()->FileName());
  }

  // ----------------------------- add some loma specific parameters
  // get also scatra stabilization sublist
  const Teuchos::ParameterList& lomadyn =
    DRT::Problem::Instance()->LOMAControlParams();
  fluidtimeparams->sublist("LOMA").set<bool>("update material",DRT::INPUT::IntegralValue<int>(lomadyn,"SGS_MATERIAL_UPDATE"));

  // ----------------------------------------------- XFEM related stuff
  {
    const Teuchos::ParameterList& xdyn = DRT::Problem::Instance()->XFEMGeneralParams();
    // Axel's unused input parameters
    fluidtimeparams->sublist("XFEM").set<bool>("DLM_condensation", DRT::INPUT::IntegralValue<int>(xdyn,"DLM_CONDENSATION")==1 );

    // interface coupling method
    fluidtimeparams->sublist("XFEM").set<int>("EMBEDDED_BOUNDARY", DRT::INPUT::IntegralValue<INPAR::XFEM::BoundaryIntegralType>(xdyn, "EMBEDDED_BOUNDARY"));

    fluidtimeparams->sublist("XFEM").set<int>("XFLUID_BOUNDARY", DRT::INPUT::IntegralValue<INPAR::XFEM::MovingBoundary>(xdyn, "XFLUID_BOUNDARY"));

    fluidtimeparams->sublist("XFEM").set<int>("INTERFACE_VEL_INITIAL", DRT::INPUT::IntegralValue<INPAR::XFEM::InterfaceInitVel>(xdyn, "INTERFACE_VEL_INITIAL"));
    fluidtimeparams->sublist("XFEM").set<int>("VEL_INIT_FUNCT_NO", xdyn.get<int>("VEL_INIT_FUNCT_NO"));
    fluidtimeparams->sublist("XFEM").set<int>("INTERFACE_VEL", DRT::INPUT::IntegralValue<INPAR::XFEM::InterfaceVel>(xdyn, "INTERFACE_VEL"));
    fluidtimeparams->sublist("XFEM").set<int>("VEL_FUNCT_NO", xdyn.get<int>("VEL_FUNCT_NO"));

    fluidtimeparams->sublist("XFEM").set<int>("INTERFACE_DISP", DRT::INPUT::IntegralValue<INPAR::XFEM::InterfaceDisplacement>(xdyn, "INTERFACE_DISP"));
    fluidtimeparams->sublist("XFEM").set<int>("DISP_FUNCT_NO", xdyn.get<int>("DISP_FUNCT_NO"));
    fluidtimeparams->sublist("XFEM").set<int>("DISP_CURVE_NO", xdyn.get<int>("DISP_CURVE_NO"));


    fluidtimeparams->sublist("XFEM").set<int>("COUPLING_STRATEGY", DRT::INPUT::IntegralValue<INPAR::XFEM::CouplingStrategy>(xdyn, "COUPLING_STRATEGY"));
    fluidtimeparams->sublist("XFEM").set<double>("Nitsche_stab", xdyn.get<double>("Nitsche_stab"));
    fluidtimeparams->sublist("XFEM").set<double>("Nitsche_stab_conv", xdyn.get<double>("Nitsche_stab_conv"));
    fluidtimeparams->sublist("XFEM").set<int>("MAX_NUM_DOFSETS", xdyn.get<int>("MAX_NUM_DOFSETS"));
    fluidtimeparams->sublist("XFEM").set<string>("VOLUME_GAUSS_POINTS_BY", xdyn.get<string>("VOLUME_GAUSS_POINTS_BY"));
    fluidtimeparams->sublist("XFEM").set<string>("BOUNDARY_GAUSS_POINTS_BY", xdyn.get<string>("BOUNDARY_GAUSS_POINTS_BY"));

    // GMSH solution output
    fluidtimeparams->sublist("XFEM").set<int>("GMSH_DEBUG_OUT",        DRT::INPUT::IntegralValue<int>(xdyn, "GMSH_DEBUG_OUT"));
    fluidtimeparams->sublist("XFEM").set<int>("GMSH_DEBUG_OUT_SCREEN", DRT::INPUT::IntegralValue<int>(xdyn, "GMSH_DEBUG_OUT_SCREEN"));
    fluidtimeparams->sublist("XFEM").set<int>("GMSH_SOL_OUT",          DRT::INPUT::IntegralValue<int>(xdyn, "GMSH_SOL_OUT"));
    fluidtimeparams->sublist("XFEM").set<int>("GMSH_DISCRET_OUT",      DRT::INPUT::IntegralValue<int>(xdyn, "GMSH_DISCRET_OUT"));
    fluidtimeparams->sublist("XFEM").set<int>("GMSH_CUT_OUT",          DRT::INPUT::IntegralValue<int>(xdyn, "GMSH_CUT_OUT"));

    fluidtimeparams->sublist("XFEM").set<int>("MONOLITHIC_XFFSI_APPROACH", DRT::INPUT::IntegralValue<INPAR::XFEM::Monolithic_xffsi_Approach>(xdyn,"MONOLITHIC_XFFSI_APPROACH"));
    fluidtimeparams->sublist("XFEM").set<int>("XFLUIDFLUID_TIMEINT", DRT::INPUT::IntegralValue<INPAR::XFEM::XFluidFluidTimeInt>(xdyn, "XFLUIDFLUID_TIMEINT"));
    fluidtimeparams->sublist("XFEM").set<int>("RELAXING_ALE", xdyn.get<int>("RELAXING_ALE"));
  }

  if( genprob.probtyp == prb_fluid_xfem2 or
      genprob.probtyp == prb_fsi_xfem      )
  {
    // check some input configurations
    INPAR::XFEM::MovingBoundary xfluid_mov_bound = DRT::INPUT::get<INPAR::XFEM::MovingBoundary>(fluidtimeparams->sublist("XFEM"), "XFLUID_BOUNDARY");
    INPAR::XFEM::InterfaceVel   interf_vel       = DRT::INPUT::get<INPAR::XFEM::InterfaceVel>(fluidtimeparams->sublist("XFEM"),"INTERFACE_VEL");
    INPAR::XFEM::InterfaceDisplacement interf_disp    = DRT::INPUT::get<INPAR::XFEM::InterfaceDisplacement>(fluidtimeparams->sublist("XFEM"),"INTERFACE_DISP");

    if(genprob.probtyp == prb_fsi_xfem and xfluid_mov_bound!= INPAR::XFEM::XFSIMovingBoundary) dserror("choose xfsi_moving_boundary!!! for prb_fsi_xfem");
    if(genprob.probtyp == prb_fluid_xfem2 and xfluid_mov_bound== INPAR::XFEM::XFSIMovingBoundary) dserror("do not choose xfsi_moving_boundary!!! for prb_fluid_xfem2");
    if(genprob.probtyp == prb_fsi_xfem and interf_disp != INPAR::XFEM::interface_disp_by_fsi) dserror("choose interface_disp_by_fsi for prb_fsi_xfem");
    if(genprob.probtyp == prb_fluid_xfem2 and interf_disp == INPAR::XFEM::interface_disp_by_fsi ) dserror("do not choose interface_disp_by_fsi for prb_fluid_xfem2");
    if(genprob.probtyp == prb_fsi_xfem and interf_vel != INPAR::XFEM::interface_vel_by_disp ) dserror("do you want to use !interface_vel_by_disp for prb_fsi_xfem?");
  }


  // --------------------------sublist for combustion-specific fluid parameters
  /* This sublist COMBUSTION FLUID contains parameters for the fluid field
   * which are only relevant for a combustion problem.                 07/08 henke */
  if (genprob.probtyp == prb_combust)
  {
    fluidtimeparams->sublist("COMBUSTION FLUID")=prbdyn.sublist("COMBUSTION FLUID");
    // parameter COMBUSTTYPE from sublist COMBUSTION FLUID is also added to sublist XFEM
    fluidtimeparams->sublist("XFEM").set<int>("combusttype", DRT::INPUT::IntegralValue<INPAR::COMBUST::CombustionType>(prbdyn.sublist("COMBUSTION FLUID"),"COMBUSTTYPE"));
    fluidtimeparams->sublist("XFEM").set<bool>("GMSH_OUTPUT", DRT::INPUT::IntegralValue<bool>(prbdyn,"GMSH_OUTPUT"));
  }

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::FLUID::TimeIntegrationScheme timeint = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn,"TIMEINTEGR");

  // sanity checks and default flags
  if (genprob.probtyp == prb_fsi or
      genprob.probtyp == prb_fsi_lung or
      genprob.probtyp == prb_gas_fsi or
      genprob.probtyp == prb_biofilm_fsi or
      genprob.probtyp == prb_thermo_fsi or
      genprob.probtyp == prb_fluid_fluid_fsi)
  {
    // in case of FSI calculations we do not want a stationary fluid solver
    if (timeint == INPAR::FLUID::timeint_stationary)
      dserror("Stationary fluid solver not allowed for FSI.");

    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

    fluidtimeparams->set<bool>("interface second order", DRT::INPUT::IntegralValue<int>(fsidyn,"SECONDORDER"));
    fluidtimeparams->set<bool>("shape derivatives",
                               DRT::INPUT::IntegralValue<int>(fsidyn,"SHAPEDERIVATIVES"));

    const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");

    if (coupling == fsi_iter_monolithicfluidsplit or
        coupling == fsi_iter_monolithicstructuresplit or
        coupling == fsi_iter_lung_monolithicstructuresplit or
        coupling == fsi_iter_lung_monolithicstructuresplit or
        coupling == fsi_iter_constr_monolithicstructuresplit or
        coupling == fsi_iter_constr_monolithicfluidsplit or
        coupling == fsi_iter_mortar_monolithicstructuresplit or
        coupling == fsi_iter_mortar_monolithicfluidsplit or
        coupling == fsi_iter_fluidfluid_monolithicstructuresplit)
    {
      // there are a couple of restrictions in monolithic FSI
      fluidtimeparams->set<bool>("do explicit predictor",false);
    }
  }
  if (genprob.probtyp == prb_freesurf)
  {
    // in case of FSI calculations we do not want a stationary fluid solver
    if (timeint == INPAR::FLUID::timeint_stationary)
      dserror("Stationary fluid solver not allowed for Freesurface problem.");

    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

    fluidtimeparams->set<bool>("interface second order", DRT::INPUT::IntegralValue<int>(fsidyn,"SECONDORDER"));
    fluidtimeparams->set<bool>("shape derivatives",
                               DRT::INPUT::IntegralValue<int>(fsidyn,"SHAPEDERIVATIVES"));

    const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling == fsi_iter_monolithicfluidsplit or
        coupling == fsi_iter_monolithicstructuresplit)
    {
      // there are a couple of restrictions in monolithic Freesurface Algorithm
      fluidtimeparams->set<bool>("do explicit predictor",false);
    }
  }
  // sanity checks and default flags
  if (genprob.probtyp == prb_fsi_xfem or
      genprob.probtyp == prb_fluid_xfem or
      genprob.probtyp == prb_fluid_xfem2)
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    fluidtimeparams->set<bool>("interface second order", DRT::INPUT::IntegralValue<int>(fsidyn,"SECONDORDER"));

    const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
    if (coupling == fsi_iter_monolithicfluidsplit or
        coupling == fsi_iter_monolithicstructuresplit)
    {
      // there are a couple of restrictions in monolithic FSI
      dserror("XFEM and monolithic FSI not tested!");
      fluidtimeparams->set<bool>("do explicit predictor",false);
    }
    else
    {
      fluidtimeparams->set<bool>("do explicit predictor",true);
    }
  }

  if (genprob.probtyp == prb_elch)
  {
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    fluidtimeparams->set<bool>("interface second order", DRT::INPUT::IntegralValue<int>(fsidyn,"SECONDORDER"));
  }

  if (genprob.probtyp == prb_poroelast)
  {
    fluidtimeparams->set<double>("initporosity",prbdyn.get<double>("INITPOROSITY"));
    fluidtimeparams->set<bool>("poroelast",true);
  }
  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  if(timeint == INPAR::FLUID::timeint_stationary or
     timeint == INPAR::FLUID::timeint_one_step_theta or
     timeint == INPAR::FLUID::timeint_bdf2 or
     timeint == INPAR::FLUID::timeint_afgenalpha or
     timeint == INPAR::FLUID::timeint_npgenalpha
    )
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for
    // one-step-theta/BDF2/af-generalized-alpha/stationary scheme
    // -----------------------------------------------------------------
    // type of time-integration (or stationary) scheme
    fluidtimeparams->set<int>("time int algo",timeint);
    // parameter theta for time-integration schemes
    fluidtimeparams->set<double>           ("theta"                    ,fdyn.get<double>("THETA"));
    // number of steps for potential start algorithm
    fluidtimeparams->set<int>              ("number of start steps"    ,fdyn.get<int>("NUMSTASTEPS"));
    // parameter theta for potential start algorithm
    fluidtimeparams->set<double>           ("start theta"              ,fdyn.get<double>("START_THETA"));
    // parameter for grid velocity interpolation
    fluidtimeparams->set<int>              ("order gridvel"            ,fdyn.get<int>("GRIDVEL"));

    fluidtimeparams->set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());
    bool dirichletcond = true;
    if (genprob.probtyp == prb_fsi or
        genprob.probtyp == prb_fsi_lung or
        genprob.probtyp == prb_gas_fsi or
        genprob.probtyp == prb_biofilm_fsi or
        genprob.probtyp == prb_thermo_fsi or
        genprob.probtyp == prb_fluid_fluid_fsi)
    {
      // FSI input parameters
      const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
      const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
      if (coupling == fsi_iter_monolithicfluidsplit or
          coupling == fsi_iter_monolithicstructuresplit or
          coupling == fsi_iter_lung_monolithicstructuresplit or
          coupling == fsi_iter_lung_monolithicfluidsplit or
          coupling == fsi_iter_constr_monolithicstructuresplit or
          coupling == fsi_iter_constr_monolithicfluidsplit or
          coupling == fsi_iter_mortar_monolithicstructuresplit or
          coupling == fsi_iter_mortar_monolithicfluidsplit or
          coupling == fsi_iter_fluidfluid_monolithicstructuresplit)
      {
        dirichletcond = false;
      }
      else
      {
        dserror("unknown coupling algorithm!");
      }
    }

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom

    if (genprob.probtyp == prb_fsi_xfem or
             genprob.probtyp == prb_fluid_xfem2 )
    {
      RCP<DRT::Discretization> soliddis = DRT::Problem::Instance()->Dis(genprob.numsf,0);

      fluid_ = rcp( new ADAPTER::XFluid2Impl( actdis, soliddis, solver, fluidtimeparams, output));
    }
    else if (genprob.probtyp == prb_combust)
    {
      fluid_ = rcp(new ADAPTER::FluidCombust(actdis, solver, fluidtimeparams, output));
    }
    else if (genprob.probtyp == prb_fluid_fluid_ale
             or genprob.probtyp == prb_fluid_fluid)
    {
      RCP<DRT::Discretization> embfluiddis  =  DRT::Problem::Instance()->Dis(genprob.numff,1);
      bool monolithicfluidfluidfsi = false;
      fluid_ = rcp(new ADAPTER::FluidFluidImpl(embfluiddis,actdis,solver,fluidtimeparams,isale,dirichletcond,monolithicfluidfluidfsi));
    }
    else if (genprob.probtyp == prb_fluid_fluid_fsi)
    {
      RCP<DRT::Discretization> bgfluiddis  =  DRT::Problem::Instance()->Dis(genprob.numff,0);
      const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
      const int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
      bool monolithicfluidfluidfsi;
      if(coupling == fsi_iter_fluidfluid_monolithicstructuresplit)
        monolithicfluidfluidfsi = true;
      else
        monolithicfluidfluidfsi = false;

      fluid_ = rcp(new ADAPTER::FluidFluidImpl(actdis,bgfluiddis,solver,fluidtimeparams,isale,dirichletcond,monolithicfluidfluidfsi));
    }
    else
    {
      // what the hell is this???
      RCP<Fluid> tmpfluid;
      tmpfluid = rcp(new ADAPTER::FluidImpl(actdis, solver, fluidtimeparams, output, isale, dirichletcond));

      if (genprob.probtyp == prb_fsi_lung)
        fluid_ = rcp(new FluidLung(rcp(new FluidWrapper(tmpfluid))));
      else if (genprob.probtyp == prb_poroelast)
        fluid_ = rcp(new ADAPTER::FluidPoro(rcp(new FluidWrapper(tmpfluid))));
      else
        fluid_ = tmpfluid;  // < default?? TODO improve me!
    }
  }
  else if (timeint == INPAR::FLUID::timeint_gen_alpha)
  {
    fluidtimeparams->set<int>("time int algo",INPAR::FLUID::timeint_gen_alpha);

    // -------------------------------------------------------------------
    // no additional parameters in list for generalized-alpha scheme
    // -------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    //------------------------------------------------------------------
    RCP<Fluid> tmpfluid;
    tmpfluid = rcp(new ADAPTER::FluidGenAlpha(actdis, solver, fluidtimeparams, output, isale , pbcmapmastertoslave));
    if (genprob.probtyp == prb_fsi_lung)
      fluid_ = rcp(new FluidLung(rcp(new FluidWrapper(tmpfluid))));
    else
      fluid_ = tmpfluid;
  }
  else
  {
    dserror("Unknown time integration for fluid\n");
  }

  // set initial field by given function
  // we do this here, since we have direct access to all necessary parameters
  if(initfield != INPAR::FLUID::initfield_zero_field)
  {
    int startfuncno = fdyn.get<int>("STARTFUNCNO");
    if (initfield != INPAR::FLUID::initfield_field_by_function and
        initfield != INPAR::FLUID::initfield_disturbed_field_from_function)
    {
      startfuncno=-1;
    }
    fluid_->SetInitialFlowField(initfield,startfuncno);
  }

  if (genprob.probtyp == prb_fluid_topopt)
    fluid_->Output();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidBaseAlgorithm::SetupInflowFluid(
  const Teuchos::ParameterList& prbdyn,
  const Teuchos::RCP<DRT::Discretization> discret,
  const Teuchos::RCP<std::map<int,std::vector<int> > > pbcmapmastertoslave)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("ADAPTER::FluidBaseAlgorithm::SetupFluid");
  Teuchos::TimeMonitor monitor(*t);

  // the inflow computation can only deal with standard fluid problems so far
  // extensions for xfluid, fsi or combust problems have to be added if necessary
  // they should not pose any additional problem
  // meshtying or xfem related parameters are not supported, yet
  if (genprob.probtyp != prb_fluid)
     dserror("Only fluid problems supported! Read comment and add your problem type!");

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!discret->HaveDofs())
  {
    dserror("FillComplete shouldn't be necessary!");
    //discret->FillComplete();
  }

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter> output = rcp(new IO::DiscretizationWriter(discret));

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  //const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  if (discret->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // get the solver number used for linear fluid solver
  const int linsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for fluid problem. Please set LINEAR_SOLVER in FLUID DYNAMIC to a valid number!");
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                           discret->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));

  discret->ComputeNullSpaceIfNecessary(solver->Params(),true);


  // -------------------------------------------------------------------
  // create a second solver for SIMPLER preconditioner if chosen from input
  // -------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<int>(fdyn,"SIMPLER"))
  {
    // add Inverse1 block for velocity dofs
    Teuchos::ParameterList& inv1 = solver->Params().sublist("Inverse1");
    inv1 = solver->Params();
    inv1.remove("SIMPLER",false); // not necessary
    inv1.remove("Inverse1",false);

    // get the solver number used for SIMPLER SOLVER
    const int linsolvernumber_simpler = fdyn.get<int>("SIMPLER_SOLVER");
    if (linsolvernumber_simpler == (-1))
      dserror("no SIMPLER_SOLVER number set for fluid problem solved with SIMPLER. Please set SIMPLER_SOLVER in FLUID DYNAMIC to a valid number!");
    // add Inverse2 block for pressure dofs
    solver->PutSolverParamsToSubParams("Inverse2", DRT::Problem::Instance()->SolverParams(linsolvernumber_simpler));
    // use CheapSIMPLE preconditioner (hardwired, change me for others)
    solver->Params().sublist("CheapSIMPLE Parameters").set("Prec Type","CheapSIMPLE");
    solver->Params().set("FLUID",true);
  }

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  RCP<ParameterList> fluidtimeparams = rcp(new ParameterList());

  // --------------------provide info about periodic boundary conditions
  fluidtimeparams->set<RCP<map<int,vector<int> > > >("periodic bc",pbcmapmastertoslave);

  fluidtimeparams->set<int>("Simple Preconditioner",DRT::INPUT::IntegralValue<int>(fdyn,"SIMPLER"));
  // fluidtimeparams->set<int>("AMG(BS) Preconditioner",DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(DRT::Problem::Instance()->FluidSolverParams(),"AZPREC")); // probably can be removed

  // -------------------------------------- number of degrees of freedom
  // number of space dimensions
  fluidtimeparams->set<int>("number of velocity degrees of freedom" ,probsize.get<int>("DIM"));

  // physical type of fluid flow (incompressible, Boussinesq Approximation, varying density, loma)
  fluidtimeparams->set<int>("Physical Type",
      DRT::INPUT::IntegralValue<INPAR::FLUID::PhysicalType>(fdyn,"PHYSICAL_TYPE"));


  // -------------------------------------------------- time integration
  // note: here, the values are taken out of the problem-dependent ParameterList prbdyn
  // (which also can be fluiddyn itself!)

  // the default time step size
  fluidtimeparams->set<double> ("time step size"      ,prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  fluidtimeparams->set<double> ("total time"          ,prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  fluidtimeparams->set<int>    ("max number timesteps",prbdyn.get<int>("NUMSTEP"));

  // -------- additional parameters in list for generalized-alpha scheme
  // parameter alpha_M
  fluidtimeparams->set<double> ("alpha_M", fdyn.get<double>("ALPHA_M"));
  // parameter alpha_F
  fluidtimeparams->set<double> ("alpha_F", fdyn.get<double>("ALPHA_F"));
  // parameter gamma
  fluidtimeparams->set<double> ("gamma",   fdyn.get<double>("GAMMA"));

  // ---------------------------------------------- nonlinear iteration
  // type of predictor
  fluidtimeparams->set<string>          ("predictor"                 ,fdyn.get<string>("PREDICTOR"));
  // set linearisation scheme
  fluidtimeparams->set<int>("Linearisation", DRT::INPUT::IntegralValue<INPAR::FLUID::LinearisationAction>(fdyn,"NONLINITER"));
  // set bool flag "Newton true or false" for combustion formulation and XFEM
  //fluidtimeparams->set<bool>("Use reaction terms for linearisation",
  //                           DRT::INPUT::IntegralValue<INPAR::FLUID::LinearisationAction>(fdyn,"NONLINITER")== INPAR::FLUID::Newton);
  // maximum number of nonlinear iteration steps
  fluidtimeparams->set<int>             ("max nonlin iter steps"     ,fdyn.get<int>("ITEMAX"));
  // maximum number of nonlinear iteration steps for initial stationary solution
  fluidtimeparams->set<int>             ("max nonlin iter steps init stat sol",fdyn.get<int>("INITSTATITEMAX"));
  // stop nonlinear iteration when both incr-norms are below this bound
  fluidtimeparams->set<double>          ("tolerance for nonlin iter" ,fdyn.get<double>("CONVTOL"));
  // set convergence check
  fluidtimeparams->set<string>          ("CONVCHECK"  ,fdyn.get<string>("CONVCHECK"));
  // set adaptoive linear solver tolerance
  fluidtimeparams->set<bool>            ("ADAPTCONV",DRT::INPUT::IntegralValue<int>(fdyn,"ADAPTCONV")==1);
  fluidtimeparams->set<double>          ("ADAPTCONV_BETTER",fdyn.get<double>("ADAPTCONV_BETTER"));
  fluidtimeparams->set<bool>            ("INFNORMSCALING", (DRT::INPUT::IntegralValue<int>(fdyn,"INFNORMSCALING")==1));

  // ----------------------------------------------- restart and output
  // restart
  fluidtimeparams->set ("write restart every", prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  fluidtimeparams->set ("write solution every", prbdyn.get<int>("UPRES"));
  // flag for writing stresses
  fluidtimeparams->set ("write stresses"  ,DRT::INPUT::IntegralValue<int>(ioflags,"FLUID_STRESS"));
  // flag for writing wall shear stress
  fluidtimeparams->set ("write wall shear stresses"  ,DRT::INPUT::IntegralValue<int>(ioflags,"FLUID_WALL_SHEAR_STRESS"));

  // ---------------------------------------------------- lift and drag
  fluidtimeparams->set<int>("liftdrag",DRT::INPUT::IntegralValue<int>(fdyn,"LIFTDRAG"));

  // -----------evaluate error for test flows with analytical solutions
  INPAR::FLUID::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn,"INITIALFIELD");
  fluidtimeparams->set<int>("eval err for analyt sol", initfield);

  // ------------------------------------------ form of convective term
  fluidtimeparams->set<string> ("form of convective term", fdyn.get<string>("CONVFORM"));

  // ------------------------------------ potential Neumann inflow terms
  fluidtimeparams->set<string> ("Neumann inflow",fdyn.get<string>("NEUMANNINFLOW"));

  //--------------------------------------mesh tying for fluid
  fluidtimeparams->set<int>("MESHTYING",
      DRT::INPUT::IntegralValue<int>(fdyn,"MESHTYING"));

  //--------------------------------------analytical error evaluation
  fluidtimeparams->set<int>("calculate error",
      Teuchos::getIntegralValue<int>(fdyn,"CALCERROR"));

  // -----------------------sublist containing stabilization parameters
  fluidtimeparams->sublist("STABILIZATION")=fdyn.sublist("STABILIZATION");

  // -----------------------------get also scatra stabilization sublist
  const Teuchos::ParameterList& scatradyn =
    DRT::Problem::Instance()->ScalarTransportDynamicParams();
  fluidtimeparams->sublist("SCATRA STABILIZATION")=scatradyn.sublist("STABILIZATION");

  // --------------------------sublist containing turbulence parameters
  {
    fluidtimeparams->sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");
    fluidtimeparams->sublist("SUBGRID VISCOSITY")=fdyn.sublist("SUBGRID VISCOSITY");
    fluidtimeparams->sublist("MULTIFRACTAL SUBGRID SCALES")=fdyn.sublist("MULTIFRACTAL SUBGRID SCALES");
    fluidtimeparams->sublist("TURBULENT INFLOW")=fdyn.sublist("TURBULENT INFLOW");

    fluidtimeparams->sublist("TURBULENCE MODEL").set<string>("statistics outfile",DRT::Problem::Instance()->OutputControlFile()->FileName());

    // overwrite canonical flow parameters by inflow type
    fluidtimeparams->sublist("TURBULENCE MODEL").set<string>("CANONICAL_FLOW",fdyn.sublist("TURBULENT INFLOW").get<string>("CANONICAL_INFLOW"));
    fluidtimeparams->sublist("TURBULENCE MODEL").set<string>("HOMDIR",fdyn.sublist("TURBULENT INFLOW").get<string>("INFLOW_HOMDIR"));
    fluidtimeparams->sublist("TURBULENCE MODEL").set<int>("DUMPING_PERIOD",fdyn.sublist("TURBULENT INFLOW").get<int>("INFLOW_DUMPING_PERIOD"));
    fluidtimeparams->sublist("TURBULENCE MODEL").set<int>("SAMPLING_START",fdyn.sublist("TURBULENT INFLOW").get<int>("INFLOW_SAMPLING_START"));
    fluidtimeparams->sublist("TURBULENCE MODEL").set<int>("SAMPLING_STOP",fdyn.sublist("TURBULENT INFLOW").get<int>("INFLOW_SAMPLING_STOP"));
  }

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::FLUID::TimeIntegrationScheme timeint = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn,"TIMEINTEGR");

  // -------------------------------------------------------------------
  // additional parameters and algorithm call depending on respective
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  if(timeint == INPAR::FLUID::timeint_stationary or
     timeint == INPAR::FLUID::timeint_one_step_theta or
     timeint == INPAR::FLUID::timeint_bdf2 or
     timeint == INPAR::FLUID::timeint_afgenalpha
    )
  {
    // -----------------------------------------------------------------
    // set additional parameters in list for
    // one-step-theta/BDF2/af-generalized-alpha/stationary scheme
    // -----------------------------------------------------------------
    // type of time-integration (or stationary) scheme
    fluidtimeparams->set<int>("time int algo",timeint);
    // parameter theta for time-integration schemes
    fluidtimeparams->set<double>           ("theta"                    ,fdyn.get<double>("THETA"));
    // number of steps for potential start algorithm
    fluidtimeparams->set<int>              ("number of start steps"    ,fdyn.get<int>("NUMSTASTEPS"));
    // parameter theta for potential start algorithm
    fluidtimeparams->set<double>           ("start theta"              ,fdyn.get<double>("START_THETA"));
    // parameter for grid velocity interpolation
    fluidtimeparams->set<int>              ("order gridvel"            ,fdyn.get<int>("GRIDVEL"));

    fluidtimeparams->set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

    bool dirichletcond = true;
    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    fluid_ = rcp(new ADAPTER::FluidImpl(discret, solver, fluidtimeparams, output, false, dirichletcond));

  }
  else if (timeint == INPAR::FLUID::timeint_gen_alpha)
  {
    fluidtimeparams->set<int>("time int algo",timeint);

    // -------------------------------------------------------------------
    // no additional parameters in list for generalized-alpha scheme
    // -------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    //------------------------------------------------------------------
    RCP<Fluid> tmpfluid; // what the hell is this???
    tmpfluid = rcp(new ADAPTER::FluidGenAlpha(discret, solver, fluidtimeparams, output, false , pbcmapmastertoslave));
    fluid_ = tmpfluid;
  }
  else
  {
    dserror("Unknown time integration for fluid\n");
  }

  // set initial field for inflow section by given function
  // we do this here, since we have direct access to all necessary parameters
  initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn.sublist("TURBULENT INFLOW"),"INITIALINFLOWFIELD");
  if(initfield != INPAR::FLUID::initfield_zero_field)
  {
    int startfuncno = fdyn.sublist("TURBULENT INFLOW").get<int>("INFLOWFUNC");
    if (initfield != INPAR::FLUID::initfield_field_by_function and
        initfield != INPAR::FLUID::initfield_disturbed_field_from_function)
    {
      startfuncno=-1;
    }
    fluid_->SetInitialFlowField(initfield,startfuncno);
  }

  return;
}
