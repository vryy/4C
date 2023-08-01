/*----------------------------------------------------------------------*/
/*! \file

\brief scalar transport field base algorithm

\level 1


*/
/*----------------------------------------------------------------------*/

#include "baci_adapter_scatra_base_algorithm.H"

#include "baci_inpar_ssi.H"
#include "baci_inpar_ssti.H"
#include "baci_inpar_sti.H"
#include "baci_io.H"
#include "baci_io_control.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_linear_solver_method_linalg.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// general time integration schemes
#include "baci_scatra_timint_bdf2.H"
#include "baci_scatra_timint_genalpha.H"
#include "baci_scatra_timint_ost.H"
#include "baci_scatra_timint_stat.H"

// HDG time integration schemes
#include "baci_scatra_resulttest_hdg.H"
#include "baci_scatra_timint_cardiac_monodomain_scheme_hdg.H"
#include "baci_scatra_timint_stat_hdg.H"

// loma specific files
#include "baci_scatra_timint_loma_genalpha.H"

// elch specific files
#include "baci_scatra_timint_elch_scheme.H"

// level set specific files
#include "baci_levelset_timint_ost.H"
#include "baci_levelset_timint_stat.H"

// cardiac monodomain specific files
#include "baci_scatra_timint_cardiac_monodomain_scheme.H"

// poro multiphase files
#include "baci_scatra_timint_poromulti.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraBaseAlgorithm::ScaTraBaseAlgorithm()
    : scatra_(Teuchos::null), issetup_(false), isinit_(false)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraBaseAlgorithm::Init(
    const Teuchos::ParameterList& prbdyn,        ///< parameter list for global problem
    const Teuchos::ParameterList& scatradyn,     ///< parameter list for scalar transport subproblem
    const Teuchos::ParameterList& solverparams,  ///< parameter list for scalar transport solver
    const std::string& disname,                  ///< name of scalar transport discretization
    const bool isale                             ///< ALE flag
)
{
  SetIsSetup(false);

  // setup scalar transport algorithm (overriding some dynamic parameters
  // with values specified in given problem-dependent ParameterList prbdyn)

  // -------------------------------------------------------------------
  // what's the current problem type?
  // -------------------------------------------------------------------
  auto probtype = DRT::Problem::Instance()->GetProblemType();

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  auto discret = DRT::Problem::Instance()->GetDis(disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!discret->Filled() or !discret->HaveDofs()) discret->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  auto output = discret->Writer();
  if (discret->NumGlobalElements() == 0)
    dserror("No elements in discretization %s", discret->Name().c_str());
  output->WriteMesh(0, 0.0);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // TODO: TAW use of solverparams???
  // change input parameter to solver number instead of parameter list?
  // -> no default paramter possible any more
  auto solver = Teuchos::rcp(new CORE::LINALG::Solver(
      solverparams, discret->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  // make a copy (inside an Teuchos::rcp) containing also all sublists
  auto scatratimeparams = Teuchos::rcp(new Teuchos::ParameterList(scatradyn));
  if (scatratimeparams == Teuchos::null) dserror("Instantiation of Teuchos::ParameterList failed!");

  // -------------------------------------------------------------------
  // overrule certain parameters for coupled problems
  // -------------------------------------------------------------------
  // the default time step size
  scatratimeparams->set<double>("TIMESTEP", prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  scatratimeparams->set<double>("MAXTIME", prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  scatratimeparams->set<int>("NUMSTEP", prbdyn.get<int>("NUMSTEP"));
  // restart
  scatratimeparams->set<int>("RESTARTEVRY", prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  scatratimeparams->set<int>("RESULTSEVRY", prbdyn.get<int>("RESULTSEVRY"));

  // -------------------------------------------------------------------
  // overrule flags for solid-based scalar transport!
  // (assumed disname = "scatra2" for solid-based scalar transport)
  // -------------------------------------------------------------------
  if (probtype == ProblemType::ac_fsi or probtype == ProblemType::biofilm_fsi or
      probtype == ProblemType::gas_fsi or probtype == ProblemType::fps3i or
      probtype == ProblemType::thermo_fsi)
  {
    // scatra1 (=fluid scalar) get's inputs from SCALAR TRANSPORT DYNAMIC/STABILIZATION, hence
    // nothing is to do here
    //    if (disname== "scatra1") //get's inputs from SCALAR TRANSPORT DYNAMIC/STABILIZATION

    if (disname == "scatra2")  // structure_scatra discretisation
    {
      // scatra2 (=structure scalar) get's inputs from FS3I DYNAMIC/STRUCTURE SCALAR STABILIZATION,
      // hence we have to replace it
      scatratimeparams->sublist("STABILIZATION") = prbdyn.sublist("STRUCTURE SCALAR STABILIZATION");
      scatratimeparams->set<std::string>(
          "CONVFORM", prbdyn.get<std::string>("STRUCTSCAL_CONVFORM"));

      // scatra2 get's in initial functions from FS3I DYNAMICS
      switch (
          DRT::INPUT::IntegralValue<INPAR::SCATRA::InitialField>(prbdyn, "STRUCTSCAL_INITIALFIELD"))
      {
        case INPAR::SCATRA::initfield_zero_field:
          scatratimeparams->set<std::string>("INITIALFIELD",
              "zero_field");  // we want zero initial conditions for the structure scalar
          scatratimeparams->set<int>("INITFUNCNO", -1);
          break;
        case INPAR::SCATRA::initfield_field_by_function:
          scatratimeparams->set<std::string>(
              "INITIALFIELD", "field_by_function");  // we want the same initial conditions for
                                                     // structure scalar as for the fluid scalar
          scatratimeparams->set<int>("INITFUNCNO", prbdyn.get<int>("STRUCTSCAL_INITFUNCNO"));
          break;
        default:
          dserror("Your STRUCTSCAL_INITIALFIELD type is not supported!");
          break;
      }

      // structure scatra does not require any Neumann inflow boundary conditions
      scatratimeparams->set<std::string>("NEUMANNINFLOW", "no");
    }
    else if (disname == "scatra1")  // fluid_scatra discretisation
    {
      // fluid scatra does not require any convective heat transfer boundary conditions
      scatratimeparams->set<std::string>("CONV_HEAT_TRANS", "no");
    }
  }

  // -------------------------------------------------------------------
  // list for extra parameters
  // (put here everything that is not available in scatradyn or its sublists)
  // -------------------------------------------------------------------
  auto extraparams = Teuchos::rcp(new Teuchos::ParameterList());

  // ------------------------------pointer to the error file (for output)
  extraparams->set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());

  // ----------------Eulerian or ALE formulation of transport equation(s)
  extraparams->set<bool>("isale", isale);

  // ------------------------------------get also fluid turbulence sublist
  const auto& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  extraparams->sublist("TURBULENCE MODEL") = fdyn.sublist("TURBULENCE MODEL");
  extraparams->sublist("SUBGRID VISCOSITY") = fdyn.sublist("SUBGRID VISCOSITY");
  extraparams->sublist("MULTIFRACTAL SUBGRID SCALES") = fdyn.sublist("MULTIFRACTAL SUBGRID SCALES");
  extraparams->sublist("TURBULENT INFLOW") = fdyn.sublist("TURBULENT INFLOW");

  // ------------------------------------get electromagnetic parameters
  extraparams->set<bool>("ELECTROMAGNETICDIFFUSION",
      DRT::INPUT::IntegralValue<int>(scatradyn, "ELECTROMAGNETICDIFFUSION"));
  extraparams->set<int>("EMDSOURCE", scatradyn.get<int>("EMDSOURCE"));

  // -------------------------------------------------------------------
  // algorithm construction depending on problem type and
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  auto timintscheme =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(scatradyn, "TIMEINTEGR");

  // low Mach number flow
  if (probtype == ProblemType::loma or probtype == ProblemType::thermo_fsi)
  {
    auto lomaparams =
        Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->LOMAControlParams()));
    switch (timintscheme)
    {
      case INPAR::SCATRA::timeint_gen_alpha:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::TimIntLomaGenAlpha(
            discret, solver, lomaparams, scatratimeparams, extraparams, output));
        break;
      }
      default:
        dserror("Unknown time integration scheme for loMa!");
        break;
    }
  }

  // electrochemistry
  else if (probtype == ProblemType::elch or
           ((probtype == ProblemType::ssi and
                Teuchos::getIntegralValue<INPAR::SSI::ScaTraTimIntType>(
                    DRT::Problem::Instance()->SSIControlParams(), "SCATRATIMINTTYPE") ==
                    INPAR::SSI::ScaTraTimIntType::elch) or
               (disname == "scatra" and
                   ((probtype == ProblemType::ssti and
                        Teuchos::getIntegralValue<INPAR::SSTI::ScaTraTimIntType>(
                            DRT::Problem::Instance()->SSTIControlParams(), "SCATRATIMINTTYPE") ==
                            INPAR::SSTI::ScaTraTimIntType::elch) or
                       (probtype == ProblemType::sti and
                           Teuchos::getIntegralValue<INPAR::STI::ScaTraTimIntType>(
                               DRT::Problem::Instance()->STIDynamicParams(), "SCATRATIMINTTYPE") ==
                               INPAR::STI::ScaTraTimIntType::elch)))))
  {
    auto elchparams =
        Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->ELCHControlParams()));

    switch (timintscheme)
    {
      case INPAR::SCATRA::timeint_one_step_theta:
      {
        if (DRT::INPUT::IntegralValue<bool>(elchparams->sublist("SCL"), "ADD_MICRO_MACRO_COUPLING"))
        {
          if (disname == "scatra")
          {
            scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntElchSCLOST(
                discret, solver, elchparams, scatratimeparams, extraparams, output));
          }
          else if (disname == "scatra_micro")
          {
            scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntElchOST(
                discret, solver, elchparams, scatratimeparams, extraparams, output));
          }
          else
            dserror("not identified");
        }
        else
        {
          scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntElchOST(
              discret, solver, elchparams, scatratimeparams, extraparams, output));
        }

        break;
      }
      case INPAR::SCATRA::timeint_bdf2:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntElchBDF2(
            discret, solver, elchparams, scatratimeparams, extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_gen_alpha:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntElchGenAlpha(
            discret, solver, elchparams, scatratimeparams, extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_stationary:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntElchStationary(
            discret, solver, elchparams, scatratimeparams, extraparams, output));
        break;
      }
      default:
        dserror("Unknown time integration scheme for electrochemistry!");
        break;
    }
  }

  // levelset and two phase flow
  else if (probtype == ProblemType::level_set or probtype == ProblemType::two_phase_flow or
           probtype == ProblemType::fluid_xfem_ls)
  {
    Teuchos::RCP<Teuchos::ParameterList> lsparams = Teuchos::null;
    switch (probtype)
    {
      case ProblemType::level_set:
        lsparams = Teuchos::rcp(new Teuchos::ParameterList(prbdyn));
        break;
      case ProblemType::two_phase_flow:
      {
        // Give access to smoothing parameter for levelset calculations.
        lsparams =
            Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->LevelSetControl()));
        lsparams->set<double>("INTERFACE_THICKNESS_TPF",
            prbdyn.sublist("SMEARED").get<double>("INTERFACE_THICKNESS"));
        // !!! no break !!!
      }
      default:
      {
        if (lsparams.is_null())
          lsparams =
              Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->LevelSetControl()));
        // overrule certain parameters for coupled problems
        // this has already been ensured for scatratimeparams, but has also been ensured for the
        // level-set parameter in a hybrid approach time step size
        lsparams->set<double>("TIMESTEP", prbdyn.get<double>("TIMESTEP"));
        // maximum simulation time
        lsparams->set<double>("MAXTIME", prbdyn.get<double>("MAXTIME"));
        // maximum number of timesteps
        lsparams->set<int>("NUMSTEP", prbdyn.get<int>("NUMSTEP"));
        // restart
        lsparams->set<int>("RESTARTEVRY", prbdyn.get<int>("RESTARTEVRY"));
        // solution output
        lsparams->set<int>("RESULTSEVRY", prbdyn.get<int>("RESULTSEVRY"));

        break;
      }
    }

    switch (timintscheme)
    {
      case INPAR::SCATRA::timeint_one_step_theta:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::LevelSetTimIntOneStepTheta(
            discret, solver, lsparams, scatratimeparams, extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_stationary:
      {
        // create instance of time integration class (call the constructor)
        switch (probtype)
        {
          case ProblemType::level_set:
          {
            dserror(
                "Stationary time integration scheme only supported for a selection of coupled "
                "level-set problems!");
            exit(EXIT_FAILURE);
          }
          default:
          {
            scatra_ = Teuchos::rcp(new SCATRA::LevelSetTimIntStationary(
                discret, solver, lsparams, scatratimeparams, extraparams, output));
            break;
          }
        }
        break;
      }
      case INPAR::SCATRA::timeint_gen_alpha:
      {
        switch (probtype)
        {
          case ProblemType::two_phase_flow:
          {
            std::cout << "\n\n\n WARNING: Level set algorithm does not yet support gen-alpha. You "
                         "thus get a standard Scatra!\n\n\n"
                      << std::endl;
            // create instance of time integration class (call the constructor)
            scatra_ = Teuchos::rcp(
                new SCATRA::TimIntGenAlpha(discret, solver, scatratimeparams, extraparams, output));
            break;
          }
          default:
            dserror("Unknown time-integration scheme for level-set problem");
            exit(EXIT_FAILURE);
        }

        break;
      }
      default:
        dserror("Unknown time-integration scheme for level-set problem");
        break;
    }  // switch(timintscheme)
  }

  // cardiac monodomain
  else if (probtype == ProblemType::cardiac_monodomain or
           (probtype == ProblemType::ssi and
               Teuchos::getIntegralValue<INPAR::SSI::ScaTraTimIntType>(
                   DRT::Problem::Instance()->SSIControlParams(), "SCATRATIMINTTYPE") ==
                   INPAR::SSI::ScaTraTimIntType::cardiac_monodomain))
  {
    auto cmonoparams = rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->EPControlParams()));

    // HDG implements all time stepping schemes within gen-alpha
    if (DRT::Problem::Instance()->SpatialApproximationType() ==
        ShapeFunctionType::shapefunction_hdg)
    {
      scatra_ = Teuchos::rcp(new SCATRA::TimIntCardiacMonodomainHDG(
          discret, solver, cmonoparams, scatratimeparams, extraparams, output));
    }
    else
    {
      switch (timintscheme)
      {
        case INPAR::SCATRA::timeint_gen_alpha:
        {
          // create instance of time integration class (call the constructor)
          scatra_ = Teuchos::rcp(new SCATRA::TimIntCardiacMonodomainGenAlpha(
              discret, solver, cmonoparams, scatratimeparams, extraparams, output));
          break;
        }
        case INPAR::SCATRA::timeint_one_step_theta:
        {
          // create instance of time integration class (call the constructor)
          scatra_ = Teuchos::rcp(new SCATRA::TimIntCardiacMonodomainOST(
              discret, solver, cmonoparams, scatratimeparams, extraparams, output));
          break;
        }
        case INPAR::SCATRA::timeint_bdf2:
        {
          // create instance of time integration class (call the constructor)
          scatra_ = Teuchos::rcp(new SCATRA::TimIntCardiacMonodomainBDF2(
              discret, solver, cmonoparams, scatratimeparams, extraparams, output));
          break;
        }
        default:
          dserror("Unknown time integration scheme for cardiac monodomain problem!");
          break;
      }  // switch(timintscheme)
    }
  }
  else if (probtype == ProblemType::poromultiphasescatra)
  {
    switch (timintscheme)
    {
      case INPAR::SCATRA::timeint_gen_alpha:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntPoroMultiGenAlpha(
            discret, solver, Teuchos::null, scatratimeparams, extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_one_step_theta:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntPoroMultiOST(
            discret, solver, Teuchos::null, scatratimeparams, extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_bdf2:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntPoroMultiBDF2(
            discret, solver, Teuchos::null, scatratimeparams, extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_stationary:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntPoroMultiStationary(
            discret, solver, Teuchos::null, scatratimeparams, extraparams, output));
        break;
      }
      default:
        dserror("Unknown time integration scheme for porous medium multiphase problem!");
        break;
    }  // switch(timintscheme)
  }
  // everything else
  else
  {
    // HDG implements all time stepping schemes within gen-alpha
    if (DRT::Problem::Instance()->SpatialApproximationType() ==
        ShapeFunctionType::shapefunction_hdg)
    {
      switch (timintscheme)
      {
        case INPAR::SCATRA::timeint_one_step_theta:
        case INPAR::SCATRA::timeint_bdf2:
        case INPAR::SCATRA::timeint_gen_alpha:
        {
          scatra_ = Teuchos::rcp(
              new SCATRA::TimIntHDG(discret, solver, scatratimeparams, extraparams, output));
          break;
        }
        case INPAR::SCATRA::timeint_stationary:
        {
          scatra_ = Teuchos::rcp(new SCATRA::TimIntStationaryHDG(
              discret, solver, scatratimeparams, extraparams, output));
          break;
        }
        default:
        {
          dserror("Unknown time-integration scheme for HDG scalar transport problem");
          break;
        }
      }
    }
    else
    {
      switch (timintscheme)
      {
        case INPAR::SCATRA::timeint_stationary:
        {
          // create instance of time integration class (call the constructor)
          scatra_ = Teuchos::rcp(
              new SCATRA::TimIntStationary(discret, solver, scatratimeparams, extraparams, output));
          break;
        }
        case INPAR::SCATRA::timeint_one_step_theta:
        {
          // create instance of time integration class (call the constructor)
          scatra_ = Teuchos::rcp(new SCATRA::TimIntOneStepTheta(
              discret, solver, scatratimeparams, extraparams, output));
          break;
        }
        case INPAR::SCATRA::timeint_bdf2:
        {
          // create instance of time integration class (call the constructor)
          scatra_ = Teuchos::rcp(
              new SCATRA::TimIntBDF2(discret, solver, scatratimeparams, extraparams, output));
          break;
        }
        case INPAR::SCATRA::timeint_gen_alpha:
        {
          // create instance of time integration class (call the constructor)
          scatra_ = Teuchos::rcp(
              new SCATRA::TimIntGenAlpha(discret, solver, scatratimeparams, extraparams, output));
          break;
        }
        default:
          dserror("Unknown time-integration scheme for scalar transport problem");
          break;
      }  // switch(timintscheme)
    }
  }

  // initialize scatra time integrator
  scatra_->Init();

  SetIsInit(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraBaseAlgorithm::Setup()
{
  CheckIsInit();

  // setup the time integrator
  scatra_->Setup();

  // get the parameter list
  auto scatradyn = scatra_->ScatraParameterList();
  // get the discretization
  auto discret = scatra_->Discretization();

  // -------------------------------------------------------------------
  // what's the current problem type?
  // -------------------------------------------------------------------
  auto probtype = DRT::Problem::Instance()->GetProblemType();

  // prepare fixing the null space for electrochemistry and sti
  if (probtype == ProblemType::elch or
      (probtype == ProblemType::sti and discret->Name() == "scatra" and
          Teuchos::getIntegralValue<INPAR::STI::ScaTraTimIntType>(
              DRT::Problem::Instance()->STIDynamicParams(), "SCATRATIMINTTYPE") ==
              INPAR::STI::ScaTraTimIntType::elch))
  {
    auto elchparams =
        Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->ELCHControlParams()));

    // create a 2nd solver for block-preconditioning if chosen from input
    if (DRT::INPUT::IntegralValue<int>(*elchparams, "BLOCKPRECOND"))
    {
      const auto& solver = scatra_->Solver();

      const int linsolvernumber = scatradyn->get<int>("LINEAR_SOLVER");
      const auto prec = Teuchos::getIntegralValue<INPAR::SOLVER::PreconditionerType>(
          DRT::Problem::Instance()->SolverParams(linsolvernumber), "AZPREC");
      if (prec != INPAR::SOLVER::PreconditionerType::cheap_simple)  // TODO adapt error message
      {
        dserror(
            "If SIMPLER flag is set to YES you can only use CheapSIMPLE as preconditioner in your "
            "fluid solver. Choose CheapSIMPLE in the SOLVER %i block in your dat file.",
            linsolvernumber);
      }

      solver->Params().sublist("CheapSIMPLE Parameters").set("Prec Type", "CheapSIMPLE");
      solver->Params().set(
          "ELCH", true);  // internal CheapSIMPLE modus for ML null space computation

      // add Inverse1 block for velocity dofs
      // tell Inverse1 block about NodalBlockInformation
      // In contrary to contact/meshtying problems this is necessary here, since we originally have
      // built the null space for the whole problem (velocity and pressure dofs). However, if we
      // split the matrix into velocity and pressure block, we have to adapt the null space
      // information for the subblocks. Therefore we need the nodal block information in the first
      // subblock for the velocities. The pressure null space is trivial to be built using a
      // constant vector
      auto& inv1 = solver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1");
      inv1.sublist("NodalBlockInformation") = solver->Params().sublist("NodalBlockInformation");
    }
  }

  SetIsSetup(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::ScaTraBaseAlgorithm::CreateScaTraFieldTest()
{
  return scatra_->CreateScaTraFieldTest();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraBaseAlgorithm::CheckIsSetup() const
{
  if (not IsSetup()) dserror("Setup() was not called.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraBaseAlgorithm::CheckIsInit() const
{
  if (not IsInit()) dserror("Init(...) was not called.");
}
