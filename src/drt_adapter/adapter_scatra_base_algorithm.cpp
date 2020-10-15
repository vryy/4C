/*----------------------------------------------------------------------*/
/*! \file

\brief scalar transport field base algorithm

\level 1


*/
/*----------------------------------------------------------------------*/

#include "adapter_scatra_base_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../linalg/linalg_solver.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_elch.H"
#include "../drt_inpar/inpar_cardiac_monodomain.H"
#include "../drt_inpar/inpar_ssi.H"
#include "../drt_inpar/inpar_ssti.H"
#include "../drt_inpar/inpar_sti.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// general time integration schemes
#include "../drt_scatra/scatra_timint_stat.H"
#include "../drt_scatra/scatra_timint_ost.H"
#include "../drt_scatra/scatra_timint_bdf2.H"
#include "../drt_scatra/scatra_timint_genalpha.H"

// HDG time integration schemes
#include "../drt_scatra/scatra_timint_cardiac_monodomain_scheme_hdg.H"
#include "../drt_scatra/scatra_timint_stat_hdg.H"
#include "../drt_scatra/scatra_resulttest_hdg.H"

// loma specific files
#include "../drt_scatra/scatra_timint_loma_genalpha.H"
#include "../drt_scatra/scatra_timint_loma_ost.H"
#include "../drt_scatra/scatra_timint_loma_bdf2.H"

// elch specific files
#include "../drt_scatra/scatra_resulttest_elch.H"
#include "../drt_scatra/scatra_timint_elch_scheme.H"

// level set specific files
#include "../drt_levelset/levelset_timint_ost.H"
#include "../drt_levelset/levelset_timint_stat.H"

// xcontact level set specific files
#include "../drt_contact_xcontact/xcontact_levelset_timint_ost.H"

// cardiac monodomain specific files
#include "../drt_scatra/scatra_timint_cardiac_monodomain_scheme.H"

// poro multiphase files
#include "../drt_scatra/scatra_timint_poromulti.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraBaseAlgorithm::ScaTraBaseAlgorithm()
    : scatra_(Teuchos::null), issetup_(false), isinit_(false)
{
  // Keep constructor empty
  return;
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
  ProblemType probtype = DRT::Problem::Instance()->GetProblemType();

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> discret = DRT::Problem::Instance()->GetDis(disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!discret->Filled() or !discret->HaveDofs()) discret->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = discret->Writer();
  if (discret->NumGlobalElements() == 0)
    dserror("No elements in discretization %s", discret->Name().c_str());
  output->WriteMesh(0, 0.0);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // TODO: TAW use of solverparams???
  // change input parameter to solver number instead of parameter list?
  // -> no default paramter possible any more
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(
      solverparams, discret->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  // make a copy (inside an Teuchos::rcp) containing also all sublists
  Teuchos::RCP<Teuchos::ParameterList> scatratimeparams =
      Teuchos::rcp(new Teuchos::ParameterList(scatradyn));
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
  if (probtype == prb_ac_fsi or probtype == prb_biofilm_fsi or probtype == prb_gas_fsi or
      probtype == prb_fps3i or probtype == prb_thermo_fsi)
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
  Teuchos::RCP<Teuchos::ParameterList> extraparams = Teuchos::rcp(new Teuchos::ParameterList());

  // ------------------------------pointer to the error file (for output)
  extraparams->set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());

  // ----------------Eulerian or ALE formulation of transport equation(s)
  extraparams->set<bool>("isale", isale);

  // ------------------------------------get also fluid turbulence sublist
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
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
  INPAR::SCATRA::TimeIntegrationScheme timintscheme =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(scatradyn, "TIMEINTEGR");

  // low Mach number flow
  if (probtype == prb_loma or probtype == prb_thermo_fsi)
  {
    Teuchos::RCP<Teuchos::ParameterList> lomaparams =
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
      case INPAR::SCATRA::timeint_one_step_theta:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::TimIntLomaOST(
            discret, solver, lomaparams, scatratimeparams, extraparams, output));
        break;
      }
      //      case INPAR::SCATRA::timeint_bdf2:
      //      {
      //        // create instance of time integration class (call the constructor)
      //        scatra_ = Teuchos::rcp(new SCATRA::TimIntLomaBDF2(discret, solver, lomaparams,
      //        scatratimeparams,extraparams, output)); break;
      //      }
      default:
        dserror("Unknown time integration scheme for loMa!");
        break;
    }
  }

  // electrochemistry
  else if (probtype == prb_elch or
           (probtype == prb_ssi and DRT::INPUT::IntegralValue<INPAR::SSI::ScaTraTimIntType>(
                                        DRT::Problem::Instance()->SSIControlParams(),
                                        "SCATRATIMINTTYPE") == INPAR::SSI::scatratiminttype_elch) or
           (probtype == prb_ssti and disname == "scatra" and
               Teuchos::getIntegralValue<INPAR::SSTI::ScaTraTimIntType>(
                   DRT::Problem::Instance()->SSTIControlParams(), "SCATRATIMINTTYPE") ==
                   INPAR::SSTI::ScaTraTimIntType::elch) or
           (probtype == prb_sti and disname == "scatra" and
               Teuchos::getIntegralValue<INPAR::STI::ScaTraTimIntType>(
                   DRT::Problem::Instance()->STIDynamicParams(), "SCATRATIMINTTYPE") ==
                   INPAR::STI::ScaTraTimIntType::elch))
  {
    Teuchos::RCP<Teuchos::ParameterList> elchparams =
        Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->ELCHControlParams()));

    switch (timintscheme)
    {
      case INPAR::SCATRA::timeint_one_step_theta:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntElchOST(
            discret, solver, elchparams, scatratimeparams, extraparams, output));
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
  else if (probtype == prb_level_set or probtype == prb_two_phase_flow or
           probtype == prb_fluid_xfem_ls or probtype == prb_xcontact)
  {
    Teuchos::RCP<Teuchos::ParameterList> lsparams = Teuchos::null;
    switch (probtype)
    {
      case prb_level_set:
        lsparams = Teuchos::rcp(new Teuchos::ParameterList(prbdyn));
        break;
      case prb_two_phase_flow:
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
        switch (probtype)
        {
          case prb_xcontact:
            // create instance of time integration class (call the constructor)
            scatra_ = Teuchos::rcp(new XCONTACT::LEVELSET::TIMINT::OneStepTheta(
                discret, solver, lsparams, scatratimeparams, extraparams, output));
            break;
          default:
          {
            // create instance of time integration class (call the constructor)
            scatra_ = Teuchos::rcp(new SCATRA::LevelSetTimIntOneStepTheta(
                discret, solver, lsparams, scatratimeparams, extraparams, output));
            break;
          }
        }
        break;
      }
      case INPAR::SCATRA::timeint_stationary:
      {
        // create instance of time integration class (call the constructor)
        switch (probtype)
        {
          case prb_xcontact:
          case prb_level_set:
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
          case prb_two_phase_flow:
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
  else if (probtype == prb_cardiac_monodomain or
           (probtype == prb_ssi and
               DRT::INPUT::IntegralValue<INPAR::SSI::ScaTraTimIntType>(
                   DRT::Problem::Instance()->SSIControlParams(), "SCATRATIMINTTYPE") ==
                   INPAR::SSI::scatratiminttype_cardiac_monodomain))
  {
    Teuchos::RCP<Teuchos::ParameterList> cmonoparams =
        rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->EPControlParams()));

    // HDG implements all time stepping schemes within gen-alpha
    if (DRT::Problem::Instance()->SpatialApproximationType() ==
        ShapeFunctionType::shapefunction_hdg)
      scatra_ = Teuchos::rcp(new SCATRA::TimIntCardiacMonodomainHDG(
          discret, solver, cmonoparams, scatratimeparams, extraparams, output));
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
  else if (probtype == prb_poromultiphasescatra)
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
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::ScaTraBaseAlgorithm::Setup()
{
  CheckIsInit();

  // setup the time integrator
  scatra_->Setup();

  // get the parameter list
  Teuchos::RCP<Teuchos::ParameterList> scatradyn = scatra_->ScatraParameterList();
  // get the discretization
  Teuchos::RCP<DRT::Discretization> discret = scatra_->Discretization();

  // -------------------------------------------------------------------
  // what's the current problem type?
  // -------------------------------------------------------------------
  ProblemType probtype = DRT::Problem::Instance()->GetProblemType();

  // prepare fixing the null space for electrochemistry and sti
  if (probtype == prb_elch or (probtype == prb_sti and discret->Name() == "scatra" and
                                  Teuchos::getIntegralValue<INPAR::STI::ScaTraTimIntType>(
                                      DRT::Problem::Instance()->STIDynamicParams(),
                                      "SCATRATIMINTTYPE") == INPAR::STI::ScaTraTimIntType::elch))
  {
    Teuchos::RCP<Teuchos::ParameterList> elchparams =
        Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->ELCHControlParams()));

    // create a 2nd solver for block-preconditioning if chosen from input
    if (DRT::INPUT::IntegralValue<int>(*elchparams, "BLOCKPRECOND"))
    {
      const Teuchos::RCP<LINALG::Solver>& solver = scatra_->Solver();

      const int linsolvernumber = scatradyn->get<int>("LINEAR_SOLVER");
      INPAR::SOLVER::AzPrecType prec = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(
          DRT::Problem::Instance()->SolverParams(linsolvernumber), "AZPREC");
      if (prec != INPAR::SOLVER::azprec_CheapSIMPLE &&
          prec != INPAR::SOLVER::azprec_TekoSIMPLE)  // TODO adapt error message
        dserror(
            "If SIMPLER flag is set to YES you can only use CheapSIMPLE or TekoSIMPLE as "
            "preconditioners in your fluid solver. Choose CheapSIMPLE or TekoSIMPLE in the SOLVER "
            "%i block in your dat file.",
            linsolvernumber);

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
      Teuchos::ParameterList& inv1 =
          solver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1");
      inv1.sublist("NodalBlockInformation") = solver->Params().sublist("NodalBlockInformation");
    }
  }

  SetIsSetup(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::ScaTraBaseAlgorithm::CreateScaTraFieldTest()
{
  if (DRT::Problem::Instance()->SpatialApproximationType() == ShapeFunctionType::shapefunction_hdg)
    return Teuchos::rcp(new SCATRA::HDGResultTest(scatra_));
  else if (DRT::Problem::Instance()->GetProblemType() == prb_elch or
           (DRT::Problem::Instance()->GetProblemType() == prb_ssi and
               DRT::INPUT::IntegralValue<INPAR::SSI::ScaTraTimIntType>(
                   DRT::Problem::Instance()->SSIControlParams(), "SCATRATIMINTTYPE") ==
                   INPAR::SSI::scatratiminttype_elch) or
           (DRT::Problem::Instance()->GetProblemType() == prb_ssti and
               Teuchos::getIntegralValue<INPAR::SSTI::ScaTraTimIntType>(
                   DRT::Problem::Instance()->SSTIControlParams(), "SCATRATIMINTTYPE") ==
                   INPAR::SSTI::ScaTraTimIntType::elch and
               scatra_->Discretization()->Name() == "scatra"))
  {
    return Teuchos::rcp(
        new SCATRA::ElchResultTest(Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntElch>(scatra_)));
  }
  else
    return Teuchos::rcp(new SCATRA::ScaTraResultTest(scatra_));
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
