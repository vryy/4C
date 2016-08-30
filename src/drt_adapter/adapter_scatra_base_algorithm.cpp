/*----------------------------------------------------------------------*/
/*!
\file adapter_scatra_base_algorithm.cpp

\brief scalar transport field base algorithm

<pre>
\level 1

\maintainer Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
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

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_scatra/scatra_resulttest.H"

// general time integration schemes
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_stat.H"
#include "../drt_scatra/scatra_timint_ost.H"
#include "../drt_scatra/scatra_timint_bdf2.H"
#include "../drt_scatra/scatra_timint_genalpha.H"

// HDG time integration schemes
#include "../drt_scatra/scatra_timint_cardiac_monodomain_scheme_hdg.H"
#include "../drt_scatra/scatra_resulttest_hdg.H"

// loma specific files
#include "../drt_scatra/scatra_timint_loma_genalpha.H"
#include "../drt_scatra/scatra_timint_loma_ost.H"
#include "../drt_scatra/scatra_timint_loma_bdf2.H"

// elch specific files
#include "../drt_scatra/scatra_timint_elch_scheme.H"

// level set specific files
#include "../drt_levelset/levelset_timint_ost.H"
#include "../drt_levelset/levelset_timint_stat.H"

// cardiac monodomain specific files
#include "../drt_scatra/scatra_timint_cardiac_monodomain_scheme.H"

// cell migration specific files
#include "../drt_immersed_problem/scatra_timint_ost_endoexocytosis.H"

#include "../drt_scatra/scatra_timint_poromulti.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::ScaTraBaseAlgorithm::ScaTraBaseAlgorithm(
    const Teuchos::ParameterList&   prbdyn,         ///< parameter list for global problem
    const Teuchos::ParameterList&   scatradyn,      ///< parameter list for scalar transport subproblem
    const Teuchos::ParameterList&   solverparams,   ///< parameter list for scalar transport solver
    const std::string&              disname,        ///< name of scalar transport discretization
    const bool                      isale           ///< ALE flag
    )
{
  // setup scalar transport algorithm (overriding some dynamic parameters
  // with values specified in given problem-dependent ParameterList prbdyn)

  // -------------------------------------------------------------------
  // what's the current problem type?
  // -------------------------------------------------------------------
  PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis(disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // TODO: TAW use of solverparams??? change input parameter to solver number instead of parameter list? -> no default paramter possible any more
  Teuchos::RCP<LINALG::Solver> solver =
    Teuchos::rcp(new LINALG::Solver(solverparams,
                           actdis->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  // make a copy (inside an Teuchos::rcp) containing also all sublists
  Teuchos::RCP<Teuchos::ParameterList> scatratimeparams = Teuchos::rcp(new Teuchos::ParameterList(scatradyn));

  // -------------------------------------------------------------------
  // overrule certain parameters for coupled problems
  // -------------------------------------------------------------------
  // the default time step size
  scatratimeparams->set<double>   ("TIMESTEP"    ,prbdyn.get<double>("TIMESTEP"));
  // maximum simulation time
  scatratimeparams->set<double>   ("MAXTIME"     ,prbdyn.get<double>("MAXTIME"));
  // maximum number of timesteps
  scatratimeparams->set<int>      ("NUMSTEP"     ,prbdyn.get<int>("NUMSTEP"));
  // restart
  scatratimeparams->set           ("RESTARTEVRY" ,prbdyn.get<int>("RESTARTEVRY"));
  // solution output
  scatratimeparams->set           ("RESULTSEVRY"       ,prbdyn.get<int>("RESULTSEVRY"));

  // -------------------------------------------------------------------
  // overrule flags for solid-based scalar transport!
  // (assumed disname = "scatra2" for solid-based scalar transport)
  // -------------------------------------------------------------------
  if (probtype == prb_ac_fsi or probtype == prb_biofilm_fsi or probtype == prb_gas_fsi or probtype ==prb_fps3i)
  {
    //scatra1 (=fluid scalar) get's inputs from SCALAR TRANSPORT DYNAMIC/STABILIZATION, hence nothing is to do here
//    if (disname== "scatra1") //get's inputs from SCALAR TRANSPORT DYNAMIC/STABILIZATION

    if (disname == "scatra2") //structure_scatra discretisation
    {
      //scatra2 (=structure scalar) get's inputs from FS3I DYNAMIC/STRUCTURE SCALAR STABILIZATION, hence we have to replace it
      scatratimeparams->sublist("STABILIZATION")=prbdyn.sublist("STRUCTURE SCALAR STABILIZATION");
      scatratimeparams->set<std::string>("CONVFORM",prbdyn.get<std::string>("STRUCTSCAL_CONVFORM"));

      //scatra2 get's in initial functions from FS3I DYNAMICS
      switch ( DRT::INPUT::IntegralValue<INPAR::SCATRA::InitialField>(prbdyn,"STRUCTSCAL_INITIALFIELD") )
      {
      case INPAR::SCATRA::initfield_zero_field:
        scatratimeparams->set<std::string>("INITIALFIELD","zero_field"); //we want zero initial conditions for the structure scalar
        scatratimeparams->set<int>("INITFUNCNO",-1);
        break;
      case INPAR::SCATRA::initfield_field_by_function:
        scatratimeparams->set<std::string>("INITIALFIELD","field_by_function"); //we want the same initial conditions for structure scalar as for the fluid scalar
        scatratimeparams->set<int>("INITFUNCNO", prbdyn.get<int>("STRUCTSCAL_INITFUNCNO") );
        break;
      default:
        dserror("Your STRUCTSCAL_INITIALFIELD type is not supported!");
        break;
      }
    }
  }

  // -------------------------------------------------------------------
  // list for extra parameters
  // (put here everything that is not available in scatradyn or its sublists)
  // -------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> extraparams
    = Teuchos::rcp(new Teuchos::ParameterList());

  // ------------------------------pointer to the error file (for output)
  extraparams->set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

  // ----------------Eulerian or ALE formulation of transport equation(s)
  extraparams->set<bool>("isale",isale);

  // ------------------------------------get also fluid turbulence sublist
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  extraparams->sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");
  extraparams->sublist("SUBGRID VISCOSITY")=fdyn.sublist("SUBGRID VISCOSITY");
  extraparams->sublist("MULTIFRACTAL SUBGRID SCALES")=fdyn.sublist("MULTIFRACTAL SUBGRID SCALES");
  extraparams->sublist("TURBULENT INFLOW")=fdyn.sublist("TURBULENT INFLOW");

  // -------------------------------------------------------------------
  // algorithm construction depending on problem type and
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::SCATRA::TimeIntegrationScheme timintscheme =
    DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(scatradyn,"TIMEINTEGR");

  // low Mach number flow
  if (probtype == prb_loma)
  {
    Teuchos::RCP<Teuchos::ParameterList> lomaparams = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->LOMAControlParams()));
    switch(timintscheme)
    {
      case INPAR::SCATRA::timeint_gen_alpha:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::TimIntLomaGenAlpha(actdis, solver, lomaparams, scatratimeparams,extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_one_step_theta:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::TimIntLomaOST(actdis, solver, lomaparams, scatratimeparams,extraparams, output));
        break;
      }
//      case INPAR::SCATRA::timeint_bdf2:
//      {
//        // create instance of time integration class (call the constructor)
//        scatra_ = Teuchos::rcp(new SCATRA::TimIntLomaBDF2(actdis, solver, lomaparams, scatratimeparams,extraparams, output));
//        break;
//      }
      default:
        dserror("Unknown time integration scheme for loMa!");
        break;
    }
  }

  // electrochemistry
  else if (probtype == prb_elch or (probtype == prb_sti and disname == "scatra" and DRT::INPUT::IntegralValue<INPAR::SCATRA::ImplType>(DRT::Problem::Instance()->STIDynamicParams(),"SCATRATYPE") == INPAR::SCATRA::impltype_elch_diffcond))
  {
    Teuchos::RCP<Teuchos::ParameterList> elchparams = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->ELCHControlParams()));

    // create a 2nd solver for block-preconditioning if chosen from input
    if (DRT::INPUT::IntegralValue<int>(*elchparams,"BLOCKPRECOND"))
    {
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      INPAR::SOLVER::AzPrecType prec = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(DRT::Problem::Instance()->SolverParams(linsolvernumber),"AZPREC");
      if (prec != INPAR::SOLVER::azprec_CheapSIMPLE &&
          prec != INPAR::SOLVER::azprec_TekoSIMPLE)  // TODO adapt error message
        dserror("If SIMPLER flag is set to YES you can only use CheapSIMPLE or TekoSIMPLE as preconditioners in your fluid solver. Choose CheapSIMPLE or TekoSIMPLE in the SOLVER %i block in your dat file.",linsolvernumber);

      solver->Params().sublist("CheapSIMPLE Parameters").set("Prec Type","CheapSIMPLE");
      solver->Params().set("ELCH",true); // internal CheapSIMPLE modus for ML null space computation

      // add Inverse1 block for velocity dofs
      // tell Inverse1 block about NodalBlockInformation
      // In contrary to contact/meshtying problems this is necessary here, since we originally have built the
      // null space for the whole problem (velocity and pressure dofs). However, if we split the matrix into
      // velocity and pressure block, we have to adapt the null space information for the subblocks. Therefore
      // we need the nodal block information in the first subblock for the velocities. The pressure null space
      // is trivial to be built using a constant vector
      Teuchos::ParameterList& inv1 = solver->Params().sublist("CheapSIMPLE Parameters").sublist("Inverse1");
      inv1.sublist("NodalBlockInformation") = solver->Params().sublist("NodalBlockInformation");
    }

    switch(timintscheme)
    {
      case INPAR::SCATRA::timeint_one_step_theta:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntElchOST(actdis, solver, elchparams, scatratimeparams, extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_bdf2:
      {
         // create instance of time integration class (call the constructor)
         scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntElchBDF2(actdis, solver, elchparams, scatratimeparams, extraparams, output));
         break;
      }
      case INPAR::SCATRA::timeint_gen_alpha:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntElchGenAlpha(actdis, solver, elchparams, scatratimeparams, extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_stationary:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntElchStationary(actdis, solver, elchparams,  scatratimeparams, extraparams, output));
        break;
      }
      default:
        dserror("Unknown time integration scheme for electrochemistry!");
        break;
    }
  }

  // levelset, combustion, and two phase flow
  else if (probtype == prb_level_set or probtype == prb_combust or probtype == prb_two_phase_flow or probtype == prb_fluid_xfem_ls)
  {
    Teuchos::RCP<Teuchos::ParameterList> lsparams = Teuchos::null;
    if (probtype == prb_level_set)
      lsparams = Teuchos::rcp(new Teuchos::ParameterList(prbdyn));
    else
    {
      lsparams = Teuchos::rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->LevelSetControl()));
      // overrule certain parameters for coupled problems
      // this has already been ensured for scatratimeparams, but has also been ensured for the level-set
      // parameter which are potentially handled down to the particle algorithm in a hybrid approach
      // time step size
      lsparams->set<double>   ("TIMESTEP"    ,prbdyn.get<double>("TIMESTEP"));
      // maximum simulation time
      lsparams->set<double>   ("MAXTIME"     ,prbdyn.get<double>("MAXTIME"));
      // maximum number of timesteps
      lsparams->set<int>      ("NUMSTEP"     ,prbdyn.get<int>("NUMSTEP"));
      // restart
      lsparams->set           ("RESTARTEVRY" ,prbdyn.get<int>("RESTARTEVRY"));
      // solution output
      lsparams->set           ("RESULTSEVRY"       ,prbdyn.get<int>("RESULTSEVRY"));

      if (probtype == prb_two_phase_flow)
      {
        //Give access to smoothing parameter for levelset calculations.
        lsparams->set<double> ("INTERFACE_THICKNESS_TPF",prbdyn.sublist("SMEARED").get<double>("INTERFACE_THICKNESS"));
      }
    }

    switch(timintscheme)
    {
      case INPAR::SCATRA::timeint_one_step_theta:
      {
        // create instance of time integration class (call the constructor)
       scatra_ = Teuchos::rcp(new SCATRA::LevelSetTimIntOneStepTheta(actdis, solver, lsparams, scatratimeparams, extraparams,output));
        break;
      }
      case INPAR::SCATRA::timeint_stationary:
      {
        // create instance of time integration class (call the constructor)
        if (probtype != prb_level_set)
          scatra_ = Teuchos::rcp(new SCATRA::LevelSetTimIntStationary(actdis, solver, lsparams, scatratimeparams,extraparams, output));
        else
          dserror("Stationary time integration scheme only supported for coupled level-set problems!");
        break;
      }
      case INPAR::SCATRA::timeint_gen_alpha:
      {
        if (probtype == prb_combust)
        {
          std::cout << "\n\n\n WARNING: Level set algorithm does not yet support gen-alpha. You thus get a standard Scatra!\n\n\n" << std::endl;
          // create instance of time integration class (call the constructor)
          scatra_ = Teuchos::rcp(new SCATRA::TimIntGenAlpha(actdis, solver, scatratimeparams,extraparams, output));
        }
        else if (probtype == prb_two_phase_flow)
        {
          std::cout << "\n\n\n WARNING: Level set algorithm does not yet support gen-alpha. You thus get a standard Scatra!\n\n\n" << std::endl;
          // create instance of time integration class (call the constructor)
          scatra_ = Teuchos::rcp(new SCATRA::TimIntGenAlpha(actdis, solver, scatratimeparams,extraparams, output));
        }
        else
          dserror("Unknown time-integration scheme for level-set problem");

        break;
      }
      default:
        dserror("Unknown time-integration scheme for level-set problem");
        break;
    }// switch(timintscheme)
  }

  // cardiac monodomain
  else if (probtype == prb_cardiac_monodomain or (probtype == prb_ssi and DRT::INPUT::IntegralValue<INPAR::SCATRA::ImplType>(DRT::Problem::Instance()->SSIControlParams(),"SCATRATYPE") == INPAR::SCATRA::impltype_cardiac_monodomain))
  {
    Teuchos::RCP<Teuchos::ParameterList> cmonoparams = rcp(new Teuchos::ParameterList(DRT::Problem::Instance()->EPControlParams()));

    // HDG implements all time stepping schemes within gen-alpha
    if (DRT::Problem::Instance()->SpatialApproximation() == "HDG")
      scatra_ = Teuchos::rcp(new SCATRA::TimIntCardiacMonodomainHDG(actdis, solver, cmonoparams, scatratimeparams, extraparams, output));
    else
    {
      switch(timintscheme)
      {
        case INPAR::SCATRA::timeint_gen_alpha:
        {
          // create instance of time integration class (call the constructor)
          scatra_ = Teuchos::rcp(new SCATRA::TimIntCardiacMonodomainGenAlpha(actdis, solver, cmonoparams, scatratimeparams,extraparams, output));
          break;
        }
        case INPAR::SCATRA::timeint_one_step_theta:
        {
          // create instance of time integration class (call the constructor)
          scatra_ = Teuchos::rcp(new SCATRA::TimIntCardiacMonodomainOST(actdis, solver, cmonoparams, scatratimeparams,extraparams, output));
          break;
        }
        case INPAR::SCATRA::timeint_bdf2:
      {
          // create instance of time integration class (call the constructor)
          scatra_ = Teuchos::rcp(new SCATRA::TimIntCardiacMonodomainBDF2(actdis, solver, cmonoparams, scatratimeparams,extraparams, output));
          break;
        }
        default:
          dserror("Unknown time integration scheme for cardiac monodomain problem!");
          break;
      }// switch(timintscheme)
    }
  }
  // cell migration endo-/exocytosis
  else if (probtype == prb_scatra_endoexocytosis)
  {
    switch(timintscheme)
    {
    case INPAR::SCATRA::timeint_one_step_theta:
    {
      // create instance of time integration class (call the constructor)
      scatra_ = Teuchos::rcp(new SCATRA::TimIntOneStepThetaEndoExocytosis(actdis,solver,scatratimeparams,extraparams,output,0));
      break;
    }
    default:
      dserror("Unknown time integration scheme for cell migration problem! Use One-Step-Theta Scheme!");
      break;
    }
  }
  else if (probtype == prb_poromultiphasescatra)
  {

    switch(timintscheme)
    {
      case INPAR::SCATRA::timeint_gen_alpha:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntPoroMultiGenAlpha(actdis, solver, Teuchos::null, scatratimeparams,extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_one_step_theta:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntPoroMultiOST(actdis, solver, Teuchos::null, scatratimeparams,extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_bdf2:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntPoroMultiBDF2(actdis, solver, Teuchos::null, scatratimeparams,extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_stationary:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::ScaTraTimIntPoroMultiStationary(actdis, solver, Teuchos::null, scatratimeparams,extraparams, output));
        break;
      }
      default:
        dserror("Unknown time integration scheme for porous medium multiphase problem!");
        break;
    }// switch(timintscheme)
  }
  // everything else
  else
  {
    // HDG implements all time stepping schemes within gen-alpha
    if (DRT::Problem::Instance()->SpatialApproximation() == "HDG")
      scatra_ = Teuchos::rcp(new SCATRA::TimIntHDG(actdis, solver, scatratimeparams, extraparams, output));
    else
    {
      switch(timintscheme)
      {
      case INPAR::SCATRA::timeint_stationary:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::TimIntStationary(actdis, solver, scatratimeparams, extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_one_step_theta:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::TimIntOneStepTheta(actdis, solver, scatratimeparams, extraparams,output));
        break;
      }
      case INPAR::SCATRA::timeint_bdf2:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::TimIntBDF2(actdis, solver, scatratimeparams,extraparams, output));
        break;
      }
      case INPAR::SCATRA::timeint_gen_alpha:
      {
        // create instance of time integration class (call the constructor)
        scatra_ = Teuchos::rcp(new SCATRA::TimIntGenAlpha(actdis, solver, scatratimeparams,extraparams, output));
        break;
      }
      default:
        dserror("Unknown time-integration scheme for scalar transport problem");
        break;
      } // switch(timintscheme)
    }
  }

  // initialize algorithm for specific time-integration scheme
  scatra_->Init();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> ADAPTER::ScaTraBaseAlgorithm::CreateScaTraFieldTest()
{
  if (DRT::Problem::Instance()->SpatialApproximation() == "HDG")
    return Teuchos::rcp(new SCATRA::HDGResultTest(scatra_));
  else
    return Teuchos::rcp(new SCATRA::ScaTraResultTest(scatra_));
}


