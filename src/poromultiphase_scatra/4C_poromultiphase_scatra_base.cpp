/*----------------------------------------------------------------------*/
/*! \file
 \brief base algorithm for scalar transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/

#include "4C_poromultiphase_scatra_base.hpp"

#include "4C_adapter_art_net.hpp"
#include "4C_adapter_porofluidmultiphase_wrapper.hpp"
#include "4C_adapter_poromultiphase.hpp"
#include "4C_adapter_scatra_base_algorithm.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_fluidporo_multiphase.hpp"
#include "4C_mat_scatra_multiporo.hpp"
#include "4C_poromultiphase_scatra_utils.hpp"
#include "4C_poromultiphase_utils.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_timint_meshtying_strategy_artery.hpp"
#include "4C_scatra_timint_poromulti.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::PoroMultiPhaseScaTraBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : AlgorithmBase(comm, globaltimeparams),
      poromulti_(Teuchos::null),
      scatra_(Teuchos::null),
      fluxreconmethod_(Inpar::POROFLUIDMULTIPHASE::gradreco_none),
      ndsporofluid_scatra_(-1),
      timertimestep_("PoroMultiPhaseScaTraBase", true),
      dttimestep_(0.0),
      divcontype_(Core::UTILS::IntegralValue<Inpar::PoroMultiPhaseScaTra::DivContAct>(
          globaltimeparams, "DIVERCONT")),
      artery_coupl_(Core::UTILS::IntegralValue<int>(globaltimeparams, "ARTERY_COUPLING"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::init(
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& algoparams,
    const Teuchos::ParameterList& poroparams, const Teuchos::ParameterList& structparams,
    const Teuchos::ParameterList& fluidparams, const Teuchos::ParameterList& scatraparams,
    const std::string& struct_disname, const std::string& fluid_disname,
    const std::string& scatra_disname, bool isale, int nds_disp, int nds_vel, int nds_solidpressure,
    int ndsporofluid_scatra, const std::map<int, std::set<int>>* nearbyelepairs)
{
  // save the dofset number of the scatra on the fluid dis
  ndsporofluid_scatra_ = ndsporofluid_scatra;

  // access the global problem
  Global::Problem* problem = Global::Problem::Instance();

  // Create the two uncoupled subproblems.

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // coupling scheme
  // -------------------------------------------------------------------
  // first of all check for possible couplings
  Inpar::POROMULTIPHASE::SolutionSchemeOverFields solschemeporo =
      Core::UTILS::IntegralValue<Inpar::POROMULTIPHASE::SolutionSchemeOverFields>(
          poroparams, "COUPALGO");
  Inpar::PoroMultiPhaseScaTra::SolutionSchemeOverFields solschemescatraporo =
      Core::UTILS::IntegralValue<Inpar::PoroMultiPhaseScaTra::SolutionSchemeOverFields>(
          algoparams, "COUPALGO");

  // partitioned -- monolithic not possible --> error
  if (solschemeporo !=
          Inpar::POROMULTIPHASE::SolutionSchemeOverFields::solscheme_twoway_monolithic &&
      solschemescatraporo ==
          Inpar::PoroMultiPhaseScaTra::SolutionSchemeOverFields::solscheme_twoway_monolithic)
    FOUR_C_THROW(
        "Your requested coupling is not available: possible couplings are:\n"
        "(STRUCTURE <--> FLUID) <--> SCATRA: partitioned -- partitioned_nested\n"
        "                                    monolithic  -- partitioned_nested\n"
        "                                    monolithic  -- monolithic\n"
        "YOUR CHOICE                       : partitioned -- monolithic");

  // monolithic -- partitioned sequential not possible
  if (solschemeporo ==
          Inpar::POROMULTIPHASE::SolutionSchemeOverFields::solscheme_twoway_monolithic &&
      solschemescatraporo == Inpar::PoroMultiPhaseScaTra::SolutionSchemeOverFields::
                                 solscheme_twoway_partitioned_sequential)
    FOUR_C_THROW(
        "Your requested coupling is not available: possible couplings are:\n"
        "(STRUCTURE <--> FLUID) <--> SCATRA: partitioned -- partitioned_nested\n"
        "                                    monolithic  -- partitioned_nested\n"
        "                                    monolithic  -- monolithic\n"
        "YOUR CHOICE                       : monolithic  -- partitioned_sequential");

  fluxreconmethod_ =
      Core::UTILS::IntegralValue<Inpar::POROFLUIDMULTIPHASE::FluxReconstructionMethod>(
          fluidparams, "FLUX_PROJ_METHOD");

  if (solschemescatraporo ==
          Inpar::PoroMultiPhaseScaTra::SolutionSchemeOverFields::solscheme_twoway_monolithic &&
      fluxreconmethod_ == Inpar::POROFLUIDMULTIPHASE::FluxReconstructionMethod::gradreco_l2)
  {
    FOUR_C_THROW(
        "Monolithic porofluidmultiphase-scatra coupling does not work with L2-projection!\n"
        "Set FLUX_PROJ_METHOD to none if you want to use monolithic coupling or use partitioned "
        "approach instead.");
  }

  poromulti_ =
      POROMULTIPHASE::UTILS::CreatePoroMultiPhaseAlgorithm(solschemeporo, globaltimeparams, Comm());

  // initialize
  poromulti_->init(globaltimeparams, poroparams, structparams, fluidparams, struct_disname,
      fluid_disname, isale, nds_disp, nds_vel, nds_solidpressure, ndsporofluid_scatra,
      nearbyelepairs);

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = scatraparams.get<int>("LINEAR_SOLVER");

  // scatra problem
  scatra_ = Teuchos::rcp(new Adapter::ScaTraBaseAlgorithm(globaltimeparams, scatraparams,
      problem->SolverParams(linsolvernumber), scatra_disname, true));

  // initialize the base algo.
  // scatra time integrator is constructed and initialized inside.
  scatra_->init();
  scatra_->ScaTraField()->set_number_of_dof_set_displacement(1);
  scatra_->ScaTraField()->set_number_of_dof_set_velocity(1);
  scatra_->ScaTraField()->set_number_of_dof_set_pressure(2);

  // do we perform coupling with 1D artery
  if (artery_coupl_)
  {
    // get mesh tying strategy
    scatramsht_ = Teuchos::rcp_dynamic_cast<ScaTra::MeshtyingStrategyArtery>(
        scatra_->ScaTraField()->Strategy());
    if (scatramsht_ == Teuchos::null) FOUR_C_THROW("cast to Meshtying strategy failed!");

    scatramsht_->set_artery_time_integrator(poro_field()->fluid_field()->ArtNetTimInt());
    scatramsht_->SetNearbyElePairs(nearbyelepairs);
  }

  // only now we must call setup() on the scatra time integrator.
  // all objects relying on the parallel distribution are
  // created and pointers are set.
  // calls setup() on the scatra time integrator inside.
  scatra_->ScaTraField()->setup();

  // do we perform coupling with 1D artery
  if (artery_coupl_)
  {
    // this check can only be performed after calling setup
    scatramsht_->CheckInitialFields();
  }

  std::vector<int> mydirichdofs(0);
  add_dirichmaps_volfrac_spec_ = Teuchos::rcp(new Epetra_Map(
      -1, 0, mydirichdofs.data(), 0, ScatraAlgo()->ScaTraField()->discretization()->Comm()));

  // done.
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::read_restart(int restart)
{
  if (restart)
  {
    // read restart data for structure field (will set time and step internally)
    poromulti_->read_restart(restart);

    // read restart data for scatra field (will set time and step internally)
    scatra_->ScaTraField()->read_restart(restart);
    if (artery_coupl_) scatramsht_->ArtScatraField()->read_restart(restart);

    // reset time and step for the global algorithm
    SetTimeStep(scatra_->ScaTraField()->Time(), restart);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::Timeloop()
{
  prepare_time_loop();

  while (NotFinished())
  {
    prepare_time_step();

    // reset timer
    timertimestep_.reset();
    // *********** time measurement ***********
    double dtcpu = timertimestep_.wallTime();
    // *********** time measurement ***********
    TimeStep();
    // *********** time measurement ***********
    double mydttimestep = timertimestep_.wallTime() - dtcpu;
    Comm().MaxAll(&mydttimestep, &dttimestep_, 1);
    // *********** time measurement ***********

    update_and_output();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::prepare_time_step(bool printheader)
{
  // the global control routine has its own time_ and step_ variables, as well as the single fields
  // keep them in sync!
  increment_time_and_step();

  if (printheader) print_header();

  SetPoroSolution();
  scatra_->ScaTraField()->prepare_time_step();
  if (artery_coupl_) scatramsht_->ArtScatraField()->prepare_time_step();
  // set structure-based scalar transport values
  SetScatraSolution();

  poromulti_->prepare_time_step();
  SetPoroSolution();
  apply_additional_dbc_for_vol_frac_species();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::prepare_time_loop()
{
  // set structure-based scalar transport values
  SetScatraSolution();
  poromulti_->prepare_time_loop();
  // initial output for scatra field
  SetPoroSolution();
  if (scatra_->ScaTraField()->HasExternalForce()) scatra_->ScaTraField()->SetExternalForce();
  scatra_->ScaTraField()->check_and_write_output_and_restart();
  if (artery_coupl_) scatramsht_->ArtScatraField()->check_and_write_output_and_restart();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::update_and_output()
{
  // set scatra on fluid (necessary for possible domain integrals)
  SetScatraSolution();
  poromulti_->update_and_output();

  // scatra field
  scatra_->ScaTraField()->update();
  scatra_->ScaTraField()->evaluate_error_compared_to_analytical_sol();
  scatra_->ScaTraField()->check_and_write_output_and_restart();
  // artery scatra field
  if (artery_coupl_)
  {
    scatramsht_->ArtScatraField()->update();
    scatramsht_->ArtScatraField()->evaluate_error_compared_to_analytical_sol();
    scatramsht_->ArtScatraField()->check_and_write_output_and_restart();
  }
  if (Comm().MyPID() == 0)
  {
    std::cout << "Finished POROMULTIPHASESCATRA STEP " << std::setw(5) << std::setprecision(4)
              << std::scientific << Step() << "/" << std::setw(5) << std::setprecision(4)
              << std::scientific << n_step() << ": dtstep = " << dttimestep_ << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::CreateFieldTest()
{
  Global::Problem* problem = Global::Problem::Instance();

  poromulti_->CreateFieldTest();
  problem->AddFieldTest(scatra_->create_sca_tra_field_test());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::SetPoroSolution()
{
  // safety check
  Teuchos::RCP<ScaTra::ScaTraTimIntPoroMulti> poroscatra =
      Teuchos::rcp_dynamic_cast<ScaTra::ScaTraTimIntPoroMulti>(scatra_->ScaTraField());
  if (poroscatra == Teuchos::null) FOUR_C_THROW("cast to ScaTraTimIntPoroMulti failed!");

  // set displacements
  poroscatra->ApplyMeshMovement(poromulti_->StructDispnp());

  // set the fluid solution
  poroscatra->set_solution_field_of_multi_fluid(
      poromulti_->RelaxedFluidPhinp(), poromulti_->fluid_field()->Phin());

  // additionally, set nodal flux if L2-projection is desired
  if (fluxreconmethod_ == Inpar::POROFLUIDMULTIPHASE::FluxReconstructionMethod::gradreco_l2)
    poroscatra->set_l2_flux_of_multi_fluid(poromulti_->FluidFlux());

  if (artery_coupl_)
  {
    scatramsht_->SetArteryPressure();
    scatramsht_->ApplyMeshMovement();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::apply_additional_dbc_for_vol_frac_species()
{
  // remove the old one
  ScatraAlgo()->ScaTraField()->remove_dirich_cond(add_dirichmaps_volfrac_spec_);

  std::vector<int> mydirichdofs(0);

  // get map and validdof-vector
  const Epetra_Map* elecolmap = ScatraAlgo()->ScaTraField()->discretization()->ElementColMap();
  Teuchos::RCP<const Epetra_Vector> valid_volfracspec_dofs =
      poro_field()->fluid_field()->valid_vol_frac_spec_dofs();

  // we identify the volume fraction species dofs which do not have a physical meaning and set a
  // DBC on them
  for (int iele = 0; iele < elecolmap->NumMyElements(); ++iele)
  {
    // dynamic_cast necessary because virtual inheritance needs runtime information
    Discret::ELEMENTS::Transport* myele = dynamic_cast<Discret::ELEMENTS::Transport*>(
        ScatraAlgo()->ScaTraField()->discretization()->gElement(elecolmap->GID(iele)));

    const Core::Mat::Material& material2 = *(myele->Material(2));

    // check the material
    if (material2.MaterialType() != Core::Materials::m_fluidporo_multiphase and
        material2.MaterialType() != Core::Materials::m_fluidporo_multiphase_reactions)
      FOUR_C_THROW("only poro multiphase and poro multiphase reactions material valid");

    // cast fluid material
    const Mat::FluidPoroMultiPhase& multiphasemat =
        static_cast<const Mat::FluidPoroMultiPhase&>(material2);

    const int numfluidphases = multiphasemat.NumFluidPhases();
    const int numvolfrac = multiphasemat.NumVolFrac();
    const int numfluidmat = multiphasemat.NumMat();

    // this is only necessary if we have volume fractions present
    // TODO: this works only if we have the same number of phases in every element
    if (numfluidmat == numfluidphases) return;

    const Core::Mat::Material& material = *(myele->Material());

    // cast scatra material
    const Mat::MatList& scatramat = static_cast<const Mat::MatList&>(material);

    if (not(scatramat.MaterialType() == Core::Materials::m_matlist or
            scatramat.MaterialType() == Core::Materials::m_matlist_reactions))
      FOUR_C_THROW("wrong type of ScaTra-Material");

    const int numscatramat = scatramat.NumMat();

    Core::Nodes::Node** nodes = myele->Nodes();

    for (int inode = 0; inode < (myele->num_node()); inode++)
    {
      if (nodes[inode]->Owner() == ScatraAlgo()->ScaTraField()->discretization()->Comm().MyPID())
      {
        std::vector<int> scatradofs =
            ScatraAlgo()->ScaTraField()->discretization()->Dof(0, nodes[inode]);
        std::vector<int> fluiddofs =
            poro_field()->fluid_field()->discretization()->Dof(0, nodes[inode]);

        for (int idof = 0; idof < numscatramat; ++idof)
        {
          int matid = scatramat.MatID(idof);
          Teuchos::RCP<Core::Mat::Material> singlemat = scatramat.MaterialById(matid);
          if (singlemat->MaterialType() == Core::Materials::m_scatra_multiporo_fluid ||
              singlemat->MaterialType() == Core::Materials::m_scatra_multiporo_solid ||
              singlemat->MaterialType() == Core::Materials::m_scatra_multiporo_temperature)
          {
            // do nothing
          }
          else if (singlemat->MaterialType() == Core::Materials::m_scatra_multiporo_volfrac)
          {
            const Teuchos::RCP<const Mat::ScatraMatMultiPoroVolFrac>& scatravolfracmat =
                Teuchos::rcp_dynamic_cast<const Mat::ScatraMatMultiPoroVolFrac>(singlemat);

            const int scalartophaseid = scatravolfracmat->PhaseID();
            // if not already in original dirich map     &&   if it is not a valid volume fraction
            // species dof identified with < 1
            if (ScatraAlgo()->ScaTraField()->DirichMaps()->CondMap()->LID(scatradofs[idof]) == -1 &&
                (int)(*valid_volfracspec_dofs)
                        [poro_field()->fluid_field()->discretization()->dof_row_map()->LID(
                            fluiddofs[scalartophaseid + numvolfrac])] < 1)
            {
              mydirichdofs.push_back(scatradofs[idof]);
              ScatraAlgo()->ScaTraField()->Phinp()->ReplaceGlobalValue(scatradofs[idof], 0, 0.0);
            }
          }
          else
            FOUR_C_THROW("only MAT_scatra_multiporo_(fluid,volfrac,solid,temperature) valid here");
        }
      }
    }
  }

  // build map
  int nummydirichvals = mydirichdofs.size();
  add_dirichmaps_volfrac_spec_ = Teuchos::rcp(new Epetra_Map(-1, nummydirichvals,
      mydirichdofs.data(), 0, ScatraAlgo()->ScaTraField()->discretization()->Comm()));

  // add the condition
  ScatraAlgo()->ScaTraField()->add_dirich_cond(add_dirichmaps_volfrac_spec_);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::SetScatraSolution()
{
  poromulti_->SetScatraSolution(ndsporofluid_scatra_, scatra_->ScaTraField()->Phinp());
  if (artery_coupl_)
    poromulti_->fluid_field()->ArtNetTimInt()->discretization()->set_state(
        2, "one_d_artery_phinp", scatramsht_->ArtScatraField()->Phinp());
}

/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::ScatraDofRowMap()
    const
{
  return scatra_->ScaTraField()->dof_row_map();
}

/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
void PoroMultiPhaseScaTra::PoroMultiPhaseScaTraBase::HandleDivergence() const
{
  switch (divcontype_)
  {
    case Inpar::PoroMultiPhaseScaTra::divcont_continue:
    {
      // warn if itemax is reached without convergence, but proceed to
      // next timestep...
      if (Comm().MyPID() == 0)
      {
        std::cout << "+---------------------------------------------------------------+"
                  << std::endl;
        std::cout << "|            >>>>>> continuing to next time step!               |"
                  << std::endl;
        std::cout << "+---------------------------------------------------------------+"
                  << std::endl
                  << std::endl;
      }
      break;
    }
    case Inpar::PoroMultiPhaseScaTra::divcont_stop:
    {
      FOUR_C_THROW("POROMULTIPHASESCATRA nonlinear solver not converged in ITEMAX steps!");
      break;
    }
    default:
      FOUR_C_THROW("unknown divercont action!");
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
