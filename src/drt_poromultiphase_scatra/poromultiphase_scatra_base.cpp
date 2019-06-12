/*----------------------------------------------------------------------*/
/*!
 \brief base algorithm for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#include "poromultiphase_scatra_base.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_adapter/ad_poromultiphase.H"
#include "../drt_poromultiphase/poromultiphase_utils.H"
#include "poromultiphase_scatra_utils.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_poromulti.H"

#include "../drt_adapter/ad_art_net.H"
#include "../drt_adapter/ad_porofluidmultiphase_wrapper.H"

#include "../drt_scatra/scatra_timint_meshtying_strategy_artery.H"
#include "../drt_scatra_ele/scatra_ele.H"
#include "../drt_mat/fluidporo_multiphase.H"
#include "../drt_mat/scatra_mat_multiporo.H"

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16  |
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::PoroMultiPhaseScaTraBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : AlgorithmBase(comm, globaltimeparams),
      poromulti_(Teuchos::null),
      scatra_(Teuchos::null),
      fluxreconmethod_(INPAR::POROFLUIDMULTIPHASE::gradreco_none),
      ndsporofluid_scatra_(-1),
      timertimestep_(comm),
      dttimestep_(0.0),
      divcontype_(DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASESCATRA::DivContAct>(
          globaltimeparams, "DIVERCONT")),
      artery_coupl_(DRT::INPUT::IntegralValue<int>(globaltimeparams, "ARTERY_COUPLING"))
{
}


/*----------------------------------------------------------------------*
 | initialize algorithm                                    vuong 08/16  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::Init(
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
  DRT::Problem* problem = DRT::Problem::Instance();

  // Create the two uncoupled subproblems.

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // coupling scheme
  // -------------------------------------------------------------------
  // first of all check for possible couplings
  INPAR::POROMULTIPHASE::SolutionSchemeOverFields solschemeporo =
      DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASE::SolutionSchemeOverFields>(
          poroparams, "COUPALGO");
  INPAR::POROMULTIPHASESCATRA::SolutionSchemeOverFields solschemescatraporo =
      DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASESCATRA::SolutionSchemeOverFields>(
          algoparams, "COUPALGO");

  // partitioned -- monolithic not possible --> error
  if (solschemeporo !=
          INPAR::POROMULTIPHASE::SolutionSchemeOverFields::solscheme_twoway_monolithic &&
      solschemescatraporo ==
          INPAR::POROMULTIPHASESCATRA::SolutionSchemeOverFields::solscheme_twoway_monolithic)
    dserror(
        "Your requested coupling is not available: possible couplings are:\n"
        "(STRUCTURE <--> FLUID) <--> SCATRA: partitioned -- partitioned_nested\n"
        "                                    monolithic  -- partitioned_nested\n"
        "                                    monolithic  -- monolithic\n"
        "YOUR CHOICE                       : partitioned -- monolithic");

  // monolithic -- partitioned sequential not possible
  if (solschemeporo ==
          INPAR::POROMULTIPHASE::SolutionSchemeOverFields::solscheme_twoway_monolithic &&
      solschemescatraporo == INPAR::POROMULTIPHASESCATRA::SolutionSchemeOverFields::
                                 solscheme_twoway_partitioned_sequential)
    dserror(
        "Your requested coupling is not available: possible couplings are:\n"
        "(STRUCTURE <--> FLUID) <--> SCATRA: partitioned -- partitioned_nested\n"
        "                                    monolithic  -- partitioned_nested\n"
        "                                    monolithic  -- monolithic\n"
        "YOUR CHOICE                       : monolithic  -- partitioned_sequential");

  fluxreconmethod_ =
      DRT::INPUT::IntegralValue<INPAR::POROFLUIDMULTIPHASE::FluxReconstructionMethod>(
          fluidparams, "FLUX_PROJ_METHOD");

  if (solschemescatraporo ==
          INPAR::POROMULTIPHASESCATRA::SolutionSchemeOverFields::solscheme_twoway_monolithic &&
      fluxreconmethod_ == INPAR::POROFLUIDMULTIPHASE::FluxReconstructionMethod::gradreco_l2)
  {
    dserror(
        "Monolithic porofluidmultiphase-scatra coupling does not work with L2-projection!\n"
        "Set FLUX_PROJ_METHOD to none if you want to use monolithic coupling or use partitioned "
        "approach instead.");
  }

  poromulti_ =
      POROMULTIPHASE::UTILS::CreatePoroMultiPhaseAlgorithm(solschemeporo, globaltimeparams, Comm());

  // initialize
  poromulti_->Init(globaltimeparams, poroparams, structparams, fluidparams, struct_disname,
      fluid_disname, isale, nds_disp, nds_vel, nds_solidpressure, ndsporofluid_scatra,
      nearbyelepairs);

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = scatraparams.get<int>("LINEAR_SOLVER");

  // scatra problem
  scatra_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm());

  // initialize the base algo.
  // scatra time integrator is constructed and initialized inside.
  scatra_->Init(
      globaltimeparams, scatraparams, problem->SolverParams(linsolvernumber), scatra_disname, true);

  // do we perform coupling with 1D artery
  if (artery_coupl_)
  {
    // get mesh tying strategy
    scatramsht_ = Teuchos::rcp_dynamic_cast<SCATRA::MeshtyingStrategyArtery>(
        scatra_->ScaTraField()->Strategy());
    if (scatramsht_ == Teuchos::null) dserror("cast to Meshtying strategy failed!");

    scatramsht_->SetArteryTimeIntegrator(PoroField()->FluidField()->ArtNetTimInt());
    scatramsht_->SetNearbyElePairs(nearbyelepairs);
  }

  // only now we must call Setup() on the scatra time integrator.
  // all objects relying on the parallel distribution are
  // created and pointers are set.
  // calls Setup() on the scatra time integrator inside.
  scatra_->ScaTraField()->Setup();

  // do we perform coupling with 1D artery
  if (artery_coupl_)
  {
    // this check can only be performed after calling setup
    scatramsht_->CheckInitialFields();
  }

  std::vector<int> mydirichdofs(0);
  add_dirichmaps_volfrac_spec_ = Teuchos::rcp(new Epetra_Map(
      -1, 0, &(mydirichdofs[0]), 0, ScatraAlgo()->ScaTraField()->Discretization()->Comm()));

  // done.
  return;
}

/*----------------------------------------------------------------------*
 | read restart information for given time step (public)   vuong 08/16  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::ReadRestart(int restart)
{
  if (restart)
  {
    // read restart data for structure field (will set time and step internally)
    poromulti_->ReadRestart(restart);

    // read restart data for scatra field (will set time and step internally)
    scatra_->ScaTraField()->ReadRestart(restart);
    if (artery_coupl_) scatramsht_->ArtScatraField()->ReadRestart(restart);

    // reset time and step for the global algorithm
    SetTimeStep(scatra_->ScaTraField()->Time(), restart);
  }

  return;
}

/*----------------------------------------------------------------------*
 | time loop                                                 vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::Timeloop()
{
  PrepareTimeLoop();

  while (NotFinished())
  {
    PrepareTimeStep();

    // reset timer
    timertimestep_.ResetStartTime();
    // *********** time measurement ***********
    double dtcpu = timertimestep_.WallTime();
    // *********** time measurement ***********
    TimeStep();
    // *********** time measurement ***********
    double mydttimestep = timertimestep_.WallTime() - dtcpu;
    Comm().MaxAll(&mydttimestep, &dttimestep_, 1);
    // *********** time measurement ***********

    UpdateAndOutput();
  }
  return;
}

/*----------------------------------------------------------------------*
 | prepare one time step                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::PrepareTimeStep(bool printheader)
{
  // the global control routine has its own time_ and step_ variables, as well as the single fields
  // keep them in sync!
  IncrementTimeAndStep();

  if (printheader) PrintHeader();

  SetPoroSolution();
  scatra_->ScaTraField()->PrepareTimeStep();
  if (artery_coupl_) scatramsht_->ArtScatraField()->PrepareTimeStep();
  // set structure-based scalar transport values
  SetScatraSolution();

  poromulti_->PrepareTimeStep();
  SetPoroSolution();
  ApplyAdditionalDBCForVolFracSpecies();

  return;
}

/*----------------------------------------------------------------------*
 | prepare the time loop                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::PrepareTimeLoop()
{
  // set structure-based scalar transport values
  SetScatraSolution();
  poromulti_->PrepareTimeLoop();
  // initial output for scatra field
  SetPoroSolution();
  scatra_->ScaTraField()->Output();
  if (artery_coupl_) scatramsht_->ArtScatraField()->Output();

  return;
}

/*----------------------------------------------------------------------*
 | update fields and output results                         vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::UpdateAndOutput()
{
  // set scatra on fluid (necessary for possible domain integrals)
  SetScatraSolution();
  poromulti_->UpdateAndOutput();

  // scatra field
  scatra_->ScaTraField()->Update();
  scatra_->ScaTraField()->EvaluateErrorComparedToAnalyticalSol();
  scatra_->ScaTraField()->Output();
  // artery scatra field
  if (artery_coupl_)
  {
    scatramsht_->ArtScatraField()->Update();
    scatramsht_->ArtScatraField()->EvaluateErrorComparedToAnalyticalSol();
    scatramsht_->ArtScatraField()->Output();
  }
  if (Comm().MyPID() == 0)
  {
    std::cout << "Finished POROMULTIPHASESCATRA STEP " << std::setw(5) << std::setprecision(4)
              << std::scientific << Step() << "/" << std::setw(5) << std::setprecision(4)
              << std::scientific << NStep() << ": dtstep = " << dttimestep_ << std::endl;
  }
}

/*----------------------------------------------------------------------*
 | Test the results of all subproblems                       vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::CreateFieldTest()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  poromulti_->CreateFieldTest();
  problem->AddFieldTest(scatra_->CreateScaTraFieldTest());
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::SetPoroSolution()
{
  // safety check
  Teuchos::RCP<SCATRA::ScaTraTimIntPoroMulti> poroscatra =
      Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntPoroMulti>(scatra_->ScaTraField());
  if (poroscatra == Teuchos::null) dserror("cast to ScaTraTimIntPoroMulti failed!");

  // set displacements
  poroscatra->ApplyMeshMovement(poromulti_->StructDispnp(), 1);

  // set the fluid solution
  poroscatra->SetSolutionFieldOfMultiFluid(
      poromulti_->RelaxedFluidPhinp(), poromulti_->FluidField()->Phin(), 2);

  // additionally, set nodal flux if L2-projection is desired
  if (fluxreconmethod_ == INPAR::POROFLUIDMULTIPHASE::FluxReconstructionMethod::gradreco_l2)
    poroscatra->SetL2FluxOfMultiFluid(poromulti_->FluidFlux(), 1);

  if (artery_coupl_)
  {
    scatramsht_->SetArteryPressure();
    scatramsht_->ApplyMeshMovement();
  }
}

/*----------------------------------------------------------------------*
 |                                                     kremheller 06/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::ApplyAdditionalDBCForVolFracSpecies()
{
  // remove the old one
  ScatraAlgo()->ScaTraField()->RemoveDirichCond(add_dirichmaps_volfrac_spec_);

  std::vector<int> mydirichdofs(0);

  // get map and validdof-vector
  const Epetra_Map* elecolmap = ScatraAlgo()->ScaTraField()->Discretization()->ElementColMap();
  Teuchos::RCP<const Epetra_Vector> valid_volfracspec_dofs =
      PoroField()->FluidField()->ValidVolFracSpecDofs();

  // we identify the volume fraction species dofs which do not have a physical meaning and set a
  // DBC on them
  for (int iele = 0; iele < elecolmap->NumMyElements(); ++iele)
  {
    // dynamic_cast necessary because virtual inheritance needs runtime information
    DRT::ELEMENTS::Transport* myele = dynamic_cast<DRT::ELEMENTS::Transport*>(
        ScatraAlgo()->ScaTraField()->Discretization()->gElement(elecolmap->GID(iele)));

    const MAT::Material& material2 = *(myele->Material(2));

    // check the material
    if (material2.MaterialType() != INPAR::MAT::m_fluidporo_multiphase and
        material2.MaterialType() != INPAR::MAT::m_fluidporo_multiphase_reactions)
      dserror("only poro multiphase and poro multiphase reactions material valid");

    // cast fluid material
    const MAT::FluidPoroMultiPhase& multiphasemat =
        static_cast<const MAT::FluidPoroMultiPhase&>(material2);

    const int numfluidphases = multiphasemat.NumFluidPhases();
    const int numvolfrac = multiphasemat.NumVolFrac();
    const int numfluidmat = multiphasemat.NumMat();

    // this is only necessary if we have volume fractions present
    // TODO: this works only if we have the same number of phases in every element
    if (numfluidmat == numfluidphases) return;

    const MAT::Material& material = *(myele->Material());

    // cast scatra material
    const MAT::MatList& scatramat = static_cast<const MAT::MatList&>(material);

    if (not(scatramat.MaterialType() == INPAR::MAT::m_matlist or
            scatramat.MaterialType() == INPAR::MAT::m_matlist_reactions))
      dserror("wrong type of ScaTra-Material");

    const int numscatramat = scatramat.NumMat();

    DRT::Node** nodes = myele->Nodes();

    for (int inode = 0; inode < (myele->NumNode()); inode++)
    {
      if (nodes[inode]->Owner() == ScatraAlgo()->ScaTraField()->Discretization()->Comm().MyPID())
      {
        std::vector<int> scatradofs =
            ScatraAlgo()->ScaTraField()->Discretization()->Dof(0, nodes[inode]);
        std::vector<int> fluiddofs =
            PoroField()->FluidField()->Discretization()->Dof(0, nodes[inode]);

        for (int idof = 0; idof < numscatramat; ++idof)
        {
          int matid = scatramat.MatID(idof);
          Teuchos::RCP<MAT::Material> singlemat = scatramat.MaterialById(matid);
          if (singlemat->MaterialType() == INPAR::MAT::m_scatra_multiporo_fluid ||
              singlemat->MaterialType() == INPAR::MAT::m_scatra_multiporo_solid ||
              singlemat->MaterialType() == INPAR::MAT::m_scatra_multiporo_temperature)
          {
            // do nothing
          }
          else if (singlemat->MaterialType() == INPAR::MAT::m_scatra_multiporo_volfrac)
          {
            const Teuchos::RCP<const MAT::ScatraMatMultiPoroVolFrac>& scatravolfracmat =
                Teuchos::rcp_dynamic_cast<const MAT::ScatraMatMultiPoroVolFrac>(singlemat);

            const int scalartophaseid = scatravolfracmat->PhaseID();
            // if not already in original dirich map     &&   if it is not a valid volume fraction
            // species dof identified with < 1
            if (ScatraAlgo()->ScaTraField()->DirichMaps()->CondMap()->LID(scatradofs[idof]) == -1 &&
                (int)(*valid_volfracspec_dofs)
                        [PoroField()->FluidField()->Discretization()->DofRowMap()->LID(
                            fluiddofs[scalartophaseid + numvolfrac])] < 1)
            {
              mydirichdofs.push_back(scatradofs[idof]);
              ScatraAlgo()->ScaTraField()->Phinp()->ReplaceGlobalValue(scatradofs[idof], 0, 0.0);
            }
          }
          else
            dserror("only MAT_scatra_multiporo_(fluid,volfrac,solid,temperature) valid here");
        }
      }
    }
  }

  // build map
  int nummydirichvals = mydirichdofs.size();
  add_dirichmaps_volfrac_spec_ = Teuchos::rcp(new Epetra_Map(-1, nummydirichvals,
      &(mydirichdofs[0]), 0, ScatraAlgo()->ScaTraField()->Discretization()->Comm()));

  // add the condition
  ScatraAlgo()->ScaTraField()->AddDirichCond(add_dirichmaps_volfrac_spec_);

  return;
}

/*----------------------------------------------------------------------*
 |                                                         vuong 05/13  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::SetScatraSolution()
{
  poromulti_->SetScatraSolution(ndsporofluid_scatra_, scatra_->ScaTraField()->Phinp());
  if (artery_coupl_)
    poromulti_->FluidField()->ArtNetTimInt()->Discretization()->SetState(
        2, "one_d_artery_phinp", scatramsht_->ArtScatraField()->Phinp());
  return;
}

/*------------------------------------------------------------------------*
 | dof map of vector of unknowns of scatra field        kremheller 06/17  |
 *------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::ScatraDofRowMap()
    const
{
  return scatra_->ScaTraField()->DofRowMap();
}

/*------------------------------------------------------------------------*
 | handle divergence of nonlinear solver                kremheller 10/18  |
 *------------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraBase::HandleDivergence() const
{
  switch (divcontype_)
  {
    case INPAR::POROMULTIPHASESCATRA::divcont_continue:
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
    case INPAR::POROMULTIPHASESCATRA::divcont_stop:
    {
      dserror("POROMULTIPHASESCATRA nonlinear solver not converged in ITEMAX steps!");
      break;
    }
    default:
      dserror("unknown divercont action!");
      break;
  }
  return;
}
