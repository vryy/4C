/*----------------------------------------------------------------------*/
/*! \file

\brief scatra time integration for elch

\level 2


 *------------------------------------------------------------------------------------------------*/
#include "../drt_fluid/fluid_utils.H"  // for splitter

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/ion.H"
#include "../drt_mat/matlist.H"

#include "../drt_nurbs_discret/drt_nurbs_discret.H"

#include "../drt_scatra/scatra_timint_meshtying_strategy_fluid_elch.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i_elch.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_std_elch.H"

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_krylov_projector.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"

#include "scatra_timint_elch.H"
#include "scatra_timint_elch_service.H"

/*----------------------------------------------------------------------*
 | constructor                                              ehrl  01/14 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntElch::ScaTraTimIntElch(Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(dis, solver, sctratimintparams, extraparams, output),
      elchparams_(params),
      equpot_(DRT::INPUT::IntegralValue<INPAR::ELCH::EquPot>(*elchparams_, "EQUPOT")),
      frt_(elchparams_->get<double>("FARADAY_CONSTANT") /
           (elchparams_->get<double>("GAS_CONSTANT") * elchparams_->get<double>("TEMPERATURE"))),
      gstatnumite_(0),
      gstatincrement_(0.),
      dlcapexists_(false),
      ektoggle_(Teuchos::null),
      dctoggle_(Teuchos::null),
      electrodeinitvols_(Teuchos::null),
      electrodesoc_(Teuchos::null),
      electrodecrates_(Teuchos::null),
      electrodeconc_(Teuchos::null),
      electrodeeta_(Teuchos::null),
      electrodecurr_(Teuchos::null),
      cellvoltage_(0.),
      cellvoltage_old_(-1.),
      cccv_condition_(Teuchos::null),
      cellcrate_(0.),
      cellcrate_old_(-1.0),
      cycling_timestep_(elchparams_->get<double>("CYCLING_TIMESTEP")),
      adapted_timestep_active_(false),
      dt_adapted_(-1.0),
      splitter_macro_(Teuchos::null)
{
  // safety check
  if (frt_ <= 0.) dserror("Factor F/RT is non-positive!");
}


/*----------------------------------------------------------------------*
 | initialize algorithm                                 rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::Init()
{
  // The diffusion-conduction formulation does not support all options of the Nernst-Planck
  // formulation Let's check for valid options
  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "DIFFCOND_FORMULATION"))
    ValidParameterDiffCond();

  // additional safety checks associated with adaptive time stepping for CCCV cell cycling
  if (cycling_timestep_ > 0.)
  {
    if (not DRT::INPUT::IntegralValue<bool>(*params_, "ADAPTIVE_TIMESTEPPING"))
      dserror(
          "Adaptive time stepping for CCCV cell cycling requires ADAPTIVE_TIMESTEPPING flag to be "
          "set!");
    if (not discret_->GetCondition("CCCVCycling"))
      dserror(
          "Adaptive time stepping for CCCV cell cycling requires corresponding boundary "
          "condition!");
    if (not(cycling_timestep_ < dta_))
      dserror(
          "Adaptive time stepping for CCCV cell cycling requires that the modified time step size "
          "is smaller than the original time step size!");
  }
}

/*----------------------------------------------------------------------*
 |  initialize splitter for preconditioning                deanda 10/17 |
 *----------------------------------------------------------------------*/

void SCATRA::ScaTraTimIntElch::SetupSplitter()
{
  // set up concentration-potential splitter
  SetupConcPotSplit();

  // set up concentration-potential-potential splitter for macro scale in multi-scale simulations
  if (macro_scale_) SetupConcPotPotSplit();
}

/*----------------------------------------------------------------------*
 | setup algorithm                                          rauch 09/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::Setup()
{
  // set up concentration-potential splitter
  SetupSplitter();

  // initialize time-dependent electrode kinetics variables (galvanostatic mode or double layer
  // contribution)
  ComputeTimeDerivPot0(true);

  // initialize dirichlet toggle:
  // for certain ELCH problem formulations we have to provide
  // additional flux terms / currents across Dirichlet boundaries for the standard element call
  Teuchos::RCP<Epetra_Vector> dirichones = LINALG::CreateVector(*(dbcmaps_->CondMap()), false);
  dirichones->PutScalar(1.0);
  dctoggle_ = LINALG::CreateVector(*(discret_->DofRowMap()), true);
  dbcmaps_->InsertCondVector(dirichones, dctoggle_);

  // screen output (has to come after SetInitialField)
  // a safety check for the solver type
  if ((NumScal() > 1) && (solvtype_ != INPAR::SCATRA::solvertype_nonlinear))
    dserror("Solver type has to be set to >>nonlinear<< for ion transport.");

  if (myrank_ == 0)
  {
    std::cout << "\nSetup of splitter: numscal = " << NumScal() << std::endl;
    std::cout << "Temperature value T (Kelvin)     = " << elchparams_->get<double>("TEMPERATURE")
              << std::endl;
    std::cout << "Constant F/RT                    = " << frt_ << std::endl;
  }

  // initialize vectors for states of charge and C rates of resolved electrodes
  std::vector<DRT::Condition*> electrodesocconditions;
  discret_->GetCondition("ElectrodeSOC", electrodesocconditions);
  if (electrodesocconditions.size() > 0)
  {
    if (isale_)
    {
      electrodeinitvols_ = Teuchos::rcp(new std::vector<double>);
      electrodeinitvols_->resize(electrodesocconditions.size(), -1.);
    }
    electrodesoc_ = Teuchos::rcp(new std::vector<double>);
    electrodesoc_->resize(electrodesocconditions.size(), -1.);
    electrodecrates_ = Teuchos::rcp(new std::vector<double>);
    electrodecrates_->resize(electrodesocconditions.size(), -1.);

    // safety checks
    for (int i = 0; i < static_cast<int>(electrodesocconditions.size()); ++i)
    {
      const double one_hour = electrodesocconditions[i]->GetDouble("one_hour");

      if (one_hour <= 0.0) dserror("One hour must not be negative");
      if (std::fmod(std::log10(one_hour / 3600.0), 1.0) != 0)
        dserror("This is not one hour in SI units");
      if (i > 0)  // in case of more than one condition
        for (int j = 0; j < i; ++j)
          if (one_hour != electrodesocconditions[j]->GetDouble("one_hour"))
            dserror("Different definitions of one hour in Electrode STATE OF CHARGE CONDITIONS.");
    }
  }

  // initialize vectors for mean reactant concentrations, mean electric overpotentials, and total
  // electric currents at electrode boundaries
  std::vector<DRT::Condition*> electrodedomainconditions;
  discret_->GetCondition("ElchDomainKinetics", electrodedomainconditions);
  std::vector<DRT::Condition*> electrodeboundaryconditions;
  discret_->GetCondition("ElchBoundaryKinetics", electrodeboundaryconditions);
  std::vector<DRT::Condition*> electrodeboundarypointconditions;
  discret_->GetCondition("ElchBoundaryKineticsPoint", electrodeboundarypointconditions);
  if (electrodedomainconditions.size() > 0 and
      (electrodeboundaryconditions.size() > 0 or electrodeboundarypointconditions.size() > 0))
    dserror(
        "At the moment, we cannot have electrode domain kinetics conditions and electrode boundary "
        "kinetics conditions at the same time!");
  else if (electrodeboundaryconditions.size() > 0 and electrodeboundarypointconditions.size() > 0)
    dserror(
        "At the moment, we cannot have electrode boundary kinetics line/surface conditions and "
        "electrode boundary kinetics point conditions at the same time!");
  else if (electrodedomainconditions.size() > 0 or electrodeboundaryconditions.size() > 0 or
           electrodeboundarypointconditions.size() > 0)
  {
    const unsigned ncond = electrodedomainconditions.size() + electrodeboundaryconditions.size() +
                           electrodeboundarypointconditions.size();
    electrodeconc_ = Teuchos::rcp(new std::vector<double>(ncond, -1.));
    electrodeeta_ = Teuchos::rcp(new std::vector<double>(ncond, -1.));
    electrodecurr_ = Teuchos::rcp(new std::vector<double>(ncond, -1.));
  }

  // extract constant-current constant-voltage (CCCV) cell cycling and half-cycle boundary
  // conditions
  std::vector<Teuchos::RCP<DRT::Condition>> cccvcyclingconditions;
  discret_->GetCondition("CCCVCycling", cccvcyclingconditions);
  std::vector<Teuchos::RCP<DRT::Condition>> cccvhalfcycleconditions;
  discret_->GetCondition("CCCVHalfCycle", cccvhalfcycleconditions);

  switch (cccvcyclingconditions.size())
  {
    // no cell cycling intended
    case 0:
    {
      // safety check
      if (cccvhalfcycleconditions.size())
        dserror(
            "Found constant-current constant-voltage (CCCV) half-cycle boundary conditions, but no "
            "CCCV cell cycling condition!");

      break;
    }

    // cell cycling intended
    case 1:
    {
      // extract constant-current constant-voltage (CCCV) cell cycling boundary condition
      const DRT::Condition& cccvcyclingcondition = *cccvcyclingconditions[0];

      // safety checks
      if (NumDofPerNode() != 2)
        dserror(
            "Must have exactly two degrees of freedom per node (lithium concentration and electric "
            "potential) for lithium-ion cell cycling!");
      if (!cccvhalfcycleconditions.size())
        dserror(
            "Found constant-current constant-voltage (CCCV) cell cycling boundary condition, but "
            "no CCCV half-cycle boundary conditions!");
      if (cccvcyclingcondition.GetInt("ConditionIDForCharge") < 0 or
          cccvcyclingcondition.GetInt("ConditionIDForDischarge") < 0)
        dserror(
            "Invalid ID of constant-current constant-voltage (CCCV) half-cycle boundary condition "
            "specified in CCCV cell cycling boundary condition!");

      std::vector<DRT::Condition*> cccvhalfcycleconditions;
      discret_->GetCondition("CCCVHalfCycle", cccvhalfcycleconditions);

      // new cccv condition
      cccv_condition_ =
          Teuchos::rcp(new SCATRA::CCCVCondition(cccvcyclingcondition, cccvhalfcycleconditions,
              DRT::INPUT::IntegralValue<bool>(*params_, "ADAPTIVE_TIMESTEPPING")));

      break;
    }

    // safety check
    default:
    {
      dserror(
          "More than one constant-current constant-voltage (CCCV) cell cycling boundary "
          "condition is not allowed!");
      break;
    }
  }
}


/*----------------------------------------------------------------------*
 | set up concentration-potential splitter                   fang 08/17 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::SetupConcPotSplit()
{
  // prepare sets for concentration and potential dofs
  std::set<int> conddofset;
  std::set<int> otherdofset;

  // fill sets
  for (int inode = 0; inode < discret_->NumMyRowNodes(); ++inode)
  {
    std::vector<int> dofs = discret_->Dof(0, discret_->lRowNode(inode));
    for (unsigned idof = 0; idof < dofs.size(); ++idof)
      if (idof < (unsigned)NumScal())
        otherdofset.insert(dofs[idof]);
      else
        conddofset.insert(dofs[idof]);
  }

  // transform sets to maps
  std::vector<int> conddofmapvec(conddofset.begin(), conddofset.end());
  const Teuchos::RCP<const Epetra_Map> conddofmap = Teuchos::rcp(
      new Epetra_Map(-1, conddofmapvec.size(), &conddofmapvec[0], 0, discret_->Comm()));
  std::vector<int> otherdofmapvec(otherdofset.begin(), otherdofset.end());
  const Teuchos::RCP<const Epetra_Map> otherdofmap = Teuchos::rcp(
      new Epetra_Map(-1, otherdofmapvec.size(), &otherdofmapvec[0], 0, discret_->Comm()));

  // set up concentration-potential splitter
  splitter_ =
      Teuchos::rcp(new LINALG::MapExtractor(*discret_->DofRowMap(), conddofmap, otherdofmap));
}

/*---------------------------------------------------------------------*
 *---------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::SetupConcPotPotSplit()
{
  // prepare sets for dofs associated with electrolyte concentration, electrolyte potential, and
  // electrode potential
  std::set<int> dofset_conc_el;
  std::set<int> dofset_pot_el;
  std::set<int> dofset_pot_ed;

  // fill sets
  for (int inode = 0; inode < discret_->NumMyRowNodes(); ++inode)
  {
    std::vector<int> dofs = discret_->Dof(0, discret_->lRowNode(inode));
    for (unsigned idof = 0; idof < dofs.size(); ++idof)
      if (idof < (unsigned)NumScal())
        dofset_conc_el.insert(dofs[idof]);
      else if (idof == (unsigned)NumScal())
        dofset_pot_el.insert(dofs[idof]);
      else
        dofset_pot_ed.insert(dofs[idof]);
  }

  // transform sets to maps
  std::vector<Teuchos::RCP<const Epetra_Map>> maps(3, Teuchos::null);
  std::vector<int> dofmapvec_conc_el(dofset_conc_el.begin(), dofset_conc_el.end());
  maps[0] = Teuchos::rcp(
      new Epetra_Map(-1, dofmapvec_conc_el.size(), &dofmapvec_conc_el[0], 0, discret_->Comm()));
  std::vector<int> dofmapvec_pot_el(dofset_pot_el.begin(), dofset_pot_el.end());
  maps[1] = Teuchos::rcp(
      new Epetra_Map(-1, dofmapvec_pot_el.size(), &dofmapvec_pot_el[0], 0, discret_->Comm()));
  std::vector<int> dofmapvec_pot_ed(dofset_pot_ed.begin(), dofset_pot_ed.end());
  maps[2] = Teuchos::rcp(
      new Epetra_Map(-1, dofmapvec_pot_ed.size(), &dofmapvec_pot_ed[0], 0, discret_->Comm()));

  // set up concentration-potential-potential splitter
  splitter_macro_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*discret_->DofRowMap(), maps));
}


/*----------------------------------------------------------------------*
 | set elch-specific element parameters                      fang 11/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::SetElementSpecificScaTraParameters(
    Teuchos::ParameterList& eleparams) const
{
  // overwrite action type
  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "DIFFCOND_FORMULATION"))
  {
    eleparams.set<int>("action", SCATRA::set_diffcond_scatra_parameter);

    // parameters for diffusion-conduction formulation
    eleparams.sublist("DIFFCOND") = elchparams_->sublist("DIFFCOND");
  }
  else
    eleparams.set<int>("action", SCATRA::set_elch_scatra_parameter);

  // general elch parameters
  eleparams.set<double>("faraday", elchparams_->get<double>("FARADAY_CONSTANT"));
  eleparams.set<double>("gas_constant", elchparams_->get<double>("GAS_CONSTANT"));
  eleparams.set<double>("frt", FRT());
  eleparams.set<int>("equpot", equpot_);
  eleparams.set<bool>("boundaryfluxcoupling",
      DRT::INPUT::IntegralValue<bool>(*elchparams_, "COUPLE_BOUNDARY_FLUXES"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::ComputeTimeStepSize(double& dt)
{
  // call base class routine
  ScaTraTimIntImpl::ComputeTimeStepSize(dt);

  // adaptive time stepping for CCCV if activated
  if (cycling_timestep_ > 0.0)
  {
    // adaptive time stepping for CCCV cell cycling is currently inactive -> Check if it should be
    // activated
    if (!adapted_timestep_active_)
    {
      // check, current phase allows adaptive time stepping
      if (cccv_condition_->IsAdaptiveTimeSteppingPhase())
      {
        // extrapolate step and adapt time step if needed
        double dt_new = ExtrapolateStateAdaptTimeStep(dt);

        // activate adaptive time stepping and set new time step
        if (dt_new != dt)
        {
          // CCCV half cycle was not changed since this time step adaptivity. Thus, reset observer
          // (tracks phase changes)
          cccv_condition_->ResetPhaseChangeObserver();
          adapted_timestep_active_ = true;
          dt_adapted_ = dt = dt_new;
        }
      }
    }
    else
    {
      // if time step adaptivity is enabled for more than 3 steps after last change of phase:
      // disable, otherwise keep adapted time step
      if (cccv_condition_->IsStepsFromLastPhaseChange(3, step_))
        adapted_timestep_active_ = false;
      else
        dt = dt_adapted_;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double SCATRA::ScaTraTimIntElch::ExtrapolateStateAdaptTimeStep(double dt)
{
  // new time step;
  double dt_new = dt;

  // extrapolate state depending on current phase and check if it exceeds bounds of current phase.
  // If so, adapt time step
  switch (cccv_condition_->GetCCCVHalfCyclePhase())
  {
    case INPAR::ELCH::CCCVHalfCyclePhase::initital_relaxation:
    {
      const double time_new = time_ + 2 * dt;                  // extrapolate
      if (time_new >= cccv_condition_->GetInitialRelaxTime())  // check
      {
        const double timetoend = cccv_condition_->GetInitialRelaxTime() - time_;
        const int stepstoend = std::ceil(timetoend / cycling_timestep_);
        dt_new = timetoend / stepstoend;  // adapt
      }
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::constant_current:
    {
      // initialize variable for cell voltage from previous time step
      if (cellvoltage_old_ < 0.0) cellvoltage_old_ = cellvoltage_;

      const double cellvoltage_new =
          cellvoltage_ + 2.0 * (cellvoltage_ - cellvoltage_old_);  // extrapolate
      if (cccv_condition_->IsExceedCellVoltage(cellvoltage_new))   // check
      {
        dt_new = cycling_timestep_;  // adapt
        cellvoltage_old_ = -1.0;
      }
      else
        cellvoltage_old_ = cellvoltage_;
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::constant_voltage:
    {
      if (cellcrate_old_ < 0.0) cellcrate_old_ = cellcrate_;
      const double cellcrate_new = cellcrate_ + 2.0 * (cellcrate_ - cellcrate_old_);  // extrapolate
      if (cccv_condition_->IsExceedCellCCRate(cellcrate_new))                         // check
      {
        dt_new = cycling_timestep_;  // adapt
        cellcrate_old_ = -1.0;
      }
      else
        cellcrate_old_ = cellcrate_;
      break;
    }
    case INPAR::ELCH::CCCVHalfCyclePhase::relaxation:
    {
      const double time_new = time_ + 2 * dt;              // extrapolate
      if (time_new >= cccv_condition_->GetRelaxEndTime())  // check
      {
        const double timetoend = cccv_condition_->GetRelaxEndTime() - time_;
        const int stepstoend = std::ceil(timetoend / cycling_timestep_);
        dt_new = timetoend / stepstoend;  // adapt
      }
      break;
    }
    default:
    {
      dserror("Unknown phase of half cycle.");
      break;
    }
  }

  return dt_new;
}


/*----------------------------------------------------------------------*
 | add parameters depending on the problem              rasthofer 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::AddProblemSpecificParametersAndVectors(
    Teuchos::ParameterList& params  //!< parameter list
)
{
  discret_->SetState("dctoggle", dctoggle_);
}


/*----------------------------------------------------------------------*
 | contains the elch-specific nonlinear iteration loop       ehrl 01/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::NonlinearSolve()
{
  bool stopgalvanostat(false);
  gstatnumite_ = 1;

  // galvanostatic control (ELCH)
  while (!stopgalvanostat)
  {
    ScaTraTimIntImpl::NonlinearSolve();

    stopgalvanostat = ApplyGalvanostaticControl();
  }  // end galvanostatic control
}


/*----------------------------------------------------------------------*
 | assemble global system of equations                       fang 09/15 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::AssembleMatAndRHS()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // check for zero or negative concentration values
  CheckConcentrationValues(phinp_);

  // call base class routine
  ScaTraTimIntImpl::AssembleMatAndRHS();
}  // SCATRA::ScaTraTimIntElch::AssembleMatAndRHS


/*----------------------------------------------------------------------*
 | prepare time loop                                         fang 10/15 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::PrepareTimeLoop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // call base class routine
  ScaTraTimIntImpl::PrepareTimeLoop();

  // check validity of material and element formulation
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", SCATRA::check_scatra_element_parameter);
  if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}  // SCATRA::ScaTraTimIntElch::PrepareTimeLoop


/*------------------------------------------------------------------------------*
 | initialization procedure prior to evaluation of first time step   fang 09/15 |
 *------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::PrepareFirstTimeStep()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // calculate initial electric potential field
  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "INITPOTCALC")) CalcInitialPotentialField();

  // call base class routine
  ScaTraTimIntImpl::PrepareFirstTimeStep();

  // initialize Nernst boundary conditions
  InitNernstBC();
}  // SCATRA::ScaTraTimIntElch::PrepareFirstTimeStep

/*----------------------------------------------------------------------------------------*
 | create scalar manager                                                      fang 12/14 |
 *----------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::CreateScalarHandler()
{
  scalarhandler_ = Teuchos::rcp(new ScalarHandlerElch());
}  // ScaTraTimIntImpl::CreateMeshtyingStrategy

/*----------------------------------------------------------------------*
 |  calculate error compared to analytical solution            gjb 10/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::EvaluateErrorComparedToAnalyticalSol()
{
  switch (calcerror_)
  {
    case INPAR::SCATRA::calcerror_no:  // do nothing (the usual case)
      break;
    case INPAR::SCATRA::calcerror_Kwok_Wu:
    {
      //   References:

      //   Kwok, Yue-Kuen and Wu, Charles C. K.
      //   "Fractional step algorithm for solving a multi-dimensional
      //   diffusion-migration equation"
      //   Numerical Methods for Partial Differential Equations
      //   1995, Vol 11, 389-397

      //   G. Bauer, V. Gravemeier, W.A. Wall, A 3D finite element approach for the coupled
      //   numerical simulation of electrochemical systems and fluid flow,
      //   International Journal for Numerical Methods in Engineering, 86
      //   (2011) 1339–1359. DOI: 10.1002/nme.3107

      // create the parameters for the error calculation
      Teuchos::ParameterList eleparams;
      eleparams.set<int>("action", SCATRA::calc_error);
      eleparams.set("total time", time_);
      eleparams.set<int>("calcerrorflag", calcerror_);
      // provide displacement field in case of ALE
      if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("phinp", phinp_);

      // get (squared) error values
      Teuchos::RCP<Epetra_SerialDenseVector> errors = Teuchos::rcp(new Epetra_SerialDenseVector(3));
      discret_->EvaluateScalars(eleparams, errors);
      discret_->ClearState();

      double conerr1 = 0.0;
      double conerr2 = 0.0;
      // for the L2 norm, we need the square root
      if (NumScal() == 2)
      {
        conerr1 = sqrt((*errors)[0]);
        conerr2 = sqrt((*errors)[1]);
      }
      else if (NumScal() == 1)
      {
        conerr1 = sqrt((*errors)[0]);
        conerr2 = 0.0;
      }
      else
        dserror("The analytical solution of Kwok and Wu is only defined for two species");

      double poterr = sqrt((*errors)[2]);

      if (myrank_ == 0)
      {
        printf("\nL2_err for Kwok and Wu (time = %f):\n", time_);
        printf(" concentration1 %15.8e\n concentration2 %15.8e\n potential      %15.8e\n\n",
            conerr1, conerr2, poterr);
      }
#if 0
    if (myrank_ == 0)
    {
      // append error of the last time step to the error file
      if ((step_==stepmax_) or (time_==maxtime_))// write results to file
      {
        ostd::stringstream temp;
        const std::string simulation = problem_->OutputControlFile()->FileName();
        //const std::string fname = simulation+".relerror";
        const std::string fname = "XXX_kwok_xele.relerror";

        double elelength=0.0;
        if(simulation.find("5x5")!=std::string::npos)
          elelength=0.2;
        else if(simulation.find("10x10")!=std::string::npos)
          elelength=0.1;
        else if(simulation.find("20x20")!=std::string::npos)
          elelength=0.05;
        else if(simulation.find("40x40")!=std::string::npos)
          elelength=0.025;
        else if(simulation.find("50x50")!=std::string::npos)
          elelength=0.02;
        else if(simulation.find("80x80")!=std::string::npos)
          elelength=0.0125;
        else std::cout << "Warning: file name did not allow a evaluation of the element size!!!" << std::endl;

        std::ofstream f;
        f.precision(12);
        f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
        f << "#| " << simulation << "\n";
        //f << "#| Step | Time | rel. L2-error velocity mag |  rel. L2-error pressure  |\n";
        f << "#| Step | Time | c1 abs. error L2 | c2 abs. error L2 | phi abs. error L2 |  element length  |\n";
        //f << step_ << " " << time_ << " " << velerr/velint << " " << preerr/pint << " "<<"\n";
        f << step_ << " " << time_ << " " << conerr1 << " " << conerr2 << " " <<  poterr << " "<< elelength << "" << "\n";
        f.flush();
        f.close();
      }
    }
#endif
    }
    break;
    case INPAR::SCATRA::calcerror_cylinder:
    {
      //   Reference:
      //   G. Bauer, V. Gravemeier, W.A. Wall, A 3D finite element approach for the coupled
      //   numerical simulation of electrochemical systems and fluid flow,
      //   International Journal for Numerical Methods in Engineering, 2011

      // create the parameters for the error calculation
      Teuchos::ParameterList eleparams;
      eleparams.set<int>("action", SCATRA::calc_error);
      eleparams.set("total time", time_);
      eleparams.set<int>("calcerrorflag", calcerror_);
      // provide displacement field in case of ALE
      if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("phinp", phinp_);

      // get (squared) error values
      Teuchos::RCP<Epetra_SerialDenseVector> errors = Teuchos::rcp(new Epetra_SerialDenseVector(3));
      discret_->EvaluateScalars(eleparams, errors);
      discret_->ClearState();

      // for the L2 norm, we need the square root
      double conerr1 = sqrt((*errors)[0]);
      double conerr2 = sqrt((*errors)[1]);
      double poterr = sqrt((*errors)[2]);

      if (myrank_ == 0)
      {
        printf("\nL2_err for concentric cylinders (time = %f):\n", time_);
        printf(" concentration1 %15.8e\n concentration2 %15.8e\n potential      %15.8e\n\n",
            conerr1, conerr2, poterr);
      }
    }
    break;
    case INPAR::SCATRA::calcerror_electroneutrality:
    {
      // compute L2 norm of electroneutrality condition

      // create the parameters for the error calculation
      Teuchos::ParameterList eleparams;
      eleparams.set<int>("action", SCATRA::calc_error);
      eleparams.set("total time", time_);
      eleparams.set<int>("calcerrorflag", calcerror_);
      // provide displacement field in case of ALE
      if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

      // set vector values needed by elements
      discret_->ClearState();
      discret_->SetState("phinp", phinp_);

      // get (squared) error values
      Teuchos::RCP<Epetra_SerialDenseVector> errors = Teuchos::rcp(new Epetra_SerialDenseVector(1));
      discret_->EvaluateScalars(eleparams, errors);
      discret_->ClearState();

      // for the L2 norm, we need the square root
      double err = sqrt((*errors)[0]);

      if (myrank_ == 0)
      {
        printf("\nL2_err for electroneutrality (time = %f):\n", time_);
        printf(" Deviation from ENC: %15.8e\n\n", err);
      }
    }
    break;
    default:
    {
      // call base class routine
      ScaTraTimIntImpl::EvaluateErrorComparedToAnalyticalSol();
      break;
    }
  }
}  // SCATRA::ScaTraTimIntImpl::EvaluateErrorComparedToAnalyticalSol


/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::Update(const int num)
{
  // perform update of time-dependent electrode variables
  ElectrodeKineticsTimeUpdate();
}


/*----------------------------------------------------------------------*
 | output solution and restart data to file                  fang 02/18 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::Output(const int num)
{
  // call base class routine
  ScaTraTimIntImpl::Output(num);

  // output electrode interior status information and cell voltage in every time step
  // in the presence of constant-current constant-voltage (CCCV) cell cycling boundary condition
  if (not DoOutput() and discret_->GetCondition("CCCVCycling"))
  {
    // print electrode interior status information to screen and files
    OutputElectrodeInfoInterior();

    // print cell voltage to screen and file
    OutputCellVoltage();
  }
}


/*----------------------------------------------------------------------*
 | problem-specific outputs                                  fang 01/15 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::OutputProblemSpecific()
{
  // print electrode domain and boundary status information to screen and files
  OutputElectrodeInfoDomain();
  OutputElectrodeInfoBoundary();

  // print electrode interior status information to screen and files
  OutputElectrodeInfoInterior();

  // print cell voltage to screen and file
  OutputCellVoltage();

  // for elch problems with moving boundary
  if (isale_) output_->WriteVector("trueresidual", trueresidual_);
}  // SCATRA::ScaTraTimIntElch::OutputProblemSpecific()


/*----------------------------------------------------------------------*
 | read problem-specific restart data                        fang 01/17 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::ReadRestartProblemSpecific(
    const int step, IO::DiscretizationReader& reader)
{
  if (isale_) reader.ReadVector(trueresidual_, "trueresidual");

  // read restart data associated with electrode state of charge conditions if applicable, needed
  // for correct evaluation of cell C rate at the beginning of the first time step after restart
  if (discret_->GetCondition("ElectrodeSOC"))
  {
    // read volumes of resolved electrodes
    if (isale_) reader.ReadRedundantDoubleVector(electrodeinitvols_, "electrodeinitvols");

    // read states of charge of resolved electrodes
    reader.ReadRedundantDoubleVector(electrodesoc_, "electrodesoc");
  }

  // extract constant-current constant-voltage (CCCV) cell cycling boundary condition if available
  const DRT::Condition* const cccvcyclingcondition = discret_->GetCondition("CCCVCycling");

  // read restart data associated with constant-current constant-voltage (CCCV) cell cycling if
  // applicable
  if (cccvcyclingcondition)
  {
    // extract cell voltage
    cellvoltage_ = reader.ReadDouble("cellvoltage");

    // extract cell C rate
    cellcrate_ = reader.ReadDouble("cellcrate");

    // is time step adaptivity activated?
    adapted_timestep_active_ = static_cast<bool>(reader.ReadInt("adapted_timestep_active"));

    // adapted time step
    dt_adapted_ = reader.ReadDouble("dt_adapted");

    // read restart of cccv condition
    cccv_condition_->ReadRestart(reader);
  }
}  // SCATRA::ScaTraTimIntElch::ReadRestartProblemSpecific


/*------------------------------------------------------------------------------*
 | output electrode boundary status information to screen and file   fang 09/15 |
 *------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::OutputElectrodeInfoBoundary()
{
  // extract electrode boundary kinetics conditions from discretization
  std::vector<Teuchos::RCP<DRT::Condition>> cond;
  discret_->GetCondition("ElchBoundaryKinetics", cond);
  std::vector<Teuchos::RCP<DRT::Condition>> pointcond;
  discret_->GetCondition("ElchBoundaryKineticsPoint", pointcond);

  // safety check
  if (cond.size() and pointcond.size())
    dserror(
        "Cannot have electrode boundary kinetics point conditions and electrode boundary "
        "kinetics line/surface conditions at the same time!");

  // process conditions
  else if (cond.size() or pointcond.size())
  {
    double sum(0.0);

    if (myrank_ == 0)
    {
      std::cout << "Electrode boundary status information:" << std::endl;
      std::cout << "+----+------------------+-------------------+--------------------+-------------"
                   "--------+--------------------+---------------+----------------------+"
                << std::endl;
      std::cout << "| ID | reference domain | boundary integral | mean concentration | electrode "
                   "potential | mean overpotential | total current | mean current density |"
                << std::endl;
    }

    // evaluate the conditions and separate via ConditionID
    for (unsigned icond = 0; icond < cond.size() + pointcond.size(); ++icond)
    {
      // extract condition ID
      int condid(-1);
      if (cond.size())
        condid = cond[icond]->GetInt("ConditionID");
      else
        condid = pointcond[icond]->GetInt("ConditionID");

      // result vector
      // physical meaning of vector components is described in PostProcessSingleElectrodeInfo
      // routine
      Teuchos::RCP<Epetra_SerialDenseVector> scalars;

      // electrode boundary kinetics line/surface condition
      if (cond.size()) scalars = EvaluateSingleElectrodeInfo(condid, "ElchBoundaryKinetics");

      // electrode boundary kinetics point condition
      else
        scalars = EvaluateSingleElectrodeInfoPoint(pointcond[icond]);

      // initialize unused dummy variable
      double dummy(0.);

      PostProcessSingleElectrodeInfo(
          *scalars, condid, true, sum, dummy, dummy, dummy, dummy, dummy);
    }  // loop over conditions

    if (myrank_ == 0)
    {
      std::cout << "+----+------------------+-------------------+--------------------+-------------"
                   "--------+--------------------+---------------+----------------------+"
                << std::endl;
      // print out the net total current for all indicated boundaries
      printf("Net total current over boundary: %10.3E\n\n", sum);
    }

    // clean up
    discret_->ClearState();
  }
}  // SCATRA::ScaTraTimIntElch::OutputElectrodeInfoBoundary()


/*------------------------------------------------------------------------------*
 | evaluate status information on single line or surface electrode   fang 09/15 |
 *------------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_SerialDenseVector> SCATRA::ScaTraTimIntElch::EvaluateSingleElectrodeInfo(
    const unsigned condid,        //!< ID of condition to be evaluated
    const std::string condstring  //!< name of condition to be evaluated
)
{
  // set vector values needed by elements
  discret_->ClearState();
  discret_->SetState("phinp", phinp_);
  // needed for double-layer capacity!
  discret_->SetState("phidtnp", phidtnp_);

  // create parameter list
  Teuchos::ParameterList eleparams;

  // set action for elements depending on type of condition to be evaluated
  if (condstring == "ElchDomainKinetics")
    eleparams.set<int>("action", SCATRA::calc_elch_domain_kinetics);
  else if (condstring == "ElchBoundaryKinetics")
    eleparams.set<int>("action", SCATRA::bd_calc_elch_boundary_kinetics);
  else
    dserror("Invalid action " + condstring + " for output of electrode status information!");

  eleparams.set("calc_status", true);  // just want to have a status output!

  // provide displacement field in case of ALE
  if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

  // Since we just want to have the status output for t_{n+1},
  // we have to take care for Gen.Alpha!
  // AddTimeIntegrationSpecificVectors cannot be used since we do not want
  // an evaluation for t_{n+\alpha_f} !!!

  // Warning:
  // Specific time integration parameter are set in the following function.
  // In the case of a genalpha-time integration scheme the solution vector phiaf_ at time n+af
  // is passed to the element evaluation routine. Therefore, the electrode status is evaluate at a
  // different time (n+af) than our output routine (n+1), resulting in slightly different values
  // at the electrode. A different approach is not possible (without major hacks) since the
  // time-integration scheme is necessary to perform galvanostatic simulations, for instance.
  // Think about: double layer effects for genalpha time-integration scheme

  // add element parameters according to time-integration scheme
  AddTimeIntegrationSpecificVectors();

  // initialize result vector
  // physical meaning of vector components is described in PostProcessSingleElectrodeInfo routine
  Teuchos::RCP<Epetra_SerialDenseVector> scalars = Teuchos::rcp(new Epetra_SerialDenseVector(11));

  // evaluate relevant boundary integrals
  discret_->EvaluateScalars(eleparams, scalars, condstring, condid);

  // clean up
  discret_->ClearState();

  return scalars;
}


/*--------------------------------------------------------------------*
 | evaluate status information on single point electrode   fang 09/15 |
 *--------------------------------------------------------------------*/
Teuchos::RCP<Epetra_SerialDenseVector> SCATRA::ScaTraTimIntElch::EvaluateSingleElectrodeInfoPoint(
    Teuchos::RCP<DRT::Condition> condition  //!< condition to be evaluated
)
{
  // remove state vectors from discretization
  discret_->ClearState();

  // add state vectors to discretization
  discret_->SetState("phinp", phinp_);
  discret_->SetState("phidtnp", phidtnp_);  // needed for double layer capacity

  // add state vectors according to time integration scheme
  AddTimeIntegrationSpecificVectors();

  // determine number of scalar quantities to be computed
  const int numscalars = 11;

  // initialize result vector
  // physical meaning of vector components is described in PostProcessSingleElectrodeInfo routine
  Teuchos::RCP<Epetra_SerialDenseVector> scalars =
      Teuchos::rcp(new Epetra_SerialDenseVector(numscalars));

  // extract nodal cloud of current condition
  const std::vector<int>* nodeids = condition->Nodes();

  // safety checks
  if (!nodeids) dserror("Electrode kinetics point boundary condition doesn't have nodal cloud!");
  if (nodeids->size() != 1)
    dserror(
        "Electrode kinetics point boundary condition must be associated with exactly one node!");

  // extract global ID of conditioned node
  const int nodeid = (*nodeids)[0];

  // initialize variable for number of processor owning conditioned node
  int procid(-1);

  // consider node only if it is owned by current processor
  if (discret_->NodeRowMap()->MyGID(nodeid))
  {
    // extract number of processor owning conditioned node
    procid = discret_->Comm().MyPID();

    // create parameter list
    Teuchos::ParameterList condparams;

    // set action for elements
    condparams.set<int>("action", SCATRA::calc_elch_boundary_kinetics_point);

    // set flag for evaluation of status information
    condparams.set<bool>("calc_status", true);

    // equip element parameter list with current condition
    condparams.set<Teuchos::RCP<DRT::Condition>>("condition", condition);

    // provide displacement field in case of ALE
    if (isale_) condparams.set<int>("ndsdisp", nds_disp_);

    // get node
    DRT::Node* node = discret_->gNode(nodeid);

    // safety checks
    if (node == NULL) dserror("Cannot find node with global ID %d on discretization!", nodeid);
    if (node->NumElement() != 1)
      dserror(
          "Electrode kinetics point boundary condition must be specified on boundary node with "
          "exactly one attached element!");

    // get element attached to node
    DRT::Element* element = node->Elements()[0];

    // determine location information
    DRT::Element::LocationArray la(discret_->NumDofSets());
    element->LocationVector(*discret_, la, false);

    // dummy matrix and right-hand side vector
    Epetra_SerialDenseMatrix elematrix_dummy;
    Epetra_SerialDenseVector elevector_dummy;

    // evaluate electrode kinetics point boundary conditions
    const int error = element->Evaluate(condparams, *discret_, la, elematrix_dummy, elematrix_dummy,
        *scalars, elevector_dummy, elevector_dummy);

    // safety check
    if (error)
      dserror("Element with global ID %d returned error code %d on processor %d!", element->Id(),
          error, discret_->Comm().MyPID());

    // remove state vectors from discretization
    discret_->ClearState();
  }

  // communicate number of processor owning conditioned node
  int ownerid(-1);
  discret_->Comm().MaxAll(&procid, &ownerid, 1);

  // broadcast results from processor owning conditioned node to all other processors
  discret_->Comm().Broadcast(scalars->Values(), numscalars, ownerid);

  return scalars;
}


/*----------------------------------------------------------------------*
 | post-process status information on single electrode       fang 09/15 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::PostProcessSingleElectrodeInfo(
    Epetra_SerialDenseVector&
        scalars,           //!< scalar quantities associated with electrode status information (in)
    const unsigned id,     //!< electrode ID (in)
    const bool print,      //!< flag for output to screen and file (in)
    double& currentsum,    //!< net current involving all conditions (out)
    double& currtangent,   //!< tangent of current w.r.t. electrode potential (out)
    double& currresidual,  //!< negative residual of current equation (out)
    double& electrodeint,  //!< physical dimensions of the electrode region (out)
    double& electrodepot,  //!< electrode potential on electrode side (out)
    double& meanoverpot    //!< mean overpotential (out)
)
{
  // get total integral of current
  double currentintegral = scalars(0);
  // get total integral of double layer current
  double currentdlintegral = scalars(1);
  // get total domain or boundary integral
  double boundaryint = scalars(2);
  // get total integral of electric potential
  double electpotentialint = scalars(3);
  // get total integral of electric overpotential
  double overpotentialint = scalars(4);
  // get total integral of electric potential difference
  double epdint = scalars(5);
  // get total integral of open circuit electric potential
  double ocpint = scalars(6);
  // get total integral of reactant concentration
  double cint = scalars(7);
  // get derivative of integrated current with respect to electrode potential
  double currderiv = scalars(8);
  // get negative current residual (right-hand side of galvanostatic balance equation)
  double currentresidual = scalars(9);
  // get total domain integral scaled with volumetric electrode surface area total boundary
  // integral scaled with boundary porosity
  double boundaryint_porous = scalars(10);

  // specify some return values
  currentsum += currentintegral;  // sum of currents
  currtangent = currderiv;        // tangent w.r.t. electrode potential on metal side
  currresidual = currentresidual;
  electrodeint = boundaryint;
  electrodepot = electpotentialint / boundaryint;
  meanoverpot = overpotentialint / boundaryint;

  // print out results to screen/file if desired
  if (myrank_ == 0 and print)
  {
    // print out results to screen
    printf(
        "| %2d |      total       |    %10.3E     |     %10.3E     |     %10.3E      |     "
        "%10.3E     |  %10.3E   |      %10.3E      |\n",
        id, boundaryint, cint / boundaryint, electrodepot, overpotentialint / boundaryint,
        currentintegral + currentdlintegral,
        currentintegral / boundaryint + currentdlintegral / boundaryint);
    printf(
        "| %2d |   electrolyte    |    %10.3E     |     %10.3E     |     %10.3E      |     "
        "%10.3E     |  %10.3E   |      %10.3E      |\n",
        id, boundaryint_porous, cint / boundaryint_porous, electrodepot,
        overpotentialint / boundaryint, currentintegral + currentdlintegral,
        currentintegral / boundaryint_porous + currentdlintegral / boundaryint_porous);

    // write results to file
    std::ostringstream temp;
    temp << id;
    const std::string fname =
        problem_->OutputControlFile()->FileName() + ".electrode_status_" + temp.str() + ".txt";

    std::ofstream f;
    if (Step() == 0)
    {
      f.open(fname.c_str(), std::fstream::trunc);
      f << "#ID,Step,Time,Total_current,Boundary_integral,Mean_current_density_electrode_"
           "kinetics,Mean_current_density_dl,Mean_overpotential,Mean_electrode_pot_diff,Mean_"
           "opencircuit_pot,Electrode_pot,Mean_concentration,Boundary_integral_porous\n";
    }
    else
      f.open(fname.c_str(), std::fstream::ate | std::fstream::app);

    f << id << "," << Step() << "," << Time() << "," << currentintegral + currentdlintegral << ","
      << boundaryint << "," << currentintegral / boundaryint << ","
      << currentdlintegral / boundaryint << "," << overpotentialint / boundaryint << ","
      << epdint / boundaryint << "," << ocpint / boundaryint << "," << electrodepot << ","
      << cint / boundaryint << "," << boundaryint_porous << "\n";
    f.flush();
    f.close();
  }  // if (myrank_ == 0)

  // galvanostatic simulations:
  // add the double layer current to the Butler-Volmer current
  currentsum += currentdlintegral;

  // update vectors
  (*electrodeconc_)[id] = cint / boundaryint;
  (*electrodeeta_)[id] = overpotentialint / boundaryint;
  (*electrodecurr_)[id] = currentintegral + currentdlintegral;
}  // SCATRA::ScaTraTimIntElch::OutputSingleElectrodeInfoBoundary


/*-------------------------------------------------------------------------*
 | output electrode domain status information to screen         fang 07/15 |
 *-------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::OutputElectrodeInfoDomain()
{
  // extract electrode domain kinetics conditions from discretization
  std::string condstring("ElchDomainKinetics");
  std::vector<DRT::Condition*> conditions;
  discret_->GetCondition(condstring, conditions);

  // output electrode domain status information to screen if applicable
  if (conditions.size() > 0)
  {
    // initialize variable for total current
    double currentsum(0.);

    // print header of output table
    if (myrank_ == 0)
    {
      std::cout << "Status of '" << condstring << "':" << std::endl;
      std::cout << "+----+--------------------+---------------------+------------------+-----------"
                   "-----------+--------------------+----------------+----------------+"
                << std::endl;
      std::cout << "| ID | Bound/domain ratio |    Total current    | Domain integral  | Mean "
                   "current density | Mean overpotential | Electrode pot. | Mean Concentr. |"
                << std::endl;
    }

    // evaluate electrode domain kinetics conditions
    for (const auto& condition : conditions)
    {
      // extract condition ID
      const int condid = condition->GetInt("ConditionID");

      Teuchos::RCP<Epetra_SerialDenseVector> scalars =
          EvaluateSingleElectrodeInfo(condid, condstring);

      // initialize unused dummy variable
      double dummy(0.);

      PostProcessSingleElectrodeInfo(
          *scalars, condid, true, currentsum, dummy, dummy, dummy, dummy, dummy);
    }  // loop over electrode domain kinetics conditions

    if (myrank_ == 0)
    {
      // print finish line of output table
      std::cout << "+----+--------------------+----------------------+-----------------+-----------"
                   "-----------+--------------------+----------------+----------------+"
                << std::endl
                << std::endl;

      // print net total current
      printf("Net total current: %10.3E\n\n", currentsum);
    }

    // clean discretization
    discret_->ClearState();
  }
}


/*-------------------------------------------------------------------------------*
 | output electrode interior status information to screen and files   fang 01/15 |
 *-------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::OutputElectrodeInfoInterior()
{
  // extract conditions for electrode state of charge
  std::vector<DRT::Condition*> conditions;
  discret_->GetCondition("ElectrodeSOC", conditions);

  // perform all following operations only if there is at least one condition for electrode state
  // of charge
  if (conditions.size() > 0)
  {
    // initialize variable for cell C rate
    cellcrate_ = 0.;

    // print header to screen
    if (myrank_ == 0)
    {
      std::cout << std::endl << "Electrode state of charge and related:" << std::endl;
      std::cout << "+----+-----------------+----------------+----------------+" << std::endl;
      std::cout << "| ID | state of charge |     C rate     | operation mode |" << std::endl;
    }

    // loop over conditions for electrode state of charge
    for (const auto& condition : conditions)
    {
      // extract condition ID
      const int condid = condition->GetInt("ConditionID");

      // add state vectors to discretization
      discret_->ClearState();
      discret_->SetState("phinp", phinp_);
      discret_->SetState("phidtnp", phidtnp_);

      // create parameter list
      Teuchos::ParameterList condparams;

      // action for elements
      condparams.set<int>("action", SCATRA::calc_elch_electrode_soc_and_c_rate);

      // number of dofset associated with displacement-related dofs
      if (isale_) condparams.set<int>("ndsdisp", nds_disp_);

      // initialize result vector
      // first component  = integral of concentration
      // second component = integral of time derivative of concentration
      // third component  = integral of domain
      // fourth component = integral of velocity divergence (ALE only)
      // fifth component  = integral of concentration times velocity divergence (ALE only)
      // sixth component  = integral of velocity times concentration gradient (ALE only)
      Teuchos::RCP<Epetra_SerialDenseVector> scalars =
          Teuchos::rcp(new Epetra_SerialDenseVector(isale_ ? 6 : 3));

      // evaluate current condition for electrode state of charge
      discret_->EvaluateScalars(condparams, scalars, "ElectrodeSOC", condid);
      discret_->ClearState();

      // extract integral of domain
      const double intdomain = (*scalars)(2);

      // store initial volume of current electrode
      if (isale_ and step_ == 0) (*electrodeinitvols_)[condid] = intdomain;

      // extract reference concentrations at 0% and 100% state of charge
      const double volratio = isale_ ? (*electrodeinitvols_)[condid] / intdomain : 1.;
      const double c_0 = condition->GetDouble("c_0%") * volratio;
      const double c_100 = condition->GetDouble("c_100%") * volratio;
      const double c_delta_inv = 1. / (c_100 - c_0);

      // get one hour for c_rate
      const double one_hour = condition->GetDouble("one_hour");

      // compute state of charge and C rate for current electrode
      const double c_avg = (*scalars)(0) / intdomain;
      const double soc = (c_avg - c_0) * c_delta_inv;
      double c_rate = (*scalars)(1);
      if (isale_)  // ToDo: The ALE case is still doing some weird stuff (strong temporal
                   // oscillations of C rate), so one should have a closer look at that...
        c_rate += (*scalars)(4) + (*scalars)(5) - c_avg * (*scalars)(3);
      c_rate *= c_delta_inv * one_hour / intdomain;

      // determine operation mode
      std::string mode;
      if (std::abs(c_rate) < 1.e-16)
        mode.assign(" at rest ");
      else if (c_rate < 0.)
        mode.assign("discharge");
      else
        mode.assign(" charge  ");

      // update state of charge and C rate for current electrode
      (*electrodesoc_)[condid] = soc;
      (*electrodecrates_)[condid] = c_rate;

      // update cell C rate
      cellcrate_ = std::max(std::abs(c_rate), cellcrate_);

      // print results to screen and files
      if (myrank_ == 0)
      {
        // print results to screen
        std::cout << "| " << std::setw(2) << condid << " |   " << std::setw(7)
                  << std::setprecision(2) << std::fixed << soc * 100. << " %     |     "
                  << std::setw(5) << std::abs(c_rate) << "      |   " << mode.c_str() << "    |"
                  << std::endl;

        // set file name
        std::ostringstream id;
        id << condid;
        const std::string filename(
            problem_->OutputControlFile()->FileName() + ".electrode_soc_" + id.str() + ".csv");

        // open file in appropriate mode and write header at beginning
        std::ofstream file;
        if (Step() == 0)
        {
          file.open(filename.c_str(), std::fstream::trunc);
          file << "Step,Time,SOC,CRate" << std::endl;
        }
        else
          file.open(filename.c_str(), std::fstream::app);

        // write results for current electrode to file
        file << Step() << "," << Time() << "," << std::setprecision(9) << std::fixed << soc << ","
             << c_rate << std::endl;

        // close file
        file.close();
      }  // if(myrank_ == 0)
    }    // loop over conditions

    // print finish line to screen
    if (myrank_ == 0)
      std::cout << "+----+-----------------+----------------+----------------+" << std::endl;
  }
}  // SCATRA::ScaTraTimIntElch::OutputElectrodeInfoInterior


/*----------------------------------------------------------------------*
 | output cell voltage to screen and file                    fang 01/15 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::OutputCellVoltage()
{
  // extract conditions for cell voltage
  std::vector<DRT::Condition*> conditions;
  discret_->GetCondition("CellVoltage", conditions);
  std::vector<DRT::Condition*> conditionspoint;
  discret_->GetCondition("CellVoltagePoint", conditionspoint);
  if (conditions.size() and conditionspoint.size())
    dserror(
        "Cannot have cell voltage line/surface conditions and cell voltage point conditions at "
        "the "
        "same time!");
  else if (conditionspoint.size())
    conditions = conditionspoint;

  // perform all following operations only if there is at least one condition for cell voltage
  if (conditions.size() > 0)
  {
    // safety check
    if (conditions.size() != 2)
      dserror("Must have exactly two boundary conditions for cell voltage, one per electrode!");

    // print header
    if (myrank_ == 0)
    {
      std::cout << std::endl << "Electrode potentials and cell voltage:" << std::endl;
      std::cout << "+----+-------------------------+" << std::endl;
      std::cout << "| ID | mean electric potential |" << std::endl;
    }

    // initialize vector for mean electric potentials of electrodes
    std::vector<double> potentials(2, 0.);

    // loop over both conditions for cell voltage
    for (const auto& condition : conditions)
    {
      // extract condition ID
      const int condid = condition->GetInt("ConditionID");

      // process line and surface conditions
      if (conditionspoint.size() == 0)
      {
        // add state vector to discretization
        discret_->ClearState();
        discret_->SetState("phinp", phinp_);

        // create parameter list
        Teuchos::ParameterList condparams;

        // action for elements
        condparams.set<int>("action", SCATRA::bd_calc_elch_cell_voltage);

        // number of dofset associated with displacement-related dofs
        if (isale_) condparams.set<int>("ndsdisp", nds_disp_);

        // initialize result vector
        // first component = electric potential integral, second component = domain integral
        Teuchos::RCP<Epetra_SerialDenseVector> scalars =
            Teuchos::rcp(new Epetra_SerialDenseVector(2));

        // evaluate current condition for electrode state of charge
        discret_->EvaluateScalars(condparams, scalars, "CellVoltage", condid);
        discret_->ClearState();

        // extract concentration and domain integrals
        double intpotential = (*scalars)(0);
        double intdomain = (*scalars)(1);

        // compute mean electric potential of current electrode
        potentials[condid] = intpotential / intdomain;
      }

      // process point conditions
      else
      {
        // initialize local variable for electric potential of current electrode
        double potential(0.);

        // extract nodal cloud
        const std::vector<int>* const nodeids = condition->Nodes();
        if (nodeids == NULL) dserror("Cell voltage point condition does not have nodal cloud!");
        if (nodeids->size() != 1)
          dserror("Nodal cloud of cell voltage point condition must have exactly one node!");

        // extract node ID
        const int nodeid = (*nodeids)[0];

        // process row nodes only
        if (discret_->NodeRowMap()->MyGID(nodeid))
        {
          // extract node
          DRT::Node* node = discret_->gNode(nodeid);
          if (node == NULL)
            dserror("Cannot extract node with global ID %d from scalar transport discretization!",
                nodeid);

          // extract degrees of freedom from node
          const std::vector<int> dofs = discret_->Dof(0, node);

          // extract local ID of degree of freedom associated with electrode potential
          const int lid = discret_->DofRowMap()->LID(*dofs.rbegin());
          if (lid < 0)
            dserror("Cannot extract degree of freedom with global ID %d!", *dofs.rbegin());

          // extract electrode potential
          potential = (*phinp_)[lid];
        }

        // communicate electrode potential
        discret_->Comm().SumAll(&potential, &potentials[condid], 1);
      }

      // print mean electric potential of current electrode to screen
      if (myrank_ == 0)
        std::cout << "| " << std::setw(2) << condid << " |         " << std::setw(6)
                  << std::setprecision(3) << std::fixed << potentials[condid] << "          |"
                  << std::endl;
    }  // loop over conditions

    // compute cell voltage
    cellvoltage_ = std::abs(potentials[0] - potentials[1]);

    // print cell voltage to screen and file
    if (myrank_ == 0)
    {
      // print cell voltage to screen
      std::cout << "+----+-------------------------+" << std::endl;
      std::cout << "| cell voltage: " << std::setw(6) << cellvoltage_ << "         |" << std::endl;
      std::cout << "+----+-------------------------+" << std::endl;

      // set file name
      const std::string filename(problem_->OutputControlFile()->FileName() + ".cell_voltage.csv");

      // open file in appropriate mode and write header at beginning
      std::ofstream file;
      if (Step() == 0)
      {
        file.open(filename.c_str(), std::fstream::trunc);
        file << "Step,Time,CellVoltage" << std::endl;
      }
      else
        file.open(filename.c_str(), std::fstream::app);

      // write results for current electrode to file
      file << Step() << "," << Time() << "," << std::setprecision(12) << std::fixed << cellvoltage_
           << std::endl;

      // close file
      file.close();
    }  // if(myrank_ == 0)
  }
}  // SCATRA::ScaTraTimIntElch::OutputCellVoltage


/*----------------------------------------------------------------------*
 | output restart data                                       fang 01/17 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::OutputRestart() const
{
  // output restart data associated with electrode state of charge conditions if applicable,
  // needed for correct evaluation of cell C rate at the beginning of the first time step after
  // restart
  if (discret_->GetCondition("ElectrodeSOC"))
  {
    // output volumes of resolved electrodes
    if (isale_) output_->WriteRedundantDoubleVector("electrodeinitvols", electrodeinitvols_);

    // output states of charge of resolved electrodes
    output_->WriteRedundantDoubleVector("electrodesoc", electrodesoc_);
  }

  // output restart data associated with constant-current constant-voltage (CCCV) cell cycling if
  // applicable
  if (discret_->GetCondition("CCCVCycling"))
  {
    // output number of current charge or discharge half-cycle
    output_->WriteInt("ihalfcycle", cccv_condition_->GetNumCurrentHalfCycle());

    // output cell voltage
    output_->WriteDouble("cellvoltage", cellvoltage_);

    // output cell C rate
    output_->WriteDouble("cellcrate", cellcrate_);

    // was the phase changed since last time step adaptivity?
    output_->WriteInt("phasechanged", static_cast<int>(cccv_condition_->IsPhaseChanged()));

    // are we in initial phase relaxation?
    output_->WriteInt(
        "phaseinitialrelaxation", static_cast<int>(cccv_condition_->IsPhaseInitialRelaxation()));

    // end time of current relaxation phase
    output_->WriteDouble("relaxendtime", cccv_condition_->GetRelaxEndTime());

    // current phase of half cycle
    output_->WriteInt("phase_cccv", static_cast<int>(cccv_condition_->GetCCCVHalfCyclePhase()));

    // when was the phase change last time?
    output_->WriteInt("steplastphasechange", cccv_condition_->GetStepLastPhaseChange());

    // adapted time step
    output_->WriteDouble("dt_adapted", dt_adapted_);

    // is time step adaptivity activated?
    output_->WriteInt("adapted_timestep_active", adapted_timestep_active_);
  }
}  // SCATRA::ScaTraTimIntElch::OutputRestart()


/*----------------------------------------------------------------------*
 | perform setup of natural convection                       fang 08/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::SetupNatConv()
{
  // calculate the initial mean concentration value
  if (NumScal() < 1) dserror("Error since numscal = %d. Not allowed since < 1", NumScal());
  c0_.resize(NumScal());

  discret_->ClearState();
  discret_->SetState("phinp", phinp_);

  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", SCATRA::calc_total_and_mean_scalars);
  eleparams.set("inverting", false);

  // provide displacement field in case of ALE
  if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

  // evaluate integrals of concentrations and domain
  Teuchos::RCP<Epetra_SerialDenseVector> scalars =
      Teuchos::rcp(new Epetra_SerialDenseVector(NumDofPerNode() + 1));
  discret_->EvaluateScalars(eleparams, scalars);
  discret_->ClearState();  // clean up

  // calculate mean concentration
  const double domint = (*scalars)[NumDofPerNode()];

  if (std::abs(domint) < EPS15) dserror("Division by zero!");

  for (int k = 0; k < NumScal(); k++) c0_[k] = (*scalars)[k] / domint;

  // initialization of the densification coefficient vector
  densific_.resize(NumScal());
  DRT::Element* element = discret_->lRowElement(0);
  Teuchos::RCP<MAT::Material> mat = element->Material();

  if (mat->MaterialType() == INPAR::MAT::m_matlist)
  {
    Teuchos::RCP<const MAT::MatList> actmat = Teuchos::rcp_static_cast<const MAT::MatList>(mat);

    for (int k = 0; k < NumScal(); ++k)
    {
      const int matid = actmat->MatID(k);
      Teuchos::RCP<const MAT::Material> singlemat = actmat->MaterialById(matid);

      if (singlemat->MaterialType() == INPAR::MAT::m_ion)
      {
        Teuchos::RCP<const MAT::Ion> actsinglemat =
            Teuchos::rcp_static_cast<const MAT::Ion>(singlemat);

        densific_[k] = actsinglemat->Densification();

        if (densific_[k] < 0.0) dserror("received negative densification value");
      }
      else
        dserror("Material type is not allowed!");
    }
  }

  // for a single species calculation
  else if (mat->MaterialType() == INPAR::MAT::m_ion)
  {
    Teuchos::RCP<const MAT::Ion> actmat = Teuchos::rcp_static_cast<const MAT::Ion>(mat);

    densific_[0] = actmat->Densification();

    if (densific_[0] < 0.0) dserror("received negative densification value");
    if (NumScal() > 1) dserror("Single species calculation but numscal = %d > 1", NumScal());
  }
  else
    dserror("Material type is not allowed!");
}  // ScaTraTimIntElch::SetupNatConv()


/*-------------------------------------------------------------------------*
 | valid parameters for the diffusion-conduction formulation    ehrl 12/12 |
 *-------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::ValidParameterDiffCond()
{
  if (myrank_ == 0)
  {
    if (DRT::INPUT::IntegralValue<INPAR::ELCH::ElchMovingBoundary>(
            *elchparams_, "MOVINGBOUNDARY") != INPAR::ELCH::elch_mov_bndry_no)
      dserror("Moving boundaries are not supported in the ELCH diffusion-conduction framework!!");

    if (DRT::INPUT::IntegralValue<int>(*params_, "NATURAL_CONVECTION"))
      dserror("Natural convection is not supported in the ELCH diffusion-conduction framework!!");

    if (DRT::INPUT::IntegralValue<INPAR::SCATRA::SolverType>(*params_, "SOLVERTYPE") !=
            INPAR::SCATRA::solvertype_nonlinear and
        DRT::INPUT::IntegralValue<INPAR::SCATRA::SolverType>(*params_, "SOLVERTYPE") !=
            INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro and
        DRT::INPUT::IntegralValue<INPAR::SCATRA::SolverType>(*params_, "SOLVERTYPE") !=
            INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken and
        DRT::INPUT::IntegralValue<INPAR::SCATRA::SolverType>(*params_, "SOLVERTYPE") !=
            INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit and
        DRT::INPUT::IntegralValue<INPAR::SCATRA::SolverType>(*params_, "SOLVERTYPE") !=
            INPAR::SCATRA::solvertype_nonlinear_multiscale_microtomacro)
      dserror(
          "The only solvertype supported by the ELCH diffusion-conduction framework is the "
          "non-linear solver!!");

    if (problem_->GetProblemType() != prb_ssi and
        DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(*params_, "CONVFORM") !=
            INPAR::SCATRA::convform_convective)
      dserror("Only the convective formulation is supported so far!!");

    if ((DRT::INPUT::IntegralValue<int>(*params_, "NEUMANNINFLOW")) == true)
      dserror("Neuman inflow BC's are not supported by the ELCH diffusion-conduction framework!!");

    if ((DRT::INPUT::IntegralValue<int>(*params_, "CONV_HEAT_TRANS")) == true)
      dserror(
          "Convective heat transfer BC's are not supported by the ELCH diffusion-conduction "
          "framework!!");

    if ((DRT::INPUT::IntegralValue<INPAR::SCATRA::FSSUGRDIFF>(*params_, "FSSUGRDIFF")) !=
        INPAR::SCATRA::fssugrdiff_no)
      dserror("Subgrid diffusivity is not supported by the ELCH diffusion-conduction framework!!");

    if ((DRT::INPUT::IntegralValue<int>(*elchparams_, "BLOCKPRECOND")) == true)
      dserror("Block preconditioner is not supported so far!!");

    // Parameters defined in "SCALAR TRANSPORT DYNAMIC"
    Teuchos::ParameterList& scatrastabparams = params_->sublist("STABILIZATION");

    if ((DRT::INPUT::IntegralValue<INPAR::SCATRA::StabType>(scatrastabparams, "STABTYPE")) !=
        INPAR::SCATRA::stabtype_no_stabilization)
      dserror(
          "No stabilization is necessary for solving the ELCH diffusion-conduction framework!!");

    if ((DRT::INPUT::IntegralValue<INPAR::SCATRA::TauType>(scatrastabparams, "DEFINITION_TAU")) !=
        INPAR::SCATRA::tau_zero)
      dserror(
          "No stabilization is necessary for solving the ELCH diffusion-conduction framework!!");

    if ((DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalTau>(scatrastabparams, "EVALUATION_TAU")) !=
        INPAR::SCATRA::evaltau_integration_point)
      dserror("Evaluation of stabilization parameter only at Gauss points!!");

    if ((DRT::INPUT::IntegralValue<INPAR::SCATRA::EvalMat>(scatrastabparams, "EVALUATION_MAT")) !=
        INPAR::SCATRA::evalmat_integration_point)
      dserror("Evaluation of material only at Gauss points!!");

    if ((DRT::INPUT::IntegralValue<INPAR::SCATRA::Consistency>(scatrastabparams, "CONSISTENCY")) !=
        INPAR::SCATRA::consistency_no)
      dserror("Consistence formulation is not in the ELCH diffusion-conduction framework!!");

    if ((DRT::INPUT::IntegralValue<int>(scatrastabparams, "SUGRVEL")) == true)
      dserror("Subgrid velocity is not incoperated in the ELCH diffusion-conduction framework!!");

    if ((DRT::INPUT::IntegralValue<int>(scatrastabparams, "ASSUGRDIFF")) == true)
      dserror(
          "Subgrid diffusivity is not incoperated in the ELCH diffusion-conduction framework!!");
  }
}


/*----------------------------------------------------------------------*
 | Initialize Nernst-BC                                      ehrl 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::InitNernstBC()
{
  // access electrode kinetics condition
  std::vector<DRT::Condition*> Elchcond;
  discret_->GetCondition("ElchBoundaryKinetics", Elchcond);
  if (!Elchcond.size()) discret_->GetCondition("ElchBoundaryKineticsPoint", Elchcond);

  for (unsigned icond = 0; icond < Elchcond.size(); ++icond)
  {
    // check if Nernst-BC is defined on electrode kinetics condition
    if (Elchcond[icond]->GetInt("kinetic model") == INPAR::ELCH::nernst)
    {
      // safety check
      if (!Elchcond[icond]->GeometryDescription())
        dserror("Nernst boundary conditions not implemented for one-dimensional domains yet!");

      if (DRT::INPUT::IntegralValue<int>(*elchparams_, "DIFFCOND_FORMULATION"))
      {
        if (icond == 0) ektoggle_ = LINALG::CreateVector(*(discret_->DofRowMap()), true);

        // 1.0 for electrode-kinetics toggle
        const double one = 1.0;

        // global node id's which are part of the Nernst-BC
        const std::vector<int>* nodegids = Elchcond[icond]->Nodes();

        // loop over all global nodes part of the Nernst-BC
        for (int ii = 0; ii < (int(nodegids->size())); ++ii)
        {
          if (discret_->NodeRowMap()->MyGID((*nodegids)[ii]))
          {
            // get node with global node id (*nodegids)[ii]
            DRT::Node* node = discret_->gNode((*nodegids)[ii]);

            // get global dof ids of all dof's with global node id (*nodegids)[ii]
            std::vector<int> nodedofs = discret_->Dof(0, node);

            // define electrode kinetics toggle
            // later on this toggle is used to blanck the sysmat and rhs
            ektoggle_->ReplaceGlobalValues(1, &one, &nodedofs[NumScal()]);
          }
        }
      }
      else
        dserror("Nernst BC is only available for diffusion-conduction formulation!");
    }
  }

  // At element level the Nernst condition has to be handled like a DC
  if (ektoggle_ != Teuchos::null) dctoggle_->Update(1.0, *ektoggle_, 1.0);
}  // SCATRA::ScaTraTimIntImpl::InitNernstBC


/*----------------------------------------------------------------------------------------*
 | initialize meshtying strategy (including standard case without meshtying)   fang 12/14 |
 *----------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::CreateMeshtyingStrategy()
{
  // fluid meshtying
  if (msht_ != INPAR::FLUID::no_meshtying)
    strategy_ = Teuchos::rcp(new MeshtyingStrategyFluidElch(this));

  // scatra-scatra interface coupling
  else if (IsS2IMeshtying())
    strategy_ = Teuchos::rcp(new MeshtyingStrategyS2IElch(this, *params_));

  // standard case without meshtying
  else
    strategy_ = Teuchos::rcp(new MeshtyingStrategyStdElch(this));
}  // SCATRA::ScaTraTimIntElch::CreateMeshtyingStrategy()

/*----------------------------------------------------------------------*
 | calculate initial electric potential field                fang 03/15 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::CalcInitialPotentialField()
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + calc initial potential field");

  // safety checks
  dsassert(step_ == 0, "Step counter is not zero!");
  switch (equpot_)
  {
    case INPAR::ELCH::equpot_divi:
    case INPAR::ELCH::equpot_enc_pde:
    case INPAR::ELCH::equpot_enc_pde_elim:
    {
      // These stationary closing equations for the electric potential are OK, since they
      // explicitly contain the electric potential as variable and therefore can be solved for the
      // initial electric potential.
      break;
    }
    default:
    {
      // If the stationary closing equation for the electric potential does not explicitly contain
      // the electric potential as variable, we obtain a zero block associated with the electric
      // potential on the main diagonal of the global system matrix used below. This zero block
      // makes the entire global system matrix singular! In this case, it would be possible to
      // temporarily change the type of closing equation used, e.g., from INPAR::ELCH::equpot_enc
      // to INPAR::ELCH::equpot_enc_pde. This should work, but has not been implemented yet.
      dserror(
          "Initial potential field cannot be computed for chosen closing equation for electric "
          "potential!");
      break;
    }
  }

  // screen output
  if (myrank_ == 0)
  {
    std::cout << "SCATRA: calculating initial field for electric potential" << std::endl;
    PrintTimeStepInfo();
    std::cout << "+------------+-------------------+--------------+--------------+" << std::endl;
    std::cout << "|- step/max -|- tol      [norm] -|-- pot-res ---|-- pot-inc ---|" << std::endl;
  }

  // prepare Newton-Raphson iteration
  iternum_ = 0;
  const int itermax = params_->sublist("NONLINEAR").get<int>("ITEMAX");
  const double itertol = params_->sublist("NONLINEAR").get<double>("CONVTOL");
  const double restol = params_->sublist("NONLINEAR").get<double>("ABSTOLRES");

  // start Newton-Raphson iteration
  while (true)
  {
    // update iteration counter
    iternum_++;

    // check for non-positive concentration values
    CheckConcentrationValues(phinp_);

    // assemble global system matrix and residual vector
    AssembleMatAndRHS();
    strategy_->CondenseMatAndRHS(sysmat_, residual_);

    // project residual, such that only part orthogonal to nullspace is considered
    if (projector_ != Teuchos::null) projector_->ApplyPT(*residual_);

    // apply actual Dirichlet boundary conditions to system of equations
    LINALG::ApplyDirichlettoSystem(sysmat_, increment_, residual_, zeros_, *(dbcmaps_->CondMap()));

    // apply artificial Dirichlet boundary conditions to system of equations
    // to hold initial concentrations constant when solving for initial potential field
    LINALG::ApplyDirichlettoSystem(
        sysmat_, increment_, residual_, zeros_, *(splitter_->OtherMap()));

    // compute L2 norm of electric potential state vector
    Teuchos::RCP<Epetra_Vector> pot_vector = splitter_->ExtractCondVector(phinp_);
    double pot_state_L2(0.0);
    pot_vector->Norm2(&pot_state_L2);

    // compute L2 norm of electric potential residual vector
    splitter_->ExtractCondVector(residual_, pot_vector);
    double pot_res_L2(0.);
    pot_vector->Norm2(&pot_res_L2);

    // compute L2 norm of electric potential increment vector
    splitter_->ExtractCondVector(increment_, pot_vector);
    double pot_inc_L2(0.);
    pot_vector->Norm2(&pot_inc_L2);

    // care for the case that nothing really happens in the potential field
    if (pot_state_L2 < 1e-5) pot_state_L2 = 1.0;

    // first iteration step: solution increment is not yet available
    if (iternum_ == 1)
    {
      // print first line of convergence table to screen
      if (myrank_ == 0)
        std::cout << "|  " << std::setw(3) << iternum_ << "/" << std::setw(3) << itermax << "   | "
                  << std::setw(10) << std::setprecision(3) << std::scientific << itertol
                  << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific
                  << pot_res_L2 << "   |      --      | (      --     ,te=" << std::setw(10)
                  << std::setprecision(3) << std::scientific << dtele_ << ")" << std::endl;

      // absolute tolerance for deciding if residual is already zero
      // prevents additional solver calls that will not improve the residual anymore
      if (pot_res_L2 < restol)
      {
        // print finish line of convergence table to screen
        if (myrank_ == 0)
          std::cout << "+------------+-------------------+--------------+--------------+"
                    << std::endl
                    << std::endl;

        // abort Newton-Raphson iteration
        break;
      }
    }

    // later iteration steps: solution increment can be printed
    else
    {
      // print current line of convergence table to screen
      if (myrank_ == 0)
        std::cout << "|  " << std::setw(3) << iternum_ << "/" << std::setw(3) << itermax << "   | "
                  << std::setw(10) << std::setprecision(3) << std::scientific << itertol
                  << "[L_2 ]  | " << std::setw(10) << std::setprecision(3) << std::scientific
                  << pot_res_L2 << "   | " << std::setw(10) << std::setprecision(3)
                  << std::scientific << pot_inc_L2 / pot_state_L2 << "   | (ts=" << std::setw(10)
                  << std::setprecision(3) << std::scientific << dtsolve_ << ",te=" << std::setw(10)
                  << std::setprecision(3) << std::scientific << dtele_ << ")" << std::endl;

      // convergence check
      if ((pot_res_L2 <= itertol and pot_inc_L2 / pot_state_L2 <= itertol) or pot_res_L2 < restol)
      {
        // print finish line of convergence table to screen
        if (myrank_ == 0)
          std::cout << "+------------+-------------------+--------------+--------------+"
                    << std::endl
                    << std::endl;

        // abort Newton-Raphson iteration
        break;
      }
    }

    // warn if maximum number of iterations is reached without convergence
    if (iternum_ == itermax)
    {
      if (myrank_ == 0)
      {
        std::cout << "+--------------------------------------------------------------+"
                  << std::endl;
        std::cout << "|            >>>>>> not converged!                             |"
                  << std::endl;
        std::cout << "+--------------------------------------------------------------+" << std::endl
                  << std::endl;
      }

      // abort Newton-Raphson iteration
      break;
    }

    // safety checks
    if (std::isnan(pot_inc_L2) or std::isnan(pot_state_L2) or std::isnan(pot_res_L2))
      dserror("calculated vector norm is NaN.");
    if (std::isinf(pot_inc_L2) or std::isinf(pot_state_L2) or std::isinf(pot_res_L2))
      dserror("calculated vector norm is INF.");

    // zero out increment vector
    increment_->PutScalar(0.);

    // store time before solving global system of equations
    const double time = Teuchos::Time::wallTime();

    // reprepare Krylov projection if required
    if (updateprojection_) UpdateKrylovSpaceProjection();

    // solve final system of equations incrementally
    strategy_->Solve(solver_, sysmat_, increment_, residual_, phinp_, 1, projector_);

    // determine time needed for solving global system of equations
    dtsolve_ = Teuchos::Time::wallTime() - time;

    // update electric potential degrees of freedom in initial state vector
    splitter_->AddCondVector(splitter_->ExtractCondVector(increment_), phinp_);

    // copy initial state vector
    phin_->Update(1., *phinp_, 0.);

    // update state vectors for intermediate time steps (only for generalized alpha)
    ComputeIntermediateValues();
  }  // Newton-Raphson iteration

  // reset global system matrix and its graph, since we solved a very special problem with a
  // special sparsity pattern
  sysmat_->Reset();
}  // SCATRA::ScaTraTimIntElch::CalcInitialPotentialField()


/*----------------------------------------------------------------------*
 |  calculate conductivity of electrolyte solution             gjb 07/09|
 *----------------------------------------------------------------------*/
double SCATRA::ScaTraTimIntElch::ComputeConductivity(
    Epetra_SerialDenseVector& sigma,  //! result vector
    bool effCond,                     //! flag for computation of effective conductivity
    bool specresist                   //! flag for computation of specific electrolyte resistance
)
{
  // we perform the calculation on element level hiding the material access!
  // the initial concentration distribution has to be uniform to do so!!
  double specific_resistance = 0.0;

  // create the parameters for the elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", SCATRA::calc_elch_conductivity);

  if (effCond == true)
    eleparams.set("effCond", true);
  else
    eleparams.set("effCond", false);

  if (specresist == true)
    eleparams.set("specresist", true);
  else
    eleparams.set("specresist", false);

  // provide number of dofset associated with velocity-related dofs
  eleparams.set<int>("ndsvel", nds_vel_);

  // provide displacement field in case of ALE
  if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

  // set vector values needed by elements
  discret_->ClearState();
  AddTimeIntegrationSpecificVectors();

  // evaluate integrals of scalar(s) and domain
  Teuchos::RCP<Epetra_SerialDenseVector> sigma_domint =
      Teuchos::rcp(new Epetra_SerialDenseVector(NumScal() + 2));
  discret_->EvaluateScalars(eleparams, sigma_domint);
  discret_->ClearState();  // clean up
  const double domint = (*sigma_domint)[NumScal() + 1];

  if (specresist == false)
    for (int ii = 0; ii < NumScal() + 1; ++ii) sigma[ii] = (*sigma_domint)[ii] / domint;
  else
    specific_resistance = (*sigma_domint)[NumScal()] / domint;

  return specific_resistance;
}  // ScaTraTimIntImpl::ComputeConductivity


/*----------------------------------------------------------------------*
 | apply galvanostatic control                                gjb 11/09 |
 *----------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntElch::ApplyGalvanostaticControl()
{
  // for galvanostatic ELCH applications we have to adjust the
  // applied cell voltage and continue Newton-Raphson iterations until
  // we reach the desired value for the electric current.

  if (DRT::INPUT::IntegralValue<int>(*elchparams_, "GALVANOSTATIC"))
  {
    // set time derivative parameters of applied voltage for a double layer capacitance current
    // density,
    if (dlcapexists_) ComputeTimeDerivPot0(false);

    // extract electrode domain and boundary kinetics conditions from discretization
    std::vector<Teuchos::RCP<DRT::Condition>> electrodedomainconditions;
    discret_->GetCondition("ElchDomainKinetics", electrodedomainconditions);
    std::vector<Teuchos::RCP<DRT::Condition>> electrodeboundaryconditions;
    discret_->GetCondition("ElchBoundaryKinetics", electrodeboundaryconditions);
    std::vector<Teuchos::RCP<DRT::Condition>> electrodeboundarypointconditions;
    discret_->GetCondition("ElchBoundaryKineticsPoint", electrodeboundarypointconditions);

    // safety checks
    if (electrodedomainconditions.size() > 0 and
        (electrodeboundaryconditions.size() > 0 or electrodeboundarypointconditions.size() > 0))
      dserror(
          "At the moment, we cannot have electrode domain kinetics conditions and electrode "
          "boundary kinetics conditions at the same time!");
    else if (electrodeboundaryconditions.size() > 0 and electrodeboundarypointconditions.size() > 0)
      dserror(
          "At the moment, we cannot have electrode boundary kinetics line/surface conditions and "
          "electrode boundary kinetics point conditions at the same time!");

    // determine type of electrode kinetics conditions to be evaluated
    std::vector<Teuchos::RCP<DRT::Condition>> conditions;
    std::string condstring;
    if (electrodedomainconditions.size() > 0)
    {
      conditions = electrodedomainconditions;
      condstring = "ElchDomainKinetics";
    }
    else if (electrodeboundaryconditions.size() > 0)
    {
      conditions = electrodeboundaryconditions;
      condstring = "ElchBoundaryKinetics";
    }
    else if (electrodeboundarypointconditions.size() > 0)
    {
      conditions = electrodeboundarypointconditions;
      condstring = "ElchBoundaryPointKinetics";
    }
    else
      dserror("Must have electrode kinetics conditions for galvanostatics!");

    // evaluate electrode kinetics conditions if applicable
    if (!conditions.empty())
    {
      const int condid_cathode = elchparams_->get<int>("GSTATCONDID_CATHODE");
      const int condid_anode = elchparams_->get<int>("GSTATCONDID_ANODE");
      int gstatitemax = (elchparams_->get<int>("GSTATITEMAX"));
      double gstatcurrenttol = (elchparams_->get<double>("GSTATCURTOL"));
      const int curvenum = elchparams_->get<int>("GSTATFUNCTNO");
      const double tol = elchparams_->get<double>("GSTATCONVTOL");
      const double effective_length = elchparams_->get<double>("GSTAT_LENGTH_CURRENTPATH");
      if (effective_length < 0.0) dserror("A negative effective length is not possible!");
      const INPAR::ELCH::ApproxElectResist approxelctresist =
          DRT::INPUT::IntegralValue<INPAR::ELCH::ApproxElectResist>(
              *elchparams_, "GSTAT_APPROX_ELECT_RESIST");

      // There are maximal two electrode conditions by definition
      // current flow i at electrodes
      Teuchos::RCP<std::vector<double>> actualcurrent =
          Teuchos::rcp(new std::vector<double>(2, 0.0));
      // residual at electrodes = i*timefac
      Teuchos::RCP<std::vector<double>> currresidual =
          Teuchos::rcp(new std::vector<double>(2, 0.0));
      Teuchos::RCP<std::vector<double>> currtangent = Teuchos::rcp(new std::vector<double>(2, 0.0));
      Teuchos::RCP<std::vector<double>> electrodesurface =
          Teuchos::rcp(new std::vector<double>(2, 0.0));
      Teuchos::RCP<std::vector<double>> electrodepot =
          Teuchos::rcp(new std::vector<double>(2, 0.0));
      Teuchos::RCP<std::vector<double>> meanoverpot = Teuchos::rcp(new std::vector<double>(2, 0.0));
      double meanelectrodesurface(0.0);
      // Assumption: Residual at BV1 is the negative of the value at BV2, therefore only the first
      // residual is calculated
      double residual(0.0);

      // for all time integration schemes, compute the current value for phidtnp
      // this is needed for evaluating charging currents due to double-layer capacity
      // This may only be called here and not inside OutputSingleElectrodeInfoBoundary!!!!
      // Otherwise you modify your output to file called during Output()
      ComputeTimeDerivative();

      double targetcurrent = problem_->Funct(curvenum - 1).EvaluateTime(time_);
      double timefacrhs = 1.0 / ResidualScaling();

      double currtangent_anode(0.0);
      double currtangent_cathode(0.0);
      double potinc_ohm(0.0);
      double potdiffbulk(0.0);
      double resistance = 0.0;

      if (conditions.size() > 2)
        dserror(
            "The framework may not work for geometric setups containing more than two "
            "electrodes! \n If you need it, check the framework exactly!!");

      // loop over all BV
      // degenerated to a loop over 2 (user-specified) BV conditions
      // note: only the potential at the boundary with id condid_cathode will be adjusted!
      for (unsigned icond = 0; icond < conditions.size(); ++icond)
      {
        // extract condition ID
        const int condid = conditions[icond]->GetInt("ConditionID");

        // result vector
        // physical meaning of vector components is described in PostProcessSingleElectrodeInfo
        // routine
        Teuchos::RCP<Epetra_SerialDenseVector> scalars;

        // electrode boundary kinetics line/surface condition
        if (condstring != "ElchBoundaryPointKinetics")
          scalars = EvaluateSingleElectrodeInfo(condid, condstring);

        // electrode boundary kinetics point condition
        else
          scalars = EvaluateSingleElectrodeInfoPoint(conditions[icond]);

        PostProcessSingleElectrodeInfo(*scalars, condid, false, (*actualcurrent)[condid],
            (*currtangent)[condid], (*currresidual)[condid], (*electrodesurface)[condid],
            (*electrodepot)[condid], (*meanoverpot)[condid]);

        if (conditions.size() == 2)
        {
          // In the case the actual current is zero, we assume that the first electrode is the
          // cathode
          if ((*actualcurrent)[condid] < 0.0 and condid_cathode != condid)
            dserror(
                "The defined GSTATCONDID_CATHODE does not match the actual current flow "
                "situation!!");
          else if ((*actualcurrent)[condid] > 0.0 and condid_anode != condid)
            dserror(
                "The defined GSTATCONDID_ANODE does not match the actual current flow "
                "situation!!");
        }
      }  // end loop over electrode kinetics

      if (conditions.size() == 1)
      {
        if (condid_cathode != 0 or condid_anode != 1)
          dserror(
              "The defined GSTATCONDID_CATHODE and GSTATCONDID_ANODE is wrong for a setup with "
              "only one electrode!!\n Choose: GSTATCONDID_CATHODE=0 and GSTATCONDID_ANODE=1");
      }

      // get the applied electrode potential of the cathode
      int icond_cathode(-1);
      for (unsigned icond = 0; icond < conditions.size(); ++icond)
        if (conditions[icond]->GetInt("ConditionID") == condid_cathode)
        {
          icond_cathode = icond;
          break;
        }
      const double potold = conditions[icond_cathode]->GetDouble("pot");
      double potnew = potold;

      // bulk voltage loss
      // U = V_A - V_C =  eta_A + delta phi_ohm - eta_C
      // -> delta phi_ohm  = V_A - V_C - eta_A + eta_C = V_A - eta_A - (V_C  - eta_C)

      potdiffbulk = ((*electrodepot)[condid_anode] - (*meanoverpot)[condid_anode]) -
                    ((*electrodepot)[condid_cathode] - (*meanoverpot)[condid_cathode]);

      // cell voltage loss = V_A - V_C
      // potdiffcell=(*electrodepot)[condid_anode]-(*electrodepot)[condid_cathode];
      // tanget at anode and cathode
      currtangent_anode = (*currtangent)[condid_anode];
      currtangent_cathode = (*currtangent)[condid_cathode];

      if (conditions.size() == 2)
        // mean electrode surface of the cathode and anode
        meanelectrodesurface = ((*electrodesurface)[0] + (*electrodesurface)[1]) / 2;
      else
        meanelectrodesurface = (*electrodesurface)[condid_cathode];

      // The linearization of potential increment is always based on the cathode side!!

      // Assumption: Residual at BV1 is the negative of the value at BV2, therefore only the first
      // residual is calculated residual := (I - timefacrhs *I_target) I_target is alway negative,
      // since the reference electrode is the cathode
      residual = (*currresidual)[condid_cathode] -
                 (timefacrhs * targetcurrent);  // residual is stored only from cathode!

      // convergence test
      {
        if (myrank_ == 0)
        {
          // format output
          std::cout << "Galvanostatic mode:" << std::endl;
          std::cout << "+-----------------------------------------------------------------------+"
                    << std::endl;
          std::cout << "| Convergence check:                                                    |"
                    << std::endl;
          std::cout << "+-----------------------------------------------------------------------+"
                    << std::endl;
          std::cout << "| iteration:                                " << std::setw(14) << std::right
                    << gstatnumite_ << " / " << gstatitemax << "         |" << std::endl;
          std::cout << "| actual reaction current at cathode:            " << std::setprecision(6)
                    << std::scientific << std::setw(14) << std::right
                    << (*actualcurrent)[condid_cathode] << "         |" << std::endl;
          std::cout << "| required total current at cathode:             " << std::setw(14)
                    << std::right << targetcurrent << "         |" << std::endl;
          std::cout << "| negative residual (rhs):                       " << std::setw(14)
                    << std::right << residual << "         |" << std::endl;
          std::cout << "+-----------------------------------------------------------------------+"
                    << std::endl;
        }

        if (gstatnumite_ > gstatitemax)
        {
          if (myrank_ == 0)
          {
            std::cout << "| --> converged: maximum number iterations reached. Not yet converged!  |"
                      << std::endl;
            std::cout << "+-----------------------------------------------------------------------+"
                      << std::endl
                      << std::endl;
          }
          return true;  // we proceed to next time step
        }
        else if (abs(residual) < gstatcurrenttol)
        {
          if (myrank_ == 0)
          {
            std::cout << "| --> converged: Newton-RHS-Residual is smaller than " << gstatcurrenttol
                      << "!      |" << std::endl;
            std::cout << "+-----------------------------------------------------------------------+"
                      << std::endl
                      << std::endl;
          }
          return true;  // we proceed to next time step
        }
        // electric potential increment of the last iteration
        else if ((gstatnumite_ > 1) and
                 (abs(gstatincrement_) < (1 + abs(potold)) * tol))  // < ATOL + |pot|* RTOL
        {
          if (myrank_ == 0)
          {
            std::cout << "| --> converged: |" << gstatincrement_ << "| < "
                      << (1 + abs(potold)) * tol << std::endl;
            std::cout << "+-----------------------------------------------------------------------+"
                      << std::endl
                      << std::endl;
          }
          return true;  // galvanostatic control has converged
        }

        // safety check
        if (abs((*currtangent)[condid_cathode]) < EPS13)
          dserror(
              "Tangent in galvanostatic control is near zero: %lf", (*currtangent)[condid_cathode]);
      }

      // calculate the cell potential increment due to ohmic resistance
      if (approxelctresist == INPAR::ELCH::approxelctresist_effleninitcond)
      {
        // update applied electric potential
        // potential drop ButlerVolmer conditions (surface ovepotential) and in the electrolyte
        // (ohmic overpotential) are conected in parallel:

        // 3 different versions:
        // I_0 = I_BV1 = I_ohmic = I_BV2
        // R(I_target, I) = R_BV1(I_target, I) = R_ohmic(I_target, I) = -R_BV2(I_target, I)
        // delta E_0 = delta U_BV1 + delta U_ohmic - (delta U_BV2)
        // => delta E_0 = (R_BV1(I_target, I)/J) + (R_ohmic(I_target, I)/J) - (-R_BV2(I_target,
        // I)/J) Attention: epsilon and tortuosity are missing in this framework
        //            -> use approxelctresist_efflenintegcond or approxelctresist_relpotcur

        // initialize conductivity vector
        Epetra_SerialDenseVector sigma(NumDofPerNode());

        // compute conductivity
        ComputeConductivity(sigma);

        // print conductivity
        if (myrank_ == 0)
        {
          for (int k = 0; k < NumScal(); ++k)
            std::cout << "| Electrolyte conductivity (species " << k + 1
                      << "):          " << std::setw(14) << std::right << sigma[k] << "         |"
                      << std::endl;

          if (equpot_ == INPAR::ELCH::equpot_enc_pde_elim)
          {
            double diff = sigma[0];

            for (int k = 1; k < NumScal(); ++k) diff += sigma[k];

            std::cout << "| Electrolyte conductivity (species elim) = " << sigma[NumScal()] - diff
                      << "         |" << std::endl;
          }

          std::cout << "| Electrolyte conductivity (all species):        " << std::setw(14)
                    << std::right << sigma[NumScal()] << "         |" << std::endl;
          std::cout << "+-----------------------------------------------------------------------+"
                    << std::endl;
        }

        // compute electrolyte resistance
        resistance = effective_length / (sigma[NumScal()] * meanelectrodesurface);
      }

      else if (approxelctresist == INPAR::ELCH::approxelctresist_relpotcur and
               conditions.size() == 2)
      {
        // actual potential difference is used to calculate the current path length
        // -> it is possible to compute the new ohmic potential step (porous media are
        // automatically included)
        //    without the input parameter GSTAT_LENGTH_CURRENTPATH
        // actual current < 0,  since the reference electrode is the cathode
        // potdiffbulk > 0,     always positive (see definition)
        // -1.0,                resistance has to be positive
        resistance = -1.0 * (potdiffbulk / (*actualcurrent)[condid_cathode]);
        // use of target current for the estimation of the resistance
        // resistance = -1.0*(potdiffbulk/(targetcurrent));
      }

      else if (approxelctresist == INPAR::ELCH::approxelctresist_efflenintegcond and
               conditions.size() == 2)
      {
        // dummy conductivity vector
        Epetra_SerialDenseVector sigma;

        const double specificresistance = ComputeConductivity(sigma, true, true);

        resistance = specificresistance * effective_length / meanelectrodesurface;
        // actual current < 0,  since the reference electrode is the cathode
        // potdiffbulk > 0,     always positive (see definition)
        // -1.0,                resistance has to be positive
      }

      else
        dserror(
            "The combination of the parameter GSTAT_APPROX_ELECT_RESIST %i and the number of "
            "electrodes %i\n is not valid!",
            approxelctresist, conditions.size());

      // calculate increment due to ohmic resistance
      potinc_ohm = -1.0 * resistance * residual / timefacrhs;

      // Do not update the cell potential for small currents
      if (abs((*actualcurrent)[condid_cathode]) < EPS10) potinc_ohm = 0.0;

      // the current flow at both electrodes has to be the same within the solution tolerances
      if (abs((*actualcurrent)[condid_cathode] + (*actualcurrent)[condid_anode]) > EPS8)
      {
        if (myrank_ == 0)
        {
          std::cout << "| WARNING: The difference of the current flow at anode and cathode      |"
                    << std::endl;
          std::cout << "| is "
                    << abs((*actualcurrent)[condid_cathode] + (*actualcurrent)[condid_anode])
                    << " larger than " << EPS8 << "!                             |" << std::endl;
          std::cout << "+-----------------------------------------------------------------------+"
                    << std::endl;
        }
      }

      // Newton step:  Jacobian * \Delta pot = - Residual
      const double potinc_cathode = residual / ((-1) * currtangent_cathode);
      double potinc_anode = 0.0;
      if (abs(currtangent_anode) > EPS13)  // anode surface overpotential is optional
      {
        potinc_anode = residual / ((-1) * currtangent_anode);
      }
      gstatincrement_ = (potinc_cathode + potinc_anode + potinc_ohm);
      // update electric potential
      potnew += gstatincrement_;

      if (myrank_ == 0)
      {
        std::cout << "| The ohmic potential increment is calculated based on                  |"
                  << std::endl;
        std::cout << "| the ohmic electrolyte resistance obtained from                        |"
                  << std::endl;
        ;
        if (approxelctresist == INPAR::ELCH::approxelctresist_effleninitcond)
          std::cout << "| GSTAT_LENGTH_CURRENTPATH and the averaged electrolyte conductivity.   |"
                    << std::endl;
        else if (approxelctresist == INPAR::ELCH::approxelctresist_relpotcur)
          std::cout << "| the applied potential and the resulting current flow.                 |"
                    << std::endl;
        else
          std::cout << "| GSTAT_LENGTH_CURRENTPATH and the integrated electrolyte conductivity. |"
                    << std::endl;

        std::cout << "+-----------------------------------------------------------------------+"
                  << std::endl;
        std::cout << "| Defined GSTAT_LENGTH_CURRENTPATH:              " << std::setw(14)
                  << std::right << effective_length << "         |" << std::endl;
        std::cout << "| Approximate electrolyte resistance:            " << std::setw(14)
                  << std::right << resistance << "         |" << std::endl;
        std::cout << "| New guess for:                                                        |"
                  << std::endl;
        std::cout << "| - ohmic potential increment:                   " << std::setw(14)
                  << std::right << potinc_ohm << "         |" << std::endl;
        std::cout << "| - overpotential increment cathode (condid " << condid_cathode
                  << "):  " << std::setw(14) << std::right << potinc_cathode << "         |"
                  << std::endl;
        std::cout << "| - overpotential increment anode (condid " << condid_anode
                  << "):    " << std::setw(14) << std::right << potinc_anode << "         |"
                  << std::endl;
        std::cout << "| -> total increment for potential:              " << std::setw(14)
                  << std::right << gstatincrement_ << "         |" << std::endl;
        std::cout << "+-----------------------------------------------------------------------+"
                  << std::endl;
        std::cout << "| old potential at the cathode (condid " << condid_cathode
                  << "):       " << std::setw(14) << std::right << potold << "         |"
                  << std::endl;
        std::cout << "| new potential at the cathode (condid " << condid_cathode
                  << "):       " << std::setw(14) << std::right << potnew << "         |"
                  << std::endl;
        std::cout << "+-----------------------------------------------------------------------+"
                  << std::endl
                  << std::endl;
      }

      //      // print additional information
      //      if (myrank_==0)
      //      {
      //        std::cout<< "  actualcurrent - targetcurrent = " <<
      //        ((*actualcurrent)[condid_cathode]-targetcurrent) << std::endl; std::cout<< "
      //        currtangent_cathode:  " << currtangent_cathode << std::endl; std::cout<< "
      //        currtangent_anode:    " << currtangent_anode << std::endl; std::cout<< "
      //        actualcurrent cathode:    " << (*actualcurrent)[condid_cathode] << std::endl;
      //        std::cout<< "  actualcurrent anode:  " << (*actualcurrent)[condid_anode] <<
      //        std::endl;
      //      }

      // replace potential value of the boundary condition (on all processors)
      conditions[icond_cathode]->Add("pot", potnew);
      gstatnumite_++;
      return false;  // not yet converged -> continue Newton iteration with updated potential
    }
  }
  return true;  // default

}  // end ApplyGalvanostaticControl()


/*----------------------------------------------------------------------------*
 | evaluate domain or boundary conditions for electrode kinetics   fang 06/15 |
 *----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::EvaluateElectrodeKineticsConditions(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,  //!< global system matrix
    Teuchos::RCP<Epetra_Vector> rhs,                    //!< global right-hand side vector
    const std::string condstring                        //!< name of condition to be evaluated
)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition '" + condstring + "'");

  discret_->ClearState();

  // create parameter list
  Teuchos::ParameterList condparams;

  // set action for elements depending on type of condition to be evaluated
  if (condstring == "ElchDomainKinetics")
    condparams.set<int>("action", SCATRA::calc_elch_domain_kinetics);
  else if (condstring == "ElchBoundaryKinetics")
    condparams.set<int>("action", SCATRA::bd_calc_elch_boundary_kinetics);
  else
    dserror("Illegal action for electrode kinetics evaluation!");

  if (isale_)  // provide displacement field in case of ALE
    condparams.set<int>("ndsdisp", nds_disp_);

  // add element parameters and set state vectors according to time-integration scheme
  AddTimeIntegrationSpecificVectors();

  // evaluate electrode kinetics conditions at time t_{n+1} or t_{n+alpha_F}
  discret_->EvaluateCondition(
      condparams, systemmatrix, Teuchos::null, rhs, Teuchos::null, Teuchos::null, condstring);
  discret_->ClearState();

  // add linearization of NernstCondition to system matrix
  if (ektoggle_ != Teuchos::null) LinearizationNernstCondition();
}  // ScaTraTimIntElch::EvaluateElectrodeKineticsConditions


/*----------------------------------------------------------------------------*
 | evaluate point boundary conditions for electrode kinetics       fang 07/15 |
 *----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::EvaluateElectrodeBoundaryKineticsPointConditions(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,  //!< global system matrix
    Teuchos::RCP<Epetra_Vector> rhs                     //!< global right-hand side vector
)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + evaluate condition 'ElchBoundaryKineticsPoint'");

  // remove state vectors from discretization
  discret_->ClearState();

  // create parameter list
  Teuchos::ParameterList condparams;

  // set action for elements
  condparams.set<int>("action", SCATRA::calc_elch_boundary_kinetics_point);

  // provide displacement field in case of ALE
  if (isale_) condparams.set<int>("ndsdisp", nds_disp_);

  // set state vectors according to time-integration scheme
  AddTimeIntegrationSpecificVectors();

  // extract electrode kinetics point boundary conditions from discretization
  std::vector<Teuchos::RCP<DRT::Condition>> conditions;
  discret_->GetCondition("ElchBoundaryKineticsPoint", conditions);

  // loop over all electrode kinetics point boundary conditions
  for (const auto& condition : conditions)
  {
    // extract nodal cloud of current condition
    const std::vector<int>* nodeids = condition->Nodes();

    // safety checks
    if (!nodeids) dserror("Electrode kinetics point boundary condition doesn't have nodal cloud!");
    if (nodeids->size() != 1)
      dserror(
          "Electrode kinetics point boundary condition must be associated with exactly one node!");

    // extract global ID of conditioned node
    const int nodeid = (*nodeids)[0];

    // consider node only if it is owned by current processor
    if (discret_->NodeRowMap()->MyGID(nodeid))
    {
      // equip element parameter list with current condition
      condparams.set<Teuchos::RCP<DRT::Condition>>("condition", condition);

      // get node
      DRT::Node* node = discret_->gNode(nodeid);

      // safety checks
      if (node == NULL) dserror("Cannot find node with global ID %d on discretization!", nodeid);
      if (node->NumElement() != 1)
        dserror(
            "Electrode kinetics point boundary condition must be specified on boundary node with "
            "exactly one attached element!");

      // get element attached to node
      DRT::Element* element = node->Elements()[0];

      // determine location information
      DRT::Element::LocationArray la(discret_->NumDofSets());
      element->LocationVector(*discret_, la, false);

      // initialize element matrix
      const int size = (int)la[0].lm_.size();
      Epetra_SerialDenseMatrix elematrix;
      if (elematrix.M() != size)
        elematrix.Shape(size, size);
      else
        memset(elematrix.A(), 0, size * size * sizeof(double));

      // initialize element right-hand side vector
      Epetra_SerialDenseVector elevector;
      if (elevector.Length() != size)
        elevector.Size(size);
      else
        memset(elevector.Values(), 0, size * sizeof(double));

      // dummy matrix and right-hand side vector
      Epetra_SerialDenseMatrix elematrix_dummy;
      Epetra_SerialDenseVector elevector_dummy;

      // evaluate electrode kinetics point boundary conditions
      const int error = element->Evaluate(condparams, *discret_, la, elematrix, elematrix_dummy,
          elevector, elevector_dummy, elevector_dummy);

      // safety check
      if (error)
        dserror("Element with global ID %d returned error code %d on processor %d!", element->Id(),
            error, discret_->Comm().MyPID());

      // assemble element matrix and right-hand side vector into global system of equations
      sysmat_->Assemble(element->Id(), la[0].stride_, elematrix, la[0].lm_, la[0].lmowner_);
      LINALG::Assemble(*residual_, elevector, la[0].lm_, la[0].lmowner_);
    }
  }

  // remove state vectors from discretization
  discret_->ClearState();
}  // SCATRA::ScaTraTimIntElch::EvaluateElectrodeBoundaryKineticsPointConditions


/*----------------------------------------------------------------------*
 | Add Linearization for Nernst-BC                           ehrl 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::LinearizationNernstCondition()
{
  // Blank rows with Nernst-BC (inclusive diagonal entry)
  // Nernst-BC is a additional constraint coupled to the original system of equation
  if (!sysmat_->Filled()) sysmat_->Complete();
  sysmat_->ApplyDirichlet(ektoggle_, false);
  LINALG::ApplyDirichlettoSystem(increment_, residual_, zeros_, ektoggle_);

  discret_->ClearState();

  // create an parameter list
  Teuchos::ParameterList condparams;
  // update total time for time curve actions
  AddTimeIntegrationSpecificVectors();
  // action for elements
  condparams.set<int>("action", SCATRA::bd_calc_elch_linearize_nernst);

  // add element parameters and set state vectors according to time-integration scheme
  // we need here concentration at t+np
  discret_->SetState("phinp", phinp_);

  std::string condstring("ElchBoundaryKinetics");
  // evaluate ElchBoundaryKinetics conditions at time t_{n+1} or t_{n+alpha_F}
  // phinp (view to phinp)
  discret_->EvaluateCondition(
      condparams, sysmat_, Teuchos::null, residual_, Teuchos::null, Teuchos::null, condstring);
  discret_->ClearState();
}  //  SCATRA::ScaTraTimIntImpl::LinearizationNernstCondition()


/*----------------------------------------------------------------------------*
 | evaluate solution-depending boundary and interface conditions   fang 10/14 |
 *----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::EvaluateSolutionDependingConditions(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,  //!< system matrix
    Teuchos::RCP<Epetra_Vector> rhs                     //!< rhs vector
)
{
  // evaluate domain conditions for electrode kinetics
  if (discret_->GetCondition("ElchDomainKinetics") != NULL)
    EvaluateElectrodeKineticsConditions(systemmatrix, rhs, "ElchDomainKinetics");

  // evaluate boundary conditions for electrode kinetics
  if (discret_->GetCondition("ElchBoundaryKinetics") != NULL)
    EvaluateElectrodeKineticsConditions(systemmatrix, rhs, "ElchBoundaryKinetics");

  // evaluate point boundary conditions for electrode kinetics
  if (discret_->GetCondition("ElchBoundaryKineticsPoint") != NULL)
    EvaluateElectrodeBoundaryKineticsPointConditions(systemmatrix, rhs);

  // call base class routine
  ScaTraTimIntImpl::EvaluateSolutionDependingConditions(systemmatrix, rhs);
}  // ScaTraTimIntElch::EvaluateSolutionDependingConditions


/*----------------------------------------------------------------------*
 | check for zero/negative concentration values               gjb 01/10 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::CheckConcentrationValues(Teuchos::RCP<Epetra_Vector> vec)
{
  // action only for ELCH applications

  // for NURBS discretizations we skip the following check.
  // Control points (i.e., the "nodes" and its associated dofs can be located
  // outside the domain of interest. Thus, they can have negative
  // concentration values although the concentration solution is positive
  // in the whole computational domain!
  if (dynamic_cast<DRT::NURBS::NurbsDiscretization*>(discret_.get()) != NULL) return;

  // this option can be helpful in some rare situations
  bool makepositive(false);

  std::vector<int> numfound(NumScal(), 0);
#if 0
  std::stringstream myerrormessage;
#endif
  for (int i = 0; i < discret_->NumMyRowNodes(); i++)
  {
    DRT::Node* lnode = discret_->lRowNode(i);
    std::vector<int> dofs;
    dofs = discret_->Dof(0, lnode);

    for (int k = 0; k < NumScal(); k++)
    {
      const int lid = discret_->DofRowMap()->LID(dofs[k]);
      if (((*vec)[lid]) < EPS13)
      {
        numfound[k]++;
        if (makepositive) ((*vec)[lid]) = EPS13;
#if 0
        myerrormessage<<"PROC "<<myrank_<<" dof index: "<<k<<setprecision(7)<<scientific<<
            " val: "<<((*vec)[lid])<<" node gid: "<<lnode->Id()<<
            " coord: [x] "<< lnode->X()[0]<<" [y] "<< lnode->X()[1]<<" [z] "<< lnode->X()[2]<<std::endl;
#endif
      }
    }
  }

  // print warning to screen
  for (int k = 0; k < NumScal(); k++)
  {
    if (numfound[k] > 0)
    {
      std::cout << "WARNING: PROC " << myrank_ << " has " << numfound[k]
                << " nodes with zero/neg. concentration values for species " << k;
      if (makepositive)
        std::cout << "-> were made positive (set to 1.0e-13)" << std::endl;
      else
        std::cout << std::endl;
    }
  }

#if 0
  // print detailed info to error file
  for(int p=0; p < discret_->Comm().NumProc(); p++)
  {
    if (p==myrank_) // is it my turn?
    {
      // finish error message
      myerrormessage.flush();

      // write info to error file
      if ((errfile_!=NULL) and (myerrormessage.str()!=""))
      {
        fprintf(errfile_,myerrormessage.str().c_str());
        // std::cout<<myerrormessage.str()<<std::endl;
      }
    }
    // give time to finish writing to file before going to next proc ?
    discret_->Comm().Barrier();
  }
#endif

  // so much code for a simple check!
  return;
}  // ScaTraTimIntImpl::CheckConcentrationValues


/*----------------------------------------------------------------------*
 | apply Dirichlet boundary conditions to state vectors      fang 01/17 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::ApplyDirichletBC(const double time,  //!< time
    Teuchos::RCP<Epetra_Vector> phinp,  //!< state vector of scalar transport degrees of freedom
    Teuchos::RCP<Epetra_Vector>
        phidt  //!< state vector of time derivatives of scalar transport degrees of freedom
)
{
  // call base class routine
  ScaTraTimIntImpl::ApplyDirichletBC(time, phinp, phidt);

  // evaluate Dirichlet boundary condition on electric potential arising from constant-current
  // constant-voltage (CCCV) cell cycling boundary condition during constant-voltage (CV) phase
  if (cccv_condition_ != Teuchos::null)
    if (cccv_condition_->GetCCCVHalfCyclePhase() ==
        INPAR::ELCH::CCCVHalfCyclePhase::constant_voltage)
    {
      // initialize set for global IDs of electric potential degrees of freedom affected by
      // constant-current constant-voltage (CCCV) cell cycling boundary condition
      std::set<int> dbcgids;

      // extract constant-current constant-voltage (CCCV) half-cycle boundary conditions
      std::vector<DRT::Condition*> cccvhalfcycleconditions;
      discret_->GetCondition("CCCVHalfCycle", cccvhalfcycleconditions);

      // loop over all conditions
      for (const auto& cccvhalfcyclecondition : cccvhalfcycleconditions)
      {
        // check relevance of current condition
        if (cccvhalfcyclecondition->GetInt("ConditionID") ==
            cccv_condition_->GetHalfCycleConditionID())
        {
          // extract cutoff voltage from condition and perform safety check
          const double cutoff_voltage = cccvhalfcyclecondition->GetDouble("CutoffVoltage");
          if (cutoff_voltage < 0.)
            dserror(
                "Cutoff voltage for constant-current constant-voltage (CCCV) cell cycling must "
                "not be negative!");

          // extract nodal cloud of current condition and perform safety check
          const std::vector<int>* nodegids = cccvhalfcyclecondition->Nodes();
          if (!nodegids or !nodegids->size())
            dserror(
                "Constant-current constant-voltage (CCCV) cell cycling boundary condition does "
                "not have a nodal cloud!");

          // loop over all nodes
          for (unsigned inode = 0; inode < nodegids->size(); ++inode)
          {
            // extract global ID of current node
            const int nodegid = (*nodegids)[inode];

            // consider only nodes stored by current processor
            if (discret_->HaveGlobalNode(nodegid))
            {
              // extract current node
              const DRT::Node* const node = discret_->gNode(nodegid);

              // consider only nodes owned by current processor
              if (node->Owner() == discret_->Comm().MyPID())
              {
                // extract global ID of electric potential degree of freedom carried by current
                // node
                const int gid = discret_->Dof(
                    0, node, 1);  // Do not remove the zero, i.e., the first function argument,
                                  // otherwise an error is thrown in debug mode!

                // add global ID to set
                dbcgids.insert(gid);

                // apply cutoff voltage as Dirichlet boundary condition
                phinp->ReplaceGlobalValue(gid, 0, cutoff_voltage);
              }
            }
          }  // loop over all nodes

          // leave loop after relevant condition has been processed
          break;
        }  // relevant condition
      }    // loop over all conditions

      // transform set into vector and then into Epetra map
      std::vector<int> dbcgidsvec(dbcgids.begin(), dbcgids.end());
      const Teuchos::RCP<const Epetra_Map> dbcmap = Teuchos::rcp(new Epetra_Map(
          -1, dbcgids.size(), &dbcgidsvec[0], DofRowMap()->IndexBase(), DofRowMap()->Comm()));

      // merge map with existing map for Dirichlet boundary conditions
      AddDirichCond(dbcmap);
    }
}  // SCATRA::ScaTraTimIntElch::ApplyDirichletBC



/*----------------------------------------------------------------------*
 | evaluate Neumann boundary conditions                      fang 01/17 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::ApplyNeumannBC(
    const Teuchos::RCP<Epetra_Vector>& neumann_loads  //!< Neumann loads
)
{
  // call base class routine
  ScaTraTimIntImpl::ApplyNeumannBC(neumann_loads);

  // evaluate Neumann boundary condition on electric potential arising from constant-current
  // constant-voltage (CCCV) cell cycling boundary condition during constant-current (CC) phase
  if (cccv_condition_ != Teuchos::null)
    if (cccv_condition_->GetCCCVHalfCyclePhase() ==
        INPAR::ELCH::CCCVHalfCyclePhase::constant_current)
    {
      // extract constant-current constant-voltage (CCCV) half-cycle boundary conditions
      std::vector<DRT::Condition*> cccvhalfcycleconditions;
      discret_->GetCondition("CCCVHalfCycle", cccvhalfcycleconditions);

      // loop over all conditions
      for (const auto& cccvhalfcyclecondition : cccvhalfcycleconditions)
      {
        // check relevance of current condition
        if (cccvhalfcyclecondition->GetInt("ConditionID") ==
            cccv_condition_->GetHalfCycleConditionID())
        {
          // extract condition
          DRT::Condition& condition = *cccvhalfcyclecondition;

          // To avoid code redundancy, we evaluate the condition using the element-based algorithm
          // for standard Neumann boundary conditions. For this purpose, we must provide the
          // condition with some features to make it look like a standard Neumann boundary
          // condition.
          std::vector<int> onoff(2, 0);
          std::vector<double> val(2, 0.);
          onoff[1] =
              1;  // activate Neumann boundary condition for electric potential degree of freedom
          val[1] = condition.GetDouble("Current");  // set value of Neumann boundary condition
          condition.Add("numdof", 2);
          condition.Add("funct", std::vector<int>(2, 0));
          condition.Add("onoff", onoff);
          condition.Add("val", val);

          // create parameter list for elements
          Teuchos::ParameterList params;

          // set action for elements
          params.set<int>("action", SCATRA::bd_calc_Neumann);

          // number of dofset associated with displacement-related dofs
          if (isale_) params.set<int>("ndsdisp", nds_disp_);

          // loop over all conditioned elements
          std::map<int, Teuchos::RCP<DRT::Element>>& geometry = condition.Geometry();
          std::map<int, Teuchos::RCP<DRT::Element>>::iterator iterator;
          for (iterator = geometry.begin(); iterator != geometry.end(); ++iterator)
          {
            // get location vector of current element
            std::vector<int> lm;
            std::vector<int> lmowner;
            std::vector<int> lmstride;
            iterator->second->LocationVector(*discret_, lm, lmowner, lmstride);

            // initialize element-based vector of Neumann loads
            Epetra_SerialDenseVector elevector(lm.size());

            // evaluate Neumann boundary condition
            iterator->second->EvaluateNeumann(params, *discret_, condition, lm, elevector);

            // assemble element-based vector of Neumann loads into global vector of Neumann loads
            LINALG::Assemble(*neumann_loads, elevector, lm, lmowner);
          }  // loop over all conditioned elements

          // leave loop after relevant condition has been processed
          break;
        }  // relevant condition
      }    // loop over all conditions
    }
}  // SCATRA::ScaTraTimIntElch::ApplyNeumannBC



/*---------------------------------------------------------------------------*
 | determine whether there are still time steps to be evaluated   fang 01/17 |
 *---------------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntElch::NotFinished()
{
  // call base class routine in case no cell cycling is performed
  if (cccv_condition_ == Teuchos::null) return ScaTraTimIntImpl::NotFinished();

  // control progress of simulation in case cell cycling is performed
  // note that the maximum number of time steps and the maximum simulation time are ignored in
  // this case
  else
  {
    // which mode was last converged step? Is this phase over? Is the current half cycle over?
    if (cccv_condition_->GetCCCVHalfCyclePhase() ==
        INPAR::ELCH::CCCVHalfCyclePhase::initital_relaxation)
    {
      if (cccv_condition_->IsInitialRelaxation(time_, Dt()))
      {
        // do nothing
      }
      else
        cccv_condition_->SetFirstCCCVHalfCycle(step_);
      return true;
    }
    else
      while (cccv_condition_->IsEndOfHalfCyclePhase(cellvoltage_, cellcrate_, time_))
        cccv_condition_->NextPhase(step_, time_);

    // all half cycles completed?
    return (cccv_condition_->NotFinished());
  }
}


/*---------------------------------------------------------------------------*
 | perform Aitken relaxation                                      fang 08/17 |
 *---------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::PerformAitkenRelaxation(
    Epetra_Vector& phinp,  //!< state vector to be relaxed
    const Epetra_Vector&
        phinp_inc_diff  //!< difference between current and previous state vector increments
)
{
  if (solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit)
  {
    // safety checks
    if (splitter_macro_ == Teuchos::null)
      dserror("Map extractor for macro scale has not been initialized yet!");

    // loop over all degrees of freedom
    for (int idof = 0; idof < splitter_macro_->NumMaps(); ++idof)
    {
      // extract subvectors associated with current degree of freedom
      const Teuchos::RCP<const Epetra_Vector> phinp_inc_dof =
          splitter_macro_->ExtractVector(*phinp_inc_, idof);
      const Teuchos::RCP<const Epetra_Vector> phinp_inc_diff_dof =
          splitter_macro_->ExtractVector(phinp_inc_diff, idof);

      // compute L2 norm of difference between current and previous increments of current degree
      // of freedom
      double phinp_inc_diff_L2(0.);
      phinp_inc_diff_dof->Norm2(&phinp_inc_diff_L2);

      // compute dot product between increment of current degree of freedom and difference between
      // current and previous increments of current degree of freedom
      double phinp_inc_dot_phinp_inc_diff(0.);
      if (phinp_inc_diff_dof->Dot(*phinp_inc_dof, &phinp_inc_dot_phinp_inc_diff))
        dserror("Couldn't compute dot product!");

      // compute Aitken relaxation factor for current degree of freedom
      if (iternum_outer_ > 1 and phinp_inc_diff_L2 > 1.e-12)
        omega_[idof] *= 1 - phinp_inc_dot_phinp_inc_diff / (phinp_inc_diff_L2 * phinp_inc_diff_L2);

      // perform Aitken relaxation for current degree of freedom
      splitter_macro_->AddVector(*phinp_inc_dof, idof, phinp, omega_[idof]);
    }
  }

  else
    // call base class routine
    ScaTraTimIntImpl::PerformAitkenRelaxation(phinp, phinp_inc_diff);
}


/*----------------------------------------------------------------------*
 | output flux vectors                                       fang 08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::OutputFlux(Teuchos::RCP<Epetra_MultiVector> flux,  //!< flux vector
    const std::string& fluxtype  //!< flux type ("domain" or "boundary")
)
{
  // safety check
  if (flux == Teuchos::null) dserror("Invalid flux vector!");

  if (fluxtype == "domain")
  {
    // In this case, flux output can be straightforwardly performed without additional
    // manipulation.
  }

  else if (fluxtype == "boundary")
  {
    // The closing equation for the electric potential is internally scaled by the factor 1/F for
    // better conditioning. Therefore, the associated boundary flux computed by the function
    // CalcFluxAtBoundary is also scaled by this factor. To avoid confusion, we remove the scaling
    // factor from the boundary flux before outputting it, so that the result can be physically
    // interpreted as the plain boundary current density without any scaling.
    splitter_->Scale(*flux, 1, elchparams_->get<double>("FARADAY_CONSTANT"));
  }

  else
    dserror("Unknown flux type! Must be either 'domain' or 'boundary'!");

  // perform actual flux output by calling base class routine
  ScaTraTimIntImpl::OutputFlux(flux, fluxtype);
}  // SCATRA::ScaTraTimIntElch::OutputFlux


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SCATRA::ScalarHandlerElch::ScalarHandlerElch() : numscal_() { return; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::ScalarHandlerElch::Setup(const ScaTraTimIntImpl* const scatratimint)
{
  // call base class
  ScalarHandler::Setup(scatratimint);

  // cast to electrochemistry time integrator
  const ScaTraTimIntElch* const elchtimint =
      dynamic_cast<const ScaTraTimIntElch* const>(scatratimint);
  if (elchtimint == NULL) dserror("cast to ScaTraTimIntElch failed!");

  // adapt number of transported scalars if necessary
  // current is a solution variable
  if (DRT::INPUT::IntegralValue<int>(
          elchtimint->ElchParameterList()->sublist("DIFFCOND"), "CURRENT_SOLUTION_VAR"))
  {
    // shape of local row element(0) -> number of space dimensions
    // int dim = problem_->NDim();
    int dim = DRT::UTILS::getDimension(elchtimint->Discretization()->lRowElement(0)->Shape());
    // number of concentrations transported is numdof-1-nsd
    numscal_.clear();
    numscal_.insert(NumDofPerNode() - 1 - dim);
  }

  // multi-scale case
  else if (elchtimint->MacroScale())
  {
    // number of transported scalars is 1
    numscal_.clear();
    numscal_.insert(1);
  }

  // standard case
  else
  {
    // number of transported scalars is numdof-1 (last dof = electric potential)
    numscal_.clear();
    numscal_.insert(NumDofPerNode() - 1);
  }
}

/*-------------------------------------------------------------------------*
|  Determine number of Scalars per node in given condition    vuong   04/16|
 *-------------------------------------------------------------------------*/
int SCATRA::ScalarHandlerElch::NumScalInCondition(
    const DRT::Condition& condition, const Teuchos::RCP<const DRT::Discretization>& discret) const
{
  CheckIsSetup();
  // for now only equal dof numbers are supported
  if (not equalnumdof_)
    dserror(
        "Different number of DOFs per node within ScaTra discretization! This is not supported for "
        "Elch!");

  return NumScal();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::BuildBlockMaps(
    const std::vector<Teuchos::RCP<DRT::Condition>>& partitioningconditions,
    std::vector<Teuchos::RCP<const Epetra_Map>>& blockmaps) const
{
  if (MatrixType() == INPAR::SCATRA::MatrixType::block_condition_dof)
  {
    // safety check
    if (DRT::INPUT::IntegralValue<int>(
            ElchParameterList()->sublist("DIFFCOND"), "CURRENT_SOLUTION_VAR"))
      dserror(
          "For chosen type of global block system matrix, current must not constitute solution "
          "variable!");

    // extract number of domain partitioning conditions
    const unsigned ncond = partitioningconditions.size();

    // each domain block specified by the domain partitioning conditions is again subdivided into
    // the specific dofs, therefore (ncond * NumDofPerNode()) block maps are set up below
    blockmaps.resize(ncond * NumDofPerNode(), Teuchos::null);

    // loop over all domain partitioning conditions
    for (unsigned icond = 0; icond < ncond; ++icond)
    {
      // we need to initialize as many sets as number of dofs per node, since all ids
      // corresponding to a specific dof shall be grouped into a separate set
      std::vector<std::set<int>> dofids(NumDofPerNode());

      // extract nodes associated with current domain partitioning condition
      const std::vector<int>* nodegids = partitioningconditions[icond]->Nodes();

      // loop over all nodes associated with current domain partitioning condition
      for (int nodegid : *nodegids)
      {
        // extract global ID of current node
        // consider current node only if node is owned by current processor need to make sure that
        // node is stored on current processor, otherwise cannot resolve "->Owner()"
        if (discret_->HaveGlobalNode(nodegid) and
            discret_->gNode(nodegid)->Owner() == discret_->Comm().MyPID())
        {
          // extract dof IDs associated with current node
          const std::vector<int> nodedofs = discret_->Dof(0, discret_->gNode(nodegid));

          // for each dof we add the id of the current node to its corresponding set
          for (unsigned dof = 0; dof < dofids.size(); ++dof) dofids[dof].insert(nodedofs[dof]);
        }
      }

      // transform sets for dof IDs into vectors and then into Epetra maps
      for (int iset = 0; iset < NumDofPerNode(); ++iset)
      {
        int nummyelements(0);
        int* myglobalelements(nullptr);
        std::vector<int> dofidvec;
        if (!dofids[iset].empty())
        {
          dofidvec.reserve(dofids[iset].size());
          dofidvec.assign(dofids[iset].begin(), dofids[iset].end());
          nummyelements = dofidvec.size();
          myglobalelements = &(dofidvec[0]);
        }
        blockmaps[NumDofPerNode() * icond + iset] = Teuchos::rcp(new Epetra_Map(-1, nummyelements,
            myglobalelements, discret_->DofRowMap()->IndexBase(), discret_->DofRowMap()->Comm()));
      }
    }
  }

  // call base class routine for other types of global system matrix
  else
    SCATRA::ScaTraTimIntImpl::BuildBlockMaps(partitioningconditions, blockmaps);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntElch::BuildBlockNullSpaces() const
{
  // call base class routine
  SCATRA::ScaTraTimIntImpl::BuildBlockNullSpaces();

  if (MatrixType() == INPAR::SCATRA::MatrixType::block_condition_dof)
  {
    // loop over blocks of global system matrix
    for (int iblock = 0; iblock < BlockMaps().NumMaps(); ++iblock)
    {
      // store number of current block as string, starting from 1
      std::ostringstream iblockstr;
      iblockstr << iblock + 1;

      // access parameter sublist associated with smoother for current matrix block
      Teuchos::ParameterList& mueluparams =
          Solver()->Params().sublist("Inverse" + iblockstr.str()).sublist("MueLu Parameters");

      // extract already reduced null space associated with current matrix block
      std::vector<double>& nullspace =
          *mueluparams.get<Teuchos::RCP<std::vector<double>>>("nullspace");

      // Each matrix block is associated with either concentration dofs or electric potential dofs
      // only. However, since the original full null space was computed for all degrees of freedom
      // on the discretization, the reduced null spaces still have the full dimension, i.e., the
      // full number of null space vectors equaling the total number of primary variables. Hence,
      // we need to decrease the dimension of each null space by one and remove the corresponding
      // zero null space vector from the null space.
      if (iblock % 2 == 0)
        // null space associated with concentration dofs
        // remove zero null space vector associated with electric potential dofs by truncating
        // null space
        nullspace.resize(BlockMaps().Map(iblock)->NumMyElements());

      else
        // null space associated with electric potential dofs
        // remove zero null space vector(s) associated with concentration dofs and retain only the
        // last null space vector associated with electric potential dofs
        nullspace.erase(
            nullspace.begin(), nullspace.end() - BlockMaps().Map(iblock)->NumMyElements());

      // decrease null space dimension and number of partial differential equations by one
      --mueluparams.get<int>("null space: dimension");
      --mueluparams.get<int>("PDE equations");
    }
  }
}
