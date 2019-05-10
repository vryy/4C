/*!----------------------------------------------------------------------

\brief One-Step-Theta time-integration scheme for endo-/exocytosis model

\level 2

\maintainer Jonas Eichinger

*----------------------------------------------------------------------*/

#include "scatra_timint_ost_endoexocytosis.H"
#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_cell.H"
#include "../drt_io/io.H"

/*----------------------------------------------------------------------*
 |  Constructor (public)                                    rauch 08/16 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntOneStepThetaEndoExocytosis::TimIntOneStepThetaEndoExocytosis(
    Teuchos::RCP<DRT::Discretization> actdis,          //!< discretization
    Teuchos::RCP<LINALG::Solver> solver,               //!< linear solver
    Teuchos::RCP<Teuchos::ParameterList> params,       //!< parameter list
    Teuchos::RCP<Teuchos::ParameterList> extraparams,  //!< supplementary parameter list
    Teuchos::RCP<IO::DiscretizationWriter> output,     //!< output writer
    const int probnum                                  //!< global problem number
    )
    : ScaTraTimIntImpl(actdis, solver, params, extraparams, output),
      TimIntOneStepTheta(actdis, solver, params, extraparams, output),
      endoexo_delay_steps_(-1),
      exo_domain_(0),
      exo_dof_(0),
      endo_theta_(-1.0),
      exo_domain_integral_value_(0.0),
      Delta_phi_(Teuchos::null),
      total_scalars_(0, 0.0),
      internalization_vec_(Teuchos::null),
      initial_exocytosis_(0.0),
      source_(Teuchos::null)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized

  return;
}


/*------------------------------------------------------------------------*
 |  initialize time integration                               rauch 08/16 |
 *------------------------------------------------------------------------*/
void SCATRA::TimIntOneStepThetaEndoExocytosis::Setup()
{
  // initialize base class
  TimIntOneStepTheta::Setup();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // create exocytotic source vector
  source_ = LINALG::CreateVector(*dofrowmap, true);

  // get delay between internalization and begin of exocytotic transport at exo. surface
  endoexo_delay_steps_ = DRT::Problem::Instance()
                             ->CellMigrationParams()
                             .sublist("ENDOEXOCYTOSIS MODULE")
                             .get<int>("ENDOEXO_DELAY");

  // initialize internalization vector
  internalization_vec_ = Teuchos::rcp(new std::vector<double>);
  internalization_vec_->resize(endoexo_delay_steps_, 0.0);


  // initialize vector in which the deltas (n -> n+1) are stored
  // entry 0 contains the difference at t-T
  Delta_phi_ = Teuchos::rcp(new std::vector<double>);
  Delta_phi_->resize(endoexo_delay_steps_, 0.0);

  exo_domain_ = DRT::INPUT::IntegralValue<int>(
      DRT::Problem::Instance()->CellMigrationParams().sublist("ENDOEXOCYTOSIS MODULE"),
      "EXODOMAIN");

  exo_dof_ = DRT::Problem::Instance()
                 ->CellMigrationParams()
                 .sublist("ENDOEXOCYTOSIS MODULE")
                 .get<int>("EXO_DOF");

  endo_theta_ = DRT::Problem::Instance()
                    ->CellMigrationParams()
                    .sublist("ENDOEXOCYTOSIS MODULE")
                    .get<double>("THETA");

  return;
}


/*------------------------------------------------------------------------*
 |  manipulations of scatra algorithm prior to Solve()        rauch 08/16 |
 *------------------------------------------------------------------------*/
void SCATRA::TimIntOneStepThetaEndoExocytosis::PreSolve()
{
  CheckIsInit();
  CheckIsSetup();

  // temp variables
  const int dof_id_internalized_scalar =
      discret_->GetCondition("ScaTraCellExt")->GetInt("ScalarID");

  double delta_phi = -1234.0;  //!< concentration difference of internalized scalar

  // calculate the new difference in molecule numbers from n -> n+1
  if (step_ > endoexo_delay_steps_)
  {
    // Member step_ is incremented before Solve. Thus, it starts with 1 for the first time step.
    int entry = ((step_ - 1) - endoexo_delay_steps_) % (endoexo_delay_steps_);
    delta_phi = Delta_phi_->at(entry);
  }
  else
    delta_phi = initial_exocytosis_;

  /////////////////////////////////////////////////////////////////////////////////////
  // Calc area of exocytotic domain (surface or volume).
  // The domain needs to have the condition 'ScaTraCellExt'.
  // In this domain the biomolecules internalized at t-endoexo_delay_
  // are segregated.

  // declare parameter list
  Teuchos::ParameterList eleparams;

  // reinitialize surface area
  exo_domain_integral_value_ = 0.0;

  // set action for surface area calculation
  if (exo_domain_ == INPAR::CELL::exo_surface)
    eleparams.set<int>("action", SCATRA::bd_calc_boundary_integral);
  else if (exo_domain_ == INPAR::CELL::exo_volume)
    eleparams.set<int>("action", SCATRA::calc_domain_integral);
  else
    dserror("Either surface or volume exocytosis has to be chosen.");


  // create result vector
  Teuchos::RCP<Epetra_SerialDenseVector> domain = Teuchos::rcp(new Epetra_SerialDenseVector(1));

  // evaluate over surface "ScaTraCellExt"
  discret_->EvaluateScalars(eleparams, domain, "ScaTraCellExt");

  // extract domain from result vector
  exo_domain_integral_value_ = (*domain)[0];

  // safety check
  if (std::abs(exo_domain_integral_value_) < 1E-15)
    dserror(
        "Exocytosis surface area/volume is close to zero!\n"
        "Something went wrong.\n"
        "Check definition of 'ScaTraCellExt' conditioned surface.");

  /////////////////////////////////////////////////////////////////////////////////////
  // Calc exocytotic source term at externalization entity (surface or volume).
  // The entity needs to have the condition 'ScaTraCellExt'.
  // At this entity the biomolecules internalized at t-endoexo_delay_
  // are segregated.

  // reinitialize source vector
  source_->Scale(-(1.0 / endo_theta_) + 1.0);

  // set action for source evaluation
  if (exo_domain_ == INPAR::CELL::exo_surface)
    eleparams.set<int>("action", SCATRA::bd_integrate_weighted_scalar);
  else if (exo_domain_ == INPAR::CELL::exo_volume)
    eleparams.set<int>("action", SCATRA::integrate_weighted_scalar);
  else
    dserror("Either surface or volume exocytosis has to be chosen.");

  // provide other data to evaluate routine
  eleparams.set<double>("scalar", (delta_phi / (endo_theta_ * dta_)));
  eleparams.set<double>("user defined prefac", -(1.0 / exo_domain_integral_value_));
  eleparams.set<int>("ScalarID", exo_dof_);

  // evaluate the source
  discret_->EvaluateCondition(eleparams, source_, "ScaTraCellExt");

  double sourcenorm = -1234.0;
  source_->Norm2(&sourcenorm);


  if (discret_->Comm().MyPID() == 0)
  {
    std::cout << "########################################################" << std::endl;
    std::cout << "# Externalized " << std::setprecision(7) << delta_phi << " molecules of species "
              << dof_id_internalized_scalar << std::endl;
    std::cout << "# L2-Norm of source vector = " << std::setprecision(7) << sourcenorm << std::endl;
    std::cout << "########################################################" << std::endl;
  }

  // write source_ to rhs
  GetNeumannLoadsPtr()->Update(1.0, *source_, 1.0);

  return;
}


/*------------------------------------------------------------------------*
 |  manipulations of scatra algorithm after call to Solve()   rauch 08/16 |
 *------------------------------------------------------------------------*/
void SCATRA::TimIntOneStepThetaEndoExocytosis::PostSolve()
{
  // temp variables
  const int dof_id_internalized_scalar =
      discret_->GetCondition("ScaTraCellExt")->GetInt("ScalarID");
  const int numscal = NumScalInCondition(*(discret_->GetCondition("ScaTraCellInt")));

  /////////////////////////////////////////////////////////////////////////////////////
  // clear state and set new state phinp
  discret_->ClearState();
  discret_->SetState("phinp", phinp_);

  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", SCATRA::calc_total_and_mean_scalars);
  eleparams.set<bool>("inverting", false);

  // provide displacement field in case of ALE
  if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

  // evaluate integrals of concentrations and domain
  Teuchos::RCP<Epetra_SerialDenseVector> scalars =
      Teuchos::rcp(new Epetra_SerialDenseVector(numscal + 1));
  discret_->EvaluateScalars(eleparams, scalars, "ScaTraCellInt");

  // clear all states
  discret_->ClearState();


  /////////////////////////////////////////////////////////////////////////////////////
  // calculate mean concentrations
  total_scalars_.resize(numscal);
  if (std::abs((*scalars)[numscal]) < 1E-15)
    dserror(
        "Domain has zero volume!\n"
        "Something went wrong.\n"
        "Check ---DESIGN TOTAL AND MEAN SCALAR VOL CONDITIONS in your .dat file!");

  // loop over all scalars and store mean scalars
  for (int k = 0; k < numscal; ++k) total_scalars_[k] = (*scalars)[k];

  // entry associated with current time step 'step_'
  int entry = -1234;
  if (step_ > endoexo_delay_steps_)
  {
    // Member step_ is incremented before Solve. Thus, it starts with 1 for the first time step.
    entry = ((step_ - 1) - endoexo_delay_steps_) % (endoexo_delay_steps_);
  }
  else
  {
    // Member step_ is incremented before Solve. Thus, it starts with 1 for the first time step.
    entry = ((step_ - 1)) % (endoexo_delay_steps_);
  }

  // in case we write to the first vector entry, the previous value is the last vector entry
  double previousvectorentry = -1234.0;
  if (entry == 0)
    previousvectorentry = internalization_vec_->at(endoexo_delay_steps_ - 1);
  else
    previousvectorentry = internalization_vec_->at(entry - 1);

  // save current amount of internalized scalar at correct position in vector belonging to
  // time-delay
  internalization_vec_->at(entry) = total_scalars_[dof_id_internalized_scalar];

  // Do in first time step
  if (step_ == 1)
  {
    // The first delta equals the amount of internalized biomolecules, since
    // initially we assume zero internalized biomolcules.
    Delta_phi_->at(0) = internalization_vec_->at(0);
  }
  else
  {
    // calc the difference between the internalized biomolecules in the current and the previous
    // step
    Delta_phi_->at(entry) = (internalization_vec_->at(entry) - previousvectorentry);
  }

  if (discret_->Comm().MyPID() == 0)
  {
    std::cout << "########################################################" << std::endl;
    std::cout << "# Internalized " << std::setprecision(7) << Delta_phi_->at(entry)
              << " molecules of membrane species " << dof_id_internalized_scalar << std::endl;
    std::cout << "########################################################" << std::endl;
  }

  return;
}


/*----------------------------------------------------------------------*
 | write additional data required for restart               rauch 09/17 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepThetaEndoExocytosis::OutputRestart() const
{
  // call base class routine
  TimIntOneStepTheta::OutputRestart();

  // additional state vectors that are needed for endo-/exocytosis restart
  output_->WriteVector("exo_source", source_);
  output_->WriteRedundantDoubleVector("Delta_phi", Delta_phi_);
  output_->WriteRedundantDoubleVector("internalization_vec", internalization_vec_);

  return;
}


/*----------------------------------------------------------------------*
 |                                                          rauch 09/17 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntOneStepThetaEndoExocytosis::ReadRestart(
    const int step, Teuchos::RCP<IO::InputControl> input)
{
  // call base class routine
  TimIntOneStepTheta::ReadRestart(step, input);

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  if (myrank_ == 0)
    std::cout << "Reading Endo-/Exocytosis restart data (time=" << time_ << " ; step=" << step_
              << ")" << std::endl;

  // read state vectors that are needed for endo-/exocytosis restart
  reader->ReadVector(source_, "exo_source");
  reader->ReadRedundantDoubleVector(Delta_phi_, "Delta_phi");
  reader->ReadRedundantDoubleVector(internalization_vec_, "internalization_vec");

  return;
}


/*----------------------------------------------------------------------*
 | Destructor dtor                                 (public) rauch 08/16 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntOneStepThetaEndoExocytosis::~TimIntOneStepThetaEndoExocytosis() { return; }
