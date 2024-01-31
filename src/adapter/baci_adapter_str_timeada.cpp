/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration which gives fine grained
       access in the adaptive time marching loop


\level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_adapter_str_timeada.H"

#include "baci_adapter_str_timeada_joint.H"
#include "baci_adapter_str_timeada_zienxie.H"
#include "baci_global_data.H"
#include "baci_io.H"
#include "baci_io_pstream.H"
#include "baci_lib_discret.H"
#include "baci_linalg_utils_sparse_algebra_create.H"
#include "baci_structure_new_dbc.H"
#include "baci_structure_new_timint_base.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureTimeAda::StructureTimeAda(Teuchos::RCP<Structure> structure)
    : StructureWrapper(structure)
{
  stm_ = Teuchos::rcp_dynamic_cast<STR::TIMINT::Base>(structure_);

  if (stm_ == Teuchos::null) dserror("cast from ADAPTER::Structure to STR::TIMINT::Base failed");

  // call the setup once if stm_ has been setup
  if (stm_->IsSetup()) SetupTimeAda();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::Structure> ADAPTER::StructureTimeAda::Create(
    const Teuchos::ParameterList& taflags,  //!< adaptive input flags
    Teuchos::RCP<STR::TIMINT::Base> ti_strategy)
{
  auto kind = INPUT::IntegralValue<INPAR::STR::TimAdaKind>(taflags, "KIND");
  switch (kind)
  {
    case INPAR::STR::timada_kind_zienxie:
      // Adaptive time integration with Zienkiewicz-Xie error indicator
      return Teuchos::rcp(new ADAPTER::StructureTimeAdaZienXie(ti_strategy));

    case INPAR::STR::timada_kind_joint_explicit:
      // Adaptive time integration using auxiliary time integrator
      return Teuchos::rcp(new ADAPTER::StructureTimeAdaJoint(ti_strategy));

    default:
      // Unknown adaptive time integration
      return Teuchos::null;
  }

  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimeAda::Setup()
{
  // call the wrapper setup
  StructureWrapper::Setup();

  // self setup
  SetupTimeAda();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimeAda::SetupTimeAda()
{
  const Teuchos::ParameterList& sdynparams = GLOBAL::Problem::Instance()->StructuralDynamicParams();

  // initialize the local variables
  timeinitial_ = 0.0;
  timefinal_ = sdynparams.get<double>("MAXTIME");
  if (timefinal_ <= timeinitial_) dserror("MAXTIME is not positive. It is invalid.");
  timedirect_ = 1.0;
  timestepinitial_ = 0;
  timestepfinal_ = sdynparams.get<int>("NUMSTEP");
  stepsizeinitial_ = sdynparams.get<double>("TIMESTEP");

  const Teuchos::ParameterList& tap = sdynparams.sublist("TIMEADAPTIVITY");

  stepsizemax_ = tap.get<double>("STEPSIZEMAX");
  stepsizemin_ = tap.get<double>("STEPSIZEMIN");
  sizeratiomax_ = tap.get<double>("SIZERATIOMAX");
  sizeratiomin_ = tap.get<double>("SIZERATIOMIN");
  sizeratioscale_ = tap.get<double>("SIZERATIOSCALE");
  errctrl_ = ctrl_dis;  // PROVIDE INPUT PARAMETER
  errnorm_ = INPUT::IntegralValue<INPAR::STR::VectorNorm>(tap, "LOCERRNORM");
  errtol_ = tap.get<double>("LOCERRTOL");
  errorder_ = 1;  // CHANGE THIS CONSTANT
  adaptstepmax_ = tap.get<int>("ADAPTSTEPMAX");

  time_ = timeinitial_;
  timestep_ = 0;
  stepsizepre_ = stepsizeinitial_;
  stepsize_ = sdynparams.get<double>("TIMESTEP");
  adaptstep_ = 0;

  outsys_ = false;
  outstr_ = false;
  outene_ = false;
  outrest_ = false;
  outsysperiod_ = tap.get<double>("OUTSYSPERIOD");
  outstrperiod_ = tap.get<double>("OUTSTRPERIOD");
  outeneperiod_ = tap.get<double>("OUTENEPERIOD");
  outrestperiod_ = tap.get<double>("OUTRESTPERIOD");
  outsizeevery_ = tap.get<int>("OUTSIZEEVERY");
  outsystime_ = timeinitial_ + outsysperiod_;
  outstrtime_ = timeinitial_ + outstrperiod_;
  outenetime_ = timeinitial_ + outeneperiod_;
  outresttime_ = timeinitial_ + outrestperiod_;
  outsizefile_ = Teuchos::null;

  // allocate displacement local error vector
  locerrdisn_ = CORE::LINALG::CreateVector(*(stm_->DofRowMap()), true);

  // enable restart for adaptive timestepping
  const int restart = GLOBAL::Problem::Instance()->Restart();
  if (restart)
  {
    // read restart of marching time-integrator and reset initial time and step for adaptive loop
    stm_->ReadRestart(restart);
    timeinitial_ = stm_->TimeOld();
    timestepinitial_ = stm_->StepOld();
    IO::DiscretizationReader ioreader(stm_->Discretization(), timestepinitial_);
    stepsizepre_ = ioreader.ReadDouble("next_delta_time");
    time_ = timeinitial_;

    // update variables which depend on initial time and step
    timedirect_ = timefinal_ > timeinitial_ ? 1.0 : -1.0;
    outsystime_ = timeinitial_ + outsysperiod_;
    outstrtime_ = timeinitial_ + outstrperiod_;
    outenetime_ = timeinitial_ + outeneperiod_;
    outresttime_ = timeinitial_ + outrestperiod_;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureTimeAda::ReadRestart(int step)
{
  SetupTimeAda();
  SetupAuxiliar();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::StructureTimeAda::Integrate()
{
  // error checking variables
  INPAR::STR::ConvergenceStatus convergencestatus = INPAR::STR::conv_success;

  int myrank = stm_->Discretization()->Comm().MyPID();

  // finalize initialization
  // (only relevant if an auxiliary time integrator is used)
  // buih STR:TIMINT::Base should be initialized outside
  // stm_->Init();

  // Richardson extrapolation to no avail
  if (MethodAdaptDis() == ada_ident)
    dserror(
        "This combination is not implemented ... Richardson's extrapolation ... Yoshida technique "
        "...");

  // initialise time loop
  time_ = timeinitial_;
  timestep_ = timestepinitial_;
  stepsize_ = stepsizepre_;

  // time loop
  while ((time_ < timefinal_) and (timestep_ < timestepfinal_))
  {
    // time step size adapting loop
    adaptstep_ = 0;
    bool accepted = false;
    double stpsiznew;
    while ((not accepted) and (adaptstep_ < adaptstepmax_))
    {
      // modify step-size #stepsize_ according to output period
      // and store output type on #outstep_
      SizeForOutput();

      // set current step size
      stm_->SetDeltaTime(stepsize_);
      stm_->SetTimeNp(time_ + stepsize_);

      // integrate system with auxiliary TIS
      // we hold \f$D_{n+1}^{AUX}\f$ on #locdiserrn_
      // and \f$V_{n+1}^{AUX}\f$ on #locvelerrn_
      IntegrateStepAuxiliar();

      // call the predictor
      PrePredict();
      PrepareTimeStep();

      // integrate system with marching TIS and
      // stm_->IntegrateStep();
      PreSolve();
      convergencestatus = Solve();

      if (convergencestatus != INPAR::STR::conv_success)
      {
        // if not converged, then we have to restart the step over
        accepted = false;

        // get the divergence action
        enum INPAR::STR::DivContAct div_action = stm_->DataSDyn().GetDivergenceAction();

        convergencestatus = PerformErrorAction(div_action, stpsiznew);
      }

      if (convergencestatus == INPAR::STR::conv_success)
      {
        // get local error vector on locerrdisn_
        EvaluateLocalErrorDis();

        // check whether step passes
        Indicate(accepted, stpsiznew);
      }

      // adjust step-size and prepare repetition of current step
      if (not accepted)
      {
        std::cout << "Repeating step " << timestep_ + 1 << "/" << timestepfinal_ << " at time "
                  << time_ << " with stepsize = " << stpsiznew << std::endl;
        std::cout << "- - - - - - - - - - - - - - - - - - - - - - - - -"
                  << " - - - - - - - - - - - - - - -" << std::endl;

        stepsize_ = stpsiznew;

        ResetStep();
      }

      // increment number of adapted step sizes in a row
      adaptstep_ += 1;
    }

    // update or break
    if (accepted)
    {
      if (myrank == 0) std::cout << "Step size accepted" << std::endl;
    }
    else if (adaptstep_ >= adaptstepmax_)
    {
      if (myrank == 0)
        std::cout << "Could not find acceptable time step size"
                  << " ... continuing" << std::endl;
    }
    else
    {
      dserror("Do not know what to do");
    }

    // calculate stresses, strains and energies
    // note: this has to be done before the update since otherwise a potential
    // material history is overwritten
    constexpr bool force_prepare = false;
    PrepareOutput(force_prepare);

    // update displacements, velocities, accelerations
    // after this call we will have disn_==dis_, etc
    // update time and step
    // update everything on the element level
    PreUpdate();

    Update();

    PostUpdate();

    stepsizepre_ = stepsize_;
    stepsize_ = stpsiznew;

    // write output
    Output();
    PostOutput();

    // print info about finished time step
    PrintStep();

    // update
    timestep_ += 1;
    time_ += stepsizepre_;
    stm_->SetStepN(timestep_);
    stm_->SetTimeN(time_);
    stm_->SetDeltaTime(stepsize_);
    stm_->SetTimeNp(time_ + stepsize_);

    UpdatePeriod();
    outrest_ = outsys_ = outstr_ = outene_ = false;

    UpdateAuxiliar();

    // the user reads but rarely listens
    if (myrank == 0)
    {
      std::cout << "Step " << timestep_ + 1 << ", Time " << time_ << ", new StepSize " << stepsize_
                << std::endl;
    }
  }

  // force write output
  Output(true);

  // that's it say what went wrong
  return convergencestatus;
}

/*----------------------------------------------------------------------*/
/*  Modify step size to hit precisely output period */
void ADAPTER::StructureTimeAda::SizeForOutput()
{
  // check output of restart data first
  if ((fabs(time_ + stepsize_) >= fabs(outresttime_)) and (outrestperiod_ != 0.0))

  {
    stepsize_ = outresttime_ - time_;
    outrest_ = true;
  }

  // check output of system vectors
  if ((fabs(time_ + stepsize_) >= fabs(outsystime_)) and (outsysperiod_ != 0.0))
  {
    stepsize_ = outsystime_ - time_;
    outsys_ = true;
    if (fabs(outsystime_) < fabs(outresttime_)) outrest_ = false;
  }

  // check output of stress/strain
  if ((fabs(time_ + stepsize_) >= fabs(outstrtime_)) and (outstrperiod_ != 0.0))
  {
    stepsize_ = outstrtime_ - time_;
    outstr_ = true;
    if (fabs(outstrtime_) < fabs(outresttime_)) outrest_ = false;
    if (fabs(outstrtime_) < fabs(outsystime_)) outsys_ = false;
  }

  // check output of energy
  if ((fabs(time_ + stepsize_) >= fabs(outenetime_)) and (outeneperiod_ != 0.0))
  {
    stepsize_ = outenetime_ - time_;
    outene_ = true;
    if (fabs(outenetime_) < fabs(outresttime_)) outrest_ = false;
    if (fabs(outenetime_) < fabs(outsystime_)) outsys_ = false;
    if (fabs(outenetime_) < fabs(outstrtime_)) outstr_ = false;
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Output action */
void ADAPTER::StructureTimeAda::Output(bool forced_writerestart)
{
  STR::TIMINT::BaseDataIO& dataio = stm_->DataIO();
  Teuchos::RCP<IO::DiscretizationWriter> output_ptr = dataio.GetOutputPtr();

  StructureWrapper::Output(forced_writerestart);
  output_ptr->WriteDouble("next_delta_time", stepsize_);
}

/*----------------------------------------------------------------------*/
/* Evaluate local error vector */
void ADAPTER::StructureTimeAda::EvaluateLocalErrorDis()
{
  const STR::TIMINT::Base& sti = *stm_;
  const auto& gstate = sti.DataGlobalState();

  if (MethodAdaptDis() == ada_orderequal)
  {
    const double coeffmarch = sti.MethodLinErrCoeffDis();
    const double coeffaux = MethodLinErrCoeffDis();
    locerrdisn_->Update(-1.0, *(gstate.GetDisNp()), 1.0);
    locerrdisn_->Scale(coeffmarch / (coeffaux - coeffmarch));
  }
  else
  {
    // schemes do not have the same order of accuracy
    locerrdisn_->Update(-1.0, *(gstate.GetDisNp()), 1.0);
  }

  // blank Dirichlet DOFs since they always carry the exact solution
  sti.GetDBC().ApplyDirichletToVector(locerrdisn_);
}

/*----------------------------------------------------------------------*/
/* Indicate error and determine new step size */
void ADAPTER::StructureTimeAda::Indicate(bool& accepted, double& stpsiznew)
{
  // norm of local discretisation error vector
  const int numneglect = stm_->GetDBCMapExtractor()->CondMap()->NumGlobalElements();
  const double norm = CalculateVectorNorm(errnorm_, locerrdisn_, numneglect);

  // check if acceptable
  accepted = (norm < errtol_);

  // debug
  int myrank = stm_->Discretization()->Comm().MyPID();
  if (myrank == 0)
  {
    std::cout << "LocErrNorm " << std::scientific << norm << ", LocErrTol " << errtol_
              << ", Accept " << std::boolalpha << accepted << std::endl;
  }

  stpsiznew = CalculateDt(norm);
}

/*----------------------------------------------------------------------*/
/* Prepare repetition of current time step */
void ADAPTER::StructureTimeAda::ResetStep()
{
  outrest_ = outsys_ = outstr_ = outene_ = false;
  // set current step size
  stm_->SetDeltaTime(stepsize_);
  stm_->SetTimeNp(time_ + stepsize_);
  // reset the integrator
  stm_->ResetStep();
}

/*----------------------------------------------------------------------*/
/* Indicate error and determine new step size */
double ADAPTER::StructureTimeAda::CalculateDt(const double norm)
{
  // get error order
  if (MethodAdaptDis() == ada_upward)
    errorder_ = stm_->MethodOrderOfAccuracyDis();
  else
    errorder_ = MethodOrderOfAccuracyDis();

  // optimal size ration with respect to given tolerance
  double sizrat = 1.0;
  if (not(norm == 0.0))  // do not divide by zero
    sizrat = std::pow(errtol_ / norm, 1.0 / (errorder_ + 1.0));
  else  // max increase if error norm == 0
    sizrat = sizeratiomax_ / sizeratioscale_;

  // debug
  int myrank = stm_->Discretization()->Comm().MyPID();
  if (myrank == 0)
  {
    printf("sizrat %g, stepsize %g, stepsizepre %g\n", sizrat, stepsize_, stepsizepre_);
  }

  // scaled by safety parameter
  sizrat *= sizeratioscale_;

  // optimal new step size
  double stpsiznew = sizrat * stepsize_;

  // redefine sizrat to be dt*_{n}/dt_{n-1}, ie true optimal ratio
  sizrat = stpsiznew / stepsizepre_;

  // limit #sizrat by maximum and minimum
  if (sizrat > sizeratiomax_)
  {
    stpsiznew = sizeratiomax_ * stepsizepre_;
  }
  else if (sizrat < sizeratiomin_)
  {
    stpsiznew = sizeratiomin_ * stepsizepre_;
  }

  // new step size subject to safety measurements
  if (stpsiznew > stepsizemax_)
  {
    stpsiznew = stepsizemax_;
  }
  else if (stpsiznew < stepsizemin_)
  {
    stpsiznew = stepsizemin_;
  }

  return stpsiznew;
}

/*----------------------------------------------------------------------*/
/* Calculate vector norm */
double ADAPTER::StructureTimeAda::CalculateVectorNorm(const enum INPAR::STR::VectorNorm norm,
    const Teuchos::RCP<Epetra_Vector> vect, const int numneglect)
{
  // L1 norm
  if (norm == INPAR::STR::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  else if (norm == INPAR::STR::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  else if (norm == INPAR::STR::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm / sqrt((double)(vect->GlobalLength() - numneglect));
  }
  // infinity/maximum norm
  else if (norm == INPAR::STR::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  else
  {
    dserror("Cannot handle vector norm");
    return 0;
  }
}

/*----------------------------------------------------------------------*/
/* Update output periods */
void ADAPTER::StructureTimeAda::UpdatePeriod()
{
  if (outrest_) outresttime_ += outrestperiod_;
  if (outsys_) outsystime_ += outsysperiod_;
  if (outstr_) outstrtime_ += outstrperiod_;
  if (outene_) outenetime_ += outeneperiod_;

  return;
}

/*----------------------------------------------------------------------*/
INPAR::STR::ConvergenceStatus ADAPTER::StructureTimeAda::PerformErrorAction(
    const INPAR::STR::DivContAct& action, double& stepsizenew)
{
  int myrank = stm_->Discretization()->Comm().MyPID();

  // here we handle how we deal with a failed Newton-Raphson, basically:
  // + stop -> error
  // + halve_step -> reduce size
  // + continue -> estimate posteriori error assuming very big error, throw warning for
  // divergence status
  // + adapt_step -> error
  // + rand_adapt_step -> error
  // + rand_adapt_step_ele_err -> error
  // + repeat_simulation -> error
  // + adapt_penaltycontact ??
  // + adapt_3D0Dptc_ele_err ??
  switch (action)
  {
    case INPAR::STR::divcont_stop:
      // write output
      Output();

      // error and stop the simulation
      dserror("Nonlinear solver did not converge! ");
      break;

    case INPAR::STR::divcont_halve_step:
      if (myrank == 0)
      {
        IO::cout << "Nonlinear solver failed to converge at time t= " << stm_->GetTimeNp()
                 << ". Divide timestep in half. "
                 << "Old time step: " << stepsize_ << IO::endl
                 << "New time step: " << 0.5 * stepsize_ << IO::endl
                 << IO::endl;
      }

      stepsizenew = 0.5 * stepsize_;
      return INPAR::STR::conv_fail_repeat;

    case INPAR::STR::divcont_continue:
      if (myrank == 0)
      {
        IO::cout << "\n WARNING: We are continuing your simulation although the nonlinear solver\n"
                    " did not converge in the current time step. We rely on the error estimator "
                    "to \n"
                    "give a good step size."
                 << IO::endl;
      }

      return INPAR::STR::conv_success;  // Do not surprise. We enforce successful
                                        // status flag to force the error estimator
                                        // to compute new step size later on.

    case INPAR::STR::divcont_adapt_step:
    case INPAR::STR::divcont_rand_adapt_step:
    case INPAR::STR::divcont_rand_adapt_step_ele_err:
      dserror(
          "Adapt the time step is handled by the adaptive time marching integrator. Use\n"
          "DIVERCONT = continue if you want to adapt the step size.");
      break;
    case INPAR::STR::divcont_repeat_simulation:
      dserror("No use to repeat a simulation when it failed. Get a coffee instead.");
      break;
    case INPAR::STR::divcont_adapt_penaltycontact:
    case INPAR::STR::divcont_adapt_3D0Dptc_ele_err:
      dserror(
          "DIVERCONT = adapt_penaltycontact/adapt_3D0Dptc_ele_err is yet to be implemented. "
          "Stay tune.");
      break;
    default:
      dserror("I don't know what to do.");
      break;
  }
  return INPAR::STR::conv_success;  // make compiler happy
}

BACI_NAMESPACE_CLOSE
