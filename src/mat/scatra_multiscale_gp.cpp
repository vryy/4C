/*----------------------------------------------------------------------*/
/*! \file
\brief submaterial associated with macro-scale Gauss point in multi-scale simulations of scalar
transport problems

\level 2

 */
/*----------------------------------------------------------------------*/
#include "scatra_multiscale_gp.H"

#include "io.H"
#include "io_control.H"

#include "dofset_predefineddofnumber.H"
#include "globalproblem.H"
#include "utils_parameter_list.H"

#include "scatra_timint_ost.H"

#include "scatra_ele_action.H"
#include "scatra_ele_parameter_timint.H"

#include "linalg_solver.H"

#include <filesystem>

// instantiate static maps
std::map<int, Teuchos::RCP<SCATRA::TimIntOneStepTheta>>
    MAT::ScatraMultiScaleGP::microdisnum_microtimint_map_;
std::map<int, int> MAT::ScatraMultiScaleGP::microdisnum_nummacrogp_map_;

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::ScatraMultiScaleGP::ScatraMultiScaleGP(
    const int ele_id, const int gp_id, const int microdisnum, const bool is_ale)
    : gp_id_(gp_id),
      ele_id_(ele_id),
      eleowner_(DRT::Problem::Instance()->GetDis("scatra")->ElementRowMap()->MyGID(ele_id)),
      microdisnum_(microdisnum),
      step_(0),
      phin_(Teuchos::null),
      phinp_(Teuchos::null),
      phidtn_(Teuchos::null),
      phidtnp_(Teuchos::null),
      hist_(Teuchos::null),
      micro_output_(Teuchos::null),
      restartname_(""),
      detFn_(1.0),
      detFnp_(1.0),
      ddetFdtn_(0.0),
      ddetFdtnp_(0.0),
      is_ale_(is_ale)

{
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
MAT::ScatraMultiScaleGP::~ScatraMultiScaleGP()
{
  // decrement number of macro-scale Gauss points associated with micro-scale time integrator
  --microdisnum_nummacrogp_map_[microdisnum_];

  // once all macro-scale Gauss point submaterials are removed, destruct micro-scale time integrator
  if (microdisnum_nummacrogp_map_[microdisnum_] == 0)
    microdisnum_microtimint_map_[microdisnum_] = Teuchos::null;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScaleGP::Init()
{
  // extract micro-scale problem
  DRT::Problem* microproblem = DRT::Problem::Instance(microdisnum_);

  // extract micro-scale discretization
  std::stringstream microdisname;
  microdisname << "scatra_multiscale_" << microdisnum_;
  Teuchos::RCP<DRT::Discretization> microdis = microproblem->GetDis(microdisname.str());

  // instantiate and initialize micro-scale state vectors
  phin_ = LINALG::CreateVector(*microdis->DofRowMap(), true);
  phinp_ = LINALG::CreateVector(*microdis->DofRowMap(), true);
  phidtn_ = LINALG::CreateVector(*microdis->DofRowMap(), true);
  phidtnp_ = LINALG::CreateVector(*microdis->DofRowMap(), true);
  hist_ = LINALG::CreateVector(*microdis->DofRowMap(), true);

  // set up micro-scale time integrator for micro-scale problem if not already done
  if (microdisnum_microtimint_map_.find(microdisnum_) == microdisnum_microtimint_map_.end() or
      microdisnum_microtimint_map_[microdisnum_] == Teuchos::null)
  {
    // extract macro-scale parameter list
    const Teuchos::ParameterList& sdyn_macro =
        DRT::Problem::Instance()->ScalarTransportDynamicParams();

    // extract micro-scale parameter list and create deep copy
    Teuchos::RCP<Teuchos::ParameterList> sdyn_micro = Teuchos::rcp(new Teuchos::ParameterList(
        DRT::Problem::Instance(microdisnum_)->ScalarTransportDynamicParams()));

    // preliminary safety check
    if (DRT::Problem::Instance(microdisnum_)->NDim() != 1)
    {
      dserror(
          "Must have one-dimensional micro scale in multi-scale simulations of scalar transport "
          "problems!");
    }
    if (DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(sdyn_macro, "TIMEINTEGR") !=
            INPAR::SCATRA::timeint_one_step_theta or
        DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(
            *sdyn_micro, "TIMEINTEGR") != INPAR::SCATRA::timeint_one_step_theta)
    {
      dserror(
          "Multi-scale calculations for scalar transport only implemented for one-step-theta time "
          "integration scheme!");
    }
    if (DRT::INPUT::IntegralValue<bool>(sdyn_macro, "SKIPINITDER") !=
        DRT::INPUT::IntegralValue<bool>(*sdyn_micro, "SKIPINITDER"))
      dserror("Flag SKIPINITDER in input file must be equal on macro and micro scales!");
    if (sdyn_macro.get<double>("TIMESTEP") != sdyn_micro->get<double>("TIMESTEP"))
      dserror("Must have identical time step size on macro and micro scales!");
    if (sdyn_macro.get<int>("NUMSTEP") != sdyn_micro->get<int>("NUMSTEP"))
      dserror("Must have identical number of time steps on macro and micro scales!");
    if (sdyn_macro.get<double>("THETA") != sdyn_micro->get<double>("THETA"))
      dserror(
          "Must have identical one-step-theta time integration factor on macro and micro scales!");
    if (microdis->NumGlobalElements() == 0)
      dserror("No elements in TRANSPORT ELEMENTS section of micro-scale input file!");
    if (microdis->gNode(0)->X()[0] != 0.0)
    {
      dserror(
          "Micro-scale domain must have one end at coordinate 0 and the other end at a coordinate "
          "> 0!");
    }

    // extract multi-scale coupling conditions from micro-scale discretization
    std::vector<Teuchos::RCP<DRT::Condition>> conditions;
    microdis->GetCondition("ScatraMultiScaleCoupling", conditions);

    // safety check
    if (conditions.size() == 0)
      dserror("Couldn't extract multi-scale coupling condition from micro-scale discretization!");

    // loop over all multi-scale coupling conditions
    for (auto& condition : conditions)
    {
      // extract nodal cloud
      const std::vector<int>* const nodeids = condition->Nodes();
      if (nodeids == nullptr) dserror("Multi-scale coupling condition does not have nodal cloud!");

      // loop over all nodes in nodal cloud
      for (int inode : *nodeids)
      {
        if (microdis->NodeRowMap()->MyGID(inode))
        {
          // extract node from micro-scale discretization
          DRT::Node* node = microdis->gNode(inode);

          // safety checks
          if (node == nullptr)
          {
            dserror(
                "Cannot extract node with global ID %d from micro-scale discretization!", inode);
          }
          else if (node->X()[0] <= 0.0)
            dserror(
                "Multi-scale coupling condition must be enforced on a node with coordinate > 0!");
        }
      }
    }

    // add proxy of velocity related degrees of freedom to scatra discretization
    Teuchos::RCP<DRT::DofSetInterface> dofsetaux = Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(
        DRT::Problem::Instance(microdisnum_)->NDim() + 1, 0, 0, true));
    if (microdis->AddDofSet(dofsetaux) != 1)
      dserror("Micro-scale discretization has illegal number of dofsets!");

    // finalize discretization
    microdis->FillComplete(true, false, false);

    // get solver number
    const int linsolvernumber = sdyn_micro->get<int>("LINEAR_SOLVER");

    // check solver number
    if (linsolvernumber < 0)
    {
      dserror(
          "No linear solver defined for scalar field in input file for micro scale! Please set "
          "LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");
    }

    // create solver
    Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(
        new LINALG::Solver(DRT::Problem::Instance(microdisnum_)->SolverParams(linsolvernumber),
            microdis->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

    // provide solver with null space information if necessary
    microdis->ComputeNullSpaceIfNecessary(solver->Params());

    // supplementary parameter list
    Teuchos::RCP<Teuchos::ParameterList> extraparams = Teuchos::rcp(new Teuchos::ParameterList());
    extraparams->set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());
    extraparams->set<bool>("isale", false);
    extraparams->sublist("TURBULENT INFLOW") =
        DRT::Problem::Instance(microdisnum_)->FluidDynamicParams().sublist("TURBULENT INFLOW");
    extraparams->sublist("TURBULENCE MODEL") =
        DRT::Problem::Instance(microdisnum_)->FluidDynamicParams().sublist("TURBULENCE MODEL");

    // instantiate and initialize micro-scale time integrator
    microdisnum_microtimint_map_[microdisnum_] = Teuchos::rcp(new SCATRA::TimIntOneStepTheta(
        microdis, solver, sdyn_micro, extraparams, Teuchos::null, microdisnum_));
    microdisnum_microtimint_map_[microdisnum_]->Init();
    microdisnum_microtimint_map_[microdisnum_]->SetNumberOfDofSetVelocity(1);
    microdisnum_microtimint_map_[microdisnum_]->Setup();

    // set initial velocity field
    microdisnum_microtimint_map_[microdisnum_]->SetVelocityField();

    // create counter for number of macro-scale Gauss points associated with micro-scale time
    // integrator
    microdisnum_nummacrogp_map_[microdisnum_] = 0;
  }

  // increment counter
  ++microdisnum_nummacrogp_map_[microdisnum_];

  // extract initial state vectors from micro-scale time integrator
  phin_->Scale(1., *microdisnum_microtimint_map_[microdisnum_]->Phin());
  phinp_->Scale(1., *microdisnum_microtimint_map_[microdisnum_]->Phinp());

  // create new result file
  NewResultFile();
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScaleGP::PrepareTimeStep(const std::vector<double>& phinp_macro)
{
  // extract micro-scale time integrator
  const Teuchos::RCP<SCATRA::TimIntOneStepTheta>& microtimint =
      microdisnum_microtimint_map_[microdisnum_];

  // set current state in micro-scale time integrator
  microtimint->SetState(phin_, phinp_, phidtn_, phidtnp_, hist_, micro_output_, phinp_macro, step_,
      DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->Time());

  // prepare time step
  microtimint->PrepareTimeStep();

  // clear state in micro-scale time integrator
  microtimint->ClearState();

  // increment time step
  ++step_;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScaleGP::Evaluate(const std::vector<double>& phinp_macro, double& q_micro,
    std::vector<double>& dq_dphi_micro, const double detFnp, const bool solve)
{
  // extract micro-scale time integrator
  const Teuchos::RCP<SCATRA::TimIntOneStepTheta>& microtimint =
      microdisnum_microtimint_map_[microdisnum_];

  // set current state in micro-scale time integrator
  microtimint->SetState(phin_, phinp_, phidtn_, phidtnp_, hist_, micro_output_, phinp_macro, step_,
      DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->Time());

  if (is_ale_)
  {
    // update determinant of deformation gradient
    detFnp_ = detFnp;

    // calculate time derivative and pass to micro time integration as reaction coefficient
    CalculateDdetFDt(microtimint);
    microtimint->SetMacroMicroReaCoeff(ddetFdtnp_);
  }

  if (step_ == 0 or !solve)
  {
    // only evaluate the micro-scale coupling quantities without solving the entire micro-scale
    // problem relevant for truly partitioned multi-scale simulations or for calculation of initial
    // time derivative of macro-scale state vector
    microtimint->EvaluateMacroMicroCoupling();
  }
  else
  {
    // solve micro-scale problem
    // note that it is not necessary to transfer the final micro-scale state vectors back to the
    // Gauss-point submaterial due to RCP usage
    microtimint->Solve();
  }

  // transfer micro-scale coupling quantities to macro scale
  q_micro = -microtimint->Q();
  dq_dphi_micro = microtimint->DqDphi();
  for (double& dq_dphi_micro_component : dq_dphi_micro) dq_dphi_micro_component *= -1.0;

  // clear state in micro-scale time integrator
  microtimint->ClearState();
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
double MAT::ScatraMultiScaleGP::EvaluateMeanConcentration() const
{
  // extract micro-scale discretization
  DRT::Discretization& discret = *microdisnum_microtimint_map_[microdisnum_]->Discretization();

  // set micro-scale state vector
  discret.ClearState();
  discret.SetState("phinp", phinp_);

  // set parameters for micro-scale elements
  Teuchos::ParameterList eleparams;
  DRT::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::calc_total_and_mean_scalars, eleparams);
  eleparams.set("inverting", false);
  eleparams.set("calc_grad_phi", false);

  // initialize result vector: first component = concentration integral, second component = domain
  // integral
  const Teuchos::RCP<Epetra_SerialDenseVector> integrals =
      Teuchos::rcp(new Epetra_SerialDenseVector(2));

  // evaluate concentration and domain integrals on micro scale
  discret.EvaluateScalars(eleparams, integrals);

  // clear discretization
  discret.ClearState();

  // compute and return mean concentration on micro scale
  return (*integrals)[0] / (*integrals)[1];
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
double MAT::ScatraMultiScaleGP::EvaluateMeanConcentrationTimeDerivative() const
{
  // extract micro-scale discretization
  DRT::Discretization& discret = *microdisnum_microtimint_map_[microdisnum_]->Discretization();

  // set micro-scale state vector
  discret.ClearState();
  discret.SetState("phidtnp", phidtnp_);

  // set parameters for micro-scale elements
  Teuchos::ParameterList eleparams;
  DRT::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::calc_mean_scalar_time_derivatives, eleparams);

  // initialize result vector: first component = integral of concentration time derivative, second
  // component = integral of domain
  const Teuchos::RCP<Epetra_SerialDenseVector> integrals =
      Teuchos::rcp(new Epetra_SerialDenseVector(2));

  // evaluate integrals of domain and time derivative of concentration on micro scale
  discret.EvaluateScalars(eleparams, integrals);

  // clear discretization
  discret.ClearState();

  // compute and return mean concentration time derivative on micro scale
  return (*integrals)[0] / (*integrals)[1];
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
void MAT::ScatraMultiScaleGP::Update()
{
  if (is_ale_)
  {
    // Update detF
    detFn_ = detFnp_;
    ddetFdtn_ = ddetFdtnp_;
  }

  // extract micro-scale time integrator
  Teuchos::RCP<SCATRA::TimIntOneStepTheta> microtimint = microdisnum_microtimint_map_[microdisnum_];

  // set current state in micro-scale time integrator
  microtimint->SetState(phin_, phinp_, phidtn_, phidtnp_, hist_, micro_output_,
      std::vector<double>(0, 0.), step_,
      DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->Time());

  // update micro-scale time integrator
  microtimint->Update();

  // clear state in micro-scale time integrator
  microtimint->ClearState();
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScaleGP::NewResultFile()
{
  // get properties from macro scale
  Teuchos::RCP<IO::OutputControl> macrocontrol = DRT::Problem::Instance()->OutputControlFile();
  std::string microprefix = macrocontrol->RestartName();
  std::string micronewprefix = macrocontrol->NewOutputFileName();

  // extract micro-scale problem and discretization
  DRT::Problem* microproblem = DRT::Problem::Instance(microdisnum_);
  std::stringstream microdisname;
  microdisname << "scatra_multiscale_" << microdisnum_;
  Teuchos::RCP<DRT::Discretization> microdis = microproblem->GetDis(microdisname.str());

  // figure out prefix of micro-scale restart files
  restartname_ = NewResultFilePath(microprefix);

  // figure out new prefix for micro-scale output files
  const std::string newfilename = NewResultFilePath(micronewprefix);

  if (eleowner_)
  {
    const int ndim = microproblem->NDim();
    const int restart = DRT::Problem::Instance()->Restart();
    bool adaptname = true;

    // in case of restart, the new output file name has already been adapted
    if (restart) adaptname = false;

    Teuchos::RCP<IO::OutputControl> microcontrol = Teuchos::rcp(new IO::OutputControl(
        microdis->Comm(), "Scalar_Transport", microproblem->SpatialApproximationType(),
        "micro-input-file-not-known", restartname_, newfilename, ndim, restart,
        DRT::Problem::Instance(microdisnum_)->IOParams().get<int>("FILESTEPS"),
        DRT::INPUT::IntegralValue<bool>(
            DRT::Problem::Instance(microdisnum_)->IOParams(), "OUTPUT_BIN"),
        adaptname));

    micro_output_ = Teuchos::rcp(new IO::DiscretizationWriter(microdis));
    micro_output_->SetOutput(microcontrol);
    micro_output_->WriteMesh(
        step_, DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->Time());
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
const std::string MAT::ScatraMultiScaleGP::NewResultFilePath(const std::string& newprefix)
{
  std::string newfilename;

  // create path from string to extract only filename prefix
  const std::filesystem::path path(newprefix);
  const std::string newfileprefix = path.filename().string();

  const size_t posn = newfileprefix.rfind('-');
  if (posn != std::string::npos)
  {
    const std::string number = newfileprefix.substr(posn + 1);
    const std::string prefix = newfileprefix.substr(0, posn);

    // recombine path and file
    const std::filesystem::path parent_path(path.parent_path());
    const std::filesystem::path filen_name(prefix);
    const std::filesystem::path recombined_path = parent_path / filen_name;

    std::ostringstream s;
    s << recombined_path.string() << "_microdis" << microdisnum_ << "_el" << ele_id_ << "_gp"
      << gp_id_ << "-" << number;
    newfilename = s.str();
  }
  else
  {
    std::ostringstream s;
    s << newprefix << "_microdis" << microdisnum_ << "_el" << ele_id_ << "_gp" << gp_id_;
    newfilename = s.str();
  }

  return newfilename;
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScaleGP::Output()
{
  // skip ghosted macro-scale elements
  if (eleowner_)
  {
    // extract micro-scale time integrator
    Teuchos::RCP<SCATRA::TimIntOneStepTheta> microtimint =
        microdisnum_microtimint_map_[microdisnum_];

    // set current state in micro-scale time integrator
    microtimint->SetState(phin_, phinp_, phidtn_, phidtnp_, hist_, micro_output_,
        std::vector<double>(0, 0.), step_,
        DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->Time());

    // output micro-scale quantities
    microtimint->Output();

    if (microtimint->DoOutputRestart() and is_ale_)
    {
      microtimint->DiscWriter()->WriteDouble("detFn", detFn_);
      microtimint->DiscWriter()->WriteDouble("ddetFdtn", ddetFdtn_);
    }

    // clear state in micro-scale time integrator
    microtimint->ClearState();
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScaleGP::ReadRestart()
{
  // extract micro-scale time integrator
  Teuchos::RCP<SCATRA::TimIntOneStepTheta> microtimint = microdisnum_microtimint_map_[microdisnum_];

  // extract restart step
  step_ = DRT::Problem::Instance()->Restart();

  // set current state in micro-scale time integrator
  microtimint->SetState(phin_, phinp_, phidtn_, phidtnp_, hist_, micro_output_,
      std::vector<double>(0, 0.), step_,
      DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->Time());

  // read restart on micro scale
  auto inputcontrol = Teuchos::rcp(new IO::InputControl(restartname_, true));
  microtimint->ReadRestart(step_, inputcontrol);

  // safety check
  if (microtimint->Step() != step_) dserror("Time step mismatch!");

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (inputcontrol == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(microtimint->Discretization(), step_));
  else
    reader = Teuchos::rcp(
        new IO::DiscretizationReader(microtimint->Discretization(), inputcontrol, step_));

  if (is_ale_)
  {
    detFn_ = reader->ReadDouble("detFn");
    ddetFdtn_ = reader->ReadDouble("ddetFdtn");
  }

  // clear state in micro-scale time integrator
  microtimint->ClearState();
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScaleGP::CalculateDdetFDt(Teuchos::RCP<SCATRA::TimIntOneStepTheta> microtimint)
{
  const double dt = microtimint->Dt();

  switch (microtimint->MethodName())
  {
    case INPAR::SCATRA::TimeIntegrationScheme::timeint_one_step_theta:
    {
      const double theta = microtimint->ScatraParameterList()->get<double>("THETA");

      const double part1 = (detFnp_ - detFn_) / dt;
      const double part2 = (1.0 - theta) * ddetFdtn_;
      ddetFdtnp_ = 1.0 / theta * (part1 - part2);

      break;
    }
    default:
    {
      dserror("time integration scheme not supported to calculate d detF / d t.");
      break;
    }
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScaleGP::SetTimeStepping(const double dt, const double time, const int step)
{
#ifdef DEBUG
  dsassert(dt > 0.0, "Time step for micro scale must be positive.");
  dsassert(time >= 0.0, "Time for micro scale must be positive.");
  dsassert(step >= 0, "Number of step for micro scale must be positive.");
#endif

  Teuchos::RCP<SCATRA::TimIntOneStepTheta> microtimint = microdisnum_microtimint_map_[microdisnum_];
  microtimint->SetDt(dt);
  microtimint->SetTimeStep(time, step);
}