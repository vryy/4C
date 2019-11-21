/*----------------------------------------------------------------------*/
/*! \file
\brief submaterial associated with macro-scale Gauss point in multi-scale simulations of scalar
transport problems

\level 2

\maintainer Christoph Schmidt
 */
/*----------------------------------------------------------------------*/
#include "scatra_multiscale_gp.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

#include "../drt_lib/drt_dofset_predefineddofnumber.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_scatra/scatra_timint_ost.H"

#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"

#include "../linalg/linalg_solver.H"

// instantiate static maps
std::map<int, Teuchos::RCP<SCATRA::TimIntOneStepTheta>>
    MAT::ScatraMultiScaleGP::microdisnum_microtimint_map_;
std::map<int, int> MAT::ScatraMultiScaleGP::microdisnum_nummacrogp_map_;

/*--------------------------------------------------------------------*
 | constructor                                             fang 11/15 |
 *--------------------------------------------------------------------*/
MAT::ScatraMultiScaleGP::ScatraMultiScaleGP(const int ele_id,  //!< macro-scale element ID
    const int gp_id,                                           //!< macro-scale Gauss point ID
    const int microdisnum  //!< number of micro-scale discretization
    )
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
      restartname_("")
{
  return;
}  // MAT::ScatraMultiScaleGP::ScatraMultiScaleGP


/*--------------------------------------------------------------------*
 | destructor                                              fang 11/15 |
 *--------------------------------------------------------------------*/
MAT::ScatraMultiScaleGP::~ScatraMultiScaleGP()
{
  // decrement number of macro-scale Gauss points associated with micro-scale time integrator
  --microdisnum_nummacrogp_map_[microdisnum_];

  // once all macro-scale Gauss point submaterials are removed, destruct micro-scale time integrator
  if (microdisnum_nummacrogp_map_[microdisnum_] == 0)
    microdisnum_microtimint_map_[microdisnum_] = Teuchos::null;

  return;
}  // MAT::ScatraMultiScaleGP::~ScatraMultiScaleGP


/*--------------------------------------------------------------------*
 | perform initializations                                 fang 02/16 |
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
      dserror(
          "Must have one-dimensional micro scale in multi-scale simulations of scalar transport "
          "problems!");
    if (DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(sdyn_macro, "TIMEINTEGR") !=
            INPAR::SCATRA::timeint_one_step_theta or
        DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(
            *sdyn_micro, "TIMEINTEGR") != INPAR::SCATRA::timeint_one_step_theta)
      dserror(
          "Multi-scale calculations for scalar transport only implemented for one-step-theta time "
          "integration scheme!");
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
    if (microdis->gNode(0)->X()[0] != 0.)
      dserror(
          "Micro-scale domain must have one end at coordinate 0 and the other end at a coordinate "
          "> 0!");

    // extract multi-scale coupling conditions from micro-scale discretization
    std::vector<Teuchos::RCP<DRT::Condition>> conditions;
    microdis->GetCondition("ScatraMultiScaleCoupling", conditions);

    // safety check
    if (conditions.size() == 0)
      dserror("Couldn't extract multi-scale coupling condition from micro-scale discretization!");

    // loop over all multi-scale coupling conditions
    for (unsigned icond = 0; icond < conditions.size(); ++icond)
    {
      // extract nodal cloud
      const std::vector<int>* const nodeids = conditions[icond]->Nodes();
      if (nodeids == NULL) dserror("Multi-scale coupling condition does not have nodal cloud!");

      // loop over all nodes in nodal cloud
      for (unsigned inode = 0; inode < (*nodeids).size(); ++inode)
      {
        if (microdis->NodeRowMap()->MyGID((*nodeids)[inode]))
        {
          // extract node from micro-scale discretization
          DRT::Node* node = microdis->gNode((*nodeids)[inode]);

          // safety checks
          if (node == NULL)
            dserror("Cannot extract node with global ID %d from micro-scale discretization!",
                (*nodeids)[inode]);
          if (node->X()[0] <= 0.)
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
      dserror(
          "No linear solver defined for scalar field in input file for micro scale! Please set "
          "LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

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
    microdisnum_microtimint_map_[microdisnum_]->Setup();

    // set initial velocity field
    microdisnum_microtimint_map_[microdisnum_]->SetVelocityField(1);

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

  return;
}  // MAT::ScatraMultiScaleGP::Init


/*--------------------------------------------------------------------*
 | prepare time step                                       fang 11/15 |
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScaleGP::PrepareTimeStep(
    const std::vector<double>& phinp_macro  //!< macro-scale state variables
)
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

  return;
}  // MAT::ScatraMultiScaleGP::PrepareTimeStep


/*--------------------------------------------------------------------*
 | evaluate micro scale                                    fang 11/15 |
 *--------------------------------------------------------------------*/
void MAT::ScatraMultiScaleGP::Evaluate(
    const std::vector<double>& phinp_macro,  //!< macro-scale state variables
    double& q_micro,                         //!< micro-scale coupling flux
    std::vector<double>& dq_dphi_micro,      //!< derivatives of micro-scale coupling flux w.r.t.
                                             //!< macro-scale state variables
    const bool solve  //!< flag indicating whether micro-scale problem should be solved
)
{
  // extract micro-scale time integrator
  const Teuchos::RCP<SCATRA::TimIntOneStepTheta>& microtimint =
      microdisnum_microtimint_map_[microdisnum_];

  // set current state in micro-scale time integrator
  microtimint->SetState(phin_, phinp_, phidtn_, phidtnp_, hist_, micro_output_, phinp_macro, step_,
      DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance("scatra")->Time());

  if (step_ == 0 or !solve)
    // only evaluate the micro-scale coupling quantities without solving the entire micro-scale
    // problem relevant for truly partitioned multi-scale simulations or for calculation of initial
    // time derivative of macro-scale state vector
    microtimint->EvaluateMacroMicroCoupling();
  else
    // solve micro-scale problem
    // note that it is not necessary to transfer the final micro-scale state vectors back to the
    // Gauss-point submaterial due to RCP usage
    microtimint->Solve();

  // transfer micro-scale coupling quantities to macro scale
  q_micro = -microtimint->Q();
  dq_dphi_micro = microtimint->DqDphi();
  for (unsigned i = 0; i < dq_dphi_micro.size(); ++i) dq_dphi_micro[i] *= -1.;

  // clear state in micro-scale time integrator
  microtimint->ClearState();

  return;
}  // MAT::ScatraMultiScaleGP::Evaluate


/*--------------------------------------------------------------------*
 | evaluate mean concentration on micro scale              fang 08/17 |
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
  eleparams.set<int>("action", SCATRA::calc_total_and_mean_scalars);
  eleparams.set("inverting", false);

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
}  // MAT::ScatraMultiScaleGP::EvaluateMeanConcentration


/*-------------------------------------------------------------------------*
 | evaluate mean concentration time derivative on micro scale   fang 03/18 |
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
  eleparams.set<int>("action", SCATRA::calc_mean_scalar_time_derivatives);

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
}  // MAT::ScatraMultiScaleGP::EvaluateMeanConcentrationTimeDerivative


/*------------------------------------------------------------------------------*
 | update micro-scale time integrator at the end of each time step   fang 12/15 |
 *------------------------------------------------------------------------------*/
void MAT::ScatraMultiScaleGP::Update()
{
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

  return;
}  // MAT::ScatraMultiScaleGP::Update


/*--------------------------------------------------------------------*
 | create new result file                                  fang 02/16 |
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
  size_t pos = microprefix.rfind('-');
  if (pos != std::string::npos)
  {
    std::string number = microprefix.substr(pos + 1);
    std::string prefix = microprefix.substr(0, pos);
    std::ostringstream s;
    s << prefix << "_microdis" << microdisnum_ << "_el" << ele_id_ << "_gp" << gp_id_ << "-"
      << number;
    restartname_ = s.str();
  }
  else
  {
    std::ostringstream s;
    s << microprefix << "_microdis" << microdisnum_ << "_el" << ele_id_ << "_gp" << gp_id_;
    restartname_ = s.str();
  }

  // figure out new prefix for micro-scale output files
  // note that the trailing number must be the same as for the macro scale
  std::string newfilename;
  size_t posn = micronewprefix.rfind('-');
  if (posn != std::string::npos)
  {
    std::string number = micronewprefix.substr(posn + 1);
    std::string prefix = micronewprefix.substr(0, posn);
    std::ostringstream s;
    s << prefix << "_microdis" << microdisnum_ << "_el" << ele_id_ << "_gp" << gp_id_ << "-"
      << number;
    newfilename = s.str();
  }
  else
  {
    std::ostringstream s;
    s << micronewprefix << "_microdis" << microdisnum_ << "_el" << ele_id_ << "_gp" << gp_id_;
    newfilename = s.str();
  }

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

  return;
}  // MAT::ScatraMultiScaleGP::NewResultFile


/*--------------------------------------------------------------------*
 | output micro-scale quantities                           fang 02/16 |
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

    // clear state in micro-scale time integrator
    microtimint->ClearState();
  }

  return;
}  // MAT::ScatraMultiScaleGP::Output


/*--------------------------------------------------------------------*
 | read restart on micro scale                             fang 03/16 |
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
  Teuchos::RCP<IO::InputControl> inputcontrol =
      Teuchos::rcp(new IO::InputControl(restartname_, true));
  microtimint->ReadRestart(step_, inputcontrol);

  // safety check
  if (microtimint->Step() != step_) dserror("Time step mismatch!");

  // clear state in micro-scale time integrator
  microtimint->ClearState();

  return;
}  // MAT::ScatraMultiScaleGP::ReadRestart
