/*----------------------------------------------------------------------------*/
/*!
\file fsi_monolithic_timintada.cpp

\brief Provide time adaptivity functionalities within the Monolithic class

\level 3

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "fsi_monolithic.H"

#include "../drt_inpar/inpar_fsi.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_timestepping/timintmstep.H"

#include "../linalg/linalg_utils.H"

#include "../drt_adapter/ad_ale_fsi.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"
#include "../drt_adapter/ad_str_fsi_timint_adaptive.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::InitTimIntAda(const Teuchos::ParameterList& fsidyn)
{
  // access to structural time adaptivity parameters
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& sada = sdyn.sublist("TIMEADAPTIVITY");

  // access to FSI time adaptivity parameters
  const Teuchos::ParameterList& fsiada = fsidyn.sublist("TIMEADAPTIVITY");

  // initialize member variables
  adaptstep_ = 0;
  accepted_ = false;

  adareason_ = "none";

  // L2-norms of estimation of temporal discretization errors
  strnorm_ = 0.0;
  flnorm_ = 0.0;
  strfsinorm_ = 0.0;
  flfsinorm_ = 0.0;
  strinnernorm_ = 0.0;
  flinnernorm_ = 0.0;

  // L-inf-norms of estimation of temporal discretization errors
  strinfnorm_ = 0.0;
  flinfnorm_ = 0.0;
  strinffsinorm_ = 0.0;
  flinffsinorm_ = 0.0;
  strinfinnernorm_ = 0.0;
  flinfinnernorm_ = 0.0;

  // time step sizes calculated according to the 6 available L2-norms
  dtstr_ = 0.0;
  dtfl_ = 0.0;
  dtstrfsi_ = 0.0;
  dtflfsi_ = 0.0;
  dtstrinner_ = 0.0;
  dtflinner_ = 0.0;
  dtnonlinsolver_ = 0.0;

  dtminused_ = false;

  numincreasesteps_ = fsiada.get<int>("NUMINCREASESTEPS");

  // read weights for averaging strategy in case of increasing time step size
  {
    double word;
    std::istringstream avgweightstream(Teuchos::getNumericStringParameter(fsiada, "AVERAGINGDT"));
    while (avgweightstream >> word) avgweights_.push_back(word);

    // check if weights for averaging make sense
    double sum = 0.0;
    for (std::vector<double>::iterator it = avgweights_.begin(); it < avgweights_.end(); ++it)
      sum += (*it);
    if (sum < 0.99)
    {
      dserror("Sum of weights for dt: %f < 1.0", sum);
    }
    if (sum > 1.01)
    {
      dserror("Sum of weights for dt: %f > 1.0", sum);
    }
  }

  dt_ = Teuchos::rcp(new TIMINT::TimIntMStep<double>(-avgweights_.size(), 1, 0.0));
  dt_->SetStep(1, Dt());

  //----------------------------------------------------------------------------
  // write adaptivity file
  //----------------------------------------------------------------------------
  std::string fileada = DRT::Problem::Instance()->OutputControlFile()->FileName();
  fileada.append(".adaptivity");
  logada_ = Teuchos::rcp(new std::ofstream(fileada.c_str()));

  //----------------------------------------------------------------------------
  // set algorithmic parameters
  //----------------------------------------------------------------------------
  dtmax_ = fsiada.get<double>("DTMAX");
  dtmin_ = fsiada.get<double>("DTMIN");

  errtolstr_ = sada.get<double>("LOCERRTOL");
  errtolfl_ = fsiada.get<double>("LOCERRTOLFLUID");

  //----------------------------------------------------------------------------
  // check on which fields time adaptivity should be based on
  //----------------------------------------------------------------------------
  if (not(DRT::INPUT::IntegralValue<INPAR::STR::TimAdaKind>(sada, "KIND") ==
          INPAR::STR::timada_kind_none))
    isadastructure_ = true;

  if (not(DRT::INPUT::IntegralValue<int>(fsiada, "AUXINTEGRATORFLUID") ==
          INPAR::FSI::timada_fld_none))
    isadafluid_ = true;

  // get error action strategy from input file
  const int erroractionstrategy = DRT::INPUT::IntegralValue<int>(fsiada, "DIVERCONT");
  if (erroractionstrategy != INPAR::FSI::divcont_stop and
      erroractionstrategy != INPAR::FSI::divcont_continue)
    isadasolver_ = true;

  flmethod_ = FluidField()->GetTimAdaMethodName();

  //----------------------------------------------------------------------------
  // Handling of Dirichlet BCs in error estimation
  //----------------------------------------------------------------------------
  // Create intersection of fluid DOFs that hold a Dirichlet boundary condition
  // and are located at the FSI interface.
  std::vector<Teuchos::RCP<const Epetra_Map>> intersectionmapsfluid;
  intersectionmapsfluid.push_back(FluidField()->GetDBCMapExtractor()->CondMap());
  intersectionmapsfluid.push_back(FluidField()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmapfluid =
      LINALG::MultiMapExtractor::IntersectMaps(intersectionmapsfluid);

  // store number of interface DOFs subject to Dirichlet BCs on structure and fluid side of the
  // interface
  numflfsidbcdofs_ = intersectionmapfluid->NumGlobalElements();

  //----------------------------------------------------------------------------
  // Check whether input parameters make sense
  //----------------------------------------------------------------------------
  const double safetyfac = fsiada.get<double>("SAFETYFACTOR");
  if (safetyfac > 1.0)
    dserror(
        "SAFETYFACTOR in FSI DYNAMIC/TIMEADAPTIVITY is %f > 1.0 and, "
        "thus, too large.",
        safetyfac);

  if (dtmax_ < dtmin_)
    dserror("DTMAX = %f is not allowed to be smaller than DTMIN = %f.", dtmax_, dtmin_);

  /* safety check for BDF2 time integration in fluid
   * (cf. PhD Thesis [Bornemann, 2003, p. 61, eq. (3.40)]) */
  if (FluidField()->TimIntScheme() == INPAR::FLUID::timeint_bdf2 and
      fsiada.get<double>("SIZERATIOMAX") >= 1.0 + sqrt(2))
  {
    dserror(
        "In case of BDF2 time integration, the maximum increase of "
        "two subsequent time step sizes is limited to be less than '1+sqrt(2)'"
        "due to stability reasons "
        "(cf. PhD Thesis [Bornemann, 2003, p. 61, eq. (3.40)]). "
        "Choose an appropriate value!");
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::TimeloopAdaDt(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  /*--------------------------------------------------------------------------*/
  /* Initialize some parameters                                               */
  /*--------------------------------------------------------------------------*/
  // get the FSI parameter list
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  // number of maximum allowed time step repetitions
  const int adaptstepmax = fsidyn.sublist("TIMEADAPTIVITY").get<int>("ADAPTSTEPMAX");
  /*--------------------------------------------------------------------------*/

  // resize MStep objects of structure needed for AB2 structural extrapolation of displacements
  StructureField()->ResizeMStepTimAda();

#ifdef DEBUG
  // check whether all fields have the same time step size
  CheckIfDtsSame();
#endif

  PrepareTimeloop();

  // the time loop
  while (NotFinished())
  {
    PrepareAdaptiveTimeStep();

    // time step adaptivity loop
    while (StepNotAccepted() and (adaptstep_ < adaptstepmax) and (not dtminused_))
    {
      PrintHeaderRepeatedStep();
      PrepareTimeStep();

      // Do the auxiliary step first
      TimeStepAuxiliar();

      // Do the time step with the marching time integrator
      TimeStep(interface);

      // adjust the time step size
      AdaptTimeStepSize();

      /* continue with the next time step if the step has been repeated with the minimum
       * time step size --> no more refinement is possible anyway!*/
      if (dtminused_)
      {
        if (Comm().MyPID() == 0)
        {
          IO::cout << "Time Step has already been calculated with minimum step size"
                   << " --> continue with next time step!"
                   << "\n";
        }

        accepted_ = true;
      }

      if (StepNotAccepted())
      {
        if (adaptstep_ >= adaptstepmax)
        {
          if (not IsAdaSolver())
          {
            if (Comm().MyPID() == 0)
            {
              IO::cout << "adaptstep_ = " << adaptstepmax
                       << " --> Max. number of adaption iterations is reached:"
                       << " continuing with next time step! "
                       << "\n";
            }
          }
          else if (IsAdaSolver() and not(erroraction_ == FSI::Monolithic::erroraction_none or
                                         erroraction_ == FSI::Monolithic::erroraction_continue))
          {
            dserror(
                "Nonlinear solver did not converge after %i repetitions of "
                "time step %i.",
                adaptstepmax, Step());
          }
        }
        else
        {
          ResetStep();
          ResetTime();

          if (Comm().MyPID() == 0)
          {
            IO::cout << "Repeat current time step with dt = " << Dt() << " based on " << adareason_
                     << ".\n";
          }
        }
      }

#ifdef DEBUG
      // check whether all fields have the same time step size
      CheckIfDtsSame();
#endif

      PrintAdaptivitySummary();
    }

    PrepareOutput();
    Update();            // Update in single fields
    UpdateDtPast(Dt());  // Update of adaptive time step sizes
    Output();

    WriteAdaFile();
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::PrepareAdaptiveTimeStep()
{
  // reset from last time step
  dtminused_ = false;
  adaptstep_ = 0;
  accepted_ = false;

  if (Comm().MyPID() == 0)
  {
    IO::cout << "\n"
             << "+++++++++++++++++++++++++++++NEW TIME STEP+++++++++++++++++++++++++++++";
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::PrintHeaderRepeatedStep() const
{
  if (adaptstep_ != 0 and Comm().MyPID() == 0)
  {
    IO::cout << IO::endl
             << "__________REAPEATING TIME STEP " << Step() << " WITH DT = " << Dt() << " FOR THE "
             << adaptstep_ << ". TIME__________" << IO::endl;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::WriteAdaFileHeader() const
{
  // write to adaptivity file
  if (Comm().MyPID() == 0 and (not logada_.is_null()))
  {
    // get string of type of auxiliary time integration scheme in structure field
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    const Teuchos::ParameterList& sada = sdyn.sublist("TIMEADAPTIVITY");
    std::string strmethod = sada.get<std::string>("KIND");

    // print the actual header
    (*logada_) << "Time Adaptivity in monolithic Fluid-Structure-Interaction:"
               << "\n"
               << " - Error estimation method in fluid field: " << flmethod_ << "\n"
               << " - Error estimation method in structure field: " << strmethod << "\n"
               << " - Error tolerance in fluid field: " << errtolfl_ << "\n"
               << " - Error tolerance in structure field: " << errtolstr_ << "\n"
               << " - Minimum allowed time step size DTMIN = " << dtmin_ << "\n"
               << " - Maximum allowed time step size DTMAX = " << dtmax_ << "\n \n"
               << std::right << std::setw(9) << "step" << std::right << std::setw(16) << "time"
               << std::right << std::setw(16) << "dt" << std::right << std::setw(16) << "#adaptiter"
               << std::right << std::setw(16) << "dt of field" << std::right << std::setw(16)
               << "L2-err str field" << std::right << std::setw(16) << "L2-err str inner"
               << std::right << std::setw(16) << "L2-err str fsi" << std::right << std::setw(16)
               << "L2-err fl field" << std::right << std::setw(16) << "L2-err fl inner"
               << std::right << std::setw(16) << "L2-err fl fsi" << std::right << std::setw(16)
               << "Linf str field" << std::right << std::setw(16) << "Linf str inner" << std::right
               << std::setw(16) << "Linf str fsi" << std::right << std::setw(16) << "Linf fl field"
               << std::right << std::setw(16) << "Linf fl inner" << std::right << std::setw(16)
               << "Linf fl fsi"
               << "\n\n";
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::WriteAdaFile() const
{
  if (logada_.is_null()) dserror("No access to adaptivity file!");

  if (Comm().MyPID() == 0)
  {
    (*logada_) << std::right << std::setw(9) << Step() << std::right << std::setw(16) << Time()
               << std::right << std::setw(16) << DtPast(1) << std::right << std::setw(16)
               << adaptstep_;

    // Who is responsible for the new time step size?
    (*logada_) << std::right << std::setw(16) << adareason_;

    // print L2-norms of the temporal discretization error indications
    (*logada_) << std::right << std::setw(16) << strnorm_ << std::right << std::setw(16)
               << strinnernorm_ << std::right << std::setw(16) << strfsinorm_ << std::right
               << std::setw(16) << flnorm_ << std::right << std::setw(16) << flinnernorm_
               << std::right << std::setw(16) << flfsinorm_;

    // print L-inf-norms of the temporal discretization error indications
    (*logada_) << std::right << std::setw(16) << strinfnorm_ << std::right << std::setw(16)
               << strinfinnernorm_ << std::right << std::setw(16) << strinffsinorm_ << std::right
               << std::setw(16) << flinfnorm_ << std::right << std::setw(16) << flinfinnernorm_
               << std::right << std::setw(16) << flinffsinorm_ << std::endl;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::PrintAdaptivitySummary() const
{
  if (Comm().MyPID() == 0)
  {
    if (Dt() != DtPast(1))  // only if time step size has changed
    {
      IO::cout << "\n"
               << "New time step size " << Dt() << " is based on " << adareason_ << ".\n";
    }

    if (dtminused_)
    {
      IO::cout << "Time step " << Step() << " has been done with minimum time step size."
               << " No further refinement possible. Go to next time step.\n";
    }

    if (not StepNotAccepted())
    {
      IO::cout << "Time step " << Step() << " has been accepted after " << adaptstep_ - 1
               << " repetitions."
               << "\n";
    }
    else
    {
      IO::cout << "Time step " << Step() << " will be repeated with dt = " << Dt() << ".\n";
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::TimeStepAuxiliar()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::TimeStepAuxiliar");

  // ---------------------------------------------------------------------------
  // Structure field
  // ---------------------------------------------------------------------------
  if (IsAdaStructure())
  {
    Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)
        ->TimeStepAuxiliar();
  }
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Fluid Field
  // ---------------------------------------------------------------------------
  if (IsAdaFluid())
  {
    FluidField()->TimeStepAuxiliar();
  }
  // ---------------------------------------------------------------------------

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::AdaptTimeStepSize()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::AdaptTimeStepSize");

  // Increment counter for repetition of time steps
  adaptstep_++;

  // prepare new time step size by copying the current one
  double dtnew = Dt();

  // ---------------------------------------------------------------------------
  // compute time step size suggestions based on structure field
  // ---------------------------------------------------------------------------
  if (IsAdaStructure())
  {
    Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)
        ->IndicateErrorNorms(
            strnorm_, strfsinorm_, strinnernorm_, strinfnorm_, strinffsinorm_, strinfinnernorm_);

    // structure based time step size suggestions
    dtstr_ = Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)
                 ->CalculateDt(strnorm_);
    dtstrfsi_ = Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)
                    ->CalculateDt(strfsinorm_);
    dtstrinner_ = Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)
                      ->CalculateDt(strinnernorm_);
  }
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // compute time step size suggestions based on fluid field
  // ---------------------------------------------------------------------------
  if (IsAdaFluid())
  {
    FluidField()->IndicateErrorNorms(
        flnorm_, flfsinorm_, flinnernorm_, flinfnorm_, flinffsinorm_, flinfinnernorm_);

    // error order
    const double order = FluidField()->GetTimAdaErrOrder();

    // calculate time step sizes resulting from errors in the fluid field
    dtfl_ = CalculateTimeStepSize(flnorm_, errtolfl_, order);
    dtflfsi_ = CalculateTimeStepSize(flfsinorm_, errtolfl_, order);
    dtflinner_ = CalculateTimeStepSize(flinnernorm_, errtolfl_, order);
  }
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // take care of possibly non-converged nonlinear solver
  // ---------------------------------------------------------------------------
  if (IsAdaSolver())
  {
    switch (erroraction_)
    {
      case FSI::Monolithic::erroraction_halve_step:
      {
        dtnonlinsolver_ = std::max(0.5 * Dt(), dtmin_);
        break;
      }
      case FSI::Monolithic::erroraction_revert_dt:
      {
        dtnonlinsolver_ = std::max(std::pow(0.95, (double)adaptstep_ - 1.0) * DtPast(0), dtmin_);
        break;
      }
      case FSI::Monolithic::erroraction_none:
      {
        dtnonlinsolver_ = std::min(std::max(1.1 * Dt(), dtmin_), dtmax_);
        break;
      }
      default:
        break;
    }
  }
  else
  {
    dtnonlinsolver_ = Dt();
  }
  // ---------------------------------------------------------------------------

  // Select new time step size
  dtnew = SelectDt();

  // optional averaging of time step size in case of increasing time step size
  if (dtnew > DtPast(1))
  {
    dtnew *= avgweights_[0];
    for (unsigned int i = 1; i < avgweights_.size(); ++i) dtnew += avgweights_[i] * DtPast(-i + 1);

    if (dtnew > dtmax_)
      dtnew = dtmax_;
    else if (dtnew < dtmin_)
      dtnew = dtmin_;
  }

  if (dtnew <= dtmin_ and DtPast(1) <= dtmin_) dtminused_ = true;

  // Who is responsible for changing the time step size?
  DetermineAdaReason(dtnew);

  // Check whether the step can be accepted, now
  accepted_ = SetAccepted();

  // take care of possibly non-converged nonlinear solver
  if (IsAdaSolver() and (erroraction_ == FSI::Monolithic::erroraction_halve_step or
                            erroraction_ == FSI::Monolithic::erroraction_revert_dt))
  {
    accepted_ = false;
  }

  // Finally, distribute new time step size to all participating fields/algorithms
  if (dtnew <= DtPast(1))  // time step size reduction is always done immediately
  {
    SetDt(dtnew);
  }
  else if (dtnew > DtPast(1))  // time step size increase may be delayed based on user's will
  {
    if (numincreasesteps_ <= 0)
    {
      // Finally, distribute new time step size to all participating fields/algorithms
      SetDt(dtnew);

      // reset
      numincreasesteps_ = DRT::Problem::Instance()
                              ->FSIDynamicParams()
                              .sublist("TIMEADAPTIVITY")
                              .get<int>("NUMINCREASESTEPS");
    }
    else if (numincreasesteps_ > 0)
    {
      numincreasesteps_--;
    }
  }
  else
  {
    dserror("Strange time step size!");
  }

  // ---------------------------------------------------------------------------
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::DetermineAdaReason(const double dt)
{
  const std::string oldreason = adareason_;

  if (IsAdaStructure() and
      (dt == dtstr_ or dt == dtstrfsi_ or dt == dtstrinner_))  // structure field?
    adareason_ = "Structure";
  else if (IsAdaFluid() and (dt == dtfl_ or dt == dtflfsi_ or dt == dtflinner_))  // fluid field?
    adareason_ = "Fluid";
  else if (IsAdaSolver() and dt == dtnonlinsolver_)  // Unconverged nonlinear solver?
    adareason_ = "Newton";

  // no change in time step size
  if (dt == DtPast(1)) adareason_ = oldreason;

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double FSI::Monolithic::CalculateTimeStepSize(
    const double errnorm, const double errtol, const double estorder) const
{
  //----------------------------------------------------------------------------
  // get some parameters first
  //----------------------------------------------------------------------------
  // FSI parameter list
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  // safety factor
  const double safetyfac = fsidyn.sublist("TIMEADAPTIVITY").get<double>("SAFETYFACTOR");

  // limiting factors for change of time step size
  const double facmax = fsidyn.sublist("TIMEADAPTIVITY").get<double>("SIZERATIOMAX");
  const double facmin = fsidyn.sublist("TIMEADAPTIVITY").get<double>("SIZERATIOMIN");

  // upper and lower bound of time step size
  const double dtmax = fsidyn.sublist("TIMEADAPTIVITY").get<double>("DTMAX");
  const double dtmin = fsidyn.sublist("TIMEADAPTIVITY").get<double>("DTMIN");
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  // Calculate scaling factor
  //----------------------------------------------------------------------------
  // the optimal scaling factor
  double fac = 1.0;
  if (not(errnorm == 0.0))  // do not divide by zero
    fac = std::pow(errtol / errnorm, 1.0 / (estorder + 1.0));
  else  // max increase if error norm == 0
    fac = facmax / safetyfac;

  // limit by safety factor <= 1.0
  fac *= safetyfac;

  //----------------------------------------------------------------------------
  // Calculate new time step size
  //----------------------------------------------------------------------------
  // optimal new step size
  double dtnew = fac * Dt();

  // redefine effective scaling factor to be dt*_{n}/dt_{n-1}, ie true optimal ratio
  fac = dtnew / DtPast(1);

  // limit effective scaling factor by maximum and minimum
  if (fac > facmax)
  {
    dtnew = facmax * DtPast(1);
  }
  else if (fac < facmin)
  {
    dtnew = facmin * DtPast(1);
  }

  // new step size subject to safety measurements
  if (dtnew > dtmax)
  {
    dtnew = dtmax;
  }
  else if (dtnew < dtmin)
  {
    dtnew = dtmin;
  }

  return dtnew;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::ResetStep()
{
  if (IsAdaStructure())
    Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)->ResetStep();
  else
    StructureField()->ResetStep();

  FluidField()->ResetStep();
  AleField()->ResetStep();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::ResetTime()
{
  // Fluid and ALE
  FluidField()->ResetTime(DtPast(1));
  AleField()->ResetTime(DtPast(1));

  // FSI routine
  SetTimeStep(Time() - DtPast(1), Step() - 1);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::SetDt(const double dtnew)
{
  // single fields
  if (IsAdaStructure())
    Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)->SetDt(dtnew);
  else
    StructureField()->SetDt(dtnew);

  FluidField()->SetDt(dtnew);
  AleField()->SetDt(dtnew);

  // FSI algorithm
  if (IsAdaStructure() or IsAdaFluid() or IsAdaSolver())
    dt_->SetStep(1, Dt());  // save step size of previous run of this time step for ResetTime()

  ADAPTER::AlgorithmBase::SetDt(dtnew);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double FSI::Monolithic::SelectDt() const
{
  // First, take the new time step size based on error estimation
  double dtnew = SelectDtErrorBased();

  // Consider convergence of nonlinear solver if the user wants to.
  if (IsAdaSolver())
  {
    // select time step size based on convergence of nonlinear solver
    if (erroraction_ != erroraction_none and erroraction_ != erroraction_continue)
    {
      dtnew = std::min(dtnew, dtnonlinsolver_);
    }
    else if (not(IsAdaStructure() or IsAdaFluid()))
    {
      dtnew = dtnonlinsolver_;
    }
  }
  else
  {
    // no change in time step size based on convergence of nonlinear solver
  }

  return dtnew;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool FSI::Monolithic::CheckIfDtsSame()
{
  // get time step sizes from all fields
  const double dtfsi = Dt();
  const double dtstruct = StructureField()->Dt();
  const double dtfluid = FluidField()->Dt();
  const double dtale = AleField()->Dt();

  double dtstrada = -1.0;
  if (IsAdaStructure())
    dtstrada =
        Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)->Dt();

  // print time step sizes in all fields
  if (Comm().MyPID() == 0)
  {
    std::cout << std::endl << "Time step sizes:" << std::endl;
    std::cout << "dt in FSI      : " << std::setprecision(16) << dtfsi << std::endl;
    std::cout << "dt in structure: " << std::setprecision(16) << dtstruct << std::endl;
    std::cout << "dt in fluid    : " << std::setprecision(16) << dtfluid << std::endl;
    std::cout << "dt in ale      : " << std::setprecision(16) << dtale << std::endl;

    if (IsAdaStructure())
      std::cout << "dt in str_ada  : " << std::setprecision(16) << dtstrada << std::endl;
  }

  // check whether time step sizes are all the same
  const double tol = 1.0e-12;
  if (dtfsi - dtstruct < tol and dtfsi - dtfluid < tol and dtfsi - dtale < tol)
  {
    if (IsAdaStructure() and dtfsi - dtstrada < tol)
      return true;
    else if (not IsAdaStructure())
      return true;
    else
      return false;
  }
  else
  {
    dserror("Time step sizes do not match among the fields.");
    return false;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
double FSI::Monolithic::DtPast(const int step) const { return (*dt_)[step]; }

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void FSI::Monolithic::UpdateDtPast(const double dtnew) { dt_->UpdateSteps(dtnew); }
