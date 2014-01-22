/*----------------------------------------------------------------------*/
/*!
\file fsi_monolithic_timintada.cpp

\brief Provide time adaptivity functionalities within the Monolithic class

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "fsi_monolithic.H"

#include "../drt_inpar/inpar_fsi.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_utils.H"

#include "../drt_ale/ale.H"

#include "../drt_adapter/ad_str_fsi_timint_adaptive.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::InitTimIntAda(const Teuchos::ParameterList& fsidyn)
{
  // access to structural time adaptivity parameters
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& sada = sdyn.sublist("TIMEADAPTIVITY");

  // initialize member variables
  dtold_ = Dt();
  dtprevious_ = Dt();

  adaptstep_ = 0;
  accepted_ = false;

  adareason_ = "none";

  strnorm_ = 0.0;
  flnorm_ = 0.0;
  strfsinorm_ = 0.0;
  flfsinorm_ = 0.0;
  strinnernorm_ = 0.0;
  flinnernorm_ = 0.0;

  dtstr_ = 0.0;
  dtfl_ = 0.0;
  dtstrfsi_ = 0.0;
  dtflfsi_ = 0.0;
  dtstrinner_ = 0.0;
  dtflinner_ = 0.0;
  dtnonlinsolver_ = 0.0;

  dtminused_ = false;

  //----------------------------------------------------------------------------
  // write adaptivity file
  //----------------------------------------------------------------------------
  std::string fileada = DRT::Problem::Instance()->OutputControlFile()->FileName();
  fileada.append(".adaptivity");
  logada_ = Teuchos::rcp(new std::ofstream(fileada.c_str()));

  //----------------------------------------------------------------------------
  // set algorithmic parameters
  //----------------------------------------------------------------------------
  dtmax_ = fsidyn.sublist("TIMEADAPTIVITY").get<double>("DTMAX");
  dtmin_ = fsidyn.sublist("TIMEADAPTIVITY").get<double>("DTMIN");

  errtolstr_ = sada.get<double>("LOCERRTOL");
  errtolfl_ = fsidyn.sublist("TIMEADAPTIVITY").get<double>("LOCERRTOLFLUID");

  //----------------------------------------------------------------------------
  // check on which fields time adaptivity should be based on
  //----------------------------------------------------------------------------
  if (not DRT::INPUT::IntegralValue<INPAR::STR::TimAdaKind>(sada,"KIND") == INPAR::STR::timada_kind_none)
    isadastructure_ = true;

  if (not (DRT::INPUT::IntegralValue<int>(fsidyn.sublist("TIMEADAPTIVITY"), "AUXINTEGRATORFLUID") == INPAR::FSI::timada_fld_none))
    isadafluid_ = true;

  if (DRT::INPUT::IntegralValue<int>(fsidyn.sublist("TIMEADAPTIVITY"), "DIVERCONT") == INPAR::FSI::divcont_halve_step)
    isadasolver_ = true;

  flmethod_ = DRT::INPUT::IntegralValue<int>(fsidyn.sublist("TIMEADAPTIVITY"),"AUXINTEGRATORFLUID");

  switch (flmethod_) // ToDo: depends also on the order of accuracy of marching time integrator
  {                  //       (cf. Thesis by B. Bornemann, 2003 and structure field)
    case INPAR::FSI::timada_fld_none:
    case INPAR::FSI::timada_fld_expleuler:
    {
      estorderfl_ = 1.0;
      break;
    }
    case INPAR::FSI::timada_fld_adamsbashforth2:
    {
      estorderfl_ = 2.0;
      break;
    }
    default:
    {
      dserror("Unknown auxiliary time integrator!");
      break;
    }
  }

  //----------------------------------------------------------------------------
  // Handling of Dirichlet BCs in error estimation
  //----------------------------------------------------------------------------
  // Create intersection of fluid DOFs that hold a Dirichlet boundary condition
  // and are located at the FSI interface.
  std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmapsfluid;
  intersectionmapsfluid.push_back(FluidField().GetDBCMapExtractor()->CondMap());
  intersectionmapsfluid.push_back(FluidField().Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmapfluid = LINALG::MultiMapExtractor::IntersectMaps(intersectionmapsfluid);

  // store number of interface DOFs subject to Dirichlet BCs on structure and fluid side of the interface
  numflfsidbcdofs_ = intersectionmapfluid->NumGlobalElements();

  //----------------------------------------------------------------------------
  // Check whether input parameters make sense
  //---------------------------------------------------------------------------
  const double safetyfac = fsidyn.sublist("TIMEADAPTIVITY").get<double>("SAFETYFACTOR");
  if (safetyfac > 1.0)
    dserror("SAFETYFACTOR in FSI DYNAMIC/TIMEADAPTIVITY is %f > 1.0 and, thus, to large.", safetyfac);

  if (dtmax_ <= dtmin_)
    dserror("DTMAX = %f has to be larger than DMIN = %f.", dtmax_, dtmin_);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::TimeloopAdaDt(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  /*----------------------------------------------------------------------*/
  /* Initialize some parameters                                           */
  /*----------------------------------------------------------------------*/
  // get the FSI parameter list
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  // number of maximum allowed time step repetitions
  const int adaptstepmax = fsidyn.sublist("TIMEADAPTIVITY").get<int>("ADAPTSTEPMAX");
  /*----------------------------------------------------------------------*/

  //resize MStep objects of structure, needed for AB2 structural extrapolation of displacements
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
        IO::cout<< "Time Step has already been calculated with minimum step size"
                << " --> continue with next time step!"
                << "\n";
        accepted_ = true;
      }

      if (StepNotAccepted())
      {
        if (adaptstep_ >= adaptstepmax)
        {
          if (Comm().MyPID() == 0)
          {
            IO::cout <<"adaptstep_ = "<< adaptstepmax
                     << " --> Max. number of adaption iterations is reached: continuing with next time step! "
                     << "\n";
          }
        }
        else
        {
          ResetStep();
          ResetTime();

          if (Comm().MyPID() == 0)
          {
            IO::cout << "Repeat current time step with dt = " << Dt()
                     << " based on " << adareason_ << ".\n";
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
    Update();
    Output();

    WriteAdaFile();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::PrepareAdaptiveTimeStep()
{
  // reset from last time step
  dtminused_ = false;
  adaptstep_ = 0;
  accepted_ = false;

  // store time step size of previous time step for multi-step auxiliary time integrators
  dtprevious_ = Dt();

  if (Comm().MyPID() == 0)
  {
    IO::cout << "\n"
             << "+++++++++++++++++++++++++++++NEW TIME STEP+++++++++++++++++++++++++++++";
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::PrintHeaderRepeatedStep() const
{
  if (adaptstep_ != 0 and Comm().MyPID() == 0)
  {
    IO::cout << "__________REAPEATING TIME STEP " << Step()
             << " WITH DT = " << Dt()
             << " FOR THE " << adaptstep_
             << ". TIME__________"
             << "\n";
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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
    (*logada_) << "Time Adaptivity in monolithic Fluid-Structure-Interaction:" << "\n"
               << " - Error estimation method in fluid field: " << flmethod_  << "\n"
               << " - Error estimation method in structure field: " << strmethod << "\n"
               << "   (0 = None, 1 = Explicit Euler, 2 = Adams Bashforth 2)" << "\n"
               << " - Error tolerance in fluid field: " << errtolfl_  << "\n"
               << " - Error tolerance in structure field: " << errtolstr_  << "\n"
               << " - Minimum allowed time step size DTMIN = " << dtmin_ << "\n"
               << " - Maximum allowed time step size DTMAX = " << dtmax_ << "\n \n"
               << std::right << std::setw(9) << "step"
               << std::right << std::setw(16) << "time"
               << std::right << std::setw(16) << "dt"
               << std::right << std::setw(16) << "#adaptiter"
               << std::right << std::setw(16) << "dt of field"
               << std::right << std::setw(16) << "err str field"
               << std::right << std::setw(16) << "err str inner"
               << std::right << std::setw(16) << "err str fsi"
               << std::right << std::setw(16) << "err fl field"
               << std::right << std::setw(16) << "err fl inner"
               << std::right << std::setw(16) << "err fl fsi"
               << "\n\n"
      ;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::WriteAdaFile() const
{
  if (logada_.is_null())
    dserror("No access to adaptivity file!");

  if (Comm().MyPID() == 0)
  {
    (*logada_)  << std::right << std::setw(9) << Step()
                << std::right << std::setw(16) << Time()
                << std::right << std::setw(16) << dtprevious_ //dtold_ //ToDo Which one is the right one???
                << std::right << std::setw(16) << adaptstep_ ;

    // Who is responsible for the new time step size?
    (*logada_)  << std::right << std::setw(16) << adareason_;

    // print norms of the truncation error indications
    (*logada_)  << std::right << std::setw(16) << strnorm_
                << std::right << std::setw(16) << strinnernorm_
                << std::right << std::setw(16) << strfsinorm_
                << std::right << std::setw(16) << flnorm_
                << std::right << std::setw(16) << flinnernorm_
                << std::right << std::setw(16) << flfsinorm_
                << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::PrintAdaptivitySummary() const
{

  if (Comm().MyPID() == 0)
  {
    if (Dt() != dtold_) // only if time step size has changed
    {
      IO::cout << "\n" << "New time step size " << Dt()
               << " is based on " << adareason_ << ".\n";
    }

    if (dtminused_)
    {
      IO::cout << "Time step " << Step()
               << " has been done with minimum time step size."
               << " No further refinement possible. Go to next time step.\n";
    }

    if (not StepNotAccepted())
      IO::cout << "Time step " << Step() << " has been accepted after " << adaptstep_-1 << " repetitions." << "\n";
    else
      IO::cout << "Time step " << Step() << " will be repeated with dt = " << Dt() << ".\n";
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::TimeStepAuxiliar()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::TimeStepAuxiliar");

  // -----------------------------------------------------------------------
  // Structure field
  // -----------------------------------------------------------------------
  if (IsAdaStructure())
  {
    Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)
        ->TimeStepAuxiliar();
  }
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // Fluid Field
  // -----------------------------------------------------------------------
  if (IsAdaFluid())
  {
    FluidField().TimeStepAuxiliar();
  }
  // -----------------------------------------------------------------------

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::AdaptTimeStepSize()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::AdaptTimeStepSize");

  // Increment counter for repetition of time steps
  adaptstep_++;

  // Save time step size of previous run of this time step for ResetTime()
  dtold_ = Dt();

  // prepare new time step size by copying the current one
  double dtnew = Dt();

  // -----------------------------------------------------------------------
  // compute time step size suggestions based on structure field
  // -----------------------------------------------------------------------
  if (IsAdaStructure())
  {
    Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)
        ->IndicateErrorNorms(strnorm_, strfsinorm_, strinnernorm_);

    // structure based time step size suggestions
    dtstr_ = Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)->CalculateDt(strnorm_);
    dtstrfsi_ = Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)->CalculateDt(strfsinorm_);
    dtstrinner_ = Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)->CalculateDt(strinnernorm_);
  }
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // compute time step size suggestions based on fluid field
  // -----------------------------------------------------------------------
  if (IsAdaFluid())
  {
    FluidField().IndicateErrorNorms(flnorm_, flfsinorm_, flinnernorm_);

    //calculate time step sizes resulting from errors in the fluid field
    dtfl_ = CalculateTimeStepSize(flnorm_, errtolfl_, estorderfl_);
    dtflfsi_ = CalculateTimeStepSize(flfsinorm_, errtolfl_, estorderfl_);
    dtflinner_ = CalculateTimeStepSize(flinnernorm_, errtolfl_, estorderfl_);
  }
  // -----------------------------------------------------------------------

  // -----------------------------------------------------------------------
  // take care of possibly non-converged nonlinear solver
  // -----------------------------------------------------------------------
  if (IsAdaSolver() and erroraction_ == FSI::Monolithic::erroraction_halve_step and numhalvestep_ > 0)
  {
    dtnonlinsolver_ = std::max(Dt()/2, dtmin_);
  }
  else if (IsAdaSolver() and erroraction_ == FSI::Monolithic::erroraction_none) //and numhalvestep_ > 0 )
  {
    dtnonlinsolver_ = std::min(std::max(1.1*Dt(), dtmin_), dtmax_);
    numhalvestep_ -= 1;
  }
  else
  {
    dtnonlinsolver_ = Dt();
  }
  // -----------------------------------------------------------------------

  // Select new time step size
  dtnew = SelectTimeStepSize();

  if (dtnew <= dtmin_ and dtold_ <= dtmin_)
    dtminused_ = true;

  // Who is responsible for changing the time step size?
  DetermineAdaReason(dtnew);

  // Check whether the step can be accepted, now
  accepted_ = SetAccepted();

  // take care of possibly non-converged nonlinear solver
  if (IsAdaSolver() and erroraction_ == FSI::Monolithic::erroraction_halve_step)
    accepted_ = false;

  // Finally, distribute new time step size to all participating fields/algorithms
  SetDt(dtnew);

  // -----------------------------------------------------------------------

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::DetermineAdaReason(const double dt)
{
  const std::string oldreason = adareason_;

  if (IsAdaStructure() and (dt == dtstr_ or dt == dtstrfsi_ or dt == dtstrinner_)) // structure field?
    adareason_ = "Structure";
  else if (IsAdaFluid() and (dt == dtfl_ or dt == dtflfsi_ or dt == dtflinner_)) // fluid field?
    adareason_ = "Fluid";
  else if (IsAdaSolver() and dt == dtnonlinsolver_) // Unconverged nonlinear solver?
    adareason_ = "Newton";

  // no change in time step size
  if (dt == dtold_)
    adareason_ = oldreason;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::Monolithic::CalculateTimeStepSize(const double errnorm,
                                              const double errtol,
                                              const double estorder
                                              ) const
{
  //----------------------------------------------------------------------
  // get some parameters first
  //----------------------------------------------------------------------
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
  //----------------------------------------------------------------------

  //----------------------------------------------------------------------
  // Calculate scaling factor
  //----------------------------------------------------------------------
  // the optimal scaling factor
  double fac = 1.0;
  if (not (errnorm == 0.0)) // do not divide by zero
    fac = std::pow(errtol / errnorm, 1.0 / (estorder + 1.0));
  else // max increase if error norm == 0
    fac = facmax / safetyfac;

  // limit by safety factor <= 1.0
  fac *= safetyfac;

  //----------------------------------------------------------------------
  // Calculate new time step size
  //----------------------------------------------------------------------
  // optimal new step size
  double dtnew = fac * Dt();

  // redefine effective scaling factor to be dt*_{n}/dt_{n-1}, ie true optimal ratio
  fac = dtnew / dtprevious_;

  // limit effective scaling factor by maximum and minimum
  if (fac > facmax)
  {
    dtnew = facmax * dtprevious_;
  }
  else if (fac < facmin)
  {
    dtnew = facmin * dtprevious_;
  }

  // new step size subject to safety measurements
  if (dtnew > dtmax)
  {
    dtnew = dtmax;
  }
  else if (dtnew < dtmin)
  {
    dtnew  = dtmin;
  }

  return dtnew;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::ResetStep()
{

  if (IsAdaStructure())
    Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)->ResetStep();
  else
    StructureField()->ResetStep();

  FluidField().ResetStep();
  AleField().ResetStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::ResetTime()
{
  // Fluid and ALE
  FluidField().ResetTime(dtold_);
  AleField().ResetTime(dtold_);

  // FSI routine
  SetTimeStep(Time() - dtold_, Step() - 1);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::SetDt(const double dtnew)
{
  // single fields
  if (IsAdaStructure())
  {
    Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)->SetDt(dtnew); //UpdateStepSize(dtnew);
//    Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)->UpdateStepSize(dtnew);
  }
  else
    StructureField()->SetDt(dtnew);

  FluidField().SetDt(dtnew);
  AleField().SetDt(dtnew);

  // FSI algorithm
  ADAPTER::AlgorithmBase::SetDt(dtnew);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::Monolithic::CheckIfDtsSame()
{
  // get time step sizes from all fields
  const double dtfsi = Dt();
  const double dtstruct = StructureField()->Dt();
  const double dtfluid = FluidField().Dt();
  const double dtale = AleField().Dt();

  double dtstrada = -1.0;
  if (IsAdaStructure())
    dtstrada = Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)->Dt();

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

  //check whether time step sizes are all the same
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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::Monolithic::CheckIfTimesSame()
{
  // get time from all fields
  const double tfsi = Time();
  const double tstruct = StructureField()->GetTimeNew();
  const double tfluid = FluidField().Time();
  const double tale = AleField().Time();

  double tstrada = tstruct;
  if (IsAdaStructure())
    tstrada = Teuchos::rcp_dynamic_cast<ADAPTER::StructureFSITimIntAda>(StructureField(), true)->GetTimeNew();

  // print time in all fields
  if (Comm().MyPID() == 0)
  {
    std::cout << std::endl << "Time:" << std::endl;
    std::cout << "time in FSI      : " << std::setprecision(16) << tfsi << std::endl;
    std::cout << "time in structure: " << std::setprecision(16) << tstruct << std::endl;
    std::cout << "time in fluid    : " << std::setprecision(16) << tfluid << std::endl;
    std::cout << "time in ale      : " << std::setprecision(16) << tale << std::endl;

    if (IsAdaStructure())
      std::cout << "time in str_ada  : " << std::setprecision(16) << tstrada << std::endl;
  }

  //check whether times are all the same
  const double tol = 1.0e-12;
  if (tfsi - tstruct < tol and tfsi - tfluid < tol and tfsi - tale < tol)
  {
    if (IsAdaStructure() and tfsi - tstrada < tol)
      return true;
    else if (not IsAdaStructure())
      return true;
    else
      return false;
  }
  else
  {
    dserror("Time does not match among the fields.");
    return false;
  }
}
