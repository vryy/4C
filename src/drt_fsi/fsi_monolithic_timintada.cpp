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

#include "fsi_monolithic.H"

#include "../drt_inpar/inpar_fsi.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../linalg/linalg_utils.H"

#include "../drt_ale/ale.H"

#include "../drt_adapter/ad_str_timint_adaptive.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::InitTimIntAda(const Teuchos::ParameterList& fsidyn)
{
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

  errtolstr_ = fsidyn.sublist("TIMEADAPTIVITY").get<double>("LOCERRTOLSTRUCTURE");
  errtolfl_ = fsidyn.sublist("TIMEADAPTIVITY").get<double>("LOCERRTOLFLUID");

  strmethod_= DRT::INPUT::IntegralValue<int>(fsidyn.sublist("TIMEADAPTIVITY"),"AUXINTEGRATORSTRUCTURE");
  flmethod_= DRT::INPUT::IntegralValue<int>(fsidyn.sublist("TIMEADAPTIVITY"),"AUXINTEGRATORFLUID");

  switch (strmethod_)
  {
    case INPAR::FSI::timada_str_none:
    case INPAR::FSI::timada_str_limitdisinc:
    case INPAR::FSI::timada_str_expleuler:
    {
      estorderstr_ = 1.0;
      break;
    }
    case INPAR::FSI::timada_str_adamsbashforth2:
    case INPAR::FSI::timada_str_zienkiewicz_xie:
    {
      estorderstr_ = 2.0;
      break;
    }
    default:
    {
      dserror("Unknown auxiliary time integrator!");
      break;
    }
  }

  switch (flmethod_)
  {
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

  /* In case of Zienkiewicz-Xie error indicator in structure, we have to make
   * sure that the structure uses a Newmark-type time integrator. Otherwise, the
   * error indicator is not applicable. */
  if (strmethod_ == INPAR::FSI::timada_str_zienkiewicz_xie || strmethod_ == INPAR::FSI::timada_str_limitdisinc)
  {                                                           // ToDo LimitDisInc() should be available for all structural time integration schemes
    /* try to get the Newmark beta. This will lead to an dserror(), if the
     * structural time integrator is not of Newmark-type. */
    StructureField()->GetNewmarkBeta();
  }

  //----------------------------------------------------------------------------
  // Handling of Dirichlet BCs in error estimation
  //----------------------------------------------------------------------------
  // Create intersection of slave DOFs that hold a Dirichlet boundary condition
  // and are located at the FSI interface.
  std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmapsstruct;
  intersectionmapsstruct.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  intersectionmapsstruct.push_back(StructureField()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmapstruct = LINALG::MultiMapExtractor::IntersectMaps(intersectionmapsstruct);

  // Create intersection of fluid DOFs that hold a Dirichlet boundary condition
  // and are located at the FSI interface.
  std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmapsfluid;
  intersectionmapsfluid.push_back(FluidField().GetDBCMapExtractor()->CondMap());
  intersectionmapsfluid.push_back(FluidField().Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmapfluid = LINALG::MultiMapExtractor::IntersectMaps(intersectionmapsfluid);

  // store number of interface DOFs subject to Dirichlet BCs on structure and fluid side of the interface
  numstrfsidbcdofs_ = intersectionmapstruct->NumGlobalElements();
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

  PrepareTimeloop();

  // the time loop
  while (NotFinished())
  {
    PrepareAdaptiveTimeStep();

    // time step adaptivity loop
    while (StepNotAccepted() && adaptstep_ < adaptstepmax && (not dtminused_))
    {
      PrintHeaderRepeatedStep();
      PrepareTimeStep();
      TimeStep(interface);
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

  // if first time step, set initial time step size from input file
  if (Step() == 0)
  {
    // get initial time step size from input file
    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    const double dtinit = fsidyn.sublist("TIMEADAPTIVITY").get<double>("DTINITIAL");

    // Set initial time step size
    SetDt(dtinit);
  }

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
  if (adaptstep_ != 0 && Comm().MyPID() == 0)
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
  if (Comm().MyPID() == 0 && (not logada_.is_null()))
  {
    (*logada_) << "Time Adaptivity in monolithic Fluid-Structure-Interaction:" << "\n"
               << " - Error estimation method in fluid field: " << flmethod_  << "\n"
               << " - Error estimation method in structure field: " << strmethod_  << "\n"
               << "   (0 = None, 1 = Explicit Euler, 2 = Adams Bashforth 2)" << "\n"
               << " - Error tolerance in fluid field: " << errtolfl_  << "\n"
               << " - Error tolerance in structure field: " << errtolstr_  << "\n"
               << " - Minimum allowed time step size DTMIN = " << dtmin_ << "\n"
               << " - Maximum allowed time step size DTMAX = " << dtmax_ << "\n \n"
               << "  step | time | \t dt  | \t #adaptiter |  dt of field  |  err str field  |  err str inner |  err str fsi | err fl field  |  err fl inner |  err fl fsi \n\n"
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
    (*logada_)  << Step()
                << "\t" << Time()
                << "\t\t" << dtold_
                << "\t" << adaptstep_ ;

    // Who is responsible for the new time step size?
    (*logada_)  << "\t"<< adareason_;

    // print norms of the truncation error indications
    (*logada_)  << "\t "<< strnorm_
                << "\t" << strinnernorm_
                << "\t" << strfsinorm_
                << "\t" << flnorm_
                << "\t" << flinnernorm_
                << "\t" << flfsinorm_
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
      IO::cout << "Time step " << Step() << " has been accepted after " << adaptstep_ << " repetitions." << "\n";
    else
      IO::cout << "Time step " << Step() << " will be repeated with dt = " << Dt() << ".\n";
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::AdaptTimeStepSize()
{
  // Increment counter for repetition of time steps
  adaptstep_++;

  // Save old time step size for ResetTime()
  dtold_ = Dt();

  // prepare new time step size
  double dtnew = Dt();

  // ---------------------------------------------------------------------------

  // check whether truncation error based time adaptivity is activated for structure field
  if (not (strmethod_ == INPAR::FSI::timada_str_none || strmethod_ == INPAR::FSI::timada_str_limitdisinc))
  {
    // Indicate local truncation errors in structure
    IndicateLocalErrorNormsStructure();

    if (strmethod_ != INPAR::FSI::timada_str_none)
    {
      //calculate time step sizes resulting from errors in the structure field
      dtstr_ = CalculateTimeStepSize(strnorm_, errtolstr_, estorderstr_);
      dtstrfsi_ = CalculateTimeStepSize(strfsinorm_, errtolstr_, estorderstr_);
      dtstrinner_ = CalculateTimeStepSize(strinnernorm_, errtolstr_, estorderstr_);
    }
  }
  else if (strmethod_ == INPAR::FSI::timada_str_limitdisinc)
  {
    /* This approach adjusts the time step size such that the amount of change
     * in the displacement solution per time step is limited by a user given
     * value.
     *
     * Reference:
     * Givoli, Dan and Henigsberg, Ilan: A simple time-step control scheme,
     * Commun. Numer. Meth. Engng., Vol. 9, pp. 873-881, 1993
     *
     * It is still doubtful, if the approach is really useful and maintains accuracy.
     * Perhaps it has to be removed sometimes.
     *
     * mayr.mt 11/2013
     */

    // prepare norms
    double disnorm = 1.0e+12;
    double velnorm = 1.0e+12;
    double accnorm = 1.0e+12;

    // Calculate 'acceleration' measure
    Teuchos::RCP<Epetra_Vector> acc = Teuchos::rcp(new Epetra_Vector(*StructureField()->Accn()));
    acc->Update(2.0*StructureField()->GetNewmarkBeta(), *StructureField()->Accnp(), 1.0-2.0*StructureField()->GetNewmarkBeta());

    // compute norms
    StructureField()->Dispnp()->Norm2(&disnorm);
    StructureField()->Veln()->Norm2(&velnorm);
    acc->Norm2(&accnorm);

    velnorm *= Dt();
    accnorm *= (Dt() * Dt() / 2);

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
    // Calculate new time step size
    //----------------------------------------------------------------------
    const double facvel = errtolstr_ * disnorm / std::max(velnorm, 1.0e-10);
    const double facacc = sqrt(2.0*errtolstr_*disnorm / std::max(accnorm, 1.0e-10));

    dtnew = std::min(facvel, facacc);

    double fac = dtnew / Dt();

    // take care of not-a-number cases
    if (isnan(fac))
      fac = 1.0 / safetyfac;

    // modify scaling factor by safety factor in order to approach the desired
    // level of truncation error from below,  but not to hit the error level exactly.
    fac = fac * safetyfac;

    // limit factor by user given bounds
    fac = std::max(fac, facmin);
    fac = std::min(fac, facmax);

    // calculate the new time step size (limited by user given bounds)
    dtnew = std::max(fac * dtold_, dtmin);
    dtnew = std::min(dtnew, dtmax);
    //----------------------------------------------------------------------

    dtstr_ = dtnew;
    dtstrfsi_ = dtnew;
    dtstrinner_ = dtnew;
  }

  // check whether truncation error based time adaptivity is activated for fluid field
  if (not flmethod_ == INPAR::FSI::timada_fld_none)
  {
    // Indicate local truncation errors in fluid field
    IndicateLocalErrorNormsFluid();

    if (flmethod_ != INPAR::FSI::timada_fld_none)
    {
      //calculate time step sizes resulting from errors in the fluid field
      dtfl_ = CalculateTimeStepSize(flnorm_, errtolfl_, estorderfl_);
      dtflfsi_ = CalculateTimeStepSize(flfsinorm_, errtolfl_, estorderfl_);
      dtflinner_ = CalculateTimeStepSize(flinnernorm_, errtolfl_, estorderfl_);
    }
  }

  // ---------------------------------------------------------------------------

  // take care of possibly non-converged nonlinear solver
  if (erroraction_ == FSI::Monolithic::erroraction_halve_step && numhalvestep_ > 0)
  {
    dtnonlinsolver_ = std::max(Dt()/2, dtmin_);
  }
  else if (erroraction_ == FSI::Monolithic::erroraction_none) //&& numhalvestep_ > 0 )
  {
    dtnonlinsolver_ = std::min(std::max(1.1*Dt(), dtmin_), dtmax_);
    numhalvestep_ -= 1;
  }
  else
  {
    dtnonlinsolver_ = Dt();
  }

  // ---------------------------------------------------------------------------

  // Select new time step size
  dtnew = SelectTimeStepSize();

  if (dtnew <= dtmin_ && dtold_ <= dtmin_)
    dtminused_ = true;

  // Who is responsible for changing the time step size?
  DetermineAdaReason(dtnew);

  // Check whether the step can be accepted, now
  accepted_ = SetAccepted();

  // get safety factor from FSI parameter list
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const double safetyfac = fsidyn.sublist("TIMEADAPTIVITY").get<double>("SAFETYFACTOR");

  if (dtnew < safetyfac * Dt())
    accepted_ = false;

  // take care of possibly non-converged nonlinear solver
  if (erroraction_ == FSI::Monolithic::erroraction_halve_step)
    accepted_ = false;

  // Finally, distribute new time step size to all participating fields/algorithms
  SetDt(dtnew);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::DetermineAdaReason(const double dt)
{
  const std::string oldreason = adareason_;

  if (strmethod_ != INPAR::FSI::timada_str_none && (dt == dtstr_ || dt == dtstrfsi_ || dt == dtstrinner_)) // structure field?
    adareason_ = "Structure";
  else if (flmethod_ != INPAR::FSI::timada_fld_none && (dt == dtfl_ || dt == dtflfsi_ || dt == dtflinner_)) // fluid field?
    adareason_ = "Fluid";
  else if (dt == dtnonlinsolver_) // Unconverged nonlinear solver?
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
                                              )
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

  // prepare new time step size
  double dtnew = Dt(); //dtold_;

  // catch case that error norm is (very close to) zero
  if (errnorm < 1.0e-14)
  {
    dtnew = std::min(dtmax, facmax * dtold_);

    return dtnew;
  }

  //----------------------------------------------------------------------
  // Calculate new time step size
  //----------------------------------------------------------------------
  double fac = std::pow(errtol / errnorm, 1.0 / (estorder + 1.0));

  // take care of not-a-number cases
  if (isnan(fac))
    fac = 1.0 / safetyfac;

  // modify scaling factor by safety factor in order to approach the desired
  // level of truncation error from below,  but not to hit the error level exactly.
  fac = fac * safetyfac;

  // limit factor by user given bounds
  fac = std::min(std::max(fac, facmin), facmax);

  // calculate the new time step size (limited by user given bounds)
  dtnew = std::min(std::max(fac * Dt(), dtmin), dtmax); //dtold_
  //----------------------------------------------------------------------

  // return the new time step size
  return dtnew;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::ExplicitEuler(const Epetra_Vector& xn,
                                    const Epetra_Vector& xndot,
                                    Epetra_Vector& extrapolation
                                    ) const
{
  // Do a single explicit Euler step
  extrapolation.Update(1.0, xn, Dt(), xndot, 0.0);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::AdamsBashforth2(const Epetra_Vector& xn,
                                      const Epetra_Vector& xndot,
                                      const Epetra_Vector& xnmdot,
                                      Epetra_Vector& extrapolation
                                      ) const
{
  // time step sizes of current and previous time step
  const double dt = Dt();
  const double dto = dtprevious_;

  // Do a single Adams-Bashforth 2 step
  extrapolation.Update(1.0, xn, 0.0);
  extrapolation.Update((2.0*dt*dto+dt*dt) / (2*dto), xndot, -dt*dt / (2.0*dto), xnmdot, 1.0);

  return;
}

/*----------------------------------------------------------------------*/
void FSI::Monolithic::IndicateLocalErrorNormsFluid()
{
  // current state
  Teuchos::RCP<const Epetra_Vector> xn = Teuchos::rcp(new Epetra_Vector(*FluidField().Veln()));
  Teuchos::RCP<const Epetra_Vector> xnp = Teuchos::rcp(new Epetra_Vector(*FluidField().Velnp()));
  Teuchos::RCP<const Epetra_Vector> xndot = Teuchos::rcp(new Epetra_Vector(*FluidField().Accn()));

  // prepare vector for solution of auxiliary time step
  Teuchos::RCP<Epetra_Vector> extrapolation = Teuchos::rcp(new Epetra_Vector(*FluidField().DofRowMap(),true));

  // ---------------------------------------------------------------------------

  // calculate time step with auxiliary time integrator, i.e. the extrapolated solution
  switch (flmethod_)
  {
    case INPAR::FSI::timada_fld_none:
    {
      break;
    }
    case INPAR::FSI::timada_fld_expleuler:
    {
      ExplicitEuler(*xn, *xndot, *extrapolation);

      break;
    }
    case INPAR::FSI::timada_fld_adamsbashforth2:
    {
      if (Step() >= 1) // AdamsBashforth2 only if at least second time step
      {
        // Acceleration from previous time step
        Teuchos::RCP<Epetra_Vector> accnm = Teuchos::rcp(new Epetra_Vector(*FluidField().ExtractVelocityPart(FluidField().Accnm())));

        AdamsBashforth2(*xn, *xndot, *accnm, *extrapolation);
      }
      else // ExplicitEuler as starting algorithm
      {
        ExplicitEuler(*xn, *xndot, *extrapolation);
      }

      break;
    }
    default:
    {
      dserror("Unknown auxiliary time integration scheme for fluid field.");
      break;
    }
  }

  // ---------------------------------------------------------------------------

  // estimate the local error
  if (flmethod_ != INPAR::FSI::timada_fld_none)
  {
    Teuchos::RCP<Epetra_Vector> error = Teuchos::rcp(new Epetra_Vector(*FluidField().DofRowMap(),true));
    error->Update(1.0, *xnp, -1.0, *extrapolation, 0.0);

    /*
     * the complete vector 'error' contains all velocity and pressure DOFs
     * at the interface and in the interior domain, DBC and non-DBC DOFs.
     * For error indication we only want to use non-DBC and non-pressure DOFs,
     * i.e. only velocity DOFs which do not carry a Dirichlet boundary condition.
     */

    // set '0' on all pressure DOFs
    Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(error->Map(), true));
    LINALG::ApplyDirichlettoSystem(error, zeros, *(FluidField().PressureRowMap()));
    // TODO: Do not misuse ApplyDirichlettoSystem()...works for this purpose here: writes zeros into all pressure DoFs

    // set '0' on Dirichlet DOFs
    zeros = Teuchos::rcp(new Epetra_Vector(error->Map(), true));
    LINALG::ApplyDirichlettoSystem(error, zeros, *(FluidField().GetDBCMapExtractor()->CondMap()));

    // extrapolation of pressure not usable for error estimation --> extract velocity part of the
    //complete error vector i.e. interior and interface velocity DOFs
    Teuchos::RCP<const Epetra_Vector> errorvel
      = Teuchos::rcp(new Epetra_Vector(*FluidField().ExtractVelocityPart(error)));

    // extract error at interface velocity DOFs
    Teuchos::RCP<const Epetra_Vector> errorvelinterface
      = FluidField().Interface()->ExtractFSICondVector(error);

    // in case of fluid split: extract other part of the full error vector
    // i.e. all interior velocity DOFs and the pressure DOFs at the interface
    Teuchos::RCP<const Epetra_Vector> errorvelinterior
      = FluidField().Interface()->ExtractOtherVector(error);

    // Compute a suitable norm of the error vectors
    // Be careful: DBC and pressure DOFs carry zero, hence they do not contribute
    // to the norm and have to be neglected in the length scaling!
    errorvel->Norm2(&flnorm_);
    errorvelinterface->Norm2(&flfsinorm_);
    errorvelinterior->Norm2(&flinnernorm_);

    // -------------------------------------------------------------------------

    // Length scaling: We just need the non-DBC and the non-pressure DOFs.
    // Take care, that we do not divide by zero

    if (sqrt(errorvel->GlobalLength() - FluidField().GetDBCMapExtractor()->CondMap()->NumGlobalElements()) != 0)
      flnorm_ /= sqrt(errorvel->GlobalLength() - FluidField().GetDBCMapExtractor()->CondMap()->NumGlobalElements());
    else
      flnorm_ = 0.0;

    if (sqrt(errorvelinterface->GlobalLength() - numflfsidbcdofs_) != 0)
      flfsinorm_ /= sqrt(errorvelinterface->GlobalLength() - numflfsidbcdofs_);
    else
      flfsinorm_ = 0.0;

    if (sqrt(errorvelinterior->GlobalLength() - FluidField().PressureRowMap()->NumGlobalElements() - (FluidField().GetDBCMapExtractor()->CondMap()->NumGlobalElements() - numflfsidbcdofs_)) != 0)
      flinnernorm_ /= sqrt(errorvelinterior->GlobalLength() - FluidField().PressureRowMap()->NumGlobalElements() - (FluidField().GetDBCMapExtractor()->CondMap()->NumGlobalElements() - numflfsidbcdofs_));
    else
      flinnernorm_ = 0.0;

    // -------------------------------------------------------------------------
  }
  else
  {
    // nothing to do for the fluid field. Just set some default values.
    flnorm_ = 0.0;
    flfsinorm_ = 0.0;
    flinnernorm_ = 0.0;
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::IndicateLocalErrorNormsStructure()
{
  // current state
  Teuchos::RCP<const Epetra_Vector> xn = Teuchos::rcp(new Epetra_Vector(*StructureField()->Dispn()));
  Teuchos::RCP<const Epetra_Vector> xnp = Teuchos::rcp(new Epetra_Vector(*StructureField()->Dispnp()));
  Teuchos::RCP<const Epetra_Vector> xndot = Teuchos::rcp(new Epetra_Vector(*StructureField()->Veln()));

  // prepare vector for solution of auxiliary time step
  Teuchos::RCP<Epetra_Vector> extrapolation = Teuchos::rcp(new Epetra_Vector(*StructureField()->DofRowMap(),true));

  // ---------------------------------------------------------------------------

  // calculate time step with auxiliary time integrator, i.e. the extrapolated solution
  switch (strmethod_)
  {
    case INPAR::FSI::timada_str_none:
    case INPAR::FSI::timada_str_limitdisinc:
    {
      break;
    }
    case INPAR::FSI::timada_str_expleuler:
    {
      ExplicitEuler(*xn, *xndot, *extrapolation);

      break;
    }
    case INPAR::FSI::timada_str_adamsbashforth2:
    {
      if (Step() >= 1) // AdamsBashforth2 if at least second time step
      {
        // velocity from previous time step
        Teuchos::RCP<Epetra_Vector> velnm = Teuchos::rcp(new Epetra_Vector(*StructureField()->Velnm()));

        AdamsBashforth2(*xn, *xndot, *velnm, *extrapolation);
      }
      else // ExplicitEuler as starting algorithm
      {
        ExplicitEuler(*xn, *xndot, *extrapolation);
      }

      break;
    }
    case INPAR::FSI::timada_str_zienkiewicz_xie:
    {
      // No extrapolation needed for this error indicator
      break;
    }
    default:
    {
      dserror("Unknown auxiliary time integration scheme for structure field.");
      break;
    }
  }

  // ---------------------------------------------------------------------------

  // estimate the local error
  Teuchos::RCP<Epetra_Vector> error = Teuchos::rcp(new Epetra_Vector(*StructureField()->DofRowMap(),true));
  switch (strmethod_)
  {
    case INPAR::FSI::timada_str_none:
    case INPAR::FSI::timada_str_limitdisinc:
    {
      // do nothing
      break;
    }
    case INPAR::FSI::timada_str_expleuler:
    case INPAR::FSI::timada_str_adamsbashforth2:
    {
      error->Update(1.0, *xnp, -1.0, *extrapolation, 0.0);
      break;
    }
    case INPAR::FSI::timada_str_zienkiewicz_xie:
    {
      /* Reference:
       * Zienkiewicz, O. C. and Xie, Y. M: A simple error estimator and adaptive
       * time stepping procedure for dynamic analysis, Earthquake Engng. Struct. Dyn.
       * Vol. 20, pp. 871-887, 1991
       */
      error->Update(1.0, *StructureField()->Accnp(), -1.0, *StructureField()->Accn(), 0.0);
      error->Scale(Dt() * Dt() * (StructureField()->GetNewmarkBeta() - (1.0 / 6.0)));
      break;
    }
    default:
    {
      dserror("Unknown auxiliary time integration scheme for structure field.");
      break;
    }
  }

  // ---------------------------------------------------------------------------

  // compute error norms
  if (strmethod_ != INPAR::FSI::timada_str_none && strmethod_ != INPAR::FSI::timada_str_limitdisinc)
  {
    /*
     * the complete vector 'error' contains all displacement DOFs
     * at the interface and in the interior domain, DBC and non-DBC DOFs.
     * For error indication we only want to use non-DBC DOFs,
     * i.e. only displacement DOFs which do not carry a Dirichlet boundary condition.
     */

    // set '0' on Dirichlet DOFs
    Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(error->Map(), true));
    LINALG::ApplyDirichlettoSystem(error, zeros, *(StructureField()->GetDBCMapExtractor()->CondMap()));

    //extract the condition part of the full error vector (i.e. only interface displacement DOFs)
    Teuchos::RCP<const Epetra_Vector> errorinterface
      = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->ExtractFSICondVector(error)));

    //in case of structure split: extract the other part of the full error vector (i.e. only interior displacement DOFs)
    Teuchos::RCP<const Epetra_Vector> errorinterior
      = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->ExtractOtherVector(error)));

    // Calculate a suitable norm of the error vectors
    // DBC and pressure DOFs carry zero, hence they do not contribute to the norm.
    // But be careful with length scaling.
    error->Norm2(&strnorm_);
    errorinterface->Norm2(&strfsinorm_);
    errorinterior->Norm2(&strinnernorm_);

    // -------------------------------------------------------------------------

    // normalize the norms with the correct lengths of the error vectors
    // For length scaling we just need the non-DBC and non-pressure DOFs.
    // Take care, that we do not divide by zero
    if (sqrt(error->GlobalLength() - StructureField()->GetDBCMapExtractor()->CondMap()->NumGlobalElements()) != 0)
      strnorm_ /= sqrt(error->GlobalLength() - StructureField()->GetDBCMapExtractor()->CondMap()->NumGlobalElements());
    else
      strnorm_ = 0.0;

    if (sqrt(errorinterface->GlobalLength() - numstrfsidbcdofs_) != 0)
      strfsinorm_ /= sqrt(errorinterface->GlobalLength() - numstrfsidbcdofs_);
    else
      strfsinorm_ = 0.0;

    if (sqrt(errorinterior->GlobalLength() - (StructureField()->GetDBCMapExtractor()->CondMap()->NumGlobalElements() - numstrfsidbcdofs_)) != 0)
      strinnernorm_ /= sqrt(errorinterior->GlobalLength() - (StructureField()->GetDBCMapExtractor()->CondMap()->NumGlobalElements() - numstrfsidbcdofs_));
    else
      strinnernorm_ = 0.0;

    // -------------------------------------------------------------------------
  }
  else
  {
    // nothing to do for the structure field. Just set some default values.
    strnorm_ = 0.0;
    strfsinorm_ = 0.0;
    strinnernorm_ = 0.0;
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::ResetStep()
{
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
  StructureField()->SetDt(dtnew);
  FluidField().SetDt(dtnew);
  AleField().SetDt(dtnew);

  // FSI algorithm
  ADAPTER::AlgorithmBase::SetDt(dtnew);

  return;
}
