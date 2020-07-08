/*======================================================================*/
/*! \file
\brief Time step adaptivity front-end for structural dynamics
\level 1
*/

/*----------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------*/
/* headers */
#include <iostream>

#include "Epetra_Vector.h"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "../drt_inpar/inpar_structure.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"

#include "strtimada.H"
#include "strtimint.H"
#include "stru_aux.H"


/*----------------------------------------------------------------------*/
/* Constructor */
STR::TimAda::TimAda(const Teuchos::ParameterList& timeparams,  //!< TIS input parameters
    const Teuchos::ParameterList& tap,                         //!< adaptive input flags
    Teuchos::RCP<TimInt> tis                                   //!< marching time integrator
    )
    : sti_(tis),
      discret_(tis->Discretization()),
      myrank_(discret_->Comm().MyPID()),
      solver_(tis->Solver()),
      output_(tis->DiscWriter()),
      //
      timeinitial_(0.0),
      timefinal_(timeparams.get<double>("MAXTIME")),
      timedirect_(Sign(timefinal_ - timeinitial_)),
      timestepinitial_(0),
      timestepfinal_(timeparams.get<int>("NUMSTEP")),
      stepsizeinitial_(timeparams.get<double>("TIMESTEP")),
      //
      stepsizemax_(tap.get<double>("STEPSIZEMAX")),
      stepsizemin_(tap.get<double>("STEPSIZEMIN")),
      sizeratiomax_(tap.get<double>("SIZERATIOMAX")),
      sizeratiomin_(tap.get<double>("SIZERATIOMIN")),
      sizeratioscale_(tap.get<double>("SIZERATIOSCALE")),
      errctrl_(ctrl_dis),  // PROVIDE INPUT PARAMETER
      errnorm_(DRT::INPUT::IntegralValue<INPAR::STR::VectorNorm>(tap, "LOCERRNORM")),
      errtol_(tap.get<double>("LOCERRTOL")),
      errorder_(1),  // CHANGE THIS CONSTANT
      adaptstepmax_(tap.get<int>("ADAPTSTEPMAX")),
      //
      time_(timeinitial_),
      timestep_(0),
      stepsizepre_(stepsizeinitial_),
      stepsize_(timeparams.get<double>("TIMESTEP")),
      locerrdisn_(Teuchos::null),
      adaptstep_(0),
      //
      outsys_(false),
      outstr_(false),
      outene_(false),
      outrest_(false),
      outsysperiod_(tap.get<double>("OUTSYSPERIOD")),
      outstrperiod_(tap.get<double>("OUTSTRPERIOD")),
      outeneperiod_(tap.get<double>("OUTENEPERIOD")),
      outrestperiod_(tap.get<double>("OUTRESTPERIOD")),
      outsizeevery_(tap.get<int>("OUTSIZEEVERY")),
      outsystime_(timeinitial_ + outsysperiod_),
      outstrtime_(timeinitial_ + outstrperiod_),
      outenetime_(timeinitial_ + outeneperiod_),
      outresttime_(timeinitial_ + outrestperiod_),
      outsizefile_(Teuchos::null)
{
  // allocate displacement local error vector
  locerrdisn_ = LINALG::CreateVector(*(discret_->DofRowMap()), true);

  // check whether energyout_ file handle was attached
  if ((not sti_->AttachedEnergyFile()) and (outeneperiod_ != 0.0) and (myrank_ == 0))
  {
    sti_->AttachEnergyFile();
  }

  // check if step size file is wanted and attach
  if ((outsizeevery_ != 0) and (myrank_ == 0))
  {
    AttachFileStepSize();
  }

  // enable restart for adaptive timestepping - however initial timestep size is still read from
  // datfile! (mhv 01/2015)
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    // read restart of marching time-integrator and reset initial time and step for adaptive loop
    tis->ReadRestart(restart);
    timeinitial_ = tis->TimeOld();
    timestepinitial_ = tis->StepOld();

    // update variables which depend on initial time and step
    timedirect_ = Sign(timefinal_ - timeinitial_);
    outsystime_ = timeinitial_ + outsysperiod_;
    outstrtime_ = timeinitial_ + outstrperiod_;
    outenetime_ = timeinitial_ + outeneperiod_;
    outresttime_ = timeinitial_ + outrestperiod_;
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Integrate adaptively in time */
int STR::TimAda::Integrate()
{
  // finalize initialization
  // (only relevant if an auxiliary time integrator is used)
  Init(sti_);

  // Richardson extrapolation to no avail
  if (MethodAdaptDis() == ada_ident)
    dserror(
        "This combination is not implemented ... Richardson's extrapolation ... Yoshida technique "
        "...");

  // initialise time loop
  time_ = timeinitial_;
  timestep_ = timestepinitial_;
  stepsize_ = stepsizeinitial_;
  UpdateStepSize();

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
      sti_->dt_->SetStep(0, stepsize_);

      // integrate system with auxiliary TIS
      // we hold \f$D_{n+1}^{AUX}\f$ on #locdiserrn_
      // and \f$V_{n+1}^{AUX}\f$ on #locvelerrn_
      IntegrateStepAuxiliar();

      // integrate system with marching TIS and
      sti_->IntegrateStep();

      // get local error vector on #locerrdisn_
      EvaluateLocalErrorDis();

      // check whether step passes
      Indicate(accepted, stpsiznew);

      // adjust step-size and prepare repetition of current step
      if (not accepted)
      {
        IO::cout << "Repeating step with stepsize = " << stpsiznew << IO::endl;
        IO::cout << "- - - - - - - - - - - - - - - - - - - - - - - - -"
                 << " - - - - - - - - - - - - - - -" << IO::endl;

        stepsize_ = stpsiznew;

        ResetStep();
      }

      // increment number of adapted step sizes in a row
      adaptstep_ += 1;
    }

    // update or break
    if (accepted)
    {
      if (myrank_ == 0) std::cout << "Step size accepted" << std::endl;
    }
    else if (adaptstep_ >= adaptstepmax_)
    {
      if (myrank_ == 0)
        std::cout << "Could not find acceptable time step size"
                  << " ... continuing" << std::endl;
    }
    else
    {
      dserror("Do not know what to do");
    }

    // increment time and step in the marching time integrator
    sti_->time_->UpdateSteps(time_ + stepsize_);
    sti_->step_ = timestep_ + 1;
    sti_->dt_->UpdateSteps(stepsize_);

    // printing and output
    PrepareOutputPeriod();
    sti_->PreUpdate();
    sti_->UpdateStepState();
    sti_->UpdateStepElement();
    sti_->PostUpdate();
    sti_->PreOutput();
    OutputPeriod();
    sti_->PostOutput();
    OutputStepSize();
    sti_->PrintStep();

    // update
    //    Update();
    sti_->stepn_ = timestep_ += 1;
    sti_->timen_ = time_ += stepsize_;
    UpdateStepSize(stpsiznew);

    UpdatePeriod();
    outrest_ = outsys_ = outstr_ = outene_ = false;

    // the user reads but rarely listens
    if (myrank_ == 0)
    {
      std::cout << "Step " << timestep_ << ", Time " << time_ << ", StepSize " << stepsize_
                << std::endl;
    }
  }

  return 0;  // ToDo Provide meaningful error code here
}

/*----------------------------------------------------------------------*/
/* Evaluate local error vector */
void STR::TimAda::EvaluateLocalErrorDis()
{
  if (MethodAdaptDis() == ada_orderequal)
  {
    const double coeffmarch = sti_->MethodLinErrCoeffDis();
    const double coeffaux = MethodLinErrCoeffDis();
    locerrdisn_->Update(-1.0, *(sti_->disn_), 1.0);
    locerrdisn_->Scale(coeffmarch / (coeffaux - coeffmarch));
  }
  else
  {
    // schemes do not have the same order of accuracy
    locerrdisn_->Update(-1.0, *(sti_->disn_), 1.0);
  }

  // blank Dirichlet DOFs since they always carry the exact solution
  Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(locerrdisn_->Map(), true));
  LINALG::ApplyDirichlettoSystem(locerrdisn_, zeros, *(sti_->GetDBCMapExtractor()->CondMap()));
}

/*----------------------------------------------------------------------*/
/* Indicate error and determine new step size */
void STR::TimAda::Indicate(bool& accepted, double& stpsiznew)
{
  // norm of local discretisation error vector
  const int numneglect = sti_->GetDBCMapExtractor()->CondMap()->NumGlobalElements();
  const double norm = STR::AUX::CalculateVectorNorm(errnorm_, locerrdisn_, numneglect);

  // check if acceptable
  accepted = (norm < errtol_);

  // debug
  if (myrank_ == 0)
  {
    std::cout << "LocErrNorm " << std::scientific << norm << ", LocErrTol " << errtol_
              << ", Accept " << std::boolalpha << accepted << std::endl;
  }

  stpsiznew = CalculateDt(norm);

  return;
}

/*----------------------------------------------------------------------*/
/* Indicate error and determine new step size */
double STR::TimAda::CalculateDt(const double norm)
{
  // get error order
  if (MethodAdaptDis() == ada_upward)
    errorder_ = sti_->MethodOrderOfAccuracyDis();
  else
    errorder_ = MethodOrderOfAccuracyDis();

  // optimal size ration with respect to given tolerance
  double sizrat = 1.0;
  if (not(norm == 0.0))  // do not divide by zero
    sizrat = std::pow(errtol_ / norm, 1.0 / (errorder_ + 1.0));
  else  // max increase if error norm == 0
    sizrat = sizeratiomax_ / sizeratioscale_;

  // debug
  if (myrank_ == 0)
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
/* Prepare repetition of current time step */
void STR::TimAda::ResetStep()
{
  outrest_ = outsys_ = outstr_ = outene_ = false;
  sti_->ResetStep();

  return;
}

/*----------------------------------------------------------------------*/
/*  Modify step size to hit precisely output period */
void STR::TimAda::SizeForOutput()
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
/* Prepare output to file(s)                                            */
void STR::TimAda::PrepareOutputPeriod() { sti_->PrepareOutput(); }

/*----------------------------------------------------------------------*/
/* Output to file(s) */
void STR::TimAda::OutputPeriod()
{
  // this flag is passed along subroutines and prevents
  // repeated initialising of output writer, printing of
  // state vectors, or similar
  bool datawritten = false;

  // output restart (try this first)
  // write restart step
  if (outrest_)
  {
    sti_->OutputRestart(datawritten);
  }

  // output results (not necessary if restart in same step)
  if (outsys_ and (not datawritten))
  {
    sti_->OutputState(datawritten);
  }

  // output stress & strain
  if (outstr_)
  {
    sti_->OutputStressStrain(datawritten);
  }

  // output energy
  if (outene_)
  {
    sti_->OutputEnergy();
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Update output periods */
void STR::TimAda::UpdatePeriod()
{
  if (outrest_) outresttime_ += outrestperiod_;
  if (outsys_) outsystime_ += outsysperiod_;
  if (outstr_) outstrtime_ += outstrperiod_;
  if (outene_) outenetime_ += outeneperiod_;

  return;
}

/*----------------------------------------------------------------------*/
/* Write step size */
void STR::TimAda::OutputStepSize()
{
  if ((outsizeevery_ != 0) and (timestep_ % outsizeevery_ == 0) and (myrank_ == 0))
  {
    (*outsizefile_) << " " << std::setw(12) << timestep_ << std::scientific << std::setprecision(8)
                    << " " << time_ << " " << stepsize_ << " " << std::setw(2) << adaptstep_
                    << std::endl;
  }
}

/*----------------------------------------------------------------------*/
/* Print constants */
void STR::TimAda::PrintConstants(std::ostream& str) const
{
  str << "TimAda:  Constants" << std::endl
      << "   Initial time = " << timeinitial_ << std::endl
      << "   Final time = " << timefinal_ << std::endl
      << "   Initial Step = " << timestepinitial_ << std::endl
      << "   Final Step = " << timestepfinal_ << std::endl
      << "   Initial step size = " << stepsizeinitial_ << std::endl
      << "   Max step size = " << stepsizemax_ << std::endl
      << "   Min step size = " << stepsizemin_ << std::endl
      << "   Max size ratio = " << sizeratiomax_ << std::endl
      << "   Min size ratio = " << sizeratiomin_ << std::endl
      << "   Size ratio scale = " << sizeratioscale_ << std::endl
      << "   Error norm = " << INPAR::STR::VectorNormString(errnorm_) << std::endl
      << "   Error order = " << errorder_ << std::endl
      << "   Error tolerance = " << errtol_ << std::endl
      << "   Max adaptations = " << adaptstepmax_ << std::endl;
  return;
}

/*----------------------------------------------------------------------*/
/* Print variables */
void STR::TimAda::PrintVariables(std::ostream& str) const
{
  str << "TimAda:  Variables" << std::endl
      << "   Current time = " << time_ << std::endl
      << "   Previous step size = " << stepsizepre_ << std::endl
      << "   Current step size = " << stepsize_ << std::endl
      << "   Current adaptive step = " << adaptstep_ << std::endl;
  return;
}


/*----------------------------------------------------------------------*/
/* Print */
void STR::TimAda::Print(std::ostream& str) const
{
  str << "TimAda" << std::endl;
  PrintConstants(str);
  PrintVariables(str);

  return;
}

/*----------------------------------------------------------------------*/
/* Attach file handle for step size file #outsizefile_                  */
void STR::TimAda::AttachFileStepSize()
{
  if (outsizefile_.is_null())
  {
    std::string filename = DRT::Problem::Instance()->OutputControlFile()->FileName() + ".stepsize";
    outsizefile_ = Teuchos::rcp(new std::ofstream(filename.c_str()));
    (*outsizefile_) << "# timestep time step-size adaptations" << std::endl;
  }

  return;
}

/*----------------------------------------------------------------------*/
/* Update step size and set new time step size                          */
void STR::TimAda::UpdateStepSize(const double dtnew)
{
  UpdateStepSize();
  stepsize_ = dtnew;

  return;
}

/*----------------------------------------------------------------------*/
/* Update step size                                                     */
void STR::TimAda::UpdateStepSize()
{
  stepsizepre_ = stepsize_;

  return;
}

/*----------------------------------------------------------------------*/
/* Set new time step size                                               */
void STR::TimAda::SetDt(const double dtnew)
{
  stepsize_ = dtnew;   // in the adaptive time integrator
  sti_->SetDt(dtnew);  // in the marching time integrator
}

/*======================================================================*/
/* Out stream */
std::ostream& operator<<(std::ostream& str, const STR::TimAda& ta)
{
  ta.Print(str);

  return str;
}

/*----------------------------------------------------------------------*/
