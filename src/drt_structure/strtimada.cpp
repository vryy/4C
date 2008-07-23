/*----------------------------------------------------------------------*/
/*!
\file strtimada.cpp

\brief Time step adaptivity front-end for structural dynamics

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include <iostream>

#include "strtimada.H"


/*----------------------------------------------------------------------*/
/* Constructor */
STR::StruTimAda::StruTimAda
(
  const Teuchos::ParameterList& sdyn,  //!< TIS input parameters
  const Teuchos::ParameterList& tap,  //!< adaptive input flags
  Teuchos::RCP<StruTimInt> tis  //!< marching time integrator
)
: sti_(tis),
  discret_(tis->Discretization()),
  mypid_(discret_->Comm().MyPID()),
  solver_(tis->GetSolver()),
  output_(tis->GetDiscretizationWriter()),
  //
  timeinitial_(0.0),
  timefinal_(sdyn.get<double>("MAXTIME")),
  timestepinitial_(0),
  timestepfinal_(sdyn.get<int>("NUMSTEP")),
  stepsizeinitial_(sdyn.get<double>("TIMESTEP")),
  //
  stepsizemax_(tap.get<double>("STEPSIZEMAX")),
  stepsizemin_(tap.get<double>("STEPSIZEMIN")),
  sizeratiomax_(tap.get<double>("SIZERATIOMAX")),
  sizeratiomin_(tap.get<double>("SIZERATIOMIN")),
  sizeratioscale_(tap.get<double>("SIZERATIOSCALE")),
  errctrl_(ctrl_dis),  // PROVIDE INPUT PARAMETER
  errnorm_(StruTimIntVector::MapNormStringToEnum(tap.get<std::string>("ERRNORM"))),
  errtol_(tap.get<double>("ERRTOL")),
  errorder_(1),  // CHANGE THIS CONSTANT
  adaptstepmax_(tap.get<int>("ADAPTSTEPMAX")),
  //
  time_(timeinitial_),
  timestep_(0),
  stepsizepre_(0),
  stepsize_(sdyn.get<double>("TIMESTEP")),
  locerrdisn_(Teuchos::null),
  adaptstep_(0)
{
  // allocate displacement local error vector
  locerrdisn_ = LINALG::CreateVector(*(discret_->DofRowMap()), true);

  // hallelujah
  return;
}

/*----------------------------------------------------------------------*/
/* Integrate adaptively in time */
void STR::StruTimAda::Integrate()
{
  // initialise time loop
  double time_ = timeinitial_;
  int timestep_ = timestepinitial_;
  stepsize_ = stepsizeinitial_;
  stepsizepre_ = stepsize_;

  // time loop
  while ( (time_ < timefinal_) and (timestep_ < timestepfinal_) )
  {
    // time step size adapting loop
    adaptstep_ = 0;
    //double err = 2.0*errtol_;
    bool accepted = false;
    double stpsiznew;
    while ( (not accepted) and (adaptstep_ < adaptstepmax_) )
    {
      // set current stepsize
      sti_->dt_->SetStep(0, stepsize_);
      //*(sti_->dt_(0)) = stepsize_;

      // integrate system with auxiliar TIS
      // we hold \f$D_{n+1}^{AUX}\f$ on #locdiserrn_
      // and \f$V_{n+1}^{AUX}\f$ on #locvelerrn_
      IntegrateStepAuxiliar();

      // integrate system with marching TIS and 
      sti_->IntegrateStep();

      // get local error vector on #locerrdisn_
      EvaluateLocalErrorDis();

      // check wether step passes
      Indicate(accepted, stpsiznew);

      // modify step-size
      if (not accepted)
      {
        printf("Repeating step with stepsize = %g\n", stpsiznew);
        stepsize_ = stpsiznew;
      }
      adaptstep_ += 1;
    }

    // update or break
    if ( (mypid_ == 0) and (accepted) )
    {
      printf("Step size accepted\n");
    }
    else if ( (mypid_ == 0) and (adaptstep_ >= adaptstepmax_) )
    {
      printf("Could not find acceptable time step size ... continuing\n");
    }
    else
    {
      dserror("Do not know what to do");
    }
    
    // sti_->time_ = time_ + stepsize_;
    sti_->time_->UpdateSteps(time_ + stepsize_);
    sti_->step_ = timestep_ + 1;
    // sti_->dt_ = stepsize_;
    sti_->dt_->UpdateSteps(stepsize_);

    sti_->UpdateStep();
    sti_->PrintStep();
    sti_->OutputStep();

    sti_->stepn_ = timestep_ += 1;
    sti_->timen_ = time_ += stepsize_;
    stepsizepre_ = stepsize_;
    stepsize_ = stpsiznew;
    
    std::cout << "Step " << timestep_ 
              << ", Time " << time_ 
              << ", Step size " << stepsize_ 
              << std::endl;
  }

  // leave for good
  return;
}

/*----------------------------------------------------------------------*/
/* Evaluate local error vector */
void STR::StruTimAda::EvaluateLocalErrorDis()
{
  // assumption: schemes do not have the same order of accuracy
  locerrdisn_->Update(-1.0, *(sti_->disn_), 1.0);
}

/*----------------------------------------------------------------------*/
/* Indicate error and determine new step size */
void STR::StruTimAda::Indicate
(
  bool& accepted,
  double& stpsiznew
)
{
  // norm of local discretisation error vector
  double norm = StruTimIntVector::CalculateNorm(errnorm_, locerrdisn_); 

  // check if acceptable
  accepted = (norm < errtol_);

  // debug
  std::cout << "ErrNorm " << norm 
            << ", ErrTol " << errtol_ 
            << ", Accept " << accepted 
            << std::endl;

  // get error order
  errorder_ = MethodOrderOfAccuracy();

  // optimal size ration with respect to given tolerance
  double sizrat = pow(errtol_/norm, 1.0/(errorder_+1.0));

  // debug
  printf("sizrat %g, stepsize %g, stepsizepre %g\n",
         sizrat, stepsize_, stepsizepre_);

  // scaled by safety parameter
  sizrat *= sizeratioscale_;
  // optimal new step size
  stpsiznew = sizrat * stepsize_;
  // redefine sizrat to be dt*_{n}/dt_{n-1}, ie true optimal ratio
  sizrat = stpsiznew/stepsizepre_;

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
  
  // get away from here
  return;
}

/*----------------------------------------------------------------------*/
/* Print constants */
void STR::StruTimAda::PrintConstants
(
  std::ostream& str
) const
{
  str << "StruTimAda:  Constants" << std::endl
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
      << "   Error norm = " << StruTimIntVector::MapNormEnumToString(errnorm_) << std::endl
      << "   Error order = " << errorder_ << std::endl
      << "   Error tolerance = " << errtol_ << std::endl
      << "   Max adaptive step = " << adaptstepmax_ << std::endl;
  return;
}

/*----------------------------------------------------------------------*/
/* Print variables */
void STR::StruTimAda::PrintVariables
(
  std::ostream& str
) const
{
  str << "StruTimAda:  Variables" << endl
      << "   Current time = " << time_ << endl
      << "   Previous step size = " << stepsizepre_ << endl
      << "   Current step size = " << stepsize_ << endl
      << "   Current adaptive step = " << adaptstep_ << endl;
  return;
}


/*----------------------------------------------------------------------*/
/* Print */
void STR::StruTimAda::Print
(
  std::ostream& str
) const
{
  str << "StruTimAda" << endl;
  PrintConstants(str);
  PrintVariables(str);
  // step aside
  return;
}


/*======================================================================*/
/* Out stream */
std::ostream& operator<<
(
  std::ostream& str, 
  const STR::StruTimAda::StruTimAda& ta
)
{
  ta.Print(str);
  return str;
}

/*----------------------------------------------------------------------*/
#endif  // #ifdef CCADISCRET
