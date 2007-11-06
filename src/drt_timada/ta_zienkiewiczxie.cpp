/*======================================================================*/
/*!
\file ta_zienkiewczxie.cpp
\brief Generalized Alpha time integration for structural problems

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#ifdef CCADISCRET

#include "ta_zienkiewiczxie.H"
#include "iostream"


/*----------------------------------------------------------------------*/
/*!
\brief Constructor with parameters (public)
\author bborn
\date 10/07
*/
ZienkiewiczXie::ZienkiewiczXie
(
  double timeinitial,
  double timefinal,
  int timestepinitial,
  int timestepfinal,
  double stepsizeinitial,
  //
  double stepsizemax,
  double stepsizemin,
  double sizeratiomax,
  double sizeratiomin,
  double sizeratioscale,
  TAErrNorm errnorm,
  double errtol,
  int errorder,
  int adaptstepmax,
  //
  DRT::Discretization& dis,
  LINALG::Solver& solver,
  IO::DiscretizationWriter& output
)
: TimeAdaptivity
  (
    timeinitial,
    timefinal,
    timestepinitial,
    timestepfinal,
    stepsizeinitial,
    //
    stepsizemax,
    stepsizemin,
    sizeratiomax,
    sizeratiomin,
    sizeratioscale,
    errnorm,
    errtol,
    errorder,
    adaptstepmax,
    //
    dis,
    solver,
    output
  )
{
   return;
}

/*----------------------------------------------------------------------*/
/*!
\brief Destructor
\author bborn
\date 10/07
*/
ZienkiewiczXie::~ZienkiewiczXie()
{
   return;
}

/*----------------------------------------------------------------------*/
/*!
\brief Integrate in time
\author bborn
\date 10/07
*/
void ZienkiewiczXie::Integrate(StruGenAlpha::StruGenAlpha& timint)
{
  // --------------------------------------------------------------------
  // preparations
  // Generalised-alpha time integration scheme parameters
  double beta, gamma, alpham, alphaf;
  timint.GetTISPara(beta, gamma, alpham, alphaf);
  if (fabs(beta-1./6.) < 1.0e-6)
  {
    dserror("Generalised-alpha's beta must be non-equal to 1/6");
  }

  // ---------------------------------------------------------------------
  // initialise time loop
  double time_ = timeinitial_;
  int timestep_ = timestepinitial_;
  stepsize_ = stepsizeinitial_;
  stepsizepre_ = stepsize_;

  // ---------------------------------------------------------------------
  // time loop
  while ( (time_ < timefinal_) 
          && (timestep_ < timestepfinal_) )
  {
    // -------------------------------------------------------------------
    // time step size adapting loop
    adaptstep_ = 0;
    //double err = 2.0*errtol_;
    bool accepted = false;
    double stpsiznew;
    while ( (!accepted) && (adaptstep_ < adaptstepmax_) )
    {
      timint.SetTimeStepSize(stepsize_);
      timint.IntegrateStep();
      timint.ExtrapolateEndState();
      const RefCountPtr<Epetra_Vector>& acc = timint.GetAcc();
      const RefCountPtr<Epetra_Vector>& accn = timint.GetAccn();
      double factor = stepsize_*stepsize_*(1.0-6.0*beta)/6.0;
      locdiserrn_->Update(factor, *accn, -factor, *acc, 0.0);
      TimeAdaptivity::Indicate(accepted, stpsiznew);
      // modify step-size
      if (!accepted)
      {
        printf("Repeating step with stepsize = %g\n", stpsiznew);
        stepsize_ = stpsiznew;
      }
      adaptstep_ += 1;
    }  // end adaptive loop
    // -------------------------------------------------------------------
    // update or break
    if ( (mypid_ == 0) && (accepted) )
    {
      printf("Step size accepted\n");
    }
    else if ( (mypid_ == 0) && (adaptstep_ >= adaptstepmax_) )
    {
      printf("Could not find acceptable time step size ... continuing\n");
    }
    else
    {
      dserror("Do not know what to do");
    }
    timint.SetTime(time_);
    timint.SetTimeStep(timestep_);
    timint.SetTimeStepSize(stepsize_);
    timint.UpdateandOutput();
    timestep_ += 1;
    time_ += stepsize_;
    stepsizepre_ = stepsize_;
    stepsize_ = stpsiznew;
    cout << "Step " << timestep_ << ", Time " << time_ << ", Step size " << stepsize_ << endl;
  }  // end time loop

  // leave for good
  return;
}


#endif  // #ifdef CCADISCRET
