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
#ifdef TRILINOS_PACKAGE

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
:  TimeAdaptivity
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

   //int nstep = timint.params_.get<int>("nstep", 5);
   //int istep = timint.params_.get<int>("step", 0);

   // Generalised-alpha time integration scheme parameters
   double beta, gamma, alpham, alphaf;
   timint.GetTISPara(beta, gamma, alpham, alphaf);
   if (fabs(beta-1./6.) < 1.0e-6)
   {
      dserror("Generalised-alpha's beta must be non-equal to 1/6");
   }

   // initialise time loop
   double time_ = timeinitial_;
   int timestep_ = timestepinitial_;
   stepsize_ = stepsizeinitial_;

   // --------------------------------------------------------------------
   // time loop
   while ( (time_ < timefinal_) 
           && (timestep_ < timestepfinal_) )
   {

      // ----------
      // time step size adapting loop
      adaptstep_ = 0;
      double err = 2.0*errtol_;
      bool accepted = false;
      double stpsiznew;
      while ( (!accepted) && (adaptstep_ < adaptstepmax_) )
      {
         timint.SetTimeStepSize(stepsize_);
         timint.IntegrateStep();
         const RefCountPtr<Epetra_Vector>& acc = timint.GetAcc();
         const RefCountPtr<Epetra_Vector>& accn = timint.GetAccn();
         double factor = stepsize_*stepsize_*(1.0-6.0*beta)/6.0;
         locdiserrn_->Update(factor, *accn, -factor, *acc, 0.0);
         TimeAdaptivity::Indicate(accepted, stpsiznew);
         // modify step-size
         if (!accepted)
         {
            stepsize_ = stpsiznew;
         }
         adaptstep_ += 1;
      }
      // 
      if (accepted)
      {
         time_ += stepsize_;
         stepsize_ = stpsiznew;
      }
      else if (adaptstep_ >= adaptstepmax_)
      {
         printf("Could not find acceptable time step size ... continuing");
         time_ += stepsize_;
         stepsize_ = stpsiznew;
      }
      else
      {
         dserror("Do not know what to do");
      }
      // set total time via parameter list
      timint.UpdateandOutput();
   }

   return;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
