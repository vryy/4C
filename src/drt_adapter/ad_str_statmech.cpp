/*----------------------------------------------------------------------*/
/*!
\file ad_str_statmech.cpp

\brief Wrapper for the time loop of statmech problems

<pre>
Maintainer: Kei Müller
            mueller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>
*/
/*----------------------------------------------------------------------*/
#include <Teuchos_Time.hpp>

#include "../drt_statmech/statmech_manager.H"
#include "../drt_structure/strtimint_statmech.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include "ad_str_statmech.H"

//#define MEASURETIME

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::StructureStatMech::Integrate()
{
  Teuchos::RCP<STR::TimIntStatMech> statmechstructure = Teuchos::rcp_dynamic_cast<STR::TimIntStatMech>(structure_);
  //getting number of dimensions for diffusion coefficient calculation
  const int ndim= DRT::Problem::Instance()->NDim();

  // set statmech internal time and time step size
  double dt = 0.0;
  if(HaveStatMech())
    StatMechManager()->UpdateTimeAndStepSize(dt,TimeOld(),true);
  else if(HaveStatMechBilayer())
  {
    StatMechManagerBilayer()->UpdateTimeAndStepSize(dt,TimeOld(),true);
  }
  SetDt(dt);
  // this is necessary in case the time step size changed with the initial step (originally timen_ is set in strtimint.cpp)
  SetTimen(TimeOld()+Dt());

  // error checking variables
  INPAR::STR::ConvergenceStatus convergencestatus = INPAR::STR::conv_success;

  while(NotFinished())
  {
#ifdef MEASURETIME
    const double t0 = Teuchos::Time::wallTime();
#endif

    // preparations for statistical mechanics in this time step
    statmechstructure->StatMechPrepareStep();
#ifdef MEASURETIME
    const double t1 = Teuchos::Time::wallTime();
#endif

    //check if new random numbers are needed for a repeated time step
    bool newrandomnumbers=true;
    //redo time step in case of bad random configuration
    do
    {
      // Update of statmech specific quantities as well as new set of random numbers
      statmechstructure->StatMechUpdate(newrandomnumbers);

      // pay attention: for a constant predictor an incremental velocity update is necessary, which has
      // been deleted out of the code in order to simplify it
      // if(!Discretization()->Comm().MyPID())
      // std::cout<<"target time = "<<Time()<<", time step = "<<Dt()<<std::endl;
      statmechstructure->Predict();

      convergencestatus=Solve();

      if(convergencestatus!=INPAR::STR::conv_success)
        newrandomnumbers=statmechstructure->PerformErrorAction();
    }
    while(!statmechstructure->IsConverged());
#ifdef MEASURETIME
    const double t2 = Teuchos::Time::wallTime();
#endif

    //periodic shift of configuration at the end of the time step in order to avoid improper output
    if(HaveStatMech())
      StatMechManager()->PeriodicBoundaryShift(*WriteAccessDispnp(), ndim, Time(), Dt());

#ifdef MEASURETIME
    const double t3 = Teuchos::Time::wallTime();
#endif

    // update all that is relevant
    statmechstructure->UpdateAndOutput();

#ifdef MEASURETIME
    const double t4 = Teuchos::Time::wallTime();
#endif

    //special output for statistical mechanics
    statmechstructure->StatMechOutput();

#ifdef MEASURETIME
    if(!Discretization()->Comm().MyPID())
    {
      std::cout<<"\n=================Time  Measurement================"<<std::endl;
      std::cout<<"TimIntStatMech::Integrate"<<std::endl;
      std::cout<<"StatMechPrepareStep :\t"<<std::setprecision(4)<<t1-t0<<"\ts"<<std::endl;
      std::cout<<"Newton              :\t"<<std::setprecision(4)<<t2-t1<<"\ts"<<std::endl;
      std::cout<<"PeriodicShift       :\t"<<std::setprecision(4)<<t3-t2<<"\ts"<<std::endl;
      std::cout<<"UpdateAndOutput     :\t"<<std::setprecision(4)<<t4-t3<<"\ts"<<std::endl;
      std::cout<<"StatMechOutput      :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t4<<"\ts"<<std::endl;
      std::cout<<"=================================================="<<std::endl;
      std::cout<<"total time         :\t"<<std::setprecision(4)<<Teuchos::Time::wallTime()-t0<<"\ts"<<std::endl;
    }
#endif
  }
  return 0;
} // void STR::TimIntStatMech::Integrate()

