
#ifdef CCADISCRET

#include "elch_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
//#include "../drt_lib/drt_validparameters.H"
//#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::Algorithm(Epetra_Comm& comm)
  :  FluidBaseAlgorithm(DRT::Problem::Instance()->FluidDynamicParams()),
     comm_(comm),
     step_(0),
     time_(0.0)
{
    const Teuchos::ParameterList& fluiddyn = DRT::Problem::Instance()->FluidDynamicParams();
    // maximum simulation time
    maxtime_=fluiddyn.get<double>("MAXTIME");
    // maximum number of timesteps
    nstep_ = fluiddyn.get<int>("NUMSTEP");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::~Algorithm()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::TimeLoop()
{
  if (Comm().MyPID() == 0)
  cout<<"ELCH problemtype under development..."<<endl;

  while (step_< nstep_ and time_<maxtime_)
  {
    PrepareTimeStep();
    FluidField().NonlinearSolve();
    Update();
    Output();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;

  FluidField().PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::Update()
{
  FluidField().Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::Algorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  FluidField().Output();

  FluidField().LiftDrag();
}



#endif // CCADISCRET
