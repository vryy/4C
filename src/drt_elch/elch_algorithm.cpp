
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
  :  FluidBaseAlgorithm(DRT::Problem::Instance()->FluidDynamicParams(),false),
     comm_(comm),
     step_(0),
     time_(0.0)
{
    // taking time loop control parameters out of fluid dynamics section
    const Teuchos::ParameterList& fluiddyn = DRT::Problem::Instance()->FluidDynamicParams();
    // maximum simulation time
    maxtime_=fluiddyn.get<double>("MAXTIME");
    // maximum number of timesteps
    nstep_ = fluiddyn.get<int>("NUMSTEP");
    
    // -------------------------------------------------------------------
    // get a vector layout from the fluid discretization to construct matching
    // vectors and matrices
    //                 local <-> global dof numbering
    // -------------------------------------------------------------------
    RCP<const Epetra_Map> velpresdofrowmap = FluidField().DofRowMap();
    
    // velocities at time n+1
    velnp_        = LINALG::CreateVector(*velpresdofrowmap,true);
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
  // time loop
  while (NotFinished())
  {
    // prepare next time step
    PrepareTimeStep();

    // solve Navier-Stokes system
    FluidField().NonlinearSolve();

    // get new velocity from Navier-Stokes solver
    velnp_=FluidField().Velnp();

    // solve con-dif-mig equation

    // update all field solvers
    Update();

    // write output to file
    Output();
  } // time loop

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

  // debug IO
#if 0
  // print out velocity-pressure vector
  cout<<*velnp_<<endl;
#endif
}


#endif // CCADISCRET
