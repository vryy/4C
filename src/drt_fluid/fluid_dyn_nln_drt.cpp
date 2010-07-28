/*!----------------------------------------------------------------------
\file fluid_dyn_nln_drt.cpp
\brief Main control routine for all fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     o generalized-alpha time-integration scheme

     and stationary solver.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

 *----------------------------------------------------------------------*/
#ifdef CCADISCRET
/*
#include <ctime>
#include <cstdlib>
#include <iostream>
*/

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "fluid_dyn_nln_drt.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/adapter_fluid_base_algorithm.H"

/*----------------------------------------------------------------------*
 * Main control routine for fluid including various solvers:
 *
 *        o instationary one-step-theta
 *        o instationary BDF2
 *        o instationary generalized-alpha (two versions)
 *        o stationary
 *
 *----------------------------------------------------------------------*/
void dyn_fluid_drt()
{
  // create a communicator
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // access to some parameter lists
  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  // create instance of fluid basis algorithm
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo = rcp(new ADAPTER::FluidBaseAlgorithm(fdyn,false));

  // read the restart information, set vectors and variables
  const int restart = probtype.get<int>("RESTART");
  if (restart) fluidalgo->FluidField().ReadRestart(restart);

  // run the simulation
  fluidalgo->FluidField().TimeLoop();

  // perform result tests if required
  DRT::Problem::Instance()->AddFieldTest(fluidalgo->FluidField().CreateFieldTest());
  DRT::Problem::Instance()->TestAll(comm);

  // have fun with your results!
  return;

} // end of dyn_fluid_drt()


#endif  // #ifdef CCADISCRET
