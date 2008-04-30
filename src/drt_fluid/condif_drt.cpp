/*!----------------------------------------------------------------------
\file condif_drt.cpp
\brief Main control routine for all (in)stationary convect.-diff. solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     o generalized-alpha time-integration scheme

     and stationary solver.

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "condif_drt.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/adapter_condif_base_algorithm.H"


/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*
 * Main control routine for convection-diffusion incl. various solvers:
 *
 *        o instationary one-step-theta
 *        o instationary BDF2
 *        o instationary generalized-alpha
 *        o stationary
 *
 *----------------------------------------------------------------------*/
void dyn_condif_drt()
{
  // create instance of convection diffusion basis algorithm
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();
  Teuchos::RCP<ADAPTER::ConDifBaseAlgorithm> condifonly = rcp(new ADAPTER::ConDifBaseAlgorithm(fdyn)); 

  // set velocity field
  int veltype = Teuchos::getIntegralValue<int>(fdyn,"CD_VELOCITY");
  switch (veltype)
  {
    case 0:  // zero  (see case 1)
    case 1:  // function
      (condifonly->ConDifField()).SetVelocityField(veltype,fdyn.get<int>("CD_VELFUNCNO"));
      break;
    case 2:  // Navier_Stokes
      dserror("condif velocity: >>Navier_Stokes<< not implemented.");
    default:
      dserror("unknown velocity field type for convection-diffusion");
  }

  // solve the convection-diffusion problem
  (condifonly->ConDifField()).Integrate();

  return;

} // end of dyn_condif_drt()

#endif  // #ifdef CCADISCRET
