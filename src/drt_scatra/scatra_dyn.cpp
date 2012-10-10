/*----------------------------------------------------------------------*/
/*!
\file scatra_dyn.cpp
\brief entry point for (passive) scalar transport problems

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_dyn.H"
#include "passive_scatra_algorithm.H"
#include "scatra_utils_clonestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "scatra_resulttest.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <iostream>


/*----------------------------------------------------------------------*
 * Main control routine for scalar transport problems, incl. various solvers
 *
 *        o Laplace-/ Poisson equation (zero velocity field)
 *          (with linear and nonlinear boundary conditions)
 *        o transport of passive scalar in velocity field given by spatial function
 *        o transport of passive scalar in velocity field given by Navier-Stokes
 *          (one-way coupling)
 *
 *----------------------------------------------------------------------*/
void scatra_dyn(int restart)
{
  // access the communicator
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();

  // access the fluid discretization
  RefCountPtr<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis("fluid");
  // access the scatra discretization
  RefCountPtr<DRT::Discretization> scatradis = DRT::Problem::Instance()->GetDis("scatra");

  // ensure that all dofs are assigned in the right order; this creates dof numbers with
  //       fluid dof < scatra dof
  fluiddis->FillComplete();
  scatradis->FillComplete();

  // set velocity field
  const INPAR::SCATRA::VelocityField veltype
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn,"VELOCITYFIELD");
  switch (veltype)
  {
    case INPAR::SCATRA::velocity_zero:  // zero  (see case 1)
    case INPAR::SCATRA::velocity_function:  // function
    {
      // we directly use the elements from the scalar transport elements section
      if (scatradis->NumGlobalNodes()==0)
        dserror("No elements in the ---TRANSPORT ELEMENTS section");

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for SCALAR_TRANSPORT problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly = rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn,false,"scatra",DRT::Problem::Instance()->SolverParams(linsolvernumber)));

      // read the restart information, set vectors and variables
      if (restart) scatraonly->ScaTraField().ReadRestart(restart);

      // set velocity field
      //(this is done only once. Time-dependent velocity fields are not supported)
      (scatraonly->ScaTraField()).SetVelocityField();

      // enter time loop to solve problem with given convective velocity
      (scatraonly->ScaTraField()).TimeLoop();

      // perform the result test if required
      DRT::Problem::Instance()->AddFieldTest(scatraonly->CreateScaTraFieldTest());
      DRT::Problem::Instance()->TestAll(comm);

      break;
    }
    case INPAR::SCATRA::velocity_Navier_Stokes:  // Navier_Stokes
    {
      // we use the fluid discretization as layout for the scalar transport discretization
      if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

      // create scatra elements if the scatra discretization is empty
      if (scatradis->NumGlobalNodes()==0)
      {
        DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,scatradis);
      }
      else
        dserror("Fluid AND ScaTra discretization present. This is not supported.");

      // we need a non-const list in order to be able to add sublists below!
      Teuchos::ParameterList prbdyn(scatradyn);
      // support for turbulent flow statistics
      const Teuchos::ParameterList& fdyn = (DRT::Problem::Instance()->FluidDynamicParams());
      prbdyn.sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for SCALAR_TRANSPORT problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create an one-way coupling algorithm instance
      Teuchos::RCP<SCATRA::PassiveScaTraAlgorithm> algo = Teuchos::rcp(new SCATRA::PassiveScaTraAlgorithm(comm,prbdyn,"scatra",DRT::Problem::Instance()->SolverParams(linsolvernumber)));

      // read restart information
      // in case a inflow generation in the inflow section has been performed, there are not any
      // scatra results available and the initial field is used
      if (restart)
      {
        if ((DRT::INPUT::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")==true) and
           (restart==fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
          algo->ReadInflowRestart(restart);
        else
          algo->ReadRestart(restart);
      }
      else
        if(DRT::INPUT::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")==true)
          dserror("Turbulent inflow generation for passive scalar transport should be performed as fluid problem!");

      // solve the whole (one-way-coupled) problem
      algo->TimeLoop();

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      DRT::Problem::Instance()->AddFieldTest(algo->FluidField().CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(algo->CreateScaTraFieldTest());
      DRT::Problem::Instance()->TestAll(comm);

      break;
    } // case 2
    default:
      dserror("unknown velocity field type for transport of passive scalar");
  }

  return;

} // end of scatra_dyn()

