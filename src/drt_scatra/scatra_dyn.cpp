/*----------------------------------------------------------------------*/
/*!
\file scatra_dyn.cpp
\brief entry point for scalar transport problems

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../drt_scatra_ele/scatra_ele.H"

#include <iostream>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "scatra_algorithm.H"
#include "scatra_resulttest.H"
#include "scatra_utils_clonestrategy.H"
#include "scatra_dyn.H"


/*----------------------------------------------------------------------*
 * Main control routine for scalar transport problems, incl. various solvers
 *
 *        o Laplace-/ Poisson equation (zero velocity field)
 *          (with linear and nonlinear boundary conditions)
 *        o transport of passive scalar in velocity field given by spatial function
 *        o transport of passive scalar in velocity field given by Navier-Stokes
 *          (one-way coupling)
 *        o scalar transport in velocity field given by Navier-Stokes with natural convection
 *          (two-way coupling)
 *
 *----------------------------------------------------------------------*/
void scatra_dyn(int restart)
{
  // access the communicator
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();

  // print problem type
  if (comm.MyPID() == 0)
  {
    std::cout << "###################################################"<< std::endl;
    std::cout << "# YOUR PROBLEM TYPE: " << DRT::Problem::Instance()->ProblemName() << std::endl;
    std::cout << "###################################################" << std::endl;
  }

  // access the problem-specific parameter list
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();

  // access the fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis("fluid");
  // access the scatra discretization
  Teuchos::RCP<DRT::Discretization> scatradis = DRT::Problem::Instance()->GetDis("scatra");

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
    case INPAR::SCATRA::velocity_function_and_curve: // function and time curve
    {
      // we directly use the elements from the scalar transport elements section
      if (scatradis->NumGlobalNodes()==0)
        dserror("No elements in the ---TRANSPORT ELEMENTS section");

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for SCALAR_TRANSPORT problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn,scatradyn,DRT::Problem::Instance()->SolverParams(linsolvernumber)));

      // read the restart information, set vectors and variables
      if (restart) scatraonly->ScaTraField()->ReadRestart(restart);

      // set initial velocity field
      // note: The order ReadRestart() before SetVelocityField() is important here!!
      // for time-dependent velocity fields, SetVelocityField() is additionally called in each PrepareTimeStep()-call
      (scatraonly->ScaTraField())->SetVelocityField();

      // enter time loop to solve problem with given convective velocity
      (scatraonly->ScaTraField())->TimeLoop();

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
        // fill scatra discretization by cloning fluid discretization
        DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,scatradis);

        // set implementation type of cloned scatra elements
        for(int i=0; i<scatradis->NumMyColElements(); ++i)
        {
          DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(scatradis->lColElement(i));
          if(element == NULL)
            dserror("Invalid element type!");
          else
            element->SetImplType(INPAR::SCATRA::impltype_std);
        }
      }

      else
        dserror("Fluid AND ScaTra discretization present. This is not supported.");

      // support for turbulent flow statistics
      const Teuchos::ParameterList& fdyn = (DRT::Problem::Instance()->FluidDynamicParams());

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        dserror("no linear solver defined for SCALAR_TRANSPORT problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

      // create a scalar transport algorithm instance
      Teuchos::RCP<SCATRA::ScaTraAlgorithm> algo = Teuchos::rcp(new SCATRA::ScaTraAlgorithm(comm,scatradyn,fdyn,"scatra",DRT::Problem::Instance()->SolverParams(linsolvernumber)));

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

      // solve the whole scalar transport problem
      algo->TimeLoop();

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      DRT::Problem::Instance()->AddFieldTest(algo->FluidField()->CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(algo->CreateScaTraFieldTest());
      DRT::Problem::Instance()->TestAll(comm);

      break;
    } // case 2
    default:
    {
      dserror("unknown velocity field type for transport of passive scalar");
      break;
    }
  }

  return;

} // end of scatra_dyn()

