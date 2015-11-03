/*--------------------------------------------------------------------------*/
/*!
\file reynolds_dyn.cpp

\brief entry point for the solution of Reynolds problems

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "reynolds_dyn.H"

/*----------------------------------------------------------------------*
 * Main control routine for Reynolds problems
 *----------------------------------------------------------------------*/
void reynolds_dyn(int restart)
{
  // access the communicator
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("scatra")->Comm();

  // print problem type
  if (comm.MyPID() == 0)
  {
    std::cout << "###################################################"
        << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: "
        << DRT::Problem::Instance()->ProblemName() << std::endl;
    std::cout << "###################################################"
        << std::endl;
  }

  // access the problem-specific parameter list
  const Teuchos::ParameterList& reynoldsdyn =
      DRT::Problem::Instance()->ReynoldsDynamicParams();
  const Teuchos::ParameterList& scatradyn =
      DRT::Problem::Instance()->ScalarTransportDynamicParams();

  // access the scatra discretization
  Teuchos::RCP<DRT::Discretization> scatradis =
      DRT::Problem::Instance()->GetDis("scatra");

  scatradis->FillComplete();

  // we directly use the elements from the scalar transport elements section
  if (scatradis->NumGlobalNodes() == 0)
    dserror("No elements in the ---TRANSPORT ELEMENTS section");

  // add proxy of velocity related degrees of freedom to scatra discretization
  if (scatradis->BuildDofSetAuxProxy(DRT::Problem::Instance()->NDim() + 1, 0, 0,
      true) != 1)
    dserror("Scatra discretization has illegal number of dofsets!");

  // finalize discretization
  scatradis->FillComplete(true, false, false);

  // get linear solver id from REYNOLDS DYNAMIC
  const int linsolvernumber = reynoldsdyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for REYNOLDS problem. Please set LINEAR_SOLVER in REYNOLDS DYNAMIC to a valid number!");

  // create instance of scalar transport basis algorithm (empty fluid discretization)
  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly = Teuchos::rcp(
      new ADAPTER::ScaTraBaseAlgorithm(reynoldsdyn, scatradyn,
          DRT::Problem::Instance()->SolverParams(linsolvernumber)));

  // set initial velocity field
  // note: The order ReadRestart() before SetVelocityField() is important here!!
  // for time-dependent velocity fields, SetVelocityField() is additionally called in each PrepareTimeStep()-call
  (scatraonly->ScaTraField())->SetVelocityField(1);

  // read the restart information, set vectors and variables
  if (restart)
    scatraonly->ScaTraField()->ReadRestart(restart);

  // enter time loop to solve problem
  (scatraonly->ScaTraField())->TimeLoop();

  // perform the result test if required
  DRT::Problem::Instance()->AddFieldTest(scatraonly->CreateScaTraFieldTest());
  DRT::Problem::Instance()->TestAll(comm);

return;

} // end of reynolds_dyn()
