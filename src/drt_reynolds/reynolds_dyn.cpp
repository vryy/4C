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

#include "reynolds_dyn.H"

#include "../drt_adapter/adapter_reynolds.H"

#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_scatra/scatra_timint_implicit.H"


/*----------------------------------------------------------------------*
 * Main control routine for Reynolds problems
 *----------------------------------------------------------------------*/
void reynolds_dyn(int restart)
{
  // access the communicator
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("reynolds")->Comm();

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

  // access the reynolds discretization
  Teuchos::RCP<DRT::Discretization> reynoldsdis =
      DRT::Problem::Instance()->GetDis("reynolds");

  reynoldsdis->FillComplete();

  // we directly use the elements from the Reynolds elements section
  if (reynoldsdis->NumGlobalNodes() == 0)
    dserror("No elements in the ---REYNOLDS ELEMENTS section");

  // add proxy of velocity related degrees of freedom to Reynolds discretization
  if (reynoldsdis->BuildDofSetAuxProxy(DRT::Problem::Instance()->NDim() + 1, 0, 0,
      true) != 1)
    dserror("Reynolds discretization has illegal number of dofsets!");

  // finalize discretization
  reynoldsdis->FillComplete(true, false, false);

  // get linear solver id from REYNOLDS DYNAMIC
  const int linsolvernumber = reynoldsdyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for REYNOLDS problem. Please set LINEAR_SOLVER in REYNOLDS DYNAMIC to a valid number!");

  // create instance of Reynolds basis algorithm
  Teuchos::RCP<ADAPTER::ReynoldsBaseAlgorithm> reynoldsonly = Teuchos::rcp(
      new ADAPTER::ReynoldsBaseAlgorithm());

  // setup Reynolds basis algorithm
  reynoldsonly->Setup(reynoldsdyn, scatradyn,
          DRT::Problem::Instance()->SolverParams(linsolvernumber));

  // set initial velocity field
  // note: The order ReadRestart() before SetVelocityField() is important here!!
  // for time-dependent velocity fields, SetVelocityField() is additionally called in each PrepareTimeStep()-call
  (reynoldsonly->ReynoldsField())->SetVelocityField(1);

  // read the restart information, set vectors and variables
  if (restart)
    reynoldsonly->ReynoldsField()->ReadRestart(restart);

  // enter time loop to solve problem
  (reynoldsonly->ReynoldsField())->TimeLoop();

  // perform the result test if required
  DRT::Problem::Instance()->AddFieldTest(reynoldsonly->CreateReynoldsFieldTest());
  DRT::Problem::Instance()->TestAll(comm);

return;

} // end of reynolds_dyn()
